# This file is part of MiniDOAS.jl (Differential Optical Absorption Spectroscopy in Julia).
# Copyright (c) 2026 Sunandan Mahant <sunandanmahant@outlook.com>; BSD-3-Clause, see LICENSE.

# =============================================================================
# Local web GUI: open a QDOAS project, browse its spectra, adjust the analysis
# windows, fit, map the fit window, and look at and export the results - in a web
# browser, served by this Julia session (gui/index.html is the page).
# =============================================================================

module MiniDOASGUIExt

using MiniDOAS
using MiniDOAS: ProjectSpec, apply_edits, save_edited_project, glob_regex, fmt_date, fmt_time, VERSION_STRING
using HTTP
using JSON3
using LinearAlgebra
using Printf
using Sockets

const GUI_PAGE = joinpath(@__DIR__, "..", "gui", "index.html")

"A background job (a batch of spectra, or a fit-window map), polled by the page."
mutable struct GuiJob
    running::Bool
    cancel::Bool
    done::Int
    total::Int
    started::Float64
    seconds::Float64
    result::Any
    error::String
end
GuiJob() = GuiJob(false, false, 0, 0, 0.0, 0.0, nothing, "")

mutable struct GuiState
    xml::String
    spec::Union{Nothing,ProjectSpec}
    jobs::Dict{String,GuiJob}
    # analysed spectra, per analysis settings (the page's window edits):
    # settings key => spectrum path => (file, result, error, seconds)
    results::Dict{String,Dict{String,Any}}
    result_keys::Vector{String}          # least recently used first
end
const GUI = GuiState("", nothing, Dict("batch" => GuiJob(), "map" => GuiJob()), Dict(), String[])
const GUI_LOCK = ReentrantLock()
const GUI_MAX_SETTINGS = 8

"Values for JSON: rounded to 6 significant digits (plenty for display), non-finite -> null."
jnum(x::Real) = isfinite(x) ? round(Float64(x); sigdigits=6) : nothing
jnum(v::AbstractVector) = [jnum(x) for x in v]

# ---------------------------------------------------------------------------
# analysis settings (the page's window edits)
# ---------------------------------------------------------------------------

"The edits as a canonical string (sorted keys), to recognise the same analysis settings."
function settings_key(edits)
    canon(x::AbstractDict) = "{" * join(("\"$k\":" * canon(x[k]) for k in sort!(collect(String.(keys(x))))), ",") * "}"
    canon(x::AbstractVector) = "[" * join(sort!(canon.(x)), ",") * "]"
    canon(x::Real) = string(Float64(x))
    canon(x) = JSON3.write(x)
    (edits === nothing || isempty(edits)) ? "{}" : canon(edits)
end

"Remembers analysed spectra under their analysis settings (a few settings are kept)."
function store_results!(key::String, items)
    lock(GUI_LOCK) do
        d = get!(() -> Dict{String,Any}(), GUI.results, key)
        filter!(!=(key), GUI.result_keys); push!(GUI.result_keys, key)
        while length(GUI.result_keys) > GUI_MAX_SETTINGS
            delete!(GUI.results, popfirst!(GUI.result_keys))
        end
        for it in items
            d[it.file] = it
        end
    end
end

"Stored results for the settings `key`, in spectrum order (optionally only `files`)."
function stored_results(key::String, files=nothing)
    lock(GUI_LOCK) do
        d = get(GUI.results, key, Dict{String,Any}())
        sel = files === nothing ? collect(keys(d)) : [f for f in files if haskey(d, f)]
        [d[f] for f in sort(sel)]
    end
end

# ---------------------------------------------------------------------------
# JSON views
# ---------------------------------------------------------------------------

function project_json(p::ProjectSpec, xml)
    shiftdesc(g) = Dict("symbols" => g.symbols,
                        "shift" => g.sh_fit ? "fitted" : @sprintf("fixed %.4g nm", g.sh_init),
                        "stretch" => g.st_type == 0 ? "none" : g.st_type == 1 ? "1st order, fitted" : "2nd order, fitted")
    offdesc(w) = join([(w.off[k][1] ? "order $(k-1) fitted" : "order $(k-1) = $(w.off[k][2])")
                       for k in 1:3 if w.off[k][1] || w.off[k][2] != 0], ", ")
    Dict("path" => xml, "name" => p.project_name, "interpolation" => string(p.interpolation),
         "convergence" => p.convergence, "pixels" => p.n_det, "calibration" => p.calib,
         "dark" => p.dark, "offset" => p.offset,
         "windows" => [Dict("name" => w.name, "lmin" => w.lmin, "lmax" => w.lmax, "lambda0" => w.lambda0,
                            "poly_order" => w.poly_order, "linear_offset" => w.offl_order,
                            "reference" => w.refone, "reference_ok" => isfile(w.refone),
                            "xs" => [Dict("sym" => c.name, "file" => c.file, "ok" => isfile(c.file),
                                          "use" => c.action === :none ? "as is" : "interpolated") for c in w.xs],
                            "shifts" => [shiftdesc(g) for g in w.shifts],
                            "offset" => offdesc(w), "gaps" => [collect(g) for g in w.gaps]) for w in p.windows],
         "spectra" => p.spectra)
end

function window_json(wr)
    cv = wr.curves
    poly = zeros(length(cv.lambda)); haspoly = false
    offl = zeros(length(cv.lambda)); hasoffl = false
    comps = Any[]
    for c in cv.components
        if c.kind === :poly
            poly .+= c.value; haspoly = true
        elseif c.kind === :offl
            offl .+= c.value; hasoffl = true
        else
            push!(comps, Dict("name" => c.name, "value" => jnum(c.value), "scd" => jnum(c.scd), "err" => jnum(c.err)))
        end
    end
    Dict("name" => wr.name, "rms" => jnum(wr.rms), "chi2" => jnum(wr.chisqr), "iterations" => wr.niter,
         "terms" => terms_json(wr),
         "shifts" => [Dict("sym" => s, "shift" => jnum(wr.shift[s]), "err" => jnum(wr.shift_err[s]),
                           "stretch" => jnum(get(wr.stretch, s, NaN)), "stretch_err" => jnum(get(wr.stretch_err, s, NaN)))
                      for s in sort(collect(keys(wr.shift)))],
         "offsets" => [Dict("name" => k, "value" => jnum(v), "err" => jnum(wr.param_err[k])) for (k, v) in sort(collect(wr.params))],
         "curves" => Dict("lambda" => jnum(cv.lambda), "od" => jnum(cv.od), "model" => jnum(cv.model),
                          "residual" => jnum(cv.residual), "fit_range" => collect(cv.fit_range),
                          "polynomial" => haspoly ? jnum(poly) : nothing,
                          "offset" => hasoffl ? jnum(offl) : nothing, "components" => comps))
end

terms_json(wr) = [Dict("term" => c, "scd" => jnum(wr.slcol[c]), "err" => jnum(wr.slerr[c])) for c in wr.columns]

function record_json(r)
    Dict("name" => r.name, "date" => fmt_date(r.date), "time" => fmt_time(r.time), "scans" => r.scans,
         "exposure" => jnum(Float64(r.exposure)), "elevation" => jnum(r.elevation), "azimuth" => jnum(r.azimuth),
         "latitude" => jnum(r.latitude), "longitude" => jnum(r.longitude))
end

function result_json(res, seconds)
    r = res.record
    Dict("spectrum" => res.spectrum, "name" => r.name, "date" => fmt_date(r.date), "time" => fmt_time(r.time),
         "record" => record_json(r), "seconds" => seconds,
         "observed" => Dict("lambda" => jnum(res.lambda), "counts" => jnum(res.raw)),
         "windows" => [window_json(wr) for wr in res.windows])
end

function table_json(items)
    header, rows = results_table(items)
    Dict{String,Any}("header" => header,
         "rows" => [[x isa Real ? jnum(x) : x === missing ? nothing : string(x) for x in r] for r in rows])
end

# ---------------------------------------------------------------------------
# HTTP plumbing
# ---------------------------------------------------------------------------

jresp(x; status=200) = HTTP.Response(status, ["Content-Type" => "application/json"], JSON3.write(x))
jerror(msg; status=400) = jresp(Dict("error" => msg); status=status)
const GUI_ALLOW_ANY_HOST = Ref(false)

"""
Requests must come from the GUI page: a Host of localhost/127.0.0.1 and, for POST, a
JSON body. This refuses cross-site form posts and DNS-rebinding tricks against the
local server. The work runs on a worker thread, off the HTTP thread.
"""
function guarded(f)
    function (req::HTTP.Request)
        host = lowercase(HTTP.header(req, "Host", ""))
        ok = GUI_ALLOW_ANY_HOST[] || any(h -> startswith(host, h), ("127.0.0.1", "localhost", "[::1]"))
        ok || return jerror("forbidden host"; status=403)
        if req.method == "POST" && !occursin("application/json", lowercase(HTTP.header(req, "Content-Type", "")))
            return jerror("expected application/json"; status=415)
        end
        try
            body = isempty(req.body) ? Dict{String,Any}() : JSON3.read(String(req.body), Dict{String,Any})
            fetch(Threads.@spawn f(body))
        catch e
            e = e isa TaskFailedException ? e.task.exception : e
            jerror(sprint(showerror, e))
        end
    end
end

current_project() = (GUI.spec === nothing && throw(ArgumentError("open a project first")); GUI.spec)
edits_of(b) = (e = get(b, "edits", nothing); e === nothing || isempty(e) ? nothing : e)

"Drives on Windows (those that can be read: an empty card reader errors), / elsewhere."
file_roots() = Sys.iswindows() ? [string(c, ":\\") for c in 'A':'Z' if readable_dir(string(c, ":\\"))] : ["/"]
readable_dir(p) = try isdir(p) catch; false end

function api_browse(b)
    path = expanduser(String(get(b, "path", "")))
    isempty(path) && (path = isempty(GUI.xml) ? homedir() : dirname(GUI.xml))
    isfile(path) && (path = dirname(path))
    isdir(path) || throw(ArgumentError("not a folder: $path"))
    path = abspath(path)
    dirs = String[]; files = Any[]
    for f in readdir(path)
        startswith(f, ".") && continue
        full = joinpath(path, f)
        try
            isdir(full) ? push!(dirs, f) : push!(files, Dict("name" => f, "size" => filesize(full)))
        catch
        end
    end
    parent = dirname(rstrip(path, ['/', '\\']))
    jresp(Dict("path" => path, "parent" => isempty(parent) ? path : parent, "dirs" => dirs, "files" => files,
               "roots" => file_roots(), "home" => homedir(), "sep" => Sys.iswindows() ? "\\" : "/"))
end

function api_project(b)
    xml = abspath(expanduser(String(b["path"])))
    isfile(xml) || throw(ArgumentError("no such file: $xml"))
    p = parse_config(xml)
    lock(GUI_LOCK) do
        GUI.xml = xml; GUI.spec = p
        empty!(GUI.results); empty!(GUI.result_keys)
    end
    jresp(project_json(p, xml))
end

function api_folder(b)
    dir = expanduser(String(b["dir"]))
    rx = glob_regex(String(get(b, "filter", "*")))
    isdir(dir) || throw(ArgumentError("not a folder: $dir"))
    found = String[]
    if get(b, "recursive", false) == true
        for (root, _, fs) in walkdir(dir), f in fs
            occursin(rx, f) && push!(found, joinpath(root, f))
        end
    else
        append!(found, [joinpath(dir, f) for f in readdir(dir) if occursin(rx, f) && isfile(joinpath(dir, f))])
    end
    jresp(Dict("spectra" => sort(found)))
end

function api_analyze(b)
    edits = edits_of(b)
    p = apply_edits(current_project(), edits)
    spectrum = String(b["spectrum"])
    t0 = time()
    res = try
        analyze(p; spectrum=spectrum, curves=true)
    catch e
        store_results!(settings_key(edits), [(file=spectrum, result=nothing, error=sprint(showerror, e), seconds=time() - t0)])
        rethrow()
    end
    t = time() - t0
    store_results!(settings_key(edits), [(file=spectrum, result=res, error="", seconds=t)])
    d = result_json(res, t)
    # with edited settings, also the fit with the project's own settings, for comparison
    if edits !== nothing && get(b, "baseline", false) == true
        try
            r0 = analyze(current_project(); spectrum=spectrum)
            d["baseline"] = Dict(wr.name => Dict("rms" => jnum(wr.rms), "terms" => terms_json(wr)) for wr in r0.windows)
        catch
        end
    end
    jresp(d)
end

function start_job(f, name, total)
    job = GUI.jobs[name]
    lock(GUI_LOCK) do
        job.running && throw(ArgumentError("a $name is already running"))
        job.running = true; job.cancel = false; job.done = 0; job.total = total
        job.started = time(); job.result = nothing; job.error = ""
    end
    Threads.@spawn try
        job.result = f(job)
    catch e
        job.error = e isa InterruptException ? "cancelled" : sprint(showerror, e)
    finally
        job.seconds = time() - job.started
        job.running = false
    end
    jresp(Dict("started" => true, "total" => total))
end

job_status(name) = (j = GUI.jobs[name];
    Dict{String,Any}("running" => j.running, "done" => j.done, "total" => j.total, "error" => j.error,
                     "seconds" => j.running ? time() - j.started : j.seconds))

function api_batch_start(b)
    edits = edits_of(b)
    p = apply_edits(current_project(), edits)
    files = String.(collect(get(b, "spectra", p.spectra)))
    isempty(files) && throw(ArgumentError("no spectra to analyse"))
    key = settings_key(edits)
    start_job("batch", length(files)) do job
        # in chunks, so that the job can be cancelled and results appear as they come
        for chunk in Iterators.partition(files, max(16, 4 * Threads.nthreads()))
            job.cancel && throw(InterruptException())
            bt = analyze_spectra(p, collect(chunk))
            store_results!(key, bt)
            job.done += length(chunk)
        end
        length(files)
    end
end

api_batch_status(_) = jresp(job_status("batch"))

function api_cancel(b)
    GUI.jobs[String(b["job"])].cancel = true
    jresp(Dict("ok" => true))
end

function api_results(b)
    files = haskey(b, "spectra") ? String.(collect(b["spectra"])) : nothing
    items = stored_results(settings_key(edits_of(b)), files)
    d = table_json(items)
    d["seconds"] = sum((x.seconds for x in items); init=0.0)
    jresp(d)
end

function api_results_csv(b)
    files = haskey(b, "spectra") ? String.(collect(b["spectra"])) : nothing
    items = stored_results(settings_key(edits_of(b)), files)
    isempty(items) && throw(ArgumentError("nothing analysed yet with the current settings"))
    io = IOBuffer()
    path = tempname() * ".csv"
    try
        write_csv(path, items)
        write(io, read(path))
    finally
        rm(path; force=true)
    end
    name = replace(current_project().project_name, r"[^\w.-]+" => "_")
    HTTP.Response(200, ["Content-Type" => "text/csv",
                        "Content-Disposition" => "attachment; filename=\"$(name)_results.csv\""], take!(io))
end

"Retrieval interval mapping (`fit_window_map`) as a background job."
function api_map_start(b)
    edits = something(edits_of(b), Dict{String,Any}())
    base = current_project()
    wname = String(b["window"])
    any(w -> w.name == wname, base.windows) || throw(ArgumentError("no analysis window $wname"))
    spectrum = String(b["spectrum"])
    starts = collect(Float64(b["start_min"]):Float64(b["step"]):Float64(b["start_max"]) + 1e-9)
    ends = collect(Float64(b["end_min"]):Float64(b["step"]):Float64(b["end_max"]) + 1e-9)
    minwidth = Float64(get(b, "min_width", 5.0))
    length(starts) * length(ends) <= 10_000 || throw(ArgumentError("too many cells; use a larger step"))
    any(ends[j] - starts[i] >= minwidth for i in eachindex(starts), j in eachindex(ends)) ||
        throw(ArgumentError("no window of at least $minwidth nm in these ranges"))
    # only this window, with the page's other settings for it
    only_w = Dict{String,Any}(w.name => Dict{String,Any}("enabled" => false) for w in base.windows if w.name != wname)
    only_w[wname] = Dict{String,Any}(get(edits, wname, Dict{String,Any}()))
    p = apply_edits(base, only_w)
    total = length(starts) * length(ends)
    start_job("map", total) do job
        m = fit_window_map(p, spectrum; window=wname, starts=starts, ends=ends, min_width=minwidth,
                           on_progress=(k, _) -> (job.done = k; !job.cancel))
        job.cancel && throw(InterruptException())
        cell(i, j) = ends[j] - starts[i] < minwidth ? nothing :
            (wr = m.results[i, j]; wr === nothing ? Dict("error" => m.errors[i, j]) :
             Dict("rms" => jnum(wr.rms), "chi2" => jnum(wr.chisqr), "iterations" => wr.niter,
                  "scd" => Dict(c => jnum(wr.slcol[c]) for c in wr.columns),
                  "err" => Dict(c => jnum(wr.slerr[c]) for c in wr.columns)))
        Dict("window" => wname, "spectrum" => spectrum, "starts" => starts, "ends" => ends,
             "cells" => [[cell(i, j) for j in eachindex(ends)] for i in eachindex(starts)])
    end
end

function api_map_status(_)
    d = job_status("map")
    j = GUI.jobs["map"]
    !j.running && j.result !== nothing && (d["map"] = j.result)
    jresp(d)
end

function api_netcdf(b)
    p = apply_edits(current_project(), edits_of(b))
    spectrum = String(b["spectrum"])
    res = analyze(p; spectrum=spectrum)
    path = tempname() * ".nc"
    write_netcdf(path, res; config_path=GUI.xml)
    body = read(path); rm(path; force=true)
    HTTP.Response(200, ["Content-Type" => "application/x-netcdf",
                        "Content-Disposition" => "attachment; filename=\"$(splitext(basename(spectrum))[1]).nc\""], body)
end

function api_save(b)
    isempty(GUI.xml) && throw(ArgumentError("open a project first"))
    dst = abspath(expanduser(String(b["path"])))
    dst == GUI.xml && throw(ArgumentError("choose a new file name; the original project is kept as it is"))
    endswith(lowercase(dst), ".xml") || (dst *= ".xml")
    isdir(dirname(dst)) || throw(ArgumentError("no such folder: $(dirname(dst))"))
    isfile(dst) && get(b, "overwrite", false) != true && return jresp(Dict("exists" => dst))
    save_edited_project(GUI.xml, dst, edits_of(b))
    parse_config(dst)                                   # the saved file must read back
    jresp(Dict("saved" => dst))
end

"Opens `url` in the default browser (best effort)."
function open_in_browser(url)
    try
        if Sys.iswindows()
            run(`cmd /c start "" $url`; wait=false)
        elseif Sys.isapple()
            run(`open $url`; wait=false)
        else
            run(pipeline(`xdg-open $url`; stdout=devnull, stderr=devnull); wait=false)
        end
    catch
    end
end

"The first port from `port` on that is free on `host`."
function free_port(host, port)
    for p in port:port+20
        try
            s = Sockets.listen(Sockets.getaddrinfo(host), p)
            close(s)
            return p
        catch
        end
    end
    port
end

"""
    serve_gui(; project="", port=8765, host="127.0.0.1", open_browser=true)

Starts the web GUI (gui/index.html) and serves it until interrupted with Ctrl-C. If
`port` is taken, the next free port is used. Start Julia with several threads
(`julia -t auto`) so that batches and fit-window maps run in parallel.

The GUI can read any file this Julia process can read and writes the files you ask
it to, so by default it only listens on this computer (127.0.0.1). Pass
`host="0.0.0.0"` only on a network you trust.
"""
function serve_gui(; project::AbstractString="", port::Integer=8765, host::AbstractString="127.0.0.1",
                   open_browser::Bool=true)
    isfile(GUI_PAGE) || error("GUI page not found: $GUI_PAGE")
    GUI_ALLOW_ANY_HOST[] = !(host in ("127.0.0.1", "localhost", "::1"))
    isempty(project) || api_project(Dict{String,Any}("path" => project))
    # compile the fit and its JSON once in the background, so that the first fit in the
    # page is as fast as the others
    GUI.spec !== nothing && !isempty(GUI.spec.spectra) && Threads.@spawn try
        p = GUI.spec
        w = first(p.windows)
        e = Dict{String,Any}(w.name => Dict{String,Any}("lmin" => w.lmin + 0.1))
        JSON3.write(result_json(analyze(apply_edits(p, e); spectrum=first(p.spectra), curves=true), 0.0))
        table_json([(file=first(p.spectra), result=analyze(p; spectrum=first(p.spectra)), error="", seconds=0.0)])
    catch
    end
    r = HTTP.Router()
    HTTP.register!(r, "GET", "/", _ -> HTTP.Response(200, ["Content-Type" => "text/html; charset=utf-8",
                                                          "Cache-Control" => "no-store"], read(GUI_PAGE)))
    HTTP.register!(r, "GET", "/api/state", guarded(_ -> jresp(Dict("version" => VERSION_STRING,
        "threads" => Threads.nthreads(), "home" => homedir(), "julia" => string(VERSION),
        "examples" => abspath(joinpath(@__DIR__, "..", "examples")),
        "example_projects" => filter(isfile, [abspath(joinpath(@__DIR__, "..", "examples", f)) for f in
            ("ace_maxdoas/work/ace_maxdoas.xml", "synthetic/synthetic_project.xml")]),
        "project" => GUI.spec === nothing ? nothing : project_json(GUI.spec, GUI.xml)))))
    HTTP.register!(r, "POST", "/api/browse", guarded(api_browse))
    HTTP.register!(r, "POST", "/api/project", guarded(api_project))
    HTTP.register!(r, "POST", "/api/folder", guarded(api_folder))
    HTTP.register!(r, "POST", "/api/analyze", guarded(api_analyze))
    HTTP.register!(r, "POST", "/api/batch", guarded(api_batch_start))
    HTTP.register!(r, "GET", "/api/batch", guarded(api_batch_status))
    HTTP.register!(r, "POST", "/api/cancel", guarded(api_cancel))
    HTTP.register!(r, "POST", "/api/results", guarded(api_results))
    HTTP.register!(r, "POST", "/api/results.csv", guarded(api_results_csv))
    HTTP.register!(r, "POST", "/api/map", guarded(api_map_start))
    HTTP.register!(r, "GET", "/api/map", guarded(api_map_status))
    HTTP.register!(r, "POST", "/api/netcdf", guarded(api_netcdf))
    HTTP.register!(r, "POST", "/api/save", guarded(api_save))
    port = free_port(host, port)
    shown = host in ("0.0.0.0", "::") ? "127.0.0.1" : host
    url = "http://$shown:$port/"
    println("MiniDOAS.jl $VERSION_STRING GUI at $url  ($(Threads.nthreads()) threads; Ctrl-C to stop)")
    flush(stdout)
    open_browser && open_in_browser(url)
    HTTP.serve(r, host, port)
end

end # module
