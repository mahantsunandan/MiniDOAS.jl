# This file is part of MiniDOAS.jl (Differential Optical Absorption Spectroscopy in Julia).
# Copyright (c) 2026 Sunandan Mahant <sunandanmahant@outlook.com>
# Derived from QDOAS, Copyright (C) 1994-2025 BIRA-IASB and S[&]T (BSD-3-Clause).
# See LICENSE and NOTICE.md.

# =============================================================================
# A fit set up once (windows with their cross sections on the analysis grid and their
# references) and applied to any number of spectra: ANALYSE_Spectrum over all windows.
# =============================================================================

"""
    DOASResult

The analysis of one spectrum, from `doas_fit` or `analyze`:

  * `windows`   one `WindowResult` per analysis window (`r["NO2"]` or `r[1]` picks one)
  * `lambda`    the wavelength grid, nm
  * `raw`       the spectrum analysed (after dark and offset correction for MFC files)
  * `record`    the `MFCSpectrum` it came from (name, date, time, angles), or `nothing`
  * `spectrum`  its file, or the name given to `doas_fit`
  * `spec`      the settings used (`ProjectSpec`)
"""
struct DOASResult
    spec::ProjectSpec
    spectrum::String
    lambda::Vector{Float64}
    windows::Vector{WindowResult}
    record::Union{Nothing,MFCSpectrum}
    raw::Vector{Float64}
end

function Base.getindex(r::DOASResult, name::AbstractString)
    i = findfirst(w -> w.name == name, r.windows)
    i === nothing && throw(KeyError(name))
    r.windows[i]
end
Base.getindex(r::DOASResult, i::Integer) = r.windows[i]
Base.keys(r::DOASResult) = [w.name for w in r.windows]

function Base.show(io::IO, r::DOASResult)
    print(io, "DOASResult(", isempty(r.spectrum) ? "spectrum" : repr(basename(r.spectrum)), ", ",
          join((@sprintf("%s RMS %.2e", w.name, w.rms) for w in r.windows), ", "), ")")
end
Base.show(io::IO, ::MIME"text/plain", r::DOASResult) = print_summary(io, r)

"""
    DOASModel(lambda, reference, windows; interpolation=:spline, security_gap=10,
              convergence=1e-4, max_iterations=0, name="MiniDOAS.jl")
    DOASModel(p::ProjectSpec)

A DOAS fit set up once for any number of spectra: the analysis windows (`FitWindow`s), with
their cross sections interpolated onto the wavelength grid `lambda` (nm, one value per
detector pixel), and the Fraunhofer `reference` spectrum I₀ on that grid. Fit spectra with
`doas_fit(model, spectrum)`.

  * `interpolation`: `:spline` (cubic) or `:linear`, for cross sections and shifted vectors;
  * `security_gap`: pixels kept on both sides of each window for shifting vectors;
  * `convergence`: relative change of χ² at which the Levenberg-Marquardt fit stops;
  * `max_iterations`: limit on the fit's iterations (0: none).

From a QDOAS project (`parse_config`), each window gets the reference and wavelengths of
its own reference file.

A model can be shared by threads: `doas_fit` keeps its working memory per call.
"""
struct DOASModel
    spec::ProjectSpec
    lambda::Vector{Float64}
    windows::Vector{Window}
    refs::Vector{Vector{Float64}}      # normalised references, per window
    ref_norms::Vector{Float64}
    xs_cache::Dict{Any,Any}            # in-memory cross sections on the grid
    reference::Vector{Float64}         # the reference as given (empty for a QDOAS project)
end

"The output fields a model made from Julia writes to QDOAS-layout NetCDF."
const DEFAULT_FIELDS = ["name", "date", "time", "rms", "chi", "iter_number", "error_flag", "residual_spectrum"]

function DOASModel(p::ProjectSpec)
    n = p.n_det
    windows = Window[]; refs = Vector{Float64}[]; norms = Float64[]
    xs_cache = Dict{Any,Any}()
    lam_ref = Dict{String,Tuple{Vector{Float64},Vector{Float64},Float64}}()
    lam = Float64[]
    for ws in p.windows
        lr, ref, ref_norm = get!(lam_ref, ws.refone) do
            cached((ws.refone, :reference, n)) do
                l, v = load_two_columns(ws.refone, n)
                v2 = copy(v)
                nrm = normalize!(v2)
                (l, v2, nrm)
            end
        end
        lam = lr
        push!(windows, build_window(ws, p, lr, xs_cache))
        push!(refs, ref); push!(norms, ref_norm)
    end
    DOASModel(p, lam, windows, refs, norms, xs_cache, Float64[])
end

function DOASModel(lambda::AbstractVector{<:Real}, reference::AbstractVector{<:Real}, windows;
                   interpolation::Symbol=:spline, security_gap::Integer=10, convergence::Real=1e-4,
                   max_iterations::Integer=0, name::AbstractString="MiniDOAS.jl")
    lam = Vector{Float64}(lambda)
    n = length(lam)
    length(reference) == n || throw(DimensionMismatch("$n wavelengths, reference of $(length(reference))"))
    n >= 2 && all(>(0), diff(lam)) || throw(ArgumentError("wavelengths must increase"))
    interpolation in (:spline, :linear) || throw(ArgumentError("interpolation must be :spline or :linear"))
    security_gap >= 0 || throw(ArgumentError("security_gap must be >= 0"))
    ws = windows isa FitWindow ? [windows] : collect(FitWindow, windows)
    isempty(ws) && throw(ArgumentError("no analysis window"))
    allunique([w.name for w in ws]) || throw(ArgumentError("analysis window names must be unique"))
    p = ProjectSpec(interpolation, security_gap, convergence, max_iterations, 999.9, n, false,
                    "", "", "", "", String[], ws, copy(DEFAULT_FIELDS), "MM/DD/YYYY", String(name), "QDOAS Results")
    raw_ref = Vector{Float64}(reference)
    ref = copy(raw_ref)
    nrm = normalize!(ref)
    xs_cache = Dict{Any,Any}()
    built = [build_window(w, p, lam, xs_cache) for w in ws]
    DOASModel(p, lam, built, [ref for _ in ws], fill(nrm, length(ws)), xs_cache, raw_ref)
end

function Base.show(io::IO, m::DOASModel)
    print(io, "DOASModel(", length(m.lambda), " pixels, ", @sprintf("%.2f-%.2f nm", first(m.lambda), last(m.lambda)),
          ", windows ", join((w.spec.name for w in m.windows), ", "), ")")
end

"""
    doas_fit(model, spectrum; curves=false) -> DOASResult
    doas_fit(model, spectra; threads=true, curves=false, on_progress=nothing) -> DOASBatch
    doas_fit(lambda, spectrum, reference, windows; kw...)

Fits spectra with every analysis window of `model`.

**One spectrum** gives a `DOASResult`. It can be

  * a file: an MFC STD spectrum, or a text file of numbers (one column of intensities, or
    wavelength and intensity: the second column is used), one value per pixel of the model;
  * an `MFCSpectrum` (from `read_mfc_std`), or a vector of intensities on the model's grid.

**Several spectra** give a `DOASBatch`: a vector of files, of `MFCSpectrum`s or of intensity
vectors; a matrix with one spectrum per column; or a string naming a folder or a pattern,
e.g. `"spectra/*.STD"` (see `spectrum_files`; `filter` and `recursive` are passed on). They
are fitted in parallel when Julia runs with several threads; a spectrum that cannot be
fitted gets `result = nothing` and its error message, and the others are unaffected.

For MFC STD files, `dark` and `offset` (files or `MFCSpectrum`s) are removed first as QDOAS
does (`correct_dark_offset`). With `curves=true` each `WindowResult` also holds the fitted
curves, for `plot_fit`. A spectrum's file (or `name`) becomes the result's `spectrum`, and an
MFC spectrum's name, date, time and angles go into the result's `record`.

The last form sets up a `DOASModel(lambda, reference, windows; kw...)` and fits `spectrum`.

Throws `EngineError` where QDOAS would fail (e.g. a non-positive intensity inside a
window, or a spectrum identical to the reference).

```julia
r = doas_fit(model, "spectra/S0000001.STD")                 # one file
b = doas_fit(model, "spectra/*.STD"; dark="dark.STD")       # every matching file
b = doas_fit(model, ["a.STD", "b.STD"])                     # a list of files
slant_column(b, "NO2")                                      # NO2 and its error, per spectrum
```
"""
function doas_fit(m::DOASModel, spectrum::AbstractVector{<:Real}; name::AbstractString="",
                  record::Union{Nothing,MFCSpectrum}=nothing, curves::Bool=false)
    p = m.spec
    n = p.n_det
    length(spectrum) == n || throw(DimensionMismatch("spectrum of $(length(spectrum)) pixels, the model has $n"))
    p = name == p.spectrum ? p : respec(p; spectrum=String(name))
    spec = Vector{Float64}(spectrum)
    raw = copy(spec)
    spe_norm = normalize!(spec)
    results = WindowResult[]
    for (w, ref, ref_norm) in zip(m.windows, m.refs, m.ref_norms)
        # ANALYSE_Spectrum skips a window whose spectrum equals its reference
        if all(i -> (ref[i] == 0 ? spec[i] == ref[i] : abs(spec[i] - ref[i]) / abs(ref[i]) <= 1e-7), 1:n)
            throw(EngineError("spectrum identical to reference in $(w.spec.name)"))
        end
        push!(results, window_result(w, p, spec, ref; spe_norm=spe_norm, ref_norm=ref_norm, curves=curves))
    end
    DOASResult(p, p.spectrum, m.lambda, results, record, raw)
end

doas_fit(m::DOASModel, s::MFCSpectrum; name::AbstractString=s.name, curves::Bool=false) =
    doas_fit(m, s.counts; name=name, record=s, curves=curves)

"A dark or offset spectrum given as a file or an `MFCSpectrum` (or nothing)."
as_mfc(x::Nothing) = nothing
as_mfc(x::MFCSpectrum) = x
as_mfc(x::AbstractString) = read_mfc_std(x)

"Reads, corrects and fits one spectrum file."
function fit_file(m::DOASModel, path::AbstractString, dark, offset, curves::Bool)
    counts, rec = read_spectrum_file(path; date_format=m.spec.date_format)
    n = m.spec.n_det
    length(counts) == n || throw(DimensionMismatch("$path has $(length(counts)) pixels, the model $n"))
    if dark !== nothing || offset !== nothing
        rec === nothing && throw(ArgumentError("dark and offset corrections need MFC STD spectra ($path)"))
        counts = correct_dark_offset(rec; dark=dark, offset=offset)
    end
    doas_fit(m, counts; name=String(path), record=rec, curves=curves)
end

function doas_fit(m::DOASModel, path::AbstractString; dark=nothing, offset=nothing, curves::Bool=false,
                  threads::Bool=true, on_progress=nothing, filter::AbstractString="*", recursive::Bool=false)
    names_many(path) && return doas_fit(m, spectrum_files(path; filter=filter, recursive=recursive); dark=dark,
                                        offset=offset, curves=curves, threads=threads, on_progress=on_progress)
    fit_file(m, path, as_mfc(dark), as_mfc(offset), curves)
end

function doas_fit(m::DOASModel, files::AbstractVector{<:AbstractString}; dark=nothing, offset=nothing,
                  threads::Bool=true, curves::Bool=false, on_progress=nothing)
    d, o = as_mfc(dark), as_mfc(offset)
    run_batch(String.(files); threads=threads, on_progress=on_progress) do i
        fit_file(m, files[i], d, o, curves)
    end
end

doas_fit(m::DOASModel, spectra::AbstractMatrix{<:Real}; kw...) = doas_fit(m, collect(eachcol(spectra)); kw...)

function doas_fit(m::DOASModel, spectra::AbstractVector; threads::Bool=true, curves::Bool=false, on_progress=nothing)
    label(s, i) = s isa MFCSpectrum && !isempty(s.name) ? s.name : s isa AbstractString ? s : "spectrum $i"
    run_batch([label(s, i) for (i, s) in enumerate(spectra)]; threads=threads, on_progress=on_progress) do i
        doas_fit(m, spectra[i]; curves=curves)
    end
end

doas_fit(lambda::AbstractVector{<:Real}, spectrum, reference::AbstractVector{<:Real}, windows;
         curves::Bool=false, kw...) = doas_fit(DOASModel(lambda, reference, windows; kw...), spectrum; curves=curves)

# =============================================================================
# Batches
# =============================================================================

"One spectrum of a batch."
const BatchItem = NamedTuple{(:file, :result, :error, :seconds),Tuple{String,Union{Nothing,DOASResult},String,Float64}}

"""
    DOASBatch

The analysis of several spectra (`analyze_spectra`, `analyze` or `doas_fit` given several
files or spectra): a vector with one `(file, result, error, seconds)` per spectrum, in the
order given. `result` is a `DOASResult`, or `nothing` - and `error` says why - for a
spectrum that could not be analysed.

`batch[i]` or `batch["S0000001.STD"]` (a file, or just its name) picks a spectrum;
`slant_column(batch, "NO2")` gives a slant column and its error per spectrum,
`measurement_time(batch)` their times; `results_columns(batch)` is the whole table (e.g.
for `DataFrame`), `write_csv(path, batch)` writes it, and `write_netcdf(folder, batch)` writes
one QDOAS-layout NetCDF file per spectrum. `batch.seconds` is the time the batch took.
"""
struct DOASBatch <: AbstractVector{BatchItem}
    items::Vector{BatchItem}
    seconds::Float64
end

Base.size(b::DOASBatch) = size(b.items)
Base.IndexStyle(::Type{DOASBatch}) = IndexLinear()
Base.getindex(b::DOASBatch, i::Int) = b.items[i]
function Base.getindex(b::DOASBatch, file::AbstractString)
    i = findfirst(x -> x.file == file, b.items)
    i === nothing && (i = findfirst(x -> basename(x.file) == file, b.items))
    i === nothing && throw(KeyError(file))
    b.items[i]
end

Base.show(io::IO, b::DOASBatch) =
    print(io, "DOASBatch(", length(b), " spectra, ", count(x -> x.result !== nothing, b.items), " analysed)")

"A one-line summary of a window: its first term with error, and the RMS."
function window_line(w::WindowResult)
    isempty(w.columns) && return @sprintf("%s: RMS %.2e", w.name, w.rms)
    c = first(w.columns)
    @sprintf("%s: %s %.3e ± %.1e, RMS %.2e", w.name, c, w.slcol[c], w.slerr[c], w.rms)
end

function Base.show(io::IO, ::MIME"text/plain", b::DOASBatch)
    nok = count(x -> x.result !== nothing, b.items)
    @printf(io, "DOASBatch: %d spectra, %d analysed, %d failed", length(b), nok, length(b) - nok)
    b.seconds > 0 && @printf(io, " in %.2f s", b.seconds)
    for (k, x) in enumerate(b.items)
        if k > 12
            print(io, "\n  ... ", length(b) - 12, " more")
            break
        end
        name = rpad(isempty(x.file) ? "spectrum" : basename(x.file), 24)
        if x.result === nothing
            print(io, "\n  ", name, " failed: ", x.error)
        else
            print(io, "\n  ", name, " ", join((window_line(w) for w in x.result.windows), " | "))
        end
    end
end

"""
    slant_column(result, window, term=window) -> (value, error)
    slant_column(batch, window, term=window) -> (values, errors)

The slant column of `term` (a cross section, or a polynomial or offset term such as `"x0"`)
fitted in the analysis window `window`, with its 1σ error. For a batch, vectors in the
batch's order, NaN where a spectrum failed. `slant_column(batch, "NO2")` is NO₂ in the NO2
window; `slant_column(batch, "O4", "NO2")` NO₂ in the O4 window.
"""
function slant_column(r::DOASResult, window::AbstractString, term::AbstractString=window)
    w = r[window]
    haskey(w.slcol, term) || throw(ArgumentError("window $window has no term $term; it has $(join(w.columns, ", "))"))
    w.slcol[term], w.slerr[term]
end

function slant_column(b::AbstractVector, window::AbstractString, term::AbstractString=window)
    v = fill(NaN, length(b)); e = fill(NaN, length(b))
    found = false
    for (i, x) in enumerate(b)
        x.result === nothing && continue
        k = findfirst(w -> w.name == window, x.result.windows)
        (k === nothing || !haskey(x.result.windows[k].slcol, term)) && continue
        w = x.result.windows[k]
        v[i], e[i] = w.slcol[term], w.slerr[term]
        found = true
    end
    found || all(x -> x.result === nothing, b) ||
        throw(ArgumentError("no analysed spectrum has the term $term in a window $window"))
    v, e
end

"""
    measurement_time(result) -> DateTime
    measurement_time(batch) -> Vector{Union{DateTime,Missing}}

The mid-exposure time of an analysed MFC STD spectrum (from its `record`), or of every
spectrum of a batch (`missing` where it failed or has no record).
"""
measurement_time(r::DOASResult) = r.record === nothing ?
    throw(ArgumentError("the result has no spectrum record (date and time)")) : measurement_time(r.record)
measurement_time(b::AbstractVector) =
    Union{DateTime,Missing}[x.result === nothing || x.result.record === nothing ? missing : measurement_time(x.result.record)
                            for x in b]

"""
Runs `f(i)` for every item, in parallel when Julia has several threads (BLAS set to one
thread meanwhile), and returns `(file=labels[i], result, error, seconds)` per item.
"""
function run_batch(f, labels::AbstractVector{<:AbstractString}; threads::Bool=true, on_progress=nothing)
    t00 = time()
    out = Vector{BatchItem}(undef, length(labels))
    done = Threads.Atomic{Int}(0)
    lk = ReentrantLock()
    function work(i)
        t0 = time()
        out[i] = try
            BatchItem((String(labels[i]), f(i), "", time() - t0))
        catch e
            e isa InterruptException && rethrow()
            BatchItem((String(labels[i]), nothing, sprint(showerror, e), time() - t0))
        end
        k = Threads.atomic_add!(done, 1) + 1
        on_progress === nothing || lock(() -> on_progress(k, length(labels), labels[i]), lk)
        nothing
    end
    parallel = threads && Threads.nthreads() > 1 && length(labels) > 1
    nblas = BLAS.get_num_threads()
    parallel && BLAS.set_num_threads(1)
    try
        if parallel
            Threads.@threads :dynamic for i in eachindex(labels)
                work(i)
            end
        else
            foreach(work, eachindex(labels))
        end
    finally
        parallel && BLAS.set_num_threads(nblas)
    end
    DOASBatch(out, time() - t00)
end

# =============================================================================
# QDOAS projects
# =============================================================================

"""
    analyze(p::ProjectSpec; spectrum=p.spectrum, curves=false) -> DOASResult
    analyze(xml_path; spectrum=..., curves=false)

Runs every enabled analysis window of a QDOAS project on one spectrum (an MFC STD file,
dark- and offset-corrected as the project says) and returns a `DOASResult`. With
`curves=true` each `WindowResult` also carries the fitted curves (`.curves`) for plotting.
Throws `UnsupportedConfig` or `EngineError` (see their docstrings).
"""
function analyze(p::ProjectSpec; spectrum::AbstractString=p.spectrum, curves::Bool=false)
    isempty(spectrum) && throw(ArgumentError("no spectrum: the project lists none, pass spectrum=\"...\""))
    p = spectrum == p.spectrum ? p : respec(p; spectrum=String(spectrum))
    spec, record = load_spectrum(p)
    doas_fit(DOASModel(p), spec; name=p.spectrum, record=record, curves=curves)
end

"""
    analyze(p, spectrum; curves=false) -> DOASResult
    analyze(p, spectra; threads=true, curves=false, on_progress=nothing) -> DOASBatch

The same with the spectra given directly: one file gives a `DOASResult`; a vector of
files, a folder or a pattern such as `"spectra/*.STD"` (see `spectrum_files`; `filter` and
`recursive` are passed on) give a `DOASBatch` (as `analyze_spectra`). `p` is a project or
the path of a project file.

```julia
r = analyze("station.xml", "spectra/S0000001.STD")     # one spectrum
b = analyze("station.xml", "spectra/*.STD")            # all of them, in parallel
```
"""
function analyze(p::ProjectSpec, spectrum::AbstractString; curves::Bool=false, threads::Bool=true,
                 on_progress=nothing, filter::AbstractString="*", recursive::Bool=false)
    names_many(spectrum) || return analyze(p; spectrum=spectrum, curves=curves)
    analyze_spectra(p, spectrum_files(spectrum; filter=filter, recursive=recursive);
                    threads=threads, curves=curves, on_progress=on_progress)
end

analyze(p::ProjectSpec, files::AbstractVector{<:AbstractString}; kw...) = analyze_spectra(p, files; kw...)
analyze(xml_path::AbstractString; kw...) = analyze(parse_config(xml_path); kw...)
analyze(xml_path::AbstractString, spectra; kw...) = analyze(parse_config(xml_path), spectra; kw...)

"""
    analyze_config(xml_path; spectrum=nothing, curves=false)

Parses a QDOAS project and analyses one spectrum: the given one, or the project's
first. Same result as `analyze`.
"""
function analyze_config(xml_path::AbstractString; spectrum=nothing, curves::Bool=false)
    p = parse_config(xml_path)
    analyze(p; spectrum=something(spectrum, p.spectrum), curves=curves)
end

"""
    analyze_spectra(p, files=project_spectra(p); threads=true, curves=false, on_progress=nothing)
        -> DOASBatch

Analyses many spectra with one project: by default those the project lists; or `files`, a
vector of files, or a folder or pattern such as `"spectra/*.STD"` (see `spectrum_files`). A
spectrum that cannot be analysed gets `result = nothing` and its error message; the others
are unaffected. Results come back in the order of `files`, as a `DOASBatch`.

When Julia runs with several threads (`julia -t auto`) and `threads=true`, spectra are
fitted in parallel; OpenBLAS is set to one thread meanwhile, since many small fits
at once are faster that way (results are unchanged to ~1e-16). `on_progress(done,
total, file)` is called after each spectrum, one call at a time.
"""
function analyze_spectra(p::ProjectSpec, files::AbstractVector{<:AbstractString}=p.spectra;
                         threads::Bool=true, curves::Bool=false, on_progress=nothing)
    run_batch(String.(files); threads=threads, on_progress=on_progress) do i
        analyze(p; spectrum=String(files[i]), curves=curves)
    end
end

analyze_spectra(p::ProjectSpec, spectra::AbstractString; filter::AbstractString="*", recursive::Bool=false, kw...) =
    analyze_spectra(p, spectrum_files(spectra; filter=filter, recursive=recursive); kw...)

analyze_spectra(xml_path::AbstractString, args...; kw...) = analyze_spectra(parse_config(xml_path), args...; kw...)
