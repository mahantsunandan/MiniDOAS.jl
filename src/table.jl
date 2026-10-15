# This file is part of MiniDOAS.jl (Differential Optical Absorption Spectroscopy in Julia).
# Copyright (c) 2026 Sunandan Mahant <sunandanmahant@outlook.com>
# Derived from QDOAS, Copyright (C) 1994-2025 BIRA-IASB and S[&]T (BSD-3-Clause).
# See LICENSE and NOTICE.md.

# =============================================================================
# Tables and summaries
# =============================================================================

"The quantities of one window result as (column name => value) pairs, in a fixed order."
function window_quantities(wr::WindowResult)
    q = Pair{String,Float64}[]
    push!(q, "RMS" => wr.rms, "Chi2" => wr.chisqr, "iterations" => wr.niter)
    for c in wr.columns
        push!(q, "SlCol($c)" => wr.slcol[c], "SlErr($c)" => wr.slerr[c])
    end
    for s in sort(collect(keys(wr.shift)))
        push!(q, "Shift($s)" => wr.shift[s], "Err Shift($s)" => wr.shift_err[s])
    end
    for s in sort(collect(keys(wr.stretch)))
        push!(q, "Stretch($s)" => wr.stretch[s], "Err Stretch($s)" => wr.stretch_err[s])
    end
    for k in sort(collect(keys(wr.params)))
        push!(q, k => wr.params[k], "Err $k" => wr.param_err[k])
    end
    q
end

fmt_date(d) = @sprintf("%04d-%02d-%02d", d...)
fmt_time(t) = @sprintf("%02d:%02d:%02d", t...)

"A single result as a one-item batch."
as_batch(r::DOASResult) = [(file=r.spectrum, result=r, error="", seconds=0.0)]
as_batch(batch::AbstractVector) = batch

const RECORD_META = (:elevation, :azimuth, :latitude, :longitude)

record_value(r, f::Symbol) = (r === nothing || r.record === nothing) ? NaN :
    f === :sza ? sza_of(r.record) : getfield(r.record, f)

"Solar zenith angle at the record's time and position (NaN without a position)."
function sza_of(rec::MFCSpectrum)
    (isnan(rec.latitude) || isnan(rec.longitude) || rec.date[1] == 0) && return NaN
    solar_position(measurement_time(rec), rec.latitude, rec.longitude).zenith
end

"""
    results_table(batch; sza=false) -> (header, rows)

One row per spectrum of a batch (`analyze_spectra`, or `doas_fit` with several spectra)
or of a single `DOASResult`: file, record name, date, mid-exposure time (UTC as in the
file), status, error, the viewing `elevation` and `azimuth` angles and the `latitude` and
`longitude` when any spectrum's header gives them (with `sza=true` also the solar zenith
angle computed from time and position), then for every window `Window.RMS`,
`Window.Chi2`, `Window.iterations`, `Window.SlCol(X)`, `Window.SlErr(X)` for each fitted
term (cross sections, then polynomial `x0..`), shifts, stretches and non-linear offsets
with their errors. SCDs are in the units of the cross sections times length (molecules
cm⁻² for cross sections in cm² molecule⁻¹); shifts in nm.
"""
function results_table(batch::AbstractVector; sza::Bool=false)
    cols = String[]; seen = Set{String}()
    for b in batch
        b.result === nothing && continue
        for wr in b.result.windows, (k, _) in window_quantities(wr)
            name = "$(wr.name).$k"
            name in seen || (push!(cols, name); push!(seen, name))
        end
    end
    meta = Symbol[f for f in RECORD_META if any(b -> !isnan(record_value(b.result, f)), batch)]
    sza && push!(meta, :sza)
    header = vcat(["file", "name", "date", "time", "status", "error"], String.(meta), cols)
    rows = Vector{Vector{Any}}()
    for b in batch
        r = b.result
        row = Any[b.file, "", "", "", r === nothing ? "failed" : "ok", b.error]
        vals = Dict{String,Float64}()
        if r !== nothing
            if r.record !== nothing
                row[2] = r.record.name
                row[3] = fmt_date(r.record.date)
                row[4] = fmt_time(r.record.time)
            end
            for wr in r.windows, (k, v) in window_quantities(wr)
                vals["$(wr.name).$k"] = v
            end
        end
        append!(row, [(v = record_value(r, f); isnan(v) ? missing : v) for f in meta])
        append!(row, [get(vals, c, missing) for c in cols])
        push!(rows, row)
    end
    header, rows
end

results_table(r::DOASResult; kw...) = results_table(as_batch(r); kw...)

"""
    results_columns(batch; sza=false) -> NamedTuple of column vectors

The `results_table` as columns, keyed by column name: a Tables.jl column table, so
`DataFrame(results_columns(batch))` (DataFrames.jl) or `CSV.write(path, ...)` work as they
are. Numbers are `Float64` (or `Missing` where a spectrum has no value); `date` and `time`
are strings.
"""
function results_columns(batch; sza::Bool=false)
    header, rows = results_table(as_batch(batch); sza=sza)
    cols = map(eachindex(header)) do j
        v = [r[j] for r in rows]
        if j <= 6
            j == 5 || j == 6 ? String.(v) : [x === "" ? missing : String(x) for x in v]
        elseif all(x -> x isa Real, v)
            Float64.(v)
        else
            Union{Missing,Float64}[x === missing ? missing : Float64(x) for x in v]
        end
    end
    NamedTuple{Tuple(Symbol.(header))}(Tuple(cols))
end

csv_cell(x::AbstractFloat) = isfinite(x) ? @sprintf("%.10g", x) : ""
csv_cell(::Missing) = ""
csv_cell(x) = (s = string(x); occursin(r"[,\"\n]", s) ? "\"" * replace(s, "\"" => "\"\"") * "\"" : s)

"""
    write_csv(path, batch_or_result; sza=false) -> path

Writes `results_table` as comma-separated values.
"""
function write_csv(path::AbstractString, batch; sza::Bool=false)
    header, rows = results_table(as_batch(batch); sza=sza)
    open(path, "w") do io
        println(io, join(csv_cell.(header), ','))
        for r in rows
            println(io, join(csv_cell.(r), ','))
        end
    end
    path
end

"""
    print_summary([io,] result)

Human-readable summary of a `DOASResult`: per window the RMS, χ², iterations,
every slant column with its error, and the shifts and stretches.
"""
function print_summary(io::IO, res::DOASResult)
    r = res.record
    if r !== nothing
        @printf(io, "%s  %s %s UTC  (%s)\n", r.name, fmt_date(r.date), fmt_time(r.time), basename(res.spectrum))
    elseif !isempty(res.spectrum)
        println(io, res.spectrum)
    end
    for wr in res.windows
        show(io, MIME("text/plain"), wr)
        println(io)
    end
end
print_summary(res::DOASResult) = print_summary(stdout, res)

function Base.show(io::IO, ::MIME"text/plain", wr::WindowResult)
    @printf(io, "  %s: RMS %.3e, chi2 %.3e, %d iterations\n", wr.name, wr.rms, wr.chisqr, wr.niter)
    for c in wr.columns
        @printf(io, "    %-22s % .4e ± %.2e\n", c, wr.slcol[c], wr.slerr[c])
    end
    # QDOAS gives a parameter that is not fitted an error of exactly 1
    for s in sort(collect(keys(wr.shift)))
        if wr.shift_err[s] == 1.0
            @printf(io, "    shift(%s) % .5f nm (fixed)", s, wr.shift[s])
        else
            @printf(io, "    shift(%s) % .5f ± %.5f nm", s, wr.shift[s], wr.shift_err[s])
        end
        get(wr.stretch_err, s, 1.0) != 1.0 && @printf(io, ", stretch % .3e ± %.1e", wr.stretch[s], wr.stretch_err[s])
        println(io)
    end
    for k in sort(collect(keys(wr.params)))
        @printf(io, "    %s % .4e ± %.2e\n", k, wr.params[k], wr.param_err[k])
    end
end
Base.show(io::IO, wr::WindowResult) = print(io, "WindowResult(", wr.name, ", RMS=", round(wr.rms; sigdigits=4), ")")
