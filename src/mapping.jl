# This file is part of MiniDOAS.jl (Differential Optical Absorption Spectroscopy in Julia).
# Copyright (c) 2026 Sunandan Mahant <sunandanmahant@outlook.com>; BSD-3-Clause, see LICENSE.

# =============================================================================
# Retrieval interval mapping (Vogel et al., AMT 6, 275, 2013): one analysis window fitted
# to one spectrum for every combination of start and end wavelength on a grid.
# =============================================================================

"""
    WindowMap

The result of `fit_window_map`: the window's name and its own fitting interval `range`, the
`starts` and `ends` of the mapped intervals (nm), and per cell (start i, end j) the fit
`results[i, j]` (a `WindowResult`, or `nothing` where the interval was too narrow or the
fit failed; `errors[i, j]` says why).

`map[quantity]` gives a matrix of one quantity (NaN where there is no fit): the name of a
fitted term for its slant column (e.g. `"NO2"`), `"err:NO2"` for its error, `"rel:NO2"`
for the relative error, or `"rms"`, `"chi2"`, `"iterations"`.
"""
struct WindowMap
    window::String
    range::Tuple{Float64,Float64}
    starts::Vector{Float64}
    ends::Vector{Float64}
    results::Matrix{Union{Nothing,WindowResult}}
    errors::Matrix{String}
end

function Base.getindex(m::WindowMap, quantity::AbstractString)
    q = String(quantity)
    get_q(wr) = q == "rms" ? wr.rms :
                q == "chi2" ? wr.chisqr :
                q == "iterations" ? Float64(wr.niter) :
                startswith(q, "err:") ? get(wr.slerr, q[5:end], NaN) :
                startswith(q, "rel:") ? get(wr.slerr, q[5:end], NaN) / abs(get(wr.slcol, q[5:end], NaN)) :
                get(wr.slcol, q, NaN)
    [r === nothing ? NaN : get_q(r) for r in m.results]
end

"A quantity of `m[quantity]` in words, for labels."
describe_quantity(q::AbstractString) =
    q == "rms" ? "RMS" : q == "chi2" ? "χ²" : q == "iterations" ? "iterations" :
    startswith(q, "err:") ? "error of $(q[5:end])" : startswith(q, "rel:") ? "relative error of $(q[5:end])" :
    "slant column of $q"

function Base.show(io::IO, m::WindowMap)
    nok = count(!isnothing, m.results)
    @printf(io, "WindowMap(\"%s\", starts %.2f-%.2f nm, ends %.2f-%.2f nm, %d fits)", m.window,
            first(m.starts), last(m.starts), first(m.ends), last(m.ends), nok)
end

"""
    fit_window_map(model::DOASModel, spectrum; window=1, starts, ends, min_width=5.0, threads=true)
    fit_window_map(p::ProjectSpec, spectrum_file; window=1, starts, ends, min_width=5.0, threads=true)

Fits one analysis window (its name or number) to one spectrum for every fitting interval
from `starts[i]` to `ends[j]` (nm) at least `min_width` nm wide, with all other settings of
the window unchanged (λ₀ follows the window's centre unless it was set). A good fit window
is one where the slant columns change little when the window moves a little (*retrieval
interval mapping*, Vogel et al., 2013). Returns a `WindowMap`; cells are fitted in
parallel when Julia has several threads.

`on_progress(done, total)` is called after each cell (one call at a time); when it returns
`false`, the remaining cells are skipped (their error is "cancelled").
"""
function fit_window_map(m::DOASModel, spectrum; window=1, starts::AbstractVector{<:Real},
                        ends::AbstractVector{<:Real}, min_width::Real=5.0, threads::Bool=true,
                        on_progress=nothing)
    k = window isa Integer ? window : findfirst(w -> w.spec.name == window, m.windows)
    (k === nothing || !(1 <= k <= length(m.windows))) && throw(ArgumentError("no analysis window $window"))
    raw = spectrum isa MFCSpectrum ? spectrum.counts : spectrum
    p = m.spec
    length(raw) == p.n_det || throw(DimensionMismatch("spectrum of $(length(raw)) pixels, the model has $(p.n_det)"))
    starts = collect(Float64, starts); ends = collect(Float64, ends)
    length(starts) * length(ends) <= 100_000 || throw(ArgumentError("too many cells; use a coarser grid"))
    spec = Vector{Float64}(raw)
    spe_norm = normalize!(spec)
    w0 = m.windows[k]
    ref, ref_norm = m.refs[k], m.ref_norms[k]
    results = Matrix{Union{Nothing,WindowResult}}(nothing, length(starts), length(ends))
    errors = fill("", length(starts), length(ends))
    cells = [(i, j) for i in eachindex(starts), j in eachindex(ends)]
    done = Threads.Atomic{Int}(0)
    go = Threads.Atomic{Bool}(true)
    lk = ReentrantLock()
    function cell(c)
        i, j = c
        if ends[j] - starts[i] < min_width
            errors[i, j] = "narrower than $min_width nm"
        elseif !go[]
            errors[i, j] = "cancelled"
        else
            try
                ws = with_range(w0.spec, starts[i], ends[j])
                w = build_window(ws, p, w0.lam, m.xs_cache)
                results[i, j] = window_result(w, p, spec, ref; spe_norm=spe_norm, ref_norm=ref_norm)
            catch e
                e isa InterruptException && rethrow()
                errors[i, j] = sprint(showerror, e)
            end
        end
        k = Threads.atomic_add!(done, 1) + 1
        on_progress === nothing || lock(lk) do
            on_progress(k, length(cells)) === false && (go[] = false)
        end
        nothing
    end
    parallel = threads && Threads.nthreads() > 1
    nblas = BLAS.get_num_threads()
    parallel && BLAS.set_num_threads(1)
    try
        if parallel
            Threads.@threads :dynamic for c in cells
                cell(c)
            end
        else
            foreach(cell, cells)
        end
    finally
        parallel && BLAS.set_num_threads(nblas)
    end
    WindowMap(w0.spec.name, (w0.spec.lmin, w0.spec.lmax), starts, ends, results, errors)
end

function fit_window_map(p::ProjectSpec, spectrum::AbstractString; window=1, kw...)
    k = window isa Integer ? window : findfirst(w -> w.name == window, p.windows)
    (k === nothing || !(1 <= k <= length(p.windows))) && throw(ArgumentError("no analysis window $window"))
    q = respec(p; spectrum=String(spectrum), windows=[p.windows[k]])
    raw, _ = load_spectrum(q)
    fit_window_map(DOASModel(q), raw; window=1, kw...)
end
