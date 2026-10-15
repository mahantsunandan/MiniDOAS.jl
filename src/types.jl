# This file is part of MiniDOAS.jl (Differential Optical Absorption Spectroscopy in Julia).
# Copyright (c) 2026 Sunandan Mahant <sunandanmahant@outlook.com>
# Derived from QDOAS, Copyright (C) 1994-2025 BIRA-IASB and S[&]T (BSD-3-Clause).
# See LICENSE and NOTICE.md.

# =============================================================================
# What a fit is made of: cross sections, shift/stretch groups, analysis windows,
# and the settings shared by all windows (a QDOAS "project").
# =============================================================================

"""
    CrossSection(name, file; action=:interp)
    CrossSection(name, lambda, sigma; action=:interp)

An absorption cross section (or any other spectrum fitted with a linear coefficient, such
as a Ring spectrum): its `name`, the wavelengths `lambda` (nm, increasing; a decreasing
column is reversed) and the values `sigma` (e.g. cm² molecule⁻¹, so that the fitted
coefficient is a slant column in molecules cm⁻²).

From a file, the first two numeric columns are read (lines starting with `;`, `#` or `*`
are comments), as QDOAS reads them. `action` is `:interp` to interpolate onto the
analysis grid (the usual case) or `:none` when the file is already on the grid of the
spectra (one line per detector pixel, used as it is).

A cross section that a QDOAS project names (`parse_config`) keeps only its `file`; it is
read when a fit is set up, and read again when the file changes.
"""
struct CrossSection
    name::String
    file::String                  # where it comes from ("" when given as data)
    action::Symbol                # :interp, or :none (already on the calibration grid)
    lambda::Vector{Float64}       # empty: read from `file` when a fit is set up
    sigma::Vector{Float64}
end

function CrossSection(name::AbstractString, file::AbstractString; action::Symbol=:interp)
    action in (:interp, :none) || throw(ArgumentError("action must be :interp or :none"))
    lam, val, _ = load_xs(file)
    CrossSection(String(name), String(file), action, lam, val)
end

function CrossSection(name::AbstractString, lambda::AbstractVector{<:Real}, sigma::AbstractVector{<:Real};
                      action::Symbol=:interp)
    action in (:interp, :none) || throw(ArgumentError("action must be :interp or :none"))
    length(lambda) == length(sigma) || throw(DimensionMismatch("$(length(lambda)) wavelengths, $(length(sigma)) values"))
    length(lambda) >= 2 || throw(ArgumentError("a cross section needs at least two points"))
    lam = Vector{Float64}(lambda); val = Vector{Float64}(sigma)
    if lam[1] > lam[2]                            # as MATRIX_Load: a decreasing column is reversed
        reverse!(lam); reverse!(val)
    end
    CrossSection(String(name), "", action, lam, val)
end

"Whether the cross section's values are held in memory (else they are read from `file`)."
isloaded(c::CrossSection) = !isempty(c.lambda)

function Base.show(io::IO, c::CrossSection)
    print(io, "CrossSection(\"", c.name, "\"")
    if isloaded(c)
        @printf(io, ", %.2f-%.2f nm, %d points", first(c.lambda), last(c.lambda), length(c.lambda))
    end
    isempty(c.file) || print(io, ", ", basename(c.file))
    print(io, ")")
end

"""
    ShiftStretch(targets...; shift=true, stretch=0, shift_init=0.0, stretch_init=0.0,
                 stretch2_init=0.0, shift_step=1e-3, stretch_step=1e-3, stretch2_step=1e-3,
                 shift_limits=(0.0, 0.0), store=true)

A wavelength shift (and optionally a stretch) shared by `targets`: `"Spectrum"` (the
measured spectrum), `"Ref"` (the reference) and/or names of cross sections. The vectors are
evaluated at λ - Δ - s₁(λ-λ₀) - s₂(λ-λ₀)².

  * `shift`: fit the shift (else it is fixed at `shift_init`, nm);
  * `stretch`: 0 none, 1 first order, 2 second order (fitted, starting from the `*_init`
    values);
  * `*_step`: the finite-difference steps of the Levenberg-Marquardt fit;
  * `shift_limits`: (min, max) for the shift; (0, 0) means unlimited;
  * `store`: report the shift and stretch in QDOAS-layout outputs.
"""
struct ShiftStretch
    symbols::Vector{String}
    sh_fit::Bool
    st_type::Int          # 0 none, 1 first order, 2 second order
    sh_init::Float64
    st_init::Float64
    st_init2::Float64
    sh_delta::Float64
    st_delta::Float64
    st_delta2::Float64
    sh_min::Float64
    sh_max::Float64
    sh_store::Bool
    st_store::Bool
    err_store::Bool
end

ShiftStretch(targets::AbstractString...; kw...) = ShiftStretch(collect(String, targets); kw...)

function ShiftStretch(targets::AbstractVector{<:AbstractString}; shift::Bool=true, stretch::Union{Bool,Integer}=0,
                      shift_init::Real=0.0, stretch_init::Real=0.0, stretch2_init::Real=0.0,
                      shift_step::Real=1e-3, stretch_step::Real=1e-3, stretch2_step::Real=1e-3,
                      shift_limits=(0.0, 0.0), store::Bool=true)
    isempty(targets) && throw(ArgumentError("ShiftStretch needs at least one target"))
    st = Int(stretch)
    st in (0, 1, 2) || throw(ArgumentError("stretch must be 0, 1 or 2"))
    for d in (shift_step, stretch_step, stretch2_step)
        abs(d) >= EPSILON || throw(ArgumentError("finite-difference steps must be at least $EPSILON"))
    end
    ShiftStretch(collect(String, targets), shift, st, shift_init, stretch_init, stretch2_init,
                 shift_step, stretch_step, stretch2_step, Float64(shift_limits[1]), Float64(shift_limits[2]),
                 store, store, store)
end

function Base.show(io::IO, g::ShiftStretch)
    print(io, "ShiftStretch(", join(repr.(g.symbols), ", "), "; shift=", g.sh_fit)
    g.sh_fit || @printf(io, " (fixed at %g nm)", g.sh_init)
    g.st_type > 0 && print(io, ", stretch=", g.st_type)
    print(io, ")")
end

"""
    FitWindow(name, lmin, lmax, cross_sections=CrossSection[];
              polynomial=3, linear_offset=-1, linear_offset_ref=false, offset=-1,
              shifts=ShiftStretch[], gaps=[], lambda0=nothing, reference="")

An analysis window: the spectrum between `lmin` and `lmax` nm (minus `gaps`, a list of
(from, to) intervals) is fitted with

    ln(I₀/I) = Σᵢ Sᵢ σᵢ(λ) + P(λ - λ₀) [+ linear offset] + residual

  * `cross_sections`: the `CrossSection`s fitted (their coefficients Sᵢ are the slant
    columns);
  * `polynomial`: order of P (0-8), or -1 for none;
  * `linear_offset`: order (0-8) of a linear intensity-offset term, or -1 for none;
    normalised by the radiance, or by the reference when `linear_offset_ref`;
  * `offset`: order (0-2) of a non-linear intensity offset fitted inside the logarithm,
    or -1 for none (QDOAS does not allow both kinds of offset in one window);
  * `shifts`: `ShiftStretch` groups (e.g. `ShiftStretch("Spectrum"; stretch=1)`);
  * `lambda0`: the centre of the polynomial and stretch, by default the middle of the
    window;
  * `reference`: a reference file for this window; only used by `write_qdoas_project`
    and by QDOAS projects (a `DOASModel` is given its reference).

QDOAS project files describe the same windows (`parse_config(...).windows`).
"""
struct FitWindow
    name::String
    lmin::Float64
    lmax::Float64
    lambda0::Float64
    refone::String
    xs::Vector{CrossSection}
    poly_order::Int       # -1 = no polynomial
    offl_order::Int       # linear offset order, -1 = none
    offl_ref::Bool        # linear offset normalised by the reference (else by the radiance)
    # non-linear offset: (fit, init, delta) for orders 0,1,2
    off::NTuple{3,Tuple{Bool,Float64,Float64}}
    shifts::Vector{ShiftStretch}
    gaps::Vector{Tuple{Float64,Float64}}
    outputs::Dict{String,Tuple{Bool,Bool,Float64}}  # sym => (store SlCol, store SlErr, sfact)
    poly_store::Tuple{Bool,Bool}                # xfit, xerr
    offl_store::Tuple{Bool,Bool}                # offfit, offerr
    off_store::NTuple{3,Tuple{Bool,Bool}}       # o?fstr, o?estr
    save_residuals::Bool
end

function FitWindow(name::AbstractString, lmin::Real, lmax::Real, cross_sections::AbstractVector=CrossSection[];
                   polynomial::Integer=3, linear_offset::Integer=-1, linear_offset_ref::Bool=false,
                   offset::Integer=-1, shifts=ShiftStretch[], gaps=Tuple{Float64,Float64}[],
                   lambda0::Union{Nothing,Real}=nothing, reference::AbstractString="")
    lo, hi = Float64(min(lmin, lmax)), Float64(max(lmin, lmax))
    hi > lo || throw(ArgumentError("$name: empty fitting interval $lmin-$lmax nm"))
    -1 <= polynomial <= 8 || throw(ArgumentError("polynomial order must be -1 (none) to 8"))
    -1 <= linear_offset <= 8 || throw(ArgumentError("linear offset order must be -1 (none) to 8"))
    -1 <= offset <= 2 || throw(ArgumentError("non-linear offset order must be -1 (none) to 2"))
    offset >= 0 && linear_offset >= 0 &&
        throw(ArgumentError("$name: a window has a linear or a non-linear offset, not both (as in QDOAS)"))
    xs = collect(CrossSection, cross_sections)
    names = [c.name for c in xs]
    allunique(lowercase.(names)) || throw(ArgumentError("$name: cross-section names must be unique"))
    groups = shifts isa ShiftStretch ? [shifts] : collect(ShiftStretch, shifts)
    seen = Set{String}()
    for g in groups, s in g.symbols
        ls = lowercase(s)
        ls in ("spectrum", "ref") || ls in lowercase.(names) ||
            throw(ArgumentError("$name: shift target \"$s\" is not Spectrum, Ref or a cross section"))
        ls in seen && throw(ArgumentError("$name: \"$s\" appears in two shift groups"))
        push!(seen, ls)
    end
    off = ntuple(k -> (k - 1 <= offset, 0.0, 1.0e-3), 3)
    l0 = lambda0 === nothing ? 0.5 * (lo + hi) : Float64(lambda0)
    FitWindow(String(name), lo, hi, l0, String(reference), xs, polynomial, linear_offset, linear_offset_ref,
              off, groups, [(Float64(min(a, b)), Float64(max(a, b))) for (a, b) in gaps],
              Dict(c.name => (true, true, 1.0) for c in xs), (true, true), (true, true),
              ntuple(k -> (off[k][1], off[k][1]), 3), true)
end

function Base.show(io::IO, w::FitWindow)
    @printf(io, "FitWindow(\"%s\", %.2f-%.2f nm, %s, polynomial %d", w.name, w.lmin, w.lmax,
            isempty(w.xs) ? "no cross sections" : join((c.name for c in w.xs), ", "), w.poly_order)
    w.offl_order >= 0 && print(io, ", linear offset ", w.offl_order)
    any(o -> o[1], w.off) && print(io, ", offset ", findlast(o -> o[1], w.off) - 1)
    isempty(w.shifts) || print(io, ", shifts of ", join((join(g.symbols, "+") for g in w.shifts), ", "))
    print(io, ")")
end

"""
Settings shared by all analysis windows, as in a QDOAS project: interpolation, the
fit's convergence criteria, the instrument (MFC STD with `n_det` pixels, dark and offset
files), the spectra, the windows and the outputs. `parse_config` reads one from a QDOAS
project file; a `DOASModel` made from Julia carries one too.
"""
struct ProjectSpec
    interpolation::Symbol     # :spline or :linear
    security_gap::Int
    convergence::Float64
    max_iterations::Int
    spike_tolerance::Float64
    n_det::Int
    revert::Bool
    calib::String
    dark::String
    offset::String
    spectrum::String                            # the spectrum `analyze` uses by default
    spectra::Vector{String}                     # every spectrum the project lists
    windows::Vector{FitWindow}
    fields::Vector{String}                      # project <output> fields, in order
    date_format::String                         # mfcstd date layout
    project_name::String
    swath_name::String                          # NetCDF results group (output swathName)
end

function Base.show(io::IO, p::ProjectSpec)
    print(io, "ProjectSpec(\"", p.project_name, "\", ", p.n_det, " pixels, windows ",
          join((w.name for w in p.windows), ", "), ", ", length(p.spectra), " spectra)")
end

"Project output fields reproduced by write_netcdf; any other field is refused."
const SUPPORTED_FIELDS = ("name", "date", "time", "rms", "chi", "iter_number", "error_flag", "residual_spectrum")

"A copy of `p` with some fields replaced, e.g. `respec(p; spectrum=path)`."
respec(p::ProjectSpec; kw...) =
    ProjectSpec((haskey(kw, f) ? kw[f] : getfield(p, f) for f in fieldnames(ProjectSpec))...)

"A copy of analysis window `w` with some fields replaced, e.g. `rewindow(w; lmin=300.0)`."
rewindow(w::FitWindow; kw...) =
    FitWindow((haskey(kw, f) ? kw[f] : getfield(w, f) for f in fieldnames(FitWindow))...)

"""
`w` with the fitting interval `lmin`-`lmax`. λ₀ moves with the window when it was at the
window's centre (the default), and stays where it was set otherwise.
"""
function with_range(w::FitWindow, lmin::Real, lmax::Real)
    mid_before = abs(w.lambda0 - 0.5 * (w.lmin + w.lmax)) < 1e-9
    rewindow(w; lmin=Float64(lmin), lmax=Float64(lmax), lambda0=mid_before ? 0.5 * (lmin + lmax) : w.lambda0)
end
