# This file is part of MiniDOAS.jl (Differential Optical Absorption Spectroscopy in Julia).
# Copyright (c) 2026 Sunandan Mahant <sunandanmahant@outlook.com>; BSD-3-Clause, see LICENSE.
# The Ring algorithm (`ring_spectrum`) is QDOAS's (engine/ring.c, mediator/mediate_xsconv.c,
# Copyright (C) 1994-2025 BIRA-IASB and S[&]T, BSD-3-Clause, see NOTICE.md).

# =============================================================================
# Preparing what a fit needs: slit functions and convolution, I0-corrected cross
# sections, Ring spectra, air and vacuum wavelengths, a wavelength calibration against
# the solar spectrum; and simple spectral filters and simulated spectra.
# =============================================================================

"""
    SlitFunction

An instrument's slit function (line shape). Implementations: `GaussianSlit`,
`BoxcarSlit`, `TriangularSlit`, `LorentzianSlit`, `TabulatedSlit`. Widths are in nm and
may be numbers or functions of the wavelength (e.g. `GaussianSlit(λ -> 0.35 + 0.001(λ - 330))`).
"""
abstract type SlitFunction end

"""
    GaussianSlit(fwhm)

A Gaussian slit function of full width at half maximum `fwhm` (nm; a number or a function
of the wavelength), cut at ±5 standard deviations.
"""
struct GaussianSlit{F} <: SlitFunction
    fwhm::F
end

"""
    BoxcarSlit(width)

A rectangular slit function of full width `width` (nm).
"""
struct BoxcarSlit{F} <: SlitFunction
    width::F
end

"""
    TriangularSlit(fwhm)

A triangular slit function of full width at half maximum `fwhm` (nm; base 2·fwhm).
"""
struct TriangularSlit{F} <: SlitFunction
    fwhm::F
end

"""
    LorentzianSlit(fwhm)

A Lorentzian slit function of full width at half maximum `fwhm` (nm), cut at ±50 FWHM.
"""
struct LorentzianSlit{F} <: SlitFunction
    fwhm::F
end

"""
    TabulatedSlit(offset, value)

A measured slit function: the response of a pixel centred at λ to light at λ + `offset`
(nm), interpolated linearly and zero outside the table.
"""
struct TabulatedSlit <: SlitFunction
    offset::Vector{Float64}
    value::Vector{Float64}
    function TabulatedSlit(offset::AbstractVector{<:Real}, value::AbstractVector{<:Real})
        length(offset) == length(value) >= 2 || throw(ArgumentError("offset and value of equal length >= 2"))
        p = sortperm(offset)
        new(Float64.(offset[p]), Float64.(value[p]))
    end
end

widthat(f::Real, λ) = Float64(f)
widthat(f, λ) = Float64(f(λ))

const FWHM_TO_SIGMA = 1 / (2 * sqrt(2 * log(2)))

"Half-width (nm) beyond which the slit function `s` at wavelength λ is zero."
support(s::GaussianSlit, λ) = 5 * (widthat(s.fwhm, λ) / (2 * sqrt(2 * log(2))))
support(s::BoxcarSlit, λ) = widthat(s.width, λ) / 2
support(s::TriangularSlit, λ) = widthat(s.fwhm, λ)
support(s::LorentzianSlit, λ) = 50 * widthat(s.fwhm, λ)
support(s::TabulatedSlit, λ) = max(abs(first(s.offset)), abs(last(s.offset)))

"A weight function Δ -> slit(Δ) for the slit function at wavelength λ (not normalised)."
function kernel(s::GaussianSlit, λ)
    σ = widthat(s.fwhm, λ) / (2 * sqrt(2 * log(2)))
    Δ -> exp(-0.5 * (Δ / σ)^2)
end
kernel(s::BoxcarSlit, λ) = (h = widthat(s.width, λ) / 2; Δ -> abs(Δ) <= h ? 1.0 : 0.0)
kernel(s::TriangularSlit, λ) = (f = widthat(s.fwhm, λ); Δ -> max(0.0, 1 - abs(Δ) / f))
kernel(s::LorentzianSlit, λ) = (γ = widthat(s.fwhm, λ) / 2; Δ -> 1 / (1 + (Δ / γ)^2))
function kernel(s::TabulatedSlit, λ)
    x, v = s.offset, s.value
    Δ -> begin
        (Δ < x[1] || Δ > x[end]) && return 0.0
        k = clamp(searchsortedlast(x, Δ), 1, length(x) - 1)
        t = (Δ - x[k]) / (x[k+1] - x[k])
        (1 - t) * v[k] + t * v[k+1]
    end
end

as_slit(s::SlitFunction) = s
as_slit(f) = GaussianSlit(f)

"""
    convolve(x, y, at, slit) -> Vector
    convolve(xs::CrossSection, at, slit) -> CrossSection

The high-resolution spectrum `(x, y)` (x increasing, nm) convolved with the slit function
`slit` (a `SlitFunction`, or a FWHM in nm for a Gaussian) and evaluated at the wavelengths
`at`: Σⱼ s(xⱼ-λ) yⱼ wⱼ / Σⱼ s(xⱼ-λ) wⱼ with trapezoidal weights wⱼ. The input must sample
the slit function finely; where it does not cover the slit at all the result is NaN.
"""
function convolve(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, at::AbstractVector{<:Real}, slit)
    length(x) == length(y) || throw(DimensionMismatch("x and y of different lengths"))
    issorted(x) || throw(ArgumentError("the wavelengths must increase"))
    convolve_sorted(x, y, at, as_slit(slit))
end

"`convolve` without checks (x increasing, `s` a SlitFunction)."
function convolve_sorted(x, y, at, s::SlitFunction)
    out = Vector{Float64}(undef, length(at))
    for (i, λ) in enumerate(at)
        h = support(s, λ)
        k = kernel(s, λ)
        lo = searchsortedfirst(x, λ - h); hi = searchsortedlast(x, λ + h)
        num = 0.0; den = 0.0
        for j in max(lo, 2):min(hi, length(x) - 1)
            g = k(x[j] - λ) * (x[j+1] - x[j-1])
            num += g * y[j]; den += g
        end
        out[i] = den > 0 ? num / den : NaN
    end
    out
end

function convolve(xs::CrossSection, at::AbstractVector{<:Real}, slit)
    isloaded(xs) || (xs = CrossSection(xs.name, xs.file; action=xs.action))
    CrossSection(xs.name, at, convolve(xs.lambda, xs.sigma, at, slit))
end

"Values of `xs` at the wavelengths `x` (cubic spline; the end values outside its range)."
function xs_at(xs::CrossSection, x::AbstractVector{<:Real})
    isloaded(xs) || (xs = CrossSection(xs.name, xs.file; action=xs.action))
    d2 = spline_deriv2(xs.lambda, xs.sigma)
    [spline_at(xs.lambda, xs.sigma, d2, Float64(v), true) for v in x]
end

"""
    convolve_i0(xs::CrossSection, solar_lambda, solar, at, slit; scd=1e17) -> CrossSection

A cross section convolved with the I₀ correction (Aliwell et al., 2002), for absorbers whose
strong, narrow bands are not resolved by the instrument:

    σ_eff(λ) = ln( (I₀ ⊗ s)(λ) / ((I₀ e^{-σ·SCD}) ⊗ s)(λ) ) / SCD

with I₀ the high-resolution solar spectrum (`solar_lambda`, `solar`), `s` the slit and `SCD`
a typical slant column of the absorber (molecules cm⁻² for σ in cm² molecule⁻¹). For small
`scd` this tends to the I₀-weighted convolution ((I₀σ) ⊗ s)/(I₀ ⊗ s); saturation of the
strong bands makes it smaller for larger `scd`. The calculation uses the solar spectrum's
wavelength grid within the cross section's range.
"""
function convolve_i0(xs::CrossSection, solar_lambda::AbstractVector{<:Real}, solar::AbstractVector{<:Real},
                     at::AbstractVector{<:Real}, slit; scd::Real=1e17)
    isloaded(xs) || (xs = CrossSection(xs.name, xs.file; action=xs.action))
    keep = findall(l -> xs.lambda[1] <= l <= xs.lambda[end], solar_lambda)
    length(keep) >= 3 || throw(ArgumentError("the solar spectrum does not cover the cross section"))
    λh = Float64.(solar_lambda[keep]); Ih = Float64.(solar[keep])
    σh = xs_at(xs, λh)
    c0 = convolve(λh, Ih, at, slit)
    c1 = convolve(λh, Ih .* exp.(-σh .* scd), at, slit)
    CrossSection(xs.name, at, log.(c0 ./ c1) ./ scd)
end

# -----------------------------------------------------------------------------
# Ring spectrum: QDOAS's Ring tool (mediate_xsconv.c mediateRingCalculate and
# ring.c raman_convolution, after Chance and Spurr 1997)
# -----------------------------------------------------------------------------

"Rotational Raman line strengths of N2 and O2 at temperature `temp` (K)."
function raman_lines(temp)
    c2 = 1.438769; emult = -c2 / temp
    qn2 = sum((N2STAT_1[i] * 2 + 1) * N2STAT_2[i] * exp(emult * N2STAT_3[i]) for i in eachindex(N2STAT_1))
    n2 = [256 * π^5 * 0.79 / 27 * (N2DEG[i] * 2 + 1) * N2NUC[i] * exp(emult * N2TERM[i]) / qn2 * N2PLACTEL[i]
          for i in eachindex(N2TERM)]
    qo2 = sum((O2STAT_1[i] * 2 + 1) * exp(emult * O2STAT_2[i]) for i in eachindex(O2STAT_1))
    o2 = [256 * π^5 * 0.21 / 27 * (O2DEG[i] * 2 + 1) * exp(emult * O2TERM[i]) / qo2 * O2PLACTEL[i]
          for i in eachindex(O2TERM)]
    n2, o2
end

"""
    ring_spectrum(solar_lambda, solar, at, slit; temperature=250.0) -> Vector

A Ring spectrum on the wavelengths `at`, computed with QDOAS's Ring algorithm (Chance and
Spurr, 1997): the high-resolution solar spectrum is convolved with the slit, then with
the rotational Raman lines of N₂ and O₂ at `temperature` (K), and divided by the convolved
solar spectrum. Fit it like a cross section (`CrossSection("Ring", at, ring)`).
"""
function ring_spectrum(solar_lambda::AbstractVector{<:Real}, solar::AbstractVector{<:Real},
                       at::AbstractVector{<:Real}, slit; temperature::Real=250.0)
    sλ = Float64.(solar_lambda)
    keep = (sλ .> at[1] - 5) .& (sλ .< at[end] + 5)
    λ = sλ[keep]
    length(λ) >= 4 || throw(ArgumentError("the solar spectrum does not cover $(at[1])-$(at[end]) nm"))
    sol = convolve(sλ, solar, λ, slit)              # solar spectrum at the instrument's resolution
    d2 = spline_deriv2(λ, sol)
    n2, o2 = raman_lines(temperature)
    raman = similar(λ)
    for (i, l) in enumerate(λ)
        ν = 1e7 / l
        σsq = 1e6 / l^2
        γn2 = (-0.601466 + 238.557 / (186.099 - σsq))^2
        γo2 = (0.07149 + 45.9364 / (48.2716 - σsq))^2
        acc = 0.0; tot = 0.0
        for (w, pos, γ) in ((n2, N2POS, γn2), (o2, O2POS, γo2))
            for j in eachindex(w)
                νp = ν + pos[j]
                xsec = w[j] * νp^4 * γ
                tot += xsec
                acc += spline_at(λ, sol, d2, 1e7 / νp, true) * xsec
            end
        end
        raman[i] = acc / tot
    end
    r2 = spline_deriv2(λ, raman)
    [(r = spline_at(λ, raman, r2, Float64(a), true); s = spline_at(λ, sol, d2, Float64(a), true);
      r > 0 && s > 0 ? r / s : 0.0) for a in at]
end

# -----------------------------------------------------------------------------
# air and vacuum wavelengths
# -----------------------------------------------------------------------------

"Refractive index of standard air (Edlén 1966, revised by Birch and Downs 1994); λ in nm (vacuum)."
air_index(λ) = (σ2 = (1e3 / λ)^2; 1 + 1e-8 * (8342.54 + 2406147 / (130 - σ2) + 15998 / (38.9 - σ2)))

"""
    air_to_vacuum(λ) / vacuum_to_air(λ)

Converts a wavelength (nm) between standard air and vacuum with the refractive index of
Edlén (1966) as revised by Birch and Downs (1994). Broadcast over vectors:
`air_to_vacuum.(λ)`. Solar spectra and most recent cross sections are given in vacuum;
older data sets often in air.
"""
air_to_vacuum(λ::Real) = (v = Float64(λ); for _ in 1:5; v = λ * air_index(v); end; v)
vacuum_to_air(λ::Real) = λ / air_index(λ)
@doc (@doc air_to_vacuum) vacuum_to_air

# -----------------------------------------------------------------------------
# polynomials, smoothing, differential spectra
# -----------------------------------------------------------------------------

"Least-squares polynomial of `order` in (x - mean(x)), returned as a function."
function polyfit(x, y, order)
    x0 = sum(x) / length(x)
    c = hcat([(x .- x0) .^ k for k in 0:order]...) \ y
    v -> sum(c[k+1] * (v - x0)^k for k in 0:order)
end

"""
    differential(lambda, y; order=3) -> Vector

`y` minus its least-squares polynomial of `order` in wavelength: the differential
("narrow-band") part of a cross section or optical depth, as the DOAS fit sees it.
"""
function differential(lambda::AbstractVector{<:Real}, y::AbstractVector{<:Real}; order::Integer=3)
    length(lambda) == length(y) || throw(DimensionMismatch("lambda and y of different lengths"))
    order >= 0 || throw(ArgumentError("order must be >= 0"))
    x = Float64.(lambda)
    x0 = sum(x) / length(x); h = max(maximum(abs, x .- x0), eps())
    A = hcat([((x .- x0) ./ h) .^ k for k in 0:order]...)
    y .- A * (A \ Float64.(y))
end

"""
    smooth(y, n; kernel=:boxcar) -> Vector

`y` smoothed over `n` points (odd) with a `:boxcar` (running mean), `:triangular`,
`:binomial` or `:gaussian` (FWHM of `n` points) kernel. At the ends the kernel is cut and
renormalised.
"""
function smooth(y::AbstractVector{<:Real}, n::Integer; kernel::Symbol=:boxcar)
    n >= 1 || throw(ArgumentError("n must be >= 1"))
    w = if kernel === :boxcar
        isodd(n) || throw(ArgumentError("n must be odd"))
        ones(n)
    elseif kernel === :triangular
        isodd(n) || throw(ArgumentError("n must be odd"))
        h = (n + 1) ÷ 2; [Float64(min(k, n + 1 - k)) for k in 1:n] ./ h
    elseif kernel === :binomial
        isodd(n) || throw(ArgumentError("n must be odd"))
        [Float64(binomial(big(n - 1), k)) for k in 0:n-1]
    elseif kernel === :gaussian
        σ = n * FWHM_TO_SIGMA; h = ceil(Int, 3n)
        [exp(-0.5 * (k / σ)^2) for k in -h:h]
    else
        throw(ArgumentError("unknown kernel $kernel"))
    end
    h = length(w) ÷ 2
    m = length(y)
    out = Vector{Float64}(undef, m)
    for i in 1:m
        num = 0.0; den = 0.0
        for k in max(1, i - h):min(m, i + h)
            g = w[k-i+h+1]
            num += g * y[k]; den += g
        end
        out[i] = num / den
    end
    out
end

# -----------------------------------------------------------------------------
# wavelength calibration against the solar spectrum
# -----------------------------------------------------------------------------

"""
    calibrate_wavelength(lambda, intensity, solar_lambda, solar; windows, absorbers=[],
                         polynomial=3, shifts=-0.4:0.02:0.4, fwhms=0.30:0.05:1.20,
                         shift_order=2, fwhm_order=1)

Calibrates a spectrometer against a high-resolution solar spectrum, in the spirit of
QDOAS's "Kurucz" calibration: in each sub-window `(a, b)` of `windows` (nm) it finds the
wavelength shift δ and Gaussian slit FWHM that best fit

    ln I(λ) = ln (solar ⊗ slit)(λ + δ) + polynomial + Σ cᵢ absorberᵢ(λ)

(`absorbers`: vectors on `lambda`, e.g. an O₃ cross section and a Ring spectrum at the
instrument's resolution), by a grid search over `shifts` and `fwhms` refined three
times. Smooth polynomials in λ are then fitted through the shifts and widths.

Returns `(centres, shift, fwhm, shift_fit, fwhm_fit, lambda)`: per sub-window the centre,
shift and FWHM; the smooth shift and FWHM as functions of λ; and the calibrated grid
`lambda .+ shift_fit.(lambda)`. Use `GaussianSlit(result.fwhm_fit)` to convolve cross
sections onto that grid.
"""
function calibrate_wavelength(lambda::AbstractVector{<:Real}, intensity::AbstractVector{<:Real},
                              solar_lambda::AbstractVector{<:Real}, solar::AbstractVector{<:Real};
                              windows, absorbers=Vector{Float64}[], polynomial::Integer=3,
                              shifts=-0.40:0.02:0.40, fwhms=0.30:0.05:1.20,
                              shift_order::Integer=2, fwhm_order::Integer=1)
    λ = Float64.(lambda); counts = Float64.(intensity)
    length(λ) == length(counts) || throw(DimensionMismatch("lambda and intensity of different lengths"))
    all(a -> length(a) == length(λ), absorbers) || throw(DimensionMismatch("absorbers must be on lambda"))
    sλ = Float64.(solar_lambda); sy = Float64.(solar)
    issorted(sλ) || throw(ArgumentError("the solar wavelengths must increase"))
    centres = Float64[]; dshift = Float64[]; dfwhm = Float64[]
    for (a, b) in windows
        idx = findall(x -> a <= x <= b, λ)
        length(idx) > polynomial + 2 + length(absorbers) || throw(ArgumentError("too few pixels in $a-$b nm"))
        any(<=(0), counts[idx]) && throw(ArgumentError("non-positive intensity in $a-$b nm"))
        x = λ[idx]; lny = log.(counts[idx])
        t = (x .- (a + b) / 2) ./ ((b - a) / 2)
        extra = [Float64.(v[idx]) for v in absorbers]
        A = hcat([t .^ k for k in 0:polynomial]..., extra...)
        function misfit(δ, w)
            m = log.(convolve_sorted(sλ, sy, x .+ δ, GaussianSlit(w)))
            any(!isfinite, m) && return Inf
            r = (lny .- m) - A * (A \ (lny .- m))
            sum(abs2, r)
        end
        best = (Inf, 0.0, 0.5)
        for w in fwhms, δ in shifts
            v = misfit(δ, w); v < best[1] && (best = (v, δ, w))
        end
        gridstep(r, d) = r isa AbstractRange ? Float64(step(r)) : length(r) > 1 ? Float64(r[2] - r[1]) : d
        dδ0 = gridstep(shifts, 0.02); dw0 = gridstep(fwhms, 0.05)
        for (dδ, dw) in ((dδ0 / 4, dw0 / 5), (dδ0 / 20, dw0 / 25), (dδ0 / 100, dw0 / 100))
            _, δ0, w0 = best
            for w in w0-4dw:dw:w0+4dw, δ in δ0-4dδ:dδ:δ0+4dδ
                w > 0.05 || continue
                v = misfit(δ, w); v < best[1] && (best = (v, δ, w))
            end
        end
        isfinite(best[1]) || throw(ArgumentError("no fit in $a-$b nm (does the solar spectrum cover it?)"))
        push!(centres, (a + b) / 2); push!(dshift, best[2]); push!(dfwhm, best[3])
    end
    so = min(shift_order, length(centres) - 1); fo = min(fwhm_order, length(centres) - 1)
    shift_fit = polyfit(centres, dshift, so)
    fwhm_fit = polyfit(centres, dfwhm, fo)
    (centres=centres, shift=dshift, fwhm=dfwhm, shift_fit=shift_fit, fwhm_fit=fwhm_fit,
     lambda=λ .+ shift_fit.(λ))
end

# -----------------------------------------------------------------------------
# simulated spectra
# -----------------------------------------------------------------------------

"""
    simulate_spectrum(lambda, reference, cross_sections, scds; shift=0.0, polynomial=Float64[],
                      offset=0.0, noise=:none, rng=Random.default_rng()) -> Vector

A synthetic measured spectrum on the grid `lambda` from the reference I₀ (on `lambda`) and
the slant columns `scds` of the `cross_sections`:

    I(λ) = I₀(λ - Δ) · exp(-Σᵢ Sᵢ σᵢ(λ - Δ)) · exp(Σₖ aₖ (λ - λ̄)ᵏ) + offset

where Δ = `shift` (nm) moves everything the detector records (so a fit should find a
spectrum shift of -Δ), `polynomial` holds the aₖ of a broadband term around the grid's
mean wavelength λ̄, and `offset` adds a constant intensity. `noise` is `:none`, `:photon`
(Gaussian with σ = √I, counts) or a number (Gaussian with σ = noise·I). Use it to test fit
settings against known slant columns.
"""
function simulate_spectrum(lambda::AbstractVector{<:Real}, reference::AbstractVector{<:Real},
                           cross_sections::AbstractVector{CrossSection}, scds::AbstractVector{<:Real};
                           shift::Real=0.0, polynomial::AbstractVector{<:Real}=Float64[], offset::Real=0.0,
                           noise=:none, rng::AbstractRNG=Random.default_rng())
    λ = Float64.(lambda)
    length(reference) == length(λ) || throw(DimensionMismatch("reference must be on lambda"))
    length(cross_sections) == length(scds) || throw(DimensionMismatch("one slant column per cross section"))
    x = λ .- shift
    i0 = Float64.(reference)
    d2 = spline_deriv2(λ, i0)
    I = [spline_at(λ, i0, d2, v, true) for v in x]
    τ = zeros(length(λ))
    for (c, s) in zip(cross_sections, scds)
        τ .+= s .* xs_at(c, x)
    end
    I .*= exp.(-τ)
    if !isempty(polynomial)
        λm = sum(λ) / length(λ)
        I .*= exp.([sum(a * (l - λm)^(k - 1) for (k, a) in enumerate(polynomial)) for l in λ])
    end
    I .+= offset
    if noise === :photon
        I .+= sqrt.(max.(I, 0.0)) .* randn(rng, length(I))
    elseif noise isa Real
        noise > 0 && (I .+= noise .* abs.(I) .* randn(rng, length(I)))
    elseif noise !== :none
        throw(ArgumentError("noise must be :none, :photon or a number"))
    end
    I
end
