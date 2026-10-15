# Synthetic DOAS data set and QDOAS project for MiniDOAS.jl's examples and tests.
# Copyright (c) 2026 Sunandan Mahant <sunandanmahant@outlook.com>; BSD-3-Clause, see LICENSE.
#
# Everything here is made up: the cross sections are smooth analytic shapes with
# realistic magnitudes, NOT laboratory data, and must not be used for real retrievals.
# What is realistic is the processing chain: an MFC STD instrument, a dark/offset file,
# a solar-like reference with Fraunhofer lines, a Ring term, a wavelength shift, a
# broadband attenuation, photon noise, and a QDOAS project (.xml) that QDOAS itself can
# open. The true slant columns are known, so the fit can be checked against them.
#
# The data are the same as those of QDOASJulia's example, on which QDOAS 3.7.5 and the
# port agree to a median of 1.2e-12 of the fit errors.

using MiniDOAS, Printf

"Deterministic pseudo-random numbers (the same stream on every Julia version)."
mutable struct LCG
    s::UInt64
end
next!(r::LCG) = (r.s = r.s * 0x5851f42d4c957f2d + 0x14057b7ef767814f; (r.s >> 11) * (1.0 / 9007199254740992.0))
gauss!(r::LCG) = sqrt(-2log(max(next!(r), 1e-300))) * cos(2π * next!(r))

const NPIX = 1024
wavelength(i) = 295.0 + 0.05 * (i - 1) + 2.0e-6 * (i - 1)^2      # nm, i = 1..NPIX

"Solar-like reference: smooth continuum with ~160 Gaussian Fraunhofer lines."
function make_sun()
    r = LCG(0x51)
    # widths of 0.25-0.5 nm (5-10 pixels): lines as a spectrometer with ~0.5 nm resolution
    # sees them. Narrower lines would be undersampled and no interpolation could shift them.
    lines = [(290.0 + 62.0 * next!(r), 0.04 + 0.30 * next!(r), 0.25 + 0.25 * next!(r)) for _ in 1:160]
    λ -> begin
        c = 2.0e4 * (1 + (λ - 295.0) / 25.0) * exp(-((λ - 335.0) / 70.0)^2)
        for (λk, depth, w) in lines
            c *= 1 - depth * exp(-((λ - λk) / w)^2)
        end
        c
    end
end

# Synthetic cross sections (cm^2/molecule): shapes chosen to be distinct from each other.
xs_so2(λ) = 1e-19 * (0.4 + sum(a * exp(-((λ - (297.0 + 2.6k)) / 0.8)^2) for (k, a) in enumerate(
                [3.0, 4.5, 5.0, 4.0, 3.2, 2.4, 1.6, 1.0, 0.6, 0.3]))) * exp(-max(0.0, λ - 315.0) / 5.0)
xs_o3(λ) = 1e-19 * 25.0 * exp(-(λ - 295.0) / 6.5) * (1 + 0.12 * sin(2π * (λ - 295.0) / 3.3))
xs_no2(λ) = 1e-19 * 2.5 * (1 + 0.35 * sin(2π * (λ - 300.0) / 2.3) + 0.15 * sin(2π * (λ - 300.0) / 0.95)) *
            (1 / (1 + exp(-(λ - 318.0) / 3.0)))

"Ring-like term: filling-in of Fraunhofer lines, (smoothed sun - sun) / sun."
function make_ring(sun)
    grid = 288.0:0.01:352.0
    s = sun.(grid)
    k = [exp(-(x / 0.9)^2) for x in -3.0:0.01:3.0]; k ./= sum(k)
    h = length(k) ÷ 2
    sm = [sum(k[j] * s[clamp(i + j - h - 1, 1, length(s))] for j in eachindex(k)) for i in eachindex(s)]
    ring = (sm .- s) ./ s
    λ -> begin                                      # linear interpolation on the 0.01 nm grid
        x = clamp((λ - first(grid)) / 0.01 + 1, 1.0, length(grid) - 1e-9)
        j = floor(Int, x); f = x - j
        (1 - f) * ring[j] + f * ring[j+1]
    end
end

"True slant columns of spectrum k (a plume of SO2 drifting through the view)."
truth(k) = (SO2=2.0e16 + 4.0e17 * exp(-((k - 6.5) / 2.2)^2),
            O3=3.0e17 + 1.5e16 * k,
            NO2=6.0e15 + 2.0e15 * sin(k / 2),
            Ring=0.8 + 0.05 * k,
            shift=0.012)   # nm: pixel i records the light of λ_i - shift. QDOAS evaluates a
                           # shifted vector at λ - s, so a fitted "Spectrum" shift should
                           # come out as -shift and a shift of "Ref" plus cross sections as +shift

const SCANS = 200
dark_per_scan(i) = 80.0 + 5.0 * sin(i / 7.0)

# The spectra are written as QDOASJulia's example wrote them (counts with one decimal,
# `INT_TIME` in capitals), so that the data set stays the one QDOAS was compared on.
# `write_mfc_std` would write the same information.
function write_std(path, counts, name, t0, scans)
    open(path, "w") do io
        println(io, "GDBGMNUP"); println(io, 1); println(io, length(counts))
        for c in counts
            @printf(io, "%.1f\n", c)
        end
        println(io, name); println(io, "Synthetic spectrometer"); println(io, "QDOASJulia synthetic example")
        println(io, "09/29/2026")
        @printf(io, "%02d:%02d:%02d\n", t0 ÷ 3600, (t0 % 3600) ÷ 60, t0 % 60)
        t1 = t0 + 30
        @printf(io, "%02d:%02d:%02d\n", t1 ÷ 3600, (t1 % 3600) ÷ 60, t1 % 60)
        println(io, 0); println(io, 0)
        println(io, "SCANS $scans"); println(io, "INT_TIME 50.000000")
        println(io, "SITE Synthetic"); println(io, "LONGITUDE 0.000000"); println(io, "LATITUDE 0.000000")
    end
end

"""
    make_synthetic(dir; nspec=12) -> (project, spectra, truth, windows)

Writes the synthetic data set and its QDOAS project into `dir` (absolute paths inside
the project, as QDOAS expects) and returns the project path, the spectra, the true slant
columns per spectrum and the analysis windows (`FitWindow`s) of the project.
"""
function make_synthetic(dir::AbstractString; nspec::Int=12)
    dir = abspath(dir)
    for d in ("spectra", "cross", "calibration")
        mkpath(joinpath(dir, d))
    end
    sun = make_sun()
    ring = make_ring(sun)
    λ = wavelength.(1:NPIX)

    open(joinpath(dir, "calibration", "calibration.clb"), "w") do io
        foreach(l -> @printf(io, "%.8f\n", l), λ)
    end
    open(joinpath(dir, "cross", "reference.ref"), "w") do io        # I0 in counts, no noise, no dark
        for l in λ
            @printf(io, "%.8f\t%.3f\n", l, SCANS * sun(l))
        end
    end
    for (name, f) in (("so2_synthetic.xs", xs_so2), ("o3_synthetic.xs", xs_o3), ("no2_synthetic.xs", xs_no2),
                      ("ring_synthetic.xs", ring))
        open(joinpath(dir, "cross", name), "w") do io
            println(io, "; synthetic cross section for MiniDOAS.jl examples - not laboratory data")
            for l in 290.0:0.01:352.0
                @printf(io, "%.3f %.8e\n", l, f(l))
            end
        end
    end

    write_std(joinpath(dir, "spectra", "D0000000.STD"), [SCANS * dark_per_scan(i) for i in 1:NPIX], "D0000000", 43200, SCANS)
    noise = LCG(0x2026)
    files = String[]; truths = NamedTuple[]
    for k in 1:nspec
        t = truth(k)
        counts = map(1:NPIX) do i
            # a calibration drift moves everything the detector records: pixel i sees the
            # light of wavelength λ[i] - shift
            l = λ[i] - t.shift
            od = xs_so2(l) * t.SO2 + xs_o3(l) * t.O3 + xs_no2(l) * t.NO2 + ring(l) * t.Ring
            broadband = 0.85 * (1 - 0.004 * (l - 320.0))
            signal = SCANS * sun(l) * exp(-od) * broadband
            signal + sqrt(signal) * gauss!(noise) + SCANS * dark_per_scan(i)
        end
        name = @sprintf("S%07d", k)
        f = joinpath(dir, "spectra", name * ".STD")
        write_std(f, counts, name, 43200 + 60k, SCANS)
        push!(files, f); push!(truths, t)
    end
    open(joinpath(dir, "truth.csv"), "w") do io
        println(io, "spectrum,SO2,O3,NO2,Ring,shift_nm")
        for (f, t) in zip(files, truths)
            @printf(io, "%s,%.6e,%.6e,%.6e,%.6f,%.4f\n", basename(f), t.SO2, t.O3, t.NO2, t.Ring, t.shift)
        end
    end

    # The analysis, described in Julia and written as a QDOAS project
    cross(f) = joinpath(dir, "cross", f)
    so2, o3, no2, rg = (CrossSection(n, cross(f)) for (n, f) in
        (("SO2", "so2_synthetic.xs"), ("O3", "o3_synthetic.xs"), ("NO2", "no2_synthetic.xs"), ("Ring", "ring_synthetic.xs")))
    windows = [
        # one fitted shift of the measured spectrum, and a non-linear constant offset
        FitWindow("SO2", 303.0, 318.0, [so2, o3, rg]; polynomial=3, offset=0,
                  shifts=[ShiftStretch("Spectrum"), ShiftStretch("Ref"; shift=false),
                          ShiftStretch("SO2", "O3", "Ring"; shift=false)]),
        # the same drift, fitted the other way: one shift shared by the reference and the
        # cross sections
        FitWindow("NO2", 328.0, 345.0, [no2, o3, rg]; polynomial=3,
                  shifts=[ShiftStretch("Spectrum"; shift=false), ShiftStretch("Ref", "NO2", "O3", "Ring")]),
    ]
    dark = joinpath(dir, "spectra", "D0000000.STD")
    project = write_qdoas_project(joinpath(dir, "synthetic_project.xml"), windows;
        name="Synthetic", lambda=λ, reference=cross("reference.ref"),
        calibration=joinpath(dir, "calibration", "calibration.clb"),
        spectra=joinpath(dir, "spectra"), filter="S*.STD", dark=dark, offset=dark,
        output=joinpath(dir, "qdoas_results.nc"))
    (project=project, spectra=files, truth=truths, windows=windows)
end
