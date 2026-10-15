# Prepares the real-data example: MAX-DOAS spectra measured on the Antarctic
# Circumnavigation Expedition (ACE), analysed with public cross sections.
#
#   julia --project=/path/to/MiniDOAS.jl examples/ace_maxdoas/prepare.jl
#
# The spectra are in data/ (from Zenodo, CC BY 4.0; see README.md). This script
# downloads the cross sections and the solar spectrum from their public sources, and
# writes into work/:
#
#   * calibration.clb: the wavelength calibration of the instrument, and its slit width,
#     found by fitting the high-resolution solar spectrum to the noon zenith spectrum
#     (`calibrate_wavelength`, in the spirit of QDOAS's "Kurucz" calibration);
#   * reference.ref: the noon zenith spectrum on that calibration (the Fraunhofer reference);
#   * cross sections convolved with the instrument's slit function (`convolve`), and a
#     Ring spectrum computed from the solar spectrum with QDOAS's Ring algorithm
#     (`ring_spectrum`);
#   * ace_maxdoas.xml: a QDOAS project using all of it (`write_qdoas_project`), which QDOAS
#     itself also runs.
#
# Needs internet access once (about 7 MB); the downloads are kept in work/downloads/.
# Copyright (c) 2026 Sunandan Mahant <sunandanmahant@outlook.com>; BSD-3-Clause, see LICENSE.

using MiniDOAS, Downloads, Printf

const HERE = @__DIR__
const DATA = joinpath(HERE, "data")
const WORK = get(ENV, "MINIDOAS_ACE_WORK", joinpath(HERE, "work"))
const DOWNLOADS = joinpath(WORK, "downloads")

# -----------------------------------------------------------------------------
# public sources
# -----------------------------------------------------------------------------
const MAINZ = "https://uv-vis-spectral-atlas-mainz.org/uvvis_data/cross_sections/"

"""
Each cross section: its symbol, where it comes from, whether its wavelengths are in
vacuum (else in standard air) and the reference to cite. The analysis grid is in vacuum
wavelengths, because the calibration is made against the (vacuum) solar spectrum.
"""
const CROSS_SECTIONS = [
    (sym="NO2", file="no2_vandaele1998_294K.txt", vacuum=true,
     url=MAINZ * "Nitrogen%20oxides/NO2_Vandaele(1998)_294K_238-667nm.txt",
     cite="Vandaele, A. C. et al. (1998), J. Quant. Spectrosc. Radiat. Transfer 59, 171-184, doi:10.1016/S0022-4073(97)00168-4 (294 K)"),
    (sym="O3", file="o3_serdyuchenko2014_223K.txt", vacuum=true,
     url=MAINZ * "Ozone/O3_Serdyuchenko(2014)_223K_213-1100nm(2013%20version).txt",
     cite="Serdyuchenko, A. et al. (2014), Atmos. Meas. Tech. 7, 625-636, doi:10.5194/amt-7-625-2014 (223 K)"),
    (sym="O4", file="o4_thalman2013_293K.txt", vacuum=true,
     url=MAINZ * "Oxygen/O4_ThalmanVolkamer(2013)_293K_335.749-600.802nm.txt",
     cite="Thalman, R. and Volkamer, R. (2013), Phys. Chem. Chem. Phys. 15, 15371-15381, doi:10.1039/C3CP50968K (293 K)"),
    (sym="BrO", file="bro_fleischmann2004_223K.txt", vacuum=true,
     url=MAINZ * "Halogen%20oxides/Br%20oxides/BrO_Fleischmann(2004)_223K_300-385nm.txt",
     cite="Fleischmann, O. C. et al. (2004), J. Photochem. Photobiol. A 168, 117-132, doi:10.1016/j.jphotochem.2004.03.026 (223 K)"),
    (sym="HCHO", file="hcho_meller2000_298K.txt", vacuum=false,
     url=MAINZ * "Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O_MellerMoortgat(2000)_298K_224.56-376.00nm(0.01nm).txt",
     cite="Meller, R. and Moortgat, G. K. (2000), J. Geophys. Res. 105, 7089-7101, doi:10.1029/1999JD901074 (298 K)"),
]
const SOLAR = (file="sao2010.solref.converted", vacuum=true,
               url="https://www.cfa.harvard.edu/atmosphere/links/sao2010.solref.converted",
               cite="Chance, K. and Kurucz, R. L. (2010), J. Quant. Spectrosc. Radiat. Transfer 111, 1289-1295, doi:10.1016/j.jqsrt.2010.01.036")

function fetch(url, name)
    path = joinpath(DOWNLOADS, name)
    if !isfile(path)
        println("  downloading ", name)
        tmp = path * ".part"
        Downloads.download(url, tmp)
        mv(tmp, path; force=true)
    end
    path
end

"The first two numeric columns of a downloaded file (header lines are skipped), sorted."
function read_two(path)
    x = Float64[]; y = Float64[]
    for l in eachline(path)
        f = split(l)
        length(f) >= 2 || continue
        a = tryparse(Float64, f[1]); b = tryparse(Float64, f[2])
        (a === nothing || b === nothing) && continue
        push!(x, a); push!(y, b)
    end
    p = sortperm(x)
    x[p], y[p]
end

function main()
    mkpath(DOWNLOADS); mkpath(joinpath(WORK, "cross"))
    println("MiniDOAS.jl ACE MAX-DOAS example: preparing ", WORK)

    # instrument calibration (Hg-Ne lines, from the data set) and the noon zenith spectrum
    λinst = read_calibration(joinpath(DATA, "calibration", "ace_uv_instrument.clb"))
    refstd = read_mfc_std(only(filter(f -> endswith(f, ".STD"), readdir(joinpath(DATA, "reference"); join=true))))
    length(refstd.counts) == length(λinst) || error("reference and calibration sizes differ")

    sλ, sy = read_two(fetch(SOLAR.url, SOLAR.file))
    SOLAR.vacuum || (sλ = air_to_vacuum.(sλ))
    xs = Dict{String,CrossSection}()
    for c in CROSS_SECTIONS
        x, y = read_two(fetch(c.url, c.file))
        c.vacuum || (x = air_to_vacuum.(x))
        # some data sets have gaps between bands: keep the continuous part around 350 nm
        cuts = [0; findall(>(1.0), diff(x)); length(x)]
        k = findfirst(i -> x[cuts[i]+1] <= 350 <= x[cuts[i+1]], 1:length(cuts)-1)
        k === nothing && error("$(c.sym): no data around 350 nm")
        xs[c.sym] = CrossSection(c.sym, x[cuts[k]+1:cuts[k+1]], y[cuts[k]+1:cuts[k+1]])
    end

    # 1. calibration: shift and slit width in 8 sub-windows, then smooth functions of λ;
    #    a second pass adds a Ring spectrum at the first pass's resolution
    println("  calibrating the wavelengths and the slit width against the solar spectrum")
    wins = [(318.0 + 10k, 328.0 + 10k) for k in 0:7]
    o3first = convolve(xs["O3"], λinst, 0.55).sigma
    cal = calibrate_wavelength(λinst, refstd.counts, sλ, sy; windows=wins, absorbers=[o3first])
    ring0 = ring_spectrum(sλ, sy, cal.lambda, GaussianSlit(cal.fwhm_fit))
    cal = calibrate_wavelength(λinst, refstd.counts, sλ, sy; windows=wins,
                               absorbers=[convolve(xs["O3"], cal.lambda, GaussianSlit(cal.fwhm_fit)).sigma, ring0])
    for k in eachindex(cal.centres)
        @printf("    %5.1f nm: shift %+.4f nm (smooth %+.4f), slit FWHM %.3f nm (smooth %.3f)\n",
                cal.centres[k], cal.shift[k], cal.shift_fit(cal.centres[k]), cal.fwhm[k], cal.fwhm_fit(cal.centres[k]))
    end
    λ = cal.lambda
    slit = GaussianSlit(cal.fwhm_fit)

    calib = write_calibration(joinpath(WORK, "calibration.clb"), λ)
    reference = write_reference(joinpath(WORK, "reference.ref"), λ, refstd.counts)

    # 2. cross sections at the instrument's resolution, on the calibrated grid
    println("  convolving the cross sections")
    conv = Dict{String,CrossSection}()
    for c in CROSS_SECTIONS
        x = xs[c.sym].lambda
        keep = filter(v -> x[1] + 1.5 <= v <= x[end] - 1.5, λ)
        cx = convolve(xs[c.sym], keep, slit)
        file = write_cross_section(joinpath(WORK, "cross", c.sym * ".xs"), cx;
            header="$(c.sym): $(c.cite)\nconvolved with a Gaussian slit (FWHM $(round(cal.fwhm_fit(340); digits=3))-$(round(cal.fwhm_fit(390); digits=3)) nm), vacuum wavelengths")
        conv[c.sym] = CrossSection(c.sym, file)
    end
    println("  computing the Ring spectrum")
    ring = CrossSection("Ring", λ, ring_spectrum(sλ, sy, λ, slit))
    conv["Ring"] = CrossSection("Ring", write_cross_section(joinpath(WORK, "cross", "Ring.xs"), ring;
        header="Ring spectrum from the solar spectrum of $(SOLAR.cite),\nQDOAS Ring algorithm (Chance and Spurr 1997), 250 K, Gaussian slit"))

    # 3. the analysis windows, written as a QDOAS project
    project = joinpath(WORK, "ace_maxdoas.xml")
    write_qdoas_project(project, ace_windows(conv); name="ACE MAX-DOAS", lambda=λ, reference=reference,
                        calibration=calib, spectra=joinpath(DATA, "spectra"), filter="*.STD",
                        output=joinpath(WORK, "qdoas_results.nc"))
    println("  wrote ", project)
    project
end

"The analysis of the example: an NO2 and an O4 window."
function ace_windows(xs)
    no2 = FitWindow("NO2", 338.0, 370.0, [xs[s] for s in ("NO2", "O3", "O4", "BrO", "HCHO", "Ring")];
                    polynomial=5, linear_offset=1,
                    shifts=[ShiftStretch("Spectrum"; stretch=1), ShiftStretch("Ref"; shift=false),
                            ShiftStretch("NO2", "O3", "O4", "BrO", "HCHO", "Ring"; shift=false)])
    o4 = FitWindow("O4", 352.0, 384.0, [xs[s] for s in ("O4", "NO2", "O3", "BrO", "Ring")];
                   polynomial=5, linear_offset=1,
                   shifts=[ShiftStretch("Spectrum"; stretch=1), ShiftStretch("Ref"; shift=false),
                           ShiftStretch("O4", "NO2", "O3", "BrO", "Ring"; shift=false)])
    [no2, o4]
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
