# MiniDOAS.jl test suite. Run with:  julia --project -e 'using Pkg; Pkg.test()'
# Copyright (c) 2026 Sunandan Mahant <sunandanmahant@outlook.com>; BSD-3-Clause, see LICENSE.
#
# Uses the synthetic data set of examples/synthetic_data.jl (known true slant columns)
# and the real spectra of examples/ace_maxdoas/data. If the environment variable
# QDOAS_DOAS_CL points at QDOAS's doas_cl, the fit is also compared with QDOAS itself.

using Test, MiniDOAS, Dates, Random, LinearAlgebra
using NCDatasets, HTTP, JSON3
include(joinpath(@__DIR__, "..", "examples", "synthetic_data.jl"))

const DIR = mktempdir()
const SYN = make_synthetic(DIR; nspec=6)
const P = parse_config(SYN.project)
const ACE = joinpath(@__DIR__, "..", "examples", "ace_maxdoas", "data")

mean(v) = sum(v) / length(v)
cor(a, b) = (x = a .- mean(a); y = b .- mean(b); sum(x .* y) / sqrt(sum(abs2, x) * sum(abs2, y)))
"Difference of two angles in degrees, in -180..180."
angdiff(a, b) = mod(a - b + 180, 360) - 180

"Every number of two window results, compared exactly."
function identical(w1::WindowResult, w2::WindowResult)
    all(f -> isequal(getfield(w1, f), getfield(w2, f)),
        (:rms, :chisqr, :niter, :residual, :slcol, :slerr, :shift, :shift_err, :stretch, :stretch_err,
         :params, :param_err, :pixels, :columns)) &&
        isequal(last.(w1.nc_vars), last.(w2.nc_vars))
end

@testset "MiniDOAS.jl" begin

@testset "version" begin
    toml = read(joinpath(@__DIR__, "..", "Project.toml"), String)
    @test occursin("version = \"$(MiniDOAS.VERSION_STRING)\"", toml)
end

@testset "QDOAS project file" begin
    @test P.project_name == "Synthetic"
    @test P.swath_name == "QDOAS Results"
    @test [w.name for w in P.windows] == ["SO2", "NO2"]
    @test length(P.spectra) == 6                          # the dark file is excluded
    @test all(f -> occursin(r"S\d{7}\.STD$", f), P.spectra)
    @test P.windows[1].poly_order == 3
    @test [c.name for c in P.windows[1].xs] == ["SO2", "O3", "Ring"]
    @test !MiniDOAS.isloaded(P.windows[1].xs[1])              # read when a fit is set up
    @test project_spectra(SYN.project) == P.spectra
    @test sprint(show, P) == "ProjectSpec(\"Synthetic\", 1024 pixels, windows SO2, NO2, 6 spectra)"
end

@testset "recovers the true slant columns" begin
    for (f, t) in zip(SYN.spectra, SYN.truth)
        r = analyze(P; spectrum=f)
        so2, no2 = r.windows
        @test r["SO2"] === so2 && r[2] === no2 && keys(r) == ["SO2", "NO2"]
        @test abs(so2.slcol["SO2"] - t.SO2) < 4 * so2.slerr["SO2"]
        @test abs(so2.slcol["O3"] - t.O3) < 4 * so2.slerr["O3"]
        @test abs(no2.slcol["NO2"] - t.NO2) < 4 * no2.slerr["NO2"]
        @test abs(so2.shift["Spectrum"] + t.shift) < 1e-3   # spectrum shifted back by -shift
        @test abs(no2.shift["Ref"] - t.shift) < 1e-3        # reference + cross sections shifted by +shift
        @test so2.rms < 1.5e-3 && no2.rms < 1.5e-3            # photon noise, no systematic misfit
        @test so2.niter >= 1
    end
end

@testset "batch, threads and determinism" begin
    b1 = analyze_spectra(P; threads=false)
    b2 = analyze_spectra(P; threads=true)
    @test b1 isa DOASBatch && b2 isa DOASBatch
    @test all(x -> x.result !== nothing, b1)
    @test [x.file for x in b1] == P.spectra
    for (x, y) in zip(b1, b2), (w1, w2) in zip(x.result.windows, y.result.windows)
        for c in w1.columns
            @test isapprox(w1.slcol[c], w2.slcol[c]; rtol=1e-10, atol=1e-10 * abs(w1.slerr[c]))
        end
    end
    bad = analyze_spectra(P, [P.spectra[1], joinpath(DIR, "missing.STD")])
    @test bad[1].result !== nothing && bad[2].result === nothing && !isempty(bad[2].error)
end

@testset "fit curves" begin
    r = analyze(P; spectrum=P.spectra[3], curves=true)
    cv = r.windows[1].curves
    @test length(cv.lambda) == length(cv.od) == length(cv.model) == length(cv.residual)
    @test cv.od ≈ cv.model .+ cv.residual
    @test sum(c -> c.value, cv.components) ≈ cv.model
    @test all(diff(cv.lambda) .> 0)
    @test r.windows[1].rms ≈ sqrt(sum(abs2, cv.residual) / length(cv.residual)) rtol = 1e-10
    @test analyze(P; spectrum=P.spectra[3]).windows[1].curves === nothing
end

@testset "fitting from Julia" begin
    # the project's windows with their cross sections in memory, and the reference as data
    λ, I0 = read_reference(P.windows[1].refone)
    windows = [MiniDOAS.rewindow(w; xs=[read_cross_section(c.file; name=c.name) for c in w.xs]) for w in P.windows]
    @test all(c -> MiniDOAS.isloaded(c), windows[1].xs)
    m = DOASModel(λ, I0, windows)
    @test sprint(show, m) == "DOASModel(1024 pixels, 295.00-348.24 nm, windows SO2, NO2)"
    d = read_mfc_std(P.dark)
    for f in P.spectra
        s = read_mfc_std(f)
        counts = correct_dark_offset(s; dark=d, offset=d)
        @test counts == MiniDOAS.load_spectrum(MiniDOAS.respec(P; spectrum=f))[1]
        r1 = analyze(P; spectrum=f)
        r2 = doas_fit(m, counts)
        @test all(identical(a, b) for (a, b) in zip(r1.windows, r2.windows))      # bit for bit
        @test r2.record === nothing && r2.spectrum == ""
    end
    # an MFCSpectrum carries its name, date and time into the result
    s = read_mfc_std(P.spectra[2])
    r = doas_fit(m, s)
    @test r.record === s && r.spectrum == "S0000002"
    @test measurement_time(s) == DateTime(2026, 9, 29, 12, 2, 15)
    # batches: a vector of spectra or a matrix. Bit for bit between serial runs; in
    # parallel, BLAS runs on one thread, which can change errors in the last digit.
    sp = [correct_dark_offset(read_mfc_std(f); dark=d, offset=d) for f in P.spectra]
    b = doas_fit(m, sp; threads=false)
    @test length(b) == 6 && all(x -> x.result !== nothing, b) && b[1].file == "spectrum 1"
    @test identical(b[4].result["NO2"], doas_fit(m, sp[4])["NO2"])
    bm = doas_fit(m, reduce(hcat, sp); threads=false)
    @test identical(bm[6].result["SO2"], b[6].result["SO2"])
    bp = doas_fit(m, sp)
    for (x, y) in zip(b, bp), (w1, w2) in zip(x.result.windows, y.result.windows), c in w1.columns
        @test isapprox(w1.slcol[c], w2.slcol[c]; rtol=1e-12) && isapprox(w1.slerr[c], w2.slerr[c]; rtol=1e-12)
    end
    # the convenience form
    @test identical(doas_fit(λ, sp[1], I0, windows[2])["NO2"], b[1].result["NO2"])
    # errors
    @test_throws DimensionMismatch doas_fit(m, sp[1][1:100])
    @test_throws EngineError doas_fit(m, I0)                       # identical to the reference
    bad = doas_fit(m, [sp[1], fill(-1.0, 1024)])
    @test bad[1].result !== nothing && bad[2].result === nothing && occursin("non-positive", bad[2].error)
    @test_throws ArgumentError DOASModel(reverse(λ), I0, windows)
    # two identical cross sections cannot both be fitted
    twin = windows[2].xs[1]
    same = FitWindow("NO2", 328.0, 345.0, [twin, CrossSection("NO2b", twin.lambda, twin.sigma)])
    e = try doas_fit(λ, sp[1], I0, same) catch err; err end
    @test e isa EngineError && occursin("linearly dependent", e.msg)
end

@testset "single spectra and multiple files" begin
    dir = dirname(P.spectra[1])
    # which files a path names
    @test spectrum_files(P.spectra[1]) == [P.spectra[1]]
    @test spectrum_files(joinpath(dir, "S*.STD")) == P.spectra
    @test spectrum_files(joinpath(dir, "s*.std")) == P.spectra                # names match whatever their case
    @test spectrum_files(dir; filter="S*.STD") == P.spectra
    @test spectrum_files(dir; filter="S*.STD;D*.STD") == spectrum_files(dir)  # the dark file too
    @test length(spectrum_files(dir)) == 7
    @test issubset(P.spectra, spectrum_files(DIR; filter="S*.STD", recursive=true))
    @test_throws ArgumentError spectrum_files(joinpath(DIR, "nothing_here"))
    @test_throws ArgumentError spectrum_files(joinpath(DIR, "nothing_here", "*.STD"))

    # a QDOAS project: one file, a pattern, a folder, a list
    r = analyze(P, P.spectra[2])
    @test identical(r["NO2"], analyze(P; spectrum=P.spectra[2])["NO2"])
    @test identical(analyze(SYN.project, P.spectra[2])["SO2"], r["SO2"])
    b = analyze(P, joinpath(dir, "S*.STD"); threads=false)
    @test b isa DOASBatch && length(b) == 6 && [x.file for x in b] == P.spectra
    @test identical(b[2].result["NO2"], r["NO2"])
    @test analyze(P, P.spectra[1:2]; threads=false) isa DOASBatch
    @test length(analyze_spectra(P, joinpath(dir, "S*.STD"))) == 6
    @test length(analyze_spectra(SYN.project, dir; filter="S*.STD")) == 6
    @test analyze_spectra(P) isa DOASBatch

    # from Julia: the same, with the dark and offset files, bit for bit
    λ, I0 = read_reference(P.windows[1].refone)
    m = DOASModel(λ, I0, P.windows)
    r1 = doas_fit(m, P.spectra[2]; dark=P.dark, offset=P.offset)
    @test all(identical(x, y) for (x, y) in zip(r1.windows, r.windows))
    @test r1.spectrum == P.spectra[2] && r1.record.name == "S0000002"
    bj = doas_fit(m, joinpath(dir, "S*.STD"); dark=P.dark, offset=P.offset, threads=false)
    @test [x.file for x in bj] == P.spectra && identical(bj[2].result["SO2"], r["SO2"])
    @test length(doas_fit(m, P.spectra; dark=read_mfc_std(P.dark), offset=read_mfc_std(P.offset))) == 6
    @test length(doas_fit(m, dir; filter="S*.STD")) == 6

    # text files: one column of intensities, or wavelength and intensity
    counts = correct_dark_offset(read_mfc_std(P.spectra[2]); dark=read_mfc_std(P.dark), offset=read_mfc_std(P.offset))
    t2 = write_columns(joinpath(DIR, "spectrum.txt"), λ, counts)
    t1 = write_columns(joinpath(DIR, "spectrum1.dat"), counts)
    @test identical(doas_fit(m, t2)["NO2"], r["NO2"]) && identical(doas_fit(m, t1)["NO2"], r["NO2"])
    @test doas_fit(m, t2).record === nothing && doas_fit(m, t2).spectrum == t2
    @test_throws ArgumentError doas_fit(m, t2; dark=P.dark)                # dark and offset are for MFC files
    @test_throws DimensionMismatch doas_fit(m, write_columns(joinpath(DIR, "short.txt"), counts[1:100]))
    @test_throws ArgumentError doas_fit(m, joinpath(DIR, "nope.STD"))
    # an MFC STD file is recognised without the .STD extension
    cp(P.spectra[2], joinpath(DIR, "spectrum_mfc.dat"); force=true)
    @test identical(doas_fit(m, joinpath(DIR, "spectrum_mfc.dat"); dark=P.dark, offset=P.offset)["SO2"], r["SO2"])

    # the batch
    @test bj["S0000003.STD"].result.spectrum == P.spectra[3] && bj[P.spectra[3]] === bj[3]
    @test_throws KeyError bj["S9999999.STD"]
    v, e = slant_column(bj, "SO2")
    @test v[2] == r["SO2"].slcol["SO2"] && e[2] == r["SO2"].slerr["SO2"]
    @test slant_column(r, "NO2", "O3") == (r["NO2"].slcol["O3"], r["NO2"].slerr["O3"])
    @test_throws ArgumentError slant_column(bj, "SO2", "NO2")
    @test_throws ArgumentError slant_column(r, "SO2", "NO2")
    @test measurement_time(bj)[2] == DateTime(2026, 9, 29, 12, 2, 15) == measurement_time(r)
    bad = doas_fit(m, [P.spectra[1], joinpath(dir, "missing.STD")])
    @test isnan(slant_column(bad, "SO2")[1][2]) && ismissing(measurement_time(bad)[2])
    txt = sprint(show, MIME("text/plain"), bad)
    @test occursin("2 spectra, 1 analysed, 1 failed", txt) && occursin("missing.STD", txt) && occursin("SO2: SO2", txt)
    @test sprint(show, bj) == "DOASBatch(6 spectra, 6 analysed)"
    @test bj.seconds > 0
    @test results_columns(bj).file == P.spectra

    # one NetCDF file per analysed spectrum
    @test basename.(write_netcdf(joinpath(DIR, "batch_nc"), bad)) == ["S0000001.nc"]
    @test basename.(write_netcdf(joinpath(DIR, "batch_nc"), bj)) == ["S000000$k.nc" for k in 1:6]
end

@testset "fit window definitions" begin
    no2 = CrossSection("NO2", [300.0, 301.0, 302.0], [1.0, 2.0, 3.0])
    @test sprint(show, no2) == "CrossSection(\"NO2\", 300.00-302.00 nm, 3 points)"
    @test CrossSection("x", [3.0, 2.0, 1.0], [30.0, 20.0, 10.0]).lambda == [1.0, 2.0, 3.0]  # reversed
    @test_throws DimensionMismatch CrossSection("x", [1.0, 2.0], [1.0])
    w = FitWindow("W", 310, 300, [no2]; polynomial=2, linear_offset=1, gaps=[(306, 305)],
                  shifts=ShiftStretch("Spectrum"; stretch=2))
    @test (w.lmin, w.lmax, w.lambda0) == (300.0, 310.0, 305.0) && w.gaps == [(305.0, 306.0)]
    @test w.poly_order == 2 && w.offl_order == 1 && w.shifts[1].st_type == 2
    @test_throws ArgumentError FitWindow("W", 300, 310, [no2]; offset=0, linear_offset=0)
    @test_throws ArgumentError FitWindow("W", 300, 310, [no2]; shifts=ShiftStretch("O3"))
    @test_throws ArgumentError FitWindow("W", 300, 310, [no2, no2])
    @test_throws ArgumentError FitWindow("W", 300, 310, [no2]; polynomial=9)
    @test_throws ArgumentError ShiftStretch("Spectrum"; stretch=3)
    g = ShiftStretch("Ref", "NO2"; shift=false, shift_init=0.01)
    @test g.symbols == ["Ref", "NO2"] && !g.sh_fit && g.sh_init == 0.01
    @test occursin("fixed at 0.01 nm", sprint(show, g))
end

@testset "QDOAS project writing" begin
    # the synthetic project itself is written by write_qdoas_project (examples/synthetic_data.jl)
    q = P.windows[1]
    @test q.off[1][1] && !q.off[2][1] && [g.symbols for g in q.shifts] == [["Spectrum"], ["Ref"], ["SO2", "O3", "Ring"]]
    # a model from Julia, with in-memory cross sections, a gap, stretch, λ0 and linear interpolation
    λ, I0 = read_reference(P.windows[1].refone)
    xs = [CrossSection(c.name, read_cross_section(c.file).lambda, read_cross_section(c.file).sigma) for c in P.windows[2].xs]
    w = FitWindow("NO2", 328.0, 345.0, xs; polynomial=4, linear_offset=0, linear_offset_ref=true,
                  gaps=[(335.0, 335.5)], lambda0=336.0, shifts=[ShiftStretch("Ref", "NO2", "O3", "Ring"; stretch=1)])
    m = DOASModel(λ, I0, w; interpolation=:linear, security_gap=6, convergence=1e-5, name="from Julia")
    xml = write_qdoas_project(joinpath(DIR, "from_julia.xml"), m; spectra=dirname(P.spectra[1]), filter="S*.STD",
                              dark=P.dark, offset=P.offset)
    p2 = parse_config(xml)
    @test p2.project_name == "from Julia" && p2.interpolation === :linear && p2.security_gap == 6
    @test p2.convergence == 1e-5 && p2.spectra == P.spectra
    w2 = only(p2.windows)
    @test (w2.lmin, w2.lmax, w2.lambda0, w2.gaps) == (w.lmin, w.lmax, w.lambda0, w.gaps)
    @test (w2.poly_order, w2.offl_order, w2.offl_ref) == (4, 0, true)
    @test [g.symbols for g in w2.shifts] == [["Ref", "NO2", "O3", "Ring"]] && w2.shifts[1].st_type == 1
    @test all(isfile(c.file) for c in w2.xs)               # in-memory cross sections were written
    d = read_mfc_std(P.dark)
    for f in P.spectra[1:3]
        @test identical(analyze(p2; spectrum=f)["NO2"],
                        doas_fit(m, correct_dark_offset(read_mfc_std(f); dark=d, offset=d))["NO2"])
    end
    @test_throws ArgumentError write_qdoas_project(joinpath(DIR, "x.xml"), DOASModel(P))   # no reference given
end

@testset "outputs" begin
    r = analyze(P; spectrum=P.spectra[2])
    nc = joinpath(DIR, "out.nc")
    write_netcdf(nc, r; config_path=SYN.project)
    NCDataset(nc) do ds
        g = ds.group["QDOAS Results"]
        @test g["Name"][1] == "S0000002"
        @test Int.(g["Date (DD-MM-YYYY)"][:, 1]) == [2026, 9, 29]
        for wr in r.windows
            gw = g.group[wr.name]
            for c in wr.columns
                haskey(gw, "SlCol($c)") && @test gw["SlCol($c)"][1] == wr.slcol[c]
            end
            @test gw["RMS"][1] == wr.rms
            @test Int(gw["iter"][1]) == wr.niter
            @test length(gw["residual_spectrum"][:, 1]) == length(wr.residual)
        end
    end
    back = read_qdoas_netcdf(nc)
    @test back.group == "QDOAS Results" && back.fields["Name"] == ["S0000002"]
    @test back.windows["NO2"]["SlCol(NO2)"] == [r["NO2"].slcol["NO2"]]
    @test back.windows["SO2"]["iter"] == [r["SO2"].niter]
    @test back.windows["SO2"]["attributes"]["fitting window range"] == "303.000 : 318.000"
    # a fit of an in-memory spectrum has no name, date or time
    λ, I0 = read_reference(P.windows[1].refone)
    rm2 = doas_fit(DOASModel(P), r.raw)
    write_netcdf(nc, rm2)
    @test !haskey(read_qdoas_netcdf(nc).fields, "Name")

    csv = write_csv(joinpath(DIR, "out.csv"), analyze_spectra(P))
    lines = readlines(csv)
    @test length(lines) == 1 + length(P.spectra)
    @test occursin("SO2.SlCol(SO2),SO2.SlErr(SO2)", lines[1])
    header, rows = results_table(r)
    @test length(rows) == 1 && rows[1][5] == "ok"
    @test occursin("latitude,longitude", lines[1])
    @test sprint(print_summary, r) |> x -> occursin("SO2", x) && occursin("RMS", x)
    @test occursin("Offset0", sprint(print_summary, r))
    # columns for Tables.jl / DataFrames
    cols = results_columns(analyze_spectra(P); sza=true)
    @test cols.file == P.spectra && cols.status == fill("ok", 6)
    @test cols[Symbol("SO2.SlCol(SO2)")] isa Vector{Float64}
    @test cols.date[1] == "2026-09-29" && all(isfinite, cols.sza)
    @test results_columns(r).name == ["S0000002"]
    @test sprint(show, MIME("text/plain"), r) == sprint(print_summary, r)
end

@testset "command line" begin
    out = joinpath(DIR, "cli")
    @test MiniDOAS.main(["-c", SYN.project, "-o", out, "-q"]) == 0
    @test sort(filter(f -> endswith(f, ".nc"), readdir(out))) == [splitext(basename(f))[1] * ".nc" for f in P.spectra]
    @test isfile(joinpath(out, "results.csv"))
    @test MiniDOAS.main(["-c", SYN.project, "-a", "NotThisProject", "-q"]) == 2
    # a spectrum that cannot be read is reported and makes the exit status 1
    @test MiniDOAS.main(["-c", SYN.project, "-f", joinpath(DIR, "no_such.STD"), "-o", out, "--no-netcdf", "-q"]) == 1
    @test MiniDOAS.main(["-c", SYN.project, "-f", joinpath(dirname(P.spectra[1]), "S00000*.STD"), "-o", joinpath(DIR, "cli2"), "-q"]) == 0
    @test length(readlines(joinpath(DIR, "cli2", "results.csv"))) == 1 + length(P.spectra)
    @test MiniDOAS.main(["--port", "x", "--gui"]) == 2
    @test MiniDOAS.main(["-V"]) == 0
end

@testset "options that are not implemented are refused" begin
    xml = read(SYN.project, String)
    function refused(from, to)
        f = joinpath(DIR, "variant.xml")
        occursin(from, xml) || error("pattern $from not in the project")
        write(f, replace(xml, from => to; count=1))
        try
            parse_config(f); false
        catch e
            e isa UnsupportedConfig
        end
    end
    @test refused("method=\"ODF\"", "method=\"ML+SVD\"")
    @test refused("kurucz=\"none\"", "kurucz=\"ref\"")
    @test refused("refsel=\"file\"", "refsel=\"auto\"")
    @test refused("cstype=\"interp\"", "cstype=\"convolute\"")
    @test refused("format=\"mfcstd\"", "format=\"ascii\"")
    @test refused("<lowpass_filter selected=\"none\">", "<lowpass_filter selected=\"boxcar\">")
end

@testset "files" begin
    # MFC STD: write and read back
    f = joinpath(DIR, "written.STD")
    t0 = DateTime(2024, 5, 1, 11, 59, 50)
    write_mfc_std(f, [100.0, 200.5, 3.0e6]; name="zenith", start=t0, stop=t0 + Second(20), scans=7,
                  exposure=1.25, elevation=90, azimuth=45.5, latitude=50.8, longitude=4.35, extra=["Note" => "test"])
    s = read_mfc_std(f)
    @test s.counts == [100.0, 200.5, 3.0e6] && s.name == "zenith" && s.scans == 7
    @test s.exposure == 1.25f0 && s.date == (2024, 5, 1) && s.time == (12, 0, 0)
    @test (s.elevation, s.azimuth, s.latitude, s.longitude) == (90.0, 45.5, 50.8, 4.35)
    @test "Note = test" in s.header
    @test sprint(show, s) == "MFCSpectrum(\"zenith\", 3 pixels, 2024-05-01 12:00:00, 7 scans, elevation 90.0°)"
    # an MFCSpectrum is written back as it was read
    for g in (P.spectra[1], joinpath(ACE, "reference", only(readdir(joinpath(ACE, "reference")))))
        a = read_mfc_std(g)
        write_mfc_std(joinpath(DIR, "copy.STD"), a)
        b = read_mfc_std(joinpath(DIR, "copy.STD"))
        @test b.counts == a.counts && b.header == a.header && b.time == a.time
    end
    @test read_mfc_std(f; npixels=5).counts == [100.0, 200.5, 3.0e6, 0.0, 0.0]
    # dark and offset as QDOAS applies them
    off = MiniDOAS.MFCSpectrum([10.0, 10.0, 10.0], 2, 0f0, "", (0, 0, 0), (0, 0, 0), NaN, NaN, NaN, NaN, String[])
    drk = MiniDOAS.MFCSpectrum([30.0, 40.0, 50.0], 2, 2f0, "", (0, 0, 0), (0, 0, 0), NaN, NaN, NaN, NaN, String[])
    # I = R - O nR/nO - (D - O nD/nO) nR tR / (nD tD), with nR = 7, tR = 1.25
    @test correct_dark_offset(s; offset=off, dark=drk) ≈ s.counts .- 10 * 7 / 2 .- ([30, 40, 50] .- 10) .* 7 * 1.25 / 4
    @test correct_dark_offset(s) == s.counts
    # columns, references, calibrations, cross sections
    g = joinpath(DIR, "cols.txt")
    write(g, "; comment\n# another\nlambda value\n300.0, 1.5\n301.0\t2.5\n\n302.0 3.5\n")
    @test read_columns(g) == [[300.0, 301.0, 302.0], [1.5, 2.5, 3.5]]
    write(g, "1 2\n3\n")
    @test_throws ArgumentError read_columns(g)
    x = [300.1, 300.2, 300.3]; y = [0.1, 1e-20, 3.0]
    @test read_reference(write_reference(g, x, y; header="a\nb")) == (x, y)
    @test readline(g) == "; a"
    @test read_calibration(write_calibration(g, x)) == x
    xs = CrossSection("NO2", x, y)
    back = read_cross_section(write_cross_section(g, xs))
    @test back.lambda == x && back.sigma == y && back.name == "cols" && back.file == g
    @test read_cross_section(g; name="NO2").name == "NO2"
end

@testset "spectral tools" begin
    # convolution of a fine line by a Gaussian slit: area kept, width as expected
    x = collect(300.0:0.001:310.0)
    line = [abs(v - 305.0) < 0.0005 ? 1.0 : 0.0 for v in x]
    at = collect(303.0:0.01:307.0)
    c = convolve(x, line, at, 0.5)
    @test sum(c) * 0.01 ≈ 0.001 rtol = 1e-3                     # area of the line
    half = at[c .>= maximum(c) / 2]
    @test last(half) - first(half) ≈ 0.5 atol = 0.021             # FWHM, on a 0.01 nm grid
    @test convolve(x, ones(length(x)), at, GaussianSlit(λ -> 0.1 + 0.01(λ - 300))) ≈ ones(length(at))
    for s in (BoxcarSlit(0.4), TriangularSlit(0.3), LorentzianSlit(0.2), TabulatedSlit([-0.2, 0.0, 0.2], [0.0, 1.0, 0.0]))
        cs = convolve(x, line, at, s)
        @test sum(cs) * 0.01 ≈ 0.001 rtol = 0.05
        @test sum(at .* cs) / sum(cs) ≈ 305.0 atol = 0.005           # centred on the line
    end
    @test convolve(x, line, [500.0], 0.5)[1] |> isnan          # outside the data
    @test_throws ArgumentError convolve(reverse(x), line, at, 0.5)
    xs = CrossSection("X", x, sin.(x))
    cx = convolve(xs, at, 0.3)
    @test cx.name == "X" && cx.lambda == at
    @test cx.sigma ≈ sin.(at) .* exp(-(0.3 / (2sqrt(2log(2))))^2 / 2) rtol = 2e-3
    # I0 correction: for a small slant column it tends to the I0-weighted convolution
    # (the plain one under a flat I0), and saturation lowers it for a large one
    solar = [1 + 0.8 * cos(2π * v / 0.15) for v in x]
    σ = CrossSection("Y", x, [1e-19 * (1 + sin(2π * v / 0.37)) for v in x])
    weighted = convolve(x, solar .* σ.sigma, at, 0.5) ./ convolve(x, solar, at, 0.5)
    small = convolve_i0(σ, x, solar, at, 0.5; scd=1e13).sigma          # optical depths ~1e-6
    @test small ≈ weighted rtol = 1e-5
    @test convolve_i0(σ, x, ones(length(x)), at, 0.5; scd=1e13).sigma ≈ convolve(σ, at, 0.5).sigma rtol = 1e-5
    @test all(convolve_i0(σ, x, solar, at, 0.5; scd=1e19).sigma .< small)
    # air and vacuum (Edlén / Birch and Downs): 400.000 and 500.000 nm in air are
    # 400.113 and 500.139 nm in vacuum
    @test air_to_vacuum(400.0) ≈ 400.113 atol = 1e-3
    @test air_to_vacuum(500.0) ≈ 500.139 atol = 1e-3
    @test vacuum_to_air(air_to_vacuum(350.0)) ≈ 350.0 atol = 1e-9
    # smoothing and differential spectra
    y = collect(1.0:10.0)
    @test smooth(y, 3) ≈ [1.5; 2.0:9.0; 9.5]
    @test smooth(fill(2.0, 9), 5; kernel=:gaussian) ≈ fill(2.0, 9)
    @test smooth([0, 0, 1.0, 0, 0], 3; kernel=:binomial) ≈ [0, 0.25, 0.5, 0.25, 0]
    @test smooth([0, 0, 1.0, 0, 0], 3; kernel=:triangular) ≈ [0, 0.25, 0.5, 0.25, 0]
    @test_throws ArgumentError smooth(y, 4)
    lam = collect(300.0:0.1:320.0)
    @test maximum(abs, differential(lam, 2 .+ 0.3 .* (lam .- 310) .+ 0.01 .* (lam .- 310) .^ 3)) < 1e-9
    d = differential(lam, sin.(lam); order=2)
    @test abs(sum(d)) < 1e-8
end

@testset "Ring spectrum and wavelength calibration" begin
    # a solar-like spectrum at 0.01 nm and an instrument with a 0.5 nm Gaussian slit
    sun = make_sun()
    hl = collect(296.0:0.01:350.0); hs = sun.(hl)
    λ = collect(300.0:0.05:345.0)
    ring = ring_spectrum(hl, hs, λ, 0.5)
    @test length(ring) == length(λ) && all(isfinite, ring)
    @test abs(mean(ring) - 1) < 0.1                     # Raman / solar: near 1
    # Raman filling-in: the Ring spectrum peaks inside Fraunhofer lines
    i0 = convolve(hl, hs, λ, 0.5)
    @test cor(differential(λ, ring), -differential(λ, log.(i0))) > 0.5
    # calibration: the instrument records λ + 0.05 nm and has a 0.45 nm slit
    meas = convolve(hl, hs, λ .+ 0.05, 0.45)
    cal = calibrate_wavelength(λ, meas, hl, hs; windows=[(305.0 + 10k, 315.0 + 10k) for k in 0:2],
                               shift_order=1, fwhm_order=0)
    @test all(abs.(cal.shift .- 0.05) .< 0.002) && all(abs.(cal.fwhm .- 0.45) .< 0.01)
    @test maximum(abs, cal.lambda .- (λ .+ 0.05)) < 0.002
    @test cal.fwhm_fit(320.0) ≈ 0.45 atol = 0.01
end

@testset "simulated spectra" begin
    sun = make_sun()
    λ = wavelength.(1:NPIX)
    I0 = 1e4 .* sun.(λ)
    so2 = read_cross_section(P.windows[1].xs[1].file; name="SO2")
    o3 = read_cross_section(P.windows[1].xs[2].file; name="O3")
    s0 = simulate_spectrum(λ, I0, [so2, o3], [2e17, 3e17])
    w = FitWindow("SO2", 305.0, 318.0, [so2, o3]; polynomial=2)
    r = doas_fit(λ, s0, I0, w)["SO2"]
    @test r.slcol["SO2"] ≈ 2e17 rtol = 1e-6
    @test r.slcol["O3"] ≈ 3e17 rtol = 1e-6
    # with photon noise and a drift of the spectrometer, the shift of the reference and
    # cross sections is recovered, and the errors describe the scatter
    w = FitWindow("SO2", 305.0, 318.0, [so2, o3]; polynomial=2, shifts=ShiftStretch("Ref", "SO2", "O3"))
    m = DOASModel(λ, I0, w)
    pulls = Float64[]; shifts = Float64[]
    for seed in 1:12
        s = simulate_spectrum(λ, I0, [so2, o3], [2e17, 3e17]; shift=0.01, polynomial=[0.0, -0.01],
                              noise=:photon, rng=MersenneTwister(seed))
        q = doas_fit(m, s)["SO2"]
        push!(pulls, (q.slcol["SO2"] - 2e17) / q.slerr["SO2"]); push!(shifts, q.shift["Ref"])
    end
    @test abs(sum(shifts) / 12 - 0.01) < 0.002
    @test sqrt(sum(abs2, pulls) / 12) < 2.0
    @test simulate_spectrum(λ, I0, CrossSection[], Float64[]; offset=5.0) ≈ I0 .+ 5
    @test_throws ArgumentError simulate_spectrum(λ, I0, CrossSection[], Float64[]; noise=:pink)
end

@testset "solar position and air mass factors" begin
    # declination at the June solstice and March equinox of 2024
    @test solar_position(DateTime(2024, 6, 20, 20, 51), 0, 0).declination ≈ 23.44 atol = 0.01
    @test abs(solar_position(DateTime(2024, 3, 20, 3, 6), 0, 0).declination) < 0.01
    # equation of time near its extremes (about -14.2 min mid-February, +16.4 min early November)
    @test solar_position(DateTime(2025, 2, 11, 12), 0, 0).equation_of_time ≈ -14.2 atol = 0.2
    @test solar_position(DateTime(2025, 11, 3, 12), 0, 0).equation_of_time ≈ 16.4 atol = 0.2
    # at true solar noon the sun is due south (north) and the zenith angle is |φ - δ|
    for (lat, lon) in ((50.8, 4.35), (-34.6, -58.4), (40.0, -105.0))
        t = DateTime(2025, 4, 7, 12)
        sp = solar_position(t, lat, lon)
        noon = t - Millisecond(round(Int, (sp.equation_of_time + 4lon) * 60_000))
        q = solar_position(noon, lat, lon)
        @test q.zenith ≈ abs(lat - q.declination) atol = 0.02
        @test abs(angdiff(q.azimuth, lat > q.declination ? 180.0 : 0.0)) < 1.0
        @test q.elevation == 90 - q.zenith
    end
    # morning sun in the east, afternoon sun in the west
    @test 0 < solar_position(DateTime(2025, 4, 7, 7), 50.8, 4.35).azimuth < 180
    @test 180 < solar_position(DateTime(2025, 4, 7, 16), 50.8, 4.35).azimuth < 360
    # refraction lifts the sun near the horizon
    low = DateTime(2025, 4, 7, 5, 30)
    @test solar_position(low, 50.8, 4.35; refraction=true).zenith < solar_position(low, 50.8, 4.35).zenith
    # the ACE spectra: position and time from their headers. The SZA the instrument
    # software recorded is on average 0.18 degrees larger, at most 0.46.
    for f in readdir(joinpath(ACE, "spectra"); join=true)[1:12:end]
        s = read_mfc_std(f)
        hsza = parse(Float64, split(only(filter(l -> startswith(l, "SZA ="), s.header)), "=")[2])
        @test solar_position(s).zenith ≈ hsza atol = 0.6
    end
    nopos = MiniDOAS.MFCSpectrum(Float64[], 1, 0f0, "x", (2024, 1, 1), (0, 0, 0), NaN, NaN, NaN, NaN, String[])
    @test_throws ArgumentError solar_position(nopos)
    @test geometric_amf(30) ≈ 2.0 && geometric_amf(90) ≈ 1.0
    @test sza_amf(60) ≈ 2.0
    @test geometric_vcd(1e16, 30) ≈ 1e16 && isnan(geometric_vcd(1e16, 90)) && isnan(geometric_vcd(1e16, -2))
end

@testset "fit-window map" begin
    m = DOASModel(P)
    s = MiniDOAS.load_spectrum(MiniDOAS.respec(P; spectrum=P.spectra[3]))[1]
    wm = fit_window_map(m, s; window="NO2", starts=327.0:1.0:329.0, ends=343.0:1.0:345.0, threads=false)
    @test wm.starts == [327.0, 328.0, 329.0] && wm.ends == [343.0, 344.0, 345.0]
    @test count(!isnothing, wm.results) == 9
    # every cell is the window fitted with that interval: the project's 328-345 nm is cell (2, 3)
    @test identical(wm.results[2, 3], analyze(P; spectrum=P.spectra[3])["NO2"])
    @test wm["NO2"][2, 3] == wm.results[2, 3].slcol["NO2"]
    @test wm["rel:NO2"][1, 1] ≈ wm.results[1, 1].slerr["NO2"] / abs(wm.results[1, 1].slcol["NO2"])
    @test wm["iterations"][3, 3] == wm.results[3, 3].niter
    # from a project file, with a minimum width
    wp = fit_window_map(P, P.spectra[3]; window=2, starts=[328.0, 338.0], ends=[340.0, 345.0], min_width=5.0,
                        threads=false)
    @test wp.results[2, 1] === nothing && wp.errors[2, 1] == "narrower than 5.0 nm"
    @test identical(wp.results[1, 2], wm.results[2, 3])
    @test isnan(wp["rms"][2, 1])
    # cancelling
    wc = fit_window_map(m, s; window=2, starts=327.0:1.0:329.0, ends=343.0:1.0:345.0, threads=false,
                        on_progress=(k, n) -> k < 2)
    @test count(==("cancelled"), wc.errors) == 7
    @test_throws ArgumentError fit_window_map(m, s; window="O4", starts=[1.0], ends=[2.0])
end

@testset "GUI helpers" begin
    e = Dict{String,Any}("SO2" => Dict{String,Any}("exclude" => ["Ring"], "poly_order" => 2, "lmin" => 304.0))
    p2 = MiniDOAS.apply_edits(P, e)
    w = p2.windows[1]
    @test [c.name for c in w.xs] == ["SO2", "O3"] && w.poly_order == 2 && w.lmin == 304.0
    @test all(g -> !("Ring" in g.symbols), w.shifts)
    @test p2.windows[2] === P.windows[2]
    saved = MiniDOAS.save_edited_project(SYN.project, joinpath(DIR, "edited.xml"), e)
    q = parse_config(saved)
    @test [c.name for c in q.windows[1].xs] == ["SO2", "O3"] && q.windows[1].poly_order == 2 && q.windows[1].lmin == 304.0
    @test [c.name for c in q.windows[2].xs] == ["NO2", "O3", "Ring"]
    r1 = analyze(p2; spectrum=P.spectra[1]); r2 = analyze(q; spectrum=P.spectra[1])
    @test r1.windows[1].slcol == r2.windows[1].slcol          # edits in memory == edits saved to file
    @test MiniDOAS.apply_edits(P, Dict("NO2" => Dict("enabled" => false))).windows |> length == 1
end

@testset "spectrum header information" begin
    # MFC STD key = value lines that QDOAS reads; they do not enter the fit
    lines = readlines(P.spectra[1])
    f = joinpath(DIR, "with_angles.STD")
    write(f, join(vcat(lines, ["ElevationAngle = 170", "AzimuthAngle = 20", "Latitude = 40.5"]), "\n") * "\n")
    r = analyze(P; spectrum=f)
    @test r.record.elevation == 10.0 && r.record.azimuth == 200.0     # > 100 degrees looks backwards
    @test r.record.latitude == 40.5 && r.record.longitude == 0.0
    @test r.windows[1].slcol == analyze(P; spectrum=P.spectra[1]).windows[1].slcol
    header, rows = results_table(r)
    @test header[7:9] == ["elevation", "azimuth", "latitude"] && rows[1][7] == 10.0
    @test "elevation" ∉ results_table(analyze(P; spectrum=P.spectra[2]))[1]
end

@testset "GUI server" begin
    G = Base.get_extension(MiniDOAS, :MiniDOASGUIExt)
    @test G !== nothing
    json(resp) = JSON3.read(String(resp.body), Dict{String,Any})
    @test G.settings_key(Dict("B" => Dict("x" => 1, "a" => [2, 1]), "A" => Dict())) ==
          G.settings_key(Dict("A" => Dict(), "B" => Dict("a" => [1, 2], "x" => 1.0)))
    @test G.settings_key(nothing) == G.settings_key(Dict()) == "{}"
    @test_throws ArgumentError MiniDOAS.apply_edits(P, Dict("SO2" => Dict("lmin" => 310.0, "lmax" => 310.5)))
    d = json(G.api_project(Dict{String,Any}("path" => SYN.project)))
    @test d["name"] == "Synthetic" && length(d["spectra"]) == 6
    e = Dict{String,Any}("SO2" => Dict{String,Any}("lmin" => 305.0))
    a = json(G.api_analyze(Dict{String,Any}("spectrum" => P.spectra[1], "edits" => e, "baseline" => true)))
    @test [w["name"] for w in a["windows"]] == ["SO2", "NO2"] && haskey(a, "baseline")
    @test a["windows"][1]["curves"]["fit_range"][1] >= 305.0
    @test a["baseline"]["SO2"]["terms"][1]["term"] == "SO2"
    @test a["record"]["scans"] == 200
    G.api_analyze(Dict{String,Any}("spectrum" => P.spectra[2], "edits" => e))
    t = json(G.api_results(Dict{String,Any}("edits" => e)))
    @test length(t["rows"]) == 2 && t["header"][1] == "file"
    @test isempty(json(G.api_results(Dict{String,Any}("edits" => nothing)))["rows"])
    csv = G.api_results_csv(Dict{String,Any}("edits" => e))
    @test count(==('\n'), String(csv.body)) == 3
    # fit-window map: every cell is the window fitted with that start and end
    G.api_map_start(Dict{String,Any}("spectrum" => P.spectra[3], "window" => "NO2", "edits" => nothing,
        "start_min" => 327.0, "start_max" => 329.0, "end_min" => 343.0, "end_max" => 345.0, "step" => 1.0))
    t0 = time()
    while G.GUI.jobs["map"].running && time() - t0 < 120
        sleep(0.05)
    end
    st = json(G.api_map_status(nothing))
    @test st["error"] == "" && st["done"] == st["total"] == 9
    mp = st["map"]
    @test mp["starts"] == [327.0, 328.0, 329.0] && mp["ends"] == [343.0, 344.0, 345.0]
    ref = analyze(P; spectrum=P.spectra[3]).windows[2]                 # the project's 328-345 nm
    @test mp["cells"][2][3]["scd"]["NO2"] ≈ ref.slcol["NO2"] rtol = 1e-5
    # batch job
    G.api_batch_start(Dict{String,Any}("spectra" => P.spectra, "edits" => nothing))
    t0 = time()
    while G.GUI.jobs["batch"].running && time() - t0 < 120
        sleep(0.05)
    end
    @test length(json(G.api_results(Dict{String,Any}("edits" => nothing)))["rows"]) == 6
    # NetCDF export of one spectrum
    nc = G.api_netcdf(Dict{String,Any}("spectrum" => P.spectra[1], "edits" => nothing))
    @test nc.status == 200 && length(nc.body) > 1000
    # saving refuses to overwrite the original, and asks before replacing another file
    @test_throws ArgumentError G.api_save(Dict{String,Any}("path" => SYN.project, "edits" => e))
    out = joinpath(DIR, "gui_saved.xml")
    @test json(G.api_save(Dict{String,Any}("path" => out, "edits" => e)))["saved"] == out
    @test haskey(json(G.api_save(Dict{String,Any}("path" => out, "edits" => e))), "exists")
    @test parse_config(out).windows[1].lmin == 305.0
    b = json(G.api_browse(Dict{String,Any}("path" => DIR)))
    @test "spectra" in b["dirs"] && !isempty(b["roots"])
    f = json(G.api_folder(Dict{String,Any}("dir" => DIR, "filter" => "S*.STD", "recursive" => true)))
    @test length(f["spectra"]) == 6
    # the server itself: the page, and requests from elsewhere are refused
    port = G.free_port("127.0.0.1", 18765)
    server = @async G.serve_gui(port=port, open_browser=false)
    t0 = time()
    page = nothing
    while time() - t0 < 60
        page = try HTTP.get("http://127.0.0.1:$port/"; retry=false, readtimeout=5) catch; nothing end
        page === nothing || break
        sleep(0.2)
    end
    @test page !== nothing && occursin("<title>MiniDOAS.jl</title>", String(page.body))
    st = HTTP.get("http://127.0.0.1:$port/api/state"; retry=false)
    @test JSON3.read(String(st.body))["version"] == MiniDOAS.VERSION_STRING
    forbidden = HTTP.get("http://127.0.0.1:$port/api/state"; headers=["Host" => "evil.example"], status_exception=false, retry=false)
    @test forbidden.status == 403
    notjson = HTTP.post("http://127.0.0.1:$port/api/project", ["Content-Type" => "text/plain"], "{}"; status_exception=false, retry=false)
    @test notjson.status == 415
end

@testset "extensions that are not loaded" begin
    @test_throws ErrorException MiniDOAS.extension(:NoSuchExt, "NoSuchPackage")
    @test occursin("Pkg.add([\"HTTP\", \"JSON3\"])", try MiniDOAS.extension(:NoSuchExt, "HTTP, JSON3") catch e; e.msg end)
end

# Optional: compare with QDOAS itself (set QDOAS_DOAS_CL=/path/to/doas_cl)
if haskey(ENV, "QDOAS_DOAS_CL") && isfile(ENV["QDOAS_DOAS_CL"])
    @testset "agreement with QDOAS doas_cl" begin
        xml = read(SYN.project, String)
        for f in P.spectra[1:2]
            nc = joinpath(DIR, "q_" * basename(f) * ".nc")
            s = replace(xml, r"<raw_spectra>.*?</raw_spectra>"s => "<raw_spectra><directory name=\"$(dirname(f))\" filters=\"$(basename(f))\" recursive=\"false\" /></raw_spectra>")
            s = replace(s, r"(<output\s+path=\")[^\"]*\"" => SubstitutionString("\\g<1>$nc\""))
            cfg = joinpath(DIR, "q.xml"); write(cfg, s)
            @test success(pipeline(`$(ENV["QDOAS_DOAS_CL"]) -c $cfg -a Synthetic`; stdout=devnull, stderr=devnull))
            r = analyze(P; spectrum=f)
            q = read_qdoas_netcdf(nc)
            for wr in r.windows, c in wr.columns
                haskey(q.windows[wr.name], "SlCol($c)") || continue
                v = q.windows[wr.name]["SlCol($c)"][1]
                σ = get(q.windows[wr.name], "SlErr($c)", [abs(v)])[1]
                @test abs(wr.slcol[c] - v) <= 1e-6 * σ
            end
        end
    end
end

end
