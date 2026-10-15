# Tests of the Makie extension (figures). CairoMakie is large, so it is not a test
# dependency; run these in an environment that has MiniDOAS and CairoMakie:
#   julia --project=<env with MiniDOAS and CairoMakie> test/makie.jl
# Copyright (c) 2026 Sunandan Mahant <sunandanmahant@outlook.com>; BSD-3-Clause, see LICENSE.

using Test, MiniDOAS, CairoMakie
include(joinpath(@__DIR__, "..", "examples", "synthetic_data.jl"))

@testset "MiniDOAS.jl figures (Makie)" begin
    @test Base.get_extension(MiniDOAS, :MiniDOASMakieExt) !== nothing
    syn = make_synthetic(mktempdir(); nspec=3)
    p = parse_config(syn.project)
    r = analyze(p; spectrum=p.spectra[2], curves=true)
    for w in ("SO2", 2)
        fig = plot_fit(r; window=w)
        @test fig isa Figure
        f = tempname() * ".png"
        save(f, fig)
        @test filesize(f) > 10_000
    end
    @test_throws ArgumentError plot_fit(analyze(p; spectrum=p.spectra[2]))      # no curves
    m = fit_window_map(p, p.spectra[2]; window="NO2", starts=326.0:1.0:330.0, ends=343.0:1.0:347.0)
    @test m.range == (328.0, 345.0)
    for q in ("NO2", "err:NO2", "rel:NO2", "rms", "chi2", "iterations")
        @test plot_window_map(m, q) isa Figure
    end
    @test MiniDOAS.describe_quantity("rel:NO2") == "relative error of NO2"
end
