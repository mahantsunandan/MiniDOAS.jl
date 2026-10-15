# Makes the README figures (docs/img) with the Makie extension: a fit-window map of the
# synthetic example, and, once examples/ace_maxdoas/prepare.jl has been run, a fit and the O4
# time series of the real spectra.
#   julia --project=<an environment with MiniDOAS and CairoMakie> docs/make_figures.jl
# Copyright (c) 2026 Sunandan Mahant <sunandanmahant@outlook.com>; BSD-3-Clause, see LICENSE.

using MiniDOAS, CairoMakie
include(joinpath(@__DIR__, "..", "examples", "synthetic_data.jl"))

syn = make_synthetic(mktempdir(); nspec=6)
p = parse_config(syn.project)
out = mkpath(joinpath(@__DIR__, "img"))

m = fit_window_map(p, p.spectra[5]; window="NO2", starts=322.0:0.5:332.0, ends=338.0:0.5:348.0)
save(joinpath(out, "window_map.png"), plot_window_map(m, "NO2"); px_per_unit=1.5)

# real spectra, once examples/ace_maxdoas/prepare.jl has been run
ace = joinpath(@__DIR__, "..", "examples", "ace_maxdoas", "work", "ace_maxdoas.xml")
if isfile(ace)
    q = parse_config(ace)
    f = only(filter(s -> endswith(s, "20170407_180207_e+15.STD"), q.spectra))
    save(joinpath(out, "plot_fit_ace.png"), plot_fit(analyze(q; spectrum=f, curves=true); window="NO2", size=(900, 1100));
         px_per_unit=1.5)
    batch = analyze_spectra(q)
    fig = Figure(size=(900, 420))
    ax = Axis(fig[1, 1]; xlabel="time (UTC hours)", ylabel="O₄ dSCD (molecules² cm⁻⁵)",
              title="O₄ differential slant columns, 7 April 2017, by viewing elevation")
    elevs = sort(unique(b.result.record.elevation for b in batch))
    for (k, e) in enumerate(elevs)
        rs = [b.result for b in batch if b.result.record.elevation == e]
        h = [r.record.time[1] + r.record.time[2] / 60 + r.record.time[3] / 3600 for r in rs]
        scatterlines!(ax, h, [r["O4"].slcol["O4"] for r in rs]; label="$(Int(e))°",
                      color=k, colormap=:viridis, colorrange=(1, length(elevs)), markersize=6)
    end
    Legend(fig[1, 2], ax, "elevation")
    save(joinpath(out, "ace_o4.png"), fig; px_per_unit=1.5)
end
println("figures written to ", out)
