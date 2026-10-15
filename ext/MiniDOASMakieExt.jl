# This file is part of MiniDOAS.jl (Differential Optical Absorption Spectroscopy in Julia).
# Copyright (c) 2026 Sunandan Mahant <sunandanmahant@outlook.com>; BSD-3-Clause, see LICENSE.

# =============================================================================
# Figures with Makie: the fit of one window, and fit-window maps.
# =============================================================================

module MiniDOASMakieExt

using MiniDOAS
using MiniDOAS: DOASResult, WindowMap, fmt_date, fmt_time
using Makie
using Printf

sci(x) = isfinite(x) ? @sprintf("%.3e", x) : "-"

function plot_fit(r::DOASResult; window=1, size=(900, 900))
    wr = window isa Integer ? r.windows[window] : r[String(window)]
    cv = wr.curves
    cv === nothing && throw(ArgumentError("the result has no fitted curves: fit with curves=true"))
    absorbers = [c for c in cv.components if c.kind === :xs]
    fig = Figure(size=size)
    title = r.record === nothing ? (isempty(r.spectrum) ? wr.name : basename(r.spectrum)) :
            "$(r.record.name)  $(fmt_date(r.record.date)) $(fmt_time(r.record.time)) UTC"
    Label(fig[0, 1:2], "$title - window $(wr.name), RMS $(sci(wr.rms)), $(wr.niter) iterations";
          fontsize=15, tellwidth=false)

    a1 = Axis(fig[1, 1:2]; xlabel="wavelength (nm)", ylabel="intensity", title="spectrum")
    vspan!(a1, [cv.fit_range[1]], [cv.fit_range[2]]; color=(:teal, 0.15))
    lines!(a1, r.lambda, r.raw; color=:black, linewidth=1)

    a2 = Axis(fig[2, 1:2]; ylabel="optical depth", title="fit: ln(I₀/I) and model")
    scatter!(a2, cv.lambda, cv.od; color=(:black, 0.6), markersize=3)
    lines!(a2, cv.lambda, cv.model; color=:crimson, linewidth=1.5)

    a3 = Axis(fig[3, 1:2]; xlabel="wavelength (nm)", ylabel="residual", title="residual")
    hlines!(a3, [0.0]; color=:gray, linestyle=:dash)
    lines!(a3, cv.lambda, cv.residual; color=:black, linewidth=1)
    linkxaxes!(a2, a3)
    hidexdecorations!(a2; grid=false)

    for (k, c) in enumerate(absorbers)
        row, col = 4 + (k - 1) ÷ 2, 1 + (k - 1) % 2
        ax = Axis(fig[row, col]; title="$(c.name): $(sci(c.scd)) ± $(sci(c.err))", ylabel="optical depth")
        scatter!(ax, cv.lambda, c.value .+ cv.residual; color=(:black, 0.5), markersize=3)
        lines!(ax, cv.lambda, c.value; color=:crimson, linewidth=1.5)
        row == 4 + (length(absorbers) - 1) ÷ 2 && (ax.xlabel = "wavelength (nm)")
    end
    fig
end

function plot_window_map(m::WindowMap, quantity::AbstractString="rms"; size=(700, 560))
    z = m[quantity]
    what = MiniDOAS.describe_quantity(quantity)
    fig = Figure(size=size)
    ax = Axis(fig[1, 1]; xlabel="start of the fit window (nm)", ylabel="end of the fit window (nm)",
              title="window $(m.window): $what", aspect=DataAspect())
    hm = heatmap!(ax, m.starts, m.ends, z; colormap=:viridis, nan_color=:transparent)
    # the window's own interval, when it lies on the map
    a, b = m.range
    if first(m.starts) <= a <= last(m.starts) && first(m.ends) <= b <= last(m.ends)
        scatter!(ax, [a], [b]; marker=:circle, markersize=14, color=:transparent, strokecolor=:white, strokewidth=2)
    end
    Colorbar(fig[1, 2], hm; label=what)
    fig
end

end # module
