# DOAS without QDOAS: prepare cross sections, simulate spectra and fit them, all in Julia.
#   julia --project=/path/to/MiniDOAS.jl -t auto examples/fit_from_julia.jl
# Copyright (c) 2026 Sunandan Mahant <sunandanmahant@outlook.com>; BSD-3-Clause, see LICENSE.
#
# The "laboratory" data here are synthetic (smooth analytic shapes, NOT real cross
# sections); with real data you would read them with `read_cross_section` and
# `read_reference`, or `read_mfc_std` for spectra.

using MiniDOAS, Printf, Random

# --- a high-resolution "solar spectrum" and cross sections (0.01 nm) -----------------------
hires = collect(325.0:0.01:375.0)
rng = MersenneTwister(1)
lines = [(325 + 50rand(rng), 0.05 + 0.4rand(rng), 0.02 + 0.08rand(rng)) for _ in 1:400]
solar = [1e4 * prod(1 - d * exp(-((λ - c) / w)^2) for (c, d, w) in lines) for λ in hires]
no2_hr = CrossSection("NO2", hires, [2.5e-19 * (1 + 0.4sin(2π * (λ - 330) / 2.4) + 0.2sin(2π * λ / 0.9)) for λ in hires])
o4_hr = CrossSection("O4", hires, [8e-46 * exp(-((λ - 360.8) / 1.6)^2) + 2e-46 * exp(-((λ - 343.4) / 1.4)^2) for λ in hires])

# --- the instrument: 0.05 nm pixels, a 0.5 nm Gaussian slit -------------------------------
λ = collect(330.0:0.05:370.0)
slit = GaussianSlit(0.5)
I0 = convolve(hires, solar, λ, slit)                       # the Fraunhofer reference
no2 = convolve(no2_hr, λ, slit)                            # cross sections at the instrument's resolution
o4 = convolve(o4_hr, λ, slit)
ring = CrossSection("Ring", λ, ring_spectrum(hires, solar, λ, slit))

# --- simulated measurements with known slant columns --------------------------------------
# (the spectrometer drifted by 0.01 nm: a fit should shift the reference by +0.01 nm)
truth = [(NO2=2e16 * k, O4=4e43 * k, Ring=0.02k) for k in 1:8]
spectra = [simulate_spectrum(λ, I0, [no2, o4, ring], [t.NO2, t.O4, t.Ring]; shift=0.01,
                             polynomial=[log(0.8), -0.004], noise=:photon, rng=rng) for t in truth]

# --- the fit: one window; the reference and cross sections share a shift and stretch -------
# (Shifting the noisy measured spectrum instead would interpolate its noise, which smooths
# it and pulls the fitted shift towards half a pixel: shift the noise-free vectors.)
window = FitWindow("NO2", 338.0, 368.0, [no2, o4, ring]; polynomial=3,
                   shifts=ShiftStretch("Ref", "NO2", "O4", "Ring"; stretch=1))
model = DOASModel(λ, I0, window)
r = doas_fit(model, spectra[3]; curves=true)
print_summary(r)

batch = doas_fit(model, spectra)                           # all of them, in parallel
println("\n  NO2 fitted ± error        true      O4 fitted ± error        true")
for (b, t) in zip(batch, truth)
    w = b.result["NO2"]
    @printf("  %.3e ± %.1e  %.3e   %.3e ± %.1e  %.3e\n", w.slcol["NO2"], w.slerr["NO2"], t.NO2,
            w.slcol["O4"], w.slerr["O4"], t.O4)
end

# --- where should the window be? -----------------------------------------------------------
m = fit_window_map(model, spectra[3]; window="NO2", starts=334:1.0:344, ends=360:1.0:370)
println("\n", m)
@printf("NO2 over all windows: %.3e to %.3e\n", extrema(filter(isfinite, m["NO2"]))...)

# --- write the same analysis as a QDOAS project, to run it in QDOAS too ----------------------
dir = mktempdir()
xml = write_qdoas_project(joinpath(dir, "from_julia.xml"), model)
println("\nQDOAS project: ", xml)

# Tables: `results_columns(batch)` is a Tables.jl table, e.g. `DataFrame(results_columns(batch))`.
# Figures: `using CairoMakie; plot_fit(r)`.
