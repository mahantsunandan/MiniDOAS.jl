# MiniDOAS.jl quick start: make the synthetic data set, fit it, write the results.
#   julia --project=/path/to/MiniDOAS.jl -t auto examples/quickstart.jl
# Copyright (c) 2026 Sunandan Mahant <sunandanmahant@outlook.com>; BSD-3-Clause, see LICENSE.

using MiniDOAS, Printf
include(joinpath(@__DIR__, "synthetic_data.jl"))

# 1. A QDOAS project with 12 synthetic spectra (see synthetic_data.jl). With your own
#    data you would start from the .xml project you saved in the QDOAS GUI.
dir = joinpath(@__DIR__, "synthetic")
syn = make_synthetic(dir)
println("Project: ", syn.project)

# 2. Read the project once; it lists its spectra.
p = parse_config(syn.project)
println(length(p.windows), " analysis windows, ", length(p.spectra), " spectra")

# 3. One spectrum, with a printed summary.
r = analyze(p; spectrum=p.spectra[6])
print_summary(r)

# 4. All spectra (in parallel when Julia runs with several threads).
t = @elapsed batch = analyze_spectra(p)
@printf("\n%d spectra in %.3f s on %d thread(s)\n\n", length(batch), t, Threads.nthreads())

println("spectrum        SO2 fitted ± error          SO2 true     NO2 fitted ± error        NO2 true")
for (b, truth) in zip(batch, syn.truth)
    so2, no2 = b.result["SO2"], b.result["NO2"]
    @printf("%-14s %.3e ± %.1e   %.3e    %.3e ± %.1e   %.3e\n", basename(b.file),
            so2.slcol["SO2"], so2.slerr["SO2"], truth.SO2, no2.slcol["NO2"], no2.slerr["NO2"], truth.NO2)
end

# 5. Results: a CSV table (and, with NCDatasets loaded, one QDOAS-layout NetCDF per spectrum).
out = mkpath(joinpath(dir, "results"))
write_csv(joinpath(out, "results.csv"), batch)
println("\nResults written to ", out)
# using NCDatasets
# for b in batch
#     write_netcdf(joinpath(out, splitext(basename(b.file))[1] * ".nc"), b.result; config_path=syn.project)
# end
