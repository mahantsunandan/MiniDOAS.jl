# This file is part of MiniDOAS.jl (Differential Optical Absorption Spectroscopy in Julia).
# Copyright (c) 2026 Sunandan Mahant <sunandanmahant@outlook.com>; BSD-3-Clause, see LICENSE.

# =============================================================================
# Functions implemented by package extensions (ext/), loaded with other packages:
#   NCDatasets          -> write_netcdf, read_qdoas_netcdf     (MiniDOASNetCDFExt)
#   HTTP and JSON3      -> serve_gui                           (MiniDOASGUIExt)
#   Makie (CairoMakie)  -> plot_fit, plot_window_map           (MiniDOASMakieExt)
# =============================================================================

"The loaded extension `name`, or an error saying which packages to load."
function extension(name::Symbol, packages::AbstractString)
    ext = Base.get_extension(@__MODULE__, name)
    ext === nothing && error("this needs $packages: run `using $packages` first " *
                             "(install with `import Pkg; Pkg.add([$(join(repr.(split(packages, ", ")), ", "))])`)")
    ext
end

"""
    write_netcdf(path, result::DOASResult; config_path="")

Writes `result` in the NetCDF layout of QDOAS's `doas_cl`: a group named after the
project's swath name (usually "QDOAS Results") with the record fields (name, date, time),
then one group per analysis window holding its variables in QDOAS's order (Chi, RMS, iter,
processing_error, residual_spectrum, SlCol/SlErr per cross section, shifts and stretches,
polynomial, offsets), with QDOAS's types and fill values. Needs NCDatasets.jl
(`using NCDatasets`).
"""
write_netcdf(path::AbstractString, result::DOASResult; config_path::AbstractString="") =
    Base.invokelatest(extension(:MiniDOASNetCDFExt, "NCDatasets").write_netcdf, path, result; config_path=config_path)

"""
    write_netcdf(folder, batch; config_path="") -> Vector{String}

One QDOAS-layout NetCDF file per analysed spectrum of a batch, in `folder` (created if
needed), named after the spectrum files (`S0000001.STD` → `S0000001.nc`; the name of its
folder is added where two spectra share a name). Returns the paths written.
"""
function write_netcdf(folder::AbstractString, batch::AbstractVector; config_path::AbstractString="")
    ok = [x for x in batch if x.result !== nothing]
    mkpath(folder)
    [write_netcdf(joinpath(folder, n), x.result; config_path=config_path)
     for (x, n) in zip(ok, netcdf_names([x.file for x in ok]))]
end

"""
    read_qdoas_netcdf(path) -> NamedTuple

Reads a QDOAS (or MiniDOAS.jl) NetCDF result file: `(group, attributes, fields, windows)` with
the results group's name and attributes, its record fields (`"Name"`, `"Date (DD-MM-YYYY)"`,
...) as arrays with one entry per record, and per analysis window a dictionary of its
variables (`windows["NO2"]["SlCol(NO2)"]`) and `"attributes"`. Fill values are replaced by
NaN (missing for integers and strings). Needs NCDatasets.jl.
"""
read_qdoas_netcdf(path::AbstractString; kw...) =
    Base.invokelatest(extension(:MiniDOASNetCDFExt, "NCDatasets").read_qdoas_netcdf, path; kw...)

"""
    serve_gui(; project="", port=8765, host="127.0.0.1", open_browser=true)

Starts the local web GUI and serves it until interrupted with Ctrl-C: open a QDOAS project,
step through and fit its spectra, drag the fit window, map the fit window, run batches,
and export results, fits and figures. If `port` is taken, the next free port is used.
Start Julia with several threads (`julia -t auto`) so that batches and maps run in
parallel. Needs HTTP.jl and JSON3.jl (`using HTTP, JSON3`).

The GUI can read any file this Julia process can read and writes the files you ask it to,
so by default it only listens on this computer (127.0.0.1). Pass `host="0.0.0.0"` only on
a network you trust.
"""
serve_gui(; kw...) = Base.invokelatest(extension(:MiniDOASGUIExt, "HTTP, JSON3").serve_gui; kw...)

"""
    plot_fit(result::DOASResult; window=1, size=(900, 900)) -> Makie.Figure

The classic DOAS fit figure for one analysis window (name or number): the measured
spectrum with the fit window, the optical depth with the model and residual, and one panel
per fitted absorber - its fitted optical depth (line) against the measured one (absorber
plus residual, dots), with its slant column. The result needs its curves:
`doas_fit(...; curves=true)` or `analyze(...; curves=true)`. Needs a Makie backend
(`using CairoMakie`).
"""
plot_fit(result::DOASResult; kw...) =
    Base.invokelatest(extension(:MiniDOASMakieExt, "Makie").plot_fit, result; kw...)

"""
    plot_window_map(map::WindowMap, quantity="rms"; size=(700, 560)) -> Makie.Figure

A heat map of a `fit_window_map`: `quantity` (as for `map[quantity]`, e.g. `"NO2"`,
`"err:NO2"`, `"rel:NO2"`, `"rms"`) against the start and end of the fitting interval.
Needs a Makie backend (`using CairoMakie`).
"""
plot_window_map(map::WindowMap, quantity::AbstractString="rms"; kw...) =
    Base.invokelatest(extension(:MiniDOASMakieExt, "Makie").plot_window_map, map, quantity; kw...)

"Whether QDOAS-layout NetCDF output is available (NCDatasets loaded)."
netcdf_available() = Base.get_extension(@__MODULE__, :MiniDOASNetCDFExt) !== nothing
