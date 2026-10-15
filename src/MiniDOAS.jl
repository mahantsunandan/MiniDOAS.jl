# MiniDOAS.jl - Differential Optical Absorption Spectroscopy in Julia.
# Copyright (c) 2026 Sunandan Mahant <sunandanmahant@outlook.com>
# The analysis core is derived from QDOAS, Copyright (C) 1994-2025 BIRA-IASB and S[&]T
# (BSD-3-Clause). See LICENSE and NOTICE.md.

"""
    MiniDOAS

Differential Optical Absorption Spectroscopy (DOAS) in Julia: read and write DOAS files,
prepare cross sections, fit spectra, and look at the results.

**Fitting.** The fit is QDOAS's optical-density fit (BIRA-IASB and S[&]T,
<https://uv-vis.aeronomie.be/software/QDOAS/>), ported routine by routine so that the
slant columns agree with QDOAS to floating-point precision. Use it in two ways:

  * with a QDOAS project file: `parse_config`, `analyze`, `analyze_spectra`
    (anything the port does not reproduce exactly raises `UnsupportedConfig`);
  * from Julia, without QDOAS: describe fit windows with `FitWindow`, `CrossSection` and
    `ShiftStretch`, set up a `DOASModel` with a wavelength grid and a reference, then
    `doas_fit` spectra.

Both take one spectrum (a file) and return a `DOASResult`, or several (a vector of files, a
folder, or a pattern such as `"spectra/*.STD"`) and return a `DOASBatch`:

    r = analyze("station.xml", "spectra/S0000001.STD")
    batch = doas_fit(model, "spectra/*.STD"; dark="dark.STD")
    no2, err = slant_column(batch, "NO2")

**Files.** MFC STD spectra (`read_mfc_std`, `write_mfc_std`), references, cross sections
and calibrations (`read_reference`, `read_cross_section`, `read_calibration`, `read_columns`
and their `write_*` counterparts), QDOAS project files (`parse_config`,
`write_qdoas_project`), results as CSV (`write_csv`) or Tables.jl-compatible columns
(`results_columns`), and QDOAS-layout NetCDF (`write_netcdf`, `read_qdoas_netcdf`; load
NCDatasets.jl).

**Tools.** Slit functions and convolution (`convolve`, `convolve_i0`), Ring spectra
(`ring_spectrum`), air/vacuum wavelengths, a solar-spectrum wavelength calibration
(`calibrate_wavelength`), smoothing and differential spectra, simulated spectra
(`simulate_spectrum`), the solar position and geometric air mass factors
(`solar_position`, `geometric_amf`, `geometric_vcd`), and fit-window maps
(`fit_window_map`).

**Extensions.** `serve_gui` (a local web GUI; load HTTP.jl and JSON3.jl), `plot_fit` and
`plot_window_map` (load a Makie backend such as CairoMakie.jl).
"""
module MiniDOAS

using Dates
using LinearAlgebra
using Printf
using Random

export UnsupportedConfig, EngineError,
       # types
       CrossSection, ShiftStretch, FitWindow, ProjectSpec, DOASModel, DOASResult, DOASBatch, WindowResult,
       MFCSpectrum, WindowMap,
       SlitFunction, GaussianSlit, BoxcarSlit, TriangularSlit, LorentzianSlit, TabulatedSlit,
       # QDOAS projects
       parse_config, project_spectra, analyze, analyze_config, analyze_spectra, run_config,
       write_qdoas_project,
       # fitting from Julia
       doas_fit, fit_window_map, slant_column,
       # files
       spectrum_files, read_mfc_std, write_mfc_std, measurement_time, correct_dark_offset,
       read_columns, write_columns, read_reference, write_reference, read_calibration,
       write_calibration, read_cross_section, write_cross_section,
       # results
       results_table, results_columns, write_csv, print_summary, write_netcdf, read_qdoas_netcdf,
       # tools
       convolve, convolve_i0, ring_spectrum, air_to_vacuum, vacuum_to_air, calibrate_wavelength,
       smooth, differential, simulate_spectrum,
       solar_position, geometric_amf, geometric_vcd, sza_amf,
       # extensions
       serve_gui, plot_fit, plot_window_map

"Version of this package (kept in step with Project.toml)."
const VERSION_STRING = "0.1.0"

const EPSILON = 1.0e-6
const MAX_REPEAT_CURFIT = 3
const CURFIT_MAX_ITER = 100

"An option of a QDOAS project that MiniDOAS.jl does not reproduce exactly (use QDOAS for it)."
struct UnsupportedConfig <: Exception
    msg::String
end
Base.showerror(io::IO, e::UnsupportedConfig) = print(io, "UnsupportedConfig: ", e.msg)

"A condition under which QDOAS itself would report an error for this analysis window."
struct EngineError <: Exception
    msg::String
end
Base.showerror(io::IO, e::EngineError) = print(io, "EngineError: ", e.msg)

unsupported(msg) = throw(UnsupportedConfig(msg))

include("numerics.jl")       # splines, pixel search, normalisation
include("types.jl")          # CrossSection, ShiftStretch, FitWindow, ProjectSpec
include("files.jl")          # MFC STD spectra, columns, references, cross sections
include("config.jl")         # QDOAS project files (.xml) and their spectra
include("fit.jl")            # the fit: model function and Levenberg-Marquardt (Curfit)
include("results.jl")        # per-window results in QDOAS's conventions
include("model.jl")          # DOASModel, doas_fit, analyze, batches
include("mapping.jl")        # fit-window maps
include("table.jl")          # CSV, result columns, printed summaries
include("qdoas_project.jl")  # writing and editing QDOAS project files
include("raman_tables.jl")   # rotational Raman lines of N2 and O2 (QDOAS ring.c)
include("tools.jl")          # slits, convolution, Ring, air/vacuum, calibration, filters
include("solar.jl")          # solar position, geometric air mass factors
include("extensions.jl")     # functions implemented by package extensions
include("cli.jl")            # command-line interface (bin/minidoas.jl)

end # module
