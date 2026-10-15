#!/usr/bin/env julia
# minidoas - command-line DOAS analysis of QDOAS projects with MiniDOAS.jl.
# Copyright (c) 2026 Sunandan Mahant <sunandanmahant@outlook.com>; BSD-3-Clause, see LICENSE.
#
# Usage (from anywhere; -t auto fits spectra in parallel):
#   julia --project=/path/to/MiniDOAS.jl -t auto /path/to/MiniDOAS.jl/bin/minidoas.jl -c project.xml
# Run with -h for all options. NetCDF output needs NCDatasets and --gui needs HTTP and
# JSON3 in the environment; they are loaded only when they are used.

using MiniDOAS
if !("--no-netcdf" in ARGS || "--gui" in ARGS || any(a -> a in ("-h", "--help", "-V", "--version"), ARGS))
    using NCDatasets
end
if "--gui" in ARGS
    using HTTP, JSON3
    try
        @eval using NCDatasets     # optional: NetCDF export from the GUI
    catch
    end
end
exit(Base.invokelatest(MiniDOAS.main, ARGS))
