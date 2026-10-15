# Opens the MiniDOAS.jl web GUI on an example project, or on a project of your own.
#   julia --project=/path/to/MiniDOAS.jl -t auto examples/run_gui.jl            # an example
#   julia --project=/path/to/MiniDOAS.jl -t auto examples/run_gui.jl ace        # real MAX-DOAS spectra
#   julia --project=/path/to/MiniDOAS.jl -t auto examples/run_gui.jl synthetic  # synthetic spectra
#   julia --project=/path/to/MiniDOAS.jl -t auto examples/run_gui.jl my_project.xml
# The GUI needs HTTP.jl and JSON3.jl (and NCDatasets.jl for NetCDF export) in the
# environment. The real-data example is prepared on first use (it downloads public cross
# sections).
# Copyright (c) 2026 Sunandan Mahant <sunandanmahant@outlook.com>; BSD-3-Clause, see LICENSE.

using MiniDOAS, HTTP, JSON3
try
    @eval using NCDatasets        # optional: NetCDF export from the GUI
catch
end

const ACE = joinpath(@__DIR__, "ace_maxdoas", "work", "ace_maxdoas.xml")
const SYNTHETIC = joinpath(@__DIR__, "synthetic", "synthetic_project.xml")

function synthetic_project()
    isfile(SYNTHETIC) && return SYNTHETIC
    include(joinpath(@__DIR__, "synthetic_data.jl"))
    Base.invokelatest(() -> Main.make_synthetic(joinpath(@__DIR__, "synthetic"))).project
end

function ace_project()
    isfile(ACE) && return ACE
    include(joinpath(@__DIR__, "ace_maxdoas", "prepare.jl"))
    Base.invokelatest(() -> Main.main())
end

arg = isempty(ARGS) ? "" : ARGS[1]
project = arg == "ace" ? ace_project() :
          arg == "synthetic" ? synthetic_project() :
          !isempty(arg) ? arg :
          isfile(ACE) ? ACE : synthetic_project()
serve_gui(project=project)
