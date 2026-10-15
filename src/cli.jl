# This file is part of MiniDOAS.jl (Differential Optical Absorption Spectroscopy in Julia).
# Copyright (c) 2026 Sunandan Mahant <sunandanmahant@outlook.com>
# Derived from QDOAS, Copyright (C) 1994-2025 BIRA-IASB and S[&]T (BSD-3-Clause).
# See LICENSE and NOTICE.md.

# =============================================================================
# Command-line interface (bin/minidoas.jl), modelled on QDOAS's doas_cl
# =============================================================================

const CLI_USAGE = """
minidoas - DOAS analysis of QDOAS projects in Julia (MiniDOAS.jl $VERSION_STRING)

Usage:
  julia --project=<MiniDOAS.jl> -t auto <MiniDOAS.jl>/bin/minidoas.jl -c PROJECT.xml [options]
  julia --project=<MiniDOAS.jl> -t auto <MiniDOAS.jl>/bin/minidoas.jl --gui [PROJECT.xml] [--port N]

Options:
  -c, --config FILE    QDOAS project file (.xml), as saved by the QDOAS GUI      (required)
  -a, --project NAME   project name inside the file; checked, as with doas_cl -a
  -f, --file SPECTRUM  spectrum to analyse; repeat for several, wildcards allowed
                       (default: every spectrum the project's raw_spectra lists)
  -o, --output DIR     output directory (default: ./minidoas_output)
      --csv FILE       summary table (default: DIR/results.csv)
      --no-netcdf      write only the summary table, no per-spectrum NetCDF files
  -q, --quiet          no progress messages
      --gui [FILE]     open the web GUI instead (on the project FILE, if given)
      --port N         port of the web GUI (default 8765, or the next free one)
      --no-browser     with --gui: do not open a browser window
  -h, --help           this text
  -V, --version        version

Writes one NetCDF file per spectrum (DIR/<spectrum name>.nc, QDOAS layout) and a CSV
table with one row per spectrum. Run Julia with several threads (-t auto) to fit
spectra in parallel. Exit status: 0 all spectra analysed, 1 some failed, 2 usage or
project error.
"""

"Expand `*` / `?` in a file argument (Windows shells leave them to the program)."
function expand_file_arg(a::AbstractString)
    occursin(r"[*?]", a) || return [String(a)]
    dir, pat = dirname(a), basename(a)
    dir = isempty(dir) ? "." : dir
    rx = glob_regex(pat)
    isdir(dir) ? sort([joinpath(dir, f) for f in readdir(dir) if occursin(rx, f) && isfile(joinpath(dir, f))]) : String[]
end

"NetCDF file name per spectrum: the spectrum's name, disambiguated by its folder if needed."
function netcdf_names(files)
    base = [splitext(basename(f))[1] for f in files]
    counts = Dict{String,Int}()
    foreach(b -> counts[b] = get(counts, b, 0) + 1, base)
    [counts[b] > 1 ? string(basename(dirname(f)), "_", b) : b for (f, b) in zip(files, base)] .* ".nc"
end

"""
    main(args) -> exit code

Entry point of `bin/minidoas.jl`; see `MiniDOAS.CLI_USAGE`. Writing NetCDF needs NCDatasets
loaded and `--gui` needs HTTP and JSON3 (bin/minidoas.jl loads them).
"""
function main(args::Vector{String}=ARGS)
    cfg = ""; proj = ""; files = String[]; outdir = "minidoas_output"; csv = ""
    netcdf = true; quiet = false; gui = false; port = 8765; browser = true
    i = 1
    need(k) = i < length(args) ? args[i+1] : (println(stderr, "missing value for $k"); nothing)
    while i <= length(args)
        a = args[i]
        if a in ("-h", "--help")
            print(CLI_USAGE); return 0
        elseif a in ("-V", "--version")
            println("MiniDOAS.jl ", VERSION_STRING); return 0
        elseif a in ("-c", "--config", "-a", "--project", "-f", "--file", "-o", "--output", "--csv")
            v = need(a); v === nothing && return 2
            if a in ("-c", "--config"); cfg = v
            elseif a in ("-a", "--project"); proj = v
            elseif a in ("-f", "--file"); append!(files, expand_file_arg(v))
            elseif a in ("-o", "--output"); outdir = v
            else csv = v
            end
            i += 2; continue
        elseif a == "--gui"
            gui = true
            if i < length(args) && !startswith(args[i+1], "-")
                cfg = args[i+1]; i += 1
            end
        elseif a == "--port"
            v = need(a); v === nothing && return 2
            port = something(tryparse(Int, v), 0)
            port > 0 || (println(stderr, "--port needs a number"); return 2)
            i += 2; continue
        elseif a == "--no-browser"
            browser = false
        elseif a == "--no-netcdf"
            netcdf = false
        elseif a in ("-q", "--quiet")
            quiet = true
        else
            println(stderr, "unknown argument: $a\n"); print(stderr, CLI_USAGE); return 2
        end
        i += 1
    end
    if gui
        serve_gui(project=cfg, port=port, open_browser=browser)
        return 0
    end
    isempty(cfg) && (print(stderr, CLI_USAGE); return 2)
    if netcdf && !netcdf_available()
        println(stderr, "NetCDF output needs NCDatasets: load it (`using NCDatasets`) or pass --no-netcdf"); return 2
    end

    p = try
        parse_config(cfg)
    catch e
        println(stderr, "cannot use $cfg: ", sprint(showerror, e)); return 2
    end
    if !isempty(proj) && proj != p.project_name
        println(stderr, "project \"$proj\" not found in $cfg (it holds \"$(p.project_name)\")"); return 2
    end
    isempty(files) && (files = p.spectra)
    isempty(files) && (println(stderr, "no spectra: the project lists none and no -f was given"); return 2)

    mkpath(outdir)
    quiet || @printf("%s: %d spectra, %d analysis windows, %d threads\n", p.project_name, length(files),
                     length(p.windows), Threads.nthreads())
    t0 = time()
    report(k, n, f) = quiet || @printf("  [%d/%d] %s\n", k, n, basename(f))
    batch = analyze_spectra(p, files; on_progress=report)

    nfail = 0
    ncnames = netcdf_names(files)
    for (b, ncname) in zip(batch, ncnames)
        if b.result === nothing
            nfail += 1
            quiet || println("  failed: ", basename(b.file), ": ", b.error)
        elseif netcdf
            write_netcdf(joinpath(outdir, ncname), b.result; config_path=abspath(cfg))
        end
    end
    csvpath = isempty(csv) ? joinpath(outdir, "results.csv") : csv
    write_csv(csvpath, batch)
    quiet || @printf("%d analysed, %d failed in %.2f s; results in %s\n", length(files) - nfail, nfail,
                     time() - t0, abspath(outdir))
    nfail == 0 ? 0 : 1
end
