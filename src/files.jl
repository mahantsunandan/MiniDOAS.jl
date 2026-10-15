# This file is part of MiniDOAS.jl (Differential Optical Absorption Spectroscopy in Julia).
# Copyright (c) 2026 Sunandan Mahant <sunandanmahant@outlook.com>
# Derived from QDOAS, Copyright (C) 1994-2025 BIRA-IASB and S[&]T (BSD-3-Clause).
# See LICENSE and NOTICE.md.

# =============================================================================
# MFC STD spectra (mfc-read.c MFC_ReadRecordStd)
# =============================================================================

mutable struct MfcState               # QDOAS keeps these in a global record struct
    nsomme::Int
    tint::Float64
    total_acq::Float64
end
MfcState() = MfcState(0, 0.0, 0.0)

"""
    MFCSpectrum

One spectrum in the MFC STD text format, as QDOAS reads it (`read_mfc_std`):

  * `counts`     the counts, one per pixel (summed over the scans)
  * `scans`      number of scans (`SCANS n` or `NumScans = n`)
  * `exposure`   exposure time in s (`ExposureTime = ms`, else `int_TIME s`; see the README
                 on QDOAS's reading of `INT_TIME`), 0 when not given
  * `name`       spectrum name (first 20 characters of its line)
  * `date`       (year, month, day)
  * `time`       (hour, minute, second) in the middle of the exposure
  * `elevation`, `azimuth`, `latitude`, `longitude`   in degrees, NaN when not given
  * `header`     every line after the counts, as in the file

`measurement_time(s)` gives the date and time as a `DateTime`.
"""
struct MFCSpectrum
    counts::Vector{Float64}
    scans::Int
    exposure::Float32
    name::String                       # pRecord->Nom (first 20 characters)
    date::NTuple{3,Int}                # (year, month, day)
    time::NTuple{3,Int}                # mid-exposure (hour, min, sec)
    # header information that does not enter the fit (NaN when the file has none)
    elevation::Float64                 # viewing elevation angle, degrees
    azimuth::Float64                   # viewing azimuth angle, degrees
    latitude::Float64
    longitude::Float64
    header::Vector{String}
end

function Base.show(io::IO, s::MFCSpectrum)
    print(io, "MFCSpectrum(\"", s.name, "\", ", length(s.counts), " pixels, ")
    print(io, @sprintf("%04d-%02d-%02d %02d:%02d:%02d", s.date..., s.time...), ", ", s.scans, " scans")
    isnan(s.elevation) || print(io, ", elevation ", s.elevation, "°")
    print(io, ")")
end

"""
    measurement_time(s::MFCSpectrum) -> DateTime

The middle of the exposure (UTC, as the instrument recorded it).
"""
measurement_time(s::MFCSpectrum) = DateTime(s.date..., s.time...)

"Date and mid-exposure time as MFC_ReadRecordStd derives them."
function header_datetime(dline, t1, t2, date_format)
    iday = imon = iyear = 0; sep = 0
    for ch in date_format
        if ch in ('Y', 'y'); iyear = sep
        elseif ch in ('M', 'm'); imon = sep
        elseif ch in ('D', 'd'); iday = sep
        else sep += 1
        end
    end
    nums = [parse(Int, x.match) for x in eachmatch(r"\d+", dline)]
    length(nums) >= 3 || (nums = [0, 0, 0])
    day, mon, year = nums[iday+1], nums[imon+1], nums[iyear+1]
    year < 30 ? (year += 2000) : year < 130 ? (year += 1900) : year < 1930 && (year += 100)
    hms(t) = (v = [parse(Int, x.match) for x in eachmatch(r"\d+", t)]; length(v) >= 3 ? v[1:3] : [0, 0, 0])
    a, b = hms(t1), hms(t2)
    n1 = a[1] * 3600 + a[2] * 60 + a[3]
    n2 = b[1] * 3600 + b[2] * 60 + b[3]
    n2 < n1 && (n2 += 86400)
    nsec = (n1 + n2) ÷ 2
    (year, mon, day), (nsec ÷ 3600, (nsec % 3600) ÷ 60, (nsec % 3600) % 60)
end

"""
MFC_ReadRecordStd: one MFC STD file, `n_wavel` pixels (0: as many as the file holds).
`st` carries the scans and exposure between files exactly as QDOAS's global state does.
"""
function read_mfc_record(path, n_wavel, st::MfcState; date_format="MM/DD/YYYY")
    lines = readlines(path)
    length(lines) >= 3 || throw(EngineError("short MFC file $path"))
    pixfin = parse(Int, split(strip(lines[3]))[1])
    n_wavel = n_wavel > 0 ? n_wavel : pixfin
    spe = zeros(n_wavel)
    length(lines) >= 3 + pixfin || throw(EngineError("MFC file $path: $pixfin pixels announced, $(length(lines) - 3) lines"))
    li = 3
    for i in 1:pixfin
        li += 1
        v = tryparse(Float64, split(strip(lines[li]))[1])
        i <= n_wavel && v !== nothing && (spe[i] = v)
    end
    header = String.(lines[li+1:end])
    # specname, spectroname, scan_dev, date, start, end, 2 floats (8 lines)
    hdr(k) = li + k <= length(lines) ? strip(lines[li+k]) : ""
    name = first(hdr(1), 20)
    date, tmid = header_datetime(hdr(4), hdr(5), hdr(6), date_format)
    li += 8
    rest = lines[min(li + 1, length(lines) + 1):end]
    # fscanf "SCANS %d" / "int_TIME %lf" are literal, case-sensitive matches that
    # stop the sequence at the first mismatch - mirror that.
    j = 1
    if j <= length(rest) && (m = match(r"^SCANS\s+(-?\d+)", rest[j])) !== nothing
        st.nsomme = parse(Int, m[1]); j += 1
        if j <= length(rest) && (m2 = match(r"^int_TIME\s+(\S+)", rest[j])) !== nothing
            st.total_acq = parse(Float64, m2[1])
        end
    end
    elev = azim = lat = lon = NaN
    num(v) = something(tryparse(Float64, v), NaN)
    for l in rest
        if (mm = match(r"^(LONGITUDE|LATITUDE)\s+(\S+)", l)) !== nothing
            mm[1] == "LONGITUDE" ? (lon = num(mm[2])) : (lat = num(mm[2]))
            continue
        end
        occursin('=', l) || continue
        mm = match(r"^\s*(\S+)\s*=\s*(\S+)", l)
        mm === nothing && continue
        k, v = lowercase(mm[1]), mm[2]
        if k == "exposuretime"
            st.tint = something(tryparse(Float64, v), 0.0) * 0.001
        elseif k == "numscans"
            st.nsomme = something(tryparse(Int, v), 0)
        elseif k == "elevationangle"
            elev = num(v)
        elseif k == "azimuthangle"
            azim = num(v)
        elseif k == "latitude"
            lat = num(v)
        elseif k == "longitude"
            lon = num(v)
        end
    end
    # mfc-read.c: an elevation above 100 degrees looks backwards
    if elev > 100
        elev = 180 - elev
        isnan(azim) || (azim = mod(azim + 180, 360))
    end
    if st.tint < 1e-3 && st.total_acq > 1e-3
        st.tint = st.total_acq
    end
    MFCSpectrum(spe, st.nsomme, Float32(st.tint), name, date, tmid, elev, azim, lat, lon, header)
end

"""
    read_mfc_std(path; npixels=0, date_format="MM/DD/YYYY") -> MFCSpectrum

Reads one spectrum in the MFC STD text format exactly as QDOAS does: line 3 holds the
number of pixels N, lines 4..N+3 the counts, then the spectrum name, spectrometer and
device names, the date (laid out as `date_format`), the start and end times of the
exposure, two numbers, and optionally `SCANS n` / `int_TIME s` and `key = value` lines
(`ExposureTime` in ms, `NumScans`, `ElevationAngle`, `AzimuthAngle`, `Latitude`,
`Longitude`). `npixels > 0` reads that many pixels (padding with zeros).
"""
read_mfc_std(path::AbstractString; npixels::Integer=0, date_format::AbstractString="MM/DD/YYYY") =
    read_mfc_record(path, npixels, MfcState(); date_format=date_format)

"""
Whether `path` holds an MFC STD spectrum: a `.STD` file, or a text file whose first line is
not a number and whose second and third lines are integers (the third the pixel count).
"""
function is_mfc_std(path::AbstractString)
    lowercase(splitext(path)[2]) == ".std" && return true
    lines = String[]
    open(path) do io
        while length(lines) < 3 && !eof(io)
            push!(lines, readline(io))
        end
    end
    length(lines) == 3 || return false
    f1 = split(strip(lines[1]))
    (isempty(f1) || tryparse(Float64, f1[1]) === nothing) &&
        tryparse(Int, strip(lines[2])) !== nothing && something(tryparse(Int, strip(lines[3])), 0) > 0
end

"""
The intensities of a spectrum file and its header information: an MFC STD file (with its
`MFCSpectrum`), or a text file of numbers - one column (intensities) or several
(wavelength, intensity, ...: the second column is used) - with `nothing`.
"""
function read_spectrum_file(path::AbstractString; date_format::AbstractString="MM/DD/YYYY")
    isfile(path) || throw(ArgumentError("no such file: $path"))
    if is_mfc_std(path)
        s = read_mfc_std(path; date_format=date_format)
        return s.counts, s
    end
    c = read_columns(path)
    (length(c) == 1 ? c[1] : c[2]), nothing
end

"""
    spectrum_files(path; filter="*", recursive=false) -> Vector{String}

The spectrum files that `path` names, sorted:

  * a file: that file;
  * a folder: its files whose names match `filter` (shell patterns such as `"*.STD"`, several
    separated by `;`), with those in its sub-folders when `recursive`;
  * a pattern such as `"spectra/*.STD"` or `"data/S00001??.STD"` (wildcards in the file
    name): the matching files.

Names match whatever their case; hidden files (starting with `.`) are left out.
"""
function spectrum_files(path::AbstractString; filter::AbstractString="*", recursive::Bool=false)
    path = expanduser(path)
    patterns(s) = (r = [glob_regex(x) for x in split(s, r"[;,]+"; keepempty=false)]; isempty(r) ? [glob_regex("*")] : r)
    function scan(dir, rx)
        found = String[]
        if recursive
            for (root, _, fs) in walkdir(dir), f in fs
                startswith(f, ".") || !any(r -> occursin(r, f), rx) || push!(found, joinpath(root, f))
            end
        else
            for f in readdir(dir)
                p = joinpath(dir, f)
                startswith(f, ".") || !isfile(p) || !any(r -> occursin(r, f), rx) || push!(found, p)
            end
        end
        sort!(found)
    end
    if occursin(r"[*?]", basename(path))
        dir = dirname(path)
        dir = isempty(dir) ? "." : dir
        isdir(dir) || throw(ArgumentError("no such folder: $dir"))
        return scan(dir, patterns(basename(path)))
    end
    isfile(path) && return [String(path)]
    isdir(path) && return scan(path, patterns(filter))
    throw(ArgumentError("no such file or folder: $path"))
end

"Whether `path` names several spectra (a folder, or a pattern with `*` or `?`) rather than one file."
names_many(path::AbstractString) = occursin(r"[*?]", basename(path)) || isdir(path)

"Numbers as text: integers without decimals, anything else with all the digits needed."
fmtnum(x::Real) = (isinteger(x) && abs(x) < 1e15) ? string(Int64(x)) : @sprintf("%.17g", x)

"""
    write_mfc_std(path, s::MFCSpectrum)
    write_mfc_std(path, counts; name="", spectrometer="", device="", start::DateTime,
                  stop=start, scans=1, exposure=0.0, elevation=NaN, azimuth=NaN,
                  latitude=NaN, longitude=NaN, extra=Pair{String,Any}[])

Writes a spectrum in the MFC STD format. An `MFCSpectrum` is written with its header lines
as they were read. Otherwise the header is made from the keywords: the date and the start
and end of the exposure (`start`, `stop`, as `DateTime`s in UTC), `SCANS` and `int_TIME`
(the exposure in s, which QDOAS uses for the dark-current correction), the angles and
position as `key = value` lines when given, and any `extra` `key => value` pairs.
"""
function write_mfc_std(path::AbstractString, s::MFCSpectrum; counts::AbstractVector{<:Real}=s.counts)
    open(path, "w") do io
        println(io, "GDBGMNUP"); println(io, 1); println(io, length(counts))
        foreach(c -> println(io, fmtnum(c)), counts)
        foreach(l -> println(io, l), s.header)
    end
    path
end

function write_mfc_std(path::AbstractString, counts::AbstractVector{<:Real}; name::AbstractString="",
                       spectrometer::AbstractString="", device::AbstractString="", start::DateTime,
                       stop::DateTime=start, scans::Integer=1, exposure::Real=0.0,
                       elevation::Real=NaN, azimuth::Real=NaN, latitude::Real=NaN, longitude::Real=NaN,
                       extra=Pair{String,Any}[])
    stop >= start || throw(ArgumentError("the exposure ends before it starts"))
    open(path, "w") do io
        println(io, "GDBGMNUP"); println(io, 1); println(io, length(counts))
        foreach(c -> println(io, fmtnum(c)), counts)
        println(io, first(name, 20)); println(io, spectrometer); println(io, device)
        println(io, Dates.format(start, "mm/dd/yyyy"))
        println(io, Dates.format(start, "HH:MM:SS")); println(io, Dates.format(stop, "HH:MM:SS"))
        println(io, 0); println(io, 0)
        println(io, "SCANS ", scans)
        println(io, "int_TIME ", fmtnum(exposure))
        isnan(elevation) || println(io, "ElevationAngle = ", fmtnum(elevation))
        isnan(azimuth) || println(io, "AzimuthAngle = ", fmtnum(azimuth))
        isnan(latitude) || println(io, "Latitude = ", fmtnum(latitude))
        isnan(longitude) || println(io, "Longitude = ", fmtnum(longitude))
        for (k, v) in extra
            println(io, k, " = ", v isa Real ? fmtnum(v) : v)
        end
    end
    path
end

"Removes the offset and dark current from the counts `spe` in place (MFC_LoadOffset/Dark order)."
function apply_dark_offset!(spe, scans, int_time::Float32, off, drk)
    n = length(spe)
    if off !== nothing && off.scans > 0
        for i in 1:n
            spe[i] -= Float64(off.counts[i]) * scans / off.scans
        end
    end
    if drk !== nothing && drk.exposure != 0f0
        den = Float64(drk.exposure * Float32(drk.scans))
        for i in 1:n
            spe[i] -= Float64(scans) * drk.counts[i] * Float64(int_time) / den
        end
    end
    spe
end

"The dark spectrum with the offset removed (MFC_LoadDark)."
function dark_minus_offset(d::MFCSpectrum, off)
    v = copy(d.counts)
    if off !== nothing && off.scans > 0
        for i in eachindex(v)
            v[i] -= Float64(off.counts[i]) * d.scans / off.scans
        end
    end
    MFCSpectrum(v, d.scans, d.exposure, d.name, d.date, d.time, d.elevation, d.azimuth, d.latitude, d.longitude, d.header)
end

"""
    correct_dark_offset(s::MFCSpectrum; dark=nothing, offset=nothing) -> Vector{Float64}

The counts of `s` corrected as QDOAS corrects MFC spectra: the `offset` spectrum (an
`MFCSpectrum`) is subtracted scaled by the number of scans, and the `dark` spectrum (with
the offset removed from it) scaled by scans × exposure time:

    I = R - O·n_R/n_O - D·(n_R t_R)/(n_D t_D)

A dark spectrum with no exposure time is not used (QDOAS's behaviour).
"""
function correct_dark_offset(s::MFCSpectrum; dark::Union{Nothing,MFCSpectrum}=nothing,
                             offset::Union{Nothing,MFCSpectrum}=nothing)
    n = length(s.counts)
    for x in (dark, offset)
        x === nothing || length(x.counts) == n || throw(DimensionMismatch("dark/offset and spectrum sizes differ"))
    end
    drk = dark === nothing ? nothing : dark_minus_offset(dark, offset)
    apply_dark_offset!(copy(s.counts), s.scans, s.exposure, offset, drk)
end

"Spectrum after QDOAS's offset and dark corrections (engine.c / mfc-read.c order)."
function load_spectrum(p::ProjectSpec)
    n = p.n_det
    st = MfcState()
    off = isempty(p.offset) ? nothing : read_mfc_record(p.offset, n, st)       # MFC_LoadOffset: raw read
    drk = isempty(p.dark) ? nothing : dark_minus_offset(read_mfc_record(p.dark, n, st), off)
    s = read_mfc_record(p.spectrum, n, st; date_format=p.date_format)
    spe = apply_dark_offset!(copy(s.counts), s.scans, s.exposure, off, drk)
    p.revert && reverse!(spe)
    spe, s
end

# =============================================================================
# References, cross sections, calibrations: columns of numbers
# =============================================================================

"AnalyseLoadVector for a two-column reference file."
function load_two_columns(path, n)
    lam = zeros(n); v = zeros(n)
    i = 0
    for l in eachline(path)
        (occursin(';', l) || occursin('*', l) || occursin('#', l)) && continue
        f = split(strip(l))
        length(f) >= 2 || continue
        i += 1
        i > n && break
        lam[i] = parse(Float64, f[1]); v[i] = parse(Float64, f[2])
    end
    i < n && throw(EngineError("reference $path has $i lines, expected $n"))
    lam, v
end

# Parsed cross sections and references, shared across calls (a batch analyses many
# spectra with the same windows). Keyed on path + size + mtime so an edited file is
# re-read.
const FILE_CACHE = Dict{Any,Any}()
const FILE_CACHE_LOCK = ReentrantLock()
const FILE_CACHE_MAX = 4000

function cached(f, key)
    st = stat(key[1])
    k = (key..., st.size, st.mtime)
    lock(FILE_CACHE_LOCK) do
        v = get(FILE_CACHE, k, nothing)
        v === nothing || return v
        length(FILE_CACHE) >= FILE_CACHE_MAX && empty!(FILE_CACHE)
        v = f()
        FILE_CACHE[k] = v
        v
    end
end

"MATRIX_Load (no range cut): wavelengths and values of a cross-section file."
function load_xs(path)
    rows = Vector{Vector{Float64}}()
    for l in eachline(path)
        s = strip(l)
        (isempty(s) || s[1] in (';', '#', '*')) && continue
        f = split(s)
        push!(rows, parse.(Float64, f))
    end
    isempty(rows) && throw(EngineError("empty cross section $path"))
    nc = length(rows[1])
    all(length(r) == nc for r in rows) || throw(EngineError("ragged cross section $path"))
    nc >= 2 || throw(EngineError("cross section $path has one column"))
    lam = [r[1] for r in rows]
    val = [r[2] for r in rows]
    if length(lam) > 1 && lam[1] > lam[2]
        reverse!(lam); reverse!(val)
    end
    lam, val, nc
end

const COMMENT_CHARS = (';', '#', '*', '%')

"""
    read_columns(path) -> Vector{Vector{Float64}}

The numeric columns of a text file. Blank lines, comment lines (starting with `;`, `#`,
`*` or `%`) and lines that do not start with a number (headers) are skipped; separators
are spaces, tabs or commas. Every data line must have the same number of columns.
"""
function read_columns(path::AbstractString)
    cols = Vector{Vector{Float64}}()
    for (ln, l) in enumerate(eachline(path))
        s = strip(l)
        (isempty(s) || s[1] in COMMENT_CHARS) && continue
        f = split(s, r"[\s,]+"; keepempty=false)
        v = tryparse.(Float64, f)
        v[1] === nothing && continue
        any(isnothing, v) && throw(ArgumentError("$path, line $ln: not a number"))
        if isempty(cols)
            append!(cols, [Float64[] for _ in v])
        elseif length(v) != length(cols)
            throw(ArgumentError("$path, line $ln: $(length(v)) columns, expected $(length(cols))"))
        end
        for (c, x) in zip(cols, v)
            push!(c, x)
        end
    end
    isempty(cols) && throw(ArgumentError("$path: no numbers"))
    cols
end

"""
    write_columns(path, columns...; header="")

Writes vectors of equal length as columns separated by tabs, every number with the
digits needed to read it back exactly. `header` lines (without the comment character)
are written first as `;` comments.
"""
function write_columns(path::AbstractString, columns::AbstractVector{<:Real}...; header::AbstractString="")
    n = length(first(columns))
    all(c -> length(c) == n, columns) || throw(DimensionMismatch("columns of different lengths"))
    open(path, "w") do io
        isempty(header) || foreach(l -> println(io, "; ", l), split(chomp(header), '\n'))
        for i in 1:n
            println(io, join((fmtnum(c[i]) for c in columns), '\t'))
        end
    end
    path
end

"""
    read_reference(path) -> (lambda, intensity)

A reference spectrum (or any spectrum) stored as two columns, wavelength and intensity.
"""
function read_reference(path::AbstractString)
    c = read_columns(path)
    length(c) >= 2 || throw(ArgumentError("$path: expected two columns"))
    c[1], c[2]
end

"""
    write_reference(path, lambda, intensity; header="")

Writes a spectrum as two columns that QDOAS reads as a reference file.
"""
write_reference(path::AbstractString, lambda, intensity; header::AbstractString="") =
    write_columns(path, lambda, intensity; header=header)

"""
    read_calibration(path) -> Vector{Float64}

A wavelength calibration: the first column, one wavelength (nm) per detector pixel.
"""
read_calibration(path::AbstractString) = first(read_columns(path))

"""
    write_calibration(path, lambda)

Writes a wavelength calibration, one wavelength per line (a QDOAS `.clb` file).
"""
write_calibration(path::AbstractString, lambda::AbstractVector{<:Real}) = write_columns(path, lambda)

"""
    read_cross_section(path; name=<file name>, action=:interp) -> CrossSection

Reads a cross section as QDOAS does (see `CrossSection`).
"""
read_cross_section(path::AbstractString; name::AbstractString=splitext(basename(path))[1], action::Symbol=:interp) =
    CrossSection(name, path; action=action)

"""
    write_cross_section(path, xs::CrossSection; header="")

Writes a cross section as two columns, wavelength and value.
"""
function write_cross_section(path::AbstractString, xs::CrossSection; header::AbstractString="")
    isloaded(xs) || throw(ArgumentError("cross section $(xs.name) is not loaded"))
    write_columns(path, xs.lambda, xs.sigma; header=header)
end
