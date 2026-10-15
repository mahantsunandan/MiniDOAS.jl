# This file is part of MiniDOAS.jl (Differential Optical Absorption Spectroscopy in Julia).
# Copyright (c) 2026 Sunandan Mahant <sunandanmahant@outlook.com>; BSD-3-Clause, see LICENSE.

# =============================================================================
# Writing QDOAS project files (.xml, in the layout the QDOAS GUI saves), and editing them.
# Only settings that MiniDOAS.jl reproduces are written; QDOAS itself opens and runs the files.
# =============================================================================

xmlesc(s) = replace(String(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;")
xmlnum(x::Real) = @sprintf("%.6g", x)
xmlbool(b::Bool) = b ? "true" : "false"
polystr(order::Integer) = order < 0 ? "none" : string(order)

function cross_section_xml(c::CrossSection, file::AbstractString)
    """        <cross_section sym="$(xmlesc(c.name))" ortho="None" cstype="$(c.action === :none ? "none" : "interp")" amftype="none" fit="true" filter="true" cstrncc="false" ccfit="true" icc="0.000" dcc="0.001" ccio="0.000" csfile="$(xmlesc(file))" amffile="" />\n"""
end

function shift_stretch_xml(g::ShiftStretch)
    st = g.st_type == 0 ? "none" : g.st_type == 1 ? "1st" : "2nd"
    """        <shift_stretch shfit="$(xmlbool(g.sh_fit))" stfit="$st" scfit="none" shstr="$(xmlbool(g.sh_store))" ststr="$(xmlbool(g.st_store))" scstr="false" errstr="$(xmlbool(g.err_store))" shini="$(xmlnum(g.sh_init))" stini="$(xmlnum(g.st_init))" stini2="$(xmlnum(g.st_init2))" scini="0.000" scini2="0.000" shdel="$(xmlnum(g.sh_delta))" stdel="$(xmlnum(g.st_delta))" stdel2="$(xmlnum(g.st_delta2))" scdel="0.0010" scdel2="0.0010" shmin="$(xmlnum(g.sh_min))" shmax="$(xmlnum(g.sh_max))" >\n""" *
    join("          <symbol name=\"$(xmlesc(s))\" />\n" for s in g.symbols) * "        </shift_stretch>\n"
end

"One `<analysis_window>` element: `reference` is its reference file, `files[name]` the file of each cross section."
function window_xml(w::FitWindow, reference::AbstractString, files::AbstractDict)
    cross = join(cross_section_xml(c, files[c.name]) for c in w.xs)
    outs = join(begin
        scol, serr, sfact = get(w.outputs, c.name, (true, true, 1.0))
        """        <output sym="$(xmlesc(c.name))" amf="false" scol="$(xmlbool(scol))" serr="$(xmlbool(serr))" sfact="$(xmlnum(sfact))" rescol="0.000000e+000" vcol="false" verr="false" vfact="1.000" />\n"""
    end for c in w.xs)
    gaps = join("        <gap min=\"$(xmlnum(a))\" max=\"$(xmlnum(b))\" />\n" for (a, b) in w.gaps)
    l0 = abs(w.lambda0 - 0.5 * (w.lmin + w.lmax)) < 1e-9 ? "" : " lambda0=\"$(@sprintf("%.17g", w.lambda0))\""
    nl = join(begin
        fit, init, delt = w.off[k]
        fstr, estr = w.off_store[k]
        "o$(k-1)fit=\"$(xmlbool(fit))\" o$(k-1)init=\"$(xmlnum(init))\" o$(k-1)delt=\"$(xmlnum(delt))\" o$(k-1)fstr=\"$(xmlbool(fstr))\" o$(k-1)estr=\"$(xmlbool(estr))\""
    end for k in 1:3)
    """
    <analysis_window name="$(xmlesc(w.name))" disable="false" kurucz="none" refsel="file" min="$(@sprintf("%.17g", w.lmin))" max="$(@sprintf("%.17g", w.lmax))"$l0 resol_fwhm="0.500" >
      <display spectrum="true" poly="true" fits="true" residual="true" predef="true" ratio="true" />
      <files saveresiduals="$(xmlbool(w.save_residuals))" refone="$(xmlesc(reference))" reftwo="" residual="" szacenter="0.000" szadelta="0.000" scanmode="after" minlon="0.000" maxlon="0.000" minlat="0.000" maxlat="0.000" refns="1" cloudfmin="0.000" cloudfmax="1.000" maxdoasrefmode="sza" east="false" center="false" west="false" backscan="false" />
      <cross_sections>
$(cross)      </cross_sections>
      <linear xpoly="$(polystr(w.poly_order))" xbase="none" xfit="$(xmlbool(w.poly_store[1]))" xerr="$(xmlbool(w.poly_store[2]))" xinvpoly="none" xinvbase="none" xinvfit="false" xinverr="false" offpoly="$(polystr(w.offl_order))" offbase="none" offfit="$(xmlbool(w.offl_store[1]))" offerr="$(xmlbool(w.offl_store[2]))" offizero="$(xmlbool(w.offl_ref))" />
      <nonlinear solfit="false" solinit="0.000" soldelt="0.001" solfstr="false" solestr="false"
                 $nl
                 comfit="false" cominit="0.000" comdelt="0.001" comstr="false" comestr="false"
                 u1fit="false" u1init="0.000" u1delt="0.001" u1str="false" u1estr="false"
                 u2fit="false" u2init="0.000" u2delt="0.001" u2str="false" u2estr="false"
                 ramfit="false" raminit="0.000" ramdelt="0.000" ramstr="false" ramestr="false"
                 comfile="" u1file="" u2file="" ramfile="" />
      <shift_stretches>
$(join(shift_stretch_xml(g) for g in w.shifts))      </shift_stretches>
      <gaps>
$(gaps)      </gaps>
      <outputs>
$(outs)      </outputs>
    </analysis_window>
"""
end

"A whole QDOAS project file around the `<analysis_window>` elements `windows_xml`."
function project_xml(; name, npix, calib, windows_xml, raw_spectra_xml, output, symbols, fields,
                     dark="", offset="", interpolation=:spline, security_gap=10, convergence=1e-4,
                     max_iterations=0, revert=false, date_format="MM/DD/YYYY", swath_name="QDOAS Results")
    syms = join("    <symbol name=\"$(xmlesc(s))\" descr=\"\" />\n" for s in symbols)
    flds = join("      <field name=\"$(xmlesc(f))\" />\n" for f in fields)
    """
<?xml version="1.0" encoding="UTF-8"?>
<qdoas>
  <paths>
  </paths>
  <symbols>
$(syms)  </symbols>
  <project name="$(xmlesc(name))" disable="false">
    <display spectra="true" data="true" calib="true" fits="true">
      <field name="name" />
      <field name="date" />
      <field name="starttime" />
      <field name="endtime" />
    </display>
    <selection>
      <sza min="0.000" max="0.000" delta="0.000" />
      <elevation min="0.000" max="0.000" tol="0.000" />
      <reference angle="90.000" tol="10.000" />
      <record min="0" max="0" />
      <cloud min="0.000" max="1.000" />
      <geolocation selected="none">
        <circle radius="0.000" long="0.000" lat="0.000" />
        <rectangle west="0.000" east="0.000" south="0.000" north="0.000" />
        <sites radius="0.000" />
      </geolocation>
    </selection>
    <analysis method="ODF" fit="none" unit="nm" interpolation="$(interpolation)" gap="$(security_gap)" converge="$(@sprintf("%.17g", convergence))" max_iterations="$(max_iterations)" spike_tolerance="999.9" >
    </analysis>
    <lowpass_filter selected="none">
    </lowpass_filter>
    <highpass_filter selected="none">
    </highpass_filter>
    <calibration ref="" method="ODF">
      <line shape="none" lorentzorder="1" slfFile="" />
      <display spectra="true" fits="true" residual="true" shiftsfp="true" />
      <polynomial shift="3" sfp="3" />
      <window min="300.0" max="340.0" intervals="4" custom_windows="" division="contiguous" size="10.00" />
      <preshift calculate="false" min="-3.0" max="3.0" />
      <cross_sections>
      </cross_sections>
      <linear xpoly="2" xbase="none" xfit="false" xerr="false" offpoly="none" offfit="false" offerr="false" offizero="false" />
      <shift_stretches>
      </shift_stretches>
      <gaps>
      </gaps>
      <outputs>
      </outputs>
    </calibration>
    <undersampling ref="" method="file" shift="0.000000" />
    <instrumental format="mfcstd" site="No Site Specified">
      <mfcstd size="$npix" revert="$(xmlbool(revert))" straylight="false" date="$(xmlesc(date_format))" lambda_min="0" lambda_max="0" calib="$(xmlesc(calib))" instr="" dark="$(xmlesc(dark))" offset="$(xmlesc(offset))" />
    </instrumental>
    <slit ref="" fwhmcor="false">
      <slit_func type="none">
      </slit_func>
    </slit>
    <output path="$(xmlesc(output))" anlys="true" calib="false" ref="false" conf="false" dirs="true" file="true" success="false" flux="" cic=" " bandWidth="1.0" swathName="$(xmlesc(swath_name))" fileFormat=".nc">
$(flds)    </output>
$(join(windows_xml))    <raw_spectra>
$(raw_spectra_xml)    </raw_spectra>
  </project>
</qdoas>
"""
end

"""
    write_qdoas_project(path, windows; lambda, reference, spectra=String[], filter="*.STD",
                        recursive=false, name="MiniDOAS.jl", calibration=nothing, dark="", offset="",
                        output="", interpolation=:spline, security_gap=10, convergence=1e-4,
                        max_iterations=0, revert=false, date_format="MM/DD/YYYY",
                        data_dir=<path without .xml>*"_files")
    write_qdoas_project(path, model::DOASModel; kw...)

Writes a QDOAS project file (MFC STD spectra, optical-density fit) for the analysis
windows `windows` (`FitWindow`s) or those of a `DOASModel`, which QDOAS opens and runs and
`parse_config` reads back.

  * `lambda`: the wavelength calibration, one value per pixel (written to `data_dir` as the
    project's calibration file unless `calibration` names one);
  * `reference`: a reference file, or the reference intensities on `lambda` (written to
    `data_dir`); a window's own `reference` file takes precedence;
  * `spectra`: a folder (with the file pattern `filter`, and its sub-folders when
    `recursive`) or a list of spectrum files;
  * cross sections that are not read from a file are written to `data_dir`, with every
    digit, so QDOAS and MiniDOAS.jl fit exactly the same numbers.

Paths in the project are absolute, as QDOAS expects. Returns `path`.
"""
function write_qdoas_project(path::AbstractString, windows::AbstractVector{FitWindow};
                             lambda::AbstractVector{<:Real}, reference::Union{AbstractString,AbstractVector{<:Real}},
                             spectra::Union{AbstractString,AbstractVector{<:AbstractString}}=String[],
                             filter::AbstractString="*.STD", recursive::Bool=false, name::AbstractString="MiniDOAS.jl",
                             calibration::Union{Nothing,AbstractString}=nothing, dark::AbstractString="",
                             offset::AbstractString="", output::AbstractString="", interpolation::Symbol=:spline,
                             security_gap::Integer=10, convergence::Real=1e-4, max_iterations::Integer=0,
                             revert::Bool=false, date_format::AbstractString="MM/DD/YYYY",
                             fields::AbstractVector{<:AbstractString}=DEFAULT_FIELDS,
                             data_dir::AbstractString=splitext(abspath(path))[1] * "_files")
    isempty(windows) && throw(ArgumentError("no analysis window"))
    allunique([w.name for w in windows]) || throw(ArgumentError("analysis window names must be unique"))
    for f in fields
        f in SUPPORTED_FIELDS || throw(ArgumentError("output field $f is not one MiniDOAS.jl writes"))
    end
    n = length(lambda)
    data_dir = abspath(data_dir)
    needs_dir = calibration === nothing || !(reference isa AbstractString) || any(w -> any(c -> isempty(c.file), w.xs), windows)
    needs_dir && mkpath(data_dir)
    calib = calibration === nothing ? write_calibration(joinpath(data_dir, "calibration.clb"), lambda) : abspath(calibration)
    ref = if reference isa AbstractString
        abspath(reference)
    else
        length(reference) == n || throw(DimensionMismatch("$n wavelengths, reference of $(length(reference))"))
        write_reference(joinpath(data_dir, "reference.ref"), lambda, reference)
    end
    safe(s) = replace(s, r"[^\w.+-]+" => "_")
    wxml = String[]
    for w in windows
        files = Dict{String,String}()
        for c in w.xs
            files[c.name] = !isempty(c.file) ? abspath(c.file) :
                write_cross_section(joinpath(data_dir, safe("$(w.name)_$(c.name)") * ".xs"), c)
        end
        push!(wxml, window_xml(w, isempty(w.refone) ? ref : abspath(w.refone), files))
    end
    raw = if spectra isa AbstractString
        isempty(spectra) ? "" :
            "      <directory name=\"$(xmlesc(abspath(spectra)))\" filters=\"$(xmlesc(filter))\" recursive=\"$(xmlbool(recursive))\" />\n"
    else
        join("      <file name=\"$(xmlesc(abspath(f)))\" />\n" for f in spectra)
    end
    symbols = unique(c.name for w in windows for c in w.xs)
    out = isempty(output) ? splitext(abspath(path))[1] * "_results.nc" : abspath(output)
    s = project_xml(; name=name, npix=n, calib=calib, windows_xml=wxml, raw_spectra_xml=raw, output=out,
                    symbols=symbols, fields=fields, dark=isempty(dark) ? "" : abspath(dark),
                    offset=isempty(offset) ? "" : abspath(offset), interpolation=interpolation,
                    security_gap=security_gap, convergence=convergence, max_iterations=max_iterations,
                    revert=revert, date_format=date_format)
    write(path, s)
    path
end

write_qdoas_project(path::AbstractString, w::FitWindow; kw...) = write_qdoas_project(path, [w]; kw...)

function write_qdoas_project(path::AbstractString, m::DOASModel; reference=nothing, kw...)
    p = m.spec
    ref = reference !== nothing ? reference :
          !isempty(m.reference) ? m.reference :
          throw(ArgumentError("pass the reference (file or intensities) for a model made from a QDOAS project"))
    write_qdoas_project(path, p.windows; lambda=m.lambda, reference=ref, interpolation=p.interpolation,
                        security_gap=p.security_gap, convergence=p.convergence, max_iterations=p.max_iterations,
                        name=p.project_name, kw...)
end

# =============================================================================
# Editing a project: the changes the GUI makes to analysis windows
# =============================================================================

"""
Applies window edits to a project: per window name, `enabled`, `lmin`, `lmax`,
`poly_order` and `exclude` (cross sections left out of the fit).
"""
function apply_edits(p::ProjectSpec, edits)
    edits === nothing && return p
    ws = FitWindow[]
    for w in p.windows
        e = get(edits, w.name, nothing)
        if e === nothing
            push!(ws, w); continue
        end
        get(e, "enabled", true) || continue
        excl = Set(String.(collect(get(e, "exclude", String[]))))
        lmin = Float64(get(e, "lmin", w.lmin)); lmax = Float64(get(e, "lmax", w.lmax))
        lmax - lmin >= 1.0 || throw(ArgumentError("$(w.name): the fitting interval $lmin-$lmax nm is too narrow"))
        shifts = ShiftStretch[]
        for g in w.shifts
            syms = filter(s -> !(s in excl), g.symbols)
            isempty(syms) && continue
            push!(shifts, ShiftStretch((f == :symbols ? syms : getfield(g, f) for f in fieldnames(ShiftStretch))...))
        end
        w2 = with_range(w, lmin, lmax)
        push!(ws, rewindow(w2; poly_order=Int(get(e, "poly_order", w.poly_order)),
                           xs=filter(c -> !(c.name in excl), w.xs), shifts=shifts))
    end
    isempty(ws) && throw(ArgumentError("all analysis windows are switched off"))
    respec(p; windows=ws)
end

regex_escape(s) = replace(s, r"([\\^$.|?*+()\[\]{}])" => s"\\\1")

"Writes window edits (see `apply_edits`) into a copy of a project file, as QDOAS would save them."
function save_edited_project(src_xml::AbstractString, dst_xml::AbstractString, edits)
    s = read(src_xml, String)
    setattr(tag, k, v) = occursin(Regex("\\b$k=\""), tag) ?
        replace(tag, Regex("\\b$k=\"[^\"]*\"") => "$k=\"$v\"") : replace(tag, r"\s*>$" => " $k=\"$v\" >")
    for (name, e) in pairs(something(edits, Dict()))
        name = String(name)
        m = match(Regex("<analysis_window\\s+name=\"" * regex_escape(xmlesc(name)) * "\"[^>]*>.*?</analysis_window>", "s"), s)
        m === nothing && continue
        blk = m.match
        head = match(r"^<analysis_window\b[^>]*>", blk).match
        h2 = setattr(head, "disable", get(e, "enabled", true) ? "false" : "true")
        haskey(e, "lmin") && (h2 = setattr(h2, "min", @sprintf("%.3f", e["lmin"])))
        haskey(e, "lmax") && (h2 = setattr(h2, "max", @sprintf("%.3f", e["lmax"])))
        b2 = replace(blk, head => h2; count=1)
        if haskey(e, "poly_order")
            po = Int(e["poly_order"])
            b2 = replace(b2, r"(<linear\b[^>]*?\bxpoly=\")[^\"]*\"" => SubstitutionString("\\g<1>$(po < 0 ? "none" : po)\""))
        end
        for sym in String.(collect(get(e, "exclude", String[])))
            q = regex_escape(xmlesc(sym))
            b2 = replace(b2, Regex("[ \\t]*<cross_section\\s+sym=\"$q\"[^>]*/>[ \\t]*\\r?\\n?") => "")
            b2 = replace(b2, Regex("[ \\t]*<output\\s+sym=\"$q\"[^>]*/>[ \\t]*\\r?\\n?") => "")
            b2 = replace(b2, Regex("[ \\t]*<symbol\\s+name=\"$q\"\\s*/>[ \\t]*\\r?\\n?") => "")
        end
        # a shift/stretch row left without symbols is dropped, as QDOAS would not save it
        b2 = replace(b2, r"[ \t]*<shift_stretch\b[^>]*>\s*</shift_stretch>[ \t]*\r?\n?"s => "")
        s = replace(s, blk => b2; count=1)
    end
    write(dst_xml, s)
    dst_xml
end
