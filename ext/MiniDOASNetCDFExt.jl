# This file is part of MiniDOAS.jl (Differential Optical Absorption Spectroscopy in Julia).
# Copyright (c) 2026 Sunandan Mahant <sunandanmahant@outlook.com>
# Derived from QDOAS, Copyright (C) 1994-2025 BIRA-IASB and S[&]T (BSD-3-Clause).
# See LICENSE and NOTICE.md.

# =============================================================================
# NetCDF output in QDOAS's layout (output_netcdf.cpp), so that everything reading
# qdoas_results.nc works unchanged; and reading such files back.
# =============================================================================

module MiniDOASNetCDFExt

using MiniDOAS
using MiniDOAS: DOASResult, QDOAS_FILL_DOUBLE, VERSION_STRING
using NCDatasets

const QDOAS_FILL_INT = Int32(-2147483647)

# NCDatasets takes one global lock around every netCDF call. Several analyses writing
# at once then hand that lock back and forth hundreds of times per file, and a 40 ms
# write took 600 ms with six in flight. netCDF runs one call at a time regardless, so
# write whole files one at a time instead.
const NC_WRITE_LOCK = ReentrantLock()

function write_netcdf(path::AbstractString, res::DOASResult; config_path::AbstractString="")
    lock(NC_WRITE_LOCK) do
        p = res.spec
        rec = res.record
        # without a spectrum record (a fit of an in-memory spectrum) there is no name, date or time
        fields = rec === nothing ? filter(f -> !(f in ("name", "date", "time")), p.fields) : p.fields
        isfile(path) && rm(path)
        NCDataset(path, "c") do ds
            g = defGroup(ds, p.swath_name)
            g.attrib["CreationTime"] = Libc.strftime("%a %b %e %H:%M:%S %Y\n", time())
            g.attrib["Qdoas"] = "Results obtained using MiniDOAS.jl $(VERSION_STRING), a Julia port of " *
                                "the Qdoas analysis (BIRA-IASB and S[&]T); not produced by Qdoas itself"
            g.attrib["InputFile"] = p.spectrum
            g.attrib["QDOASConfigProject"] = p.project_name
            g.attrib["QDOASConfig"] = String(config_path)

            defDim(g, "n_alongtrack", 1)
            "date" in fields && defDim(g, "date", 3)
            "time" in fields && defDim(g, "time", 3)
            if any(wr -> any(kv -> startswith(first(kv), "Stretch(") || startswith(first(kv), "Err Stretch("), wr.nc_vars), res.windows)
                defDim(g, "2", 2)
            end
            for wr in res.windows
                for (k, v) in wr.nc_vars
                    k == "residual_spectrum" && defDim(g, "$(wr.name).n_datapoint", length(v))
                end
            end

            for f in fields
                if f == "name"
                    v = defVar(g, "Name", String, ("n_alongtrack",))
                    v[1] = rec.name
                elseif f == "date"
                    v = defVar(g, "Date (DD-MM-YYYY)", Int32, ("date", "n_alongtrack"); fillvalue=QDOAS_FILL_INT)
                    v[:, 1] = Int32.(collect(rec.date))
                elseif f == "time"
                    v = defVar(g, "Time (hh:mm:ss)", Int32, ("time", "n_alongtrack"); fillvalue=QDOAS_FILL_INT)
                    v[:, 1] = Int32.(collect(rec.time))
                end
            end

            for wr in res.windows
                gw = defGroup(g, wr.name)
                for (k, v) in wr.nc_attribs
                    gw.attrib[k] = v
                end
                for (k, v) in wr.nc_vars
                    if v isa Int32
                        x = defVar(gw, k, Int32, ("n_alongtrack",); fillvalue=QDOAS_FILL_INT)
                        x[1] = v
                    elseif k == "residual_spectrum"
                        x = defVar(gw, k, Float64, ("$(wr.name).n_datapoint", "n_alongtrack"); fillvalue=QDOAS_FILL_DOUBLE)
                        x[:, 1] = v
                    elseif v isa AbstractVector
                        x = defVar(gw, k, Float64, ("2", "n_alongtrack"); fillvalue=QDOAS_FILL_DOUBLE)
                        x[:, 1] = v
                    else
                        x = defVar(gw, k, Float64, ("n_alongtrack",); fillvalue=QDOAS_FILL_DOUBLE)
                        x[1] = Float64(v)
                    end
                end
            end
        end
        path
    end
end

"A variable's values with fill values as NaN (floats) or missing (others)."
function values_of(v)
    a = Array(v)                                  # NCDatasets returns missing for fill values
    T = Base.nonmissingtype(eltype(a))
    if T <: AbstractFloat
        b = Float64[x === missing ? NaN : Float64(x) for x in a]
        # QDOAS's double fill value, also where a file declares none
        b[b .>= 0.99 * QDOAS_FILL_DOUBLE] .= NaN
        reshape(b, size(a))
    else
        a
    end
end

function read_qdoas_netcdf(path::AbstractString; group::Union{Nothing,AbstractString}=nothing)
    NCDataset(path) do ds
        gname = group === nothing ? (isempty(ds.group) ? nothing : first(keys(ds.group))) : String(group)
        g = gname === nothing ? ds : ds.group[gname]
        attributes = Dict{String,Any}(k => v for (k, v) in g.attrib)
        fields = Dict{String,Any}(k => values_of(g[k]) for k in keys(g))
        windows = Dict{String,Dict{String,Any}}()
        for (wname, gw) in g.group
            d = Dict{String,Any}(k => values_of(gw[k]) for k in keys(gw))
            d["attributes"] = Dict{String,Any}(k => v for (k, v) in gw.attrib)
            windows[wname] = d
        end
        (group=gname === nothing ? "" : gname, attributes=attributes, fields=fields, windows=windows)
    end
end

end # module
