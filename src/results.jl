# This file is part of MiniDOAS.jl (Differential Optical Absorption Spectroscopy in Julia).
# Copyright (c) 2026 Sunandan Mahant <sunandanmahant@outlook.com>
# Derived from QDOAS, Copyright (C) 1994-2025 BIRA-IASB and S[&]T (BSD-3-Clause).
# See LICENSE and NOTICE.md.

# =============================================================================
# ANALYSE_CurFitMethod + ANALYSE_Spectrum (per window)
# =============================================================================

"""
    WindowResult

The fit of one analysis window to one spectrum:

  * `slcol[x]`, `slerr[x]`   slant column and its 1σ error for every linear term `x`: the
                             cross sections, the polynomial `x0`, `x1`, ... and the linear
                             offset `offl0`, ... (`columns` lists them in fit order)
  * `rms`, `chisqr`, `niter` fit quality and number of iterations
  * `shift[s]`, `stretch[s]` shift (nm) and stretch of each shifted target `s`, with
                             `shift_err`, `stretch_err` (an error of exactly 1 marks a
                             value that was not fitted, as in QDOAS)
  * `params`, `param_err`    the non-linear offsets (`Offset0`, ...)
  * `residual`               the residual (optical depth) as QDOAS stores it
  * `pixels`                 the fitted pixels
  * `curves`                 with `curves=true`: the fitted curves for plotting -
                             `lambda`, `od`, `model`, `residual` and one entry per term in
                             `components` (`name`, `kind`, `value`, `scd`, `err`)

The polynomial is reported, as by QDOAS, for ln(I/I₀) (the opposite sign of the fitted
optical depth), with the normalisation of spectrum and reference in `x0`.
"""
struct WindowResult
    name::String
    ok::Bool
    rms::Float64
    chisqr::Float64
    niter::Int
    residual::Vector{Float64}        # as QDOAS stores it: DimL values from SvdPDeb
    slcol::Dict{String,Float64}
    slerr::Dict{String,Float64}
    shift::Dict{String,Float64}      # sym -> value (fitted or fixed)
    shift_err::Dict{String,Float64}
    stretch::Dict{String,Float64}
    stretch_err::Dict{String,Float64}
    params::Dict{String,Float64}     # non-linear offsets
    param_err::Dict{String,Float64}
    pixels::Vector{Int}              # 1-based fitted pixels
    nc_vars::Vector{Pair{String,Any}}      # what QDOAS would write, in its order
    nc_attribs::Vector{Pair{String,String}}
    columns::Vector{String}          # linear parameters in fit order (cross sections, x0.., offl0..)
    curves::Any                      # fitted curves for plotting (analyze(...; curves=true)), else nothing
end

const QDOAS_FILL_DOUBLE = 9.9692099683868690e+306

function fit_window(w::Window, p::ProjectSpec, spec, ref)
    cubic = p.interpolation === :spline
    lam = w.lam
    spline_spec = cubic ? spline_deriv2(lam, spec) : zeros(length(lam))
    spline_ref = cubic ? spline_deriv2(lam, ref) : zeros(length(lam))
    wk = Work(w)
    F(A) = analyse_function(w, wk, spec, spline_spec, ref, spline_ref, A, cubic)

    nA = length(w.nl)
    st = FitState(Float64[], Float64[], 0.001, 0.0, ones(nA))
    niter = 0
    A = initial_params(w)
    if nA == 0
        st.Yfit, st.c = F(A)
        st.chisqr = fchisq(st.Yfit, w.nfree)
    else
        deltaA = [q.delta for q in w.nl]
        while true
            old = st.chisqr
            curfit!(st, w, wk, F, A, deltaA, niter)
            deltaA .*= 0.4
            niter += 1
            (st.chisqr != 0 && abs(st.chisqr - old) / st.chisqr > p.convergence &&
             (p.max_iterations == 0 || niter < p.max_iterations)) || break
        end
    end
    st, niter, wk, A
end

function window_result(w::Window, p::ProjectSpec, spec, ref; spe_norm=1.0, ref_norm=1.0, curves::Bool=false)
    st, niter, wk, A = fit_window(w, p, spec, ref)
    n = length(w.lam)

    absolu = zeros(n)
    for (k, i) in enumerate(w.pix)
        absolu[i] = st.Yfit[k] - 0.0
    end

    # Spike removal (ANALYSE_Spectrum: exclude pixels above tol*RMS and refit) is not
    # implemented. With the usual setting (999.9, off) it never fires; if it would, the
    # spectrum is refused rather than analysed differently from QDOAS.
    rms0 = sqrt(sum(abs2, st.Yfit) / length(w.pix))
    any(abs(v) > rms0 * p.spike_tolerance for v in st.Yfit) &&
        unsupported("spike removal would trigger in $(w.spec.name)")

    rms = sqrt(sum(abs2, st.Yfit) / length(w.pix))
    deb = first(first(w.ranges))
    residual = absolu[deb:min(deb + length(w.pix) - 1, n)]

    slcol = Dict{String,Float64}(); slerr = Dict{String,Float64}()
    for (j, c) in enumerate(w.cols)
        v = st.c[j]
        # ANALYSE_CurFitMethod reports the polynomial as log(spec/ref) rather than
        # log(ref/spec) - "better to compare with Matlab's polyfit" - so its sign is
        # flipped, and x0 also absorbs the two normalisation factors.
        if c.kind in (:poly, :offl) && abs(spe_norm) > EPSILON && ref_norm / spe_norm > EPSILON
            v = -v
            c.kind === :poly && c.poly_power == 0 && (v -= log(ref_norm / spe_norm))
        end
        slcol[c.name] = v
        slerr[c.name] = sqrt(wk.sigmasq[j] * st.chisqr)
    end

    shift = Dict{String,Float64}(); shift_err = Dict{String,Float64}()
    stretch = Dict{String,Float64}(); stretch_err = Dict{String,Float64}()
    function record(sym, idx, init)
        shift[sym] = idx[1] > 0 ? A[idx[1]] : init[1]
        shift_err[sym] = idx[1] > 0 ? st.sigmaA[idx[1]] : 1.0
        stretch[sym] = idx[2] > 0 ? A[idx[2]] * w.sf1 : (init[2] / w.sf1) * w.sf1
        stretch_err[sym] = idx[2] > 0 ? st.sigmaA[idx[2]] * w.sf1 : 1.0
    end
    for g in w.spec.shifts, sym in g.symbols
        if lowercase(sym) == "spectrum"
            record(sym, w.spec_shift, w.spec_init)
        elseif lowercase(sym) == "ref"
            record(sym, w.ref_shift, w.ref_init)
        else
            c = w.cols[findfirst(x -> lowercase(x.name) == lowercase(sym), w.cols)]
            record(c.name, (c.fit_shift, c.fit_stretch, c.fit_stretch2),
                   (c.init_shift, c.init_stretch, c.init_stretch2))
        end
    end

    params = Dict{String,Float64}(); perr = Dict{String,Float64}()
    facts = (1.0, w.fact1, w.fact2)
    for k in 1:3
        if w.off_idx[k] > 0
            params["Offset$(k-1)"] = A[w.off_idx[k]] / facts[k]
            perr["Offset$(k-1)"] = st.sigmaA[w.off_idx[k]] / facts[k]
        end
    end

    nc_vars, nc_attribs = output_layout(w, p, st, niter, rms, residual, A, slcol, slerr)
    WindowResult(w.spec.name, true, rms, st.chisqr, niter, residual, slcol, slerr,
                 shift, shift_err, stretch, stretch_err, params, perr, copy(w.pix),
                 nc_vars, nc_attribs, [c.name for c in w.cols],
                 curves ? fit_curves(w, p, spec, ref, A, slcol, slerr) : nothing)
end

"""
The fit as curves over the fitted pixels, for plotting: optical depth ln(I0/I) (after
the shifts, stretches and offsets the fit found), the model, the residual, and each
linear term's contribution. Re-evaluates the model once at the final parameters;
nothing here feeds back into the numbers above.
"""
function fit_curves(w::Window, p::ProjectSpec, spec, ref, A, slcol, slerr)
    cubic = p.interpolation === :spline
    wk = Work(w)
    ss = cubic ? spline_deriv2(w.lam, spec) : zeros(length(w.lam))
    sr = cubic ? spline_deriv2(w.lam, ref) : zeros(length(w.lam))
    resid, c = analyse_function(w, wk, spec, ss, ref, sr, copy(A), cubic)
    Am = wk.Amat
    model = Am * c
    comps = [(name=col.name, kind=col.kind, value=Am[:, j] .* c[j],
              scd=get(slcol, col.name, NaN), err=get(slerr, col.name, NaN)) for (j, col) in enumerate(w.cols)]
    (lambda=w.lam[w.pix], od=resid .+ model, model=model, residual=resid, components=comps,
     fit_range=(w.spec.lmin, w.spec.lmax))
end

"""
Mirror output.c register_analysis_output + register_cross_results: which variables
QDOAS writes for a window, with what values, in which order. Symbols follow the
window's cross-reference table: cross sections, polynomial, non-linear offsets,
then shift-only symbols (Spectrum, Ref) in shift-group order.
"""
function output_layout(w::Window, p::ProjectSpec, st, niter, rms, residual, A, slcol, slerr)
    ws = w.spec
    vars = Pair{String,Any}[]
    for f in p.fields
        f == "chi" && push!(vars, "Chi" => st.chisqr)
        f == "rms" && push!(vars, "RMS" => rms)
        f == "iter_number" && push!(vars, "iter" => Int32(niter))
        f == "error_flag" && push!(vars, "processing_error" => Int32(0))
        f == "residual_spectrum" && ws.save_residuals && push!(vars, "residual_spectrum" => residual)
    end

    # which symbol carries each shift group's stored outputs (the first one)
    first_of = Dict{String,ShiftStretch}()
    for g in ws.shifts
        isempty(g.symbols) || (first_of[lowercase(g.symbols[1])] = g)
    end
    fitidx(c) = (c.fit_shift, c.fit_stretch, c.fit_stretch2)
    # get_stretches tests TabCross[0] (the first entry), not the symbol itself
    tc0 = isempty(w.cols) ? (0, 0, 0) : fitidx(w.cols[1])
    item(i) = i > 0 ? i - 1 : -1                       # to QDOAS's 0-based / ITEM_NONE
    function shift_vars!(sym, idx, init)
        g = get(first_of, lowercase(sym), nothing)
        g === nothing && return
        shv = idx[1] > 0 ? A[idx[1]] : init[1]
        she = idx[1] > 0 ? st.sigmaA[idx[1]] : 1.0
        stv = [idx[2] > 0 ? A[idx[2]] * w.sf1 : (init[2] / w.sf1) * w.sf1,
               idx[3] > 0 ? A[idx[3]] * w.sf2 : (init[3] / w.sf2) * w.sf2]
        ste = [idx[2] > 0 ? st.sigmaA[idx[2]] * w.sf1 : 1.0, idx[3] > 0 ? st.sigmaA[idx[3]] * w.sf2 : 1.0]
        item(tc0[2]) == 0 && (stv[1] = QDOAS_FILL_DOUBLE)
        item(tc0[3]) == 0 && (stv[2] = QDOAS_FILL_DOUBLE)
        g.sh_store && push!(vars, "Shift($sym)" => shv)
        g.sh_store && g.err_store && push!(vars, "Err Shift($sym)" => she)
        g.st_store && push!(vars, "Stretch($sym)" => stv)
        g.st_store && g.err_store && push!(vars, "Err Stretch($sym)" => ste)
    end

    attribs = Pair{String,String}[]
    for c in w.cols
        if c.kind === :xs
            scol, serr, sfact = get(ws.outputs, c.name, (false, false, 1.0))
            scol && push!(vars, "SlCol($(c.name))" => (sfact != 0 ? slcol[c.name] / sfact : QDOAS_FILL_DOUBLE))
            serr && push!(vars, "SlErr($(c.name))" => (sfact != 0 ? slerr[c.name] / sfact : QDOAS_FILL_DOUBLE))
            shift_vars!(c.name, fitidx(c), (c.init_shift, c.init_stretch, c.init_stretch2))
            xsf = ws.xs[findfirst(x -> x.name == c.name, ws.xs)].file
            push!(attribs, c.name => xsf)
        else
            store = c.kind === :poly ? ws.poly_store : ws.offl_store
            store[1] && push!(vars, "SlCol($(c.name))" => slcol[c.name])
            store[2] && push!(vars, "SlErr($(c.name))" => slerr[c.name])
        end
    end
    names = ("Offset (Constant)", "Offset (Order 1)", "Offset (Order 2)")
    facts = (1.0, w.fact1, w.fact2)
    for k in 1:3
        (ws.off[k][1] || ws.off[k][2] != 0) || continue
        fstr, estr = ws.off_store[k]
        fitted = w.off_idx[k] > 0
        val = fitted ? A[w.off_idx[k]] / facts[k] : ws.off[k][2]
        err = fitted ? st.sigmaA[w.off_idx[k]] / facts[k] : 1.0
        fstr && push!(vars, names[k] => val)
        fstr && estr && push!(vars, "Err($(names[k]))" => err)
    end
    for g in ws.shifts, sym in g.symbols
        ls = lowercase(sym)
        if ls == "spectrum"
            shift_vars!(sym, w.spec_shift, w.spec_init)
        elseif ls == "ref"
            shift_vars!(sym, w.ref_shift, w.ref_init)
        end
    end
    push!(attribs, "fitting window range" => @sprintf("%.3f : %.3f", ws.lmin, ws.lmax))
    vars, attribs
end
