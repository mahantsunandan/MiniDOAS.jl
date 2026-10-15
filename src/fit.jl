# This file is part of MiniDOAS.jl (Differential Optical Absorption Spectroscopy in Julia).
# Copyright (c) 2026 Sunandan Mahant <sunandanmahant@outlook.com>
# Derived from QDOAS, Copyright (C) 1994-2025 BIRA-IASB and S[&]T (BSD-3-Clause).
# See LICENSE and NOTICE.md.

# =============================================================================
# Per-window setup (ANALYSE_Load*, ANALYSE_SvdInit)
# =============================================================================

mutable struct Column
    name::String
    kind::Symbol                  # :xs or :poly
    vector::Vector{Float64}
    deriv2::Vector{Float64}
    poly_power::Int
    fit_shift::Int                # index into non-linear params, 0 = none
    fit_stretch::Int
    fit_stretch2::Int
    init_shift::Float64
    init_stretch::Float64         # normalised by StretchFact1 during the fit
    init_stretch2::Float64
end

struct NLParam
    name::String
    delta::Float64
    init::Float64
    minp::Float64
    maxp::Float64
end

struct Window
    spec::FitWindow
    lam::Vector{Float64}
    ranges::Vector{UnitRange{Int}}   # 1-based pixel ranges (spectral windows minus gaps)
    pix::Vector{Int}                 # all fitted pixels, in order
    limmin::Int
    limmax::Int
    lambda0::Float64
    cols::Vector{Column}             # cross sections then polynomial terms
    nl::Vector{NLParam}
    # non-linear roles
    spec_shift::Tuple{Int,Int,Int}   # (shift, stretch, stretch2) param idx for Spectrum
    spec_init::Tuple{Float64,Float64,Float64}
    ref_shift::Tuple{Int,Int,Int}
    ref_init::Tuple{Float64,Float64,Float64}
    off_idx::NTuple{3,Int}           # param index per offset order (0 = not fitted)
    off_init::NTuple{3,Float64}
    off_active::Int                  # highest active offset order, -1 none
    fact1::Float64                   # Fact for offset order 1 / 2
    fact2::Float64
    sf1::Float64                     # StretchFact1 / StretchFact2
    sf2::Float64
    nfree::Int
    matrix_depends_on_nl::Bool       # NP != 0
end

"""
A cross section on the analysis grid `lam` with the second derivatives of its spline:
used as it is (`:none`), copied when it is already on the grid, else interpolated.
"""
function xs_on_grid(xl, xv, action, lam, interpolation, sym)
    n = length(lam)
    if action === :none
        length(xv) >= n || throw(EngineError("cross section $sym: $(length(xv)) points, need $n"))
        v = xv[1:n]
    elseif length(xl) == n && all(i -> (lam[i] == 0 ? xl[i] == lam[i] :
                                        abs(xl[i] - lam[i]) / abs(lam[i]) <= 1e-7), 1:n)
        v = copy(xv)
    else
        xd2 = spline_deriv2(xl, xv)
        cub = interpolation === :spline
        v = [spline_at(xl, xv, xd2, lam[i], cub) for i in 1:n]
    end
    (v, spline_deriv2(lam, v))
end

const XS_MEMO_LOCK = ReentrantLock()

"""
The cross section `c` on the grid `lam`: from the file cache when it is read from a file,
else remembered in `memo`. The memo is keyed on `c` itself (whose arrays compare by
identity, and which the key keeps alive) and on the grid's values.
"""
function xs_vectors(c::CrossSection, lam::Vector{Float64}, interpolation::Symbol, memo::Dict)
    if !isloaded(c)
        key = (c.file, c.action, interpolation, hash(lam))
        return cached(key) do
            xl, xv, _ = load_xs(c.file)
            xs_on_grid(xl, xv, c.action, lam, interpolation, c.name)
        end
    end
    key = (c, interpolation, lam)
    lock(XS_MEMO_LOCK) do
        get!(memo, key) do
            xs_on_grid(c.lambda, c.sigma, c.action, lam, interpolation, c.name)
        end
    end
end

function build_window(ws::FitWindow, p::ProjectSpec, lam::Vector{Float64}, xs_cache::Dict)
    n = length(lam)

    # --- ranges (ANALYSE_LoadGaps) -------------------------------------------
    fen = [(min(ws.lmin, ws.lmax), max(ws.lmin, ws.lmax))]
    for (g1, g2) in ws.gaps
        l1, l2 = min(g1, g2), max(g1, g2)
        iw = findfirst(w -> l1 > w[1] && l2 < w[2], fen)
        iw === nothing && throw(EngineError("gap $l1-$l2 outside window $(ws.name)"))
        a, b = fen[iw]
        splice!(fen, iw, [(a, l1), (l2, b)])
    end
    ranges = [fnpixel(lam, a, :after):fnpixel(lam, b, :before) for (a, b) in fen]
    pix = reduce(vcat, collect.(ranges))
    diml = length(pix)
    deb, fin = first(first(ranges)), last(last(ranges))
    limmin = max(deb - p.security_gap, 1)
    limmax = min(fin + p.security_gap, n)
    lambda0 = ws.lambda0

    # --- columns: cross sections then polynomial -----------------------------
    cols = Column[]
    for c in ws.xs
        vec, d2 = xs_vectors(c, lam, p.interpolation, xs_cache)
        push!(cols, Column(c.name, :xs, vec, d2, 0, 0, 0, 0, 0.0, 0.0, 0.0))
    end
    for k in 0:ws.poly_order
        push!(cols, Column("x$k", :poly, ones(n), zeros(n), k, 0, 0, 0, 0.0, 0.0, 0.0))
    end
    for k in 0:ws.offl_order                      # ANALYSE_LoadLinear, second entry
        push!(cols, Column("offl$k", :offl, ones(n), zeros(n), k, 0, 0, 0, 0.0, 0.0, 0.0))
    end
    ncols = length(cols)

    # --- non-linear parameters (NF order: Sol, offsets, then shift/stretch) ----
    nl = NLParam[]
    off_idx = [0, 0, 0]
    for k in 1:3
        fit, init, delt = ws.off[k]
        if fit
            push!(nl, NLParam("Offset$(k-1)", delt, init, 0.0, 0.0))
            off_idx[k] = length(nl)
        end
    end
    off_active = -1
    for k in 1:3
        (ws.off[k][1] || ws.off[k][2] != 0) && (off_active = k - 1)
    end

    spec_shift = (0, 0, 0); spec_init = (0.0, 0.0, 0.0)
    ref_shift = (0, 0, 0);  ref_init = (0.0, 0.0, 0.0)
    np = 0
    for g in ws.shifts
        ish = ist = ist2 = 0
        for (k, sym) in enumerate(g.symbols)
            if k == 1
                g.sh_fit && (push!(nl, NLParam("Shift($sym)", g.sh_delta, g.sh_init, g.sh_min, g.sh_max)); ish = length(nl))
                g.st_type >= 1 && (push!(nl, NLParam("Stretch($sym)", g.st_delta, g.st_init, 0.0, 0.0)); ist = length(nl))
                g.st_type == 2 && (push!(nl, NLParam("Stretch2($sym)", g.st_delta2, g.st_init2, 0.0, 0.0)); ist2 = length(nl))
            end
            if lowercase(sym) == "spectrum"
                spec_shift = (ish, ist, ist2); spec_init = (g.sh_init, g.st_init, g.st_init2)
            elseif lowercase(sym) == "ref"
                ref_shift = (ish, ist, ist2); ref_init = (g.sh_init, g.st_init, g.st_init2)
            else
                ic = findfirst(c -> lowercase(c.name) == lowercase(sym), cols)
                ic === nothing && unsupported("shift on unknown symbol $sym in $(ws.name)")
                col = cols[ic]
                col.fit_shift, col.fit_stretch, col.fit_stretch2 = ish, ist, ist2
                col.init_shift, col.init_stretch, col.init_stretch2 = g.sh_init, g.st_init, g.st_init2
                np += (ish > 0) + (ist > 0) + (ist2 > 0)
            end
        end
    end

    nfit = ncols + length(nl)
    nfree = diml - nfit
    nfree > 0 || throw(EngineError("no degrees of freedom in $(ws.name)"))

    # --- normalisation factors (ANALYSE_SvdInit) ------------------------------
    norm1 = norm2 = 0.0
    for i in pix
        dx = (lam[i] - lambda0) * (lam[i] - lambda0)
        norm1 += dx
        norm2 += dx * dx
    end
    sf1 = sf2 = 0.0
    for j in limmin:limmax
        dx = lam[j] - lambda0
        dx = lam[j] - lambda0 - 0.0 * dx - 0.0 * dx * dx   # Feno->Stretch = 0 for ground-based
        dx *= dx
        sf1 += dx
        sf2 += dx * dx
    end
    (norm1 <= 0 || norm2 <= 0 || sf1 <= 0 || sf2 <= 0) && throw(EngineError("sqrt argument in SvdInit"))
    sf1 = 1.0 / sqrt(sf1)
    sf2 = 1.0 / sqrt(sf2)

    Window(ws, lam, ranges, pix, limmin, limmax, lambda0, cols, nl,
           spec_shift, spec_init, ref_shift, ref_init,
           (off_idx[1], off_idx[2], off_idx[3]),
           (ws.off[1][2], ws.off[2][2], ws.off[3][2]), off_active,
           sqrt(norm1), sqrt(norm2), sf1, sf2, nfree,
           # the radiance-normalised linear offset depends on the shifted spectrum,
           # so QDOAS never caches that decomposition (Decomp stays 1)
           np > 0 || (ws.offl_order >= 0 && !ws.offl_ref))
end

"Initial non-linear parameter vector (Fitp) with QDOAS's stretch normalisation."
function initial_params(w::Window)
    A = [q.init for q in w.nl]
    for (k, q) in enumerate(w.nl)
        startswith(q.name, "Stretch(") && (A[k] = q.init / w.sf1)
        startswith(q.name, "Stretch2(") && (A[k] = q.init / w.sf2)
    end
    A
end

# =============================================================================
# ShiftVector / ANALYSE_Function
# =============================================================================

"ShiftVector: target[LimMin:LimMax] = source evaluated at the shifted grid."
function shift_vector!(target, w::Window, source, deriv, dsh, dst, dst2, cubic)
    lam = w.lam; l0 = w.lambda0
    @inbounds for j in w.limmin:w.limmax
        x0 = lam[j] - l0
        y = lam[j] - (0.0 + 0.0 * x0 + 0.0 * x0 * x0)      # second shift/stretch = 0
        x0 = (y - l0 + 0.0)
        xs = y - (dsh + dst * x0 * w.sf1 + dst2 * x0 * x0 * w.sf2)
        target[j] = spline_at(lam, source, deriv, xs, cubic)
    end
    target
end

param(A, idx, init) = idx > 0 ? A[idx] : init

const QRT = typeof(qr(zeros(2, 2), ColumnNorm()))

"""
Per-call scratch space; also carries the linear system between calls when the
design matrix does not depend on the non-linear parameters (Feno->Decomp = 0).
"""
mutable struct Work
    spec_interp::Vector{Float64}
    ref_shifted::Vector{Float64}
    xs_tmp::Vector{Float64}
    Amat::Matrix{Float64}            # DimL x DimC, unnormalised columns
    norms::Vector{Float64}
    qrf::Union{Nothing,QRT}
    sigmasq::Vector{Float64}
    decomposed::Bool
end

function Work(w::Window)
    n = length(w.lam)
    Work(zeros(n), zeros(n), zeros(n), zeros(length(w.pix), length(w.cols)),
         zeros(length(w.cols)), nothing, zeros(length(w.cols)), false)
end

"""
ANALYSE_Function (OD fitting): returns (Yfit, fitParamsC) for non-linear
parameters `A`. Yfit is the residual log(I0) - log(I) - sum(A_j c_j).
"""
function analyse_function(w::Window, wk::Work, spec, spline_spec, ref, spline_ref, A, cubic)
    pix = w.pix
    npts = length(pix)
    ymean = 0.0                                   # ANALYSE_CurFitMethod: mean of RefTrav
    for i in pix
        ymean += ref[i]
    end
    ymean /= npts

    # spectrum shift/stretch
    si = wk.spec_interp
    fill!(si, 0.0)
    sh = param(A, w.spec_shift[1], w.spec_init[1])
    st = param(A, w.spec_shift[2], w.spec_init[2] / w.sf1)
    st2 = param(A, w.spec_shift[3], w.spec_init[3] / w.sf2)
    shift_vector!(si, w, spec, spline_spec, sh, st, st2, cubic)

    # mean over fitted pixels
    xmean = 0.0
    for i in pix
        xmean += si[i]
    end
    xmean /= npts

    # non-linear offset
    if w.off_active >= 0
        o0 = param(A, w.off_idx[1], w.off_init[1])
        o1 = w.off_idx[2] > 0 ? A[w.off_idx[2]] / w.fact1 : w.off_init[2]
        o2 = w.off_idx[3] > 0 ? A[w.off_idx[3]] / w.fact2 : w.off_init[3]
        for i in w.limmin:w.limmax
            dx = w.lam[i] - w.lambda0
            off = o0
            w.off_active >= 1 && (off += o1 * dx)
            w.off_active >= 2 && (off += o2 * dx * dx)
            si[i] -= off * xmean
        end
    end

    spec_nolog = si[pix]                          # backup before the logarithm
    for i in w.limmin:w.limmax
        si[i] <= 0 && throw(EngineError("log of non-positive spectrum at pixel $i"))
        si[i] = log(si[i])
    end
    X = si[pix]

    # design matrix (rebuilt whenever it can depend on the parameters)
    if !wk.decomposed || w.matrix_depends_on_nl
        Am = wk.Amat
        prev = 0
        for (jc, c) in enumerate(w.cols)
            if c.kind === :xs
                shx = param(A, c.fit_shift, c.init_shift)
                stx = param(A, c.fit_stretch, c.init_stretch / w.sf1)
                st2x = param(A, c.fit_stretch2, c.init_stretch2 / w.sf2)
                tmp = wk.xs_tmp
                fill!(tmp, 0.0)
                shift_vector!(tmp, w, c.vector, c.deriv2, shx, stx, st2x, cubic)
                for (k, i) in enumerate(pix)
                    Am[k, jc] = tmp[i]
                end
            elseif c.kind === :offl && c.poly_power == 0
                if w.spec.offl_ref                     # LINEAR_OFFSET_REF: ymean / I0
                    for (k, i) in enumerate(pix)
                        Am[k, jc] = abs(ref[i]) > 1.0e-14 ? ymean / ref[i] : 0.0
                    end
                else                                   # LINEAR_OFFSET_RAD: -xmean / I
                    for k in 1:npts
                        Am[k, jc] = abs(spec_nolog[k]) > 1.0e-14 ? -xmean / spec_nolog[k] : 0.0
                    end
                end
            else
                if c.poly_power == 0
                    for k in 1:npts
                        Am[k, jc] = 1.0          # vector = ones
                    end
                else
                    for (k, i) in enumerate(pix)
                        Am[k, jc] = Am[k, jc-1] * (w.lam[i] - w.lambda0)
                    end
                end
            end
        end
        # LINEAR_decompose (DECOMP_EIGEN_QR): unit-norm columns, pivoted QR,
        # variances from the Cholesky inverse of A'A.
        An = similar(Am)
        for j in 1:size(Am, 2)
            s = 0.0
            for k in 1:npts
                s += Am[k, j] * Am[k, j]
            end
            s == 0 && throw(EngineError("null column $(w.cols[j].name)"))
            wk.norms[j] = sqrt(s)
            for k in 1:npts
                An[k, j] = Am[k, j] / wk.norms[j]
            end
        end
        wk.qrf = qr(An, ColumnNorm())
        ch = cholesky(Symmetric(An' * An); check=false)
        issuccess(ch) || throw(EngineError("the fitted vectors of $(w.spec.name) are linearly dependent " *
                                           "(two cross sections alike, or a term that is zero in the window?)"))
        cov = inv(ch)
        for j in 1:size(Am, 2)
            wk.sigmasq[j] = cov[j, j] / (wk.norms[j] * wk.norms[j])
        end
        wk.decomposed = true
    end

    # reference shift/stretch
    rs = wk.ref_shifted
    fill!(rs, 0.0)
    rsh = param(A, w.ref_shift[1], w.ref_init[1])
    rst = param(A, w.ref_shift[2], w.ref_init[2] / w.sf1)
    rst2 = param(A, w.ref_shift[3], w.ref_init[3] / w.sf2)
    shift_vector!(rs, w, ref, spline_ref, rsh, rst, rst2, cubic)
    for i in w.limmin:w.limmax
        rs[i] <= 0 && throw(EngineError("log of non-positive reference at pixel $i"))
        rs[i] = log(rs[i])
    end
    Y = rs[pix]

    b = Y .- X
    c = (wk.qrf::QRT) \ b
    c ./= wk.norms

    Am = wk.Amat
    for (j, cj) in enumerate(c)
        for k in 1:npts
            X[k] += Am[k, j] * cj
        end
    end
    Yfit = Y .- X
    return Yfit, c
end

# =============================================================================
# curfit.c
# =============================================================================

function fchisq(Yfit, nfree)
    nfree <= 0 && return 0.0
    chisq = 0.0
    for v in Yfit
        d = 0.0 - v
        (abs(d) <= 1e16 && abs(d) >= 1e-16) && (chisq += d * d)
    end
    chisq / nfree
end

"CurfitMatinv: Gauss-Jordan inversion with full pivoting (in place). Returns det."
function curfit_matinv!(a::Matrix{Float64})
    n = size(a, 1)
    ik = zeros(Int, n); jk = zeros(Int, n)
    det = 1.0
    for k in 1:n
        amax = 0.0
        while true
            amax = 0.0
            for i in k:n, j in k:n
                if abs(amax) <= abs(a[i, j])
                    amax = a[i, j]; ik[k] = i; jk[k] = j
                end
            end
            amax == 0.0 && return 0.0
            i = ik[k]
            if i > k
                for j in 1:n
                    save = a[k, j]; a[k, j] = a[i, j]; a[i, j] = -save
                end
            end
            j = jk[k]
            if j > k && i >= k
                for ii in 1:n
                    save = a[ii, k]; a[ii, k] = a[ii, j]; a[ii, j] = -save
                end
            end
            (ik[k] < k || jk[k] < k) || break
        end
        for i in 1:n
            i != k && (a[i, k] /= -amax)
        end
        for i in 1:n
            i == k && continue
            for j in 1:n
                j != k && (a[i, j] += a[i, k] * a[k, j])
            end
        end
        for j in 1:n
            j != k && (a[k, j] /= amax)
        end
        a[k, k] = 1.0 / amax
        det *= amax
    end
    for l in 1:n
        k = n - l + 1
        j = ik[k]
        if j > k
            for i in 1:n
                save = a[i, k]; a[i, k] = -a[i, j]; a[i, j] = save
            end
        end
        i = jk[k]
        if i > k
            for jj in 1:n
                save = a[k, jj]; a[k, jj] = -a[i, jj]; a[i, jj] = save
            end
        end
    end
    det
end

"""
One call of QDOAS's Curfit (Bevington CURFIT). `st` carries Yfit, the linear
parameters and lambda between calls exactly as the C code does.
"""
function curfit!(st, w, wk, F, A, deltaA, niter_outer)
    nA = length(A)
    nY = length(w.pix)
    niter = niter_outer
    B = copy(A)
    if niter == 0
        st.Yfit, st.c = F(A)
    end
    # numerical derivatives (forward differences), one per parameter
    deriv = zeros(nA, nY)
    for j in 1:nA
        Aj = A[j]; Dj = deltaA[j]
        A[j] = Aj + Dj
        Yfit2, _ = F(A)
        Dj == 0 && throw(EngineError("zero delta"))
        for i in 1:nY
            deriv[j, i] = (Yfit2[i] - st.Yfit[i]) / Dj
        end
        A[j] = Aj
    end
    beta = zeros(nA)
    alpha = zeros(nA, nA)
    for i in 1:nY
        for j in 1:nA
            beta[j] += 1.0 * (0.0 - st.Yfit[i]) * deriv[j, i]
            for k in 1:j
                alpha[j, k] += 1.0 * deriv[j, i] * deriv[k, i]
            end
        end
    end
    for j in 1:nA, k in 1:j
        alpha[k, j] = alpha[j, k]
    end
    chisq1 = fchisq(st.Yfit, w.nfree)

    chisqr = 0.0
    arr = zeros(nA, nA)
    while true
        oldchisq = chisqr
        for j in 1:nA
            for k in 1:nA
                alpha[j, j] * alpha[k, k] <= 0 && throw(EngineError("Curfit sqrt argument ($(w.nl[j].name))"))
                arr[j, k] = alpha[j, k] / sqrt(alpha[j, j] * alpha[k, k])
            end
            arr[j, j] = 1.0 + st.lambda
        end
        curfit_matinv!(arr) == 0.0 && throw(EngineError("Curfit matrix inversion"))
        for j in 1:nA
            B[j] = A[j]
            for k in 1:nA
                B[j] += beta[k] * arr[j, k] / sqrt(alpha[j, j] * alpha[k, k])
            end
        end
        st.Yfit, st.c = F(B)
        chisqr = fchisq(st.Yfit, w.nfree)
        chisq1 < chisqr && (st.lambda *= 10.0)
        niter += 1
        niter > CURFIT_MAX_ITER && throw(EngineError("Curfit did not converge in $CURFIT_MAX_ITER iterations"))
        ((chisq1 < chisqr) && (chisqr != oldchisq)) || break
    end

    # range check (only where min != max; none of our parameters set limits)
    out = false
    for j in 1:nA
        q = w.nl[j]
        if q.minp != q.maxp
            lo, hi = min(q.minp, q.maxp), max(q.minp, q.maxp)
            (q.minp != 0 && q.minp == q.maxp) && (lo = -q.maxp)
            B[j] > hi && (B[j] = hi; out = true)
            B[j] < lo && (B[j] = lo; out = true)
        end
    end
    if out
        st.Yfit, st.c = F(B)
        chisqr = fchisq(st.Yfit, w.nfree)
    end

    for j in 1:nA
        A[j] = B[j]
        arr[j, j] / alpha[j, j] * chisqr <= 0 && throw(EngineError("Curfit sigma ($(w.nl[j].name))"))
        st.sigmaA[j] = sqrt(arr[j, j] / alpha[j, j] * chisqr)
    end
    st.lambda *= 0.1
    st.chisqr = chisqr
    return st
end

mutable struct FitState
    Yfit::Vector{Float64}
    c::Vector{Float64}
    lambda::Float64
    chisqr::Float64
    sigmaA::Vector{Float64}
end
