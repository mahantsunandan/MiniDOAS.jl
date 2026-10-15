# This file is part of MiniDOAS.jl (Differential Optical Absorption Spectroscopy in Julia).
# Copyright (c) 2026 Sunandan Mahant <sunandanmahant@outlook.com>
# Derived from QDOAS, Copyright (C) 1994-2025 BIRA-IASB and S[&]T (BSD-3-Clause).
# See LICENSE and NOTICE.md.

# =============================================================================
# spline.c
# =============================================================================

"SPLINE_Deriv2: natural cubic spline second derivatives (NR, QDOAS formulation)."
function spline_deriv2(X::AbstractVector{Float64}, Y::AbstractVector{Float64})
    n = length(X)
    Y2 = zeros(n)
    u = zeros(n)
    Y2[1] = u[1] = 0.0
    dx_old = X[2] - X[1]
    dydx_old = (Y[2] - Y[1]) / dx_old
    i = 2
    while i <= n - 1
        X[i+1] - X[i] <= 0 && throw(EngineError("non increasing abscissae"))
        d2x = 1.0 / (X[i+1] - X[i-1])
        sig = dx_old * d2x
        dp = 1.0 / (sig * Y2[i-1] + 2.0)
        Y2[i] = (sig - 1.0) * dp
        dx = X[i+1] - X[i]
        dydx = (Y[i+1] - Y[i]) / dx
        ui = dydx - dydx_old
        ui = (6.0 * ui * d2x - sig * u[i-1]) * dp
        abs(ui) < 1.0e-300 && (ui = 0.0)
        u[i] = ui
        dx_old = dx
        dydx_old = dydx
        i += 1
    end
    # upper boundary (natural): qn = un = 0
    Y2[i] = (0.0 - 0.0 * u[i-1]) / (0.0 * Y2[i-1] + 1.0)
    for k in i-1:-1:1
        Y2[k] = Y2[k] * Y2[k+1] + u[k]
    end
    Y2
end

"SPLINE_Vector for a single abscissa."
@inline function spline_at(xa, ya, y2a, x::Float64, cubic::Bool)
    na = length(xa)
    x <= xa[1] && return ya[1]
    x >= xa[na] && return ya[na]
    k = searchsortedfirst(xa, x) - 1        # largest k with xa[k] < x
    xhi = xa[k+1]
    h = xhi - xa[k]
    a = (xhi - x) / h
    b = 1.0 - a
    y = a * ya[k] + b * ya[k+1]
    if cubic
        y += ((a * a * a - a) * y2a[k] + (b * b * b - b) * y2a[k+1]) * (h * h) / 6.0
    end
    y
end

"VECTOR_NormalizeVector: v / ||v||, returns the norm."
function normalize!(v)
    nsq = sum(abs2, v)
    nsq == 0 && throw(EngineError("cannot normalise a null vector"))
    nrm = sqrt(nsq)
    v ./= nrm
    nrm
end

"FNPixel (0-based in QDOAS; returned 1-based here)."
function fnpixel(lam, value, sel::Symbol)
    n = length(lam)
    value <= lam[1] && return 1
    value >= lam[n] && return n
    rc = (n - 1) >> 1                      # 0-based
    klo, khi = 0, n - 1
    while khi - klo > 1
        rc = (khi + klo) >> 1
        abs(lam[rc+1] - value) < EPSILON && break
        if lam[rc+1] > value
            khi = rc
        else
            klo = rc
        end
    end
    if abs(lam[rc+1] - value) > EPSILON
        if sel === :before
            (rc > 0 && lam[rc+1] > value) && (rc -= 1)
        else
            (rc < n - 1 && lam[rc+1] < value) && (rc += 1)
        end
    end
    rc + 1
end
