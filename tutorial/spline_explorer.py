#!/usr/bin/env python3
"""
Interactive explorer: what do the knots and the coefficients each actually do?

    python3 tutorial/spline_explorer.py            # sliders (needs a display)
    python3 tutorial/spline_explorer.py -K 5       # a different spline df
    python3 tutorial/spline_explorer.py --save     # write a static PNG instead

Two sets of sliders, and the difference between them is the whole lesson:

    KNOT sliders   move the joins -> the BASIS is rebuilt -> every panel changes
    GAMMA sliders  reweight a fixed basis -> panel 1 never moves

The natural-spline basis is computed here in Python and is verified against
R's splines::ns() to ~1e-14 (see `--verify`), so this is the same basis the
models in figures/ actually use — not a look-alike.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import BSpline

HERE = Path(__file__).resolve().parent
BOUNDARY = (0.0, 160.0)
BLUE = "#377eb8"
PALETTE = ["#e41a1c", "#4daf4a", "#984ea3", "#ff7f00", "#a65628", "#17becf"]
KNOT_COL = "#555555"

# knot rule of analyze_all_v2.R: quantiles of x, rounded to the nearest 10
DEFAULT_KNOTS = {2: [90.0], 3: [60.0, 110.0], 4: [50.0, 90.0, 120.0],
                 5: [40.0, 70.0, 100.0, 130.0]}


# --------------------------------------------------------------------------
def ns_basis(x, knots, boundary=BOUNDARY):
    """R's splines::ns(x, knots, Boundary.knots, intercept = FALSE).

    Cubic B-spline design matrix, then project out the two natural-boundary
    constraints (second derivative zero at each boundary knot) by QR — the
    same construction R uses, so the columns match to machine precision.
    """
    x = np.asarray(x, float)
    knots = np.sort(np.asarray(knots, float))
    lo, hi = float(boundary[0]), float(boundary[1])
    t = np.concatenate([[lo] * 4, knots, [hi] * 4])
    n = len(t) - 4

    def design(pts, der=0):
        M = np.empty((len(pts), n))
        for i in range(n):
            c = np.zeros(n)
            c[i] = 1.0
            sp = BSpline(t, c, 3, extrapolate=True)
            M[:, i] = (sp.derivative(der) if der else sp)(pts)
        return M

    B = design(x)
    const = design(np.array([lo, hi]), der=2)
    B, const = B[:, 1:], const[:, 1:]          # intercept = FALSE
    Q, _ = np.linalg.qr(const.T, mode="complete")
    return B @ Q[:, 2:]


def fitted_gammas(K):
    """Fitted difference coefficients from the worked example, if exported."""
    g = np.zeros(K + 1)
    cf = HERE / "coefs.csv"
    if not cf.exists():
        return g
    d = np.atleast_1d(np.genfromtxt(cf, delimiter=",", names=True,
                                    dtype=None, encoding="utf-8"))
    for row in d:
        if int(row["K"]) != K or np.isnan(row["estimate"]):
            continue
        j = row["gamma"]
        if not np.isnan(j) and 0 <= int(j) <= K:
            g[int(j)] = float(row["estimate"])
    return g


def verify():
    """Check the Python basis against R's exported CSVs, if present."""
    ok = True
    for K in (2, 3, 4, 5):
        f = HERE / f"basis_K{K}.csv"
        if not f.exists():
            print(f"  K={K}: no R export to compare against — run export_basis.R")
            continue
        r = np.genfromtxt(f, delimiter=",", names=True)
        R = np.column_stack([r[f"B{k}"] for k in range(1, K + 1)])
        P = ns_basis(r["x"], DEFAULT_KNOTS[K])
        err = np.abs(P - R).max()
        ok &= err < 1e-10
        print(f"  K={K}: max |python − R| = {err:.2e}")
    print("  verified" if ok else "  MISMATCH")
    return 0 if ok else 1


# --------------------------------------------------------------------------
def build(K, interactive=True):
    x = np.linspace(BOUNDARY[0], BOUNDARY[1], 600)
    kn = np.array(DEFAULT_KNOTS[K], float)
    kn0 = kn.copy()
    g0 = fitted_gammas(K)
    g = g0.copy()
    state = {"B": ns_basis(x, kn)}

    n_rows = max(K + 1, len(kn)) + 1
    fig = plt.figure(figsize=(13.5, 10.5))
    bottom = 0.06 + 0.038 * n_rows
    gs = fig.add_gridspec(3, 1, height_ratios=[1, 1, 1.15], left=0.07, right=0.98,
                          top=0.905, bottom=bottom, hspace=0.42)
    ax1, ax2, ax3 = (fig.add_subplot(gs[i]) for i in range(3))
    title = fig.suptitle("", fontsize=12, fontweight="bold")

    basis_lines = [ax1.plot([], [], color=PALETTE[k % len(PALETTE)], lw=2)[0]
                   for k in range(K)]
    basis_tags = [ax1.text(0, 0, "", ha="center", fontsize=9, fontweight="bold",
                           color=PALETTE[k % len(PALETTE)]) for k in range(K)]
    piece_lines = [ax2.plot([], [], color=PALETTE[k % len(PALETTE)], lw=2)[0]
                   for k in range(K)]
    (level_line,) = ax2.plot([], [], color="0.35", lw=2, ls="--")
    (total_line,) = ax3.plot([], [], color=BLUE, lw=3)
    (dots,) = ax3.plot([], [], "o", ms=7, mfc=BLUE, mec="white", mew=1.3, zorder=5)
    info = ax3.text(0.01, 0.05, "", transform=ax3.transAxes, fontsize=9, color="0.25")

    kn_lines, bd_lines, kn_tags = [], [], []
    for a in (ax1, ax2, ax3):
        kn_lines.append([a.axvline(v, ls="--", lw=1.0, color=KNOT_COL, zorder=0)
                         for v in kn])
        bd_lines.append([a.axvline(v, ls=":", lw=1.3, color="0.2", zorder=0)
                         for v in BOUNDARY])
        a.axhline(0, color="0.75", lw=0.8)
        a.set_xlim(BOUNDARY[0] - 5, BOUNDARY[1] + 5)
    kn_tags = [ax1.text(v, 0, "", ha="center", fontsize=7, color="0.35") for v in kn]

    ax1.set_title("1. The basis $B_k(x)$ — rebuilt when a KNOT moves, untouched by the γ sliders",
                  fontsize=10, loc="left")
    ax1.set_ylabel("$B_k(x)$")
    ax2.set_title(r"2. Weighted pieces: $\gamma_0$ (dashed) and each $\gamma_k B_k(x)$",
                  fontsize=10, loc="left")
    ax2.set_ylabel("contribution (°)")
    ax3.set_title(r"3. The curve $\gamma_0+\sum_k\gamma_k B_k(x)$ — the only interpretable object",
                  fontsize=10, loc="left")
    ax3.set_xlabel("thoracohumeral elevation (°)")
    ax3.set_ylabel("in vivo − ex vivo (°)")

    def redraw(rebuild_basis=False):
        if rebuild_basis:
            state["B"] = ns_basis(x, kn)
        B = state["B"]
        for k in range(K):
            basis_lines[k].set_data(x, B[:, k])
            i = int(np.argmax(np.abs(B[:, k])))
            basis_tags[k].set_position((x[i], B[i, k]))
            basis_tags[k].set_text(f"$B_{k+1}$")
            piece_lines[k].set_data(x, g[k + 1] * B[:, k])
        level_line.set_data(x, np.full_like(x, g[0]))
        y = g[0] + B @ g[1:]
        total_line.set_data(x, y)
        dots.set_data(kn, np.interp(kn, x, y))
        for group in kn_lines:
            for ln, v in zip(group, kn):
                ln.set_xdata([v, v])
        lo1 = ax1.get_ylim()[0]
        for t_, v in zip(kn_tags, kn):
            t_.set_position((v, lo1)); t_.set_text(f"{v:.0f}")
        info.set_text(f"mean |difference| = {np.mean(np.abs(y)):.2f}°     "
                      f"range {y.min():+.1f}° … {y.max():+.1f}°")
        title.set_text(
            f"ns(df={K}):  {len(kn)} interior knots ("
            + ", ".join(f"{v:.0f}" for v in kn)
            + f") + 2 boundary = {len(kn)+2} knots   |   {K} basis functions   |   "
            f"{K+1} parameters")
        for a in (ax1, ax2, ax3):
            a.relim(); a.autoscale_view(scalex=False)
        fig.canvas.draw_idle()

    redraw(rebuild_basis=True)
    if not interactive:
        return fig, None

    from matplotlib.widgets import Slider, Button
    gsl, ksl = [], []

    for i in range(K + 1):
        axs = fig.add_axes([0.09, bottom - 0.055 - 0.038 * i, 0.36, 0.023])
        lab = r"$\gamma_0$ level" if i == 0 else rf"$\gamma_{i}$ on $B_{i}$"
        span = max(30.0, float(np.max(np.abs(g0))) * 2.5 if np.any(g0) else 40.0)
        s = Slider(axs, lab, -span, span, valinit=g[i], valfmt="%+.1f",
                   color=("0.5" if i == 0 else PALETTE[(i - 1) % len(PALETTE)]))

        def gcb(val, idx=i):
            g[idx] = val
            redraw(rebuild_basis=False)     # basis unchanged: that is the point

        s.on_changed(gcb)
        gsl.append(s)

    for j in range(len(kn)):
        axs = fig.add_axes([0.60, bottom - 0.055 - 0.038 * j, 0.30, 0.023])
        s = Slider(axs, f"knot {j+1}", BOUNDARY[0] + 5, BOUNDARY[1] - 5,
                   valinit=kn[j], valfmt="%.0f°", color=KNOT_COL)

        def kcb(val, idx=j):
            # keep knots strictly increasing and inside the boundary knots
            lo = BOUNDARY[0] + 2 if idx == 0 else kn[idx - 1] + 2
            hi = BOUNDARY[1] - 2 if idx == len(kn) - 1 else kn[idx + 1] - 2
            kn[idx] = float(np.clip(val, lo, hi))
            redraw(rebuild_basis=True)      # the basis is a function of the knots

        s.on_changed(kcb)
        ksl.append(s)

    def reset(_):
        for s in ksl:
            s.reset()
        kn[:] = kn0
        for s in gsl:
            s.reset()
        redraw(rebuild_basis=True)

    def zero(_):
        for s in gsl:
            s.set_val(0.0)

    b1 = Button(fig.add_axes([0.60, 0.022, 0.13, 0.032]), "reset all")
    b2 = Button(fig.add_axes([0.76, 0.022, 0.13, 0.032]), "γ = 0")
    b1.on_clicked(reset)
    b2.on_clicked(zero)
    fig.text(0.09, 0.03, "left: reweight a fixed basis      right: move the joins, "
                         "rebuilding the basis", fontsize=8.5, color="0.35")
    fig._widgets = (gsl, ksl, b1, b2)
    return fig, (gsl, ksl)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-K", type=int, default=4, choices=[2, 3, 4, 5],
                    help="natural-spline df (default 4, as in analyze_all_v2.R)")
    ap.add_argument("--save", metavar="PNG", nargs="?", const="auto",
                    help="write a static figure instead of opening sliders")
    ap.add_argument("--verify", action="store_true",
                    help="check the Python basis against R's exported CSVs and exit")
    a = ap.parse_args()

    if a.verify:
        sys.exit(verify())

    if a.save:
        matplotlib.use("Agg")
        fig, _ = build(a.K, interactive=False)
        out = HERE / (f"spline_explorer_K{a.K}.png" if a.save == "auto" else a.save)
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"wrote {out}")
        return

    try:
        build(a.K, interactive=True)
        plt.show()
    except Exception as exc:
        print(f"interactive mode failed ({exc}); falling back to --save", file=sys.stderr)
        matplotlib.use("Agg")
        fig, _ = build(a.K, interactive=False)
        out = HERE / f"spline_explorer_K{a.K}.png"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
