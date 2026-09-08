# Tutorial — how the γ coefficients relate to the knots

An independent, runnable answer to the question in
[`../figures/generalized_v2/TODO.md`](../figures/generalized_v2/TODO.md):

> why only 3 [dots] for K = 4 with 5 param, which of the 5 param are the controlling dots.
> What are the role of the 2 others, can they be displayed too?

```bash
Rscript tutorial/export_basis.R          # once — exports the real basis from R
python3 tutorial/spline_explorer.py      # sliders, one per coefficient
python3 tutorial/spline_explorer.py --save -K 5     # static PNG instead
```

Needs `numpy` and `matplotlib` only. The basis is **exported from R**, not
reimplemented, so what you manipulate is this project's actual model.

---

## The short answer

**None of the coefficients is a dot.** Knots and coefficients are different kinds of
object, and there is no one-to-one map between them:

| | what it is | chosen how |
| --- | --- | --- |
| **knot** | an x-position where two cubic pieces join | **fixed before fitting** — structure |
| **basis function** $B_k(x)$ | a fixed curve, derived from the knots | **determined** by the knots |
| **coefficient** $\gamma_k$ | a weight multiplying $B_k$ | **estimated** from the data |

The chain runs one way: **knots → basis → (weighted by γ) → curve.** A knot is not a
parameter; a parameter is not attached to a knot.

## Why 3 dots but 5 parameters

Because the plates only draw the **interior** knots as dots. There are two more:

```
K = 4  ->  3 interior knots  (50, 90, 120)      drawn dashed, with dots on the curve
       +   2 boundary knots  (0, 160)           drawn dotted, no dots
       =   5 knots in total

           4 basis functions B1..B4
           1 level term  γ0   +   4 weights γ1..γ4   =   5 parameters
```

So the **two "missing" ones are the boundary knots**, and they are already on the figures —
as the dotted verticals at 0° and 160°. They are where the natural-spline constraint takes
effect: beyond them the curve is forced to be a straight line.

That gives a tidy identity, true at every df:

$$\#\text{parameters} \;=\; \#\text{interior knots} + 2 \;=\; \#\text{knots in total}$$

| df (K) | interior | boundary | total knots | basis functions | parameters |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 2 | 1 | 2 | 3 | 2 | 3 |
| 3 | 2 | 2 | 4 | 3 | 4 |
| **4** | **3** | **2** | **5** | **4** | **5** |
| 5 | 4 | 2 | 6 | 5 | 6 |

**It is a dimension count, not a correspondence.** The spline space over 5 knots happens to
be 5-dimensional; that does not pair parameter 1 with knot 1.

## The part that surprises people

You might still hope the coefficients are at least *ordered* along x — γ₁ governing the
left, γ₄ the right. They are not. Measured on the real K = 4 basis:

| basis function | peaks at |
| --- | ---: |
| $B_1$ | 86.5° |
| $B_2$ | 121.3° |
| **$B_3$** | **45.7°** ← leftmost |
| $B_4$ | 160.0° |

**$B_3$ peaks to the left of $B_1$.** `ns()` returns an orthogonalised basis chosen for
numerical stability, not a left-to-right sequence of bumps. Panel 1 of the explorer shows
this at a glance.

Two consequences follow, and they are why the analysis is set up the way it is:

- **No single γ has a location.** Each $B_k$ is nonzero over most of the range, so changing
  one coefficient changes the curve nearly everywhere. Panel 2 makes this visible: move one
  slider and one coloured contribution changes shape across the whole plot.
- **Only the sum is interpretable.** $\gamma_0 + \sum_k \gamma_k B_k(x)$ is the difference
  curve; the individual γ are bookkeeping. This is exactly why
  [`../figures/generalized_v2/`](../figures/generalized_v2/) reports a **joint** χ² test on
  the whole γ block rather than leaning on the per-coefficient Wald strip — the joint
  statistic $\gamma^\top V^{-1}\gamma$ is invariant to how the basis is parametrised, while
  individual Wald tests move around as the knots move.

## What the three panels show

1. **The basis** $B_k(x)$, with interior knots dashed and boundary knots dotted. Move any
   slider: **these never change.** They are fixed the moment you choose the knots.
2. **The weighted pieces** $\gamma_0$ (dashed, flat) and each $\gamma_k B_k(x)$. One slider
   moves one colour, and it moves across a wide span of x, not at a knot.
3. **The sum** — the difference curve itself, with the fitted values as a dotted reference
   and the mean |difference| printed live.

Things worth trying:

- **`all zero`** then raise $\gamma_0$ alone: a flat offset, the same at every elevation.
  That is the "level shift" term, and it is the only one with a location-free meaning.
- Set $\gamma_0$ to zero and raise one $\gamma_k$: the curve bends, but it also moves at
  elevations far from any knot.
- Compare `-K 2` with `-K 5`: with two basis functions the curve can barely bend; with five
  it can follow an S. That is the flexibility the df sweep in
  [`../figures/04_natural_spline_hlme/`](../figures/04_natural_spline_hlme/) is trading
  against residual behaviour.

## Files

| file | what it is |
| --- | --- |
| `export_basis.R` | run once — writes the real basis, knots and fitted γ out of R |
| `spline_explorer.py` | the interactive figure |
| `basis_K<k>.csv` | basis functions on a 600-point grid, one file per df |
| `knots_K<k>.csv` | the interior knots for that df |
| `coefs.csv` | fitted difference coefficients from the worked example, so sliders start at the real fit |
| `spline_explorer_K4.png` | a static render, for when there is no display |
