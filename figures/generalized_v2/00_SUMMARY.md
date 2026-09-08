# Generalized sweep v2 — random spline + hlme

[Iteration 04](../04_natural_spline_hlme/)'s model applied to all 72
joint x motion x DoF combinations. Produced by
[`../../analyze_all_v2.R`](../../analyze_all_v2.R).

## The model

```r
hlme(fixed   = Y ~ ns(TIME, knots) * cond,      # knots: 3 interior, per cell
     random  = ~ ns(TIME, knots),               # a random CURVE per shoulder
     subject = "IDnum", ng = 1,
     idiag   = TRUE)                            # D diagonal
```

For a shoulder $i$, observation $j$, in condition $c_i \in \{\text{ex vivo},\ \text{in vivo}\}$,
at elevation $x_{ij}$:

$$
\begin{aligned}
y_{ij} =\;& \underbrace{\beta_0 + \sum_{k=1}^{4}\beta_k B_k(x_{ij})}_{\text{ex-vivo reference}}
 + \underbrace{\Big[\gamma_0 + \sum_{k=1}^{4}\gamma_k B_k(x_{ij})\Big]\mathbb{1}(c_i=\text{in vivo})}_{\text{in-vivo departure}} \\[4pt]
&+ \underbrace{b_{i0} + \sum_{k=1}^{4} b_{ik} B_k(x_{ij})}_{\text{this shoulder's own curve}} + \varepsilon_{ij},
\qquad \mathbf{b}_i \sim \mathcal{N}_{5}(\mathbf{0}, D),
\quad \varepsilon_{ij}\sim\mathcal{N}(0,\sigma_\varepsilon^2).
\end{aligned}
$$

$D$ is **diagonal** (`idiag = TRUE`): 5 variances, no covariances. Per cell that is
**10 fixed effects** (intercept, 4 basis terms, the level shift, 4 interactions),
5 random-effect variances and 1 residual sd.

$\mathbb{1}(c_i = \text{in vivo})$ is 1 on in-vivo rows and 0 otherwise, so
$\gamma_0$ shifts the curve and $\gamma_1 \dots \gamma_{4}$ bend it; all $\gamma = 0$ would mean the
conditions are indistinguishable. That is exactly the **joint test** tabulated below.

Cells with only one condition are fitted without the $\gamma$ block and reported
descriptively.

**Per cell:** thin to <=100 points/shoulder, trim to the x-overlap, build the basis from
that cell's own quantiles with boundary knots from its untrimmed range. Curves are drawn
only where at least 2 shoulders have data.

## Results

- **compare** (both conditions): 33
- **single** (one condition):    30
- **skipped** (too few units):   9
- **did not converge**:          0
- of the compared, **25** reach BH-FDR p < 0.05 on the JOINT difference test

> **Exploratory, not causal.** Condition is perfectly confounded with source study
> (no study measured both), so a difference mixes muscle activation with protocol,
> measurement and population. The real replication is the number of **studies**
> (`st ex/in`), not rows or shoulders. Read the effect size in degrees.

> The **joint** test (all difference coefficients at once) is the one to read, not the
> per-coefficient Wald counts: the basis coefficients are correlated, so which one
> carries the signal shifts with the knots while the joint statistic does not.

## Reading the planches

Each motion has two plates — `planche_<motion>.png/.pdf` (population curves) and
`planche_diff_<motion>.png/.pdf` (the in-vivo minus ex-vivo difference). Rows are the
four joints proximal to distal, columns the three rotational DoF, and **all twelve
panels share one abscissa** so cells are directly comparable.

| glyph | meaning |
| --- | --- |
| **●** (filled dot, top right) | that in-vivo difference coefficient passed \|coef/se\| >= 1.96 |
| **○** (hollow dot, top right) | that coefficient did **not** pass |
| **\*\*\*  \*\*  \*  ns** (top right) | the **joint** chi2 test on *all* difference coefficients at once, BH-FDR adjusted: `***` p<0.001, `**` p<0.01, `*` p<0.05, `ns` not significant |
| small dot **on a curve** | the fitted value at an interior knot |
| **dashed vertical** | an interior knot — where the piecewise cubics join |
| orange / green | ex vivo / in vivo; the band is the pointwise 95% CI |
| blue band, difference plate | pointwise 95% CI of the difference; **squares along the bottom** mark where it excludes zero |
| bottom-left text | shoulders ex/in, studies ex/in, Wald counts ref/diff, mean and max \|delta\| in degrees |

### Which dot is which coefficient

The strip has **5 dots, read LEFT TO RIGHT**; on the figure the index is printed
directly under each dot. The order is fixed:

| position | symbol | model term | meaning |
| ---: | :---: | :--- | :--- |
| **1st** dot | $\gamma_0$ | `condin vivo` | the **constant level shift** between conditions |
| **2nd** dot | $\gamma_1$ | `ns1:condin vivo` | shape term on basis function $B_1$ |
| **3rd** dot | $\gamma_2$ | `ns2:condin vivo` | shape term on basis function $B_2$ |
| **4th** dot | $\gamma_3$ | `ns3:condin vivo` | shape term on basis function $B_3$ |
| **5th** dot | $\gamma_4$ | `ns4:condin vivo` | shape term on basis function $B_4$ |

So `○●○○●` reads: level shift **failed**, $\gamma_1$ **passed**, $\gamma_2$ failed,
$\gamma_3$ failed, $\gamma_4$ **passed**.

Every strip can be checked against **`00_wald_tests.csv`**, which has one row per
coefficient with its `term`, `wald` and `passed`.

The ● / ○ strip and the stars answer **different questions**, and they routinely
disagree: the strip says which individual coefficients are separable from zero, the
stars say whether the difference is nonzero *as a whole*. Because the basis
coefficients are correlated, a cell can show `○○○○○` and still be strongly
significant jointly. **Trust the stars, use the strip only as texture.**

| joint | motion | DoF | sh ex/in | st ex/in | mean abs diff | joint chi2 | joint p (FDR) | shape p (FDR) | Wald ref | Wald diff | symbols |
| --- | --- | ---: | :---: | :---: | ---: | ---: | ---: | ---: | :---: | :---: | :--- |
| glenohumeral | internal-external rotation 0 degree-abducted | 3 | 10/21 | 1/2 | 80.3° | 5.4 | 0.4 ns | 0.99 | 5/5 | 1/5 | `●○○○○` |
| glenohumeral | frontal plane elevation | 3 | 10/23 | 1/3 | 60.5° | 33.4 | 6.5e-06 *** | 0.028 | 4/5 | 2/5 | `●○●○○` |
| glenohumeral | sagittal plane elevation | 3 | 10/24 | 1/4 | 54.9° | 274.4 | 1.1e-55 *** | 4.8e-42 | 5/5 | 5/5 | `●●●●●` |
| glenohumeral | internal-external rotation 0 degree-abducted | 1 | 10/21 | 1/2 | 48.9° | 3.7 | 0.62 ns | 0.68 | 0/5 | 0/5 | `○○○○○` |
| acromioclavicular | sagittal plane elevation | 2 | 10/3 | 1/2 | 20.9° | 167.1 | 2.5e-33 *** | 2.6e-33 | 4/5 | 4/5 | `○●●●●` |
| scapulothoracic | horizontal flexion | 2 | 10/7 | 2/1 | 20.8° | 47.9 | 9.3e-09 *** | 4.5e-07 | 4/5 | 3/5 | `●●○○●` |
| glenohumeral | internal-external rotation 0 degree-abducted | 2 | 10/21 | 1/2 | 18.8° | 14.6 | 0.017 * | 0.32 | 3/5 | 1/5 | `●○○○○` |
| acromioclavicular | frontal plane elevation | 2 | 10/4 | 1/3 | 18.3° | 61.7 | 2e-11 *** | 2.1e-11 | 3/5 | 4/5 | `○●●●●` |
| glenohumeral | sagittal plane elevation | 1 | 10/24 | 1/4 | 17.1° | 265.3 | 4.8e-54 *** | 3.4e-42 | 5/5 | 5/5 | `●●●●●` |
| scapulothoracic | horizontal flexion | 1 | 10/7 | 2/1 | 12.1° | 36.8 | 1.5e-06 *** | 5.2e-07 | 5/5 | 3/5 | `○●●●○` |
| scapulothoracic | horizontal flexion | 3 | 10/7 | 2/1 | 12.0° | 23.7 | 0.00046 *** | 0.28 | 4/5 | 1/5 | `●○○○○` |
| sternoclavicular | sagittal plane elevation | 3 | 13/3 | 4/2 | 11.8° | 26.7 | 0.00013 *** | 6.1e-05 | 4/5 | 2/5 | `○●●○○` |
| sternoclavicular | frontal plane elevation | 3 | 13/4 | 4/3 | 11.7° | 55.0 | 3.6e-10 *** | 1.2e-10 | 4/5 | 3/5 | `○●●●○` |
| scapulothoracic | frontal plane elevation | 2 | 13/31 | 4/5 | 10.9° | 93.0 | 8.7e-18 *** | 2.1e-18 | 4/5 | 2/5 | `○●○○●` |
| glenohumeral | frontal plane elevation | 2 | 10/23 | 1/3 | 10.8° | 77.6 | 1.3e-14 *** | 3.4e-15 | 4/5 | 3/5 | `○●●○●` |
| scapulothoracic | internal-external rotation 0 degree-abducted | 2 | 10/21 | 1/2 | 10.4° | 7.6 | 0.22 ns | 0.89 | 2/5 | 1/5 | `●○○○○` |
| scapulothoracic | internal-external rotation 0 degree-abducted | 1 | 10/21 | 1/2 | 9.7° | 55.7 | 2.8e-10 *** | 1.1e-10 | 4/5 | 5/5 | `●●●●●` |
| scapulothoracic | sagittal plane elevation | 2 | 13/25 | 4/5 | 8.9° | 106.7 | 1.3e-20 *** | 2.5e-21 | 4/5 | 2/5 | `○●○●○` |
| glenohumeral | frontal plane elevation | 1 | 10/23 | 1/3 | 8.7° | 22.6 | 0.00069 *** | 0.0011 | 2/5 | 2/5 | `○●●○○` |
| acromioclavicular | frontal plane elevation | 3 | 10/4 | 1/3 | 7.7° | 6.9 | 0.26 ns | 0.61 | 5/5 | 0/5 | `○○○○○` |
| scapulothoracic | frontal plane elevation | 1 | 13/31 | 4/5 | 7.4° | 17.8 | 0.005 ** | 0.037 | 2/5 | 2/5 | `●○●○○` |
| acromioclavicular | sagittal plane elevation | 1 | 10/3 | 1/2 | 7.3° | 12.3 | 0.04 * | 0.073 | 4/5 | 0/5 | `○○○○○` |
| glenohumeral | sagittal plane elevation | 2 | 10/24 | 1/4 | 7.1° | 258.3 | 9.8e-53 *** | 6.5e-43 | 4/5 | 4/5 | `●●●●○` |
| scapulothoracic | internal-external rotation 0 degree-abducted | 3 | 10/21 | 1/2 | 7.0° | 7.3 | 0.24 ns | 0.4 | 1/5 | 0/5 | `○○○○○` |
| scapulothoracic | sagittal plane elevation | 1 | 13/25 | 4/5 | 6.9° | 21.9 | 0.0009 *** | 0.00059 | 3/5 | 3/5 | `○○●●●` |
| sternoclavicular | sagittal plane elevation | 2 | 13/3 | 4/2 | 6.5° | 57.2 | 1.5e-10 *** | 8.8e-11 | 5/5 | 3/5 | `○●○●●` |
| scapulothoracic | frontal plane elevation | 3 | 13/31 | 4/5 | 6.3° | 36.3 | 1.8e-06 *** | 5.7e-07 | 5/5 | 3/5 | `○●●●○` |
| sternoclavicular | frontal plane elevation | 2 | 13/4 | 4/3 | 6.2° | 17.7 | 0.0051 ** | 0.0029 | 5/5 | 1/5 | `○●○○○` |
| sternoclavicular | frontal plane elevation | 1 | 13/4 | 4/3 | 3.9° | 2.5 | 0.78 ns | 0.89 | 5/5 | 0/5 | `○○○○○` |
| sternoclavicular | sagittal plane elevation | 1 | 13/3 | 4/2 | 3.8° | 3.7 | 0.62 ns | 0.61 | 5/5 | 0/5 | `○○○○○` |
| acromioclavicular | sagittal plane elevation | 3 | 10/3 | 1/2 | 3.5° | 13.1 | 0.03 * | 0.034 | 5/5 | 1/5 | `○●○○○` |
| scapulothoracic | sagittal plane elevation | 3 | 13/25 | 4/5 | 2.3° | 76.0 | 2.4e-14 *** | 4.3e-13 | 5/5 | 5/5 | `●●●●●` |
| acromioclavicular | frontal plane elevation | 1 | 10/4 | 1/3 | 2.2° | 9.8 | 0.1 ns | 0.073 | 3/5 | 0/5 | `○○○○○` |
