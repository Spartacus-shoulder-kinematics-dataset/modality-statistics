The forest graph is a begining. There is so much to perfect. But the idea is here ! :) - Significant map is a step two

[DONE] -> regroup by dof, elevation frontal -> scapular -> sagittal, POC for scapulothoracic
dof2, then all 12 joint x DoF groups. `00_sigmap_poc_scapulothoracic_dof2` and
`00_sigmap_elevation`. Rows are joint x DoF, motions are the inner axis, one shared abscissa
so the groups are comparable end to end. The header gives NAME / DoF / motion as three
centred lines and is written once per group, never repeated per row. Second axes = mean
difference bar with its 95% CI, stacked close. Map shows the 95% significant area over the
fitted overlap, with ex-vivo and in-vivo coverage in their own colours above and below.
Far-right column = FDR stars + the gamma dot strip with its index.

[DONE] -> rotation plate `00_sigmap_rotation`, IER 0° and 90°.

**But read the README first — the two motions asked for are largely empty, and that is the
finding, not a bug:**
- **scapular plane elevation has NO compared cell anywhere in the sweep** (<= 2 ex-vivo
  shoulders everywhere against MIN_PER_COND = 3; scapulothoracic is 2 ex from 2 studies vs
  49 in). Every scapular row is a labelled placeholder.
- **IER 90° has NO ex-vivo shoulder in any joint.** Only IER 0° compares, and only at
  scapulothoracic and glenohumeral.

-> STILL OPEN on the forest: the right-hand columns overprint each other (mean |Δ|,
elevation range, Wald and joint collide into unreadable text at the header and every row)
and the motion labels are clipped off the left edge. Deliberately left alone this pass.



[DONE in /tutorial] -> still need to figure out the link between the gammas and the dot controlling the spline. why only 3 for K= 4 with 5 param which of the 5 param are the controlling dots. What are the role of the 2 others, can they be displayed too ? 

-> Find a way to store the results, to avoid regenerating the results every time. Decoupling figures from result computations, in generalized v2



