# Generalized sweep — results index

Ran over **72** combinations (joint x movement x DoF) of angular data.
- **compare** (in-vivo vs ex-vivo tested): 33
- **single** (one condition only): 30
- **skipped** (too few units): 9
- **errors**: 0
- of the compared, **32** reach BH-FDR p < 0.05 on the shape test (EXPLORATORY).

> **These contrasts are EXPLORATORY and NOT causal.** Condition is perfectly
> confounded with source study (no study measured both), so a difference mixes
> muscle activation with protocol / measurement / population. The real
> replication is the number of **studies** (`st ex/in`), not rows or shoulders,
> and p-values are pseudo-replicated — read the **effect size in degrees** and
> the study counts, not the stars. p-values are BH-FDR adjusted across the
> compared family. See `STATISTICAL_METHODS_REVIEW.md` / `_RESPONSE.md`.

**Planches** (`planche_<motion>.png`): one plate per movement, rows = joints
(sternoclavicular -> acromioclavicular -> scapulothoracic -> glenohumeral),
cols = DoF 1/2/3. See `00_master_summary.csv` for every combo and each combo's
folder for its drill-down curves.

## Compared combinations (sorted by |effect size|, largest first)

| joint | movement | DoF | shoulders ex/in | **studies ex/in** | signed Δ° | \|max\|° | level shift° | level p (FDR) | shape p (FDR) |
|---|---|---|---|---|---|---|---|---|---|
| glenohumeral | internal-external rotation 0 degree-abducted | 3 | 10/21 | **1/2** | +70.6 | 75.8 | +48.01 | 6e-05 | 0 |
| glenohumeral | internal-external rotation 0 degree-abducted | 1 | 10/21 | **1/2** | +63.5 | 70.7 | +42.85 | 1.6e-05 | 0 |
| glenohumeral | frontal plane elevation | 3 | 10/23 | **1/3** | -60.9 | 74.9 | -43.42 | 8.3e-10 | 0 |
| glenohumeral | sagittal plane elevation | 3 | 10/24 | **1/4** | -48.8 | 105.6 | -35.05 | 1.9e-13 | 0 |
| acromioclavicular | sagittal plane elevation | 2 | 10/3 | **1/2** | -21.7 | 29.9 | -15.52 | 0.00034 | 0 |
| scapulothoracic | horizontal flexion | 2 | 10/7 | **2/1** | +19.8 | 21.3 | +14.15 | 3.7e-11 | 0 |
| glenohumeral | internal-external rotation 0 degree-abducted | 2 | 10/21 | **1/2** | +18.6 | 21.5 | +12.94 | 8.8e-09 | 0 |
| acromioclavicular | frontal plane elevation | 2 | 10/4 | **1/3** | -18.0 | 20.7 | -12.70 | 0.0027 | 0 |
| glenohumeral | sagittal plane elevation | 1 | 10/24 | **1/4** | +12.5 | 75.8 | +11.15 | 0.0047 | 0 |
| scapulothoracic | horizontal flexion | 1 | 10/7 | **2/1** | -12.9 | 24.0 | -10.39 | 0.00013 | 0 |
| sternoclavicular | frontal plane elevation | 3 | 13/4 | **4/3** | +12.8 | 18.8 | +8.74 | 0.14 | 0 |
| scapulothoracic | horizontal flexion | 3 | 10/7 | **2/1** | -12.6 | 15.7 | -10.54 | 9.3e-09 | 0 |
| sternoclavicular | sagittal plane elevation | 3 | 13/3 | **4/2** | +10.9 | 17.9 | +9.20 | 0.16 | 0 |
| scapulothoracic | frontal plane elevation | 2 | 13/31 | **4/5** | -10.4 | 15.1 | -6.42 | 6e-05 | 0 |
| scapulothoracic | internal-external rotation 0 degree-abducted | 2 | 10/21 | **1/2** | -10.4 | 10.7 | -7.24 | 0.00012 | 0.13 |
| glenohumeral | frontal plane elevation | 2 | 10/23 | **1/3** | +10.0 | 15.3 | +5.65 | 0.005 | 0 |
| acromioclavicular | frontal plane elevation | 3 | 10/4 | **1/3** | -9.1 | 10.3 | -6.47 | 0.046 | 0 |
| glenohumeral | frontal plane elevation | 1 | 10/23 | **1/3** | +0.8 | 20.9 | +1.47 | 0.81 | 0 |
| scapulothoracic | sagittal plane elevation | 2 | 13/25 | **4/5** | -8.6 | 14.2 | -5.35 | 0.014 | 0 |
| acromioclavicular | sagittal plane elevation | 1 | 10/3 | **1/2** | +7.7 | 8.3 | +7.99 | 0.0095 | 0 |
| scapulothoracic | frontal plane elevation | 1 | 13/31 | **4/5** | -7.1 | 9.5 | -4.99 | 0.083 | 0 |
| scapulothoracic | internal-external rotation 0 degree-abducted | 1 | 10/21 | **1/2** | -5.8 | 14.5 | -3.22 | 0.12 | 0 |
| sternoclavicular | sagittal plane elevation | 2 | 13/3 | **4/2** | -6.6 | 8.6 | +28.70 | 8.3e-10 | 0 |
| scapulothoracic | frontal plane elevation | 3 | 13/31 | **4/5** | +6.3 | 10.7 | +4.92 | 0.084 | 0 |
| glenohumeral | sagittal plane elevation | 2 | 10/24 | **1/4** | +4.4 | 9.1 | +3.10 | 0.11 | 0 |
| scapulothoracic | internal-external rotation 0 degree-abducted | 3 | 10/21 | **1/2** | -5.9 | 6.3 | -3.97 | 0.051 | 1.8e-05 |
| sternoclavicular | sagittal plane elevation | 1 | 13/3 | **4/2** | +5.1 | 7.6 | +3.15 | 0.45 | 0 |
| scapulothoracic | sagittal plane elevation | 1 | 13/25 | **4/5** | +5.1 | 9.3 | +3.94 | 0.16 | 0 |
| sternoclavicular | frontal plane elevation | 2 | 13/4 | **4/3** | -2.5 | 7.5 | -2.55 | 0.43 | 0 |
| acromioclavicular | sagittal plane elevation | 3 | 10/3 | **1/2** | -4.1 | 6.6 | -6.28 | 0.023 | 0 |
| acromioclavicular | frontal plane elevation | 1 | 10/4 | **1/3** | +3.5 | 6.6 | +2.44 | 0.41 | 0 |
| scapulothoracic | sagittal plane elevation | 3 | 13/25 | **4/5** | +0.5 | 6.5 | +0.83 | 0.61 | 0 |
| sternoclavicular | frontal plane elevation | 1 | 13/4 | **4/3** | -2.3 | 3.4 | -1.70 | 0.54 | 0 |
