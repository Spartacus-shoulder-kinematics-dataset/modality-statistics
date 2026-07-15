#!/usr/bin/env python3
"""Export tidy long-format CSVs from the raw Spartacus data.

Writes two files:

1. monolix_st_frontal_dof2.csv  -- the single feature-example subset
   (scapulothoracic / frontal plane elevation / DoF=2 / rad), for MonolixSuite
   (which works on one dataset at a time). Columns: ID, TIME, Y, in_vivo.

2. spartacus_angles_long.csv    -- ALL angular data (unit == rad), every
   joint x humeral_motion x degree_of_freedom, for the generalized R sweep
   (analyze_all.R). Columns: joint, humeral_motion, degree_of_freedom,
   ID, TIME, Y, in_vivo.

Column roles (both files):
    ID       -> subject/unit id   (article + "_" + shoulder_id)
    TIME     -> regressor (x)      (humerothoracic_angle)
    Y        -> observation (y)    (value, the joint angle)
    in_vivo  -> categorical        (True/False)

Usage:  python3 prepare_monolix_data.py
"""

import csv
import sys

SRC          = "corrected_confident_data.csv"
DST_SUBSET   = "monolix_st_frontal_dof2.csv"
DST_ALL      = "spartacus_angles_long.csv"

SUBSET = {                       # the feature-example filter
    "joint": "scapulothoracic",
    "humeral_motion": "frontal plane elevation",
    "degree_of_freedom": "2",
    "unit": "rad",
}


def is_number(s):
    try:
        float(s)
        return True
    except (TypeError, ValueError):
        return False


def main():
    n_all = 0
    n_sub = 0
    ids_all = set()
    ids_sub = set()

    with open(SRC, newline="") as fin, \
         open(DST_ALL, "w", newline="") as fall, \
         open(DST_SUBSET, "w", newline="") as fsub:

        reader = csv.DictReader(fin)
        w_all = csv.writer(fall)
        w_sub = csv.writer(fsub)
        # NB: `study` (= article) is preserved separately from ID so that
        # study-level effects can be fitted/audited (condition is confounded
        # with study — see STATISTICAL_METHODS_REVIEW.md).
        w_all.writerow(["joint", "humeral_motion", "degree_of_freedom",
                        "study", "ID", "TIME", "Y", "in_vivo"])
        w_sub.writerow(["study", "ID", "TIME", "Y", "in_vivo"])

        for row in reader:
            if row.get("unit") != "rad":            # angular data only
                continue
            x, y = row.get("humerothoracic_angle"), row.get("value")
            if not (is_number(x) and is_number(y)):  # need finite TIME/Y
                continue

            study = row["article"]
            uid = f"{study}_{row['shoulder_id']}"
            in_vivo = row.get("in_vivo", "")

            w_all.writerow([row["joint"], row["humeral_motion"],
                            row["degree_of_freedom"], study, uid, x, y, in_vivo])
            n_all += 1
            ids_all.add(uid)

            if all(row.get(k) == v for k, v in SUBSET.items()):
                w_sub.writerow([study, uid, x, y, in_vivo])
                n_sub += 1
                ids_sub.add(uid)

    print(f"Wrote {DST_ALL}")
    print(f"  angular rows: {n_all}   distinct ID: {len(ids_all)}")
    print(f"Wrote {DST_SUBSET}  (feature subset)")
    print(f"  rows: {n_sub}   distinct ID: {len(ids_sub)}")

    if n_all == 0 or n_sub == 0:
        print("WARNING: no rows matched — check column values.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
