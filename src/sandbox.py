#!python

import numpy as np
import numba


@numba.njit
def determine_isotopic_rt_difference(
    mzs,
    rts,
    isotopic_distance,
    ppm,
    best_hit_only=False
):
    mz_order = np.argsort(mzs)
    mzs_in_mz_order = mzs[mz_order]
    rts_in_mz_order = rts[mz_order]
    if isotopic_distance > 0:
        lower_limits = np.searchsorted(
            mzs_in_mz_order,
            (mzs_in_mz_order + isotopic_distance) / (1 + ppm * 10**-6),
            "left"
        )
    else:
        lower_limits = np.arange(len(mzs)) + 1
    upper_limits = np.searchsorted(
        mzs_in_mz_order,
        (mzs_in_mz_order + isotopic_distance) * (1 + ppm * 10**-6),
        "right"
    )
    if not best_hit_only:
        first_rts = np.repeat(rts_in_mz_order, upper_limits - lower_limits)
        second_rts = np.concatenate(
            [
                rts_in_mz_order[l: h] for l, h in zip(
                    lower_limits,
                    upper_limits
                )
            ]
        )
    else:
        isotopes_found = upper_limits > lower_limits
        first_rts = rts_in_mz_order[isotopes_found]
        lower_limits = lower_limits[isotopes_found]
        upper_limits = upper_limits[isotopes_found]
        second_rts = np.array(
            [
                np.min(rts_in_mz_order[l: h]) for l, h in zip(
                    lower_limits,
                    upper_limits
                )
            ]
        )
    isotope_rt_diffs = np.abs(second_rts - first_rts)
    isotope_rt_diffs, isotope_count = np.unique(
        isotope_rt_diffs,
        return_counts=True
    )
    isotope_count = np.cumsum(isotope_count)
    return isotope_rt_diffs, isotope_count
