#!/usr/bin/env python3

import os, sys
import glob
import numpy as np
import pandas as pd
import logging
from scipy import stats
import concurrent.futures


def ks_test(a, b):
    """
    Accepts two distributions a and b
    Returns: tuple (KS Distance and p-value)
    # Significant p-value means a and b belong to different distributions
    """
    a_prime = remove_nan(a.flatten())
    b_prime = remove_nan(b.flatten())
    if len(a_prime) == 0 or len(b_prime) == 0:
        return [np.nan] * 2
    return stats.ks_2samp(a_prime, b_prime)


def ad_test(a, b):
    """
    statistic (float)
        Normalized k-sample Anderson-Darling test statistic.

    critical_values (array)
        The critical values for significance levels 25%, 10%, 5%, 2.5%, 1%, 0.5%, 0.1%.

    significance_level (float)
        An approximate significance level at which the null hypothesis for the provided samples can be rejected. The value is floored / capped at 0.1% / 25%.
    """
    a_prime = remove_nan(a.flatten())
    b_prime = remove_nan(b.flatten())
    if len(a_prime) == 0 or len(b_prime) == 0:
        return [np.nan] * 3
    return stats.anderson_ksamp([a_prime, b_prime], midrank=True)


def positional_drift(base_array, query_array, cutoff=0.002):
    """
    Description:
        Identify positions that have a low deviation from the center in the base and
        compare those positions in the query. Return the number of positions that drifted 
        in the query
    Accept: (numpy array1, numpy array2, float)
        base    diff array for the base
        query   diff array for the query
        cutoff  deviation cutoff from center
    Returns: (int1, int2, int3) 
        int1: number of positions within cutoff in base
        int2: number of positions that can be compared between base and query
        int3: number of positions (out of int1) outside cutoff in query

    Derivation with KC:
    mouse1 = np.load("Bacteroides-uniformis-ATCC-8492.W4M1_vs_W8M1.prob_diff.npy")
    nan_pos_mouse1 = ~np.isnan(mouse1)
    mouse1_no_nans = mouse1[nan_pos_mouse1]
    mouse1_no_nans[abs(mouse1_no_nans) > 0.002].shape

    mouse7 = np.load("Bacteroides-uniformis-ATCC-8492.W4M7_vs_W8M7.prob_diff.npy")
    mouse7_no_nans = mouse7[nan_pos_mouse1]
    filter_pos_mouse7 = np.array(~np.isnan(mouse7_no_nans) & (abs(mouse1_no_nans) < 0.002))
    mouse7_filter = mouse7_no_nans[filter_pos_mouse7]
    mouse7_filter[abs(mouse7_filter) > 0.002].shape
    """

    assert base_array.shape == query_array.shape, "array size mismatch!"

    base = base_array.flatten()
    query = query_array.flatten()
    # Get index of positions that have non-NaN values in Base
    base_non_nan_pos_idx = ~np.isnan(base)

    # Filter out NaN positions
    base_no_nans = base[base_non_nan_pos_idx]
    # base_len = len(base_no_nans[abs(base_no_nans) > cutoff])
    base_conserved_len = len(base_no_nans[abs(base_no_nans) < cutoff])

    # if no positions are found inside the cutoff, there is nothing to compare
    if base_conserved_len == 0:
        query_comparable_pos = np.nan
        query_drift = np.nan
    else:
        # Select positions in Query where Base has non-NAN values
        query_no_nans = query[base_non_nan_pos_idx]

        # Number of positions that can be compared (w/o NaNs)
        # query_comparable_pos = len(query_no_nans[~np.isnan(query_no_nans)])

        # Select positions that are not NaNs in Query
        # AND
        # absolute value of Base is less than cutoff
        query_filtered = query_no_nans[
            (~np.isnan(query_no_nans)) & (abs(base_no_nans) < cutoff)
        ]

        # Calc the number of common positions that have changed significantly in Query
        query_drift = len(query_filtered[abs(query_filtered) > cutoff])

    return (base_conserved_len, query_drift)


def remove_nan(prob_array):
    return prob_array[~np.isnan(prob_array)]


def parse_critical_values(crit_array):
    if np.isnan(crit_array).all():
        crit_array = [np.nan] * 7

    critical_dict = {}
    keys = ["p25", "p10", "p5", "p2_5", "p1", "p0_5", "p0_1"]
    for idx, value in enumerate(crit_array):
        critical_dict.update({keys[idx]: value})

    return critical_dict


def compare_arrays_to_base(np_base, query_array_path, cutoff=None):
    query_dict = {}
    name_split = os.path.basename(query_array_path).split(".")
    query_name = name_split[1]
    np_query = np.load(query_array_path)
    logging.info(f" ... {query_name}")

    ks_dist, ks_p_val = ks_test(np_base, np_query)
    # ad_stat, ad_critical_val, p_val_level = ad_test(np_base, np_query)
    # critical_dict = parse_critical_values(ad_critical_val)
    base_conserved_len, query_drift = positional_drift(np_base, np_query, cutoff=cutoff)
    query_dict = {
        "query_name": query_name,
        "ks_dist": ks_dist,
        "ks_p_val": ks_p_val,
        # "ad_stat": ad_stat,
        # "p_val_level": p_val_level,
        "base_conserved_pos": base_conserved_len,
        # "query_comparable_pos": query_comparable_pos,
        "query_conserved_pos_drift": query_drift,
    }

    # query_dict.update(critical_dict)
    return query_dict


def compare_all(query_array_paths, base_array_paths, cores=None, cutoff=None):
    df_dict_array = list()
    num_profiles = len(query_array_paths)
    logging.info(f"{num_profiles} profiles found.")
    processed = dict()
    for outer_idx, base_array in enumerate(base_array_paths):
        # logging.info(f"[{outer_idx+1}/{num_profiles}] Starting with {base_array}")
        name_split = os.path.basename(base_array).split(".")
        org_name = name_split[0]
        base_name = name_split[1]

        np_base = np.load(base_array)
        query_dict = {}
        logging.info(f"[{outer_idx+1}/{num_profiles}] Comparing {base_name} with ...")

        # Parallelize comparisons to base
        with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor:
            future = [
                executor.submit(compare_arrays_to_base, np_base, query_path, cutoff)
                for query_path in query_array_paths
            ]
            for f in concurrent.futures.as_completed(future):
                result = f.result()
                # logging.info(
                #     f"Saving results for {base_name} vs {result['query_name']}"
                # )
                result.update({"base_name": base_name, "org_name": org_name})
                df_dict_array.append(result)

    return df_dict_array


if __name__ == "__main__":
    logging.basicConfig(
        # filename=logfile,
        # filemode='w+',
        level=logging.INFO,
        format="%(asctime)s\t[%(levelname)s]:\t%(message)s",
    )

    prefix = sys.argv[1]
    control_profiles = sys.argv[2]
    np_folder_path = sys.argv[3].rstrip("/")
    cores = None
    cutoff = 0.002

    with open(control_profiles) as control_profile_files_list:
        control_profile_paths = [
            filepath.rstrip("\n") for filepath in control_profile_files_list
        ]

    all_array_paths = glob.glob(f"{np_folder_path}/*.npy")
    query_array_paths = set(all_array_paths) - set(control_profile_paths)
    # [print(query_path) for query_path in query_array_paths]

    df = pd.DataFrame(
        compare_all(
            query_array_paths, control_profile_paths, cores=cores, cutoff=cutoff
        )
    )
    df.to_csv(f"{prefix}.02c_compare.csv", index=False)

    logging.info("All Done. Huzzah!")
