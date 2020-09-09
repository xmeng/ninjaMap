#!/usr/bin/env python3

import concurrent.futures

# import glob
import itertools
import logging
import os
import re
import sys

import matplotlib
import numpy as np

# import numpy.ma as ma
import pandas as pd
from scipy import stats

matplotlib.use("agg")
import matplotlib.pyplot as plt


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


def cos_sim(a, b):
    """calculate cosine distance between reference and given."""
    a_prime = a.flatten()
    b_prime = b.flatten()
    if len(a_prime) == len(b_prime):
        return np.dot(a_prime, b_prime) / (
            np.linalg.norm(a_prime) * np.linalg.norm(b_prime)
        )
    else:
        raise Exception("LengthError: Cos_Sim requires arrays to be of the same size.")


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
        query_comparable_pos = len(query_no_nans[~np.isnan(query_no_nans)])

        # Select positions that are not NaNs in Query
        # AND
        # absolute value of Base is less than cutoff
        query_filtered = query_no_nans[
            (~np.isnan(query_no_nans)) & (abs(base_no_nans) < cutoff)
        ]

        # Calc the number of common positions that have changed significantly in Query
        query_drift = len(query_filtered[abs(query_filtered) > cutoff])

    return (base_conserved_len, query_comparable_pos, query_drift)


def remove_nan(prob_array):
    return prob_array[~np.isnan(prob_array)]


def summary_stats(d, cutoff=0.002):
    """
    Accept: an n-dim numpy array
    Do: flatten array and remove np.nan values
    Return: dict of summary stats of array with keys:
        min,max,mean,median,q10,q90,mad,iqr
    """
    d_prime = remove_nan(d.flatten())
    d_prime_len = len(d_prime)
    if d_prime_len > 0:
        d_stats = {
            "Total_Length": len(d.flatten()),
            "Non_NA_Length": d_prime_len,
            "Non_Zero_Length": len(d_prime[d_prime != 0]),
            "Extreme_Positions": len(d_prime[abs(d_prime) == 1]),
            "Within_Threshold": len(d_prime[abs(d_prime) < cutoff]),
            "Min": np.nanmin(d_prime),
            "Max": np.nanmax(d_prime),
            "Mean": np.nanmean(d_prime),
            "Std_Dev": np.nanstd(d_prime),
            "Variance": np.nanvar(d_prime),
            "Q10": np.nanquantile(d_prime, 0.1),
            "Median": np.nanmedian(d_prime),
            "Q90": np.nanquantile(d_prime, 0.9),
            "MAD": stats.median_absolute_deviation(d_prime, nan_policy="omit"),
            "IQR": stats.iqr(d_prime, nan_policy="omit"),
            "Skew": stats.skew(d_prime, nan_policy="omit"),
            "Kurtosis": stats.kurtosis(d_prime, nan_policy="omit"),
        }
    else:
        logging.info(
            "Not enough data points for reliable summary statistics. Adding placeholders ..."
        )
        nan_value = np.nan
        d_stats = {
            "Total_Length": len(d.flatten()),
            "Non_NA_Length": d_prime_len,
            "Non_Zero_Length": nan_value,
            "Extreme_Positions": nan_value,
            "Within_Threshold": nan_value,
            "Min": nan_value,
            "Max": nan_value,
            "Mean": nan_value,
            "Std_Dev": nan_value,
            "Variance": nan_value,
            "Q10": nan_value,
            "Median": nan_value,
            "Q90": nan_value,
            "MAD": nan_value,
            "IQR": nan_value,
            "Skew": nan_value,
            "Kurtosis": nan_value,
        }
    summ = "; ".join([f"{k}:{d_stats[k]}" for k in sorted(d_stats.keys())])
    logging.info(f"Depth summary: {summ}")
    return d_stats


def summary_stats_select(d, cutoff=0.002):
    """
    Accept: an n-dim numpy array
    Do: flatten array and remove np.nan values
    Return: dict of summary stats of array with keys:
        min,max,mean,median,q10,q90,mad,iqr
    """
    d_prime = remove_nan(d.flatten())
    d_prime_len = len(d_prime)
    if d_prime_len > 0:
        d_stats = {
            "Total_Length": len(d.flatten()),
            "Non_NA_Length": d_prime_len,
            "Non_Zero_Length": len(d_prime[d_prime != 0]),
            "Extreme_Positions": len(d_prime[abs(d_prime) == 1]),
            "Within_Threshold": len(d_prime[abs(d_prime) < cutoff]),
            "Skew": stats.skew(d_prime, nan_policy="omit"),
            "Kurtosis": stats.kurtosis(d_prime, nan_policy="omit"),
        }
    else:
        logging.info(
            "Not enough data points for reliable summary statistics. Adding placeholders ..."
        )
        nan_value = np.nan
        d_stats = {
            "Total_Length": len(d.flatten()),
            "Non_NA_Length": d_prime_len,
            "Non_Zero_Length": nan_value,
            "Extreme_Positions": nan_value,
            "Within_Threshold": nan_value,
            "Skew": nan_value,
            "Kurtosis": nan_value,
        }
    summ = "; ".join([f"{k}:{d_stats[k]}" for k in sorted(d_stats.keys())])
    logging.info(f"Depth summary: {summ}")
    return d_stats


def compare_arrays_to_base(
    np_base_complete, query_array_path, selected_pos=None, cutoff=None
):
    query_dict = {}
    name_split = os.path.basename(query_array_path).split(".")
    query_name = name_split[1]
    np_query_complete = np.load(query_array_path)
    logging.info(f" ... {query_name}")

    if selected_pos is None:
        np_base = np_base_complete
        np_query = np_query_complete
    else:
        np_base = np_base_complete[selected_pos]
        np_query = np_query_complete[selected_pos]

    ks_dist, ks_p_val = ks_test(np_base, np_query)
    cosine_sim = cos_sim(np_base, np_query)
    base_conserved_len, query_comparable_pos, query_drift = positional_drift(
        np_base, np_query, cutoff=cutoff
    )
    stats = summary_stats(np_query, cutoff=cutoff)
    query_dict = {
        "query_name": query_name,
        "ks_dist": ks_dist,
        "ks_p_val": ks_p_val,
        "base_conserved_pos": base_conserved_len,
        "query_comparable_pos": query_comparable_pos,
        "query_conserved_pos_drift": query_drift,
        "cosine_similarity": cosine_sim,
    }

    query_dict.update(stats)
    return query_dict


def compare_all(
    query_array_paths, base_array_paths, selected_pos=None, cores=None, cutoff=None
):
    df_dict_array = list()
    num_profiles = len(query_array_paths)
    logging.info(f"{num_profiles} profiles found.")
    # processed = dict()
    for outer_idx, base_array in enumerate(base_array_paths):
        # logging.info(f"[{outer_idx+1}/{num_profiles}] Starting with {base_array}")
        name_split = os.path.basename(base_array).split(".")
        org_name = name_split[0]
        base_name = name_split[1]

        np_base = np.load(base_array)
        # query_dict = {}
        logging.info(f"[{outer_idx+1}/{num_profiles}] Comparing {base_name} with ...")

        # Parallelize comparisons to base
        with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor:
            future = [
                executor.submit(
                    compare_arrays_to_base, np_base, query_path, selected_pos, cutoff
                )
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


def plot_diff_hist(ax, diff_vector, title):
    logging.info(f"Generating histogram for {title} ...")
    ax.set_yscale("log")
    ax.set_title(title)
    ax.hist(diff_vector.flatten(), bins=np.arange(-1, 1.01, 0.01))
    ax.set_xlim([-1.15, 1.15])
    ax.set_rasterized(True)


def visualizing_differences(diff_vector, selected_pos, axs):
    title = os.path.basename(diff_vector).split(".")[1]
    loaded_diff = np.load(diff_vector)
    base_name = title.split("_vs_")[0]
    query_name = title.split("_vs_")[1]
    # print(selected_pos.shape)
    # masking doesn't seem to work with matplotlib
    # masked_diff = ma.array(
    #     loaded_diff, mask=np.logical_not(selected_pos), fill_value=np.nan
    # )
    # masked_diff = np.ma.masked_where(np.logical_not(selected_pos), loaded_diff)
    # loaded_diff[selected_pos] = np.nan
    # print(loaded_diff.shape)
    # plot_diff_hist(axs, loaded_diff, title)

    masked_diff = loaded_diff[selected_pos]
    # print(masked_diff.shape)
    plot_diff_hist(axs, masked_diff, title)
    # return summary_stats(masked_diff).update({"Base": base_name, "Query": query_name})
    stat_dict = summary_stats_select(masked_diff)
    stat_dict.update({"Base": base_name, "Query": query_name, "Sample": title})
    return stat_dict


def within_limits(array, cutoff):
    """
    returns a boolean array as int(0 or 1)
    """
    logging.info(f"Loading {array} at {cutoff}")
    loaded = np.load(array)
    with np.errstate(invalid="ignore"):
        # When comparing, nan always returns False.
        # return values as 0/1
        return (abs(loaded) < cutoff).astype(int)


def select_conserved_positions(arrays, method="majority", cutoff=0.002):
    """
    Description:
        Given a list of diff prob arrays and a absolute cutoff value, 
        return a vector based on the method containing the:
        - intersection of all positions within the cutoff; or
        - union of all positions within the cutoff; or
        - most prevelant (3/5) all positions within the cutoff
    
    Return: np.ndarray (boolean)
    """
    num_arrays = len(arrays)
    logging.info(
        f"Received {num_arrays} arrays, with method {method} and cutoff {cutoff}"
    )
    ## for each array:
    ##   get an int boolean (0/1) of indices within abs(cutoff)
    ##   sum each base value at each position across all arrays
    agg_array = np.sum([within_limits(array, cutoff) for array in arrays], axis=0)

    # Make sure I have the axis right.
    # If any position's total is greater than the number of arrays,
    # Means I'm adding along the wrong axis
    assert np.sum((agg_array > num_arrays).astype(int)) == 0, "Wrong axis!"
    selected_idx = []
    if method == "majority":
        if num_arrays % 2 == 1:
            majority_at = (num_arrays - 1) / 2
        else:
            majority_at = (num_arrays) / 2
        # majority rule
        ##   if value > majority --> majority of arrays have this base at this position within the cutoff.
        selected_idx = (agg_array > majority_at).astype(bool)
    elif method == "intersection":
        # intersection index
        ##   if value == num_arrays --> all arrays have this base at this position within the cutoff.
        selected_idx = (agg_array == num_arrays).astype(bool)
    elif method == "union":
        # union index
        ##   if value > 0 --> at least one array has this base at this position within the cutoff.
        selected_idx = (agg_array > 0).astype(bool)

    selected_idx_num = np.sum(selected_idx)
    logging.info(
        f"Method: {method}, Shape:{selected_idx.shape}, Num_Positions:{selected_idx_num}"
    )

    return selected_idx


def get_conserved_positions(arr_paths, method, cutoff=0.002, cores=None):
    """take a list of N arrays
    find conserved positions based on all combinations of N-1 arrays
    assign scores to each conserved position based on it's occurance in each combination

    Args:
        arr ([type]): [description]
        method ([type]): [description]
    """
    num_arr = len(arr_paths)
    arr_combinations = list()
    logging.info(f"Found {num_arr} arrays")
    # based on all combinations of N-1 arrays
    index_combinations = list(itertools.combinations(range(num_arr), num_arr - 1))
    for idx_combo in index_combinations:
        arr_combinations.append([arr_paths[i] for i in idx_combo])
    logging.info(f"Created {len(arr_combinations)} combinations")

    # find conserved positions across N-1 arrays
    all_selected_indices = list()
    with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor:
        future = [
            executor.submit(select_conserved_positions, arr_sub, method, cutoff)
            for arr_sub in arr_combinations
        ]
        for f in concurrent.futures.as_completed(future):
            result = f.result()
            all_selected_indices.append(result.astype(int))

    selected_idx, pos_stats = pick_longest(all_selected_indices)
    summary = ";\t".join([f"{k}:{pos_stats[k]}" for k in sorted(pos_stats.keys())])
    logging.info(f"[Leave One Out]: {summary}")

    selected_idx_num = np.sum(selected_idx)
    logging.info(
        f"[Leave One Out] Method: {method}, Shape:{selected_idx.shape}, Num_Positions:{selected_idx_num}"
    )

    return (selected_idx.astype("bool"), pos_stats)


def pick_longest(arr):
    """[summary]

    Args:
        arr ([type]): [description]

    Returns:
        [type]: [description]
    """
    all_lengths = [np.sum(a) for a in arr]
    pos_stats = {
        "Min": np.nanmin(all_lengths),
        "Max": np.nanmax(all_lengths),
        "Mean": np.nanmean(all_lengths),
        "Std_Dev": np.nanstd(all_lengths),
        "Variance": np.nanvar(all_lengths),
    }
    return (arr[np.argmax(all_lengths)], pos_stats)


def read_list_file(list_file):
    with open(list_file) as files_list:
        paths = [filepath.rstrip("\n") for filepath in files_list]
    return paths


def infer_week_and_mouse(sample_name):
    valid_format = r"W\d+M\d+"
    p = re.compile(valid_format)
    if p.match(sample_name):
        week_fmt = re.compile(r"^W\d+")
        week = week_fmt.search(sample_name).group().replace("W", "")

        mouse_fmt = re.compile(r"M\d+")
        mouse = mouse_fmt.search(sample_name).group().replace("M", "")
        return week, mouse
    else:
        return np.nan, np.nan


def get_col_mean_median(df, column):
    # col_arr = df[[column]].to_numpy()
    # col_mean = np.nanmean(col_arr)
    # col_median = np.nanmedian(col_arr)

    col_mean = df.loc[:, column].mean(skipna=True)
    col_median = df.loc[:, column].median(skipna=True)

    return col_mean, col_median


def selected_pos_depths(csv_file, informative_pos):
    # Find all positions where at least one nucleotide has a value
    selected_pos_idx = np.sum(informative_pos, axis=1) > 0
    selected_pos_len = len(selected_pos_idx)

    data_frame = pd.read_csv(
        csv_file,
        header=0,
        usecols=["ref_depth", "total_depth"],
        dtype={"ref_depth": "float", "total_depth": "float",},
    )

    og_ref_mean, og_ref_median = get_col_mean_median(data_frame, "ref_depth")
    og_total_mean, og_total_median = get_col_mean_median(data_frame, "total_depth")

    genome_len, _ = data_frame.shape
    assert (
        selected_pos_len == genome_len
    ), f"Arrays do not have the same shape ({selected_pos_len} vs {genome_len}). Can't compare."

    selected_df = data_frame[selected_pos_idx]
    sel_ref_mean, sel_ref_median = get_col_mean_median(selected_df, "ref_depth")
    sel_total_mean, sel_total_median = get_col_mean_median(selected_df, "total_depth")

    logging.info(
        f"{og_ref_mean}, {og_ref_median},"
        f"{og_total_mean}, {og_total_median},"
        f"{sel_ref_mean}, {sel_ref_median},"
        f"{sel_total_mean}, {sel_total_median}"
    )
    return (
        sel_ref_mean,
        sel_ref_median,
        sel_total_mean,
        sel_total_median,
        og_ref_mean,
        og_ref_median,
        og_total_mean,
        og_total_median,
    )


def get_selected_pos_raw_data(init_profiles_file, informative_pos):
    df = pd.read_table(
        init_profiles_file, header=None, names=["Query", "S3Path"], index_col=["Query"]
    )
    (
        df["Selected_Pos_Mean_Ref_Depth"],
        df["Selected_Pos_Median_Ref_Depth"],
        df["Selected_Pos_Mean_Total_Depth"],
        df["Selected_Pos_Median_Total_Depth"],
        df["Whole_Genome_Mean_Ref_Depth"],
        df["Whole_Genome_Median_Ref_Depth"],
        df["Whole_Genome_Mean_Total_Depth"],
        df["Whole_Genome_Median_Total_Depth"],
    ) = zip(*df["S3Path"].apply(lambda x: selected_pos_depths(x, informative_pos)))

    logging.info(f"{df.head()}")

    return df


def save_predictions(list_of_df, output_file, selected_pos_depth):
    stats_df = pd.DataFrame(list_of_stats)
    # .to_csv(f"{prefix}.02d_stats.csv", index=False)
    output_cols = [
        "Organism",
        "Sample",
        "Mouse",
        "Week",
        "Extreme_Positions",
        "Selected_Pos_Mean_Ref_Depth",
        "Selected_Pos_Median_Ref_Depth",
        "Selected_Pos_Mean_Total_Depth",
        "Selected_Pos_Median_Total_Depth",
        "Whole_Genome_Mean_Ref_Depth",
        "Whole_Genome_Median_Ref_Depth",
        "Whole_Genome_Mean_Total_Depth",
        "Whole_Genome_Median_Total_Depth",
        "Prediction",
    ]
    stats_df["Organism"] = prefix
    stats_df["Week"], stats_df["Mouse"] = zip(
        *stats_df["Query"].apply(infer_week_and_mouse)
    )
    stats_df["Prediction"] = stats_df["Extreme_Positions"].apply(
        lambda x: "Invader" if x > 0 else "Input"
    )

    stats_df.join(selected_pos_depth, on="Query", how="left", rsuffix="_other").to_csv(
        output_file, index=False, columns=output_cols
    )


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s\t[%(levelname)s]:\t%(message)s",
    )

    prefix = sys.argv[1]
    method = sys.argv[2]
    control_profiles = sys.argv[3]
    sample_profiles = sys.argv[4]
    init_profiles = sys.argv[5]

    cutoff = 0.002
    cores = 10

    if method not in ["intersection", "union", "majority"]:
        logging.error(f"Method '{method}' not recognized. Exiting")
        sys.exit(1)

    control_array_paths = read_list_file(control_profiles)
    informative_pos, position_stats = get_conserved_positions(
        control_array_paths, method=method, cutoff=cutoff, cores=cores
    )
    # np.save(f"{prefix}.{method}.informative_pos.npy", informative_pos)
    pos_depth_df = get_selected_pos_raw_data(init_profiles, informative_pos)

    with open(f"{prefix}.02d_loo_pos_stats.tsv", "w") as filehandle:
        header = "\t".join([f"{k}" for k in sorted(position_stats.keys())])
        filehandle.write(f"{header}\n")

        summary = "\t".join(
            [f"{position_stats[k]}" for k in sorted(position_stats.keys())]
        )
        filehandle.write(f"{summary}\n")

    if np.sum(informative_pos) == 0:
        logging.info(
            "Cannot find any informative positions from the control samples. Exiting."
        )
        sys.exit(0)

    query_array_paths = read_list_file(sample_profiles)
    logging.info(
        f"Found {len(control_array_paths)} control profiles and {len(query_array_paths)} query profiles"
    )

    # Setup matplotlib figures
    all_array_paths = control_array_paths + query_array_paths
    num_comparisons = len(all_array_paths)
    logging.info(f"Making {num_comparisons} comparisons ...")
    plot_height = int(num_comparisons * 10)
    plt.rcParams["figure.figsize"] = (45, plot_height)

    figure, hist_axs = plt.subplots(num_comparisons, 1, sharex=True, sharey=True)
    list_of_stats = [
        visualizing_differences(diff_vector, informative_pos, hist_axs[idx])
        for idx, diff_vector in enumerate(all_array_paths)
    ]
    save_predictions(list_of_stats, f"{prefix}.02d_stats.csv", pos_depth_df)

    hist_plot = f"{prefix}.02d_hist_plot.png"
    logging.info(f"Saving histogram figure to disk as: {hist_plot}")
    figure.savefig(hist_plot, dpi=80, bbox_inches="tight")

    # sys.exit(0)

    # df = pd.DataFrame(
    #     compare_all(
    #         query_array_paths,
    #         control_array_paths,
    #         selected_pos=informative_pos,
    #         cores=cores,
    #         cutoff=cutoff,
    #     )
    # )
    # df.to_csv(f"{prefix}.02d_compare.csv", index=False)

    logging.info("All Done. Huzzah!")
