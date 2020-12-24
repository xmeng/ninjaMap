#!/usr/bin/env python3

import json
import logging
import os
import sys
from collections import defaultdict

import matplotlib
import numpy as np
import numpy.ma as ma
import pandas as pd
from scipy import stats

matplotlib.use("agg")
import matplotlib.pyplot as plt

logging.getLogger(__name__).addHandler(logging.NullHandler())


def combine_profiles(list_of_dfs, genome_size, abs_min_depth):
    """
    Across all Samples
        For each position,
            calculate probability for A, C, G and T
            calculate std error as sqrt[P * (1-P)/N], where P is probability for Nuc and N is the depth at this position.
    Returns: 3 numpy arrays, each array is the size of the genome.
        probs (dimensions = 4L x genome_size) = probability for each base at this position [[pA,pC,pG,pT] ... ]
        err = (dimensions = 4L x genome_size) = standard error at this position [[eA,eC,eG,eT] ... ]
        depth_vector (dimensions = 1L x genome_size) = total depth at this position [d0,d1,d2,d3 ... ]

    """
    num_files = len(list_of_dfs)
    logging.info(f"Found {num_files} file(s)")
    probs = np.zeros((genome_size, 4))
    prob_vector = np.zeros((genome_size, 4))
    depth_vector = np.zeros(genome_size)
    num_items = np.zeros(num_files)
    # per Sample
    for file_num, df_file in enumerate(list_of_dfs):
        # Not reading positions column because these are not reliable when multiple contigs are present.
        data_frame = pd.read_csv(
            df_file,
            header=0,
            usecols=["A", "C", "G", "T", "total_depth"],
            dtype={
                "A": "float",
                "C": "float",
                "G": "float",
                "T": "float",
                "total_depth": "int",
            },
        )
        num_items[file_num], _ = data_frame.shape
        logging.info(
            f"\t[{file_num+1}/{num_files}] Processing {os.path.basename(df_file)} with {num_items[file_num]} rows ..."
        )
        # per Position
        for row in data_frame.itertuples(index=True):
            idx = int(row.Index)
            prob_vector[idx] += np.array([row.A, row.C, row.G, row.T])
            depth_vector[idx] += row.total_depth

    # all arrays should be the same size (genome size). fail if not.
    assert (num_items == num_items[0]).all(), "Genome data are not of the same size."

    # abs_max_depth = 500
    depth_stats = summary_stats(depth_vector)
    min_depth = max(depth_stats["Q10"], int(abs_min_depth / 2))
    max_depth = max(depth_stats["Q90"], abs_min_depth * 2)
    logging.info(
        f"Calculating positional probabilities and std errors for positions between depths {min_depth} and {max_depth}"
    )
    # Aggregate across all samples
    for idx, _ in enumerate(prob_vector):
        if (
            np.isnan(depth_vector[idx])
            or (depth_vector[idx] < min_depth)
            or (depth_vector[idx] > max_depth)
        ):
            probs[idx] = np.nan
        else:
            probs[idx] = prob_vector[idx] / depth_vector[idx]

    return (probs, depth_vector, depth_stats)


def summary_stats(d):
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


def read_list_file(filename):
    file_paths_list = list()
    logging.info(f"Reading {filename} ...")
    with open(filename, "r") as file:
        file_paths_list = [line.strip() for line in file]

    return file_paths_list


def plot_diff(ax, diff_prob_vector, pos_depth, neg_depth, title):
    logging.info(f"Generating figure of differences ...")

    # Set min/max plotting range
    ax.set_ylim([-1.15, 1.15])
    ax.set_yticks(np.arange(-1, 1.15, 0.25))
    ax.set_title(title)

    # Add depth tracks
    ax.plot(pos_depth / np.nanmax(pos_depth), linewidth=1, color="b", alpha=0.75)
    ax.plot(-(neg_depth / np.nanmax(neg_depth)), linewidth=1, color="g", alpha=0.75)

    # Add prob tracks
    ax.plot(diff_prob_vector, "o", markersize=1, color="r", alpha=0.25)
    ax.set_rasterized(True)


def plot_diff_hist(ax, diff_prob_vector, title):
    logging.info(f"Generating histogram of differences ...")
    ax.set_yscale("log")
    ax.set_title(title)
    ax.hist(diff_prob_vector.flatten(), bins=np.arange(-1, 1.01, 0.01))
    ax.set_xlim([-1.15, 1.15])
    ax.set_rasterized(True)


def vector_difference(array1, array2, nan_value=np.nan):
    """
    subtract two vectors, remove nan values, 
    return vector of difference
    return vector may not be of the same length as input arrays
    """
    # print(array1.shape)
    logging.info(f"Comparing vectors ...")
    if array1.shape != array2.shape:
        return

    # genome_size = len(array1)
    diff = array1 - array2

    # If both arrays had the value 0 at the same position for the same nucl
    # make the nucl at that position a nan, so as not to inflate the number of 0s.
    # 0 would mean that both the arrays had equal contribution for that nucl at that
    # position and hence canceled each other out.
    diff[(array1 == 0) & (array2 == 0)] = np.nan

    ## if any nucl(ATGC) at a position is NA remove that position and
    ## calculate the length of the vector
    ## Replaced: Too aggressive, removes positions where other nucleotides
    ## may have legit values.
    # diff_no_na = diff[~np.isnan(diff).any(axis=1)]

    ## if all nucl(ATGC) at a position are NAs remove that position and
    ## calculate the length of the vector
    diff_no_na = diff[~np.isnan(diff).all(axis=1)]

    len_diff_no_na = len(diff_no_na)
    logging.info(f"\tShapes with NA: {diff.shape}, without NA: {diff_no_na.shape}")
    logging.info(f"\tLengths with NA: {len(diff)}, without NA: {len_diff_no_na}")

    if np.isnan(nan_value):
        return diff
    # elif len_diff_no_na == 0:
    #     return np.full_like(diff, nan_value)
    else:
        return np.nan_to_num(diff, nan=nan_value)


def remove_nan(prob_array):
    return prob_array[~np.isnan(prob_array)]


def array_properties(diff):
    array_properties_dict = dict()
    flat_diff = remove_nan(diff.flatten())
    array_properties_dict.update(summary_stats(diff))

    # Query
    q_top = len(flat_diff[(-1 <= flat_diff) & (flat_diff <= -0.75)])
    q_middle = len(flat_diff[(-0.25 <= flat_diff) & (flat_diff <= 0)])
    q_zero = len(flat_diff[(flat_diff == 0)])

    # Absolute
    abs_diff = abs(flat_diff)
    a_top = len(abs_diff[abs_diff >= 0.75])
    a_middle = len(abs_diff[(0 <= abs_diff) & (abs_diff <= 0.25)])
    a_zero = len(abs_diff[(abs_diff == 0)])

    array_properties_dict.update(
        {
            "query_side_top": q_top,
            "query_side_middle": q_middle,
            "absolute_top": a_top,
            "absolute_middle": a_middle,
            "q_zero": q_zero,
            "a_zero": a_zero,
        }
    )
    return array_properties_dict


def visualizing_differences(base_dict, query_dict, genome_size, diff_ax, hist_ax):
    # base_dict = profile_dict[base]
    # query_dict = profile_dict[query]
    plot_title = (
        f"{base_dict['sample_name']} "
        f"(Depth q10:{base_dict['depth_stats']['Q10']}, median:{base_dict['depth_stats']['Median']}, q90:{base_dict['depth_stats']['Q90']}) "
        f"vs "
        f"{query_dict['sample_name']} "
        f"(Depth q10:{query_dict['depth_stats']['Q10']}, median:{query_dict['depth_stats']['Median']}, q90:{query_dict['depth_stats']['Q90']})"
    )

    logging.info(f"Plotting {plot_title}")

    ## Isolate differences between two pruned profiles
    prob_diff = vector_difference(base_dict["prob"], query_dict["prob"])

    diff_name = f"{prefix}.{base_dict['sample_name']}_vs_{query_dict['sample_name']}.prob_diff.npy"
    np.save(diff_name, prob_diff)

    # Plot difference clouds
    try:
        plot_diff(
            diff_ax, prob_diff, base_dict["depth"], query_dict["depth"], plot_title,
        )
    except:
        logging.info(f"Not enough data to plot a diff cloud for {plot_title}.")

    # Plot difference histograms
    try:
        plot_diff_hist(hist_ax, prob_diff, plot_title)
    except:
        logging.info(f"Not enough data to plot a histogram for {plot_title}.")

    return array_properties(prob_diff)


def generate_profile(query, min_median_depth, genome_size, profiles):
    if query in profiles:
        logging.info(f"Profile already exists for {query} ...")
        return profiles[query]
    else:
        logging.info(f"Calculating profile for {query} ...")
        name = os.path.basename(query).split(".")[0].split("_vs_")[0]
        (probs, depth, depth_stats,) = combine_profiles(
            [query], genome_size=genome_size, abs_min_depth=min_median_depth
        )
    return {
        "prob": probs,
        "depth": depth,
        "sample_name": name,
        "depth_stats": depth_stats,
    }


def save_stats(profiles, genome_size, prefix):
    stats_df_list = list()
    for sample_path in profiles.keys():
        stats_dict = profiles[sample_path]["depth_stats"]
        stats_dict.update(
            {"name": profiles[sample_path]["sample_name"], "genome_size": genome_size}
        )
        stats_df_list.append(stats_dict)
    df = pd.DataFrame(stats_df_list)
    df.to_csv(f"{prefix}.02b_dist_stats.csv", index=False)


def has_sufficient_depth(median_depth, min_median_depth):
    return median_depth <= min_median_depth


if __name__ == "__main__":
    # logging.getLogger(__name__)
    logging.basicConfig(
        # filename=logfile,
        # filemode='w+',
        level=logging.INFO,
        format="%(asctime)s\t[%(levelname)s]:\t%(message)s",
    )

    prefix = sys.argv[1]
    comparison_file = sys.argv[2]
    min_median_depth = 8
    compare_df = pd.read_table(
        comparison_file, header=None, index_col=None, names=["base", "query"]
    )

    # Setup matplotlib figures
    num_comparisons, num_cols = compare_df.shape
    logging.info(f"Making {num_comparisons} comparisons ...")
    plot_height = int(num_comparisons * 10)
    plt.rcParams["figure.figsize"] = (45, plot_height)

    figure1, diff_axs = plt.subplots(num_comparisons, 1, sharex=True)
    figure2, hist_axs = plt.subplots(num_comparisons, 1, sharex=True, sharey=True)

    # Calculate the profiles for each path and store as dict values with paths as keys.
    profiles = dict()
    diff_df_list = list()
    diff_probs_list = list()
    genome_size = 0
    for row in compare_df.itertuples():
        idx = int(row.Index)
        logging.info(f"Comparison [{idx+1}/{num_comparisons}]")
        if genome_size == 0:
            logging.info("Estimating genome size ...")
            genome_size, _ = pd.read_csv(row.base, header=0, usecols=["pos"]).shape
            logging.info(f"\t{genome_size} bases.")
        # Check profile dictionary to not calclulate the profile if already calculated.
        # Base
        profiles[row.base] = generate_profile(
            row.base, min_median_depth, genome_size, profiles
        )

        # Query
        profiles[row.query] = generate_profile(
            row.query, min_median_depth, genome_size, profiles
        )

        # Plot by pairs
        pair_diff_data = visualizing_differences(
            profiles[row.base],
            profiles[row.query],
            genome_size,
            diff_axs[idx],
            hist_axs[idx],
        )

        pair_diff_data.update(
            {
                "genome_size": genome_size,
                "base": profiles[row.base]["sample_name"],
                "query": profiles[row.query]["sample_name"],
            }
        )
        diff_df_list.append(pair_diff_data)

    ## Save profiles to disk
    output_df = pd.DataFrame(diff_df_list)
    output_df.to_csv(f"{prefix}.02b_strain_profiles.csv", index=False)

    ## Save stats
    save_stats(profiles, genome_size, prefix)

    # Save Diff Plot
    # diff_plot = f"{prefix}.02b_diff_plot.pdf"
    diff_plot = f"{prefix}.02b_diff_plot.png"
    logging.info(f"Saving difference plot to disk as: {diff_plot}")
    figure1.savefig(diff_plot, dpi=80, bbox_inches="tight")

    # Save Historgams Plot
    # hist_plot = f"{prefix}.02b_hist_plot.pdf"
    hist_plot = f"{prefix}.02b_hist_plot.png"
    logging.info(f"Saving histogram figure to disk as: {hist_plot}")
    figure2.savefig(hist_plot, dpi=80, bbox_inches="tight")

    logging.info("Huzzah!")
