#!/usr/bin/env python3

import os, sys
import numpy as np
import json
import logging
import pandas as pd
from collections import defaultdict

logging.getLogger(__name__).addHandler(logging.NullHandler())


def combine_profiles(
    list_of_dfs, genome_size,
):
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
    err = np.zeros((genome_size, 4))
    prob_vector = np.zeros((genome_size, 4))
    depth_vector = np.zeros(genome_size)
    num_items = np.zeros(num_files)
    # per Sample
    for file_num, df_file in enumerate(list_of_dfs):
        # Not reading positions column because these are not reliable when multiple contigs are present.
        data_frame = pd.read_csv(
            df_file, header=0, usecols=["A", "C", "G", "T", "total_depth"], dtype=int,
        )
        num_items[file_num], ncols = data_frame.shape
        logging.info(
            f"\t[{file_num+1}/{num_files}] Processing {os.path.basename(df_file)} with {num_items[file_num]} rows ..."
        )
        # per Position
        for row in data_frame.itertuples(index=True):
            idx = int(row.Index)
            prob_vector[idx] += np.array([row.A, row.C, row.G, row.T])
            depth_vector[idx] += row.total_depth

    # all arrays should be the same size (genome size). fail if not.
    try:
        assert (num_items == num_items[0]).all()
    except AssertionError as err:
        logging.exception("Genome data are not of the same size.")
        raise err

    logging.info(f"Calculating positional probabilities and std errors.")
    # Aggregate across all samples
    for idx, arr in enumerate(prob_vector):
        if np.isnan(depth_vector[idx]).any() or (depth_vector[idx] < 1):
            probs[idx] = np.nan
            err[idx] = np.nan
        else:
            probs[idx] = prob_vector[idx] / depth_vector[idx]
            # p_best = np.max(probs[idx])
            err[idx] = np.sqrt(probs[idx] * (1 - probs[idx]) / depth_vector[idx])

    return (probs, err, depth_vector)


# def determine_thresholds(
#     probs,
#     errors,
#     depth,
#     abs_min_prob=0,
#     abs_max_error=1e-3,
#     abs_min_depth=5,
#     dynamic=False,
# ):
#     logging.info(f"Dynamic thresholds is set as: {dynamic}")
#     # Consider datapoint only if it's prob is in the upper 25% of the data
#     inferred_prob_thresh = np.nanpercentile(probs, 75, axis=0)
#     absolute_prob_thresh = np.full_like(inferred_prob_thresh, abs_min_prob)
#     if dynamic:
#         prob_thresholds = np.maximum(
#             inferred_prob_thresh, absolute_prob_thresh
#         )  # will be used as a minimum value
#     else:
#         prob_thresholds = absolute_prob_thresh
#     logging.info(f"Inferred Probability Thresholds: {inferred_prob_thresh}")
#     logging.info(f"Absolute Probability Thresholds: {absolute_prob_thresh}")
#     logging.info(f"Selected Probability Thresholds (min): {prob_thresholds}")

#     # Consider datapoint only if it's error in the lower 90% of the data
#     inferred_err_thresh = np.nanpercentile(errors, 90, axis=0)
#     absolute_err_thresh = np.full_like(inferred_err_thresh, abs_max_error)
#     if dynamic:
#         err_thresholds = np.minimum(
#             inferred_err_thresh, absolute_err_thresh
#         )  # will be used as a maximum value
#     else:
#         err_thresholds = absolute_err_thresh
#     logging.info(f"Inferred Error Thresholds: {inferred_err_thresh}")
#     logging.info(f"Absolute Error Thresholds: {absolute_err_thresh}")
#     logging.info(f"Selected Error Thresholds (max): {err_thresholds}")

#     #     depth_threshold = max(np.nanpercentile(depth, 5), abs_min_depth)
#     depth_threshold = abs_min_depth
#     logging.info(f"Selected Depth Thresholds (min): {depth_threshold}")

#     cutoffs = {
#         "min_prob": prob_thresholds,
#         "max_error": err_thresholds,
#         "min_depth": depth_threshold,
#     }
#     return cutoffs


def prune_combined_profile(combined_profile, pruned_profile, min_depth=10):
    """
    Read the expanded profile:
        for each position:
            remove if 
                total depth < min_depth,
                prob == 1
            to reduce noise when visualizing,
            reduce base depths to 0 if:
                base depth < min_ind_depth

    # return position dict {pos = nucl = depth}; or
    return dataFrame [position, nucleotide, depth]
    """
    logging.info("Pruning combined profile ...")
    nucl = ["A", "C", "G", "T"]
    positions = defaultdict(int)
    with open(combined_profile, "r") as combined, open(pruned_profile, "w") as pruned:
        pruned.write("pos,nucl_depth,base,total_depth,nucl_prob,nucl_err\n")
        for lnum, line in enumerate(combined):
            if lnum == 0:
                # pos,depth,prob_A,prob_C,prob_G,prob_T,err_A,err_C,err_G,err_T
                continue
            (
                pos,
                depth,
                prob_A,
                prob_C,
                prob_G,
                prob_T,
                err_A,
                err_C,
                err_G,
                err_T,
            ) = line.split(",")
            # ignore low depth regions
            depth = float(depth)
            if depth < min_depth:
                continue
            prob_list = [float(prob_A), float(prob_C), float(prob_G), float(prob_T)]

            # ignore if probability is nan
            if np.isnan(prob_list).any():
                continue

            # ignore if probability is 1
            if max(prob_list) == 1:
                continue

            # err_list = list(map(float, [err_A, err_C, err_G, err_T]))
            err_list = [float(err_A), float(err_C), float(err_G), float(err_T)]
            # print(f"{prob_list}, {type(prob_list)}")
            for idx, nuc_prob in enumerate(prob_list):
                pos = int(pos)
                nucl_depth = depth * nuc_prob
                if nucl_depth > min_depth:
                    pruned.write(
                        f"{pos},{nucl_depth},{nucl[idx]},{depth},{nuc_prob},{err_list[idx]}\n"
                    )
                positions[pos] += 1
    logging.info(f"Kept {len(positions)} positions.")
    return positions


# def plot_on_genome(df, y_axis, x_axis="pos", title=None):
#     if df is None:
#         print("No dataframe found")
#         return

#     if title is None:
#         title = y_axis

#     agg = ds.Canvas().points(df, x_axis, y_axis, ds.count())
#     #     agg = ds.Canvas().points(df,x_axis,y_axis, ds.any())
#     return hd.shade(hv.Image(agg), cmap=list(Viridis256)).opts(
#         height=500, title=title, responsive=True
#     )


# def plot_distribution(df, column, title=None, min_value=0):
#     if title is None:
#         title = f"{column} Distribution"

#     return hv.Distribution(df[column][df[column] >= min_value]).opts(
#         height=500, title=title, responsive=True
#     )


def profile_differences(origin_pruned_pos, challenge_pruned_pos, min_overlap=1000):
    """
    check how many positions overlap.
        for the positions that overlap:
            keep positions whose values are different
    return two dictionaries of same length
        origin = {pos : nucl}
        challenge = {pos : nucl}
    """
    from collections import defaultdict

    origin = defaultdict(int)
    challenge = defaultdict(int)
    comparable_pos = origin_pruned_pos.keys() & challenge_pruned_pos.keys()
    len_comp_pos = len(comparable_pos)
    if len_comp_pos < min_overlap:
        logging.info(
            f"The two profiles do not overlap enough ({len_comp_pos}) to make a confident call for source. Stopping ..."
        )
        return

    for pos in comparable_pos:
        if origin_pruned_pos[pos] == challenge_pruned_pos[pos]:
            continue
        origin[pos] = origin_pruned_pos[pos]
        challenge[pos] = challenge_pruned_pos[pos]

    logging.info(f"Found {len(origin)} positions with differences")
    return (origin, challenge)


def save_combined_profile(probs, errors, depth, filename):
    with open(filename, "w") as file_handle:
        file_handle.write(
            "pos,depth,prob_A,prob_C,prob_G,prob_T,err_A,err_C,err_G,err_T\n"
        )
        for idx, arr in enumerate(depth):
            total_depth = depth[idx]
            acgt_prob = ",".join(map(str, probs[idx]))
            acgt_err = ",".join(map(str, errors[idx]))
            file_handle.write(f"{idx},{total_depth},{acgt_prob},{acgt_err}\n")
        logging.info(f"Wrote {idx + 1} lines to {filename}.")


def generate_pruned_profile(stat_files_list, genome_size, prefix=None):
    # Create temporary directory
    outdir_name = prefix
    cwd = os.getcwd()
    outdir = os.path.join(cwd, outdir_name)
    os.makedirs(outdir, exist_ok=True)

    # Origin
    ## Combine all samples
    probs, errors, depth = combine_profiles(stat_files_list, genome_size=genome_size)
    ## Determine thresholds for pruned
    # thresholds = determine_thresholds(probs, errors, depth)

    # TODO @sunit: Implement a way to not have to write out the profile.
    ## Save initial combined profile
    combined_profile_name = os.path.join(outdir, f"{prefix}.combinedProfile.csv")
    save_combined_profile(probs, errors, depth, combined_profile_name)

    ## Aggregate combined profile into one pruned profile
    pruned_profile_name = os.path.join(outdir, f"{prefix}.prunedProfile.csv")
    pruned = prune_combined_profile(
        combined_profile_name, pruned_profile_name, min_depth=10
    )

    # cleanup temporary
    # os.remove(combined_profile_name)
    # os.remove(pruned_profile_name)
    # os.rmdir(outdir)

    return pruned


if __name__ == "__main__":
    # logging.getLogger(__name__)
    logging.basicConfig(
        # filename=logfile,
        # filemode='w+',
        level=logging.INFO,
        format="%(asctime)s\t[%(levelname)s]:\t%(message)s",
    )

    origin_stat_files = sys.argv[1]
    # challenge_stat_files = sys.argv[2]
    prefix = sys.argv[2]

    # Origin
    logging.info("Estimating genome size ...")
    with open(origin_stat_files, "r") as origin:
        origin_stat_files_list = [line.strip() for line in origin]

    ## Genome Size
    genome_size, ncols = pd.read_csv(
        origin_stat_files_list[0], header=0, usecols=["pos"]
    ).shape
    logging.info(f"\t{genome_size} bases.")
    origin_pruned = generate_pruned_profile(origin_stat_files_list, genome_size, prefix)

    # Challenge
    # with open(challenge_stat_files, "r") as challenge:
    #     challenge_stat_files_list = [line.strip() for line in challenge]

    # challenge_pruned = generate_pruned_profile(
    #     challenge_stat_files_list, genome_size, f"{prefix}_challenge"
    # )

    # Isolate differences between two pruned profiles
    # origin_diff, challenge_diff = profile_differences(
    #     origin_pruned, challenge_pruned
    # )

    # Save profiles to disk
    # output = f"{prefix}.pruned.json"
    # logging.info(f"Saving pruned to disk as: {output}")
    # with open(output, "w") as out_handle:
    #     json.dump([origin_diff, challenge_diff], out_handle)
    logging.info("Huzzah!")
