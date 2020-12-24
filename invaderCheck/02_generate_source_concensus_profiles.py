#!/usr/bin/env python3

import os, sys
import numpy as np
import json
import logging
import pandas as pd

logging.getLogger(__name__).addHandler(logging.NullHandler())


def generate_profile(
    list_of_dfs, genome_size,
):
    """
    Across all Samples
        For each position,
            calculate probability for A, T, G and C
            calculate std error as sqrt[P * (1-P)/N], where P is probability for Nuc and N is the depth at this position.
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
        # Not reading positions because these are not reliable when multiple contigs are present.
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
            err[idx] = np.sqrt(probs[idx] * (1 - probs[idx]) / depth_vector[idx])

    return (probs, err, depth_vector)


def determine_thresholds(
    probs,
    errors,
    depth,
    abs_min_prob=0.9,
    abs_max_error=1e-3,
    abs_min_depth=5,
    dynamic=True,
):
    logging.info(f"Dynamic thresholds is set as: {dynamic}")
    # Consider datapoint only if it's prob is in the upper 25% of the data
    inferred_prob_thresh = np.nanpercentile(probs, 75, axis=0)
    absolute_prob_thresh = np.full_like(inferred_prob_thresh, abs_min_prob)
    if dynamic:
        prob_thresholds = np.maximum(
            inferred_prob_thresh, absolute_prob_thresh
        )  # will be used as a minimum value
    else:
        prob_thresholds = absolute_prob_thresh
    logging.info(f"Inferred Probability Thresholds: {inferred_prob_thresh}")
    logging.info(f"Absolute Probability Thresholds: {absolute_prob_thresh}")
    logging.info(f"Selected Probability Thresholds (min): {prob_thresholds}")

    # Consider datapoint only if it's error in the lower 90% of the data
    inferred_err_thresh = np.nanpercentile(errors, 90, axis=0)
    absolute_err_thresh = np.full_like(inferred_err_thresh, abs_max_error)
    if dynamic:
        err_thresholds = np.minimum(
            inferred_err_thresh, absolute_err_thresh
        )  # will be used as a maximum value
    else:
        err_thresholds = absolute_err_thresh
    logging.info(f"Inferred Error Thresholds: {inferred_err_thresh}")
    logging.info(f"Absolute Error Thresholds: {absolute_err_thresh}")
    logging.info(f"Selected Error Thresholds (max): {err_thresholds}")

    #     depth_threshold = max(np.nanpercentile(depth, 5), abs_min_depth)
    depth_threshold = abs_min_depth
    logging.info(f"Selected Depth Thresholds (min): {depth_threshold}")

    cutoffs = {
        "min_prob": prob_thresholds,
        "max_error": err_thresholds,
        "min_depth": depth_threshold,
    }
    return cutoffs


def get_concensus_profile(
    combined_profile, consensus_profile, min_prob, max_error, min_depth
):
    """
    Read the expanded profile and:
        generate concensus for each position
        write concensus to concensus profile file
        return position dict {pos = nucl}
    """
    from collections import defaultdict

    nucl = ["A", "C", "G", "T"]
    positions = defaultdict(int)
    with open(combined_profile, "r") as init, open(consensus_profile, "w") as con:
        con.write(f"pos,depth,base,prob,err\n")
        for lnum, line in enumerate(init):
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
            if float(depth) < min_depth:
                continue

            prob_list = [prob_A, prob_C, prob_G, prob_T]
            prob_list = list(map(float, prob_list))

            # ignore if probability is nan
            if np.isnan(prob_list).any():
                continue

            err_list = [err_A, err_C, err_G, err_T]
            err_list = list(map(float, err_list))

            # get the index of the highest base prob
            concensus_prob = max(prob_list)
            concensus_idx = prob_list.index(concensus_prob)

            concensus_nucl = nucl[concensus_idx]
            concensus_err = err_list[concensus_idx]
            concensus_depth = float(depth) * concensus_prob
            # ignore ambiguous calls
            concensus_min_prob = min_prob[concensus_idx]
            concensus_max_error = max_error[concensus_idx]
            if (concensus_prob < concensus_min_prob) or (
                concensus_err > concensus_max_error
            ):
                continue
            pos = int(pos)
            positions[pos] = concensus_nucl
            #       Output: postion, Base, prob, depth, err
            con.write(
                f"{pos},{concensus_depth},{concensus_nucl},{concensus_prob},{concensus_err}\n"
            )
    logging.info(f"Kept {len(positions)} positions.")
    return positions


def profile_differences(
    origin_concensus_pos, challenge_concensus_pos, min_overlap=1000
):
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
    comparable_pos = origin_concensus_pos.keys() & challenge_concensus_pos.keys()
    len_comp_pos = len(comparable_pos)
    if len_comp_pos < min_overlap:
        logging.info(
            f"The two profiles do not overlap enough ({len_comp_pos}) to make a confident call for source. Stopping ..."
        )
        return

    for pos in comparable_pos:
        if origin_concensus_pos[pos] == challenge_concensus_pos[pos]:
            continue
        origin[pos] = origin_concensus_pos[pos]
        challenge[pos] = challenge_concensus_pos[pos]

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


def generate_concensus_profile(stat_files_list, genome_size, prefix=None):
    # Create temporary directory
    outdir_name = prefix
    cwd = os.getcwd()
    outdir = os.path.join(cwd, outdir_name)
    os.makedirs(outdir, exist_ok=True)

    # Origin
    ## Combine all samples
    probs, errors, depth = generate_profile(stat_files_list, genome_size=genome_size)
    ## Determine thresholds for concensus
    thresholds = determine_thresholds(probs, errors, depth)

    # TODO @sunit: Implement a way to not have to write out the profile.
    ## Save initial combined profile
    combined_profile_name = os.path.join(outdir, f"{prefix}.combinedProfile.csv")
    save_combined_profile(probs, errors, depth, combined_profile_name)

    ## Aggregate combined profile into one concensus profile
    consensus_profile_name = os.path.join(outdir, f"{prefix}.concensusProfile.csv")
    concensus = get_concensus_profile(
        combined_profile_name, consensus_profile_name, **thresholds
    )

    # cleanup temporary
    os.remove(combined_profile_name)
    os.remove(consensus_profile_name)
    os.rmdir(outdir)

    return concensus


if __name__ == "__main__":
    # logging.getLogger(__name__)
    logging.basicConfig(
        # filename=logfile,
        # filemode='w+',
        level=logging.INFO,
        format="%(asctime)s\t[%(levelname)s]:\t%(message)s",
    )

    origin_stat_files = sys.argv[1]
    challenge_stat_files = sys.argv[2]
    prefix = sys.argv[3]

    # Origin
    with open(origin_stat_files, "r") as origin:
        origin_stat_files_list = [line.strip() for line in origin]

    ## Genome Size
    genome_size, ncols = pd.read_csv(
        origin_stat_files_list[0], header=0, usecols=["pos"]
    ).shape

    origin_consensus = generate_concensus_profile(
        origin_stat_files_list, genome_size, f"{prefix}_origin"
    )

    # Challenge
    with open(challenge_stat_files, "r") as challenge:
        challenge_stat_files_list = [line.strip() for line in challenge]

    challenge_consensus = generate_concensus_profile(
        challenge_stat_files_list, genome_size, f"{prefix}_challenge"
    )

    # Isolate differences between two concensus profiles
    origin_diff, challenge_diff = profile_differences(
        origin_consensus, challenge_consensus
    )

    # Save profiles to disk
    output = f"{prefix}.concensus.json"
    logging.info(f"Saving concensus to disk as: {output}")
    with open(output, "w") as out_handle:
        json.dump([origin_diff, challenge_diff], out_handle)
    logging.info("Huzzah!")
