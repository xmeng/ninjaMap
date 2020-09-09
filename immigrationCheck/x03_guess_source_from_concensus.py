#!/usr/bin/env python3
import os, sys
import pandas as pd
import json
import logging
import numpy as np


logging.getLogger(__name__).addHandler(logging.NullHandler())


def compare_query_profile(list_of_files, origin_pos, challenge_pos, min_depth=2):
    num_files = len(list_of_files)
    logging.info(f"Found {num_files} file(s)")

    num_positions = len(origin_pos.keys())
    results = list()
    # per Sample
    for file_num, df_file in enumerate(list_of_files):
        df = dict()
        filename = os.path.basename(df_file).split(".")[0]
        logging.info(f"\t[{file_num+1}/{num_files}] Processing {filename} ...")
        origin_prob = 0
        challenge_prob = 0
        other_prob = 0
        no_depth = 0
        cum_depth = 0
        data_frame = pd.read_csv(
            df_file, header=0, usecols=["A", "C", "G", "T", "total_depth"]
        )
        nrows, ncols = data_frame.shape
        perc_cov = 100 * num_positions / nrows
        logging.info(
            f"\tAggregating data from {num_positions} bases out of {nrows} ({perc_cov}%) ..."
        )
        # per Position
        for row in data_frame.itertuples(index=True):
            pos = str(row.Index)
            if pos not in origin_pos:
                continue

            origin_nuc = origin_pos[str(pos)]
            challenge_nuc = challenge_pos[str(pos)]
            nuc_probs = {
                "A": row.A,
                "C": row.C,
                "G": row.G,
                "T": row.T,
            }
            #    raise Exception(f"Origin {origin_nuc} @ {pos} = {nuc_probs[origin_nuc]/row.total_depth};\tChallenge {challenge_nuc} @ {pos} = {nuc_probs[challenge_nuc]/row.total_depth};")
            #    Origin G @ 145 = 33;	Challenge C @ 145 = 0;
            if row.total_depth > min_depth and not (np.isnan(row.total_depth)):
                origin_prob += nuc_probs[origin_nuc] / row.total_depth
                challenge_prob += nuc_probs[challenge_nuc] / row.total_depth
                other_prob += (
                    1
                    - (nuc_probs[origin_nuc] / row.total_depth)
                    - (nuc_probs[challenge_nuc] / row.total_depth)
                )
                cum_depth += row.total_depth
            else:
                no_depth += 1
        df = {
            "Sample": filename,
            "Origin": 100 * origin_prob / num_positions,
            "Challenge": 100 * challenge_prob / num_positions,
            "Other": 100 * other_prob / num_positions,
            "No_Info": 100 * no_depth / num_positions,
            "Bases_Compared": num_positions,
            "Average_Depth": cum_depth / num_positions,
            "Percent_Coverage": perc_cov,
        }
        results.append(df)
    return pd.DataFrame(results)


def load_concensus_profile(concensus_json):
    with open(concensus_json, "r") as concensus:
        origin_pos, challenge_pos = json.load(concensus)
    # logging.info(f"Origin: {len(origin_pos)}\tChallenge:{len(challenge_pos)}")
    orLen = len(origin_pos)
    chaLen = len(challenge_pos)
    try:
        assert orLen == chaLen
    except AssertionError as err:
        logging.exception("Malformed concensus profile!")
        raise err
    return (origin_pos, challenge_pos)


if __name__ == "__main__":
    query_stat_files = sys.argv[1]
    concensus_json = sys.argv[2]
    output = sys.argv[3]

    logging.basicConfig(
        # filename=logfile,
        # filemode='w+',
        level=logging.INFO,
        format="%(asctime)s\t[%(levelname)s]:\t%(message)s",
    )

    origin_pos, challenge_pos = load_concensus_profile(concensus_json)

    with open(query_stat_files, "r") as origin:
        query_stat_files_list = [line.strip() for line in origin]
    df = compare_query_profile(query_stat_files_list, origin_pos, challenge_pos)
    logging.info(f"Writing output to: {output}")
    df.to_csv(output, index=False)
    logging.info("Huzzah!")
