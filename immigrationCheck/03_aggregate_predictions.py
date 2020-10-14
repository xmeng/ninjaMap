#!/usr/bin/env python3

import sys
import pandas as pd
import logging


def load_dataframe(file_path):
    return pd.read_csv(file_path, header=0)


def aggregate_predictions(pred_paths_file):
    with open(pred_paths_file, "r") as file_paths:
        list_of_dfs = [load_dataframe(file_path.strip()) for file_path in file_paths]

    return pd.concat(list_of_dfs).reset_index(drop=True)


def visualize_predictions(df):
    pass


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s\t[%(levelname)s]:\t%(message)s",
    )

    pred_paths_file = sys.argv[1]
    out_prefix = sys.argv[2]

    df = aggregate_predictions(pred_paths_file)
    df.to_csv(f"{out_prefix}.03_predictions.csv", index=False)

    # visualize_predictions(df)
