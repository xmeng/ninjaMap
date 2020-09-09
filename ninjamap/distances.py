#!/usr/bin/env python3

import os
import sys
from itertools import combinations

import numpy as np
from scipy import stats
from scipy.spatial import distance

array_dir = sys.argv[1]
# array_dir = "/Users/sunit.jain/Research/Alice/in_vivo/Mouse_Backfill/bowtie2/individual_genome_comparison/pre_challenge/b_ovatus/vectors"


def load_clean_array(array_path):
    prob_array = np.load(array_path)
    # remove NANs from array
    return prob_array[~np.isnan(prob_array)]
    # return normalize_array(prob_array)


def normalize_array(prob_array):
    return stats.zscore(prob_array)


def compare(two_array_tuple):
    x, y = two_array_tuple
    # y_name = os.path.splitext(os.path.basename(y))[0]
    x_name = str(os.path.basename(x)).split("_")[0]
    y_name = str(os.path.basename(y)).split("_")[0]
    x_arr = load_clean_array(x)
    y_arr = load_clean_array(y)

    dist = euc_dist(x_arr, y_arr)
    ks_stat, p_value = ks_test(x_arr, y_arr)
    x_name, y_name = sorted([x_name, y_name])
    return (x_name, y_name, str(dist), str(ks_stat), str(p_value))


def euc_dist(x_arr, y_arr):
    return distance.euclidean(x_arr, y_arr)


def ks_test(x_arr, y_arr):
    # If the K-S statistic is small or the p-value is high,
    # then we cannot reject the hypothesis that the distributions of the two samples are the same.
    # x_arr = x_arr[x_arr != 0]
    # y_arr = y_arr[y_arr != 0]
    return stats.ks_2samp(x_arr, y_arr)


prob_array_files = [
    os.path.join(array_dir, file)
    for file in os.listdir(array_dir)
    if file.endswith(".ref_prob.npy")
]

depth_array_files = [
    os.path.join(array_dir, file)
    for file in os.listdir(array_dir)
    if file.endswith(".ref_depth.npy")
]

# arr1 = np.load(f'{array_dir}/Com3_vs_db_Clostridium-hathewayi-DSM-13479-MAF-NJ35.coord_sorted.npy')
# arr2 = np.load(f'{array_dir}/Pat1-3_vs_db_Clostridium-hathewayi-DSM-13479-MAF-NJ35.coord_sorted.npy')
# a = '/Users/sunit.jain/Research/Alice/in_vivo/Mouse_Backfill/bowtie2/individual_genome_comparison/pre_challenge/c_hath/vectors/SCV1-W5M6_vs_db_Clostridium-hathewayi-DSM-13479-MAF-NJ35.coord_sorted.npy'
# b = '/Users/sunit.jain/Research/Alice/in_vivo/Mouse_Backfill/bowtie2/individual_genome_comparison/pre_challenge/c_hath/vectors/SCV1-W3M6_vs_db_Clostridium-hathewayi-DSM-13479-MAF-NJ35.coord_sorted.npy'
# dist(a,b)

# load arrays, two at a time.
pairs = combinations(prob_array_files, 2)
# compute distance.
pairwise = [compare(pair) for pair in pairs]
with open(os.path.join(array_dir, "distances.ref_prob.tsv"), "w") as f:
    f.write("x\ty\teuc_dist\tks_stat\tks_pvalue_low_same\n")
    for p in pairwise:
        line = "\t".join(p)
        f.write(f"{line}\n")

# load arrays, two at a time.
pairs = combinations(depth_array_files, 2)
# compute distance.
pairwise = [compare(pair) for pair in pairs]
with open(os.path.join(array_dir, "distances.ref_depth.tsv"), "w") as f:
    f.write("x\ty\teuc_dist\tks_stat\tks_pvalue_low_same\n")
    for p in pairwise:
        line = "\t".join(p)
        f.write(f"{line}\n")
