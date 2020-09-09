#!/usr/bin/env python3

import os
import subprocess
import sys
from random import random

import pysam
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator


class Contigs:
    genome_size = int(0)
    num_contigs = int(0)

    def __init__(self, record):
        self.name, self.sequence = record.id, str(record.seq)
        self.length = int(len(self.sequence))

        Contigs.num_contigs += 1
        Contigs.genome_size += self.length

    @classmethod
    def read_fasta(cls, fasta_file):
        return [
            Contigs(fasta_record) for fasta_record in SeqIO.parse(fasta_file, "fasta")
        ]


def sort_and_index(file_name, by="coord", cores=4, memPerCore=2):
    """
    Sorts and indexes a bam file by coordinates.
    """
    sorted_name = None
    if by.lower() == "coord":
        sorted_name = file_name.replace(".bam", "") + ".sortedByCoord.bam"
        if not os.path.exists(sorted_name):
            sort_command = f"samtools sort -@ {cores} -m {memPerCore}G -o {sorted_name} {file_name}"
            subprocess.run(sort_command, shell=True, check=True)

        subprocess.run(
            f"samtools index -@ {cores} {sorted_name}", shell=True, check=True
        )
    elif by.lower() == "name":
        sorted_name = file_name.replace(".bam", "") + ".sortedByName.bam"
        if not os.path.exists(sorted_name):
            sort_command = f"samtools sort -n -@ {cores} -m {memPerCore}G -o {sorted_name} {file_name}"
            subprocess.run(sort_command, shell=True, check=True)
    else:
        raise Exception("Bam file can only be sorted by 'coord' or 'name'.")

    # os.remove(file_name)
    return sorted_name


def subsample_fastq(fastq_file, num_reads, output=None):
    with open(fastq_file) as f:
        for num_lines, l in enumerate(f):
            pass
    num_lines += 1
    total_reads = num_lines / 4

    if num_reads >= total_reads:
        return (None, total_reads)

    fraction = num_reads / total_reads

    selected_reads = 0
    with open(fastq_file, "r") as fastq, open(output, "w") as out:
        for title, seq, qual in FastqGeneralIterator(fastq):
            if random() <= fraction:
                selected_reads += 1
                out.write(f"{title}\n")

    return (selected_reads, total_reads)


def subsample_aligned_reads(bamfile_name, num_reads, output_file, cores=4):
    # sort bam by name
    print(f"\tSorting {bamfile_name} by name")
    name_sorted_bamfile = sort_and_index(bamfile_name, by="name", cores=cores)
    # run samtools fastq to get list of reads
    print(f"\tConverting {bamfile_name} to fastq, using {cores} cores")
    subprocess.run(
        f"samtools fastq -@ {cores} -1 {bamfile_name}.fwd.fastq -2 /dev/null -0 /dev/null -s /dev/null -n {name_sorted_bamfile}",
        shell=True,
        check=True,
    )
    # subsample fastq to desired number of reads and write to file
    print(f"\tSubsampling approx {num_reads} from fastq.")
    num_selected_reads, total_reads = subsample_fastq(
        f"{bamfile_name}.fwd.fastq", num_reads, output_file
    )
    os.remove(f"{bamfile_name}.fwd.fastq")
    return num_selected_reads, total_reads


if __name__ == "__main__":
    bamfile_name = sys.argv[1]
    genome = sys.argv[2]
    downsample_depth = int(sys.argv[3])
    at_coverage_perc = float(sys.argv[4])
    cores = 10
    average_readpair_size = 280

    contigs = Contigs.read_fasta(genome)
    pairs_for_1x = (
        Contigs.genome_size * at_coverage_perc / 100
    ) / average_readpair_size
    pairs_for_desired_depth = int(pairs_for_1x * downsample_depth)

    prefix = os.path.splitext(os.path.basename(bamfile_name))[0]

    output_file = f"{prefix}_depth{downsample_depth}x_at{at_coverage_perc}Pcov.list"
    print(f"Sampling for about {pairs_for_desired_depth} reads")
    num_selected_reads, total_reads = subsample_aligned_reads(
        bamfile_name, pairs_for_desired_depth, output_file, cores=cores
    )

    if num_selected_reads is None:
        print(
            f"[ERROR] Number of reads required ({pairs_for_desired_depth})",
            f"to meet depth criteria ({downsample_depth}x) at {at_coverage_perc}% coverage,",
            f"is greater than the total number of reads ({total_reads}).",
        )
        sys.exit(1)
    else:
        print(
            f"Huzzah! Generated a list of ({num_selected_reads}/{total_reads})",
            f"which is close to our initial estimate of {pairs_for_desired_depth}",
            f"to meet depth criteria ({downsample_depth}x) at {at_coverage_perc}% coverage",
        )
