#!/usr/bin/env python3

import concurrent.futures
import os
import sys
from collections import Counter, defaultdict
from functools import partial

import numpy as np
from scipy import stats
import pysam
from Bio import SeqIO


class Contigs:
    genome_size = int(0)
    num_contigs = int(0)

    def __init__(self, record):
        self.name, self.sequence = record.id, str(record.seq)
        self.length = int(len(self.sequence))
        self.ref_prob_array = np.empty(self.length)
        self.ref_depth_array = np.empty(self.length)
        self.ref_covered_bases = float(0)
        self.total_depth = float(0)
        self.ref_depth = float(0)

        Contigs.num_contigs += 1
        Contigs.genome_size += self.length

    def __hash__(self):
        return hash(str(self.sequence))

    def __eq__(self, other):
        return self.name == other.name

    def __ne__(self, other):
        # Not strictly necessary, but to avoid having both
        # x==y and x!=y being True at the same time
        return not (self == other)

    def weighted_ref_depth(self):
        if self.ref_covered_bases == 0:
            return 0

        return Contigs.weighted_cal(self.ref_depth, self.ref_covered_bases, self.length)

    @staticmethod
    def weighted_cal(success, sample_size, total):
        return (success / float(sample_size)) * float(total)

    @classmethod
    def read_fasta(cls, fasta_file):
        return [
            Contigs(fasta_record) for fasta_record in SeqIO.parse(fasta_file, "fasta")
        ]


# def parallel_calculate_coverage(bamfile_name, contigs, min_qual, min_pid, min_paln, cores):
#     with concurrent.futures.ThreadPoolExecutor(max_workers=cores) as executor:
#         future = [executor.submit(calculate_coverage_per_contig, bamfile_name, min_qual, min_pid, min_paln, contig) for contig in contigs]
#         for f in concurrent.futures.as_completed(future):
#             result = f.result()
#         ## When you don't wish to return anything, but monitor process for exceptions.
#         jobs = concurrent.futures.wait(future, return_when='FIRST_EXCEPTION')
#         if len(jobs.not_done) > 0:
#             print(f'Some ({len(jobs.not_done)}) exceptions occurred')
#             for job in jobs.not_done:
#                 print(f'{job.exception()}')
#             sys.exit(1)


def calculate_coverage(bamfile_name, contigs, min_qual, min_pid, min_paln, cores):
    [
        calculate_coverage_per_contig(bamfile_name, min_qual, min_pid, min_paln, contig)
        for contig in contigs
    ]


def calculate_coverage_per_contig(
    bamfile_name, minQual, thresh_pid, thresh_paln, contig
):
    with pysam.AlignmentFile(bamfile_name, mode="rb") as bamfile:
        check_read = partial(keep_read, min_pid=thresh_pid, min_paln=thresh_paln)
        # https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignmentFile.count_coverage
        cov_counts = bamfile.count_coverage(
            contig.name,
            start=0,
            end=contig.length,
            read_callback=check_read,
            # quality_threshold is the minimum quality score (in phred) a base has to reach to be counted.
            quality_threshold=minQual,
        )

        for i, ref_allele in enumerate(contig.sequence):
            base_cov = defaultdict(int)
            ref_prob = 0
            ref_depth = 0
            total_depth = 0
            ref_allele = str(ref_allele).upper()
            # if refseq has N, the prob is 'NA'
            if ref_allele in ["A", "T", "C", "G"]:
                base_cov["A"] = cov_counts[0][i]
                base_cov["C"] = cov_counts[1][i]
                base_cov["G"] = cov_counts[2][i]
                base_cov["T"] = cov_counts[3][i]
                total_depth = (
                    base_cov["A"] + base_cov["C"] + base_cov["G"] + base_cov["T"]
                )
                # prob 0 if no depth
                if total_depth > 0:
                    contig.total_depth += total_depth
                    ref_depth = base_cov[ref_allele]
                    ref_prob = float(ref_depth) / float(total_depth)
                    if ref_depth > 0:
                        contig.ref_covered_bases += 1
                        contig.ref_depth += ref_depth
            else:
                ref_prob = np.nan
                ref_depth = np.nan

            contig.ref_prob_array[i] = ref_prob
            contig.ref_depth_array[i] = ref_depth
    return contig


def keep_read(aln, min_pid, min_paln):
    edit_dist = dict(aln.tags)["NM"]
    aln_len = aln.query_alignment_length
    read_length = aln.infer_read_length()

    pid = 100 * (aln_len - edit_dist) / float(aln_len)
    p_aln_len = aln_len * 100 / float(read_length)

    return (min_pid <= pid <= 100) and (min_paln <= p_aln_len <= 100)


def distance_from_reference(prob_array):
    # remove NANs from array
    prob_array = remove_nan(prob_array)
    # create a unit array of same length representing reference.
    ref_array = np.ones_like(prob_array)

    return (
        root_mean_squared_error(ref_array, prob_array),
        cos_sim(ref_array, prob_array),
    )


def root_mean_squared_error(ref_array, prob_array):
    """ Calculate mean squared error given two arrays of same length """
    if len(ref_array) == len(prob_array):
        return np.sqrt(np.square(np.subtract(ref_array, prob_array)).mean())


def cos_sim(ref_array, prob_array):
    """calculate cosine distance between reference and given."""
    if len(ref_array) == len(prob_array):
        return np.dot(ref_array, prob_array) / (
            np.linalg.norm(ref_array) * np.linalg.norm(prob_array)
        )


def remove_nan(prob_array):
    return prob_array[~np.isnan(prob_array)]


if __name__ == "__main__":
    bam_file = sys.argv[1]
    genome = sys.argv[2]

    cores = 4
    min_qual = 20  # quality_threshold is the minimum quality score (in phred) a base has to reach to be counted.
    min_pid = 99
    min_paln = 70

    output_dir = f"q{min_qual}_id{min_pid}_aln{min_paln}_vectors"
    os.makedirs(output_dir, exist_ok=True)

    file_name = os.path.splitext(os.path.basename(bam_file))[0]
    prefix = f"{output_dir}/{file_name}.q{min_qual}_id{min_pid}_aln{min_paln}"

    ref_prob_array_name = f"{prefix}.ref_prob.npy"
    ref_depth_array_name = f"{prefix}.ref_depth.npy"

    contigs = Contigs.read_fasta(genome)

    calculate_coverage(bam_file, contigs, min_qual, min_pid, min_paln, cores)
    # parallel_calculate_coverage(bam_file, contigs, min_qual, min_pid, min_paln, cores)

    depth = (
        np.sum([contig.weighted_ref_depth() for contig in contigs])
        / Contigs.genome_size
    )
    coverage = (
        100
        * np.sum([contig.ref_covered_bases for contig in contigs])
        / Contigs.genome_size
    )

    # combine all prob arrays in order of contigs
    genome_ref_prob = np.concatenate(
        [remove_nan(contig.ref_prob_array) for contig in contigs]
    )
    np.save(ref_prob_array_name, genome_ref_prob)

    genome_depth_prob = np.concatenate(
        [remove_nan(contig.ref_depth_array) for contig in contigs]
    )
    np.save(ref_depth_array_name, genome_depth_prob)

    depth_ent = stats.entropy(genome_depth_prob)
    rmse, sim = distance_from_reference(genome_ref_prob)
    print(f"{file_name}\t{depth}\t{coverage}\t{depth_ent}\t{rmse}\t{sim}")

# Calculate this vector for each sample.
# Option1: Cosine Similarity for each sample vs Reference.
# Option2: Eucledian distance for each pair. Later, could be reduced to Weekly sample vs Comm3/Pat1.
# Option3: Gaussian Mixture models (Expectation Maximization) on all vectors.
# Option4: Covariation models. Subset vector for genes, generate set of genes that covary in samples from W1-4, Comm1-3 and Pat#-1-3. If the same set covaries by same proportion in W5-8, it is one of the SC strains, else acquired from Challenge. Use saline challenge as control.
