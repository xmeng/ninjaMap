#!/usr/bin/env python3
from concurrent.futures import ThreadPoolExecutor

import argparse
import gzip
import logging
import random
import sys

def usage():
    parser = argparse.ArgumentParser()
    parser.add_argument("fwd", help="input FASTQ filename")
    parser.add_argument("rev", help="input FASTQ filename")
    parser.add_argument("prefix", help="output FASTQ filename")
    parser.add_argument("-f", "--fraction", type=float, help="fraction of reads to sample")
    parser.add_argument("-n", "--number", type=int, help="number of reads to sample")
    parser.add_argument("-s", "--sample", type=int, help="number of output files to write", default=1)
    args = parser.parse_args()

    if args.fraction and args.number:
        sys.exit("give either a fraction or a number, not both")

    if not args.fraction and not args.number:
        sys.exit("you must give either a fraction or a number")
    
    return args

def count_records(fastq):
    with gzip.open(fastq) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def write_sampled_seq(fastq, output_files, output_sequence_sets, total_records, log):
    record_number = 0
    chunk_size = int(total_records/10)
    with gzip.open(fastq, 'rt') as fq:
        for line1 in fq:
            line2 = next(fq)
            line3 = next(fq)
            line4 = next(fq)
            # foreach output file
            for i, output in enumerate(output_files):
                # if this record num is in this output set; write
                if record_number in output_sequence_sets[i]:
                    output.write(line1)
                    output.write(line2)
                    output.write(line3)
                    output.write(line4)
            record_number += 1
            if record_number % chunk_size == 0:
                log.info(f"{str((record_number / total_records) * 100)}% of {fastq} done")
    for output in output_files:
        output.close()

def run_io_tasks_in_parallel(tasks):
    with ThreadPoolExecutor() as executor:
        running_tasks = [executor.submit(task) for task in tasks]
        for running_task in running_tasks:
            running_task.result()

if __name__=='__main__':
    args = usage()
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s\t[%(levelname)s]:\t%(message)s')
    
    logging.info("counting records in fwd file ...")
    total_fwd_records = count_records(args.fwd)
    logging.info(f'\tfound {total_fwd_records}')

    logging.info("counting records in rev file ...")
    total_rev_records = count_records(args.rev)
    logging.info(f'\tfound {total_rev_records}')

    total_records = 0

    if total_fwd_records == total_rev_records:
        logging.info(f"both files have the same number of records ... {total_fwd_records}")
        total_records = total_fwd_records
    else:
        logging.critical(f'[FATAL] Paired files have different number of records: FWD={total_fwd_records} reads and REV={total_rev_records} reads')
        sys.exit(1)

    if args.fraction:
        args.number = int(total_records * args.fraction)

    logging.info("sampling " + str(args.number) + " out of " + str(total_records) + " records")

    output_fwd_files = []
    output_rev_files = []
    output_sequence_sets = []
    for i in range(args.sample):
        output_fwd_files.append(gzip.open(f'{args.prefix}.{str(i)}.fwd.fastq.gz', "wt"))
        output_rev_files.append(gzip.open(f'{args.prefix}.{str(i)}.rev.fastq.gz', "wt"))
        output_sequence_sets.append(set(random.sample(range(total_records + 1), args.number)))

    # logging.info('creating subsampled FWD file(s) ...')
    # write_sampled_seq(args.fwd, output_fwd_files, output_sequence_sets, total_records, logging)

    # logging.info('creating subsampled REV file(s) ...')
    # write_sampled_seq(args.rev, output_rev_files, output_sequence_sets, total_records, logging)
    # logging.info("done!")

    logging.info('creating subsampled FWD and REV file(s) ...')
    run_io_tasks_in_parallel([
        lambda: {'FWD': write_sampled_seq(args.fwd, output_fwd_files, output_sequence_sets, total_records, logging)},
        lambda: {'REV': write_sampled_seq(args.rev, output_rev_files, output_sequence_sets, total_records, logging)}
    ])
    
    
