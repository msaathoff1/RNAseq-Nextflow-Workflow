#!/usr/bin/env python3
# %%

import logging
import argparse
import gffutils
import os
import sqlite3
import math

# Set up logger for error reporting of the database
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler('logFile.log')
fh.setLevel(logging.INFO)
fh.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
logger.addHandler(fh)


# This function generates a file containing rpk values for each count file that is input. It also returns a sum
# of all the rpk values in the output file. --> [total_rpk, dictionary{gene_id: rpk_values} ]
# db refers to the database connection.

def count_file_to_dict(count_file):
    count_file_dict = {}
    additionalInfo = ['__no_feature', '__ambiguous', '__too_low_aQual', '__not_aligned', '__alignment_not_unique']
    try:
        with open(count_file, 'r') as f:
            for line in f:
                rowValues = line.strip().split('\t')
                gene_id = str(rowValues[0])
                count_value = float(rowValues[1])

                # Omit the last five rows of the count file that contain special information.
                if gene_id in additionalInfo:
                    continue
                count_file_dict[gene_id] = count_value
    except IOError:
        logger.error('Cannot unpack count file.\n')
        raise SystemExit('Cannot unpack count file.')
    return count_file_dict


def get_rpk(count_dict, db):
    rpk_sum = 0
    rpk_dict = {}
    gene_id_list = count_dict.keys()
    for gene_id in gene_id_list:
        gene_length = 0
        # Check to see if the gene_id is actually in the gff file
        try:
            candidates = db[gene_id]

        except gffutils.exceptions.FeatureNotFoundError:
            logger.warning(
                'Gene ID {} could not be found in the gff file. If this is unintended, check gff '
                'file and count data.\n'.format(gene_id))
            continue

        # Actually do the math
        else:
            try:
                count = count_dict[gene_id]
            except KeyError:
                logger.error("Gene {} is not in the dictionary, check count files.\n".format(gene_id))
                raise SystemExit("Gene {} is not in the dictionary, check count files.\n".format(gene_id))

            # find the gene length by summing all the non-overlapping exons for a gene.
            for exon in db.children(gene_id, level=2, featuretype='exon'):
                try:
                    exon_length = exon.end - exon.start + 1
                except ValueError:
                    logger.warning('Could not identify exon start and end for Gene ID {}\n'.format(gene_id))
                    continue

                gene_length += exon_length
                RPK = count / (gene_length / 1000)

            rpk_dict[gene_id] = RPK
            rpk_sum = rpk_sum + RPK

    return rpk_sum, rpk_dict


# This function takes rpk files and generates a tpm file based on that file and a provided scaling value.

def get_tpm(rpk_dict, scaling, count_type, direction, SampleName):
    try:
        with open(str(SampleName) + '.' + direction + "." + count_type + '.tpm.txt', 'w') as ftpm:
            tpm_sum = 0
            for gene_id, rpk_value in rpk_dict.items():
                try:
                    tpm_value = rpk_value / scaling
                    tpm_sum += tpm_value
                except ZeroDivisionError:
                    logger.error('Check to ensure scaling value is NOT zero!\n')
                    raise SystemExit('Check to ensure scaling value is NOT zero!')

                try:
                    ftpm.write("{}\t{}\n".format(gene_id, str(tpm_value)))
                except IOError:
                    logger.error('Could not write to output file.\n')
                    raise SystemExit('Could not write to output file.')

    except IOError:
        logger.error('Output TPM File could not be generated.\n')
        raise SystemExit('Output TPM File could not be generated.')
    return tpm_sum


def __main__():
    # Set up command line input
    parser = argparse.ArgumentParser(description='Calculate TPM values')
    parser.add_argument('--GffFile', required=True, help='GFF file')
    parser.add_argument('--UniqueFWD', required=True, help='HTSeq-Count file for unique_fwd reads')
    parser.add_argument('--UniqueREV', help='HTSeq-Count file for unique_rev reads')
    parser.add_argument('--NonuniqueFWD', required=True, help='HTSeq-Count file for nonunique_fwd reads')
    parser.add_argument('--NonuniqueREV', help='HTSeq-Count file for nonunique_rev reads')
    parser.add_argument('--Stranded', action='store_true', help='If true will look for 4 files (fwd/rev) instead of 2.')
    parser.add_argument('--SampleName', required=True, help="The sample name of the associated files - will be file "
                                                            "name prefix for output.")
    args = parser.parse_args()

    # Name the database
    db = args.GffFile.replace('.gff', '.db')

    # Create Database
    if not os.path.isfile(db):
        logger.info('Creating database {}...\n'.format(db))
        try:
            db = gffutils.create_db(args.GffFile, dbfn=db, force=True, keep_order=True, merge_strategy='merge',
                                    sort_attribute_values=True)
        except sqlite3.OperationalError:
            logger.error('Cannot create database {}. Check gff file.\n'.format(db))
            print('Cannot create database {}. Check gff file.\n'.format(db))
            raise SystemExit(1)
    else:
        logger.info('Connecting to existing database {}...\n'.format(db))
        try:
            db = gffutils.FeatureDB(db, keep_order=True)
        except ValueError:
            logger.error('DB {} could not be read\n'.format(db))
            print('DB {} could not be read\n'.format(db))
            raise SystemExit(1)

    if args.Stranded:
        # Check to make sure there are 4 file inputs for a stranded dataset
        if not args.UniqueREV or not args.NonuniqueREV:
            raise SystemExit('If the dataset is stranded, four files must be provided: UniqueREV, NonuniqueREV, '
                             'UniqueFWD, and NonuniqueFWD')

        # Store count values as dictionaries for use
        UniqueFWD_count_dict = count_file_to_dict(args.UniqueFWD)
        UniqueREV_count_dict = count_file_to_dict(args.UniqueREV)
        NonuniqueFWD_count_dict = count_file_to_dict(args.NonuniqueFWD)
        NonuniqueREV_count_dict = count_file_to_dict(args.NonuniqueREV)

        # Generate TRUE Nonunique counts (currently nonunique counts are unique+nonunique)
        true_nonunique_fwd_count_dict = {key: NonuniqueFWD_count_dict[key] - UniqueFWD_count_dict.get(key, 0) for key in
                                         NonuniqueFWD_count_dict}
        true_nonunique_rev_count_dict = {key: NonuniqueREV_count_dict[key] - UniqueREV_count_dict.get(key, 0) for key in
                                         NonuniqueREV_count_dict}

        # Get RPK info - function returns two values [total_rpk, dictionary{gene_id: rpk_values} ]
        UniqueFWD_rpk = get_rpk(UniqueFWD_count_dict, db)
        UniqueREV_rpk = get_rpk(UniqueREV_count_dict, db)
        NonuniqueFWD_rpk = get_rpk(true_nonunique_fwd_count_dict, db)
        NonuniqueREV_rpk = get_rpk(true_nonunique_rev_count_dict, db)

        # Calculate scaling value
        total_rpk_sum = UniqueFWD_rpk[0] + UniqueREV_rpk[0] + NonuniqueFWD_rpk[0] + NonuniqueREV_rpk[0]
        scaling_value = total_rpk_sum / 1000000

        # Calculate tpm value
        unique_fwd_tpm = get_tpm(UniqueFWD_rpk[1], scaling_value, 'unique', 'fwd', args.SampleName)
        unique_rev_tpm = get_tpm(UniqueREV_rpk[1], scaling_value, 'unique', 'rev', args.SampleName)
        nonunique_fwd_tpm = get_tpm(NonuniqueFWD_rpk[1], scaling_value, 'nonunique', 'fwd', args.SampleName)
        nonunique_rev_tpm = get_tpm(NonuniqueREV_rpk[1], scaling_value, 'nonunique', 'rev', args.SampleName)

        # Double check all is good
        total_tpm_sum = math.floor(unique_fwd_tpm + unique_rev_tpm + nonunique_fwd_tpm + nonunique_rev_tpm)
        if total_tpm_sum != 1000000:
            raise SystemExit("TPM total for {} is not one million. Check everything's good.".format(args.SampleName))

    else:
        # Store count values as dictionaries for use
        UniqueFWD_count_dict = count_file_to_dict(args.UniqueFWD)
        NonuniqueFWD_count_dict = count_file_to_dict(args.NonuniqueFWD)

        # Generate TRUE Nonunique counts (currently nonunique counts are unique+nonunique)
        true_nonunique_fwd_count_dict = {key: NonuniqueFWD_count_dict[key] - UniqueFWD_count_dict.get(key, 0) for key in
                                         NonuniqueFWD_count_dict}

        # Get RPK info - function returns two values [total_rpk, dictionary{gene_id: rpk_values} ]
        UniqueFWD_rpk = get_rpk(UniqueFWD_count_dict, db)
        NonuniqueFWD_rpk = get_rpk(true_nonunique_fwd_count_dict, db)

        # Calculate scaling value
        total_rpk_sum = UniqueFWD_rpk[0] + NonuniqueFWD_rpk[0]
        scaling_value = total_rpk_sum / 1000000

        # Calculate tpm value
        unique_fwd_tpm = get_tpm(UniqueFWD_rpk[1], scaling_value, 'unique', 'fwd', args.SampleName)
        nonunique_fwd_tpm = get_tpm(NonuniqueFWD_rpk[1], scaling_value, 'nonunique', 'fwd', args.SampleName)

        # Double check all is good
        total_tpm_sum = math.floor(unique_fwd_tpm + nonunique_fwd_tpm)
        if total_tpm_sum != 1000000:
            raise SystemExit("TPM total for {} is not one million. Check everything's good.".format(args.SampleName))


if __name__ == "__main__": __main__()


