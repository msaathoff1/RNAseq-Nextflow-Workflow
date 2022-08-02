#!/usr/bin/env python3
# %%

import logging
import argparse
import gffutils
import os
import sqlite3

# %%

import logging
import argparse
import gffutils
import os
import sqlite3

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

def get_rpk(count_file, db):
    rpk_sum = 0
    rpk_dict = {}
    additionalInfo = ['__no_feature', '__ambiguous', '__too_low_aQual', '__not_aligned', '__alignment_not_unique']
    try:
        with open(count_file, 'r') as f:
            for line in f:
                gene_length = 0
                rowValues = line.strip().split('\t')
                gene_id = str(rowValues[0])

                # Omit the last five rows of the count file that contain special information.
                if gene_id in additionalInfo:
                    continue

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
                        count = int(rowValues[1])
                    except ValueError:
                        logger.error('Count value could not be converted to integer, check count file.\n')
                        raise SystemExit('Count value could not be converted to integer, check count file.')

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
    except IOError:
        logger.error('File {} could not be opened.\n'.format(count_file))
        raise SystemExit('File {} could not be opened.'.format(count_file))
    return rpk_sum, rpk_dict


# This function takes rpk files and generates a tpm file based on that file and a provided scaling value.

def get_tpm(rpk_dict, scaling, count_type, direction, SampleName):
    try:
        with open(str(SampleName) + '.' + direction + "." + count_type + '.tpm.txt', 'w') as ftpm:
            for gene_id, rpk_value in rpk_dict.items():
                try:
                    tpm_value = rpk_value / scaling
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
    return


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

        # Get RPK info - function returns two values [total_rpk, dictionary{gene_id: rpk_values} ]
        UniqueFWD_rpk = get_rpk(args.UniqueFWD, db)
        UniqueREV_rpk = get_rpk(args.UniqueREV, db)
        NonuniqueFWD_rpk = get_rpk(args.NonuniqueFWD, db)
        NonuniqueREV_rpk = get_rpk(args.NonuniqueREV, db)

        # Calculate scaling value
        total_unique_rpk_sum = UniqueFWD_rpk[0] + UniqueREV_rpk[0]
        total_nonunique_rpk_sum = NonuniqueFWD_rpk[0] + NonuniqueREV_rpk[0]

        unique_scaling_value = total_unique_rpk_sum / 1000000
        nonunique_scaling_value = (total_nonunique_rpk_sum / 1000000) - unique_scaling_value

        # Calculate tpm value
        get_tpm(UniqueFWD_rpk[1], unique_scaling_value, 'unique', 'fwd', args.SampleName)
        get_tpm(UniqueREV_rpk[1], unique_scaling_value, 'unique', 'rev', args.SampleName)
        get_tpm(NonuniqueFWD_rpk[1], nonunique_scaling_value, 'nonunique', 'fwd', args.SampleName)
        get_tpm(NonuniqueREV_rpk[1], nonunique_scaling_value, 'nonunique', 'rev', args.SampleName)

    else:
        # Get TPM values - function returns two values [total_rpk, dictionary{gene_id: rpk_values} ]
        UniqueFWD_rpk = get_rpk(args.UniqueFWD, db)
        NonuniqueFWD_rpk = get_rpk(args.NonuniqueFWD, db)

        # Calculate the scaling value
        total_unique_rpk_sum = UniqueFWD_rpk[0]
        total_nonunique_rpk_sum = NonuniqueFWD_rpk[0]

        unique_scaling_value = total_unique_rpk_sum / 1000000
        nonunique_scaling_value = (total_nonunique_rpk_sum / 1000000) - unique_scaling_value


        # Calculate the TPM Values
        get_tpm(UniqueFWD_rpk[1], unique_scaling_value, 'unique', 'fwd', args.SampleName)
        get_tpm(NonuniqueFWD_rpk[1], nonunique_scaling_value, 'nonunique', 'fwd', args.SampleName)


if __name__ == "__main__": __main__()


