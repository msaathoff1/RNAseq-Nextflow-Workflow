# %% DESeq2 Script
import logging
import argparse
import os
import collections
import subprocess

# Set up logger for error handling and import subprocess

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler('logFile.log')
fh.setLevel(logging.INFO)
fh.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
logger.addHandler(fh)


# This function takes a directory and returns a list of the files inside. - Use if you have a lot of files and
# don't want to list them all individually as an input.

def list_of_files(countfiles_dir):
    path = countfiles_dir
    return os.listdir(path)


# The function takes a directory called countFiles and a .txt file called configFile. The config file must be a two
# column tab separated file with header "sampleName \t condition". It returns [a_list_of_the_conditions, sampleTable].
# Note the list of conditions DOES contain duplicates. The sample table is a list of lists, with each entry being:
# [Sample Name, File Name, Condition]. There is one entry for each sample/file pair.

def prepFiles(countfiles_list, configFile):
    with open(configFile, 'r') as config:
        # strip the header
        firstline = config.readline().rstrip()
        colTitle = firstline.split('\t')
        conditions = []
        sampleTable = []

        # check the config file formatting
        if colTitle[0] == "sampleName" and colTitle[1] == "condition":
            for line in config:
                try:
                    rowValues = line.strip().split('\t')
                    configSampleName = rowValues[0]
                    condition = rowValues[1]
                    conditions.append(condition)  # build a list of all conditions in the config file
                except IndexError:
                    logger.error("Index out of range. Check config file for formatting problems. Mind the tabs.")
                    raise SystemExit("Index out of range. Check config file for formatting problems. Mind the tabs.")

                try:
                    for file in countfiles_list:
                        filename = str(os.fsdecode(file))  # get the names of the files.
                        if configSampleName in filename:  # match the sample names to the files.
                            logger.info([str(configSampleName), str(filename), str(condition)])
                            sampleTable.append([str(configSampleName), str(filename), str(condition)])
                            # NOTE - SAMPLE/FILE NAMES MUST BE DISTINCT - "in" will conflate "treated, untreated"
                        else:
                            continue
                except FileNotFoundError:
                    logger.error("Cannot find count data files. Check the given directory and path.")
                    raise SystemExit("Cannot find count data files. Check the given directory and path.")
        else:
            logger.error("Please format the tab-separated config file with a first column named 'sampleName' and a "
                         "second column named 'condition'")
            raise SystemExit("Please format the tab-separated config file with a first column named 'sampleName' and "
                             "a second column named 'condition'")
    if len(sampleTable) == 0:
        logger.error("Cannot find sample information! - SampleTable is empty!")
        raise SystemExit("Cannot find sample information! - SampleTable is empty!")
    elif len(conditions) == 0:
        logger.error("Cannot find conditions!")
        raise SystemExit("Cannot find conditions!")
    else:
        return conditions, sampleTable


def perform_deseq2(conditions, sampleTable, countfiles_list):
    # check for min. of 2 replicates - ensure each condition is listed at least twice in the conditions list
    # (this one does contain duplicates; should be listed once for each replicate ergo min of 2.)
    count = collections.Counter()
    for cond in conditions:
        count[cond] += 1
    if any([True for k, v in count.items() if v < 2]):
        logger.error("You must have at least 2 replicates per condition.")
        raise SystemExit("You must have at least 2 replicates per condition.")
    else:
        # fsencode(countFile_dir) ??
        conditions_set = set(conditions)
        conditions_list = list(conditions_set)  # get rid of the duplicates

        # check there are at least 2 conditions and 4 samples provided
        if len(conditions_list) < 2:
            logger.error("Two or more conditions are required.")
            raise SystemExit("Two or more conditions are required")
        elif len(sampleTable) < 4:
            logger.error("You must provide 2 or more conditions AND 2 or more replicates. Some samples may be missing.")
            raise SystemExit("You must provide 2 or more conditions AND 2 or more replicates. Some samples may be "
                             "missing.")
        else:
            os.mkdir('countfiles_dir')
            os.mkdir('deseq2results')
            for file in countfiles_list:
                origin = file
                destination = 'countfiles_dir/' + str(os.fsdecode(file))
                os.rename(origin, destination)
            for i in range(len(conditions_list)):  # iterate through each possible combo of conditions
                for j in range(i + 1, len(conditions_list)):
                    condition_a = conditions_list[i]
                    condition_b = conditions_list[j]
                    with open("SampleTable_R.txt", 'w') as sampleTable_R:  # build the sample table for R.
                        sampleTable_R.write("sampleName\tfileName\tcondition\n")
                        for entry in sampleTable:
                            if entry[2] == condition_a or entry[2] == condition_b:
                                sampleTable_R.write("\t".join(str(i) for i in entry))
                                sampleTable_R.write("\n")
                    sampleTable_R.close()  # explicitly close so there is a sep. file for each combo
                    logger.info(["DEseq2Rbit.r", "SampleTable_R.txt", "countfiles_dir", str(condition_a), str(condition_b)])
                    subprocess.run(["DEseq2Rbit.r", "SampleTable_R.txt", "countfiles_dir", str(condition_a), str(condition_b)])
    return "done"


def __main__():
    # Set up command line input
    parser = argparse.ArgumentParser(description='Calculate TPM values')
    parser.add_argument('--countFiles', required=True, help='Directory containing output files from HTSeq-count.')
    parser.add_argument('--config', required=True, help='File matching sample to condition status.')
    parser.add_argument('--countFilesType', required=True, help="Takes either 'dir' or 'list' as a setting.")

    args = parser.parse_args()
    if args.countFilesType == 'dir':
        # Prep the files and run DESeq2
        files_list = list_of_files(args.countFiles)
        os.chdir(args.countFiles)
        get_info = prepFiles(files_list, args.config)
        all_results = perform_deseq2(get_info[0], get_info[1], files_list)
        print(all_results)
    elif args.countFilesType == 'list':
        countFiles = args.countFiles.split()
        get_info = prepFiles(countFiles, args.config)
        all_results = perform_deseq2(get_info[0], get_info[1], countFiles)
        print(all_results)


if __name__ == "__main__": __main__()

# create a mode where prep is called within deseq2? - might be nice. Then just 1 f call.
