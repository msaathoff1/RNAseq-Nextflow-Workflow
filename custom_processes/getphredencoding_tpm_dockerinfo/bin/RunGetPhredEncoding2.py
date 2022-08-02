#!/usr/bin/env python3
import argparse
from zipfile import ZipFile

parser = argparse.ArgumentParser(description="This will use a pre-specified query to examine a database and return "
                                             "information")
parser.add_argument('fqName', help='The input.zip file minus the .zip suffix. ie name.zip becomes name')
args = parser.parse_args()

def get_phred(fqname):
    with open(str(fqname) + '/fastqc_data.txt', 'r') as fq:
        for i, line in enumerate(fq):
            if i == 5:
                rowValues = line.strip().split('\t')
                encoding = rowValues[1]
                encoding = encoding.lower()
                phred = 'phred33'
                if 'solexa' in encoding:
                    phred = 'phred64'
                elif 'illumina' in encoding:
                    s1 = encoding.find('illumina')
                    if s1 != -1:
                        s2 = encoding[(s1+9):(s1+12)]
                        for version in ['1.3', '1.4', '1.5', '1.6', '1.7']:
                            if version in s2:
                                phred = 'phred64'
    return phred
phredEncoding = get_phred(args.fqName)
print(phredEncoding)
