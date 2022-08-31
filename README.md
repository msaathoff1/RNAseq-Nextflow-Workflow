# RNAseq-Nextflow-Workflow

Process RNA sequencing data sets in preparation for upload, storage, and
viewing as a part of VEuPathDB. Generate required bigwig files, TPM-
normalized count files, pre-calculated differential expression values for
each gene, and quality information.

Contents:
install.me.sh
nfscript.nf
nextflow.config : this is a template that must be filled.
custom_processes

custom_processes is a directory that contains:

deseq2_dockerinfo
getphredencoding_tpm_dockerinfo
hisat2samtools_dockerinfo
{UTR_getphredencoding_tpm_dockerinfo } This is specific to this project.

Installation:
Docker must be installed independently by the user before install.me.sh
can be run. Instructions can be found at:
https://docs.docker.com/desktop/install/ubuntu/ . 

Once this is complete,
install.me.sh can be run as a bash script. It will install Nextflow
itself as well as build/pull all docker images.

Input:
1. The inputs for this workflow are .fastq files stored in a directory
whose path is specified in the nextflow.config file.
2. The complete and accurate nextflow.config file must be in the
directory from which Nextflow is being run.

Usage:
nextflow [options] nfscript.nf
Options:
-C
  Use the specified configuration file(s) overriding any defaults
-D
  Set JVM properties
-bg
  Execute nextflow in background
-c, -config
  Add the specified file to configuration set
-d, -dockerize
  Launch nextflow via Docker (experimental)
-h
  Print this help
-log
  Set nextflow log file path
-q, -quiet
  Do not print information messages
-syslog
  Send logs to syslog server (eg. localhost:514)
-v, -version
  Print the program version

The nextflow.config file must be in the directory from where Nextflow is
being run. The parameters must all be set to their appropriate values.
The required configuration values are:

params.reads = [path/*{1,2}.fastq.gz]
params.reference = [path/reference_genome.fasta]
params.gff = [path/reference_genome.gff]
params.config = [path/sample_table.txt]
params.outdir = [path/outdir]
params.paired = [True/False]
params.stranded = [True/False]

docker {
enabled = true
}

All input read files must be named with the convention
Name_{1,2}.fastq.gz (or .fastq). Even if the input files are single-end,
a “_1” suffix must be added to the name in order for workflow naming
conventions to be satisfied. This preserves file validation within the
workflow. The reference genome must be in .fasta format. The
reference .gff file must also be provided. The config parameter is a tab
separated file containing the name of each sample and its corresponding
experimental condition. It must have the header “sampleName\tcondition”
where ‘\t’ denotes a tab character. See the DESeq2 section of methods for
an example. If differential expression analysis should not be conducted,
this file can be left empty, but it must be created. The output directory
(outdir) will default to being created in the current directory and is
named by the user. The paired and stranded parameters should be set to
true/false as appropriate. Docker is enabled so processes are run within
their respective containers.

Output:
This Nextflow workflow generates normalized read count files, bigwig
files, directories containing quality results, and a directory containing
differential expression values for each combination of experimental
conditions [deseq2results]. Note DESeq2 is performed for unique and sense
(or unstranded) read counts only.


