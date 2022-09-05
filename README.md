# RNAseq-Nextflow-Workflow

Process RNA sequencing data sets in preparation for upload, storage, and
viewing as a part of VEuPathDB. Generate required bigwig files, TPM-
normalized count files, pre-calculated differential expression values for
each gene, and quality information.

#### Contents: <br />

install.me.sh <br />
nfscript.nf <br />
nextflow.config : this is a template that must be filled. <br />
custom_processes <br />

custom_processes is a directory that contains: <br />

deseq2_dockerinfo <br />
getphredencoding_tpm_dockerinfo <br />
hisat2samtools_dockerinfo <br />
{UTR_getphredencoding_tpm_dockerinfo } This is specific to this project. <br />

#### Installation: <br />

Docker must be installed independently by the user before install.me.sh
can be run. Instructions can be found at:
https://docs.docker.com/desktop/install/ubuntu/ . 

Once this is complete,
install.me.sh can be run as a bash script. It will install Nextflow
itself as well as build/pull all docker images.

#### Input:

1. The inputs for this workflow are .fastq files stored in a directory
whose path is specified in the nextflow.config file.
2. The complete and accurate nextflow.config file must be in the
directory from which Nextflow is being run.

#### Usage: <br />

nextflow [options] nfscript.nf <br />

#### Options: <br />

-C <br />
  Use the specified configuration file(s) overriding any defaults <br />
-D <br />
  Set JVM properties <br />
-bg <br />
  Execute nextflow in background <br />
-c, -config <br />
  Add the specified file to configuration set <br />
-d, -dockerize <br />
  Launch nextflow via Docker (experimental) <br />
-h <br />
  Print this help <br />
-log <br />
  Set nextflow log file path <br />
-q, -quiet <br />
  Do not print information messages <br />
-syslog <br />
  Send logs to syslog server (eg. localhost:514) <br />
-v, -version <br />
  Print the program version <br />

The nextflow.config file must be in the directory from where Nextflow is
being run. The parameters must all be set to their appropriate values.
The required configuration values are:

params.reads = [path/*{1,2}.fastq.gz] <br />
params.reference = [path/reference_genome.fasta] <br />
params.gff = [path/reference_genome.gff] <br />
params.config = [path/sample_table.txt] <br />
params.outdir = [path/outdir] <br />
params.paired = [True/False] <br />
params.stranded = [True/False] <br />

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

#### Output: <br />

This Nextflow workflow generates normalized read count files, bigwig
files, directories containing quality results, and a directory containing
differential expression values for each combination of experimental
conditions [deseq2results]. Note DESeq2 is performed for unique and sense
(or unstranded) read counts only.


