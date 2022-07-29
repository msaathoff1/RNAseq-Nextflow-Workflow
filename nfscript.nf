#!/usr/bin/env nextflow
/*
* Nextflow workflow for processing RNA-seq Data
*
* Author: Megan Saathoff
*/

/*
* pipeline input parameters
*/

params.reads = "/home/avianeagle/bioinfo/TestSamples/Chr1Subset/*_{1,2}.fastq.gz"
params.reference = "/home/avianeagle/bioinfo/TestSamples/Chr1Subset/PlasmoDB-57_PbergheiANKA_chr1.fasta"
params.outdir = "/home/avianeagle/bioinfo/testing_results"
params.gff = "/home/avianeagle/bioinfo/TestSamples/Chr1Subset/PlasmoDB-57_PbergheiANKA_chr1.gff"
params.config = "/home/avianeagle/bioinfo/TestSamples/sampletable-nextflow.txt"
params.paired = false
params.stranded = true

log.info """\
	R N A S E Q - NF P I P E L I N E
	================================
	reference: ${params.reference}
	reads: ${params.reads}
	outdir: ${params.outdir}
	"""
	.stripIndent()

/*
* Index provided reference using reference file
*/
process index {

	container = 'biocontainers/hisat2:2.2.1--h87f3376_4'
	
	input:
	path reference from params.reference
	
	output:
	path 'referenceIndex*.ht2' into index_ch
	val 'referenceIndex' into index_value_ch
	
	script:
	"""
	hisat2-build ${reference} referenceIndex
	"""
	
}
/*
* Form read pairs if data is paired - form sample list and names if not
*/
if (params.paired) {

	Channel
		.fromFilePairs(params.reads, checkIfExists:true)
		.set {read_pairs_ch}
/*
* Pass read pairs into a channel for quality assessment
*/
	Channel
		.fromFilePairs(params.reads, checkIfExists:true)
		.set {qualityEncoding_ch}
}
else {
	Channel
		.fromPath(params.reads, checkIfExists:true)
		.map { file -> tuple(file.simpleName, [file]) }
		.set {read_pairs_ch}
	
	Channel
		.fromPath(params.reads, checkIfExists:true)
		.map { file -> tuple(file.simpleName, [file]) }
		.set {qualityEncoding_ch}
}
/*
* Run FASTQC for quality purposes
*/
	process quality {
	
		container = 'biocontainers/fastqc:v0.11.9_cv7'
		
		publishDir params.outdir, mode:'copy'
		
		input:
		tuple val(pair_id), path(reads) from qualityEncoding_ch
		
		output:
		path('fastqc_dir', type:'dir') into fastqc_ch
		tuple val(pair_id), path('fastqc_dir', type:'dir') into get_phred_ch
		
		script:
		"""
		mkdir fastqc_dir
		fastqc -o fastqc_dir --extract ${reads}
		"""
	
	}
/*
* Get phred score
*/
	process get_phred {
		
		container = 'python:3.11.0b4-alpine3.16'
		
		input:
		tuple val(pair_id), path('fastqc_dir') from get_phred_ch
		
		output:
		path("phredEncoding.txt") into (trim_phred_ch, align_phred_ch)
		
		shell:
		if (params.paired)
			"""
			python3 ~/bin/RunGetPhredEncoding2.py ${fastqc_dir}/${pair_id}_1_fastqc > phredEncoding.txt
			"""
		else
			"""
			python3 ~/bin/RunGetPhredEncoding2.py ${fastqc_dir}/${pair_id}_fastqc > phredEncoding.txt
			"""
	}
/*
* Trim reads for quality purposes
*/

	process trim{
	
		container = 'biocontainers/trim-galore:0.6.7--hdfd78af_0'
		
		input:
		tuple val(pair_id), path(reads) from read_pairs_ch
		val(phredEncoding) from trim_phred_ch.map { it.text.strip() }
		
		output:
		tuple val(pair_id), path("${pair_id}_1_val_1.fq.gz"), path("${pair_id}_2_val_2.fq.gz") into trimmed_pairs_ch
		
		script:
		if(params.paired)
			"""
			trim_galore --paired --${phredEncoding} ${reads[0]} ${reads[1]}
			"""
		else
			"""
			trim_galore --${phredEncoding} ${reads[0]}
			mv ${pair_id}_trimmed.fq.gz ${pair_id}_1_val_1.fq.gz
			touch ${pair_id}_2_val_2.fq.gz
			"""
 }

/*
* Align reads to reference
*/

	process align {
	
		container = 'hisat2samtools'
		
		publishDir params.outdir, mode:'copy'
		
		input:
		tuple val(pair_id), path("${pair_id}_1_val_1.fq.gz"), path("${pair_id}_2_val_2.fq.gz") from trimmed_pairs_ch
		path 'referenceIndex.*.ht2' from index_ch
		val 'referenceIndex' from index_value_ch
		val(phredEncoding) from align_phred_ch.map { it.text.strip() }
		
		output:
		tuple val(pair_id), path("${pair_id}.bam") into (count_bam_ch, qsplit_bam_ch)
		tuple val(pair_id), path("${pair_id}.bam.bai") into bam_index_ch
		
		script:
		if (params.paired)
		"""
		hisat2 --${phredEncoding} -x ${referenceIndex} -1 ${pair_id}_1_val_1.fq.gz -2 ${pair_id}_2_val_2.fq.gz | samtools view -bS | samtools sort > ${pair_id}.bam
		samtools index ${pair_id}.bam
		"""
		else
		"""
		hisat2 --${phredEncoding} -x ${referenceIndex} -U ${pair_id}_1_val_1.fq.gz | samtools view -bS | samtools sort > ${pair_id}.bam
		samtools index ${pair_id}.bam
		"""
 }
 /*
* Split the bam file for generating bigwig files
*/
	process split_bam {
		
		container = 'biocontainers/samtools:v1.9-4-deb_cv1'
		
		publishDir params.outdir, mode:'copy'
		
		input:
		tuple val(pair_id), path("${pair_id}.bam") from qsplit_bam_ch
		
		output:
		tuple val(pair_id), path("${pair_id}.unique.fwd.bam"), path("${pair_id}.nonunique.fwd.bam"), path("${pair_id}.unique.rev.bam"), path("${pair_id}.nonunique.rev.bam") into split_bam_ch
		tuple val(pair_id), path("${pair_id}.unique.fwd.bam.bai"), path("${pair_id}.nonunique.fwd.bam.bai"), path("${pair_id}.unique.rev.bam.bai"), path("${pair_id}.nonunique.rev.bam.bai") into split_index_ch
		
		script:
		if (params.stranded)
			if(params.paired)
				"""
				## Split unique and non-unique
				samtools view -h  ${pair_id}.bam | grep -P '^@|NH:i:1(\\s|\$)' |samtools view -h -bS - > ${pair_id}.unique.bam
				samtools index ${pair_id}.unique.bam
				samtools view -h  ${pair_id}.bam  | grep -P '^@|NH:i:([2-9]|1\\d)' |samtools view -h -bS - > ${pair_id}.nonunique.bam
				samtools index ${pair_id}.nonunique.bam
				
				## Generate unique fwd
				samtools view -b -f 163 ${pair_id}.unique.bam > unique.fwd1.bam
				samtools index unique.fwd1.bam
				
				samtools view -b -f 83 ${pair_id}.unique.bam > unique.fwd2.bam
				samtools index unique.fwd2.bam
				
				samtools merge -f ${pair_id}.unique.fwd.bam unique.fwd1.bam unique.fwd2.bam
				samtools index ${pair_id}.unique.fwd.bam
				
				## Generate non-unique fwd
				samtools view -b -f 163 ${pair_id}.nonunique.bam > nonunique.fwd1.bam
				samtools index nonunique.fwd1.bam
				
				samtools view -b -f 83 ${pair_id}.nonunique.bam > nonunique.fwd2.bam
				samtools index nonunique.fwd2.bam
				
				samtools merge -f ${pair_id}.nonunique.fwd.bam nonunique.fwd1.bam nonunique.fwd2.bam
				samtools index ${pair_id}.nonunique.fwd.bam
				
				## Generate unique rev
				samtools view -b -f 147 ${pair_id}.unique.bam > unique.rev1.bam
				samtools index unique.rev1.bam
				
				samtools view -b -f 99 ${pair_id}.unique.bam > unique.rev2.bam
				samtools index unique.rev2.bam
				
				samtools merge -f ${pair_id}.unique.rev.bam unique.rev1.bam unique.rev2.bam
				samtools index ${pair_id}.unique.rev.bam
				
				## Generate non-unique rev
				samtools view -b -f 147 ${pair_id}.nonunique.bam > nonunique.rev1.bam
				samtools index nonunique.rev1.bam
				
				samtools view -b -f 99 ${pair_id}.nonunique.bam > nonunique.rev2.bam
				samtools index nonunique.rev2.bam
				
				samtools merge -f ${pair_id}.nonunique.rev.bam nonunique.rev1.bam nonunique.rev2.bam
				samtools index ${pair_id}.nonunique.rev.bam
				"""
			else
				"""
				## Split unique and non-unique
				samtools view -h  ${pair_id}.bam | grep -P '^@|NH:i:1(\\s|\$)' |samtools view -h -bS - > ${pair_id}.unique.bam
				samtools index ${pair_id}.unique.bam
				samtools view -h  ${pair_id}.bam  | grep -P '^@|NH:i:([2-9]|1\\d)' |samtools view -h -bS - > ${pair_id}.nonunique.bam
				samtools index ${pair_id}.nonunique.bam
				
				samtools view -b -F 20 ${pair_id}.unique.bam > ${pair_id}.unique.fwd.bam
				samtools view -b -F 20 ${pair_id}.nonunique.bam > ${pair_id}.nonunique.fwd.bam
				samtools view -b -f 16 ${pair_id}.unique.bam > ${pair_id}.unique.rev.bam
				samtools view -b -f 16 ${pair_id}.nonunique.bam > ${pair_id}.nonunique.rev.bam
				
				samtools index ${pair_id}.unique.fwd.bam
				samtools index ${pair_id}.nonunique.fwd.bam
				samtools index ${pair_id}.unique.rev.bam
				samtools index ${pair_id}.nonunique.rev.bam
				"""
		else
			"""
			## Split unique and non-unique
			samtools view -h  ${pair_id}.bam | grep -P '^@|NH:i:1(\\s|\$)' |samtools view -h -bS - > ${pair_id}.unique.fwd.bam
			samtools view -h  ${pair_id}.bam  | grep -P '^@|NH:i:([2-9]|1\\d)' |samtools view -h -bS - > ${pair_id}.nonunique.fwd.bam
			
			samtools index ${pair_id}.unique.fwd.bam
			samtools index ${pair_id}.nonunique.fwd.bam
			
			touch ${pair_id}.unique.rev.bam
			touch ${pair_id}.unique.rev.bam.bai
			touch ${pair_id}.nonunique.rev.bam
			touch ${pair_id}.nonunique.rev.bam.bai
			"""
	}
/*
* Generate BigWig file
*/
	process generate_bigwig {
	
		container = 'biocontainers/deeptools:3.5.1--py_0'
	
		publishDir params.outdir, mode:'copy'
		
		input:
		tuple val(pair_id), path("${pair_id}.unique.fwd.bam"), path("${pair_id}.nonunique.fwd.bam"), path("${pair_id}.unique.rev.bam"), path("${pair_id}.nonunique.rev.bam") from split_bam_ch
		tuple val(pair_id), path("${pair_id}.unique.fwd.bam.bai"), path("${pair_id}.nonunique.fwd.bam.bai"), path("${pair_id}.unique.rev.bam.bai"), path("${pair_id}.nonunique.rev.bam.bai") from split_index_ch
		
		output:
		tuple val(pair_id), path("${pair_id}.unique.fwd.bw"), path("${pair_id}.nonunique.fwd.bw"), path("${pair_id}.unique.rev.bw"), path("${pair_id}.nonunique.rev.bw") into bigwig_ch
		
		script:
		if(params.stranded)
			"""
			bamCoverage -b ${pair_id}.unique.fwd.bam -o ${pair_id}.unique.fwd.bw
			bamCoverage -b ${pair_id}.nonunique.fwd.bam -o ${pair_id}.nonunique.fwd.bw
			bamCoverage -b ${pair_id}.unique.rev.bam -o ${pair_id}.unique.rev.bw
			bamCoverage -b ${pair_id}.nonunique.rev.bam -o ${pair_id}.nonunique.rev.bw
			"""
		else
			"""
			bamCoverage -b ${pair_id}.unique.fwd.bam -o ${pair_id}.unique.fwd.bw
			bamCoverage -b ${pair_id}.nonunique.fwd.bam -o ${pair_id}.nonunique.fwd.bw
			touch ${pair_id}.unique.rev.bw
			touch ${pair_id}.nonunique.rev.bw
			"""
 } 
 /*
* Calculate the counts
*/
	process count{
	
		container = 'biocontainers/htseq:v0.11.2-1-deb-py3_cv1'
		
		publishDir params.outdir, mode:'copy'
		
		input:
		tuple val(pair_id), path("${pair_id}.bam") from count_bam_ch
		tuple val(pair_id), path("${pair_id}.bam.bai") from bam_index_ch
		path refGFF from params.gff
		
		output:
		tuple val(pair_id), path("${pair_id}.fwd.unique.txt"), path("${pair_id}.rev.unique.txt"), path("${pair_id}.fwd.nonunique.txt"), path("${pair_id}.rev.nonunique.txt") into count_ch
		path("${pair_id}.fwd.unique.txt") into DESeq2Q_channel
		
		script:
		if (params.stranded)
			"""
			htseq-count --stranded reverse -t exon -i gene_id -f bam -r name ${pair_id}.bam $refGFF > ${pair_id}.rev.unique.txt
			htseq-count --stranded yes -t exon -i gene_id -f bam -r name ${pair_id}.bam $refGFF > ${pair_id}.fwd.unique.txt
			htseq-count --stranded reverse -t exon -i gene_id --nonunique all -f bam -r name ${pair_id}.bam $refGFF > ${pair_id}.rev.nonunique.txt
			htseq-count --stranded yes -t exon -i gene_id --nonunique all -f bam -r name ${pair_id}.bam $refGFF > ${pair_id}.fwd.nonunique.txt
			"""
		else
			"""
			htseq-count --stranded no -t exon -i gene_id -f bam -r name ${pair_id}.bam $refGFF > ${pair_id}.fwd.unique.txt
			htseq-count --stranded no --nonunique all -t exon -i gene_id -f bam -r name ${pair_id}.bam $refGFF > ${pair_id}.fwd.nonunique.txt
			touch ${pair_id}.rev.unique.txt
            touch ${pair_id}.rev.nonunique.txt
			"""
	}

 /*
* Calculate TPM values
*/
	process TPM{
	
		container = 'python:3.11.0b4-alpine3.16'
		
		publishDir params.outdir, mode:'copy'
		
		input:
		tuple val(pair_id), path("${pair_id}.fwd.unique.txt"), path("${pair_id}.rev.unique.txt"), path("${pair_id}.fwd.nonunique.txt"), path("${pair_id}.rev.nonunique.txt") from count_ch
		path gffFile from params.gff
		
		output:
		tuple val(pair_id), path("${pair_id}.fwd.unique.tpm.txt"), path("${pair_id}.fwd.nonunique.tpm.txt"), path("${pair_id}.rev.unique.tpm.txt"), path("${pair_id}.rev.nonunique.tpm.txt") into tpm_ch
		script:
		if (params.stranded)
			"""
			python3 ~/bin/TPMScript.py --Stranded --GffFile $gffFile --UniqueFWD ${pair_id}.fwd.unique.txt --NonuniqueFWD ${pair_id}.fwd.nonunique.txt --UniqueREV ${pair_id}.rev.unique.txt --NonuniqueREV ${pair_id}.rev.nonunique.txt --SampleName $pair_id
			"""
		else
			"""
			python3 ~/bin/TPMScript.py --GffFile $gffFile --UniqueFWD ${pair_id}.fwd.unique.txt --NonuniqueFWD ${pair_id}.fwd.nonunique.txt --SampleName $pair_id
			touch ${pair_id}.rev.unique.tpm.txt
            touch ${pair_id}.rev.nonunique.tpm.txt
			"""
	}

/*
* Run DESeq2 for relevant sample pairs - Only does 2 conditions; requires at least 2 replicates; only does for fwd_unique data
*/
	process DESeq2{
	
		publishDir params.outdir, mode:'copy'
		
		input: 
		path file_paths from DESeq2Q_channel.collect()
		path('config') from params.config
		
		output:
		path('deseq2results', type:'dir')
		
		script:
		"""
		python3 ~/bin/DEseq2.py --countFiles "${file_paths}" --countFilesType 'list' --config ${config}
		"""
	
	
	}



