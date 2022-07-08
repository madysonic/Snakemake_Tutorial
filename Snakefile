# Tutorial by Nathan Watson-Haigh
# SAGC workshop 6th July 2022
# A snakemake workflow for indexing, fastqc, trimming and alignment 


# Rather than install all the required tools, we could us a singularity container
# A global singularity image. Need to specify --use-singularity and have singularity available on the command line

#singularity:
#	"file:///shared/.singularity/workflows-training.simg"
#	"docker://continuumio/miniconda3:4.5.12"


# List of sample names - used as a variable in the workflow
SAMPLES = [
	"Lancer",
	"Mace"
]



################
# Pseudo-rules #
################

# By convention, the first rule should be called "all" and it's "input" defined as
# the list of ALL the files you want the workflow to create. e.g.:

# Use wildcards/variables for easier implementation 

rule all:
	input:
		expand("references/reference.fasta.gz.{ext}", ext=['amb', 'ann', 'bwt', 'pac', 'sa']),
		expand("raw_reads/{SAMPLE}_{read}_fastqc.html", SAMPLE=SAMPLES, read=['R1', 'R2']),
		expand("qc_reads/{SAMPLE}_{read}.fastq.gz", SAMPLE=SAMPLES, read=['R1', 'R2']),
		expand("mapped/{SAMPLE}.bam", SAMPLE=SAMPLES)


################
# Rules Proper #
################

# Each task to run
# Remember to specify the conda environment yaml in each rule

rule bwa_index:
	input:
		"{ref}",
	output:
		"{ref}.amb",
		"{ref}.ann",
		"{ref}.bwt",
		"{ref}.pac",
		"{ref}.sa",
	conda:
		"default.yaml",
	shell:
		"""
		bwa index \
		  -a bwtsw \
		  {input}
		"""
rule fastqc:
	input:
		"raw_reads/{prefix}.fastq.gz",
	output:
		zip = "raw_reads/{prefix}_fastqc.zip",
		html = "raw_reads/{prefix}_fastqc.html",
	conda:
		"default.yaml",
	shell:
		"""
		fastqc --threads {threads} {input}
		"""
rule trimmomatic:
	input:
		r1 = "raw_reads/{SAMPLE}_R1.fastq.gz",
		r2 = "raw_reads/{SAMPLE}_R2.fastq.gz",
		adapters = "misc/TruSeq3-PE.fa"
	output:
		r1 = "qc_reads/{SAMPLE}_R1.fastq.gz",
		r2 = "qc_reads/{SAMPLE}_R2.fastq.gz",
		r1_unpaired = "qc_reads/{SAMPLE}_R1.unpaired.fastq.gz",
		r2_unpaired = "qc_reads/{SAMPLE}_R2.unpaired.fastq.gz",
	conda:
		"default.yaml",
	shell:
		"""
		trimmomatic PE \
		-threads {threads} \
		{input.r1} {input.r2} \
		{output.r1} {output.r1_unpaired} \
		{output.r2} {output.r2_unpaired} \
		ILLUMINACLIP:{input.adapters}:2:30:10:3:true \
		LEADING:2 \
		TRAILING:2 \
		SLIDINGWINDOW:4:15 \
		MINLEN:36
		"""

rule bwa_mem:
	input:
		index = expand("references/reference.fasta.gz.{ext}", ext=["amb","ann","bwt","pac","sa"]),
		r1    = "qc_reads/{SAMPLE}_R1.fastq.gz",
		r2    = "qc_reads/{SAMPLE}_R2.fastq.gz",
	output:
		"mapped/{SAMPLE}.bam"
	params:
		#prefix = "references/reference.fasta.gz",
		prefix = lambda wildcards, input: input["index"][0][:-4],
	conda:
		"default.yaml",
	shell:
		"""
		bwa mem -t 1 \
		  {params.prefix} \
		  {input.r1} {input.r2} \
		| samtools view -b \
		> {output}
		"""
