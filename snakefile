samples = glob_wildcards("data/mydatalocal/ataqseq/subset/subset/{sample}.fastq.gz")
configfile: "config/config.yaml",
rule all:
	input:
		expand("results/fastqc_init/{sample}_fastqc.zip", sample = config["samples"]),
		expand("results/fastqc_init/{sample}_fastqc.html", sample = config["samples"]),
		expand("results/trim/{sample}_1_paired_trimmed.fastq.gz",sample = config["samples"]),
		expand("results/trim/{sample}_2_paired_trimmed.fastq.gz",sample =config["samples"]),
		expand("results/post_fastqc/{sample}_fastqc.zip", sample = config["samples"]),
		expand("results/post_fastqc/{sample}_fastqc.html", sample = config["samples"]),
		expand("results/alignment/{sample}.bam", sample= config["samples"]),
		expand("results/igvtools/{sample}.tdf", sample= config["samples"])
rule unzip:
	input:
		lambda wildcards: config["samples"][wildcards.sample]
	output:
		temp("tmp/{sample}.fastq")
	shell:
		"""
		mkdir -p tmp
		gunzip -c {input} > {output}
		"""
rule fastqc_init:
	input:
		"tmp/{sample}.fastq"
	output:
		"results/fastqc_init/{sample}_fastqc.html",
		"results/fastqc_init/{sample}_fastqc.zip"
	conda:
		"env.yaml"
	threads: 2
	shell:
		"""
		mkdir -p results/fastqc_init
		fastqc {input} -o "results/fastqc_init" -t {threads}
		"""
rule trimming:
	input:
		r1="data/mydatalocal/ataqseq/subset/subset/{sample}_1.fastq.gz",
		r2="data/mydatalocal/ataqseq/subset/subset/{sample}_2.fastq.gz"
	output:
		fwd_paired = "results/trim/{sample}_1_paired_trimmed.fastq.gz",
		fwd_unpaired = "results/trim/{sample}_1_unpaired_trimmed.fastq.gz",
		rev_paired = "results/trim/{sample}_2_paired_trimmed.fastq.gz",
		rev_unpaired = "results/trim/{sample}_2_unpaired_trimmed.fastq.gz"
	params:
		TRIMMOMATIC_JAR = "/path/to/download-software/trimmomatic-0.38-1/trimmomatic.jar"
	log: 
		"logs/trimmomatic/{sample}.PE.trimomatic.log"
	shell: 
		"""
		mkdir -p results/trim
		mkdir -p logs/trimmomatic
		{input.r1} {input.r2} \
		{output.fwd_paired} {output.fwd_unpaired} \
		{output.rev_paired} {output.rev_unpaired} ILLUMINACLIP:data/mydatalocal/ataqseq/raw/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 HEADCROP:18

		"""
rule post_fastqc:
	input:  
		fwd_paired = "results/trim/{sample}_1_paired_trimmed.fastq.gz",
		fwd_unpaired = "results/trim/{sample}_1_unpaired_trimmed.fastq.gz",
		rev_paired = "results/trim/{sample}_2_paired_trimmed.fastq.gz",
		rev_unpaired = "results/trim/{sample}_2_unpaired_trimmed.fastq.gz"
	output: 
		html="results/post_fastqc/{sample}_fastqc.html",
		zip="results/post_fastqc/{sample}_fastqc.zip"
	shell: 
		"""
		mkdir -p results/post_fastqc
		fastqc {input} -o "results/post_fastqc" -t 8
		"""
rule mapping:
	input: 
		fwd = "results/trim/{sample}_1_paired_trimmed.fastq.gz",
		rev = "results/trim/{sample}_2_paired_trimmed.fastq.gz"
	output: 
		"results/alignment/{sample}.bam"
	params: 
		threads=8
	
	log:
		"logs/alignment/{sample}.log"
	shell: 
		"""
		mkdir -p results/alignment
		mkdir -p logs/bwa
		"bwa mem {input.fwd} {input.rev} {params.threads} | \
		samtools view -b  - > {output} 2>{log}
		"""
rule bam_to_tdf:
	input: 
		"results/alignment/{sample}.bam"
	output: 
		"results/igvtools/{sample}.tdf"
	params:
		igv_genome = "mm10/GRC38"
	log:
		"logs/bam_to_tdf.{sample}.log"
	conda:
		"envs/igvtools.yaml"
	shell:
		"""
		igvtools count \
		--minMapQuality 20 \
		-z 5 \
		-w 25 \
		-e 225 \
		{input} \
		{output} \
		{params.igv_genome}
		"""
rule deeptools_summary:
	input:
		expand("results/alignment/{sample}.bam", sample= config["samples"])
	output:
		sum     = "results/deeptools/multibamsum.npz",
		counts  = "results/deeptools/multibamsum.tab"
	params:
		labels=config["samples"].tolist()
	threads: 8
	log:
		"logs/deeptools_summary.log"
	conda:
		"envs/deeptools.yaml"
	shell:
		"""
		multiBamSummary bins \
		-p {threads} \
		-b {input} \
		--minMappingQuality 20 \
		-e 225 \
		--ignoreDuplicates \
		--centerReads \
		--labels {params.labels} \
		-out {output.sum} \
		--outRawCounts {output.counts}
		"""
rule deeptools_correlation:
	input: 
		"results/deeptools/multibamsum.npz"
	output:
		fig     = "results/deeptools/pearsoncor_multibamsum.png",
		matrix  = "results/deeptools/pearsoncor_multibamsum_matrix.txt"
	conda:
		"envs/deeptools.yaml"
	log:
		"logs/deeptools_correlation.log"
	shell:
		"""
		plotCorrelation \
		--corData {input} \
		--plotFile {output.fig} \
		--outFileCorMatrix {output.matrix} \
		--corMethod pearson \
		--whatToPlot heatmap \
		--skipZeros \
		--plotTitle "Pearson Correlations of BWA Alignments" \
		--plotNumbers \
		--colorMap RdYlBu
		"""
rule deeptools_coverage:
	input: 
		expand("results/alignment/{sample}.bam", sample= config["samples"])
	output:
		fig     = "results/deeptools/multibamsum_cov.png",
		counts  = "results/deeptools/multibamsum_cov_counts.txt"
	threads: 8
	params:
		labels=config["samples"].tolist()
	conda:
		"envs/deeptools.yaml"
	log:
		"logs/deeptools_coverage.log"
	shell:
		"""
		plotCoverage \
		-p {threads} \
		-b {input} \
		--plotFile {output.fig} \
		--output.Counts {output.counts} \
		-e 225 \
		--plotTitle "Coverage of BWA Alignments" \
		--labels {params.labels} \
		--minMappingQuality 20 \
		--ignoreDuplicates \
		--skipZeros \
		--centerReads
		"""
rule macs2_callpeaks:
	input:
		"results/alignment/{id}Input.bam"
	output:
		"results/macs2/{id}_peaks.narrowPeak",
		"results/macs2/{id}_peaks.xls",
		"results/macs2/{id}_summits.bed"
	params:
		id          = "{id}",
		organism    = config["macs2"]["gsize"]
	conda:
		"envs/macs2.yaml"
	log:
		"logs/macs2_callpeaks.{id}.log"
	shell:
		"""
		macs2 callpeak \
		--name {params.id} \
		--treatment {input} \
		--control {input.input} \
		--gsize {params.organism} \
		--outdir results/macs2
		"""
