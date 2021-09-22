nextflow.enable.dsl = 2

process VIRSORTER2 {
	cpus 30
    publishDir "${params.outdir}/virsorter2", mode: params.publish_dir_mode

    container "quay.io/biocontainers/virsorter=2.2.3--pyhdfd78af_1"
    conda (params.enable_conda ? "bioconda::virsorter=2.2.3" : null)
	
	input:
	path fasta
	path db

	output:
	path "for-dramv/final-viral-combined-for-dramv.fa", emit: fasta
	path "for-dramv/viral-affi-contigs-for-dramv.tab", emit: tab

	script:
	"""
	virsorter run \
		--min-length 0 --min-score 0 -j $task.cpus \
		--seqname-suffix-off --viral-gene-enrich-off \
		--prep-for-dramv \
		-i $fasta -w ./ --db-dir $db \
		all
	"""
}

process DRAMV {
	cpus 30
    publishDir "${params.outdir}/DRAMv", mode: params.publish_dir_mode

    container "quay.io/biocontainers/dram=1.2.4--pyhdfd78af_1"
    conda (params.enable_conda ? "bioconda::dram=1.2.4" : null)
	
	input:
	path fasta
	path tab
	path db

	output:
	path "proteins.txt", emit: txt
	path "hits.tsv", emit: tsv

	script:
	"""
	DRAM-v.py annotate \
		-i $fasta -v $tab \
		-o . \
		--skip_trnascan --threads $task.cpus --min_contig_size 0

	cut -f3,11 annotations.tsv > hits.tsv
	cut -f11 annotations.tsv | grep -v "^\$\$" | sort | uniq  > proteins.txt
	"""
}

process ACC_TO_TAXID{
	publishDir "${params.outdir}/esearch", mode: params.publish_dir_mode

	container "quay.io/biocontainers/entrez-direct=15.6--he881be0_1"
    conda (params.enable_conda ? "bioconda::dram=1.2.4" : null)	
	
    input:
    val ids

    output:
	path "taxids.tsv"

    script:
    """
    #!/usr/bin/env bash

    esearch -db protein -query "${ids.replaceAll('\n', ' ')}" | 
		esummary |
        xtract -pattern DocumentSummary -element Caption,TaxId \
        > taxids.tsv
    """
}

process TAXONKIT_LINEAGE {
    publishDir "${params.outdir}/taxonkit", mode: params.publish_dir_mode

    container "quay.io/biocontainers/taxonkit=0.8.0--h9ee0642_0"
    conda (params.enable_conda ? "bioconda::taxonkit=0.8.0" : null)
	
	input:
	path taxids

	output:
	path "lineages.tsv"

	script:
	"""
	cat <(echo "acc taxid lineage ranks" | tr " " "\\t") \
		<(taxonkit lineage -R -i 2 $taxids) \
		> lineages.tsv
	"""
}

process VOTE {
    publishDir "${params.outdir}/consensus", mode: params.publish_dir_mode

	container "nakor:metaflowmics-python=0.0.1"
	conda (params.enable_conda ? "conda-forge::pandas" : null)
	
	input:
	path mapping
	path lineages

	output:
	path "taxonomy.csv"

	script:
	"""
	python bin/consensus-taxonomy.py 
	"""	
}

workflow lineage {

	take:
	contigs
	vs2_db
	dram_db
	

	main:
	VIRSORTER2(
		contigs,
		vs2_db
	)
	DRAMV(
		VIRSORTER2.out.fasta,
		VIRSORTER2.out.tab,
		dram_db
	)
	
	ids_split = Channel.fromPath(
		DRAMV.out.proteins
	).splitText( by: params.chunks )
	
	ACC_TO_TAXID(
		ids_split
	)

	TAXONKIT_LINEAGE(
		ACC_TO_TAXID.out.collectFile(name: "taxids.tsv")
	)
	
	VOTE(
		DRAMV.out.tsv,
		TAXONKIT_LINEAGE.out
	)
}

workflow {
	lineage(
		file("$params.contigs", checkIfExists: true),
		Channel.fromPath("$params.vs2_db/*").collect(),
		file("$params.dramv_db", checkIfExists: true),
	)
}
