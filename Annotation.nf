#!/usr/bin/env nextflow

/*
* AUTHOR: CERC Genomic Medicine,  Daniel Taliun, PhD <daniel.taliun@mcgill.ca>
* VERSION: 3.0
* YEAR: 2022
*/


process vcf_by_chrom {
	label "VEP"

	cache "lenient"
	executor "local"
	cpus 1

	input:
	tuple path(vcf), path(vcf_index)

	output:
	tuple stdout, path(vcf), path(vcf_index)

	"""
	tabix -l ${vcf} | tr "\n" ","
	"""
}


process annotate_chunks {
	label "VEP"

	cache "lenient"
        scratch true
        errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return "retry" }
	maxRetries 3
	cpus 1

	containerOptions "-B ${params.vep_cache}:/opt/vep/.vep"

	input:
	tuple val(chrom), path(vcf), path(vcf_index), val(start), val(stop)

	output:
	tuple(val(chrom), val("${vcf.getSimpleName()}"), path("${vcf.getSimpleName()}.${chrom}_${start}_${stop}.vep.vcf.gz"), emit: annotations)
	path "*.vep.log", emit: logs

        publishDir "results/logs", pattern: "*.vep.log", mode: "copy"

	script:
	if (params.assembly == "GRCh38")
		"""
		export PERL5LIB=/opt/vep/.vep/Plugins/:\$PERL5LIB
		loftee_args=human_ancestor_fa:/opt/vep/.vep/loftee_db_${params.assembly}/human_ancestor.fa.gz,gerp_bigwig:/opt/vep/.vep/loftee_db_${params.assembly}/gerp_conservation_scores.homo_sapiens.GRCh38.bw,conservation_file:/opt/vep/.vep/loftee_db_${params.assembly}/loftee.sql${params.loftee_flags}
		bcftools view ${params.drop_genotypes} ${vcf} ${chrom}:${start}-${stop} | vep --cache --offline --assembly ${params.assembly} --format vcf --vcf --compress_output bgzip --force_overwrite --dir_cache /opt/vep/.vep/ --plugin LoF,loftee_path:/opt/vep/.vep/loftee_${params.assembly},\${loftee_args} --dir_plugins /opt/vep/.vep/loftee_${params.assembly} --plugin CADD,/opt/vep/.vep/CADD_${params.assembly}/whole_genome_SNVs.tsv.gz,/opt/vep/.vep/CADD_${params.assembly}/InDels.tsv.gz --plugin CONTEXT ${params.vep_flags} --warning_file STDERR --output_file STDOUT > ${vcf.getSimpleName()}.${chrom}_${start}_${stop}.vep.vcf.gz 2> ${vcf.getSimpleName()}.${chrom}_${start}_${stop}.vep.log
		sanity_check.py -a ${vcf.getSimpleName()}.${chrom}_${start}_${stop}.vep.vcf.gz
	 	"""
	else if (params.assembly == "GRCh37")
		"""
		export PERL5LIB=/opt/vep/.vep/Plugins/:\$PERL5LIB
		loftee_args=human_ancestor_fa:/opt/vep/.vep/loftee_db_${params.assembly}/human_ancestor.fa.gz,conservation_file:/opt/vep/.vep/loftee_db_${params.assembly}/phylocsf_gerp.sql${params.loftee_flags}
		bcftools view ${params.drop_genotypes} ${vcf} ${chrom}:${start}-${stop} | vep --cache --offline --assembly ${params.assembly} --format vcf --vcf --compress_output bgzip --force_overwrite --dir_cache /opt/vep/.vep/ --plugin LoF,loftee_path:/opt/vep/.vep/loftee_${params.assembly},\${loftee_args} --dir_plugins /opt/vep/.vep/loftee_${params.assembly} --plugin CADD,/opt/vep/.vep/CADD_${params.assembly}/whole_genome_SNVs.tsv.gz,/opt/vep/.vep/CADD_${params.assembly}/InDels.tsv.gz --plugin CONTEXT ${params.vep_flags} --warning_file STDERR --output_file STDOUT > ${vcf.getSimpleName()}.${chrom}_${start}_${stop}.vep.vcf.gz 2> ${vcf.getSimpleName()}.${chrom}_${start}_${stop}.vep.log
		sanity_check.py -a ${vcf.getSimpleName()}.${chrom}_${start}_${stop}.vep.vcf.gz
		"""
	else
		error "Invalid assembly name: ${params.assembly}"
}


process concatenate_chunks {
	label "VEP"

	cache "lenient"
	errorStrategy "retry"
	maxRetries 3
	cpus 1

	input:
	tuple val(chrom), val(name), path(annotated_vcfs)

        output:
        tuple val(chrom), val(name), path("${name}.${chrom}.vep.vcf.gz"), path("${name}.${chrom}.vep.vcf.gz.csi")

        publishDir "results/vcf", mode: "copy"
	
	"""
	for f in ${annotated_vcfs}; do bcftools index \${f}; done
	for f in ${annotated_vcfs}; do echo "\${f}"; done | sort -V > files.txt
	bcftools concat -a -f files.txt -Oz -o ${name}.${chrom}.vep.vcf.gz
        bcftools index ${name}.${chrom}.vep.vcf.gz
	"""
}


process summarize {
	label "SUMMARY"

	cache "lenient"
 	errorStrategy "retry"
	maxRetries 3
        cpus 1

        input:
        tuple val(chrom), val(name), path(vcf), path(vcf_index)

        output:
        tuple val(name), path("*.summary.txt")

        publishDir "results/summary", pattern: "*.summary.txt", mode: "copy"

        """
	variant_summary.py -a ${vcf} -o ${name}.${chrom}.summary.txt
        """
}


process concatenate_chromosomes {
        label "VEP"

        cache "lenient"
        cpus 1

        input:
        tuple val(name), path(vcfs), path(vcf_indices)

        output:
        tuple path("${name}.vep.vcf.gz"), path("${name}.vep.vcf.gz.csi")
        
        publishDir "results/vcf", mode: "copy"
        
        """ 
        for f in ${vcfs}; do echo "\${f}"; done | sort -V > files.txt
        bcftools concat -a -f files.txt -Oz -o ${name}.vep.vcf.gz
        bcftools index ${name}.vep.vcf.gz
        """
}


process notebooks {
	label "SUMMARY"

	cache "lenient"
	executor "local"
	cpus 1

	input:
	tuple val(name), path(summaries)

	output:
	path "*.html"

	publishDir "results/html", pattern: "*.html", mode: "copy"

	"""
	cp $workflow.projectDir/JNotebooks/*.ipynb .
	jupyter nbconvert --to html autosomal_pass.ipynb --output=${name}.auto_pass.html --execute --no-input --ExecutePreprocessor.timeout=1200
	jupyter nbconvert --to html autosomal_fail.ipynb --output=${name}.auto_fail.html --execute --no-input --ExecutePreprocessor.timeout=1200
	jupyter nbconvert --to html all_pass.ipynb --output=${name}.all_pass.html --execute --no-input --ExecutePreprocessor.timeout=1200
	jupyter nbconvert --to html all_fail.ipynb --output=${name}.all_fail.html --execute --no-input --ExecutePreprocessor.timeout=1200
	"""
}

workflow {
	range = 1..300000000
	step = params.window_size
	windows = Channel.from(range.by(step)).map { w -> [ w, w + step - 1 ] }

	vcfs = Channel.fromPath(params.vcfs).map{ vcf -> [ vcf, vcf + ".tbi" ] }
	vcfs_chunks = vcf_by_chrom(vcfs).flatMap({ chroms, vcf, vcf_index -> chroms.split(',').collect { [it, vcf, vcf_index] }}).toSortedList({a, b -> a[0] <=> b[0]}).flatMap({it}).combine(windows)

	vcfs_chunks_annotated = annotate_chunks(vcfs_chunks).annotations
	vcfs_annotated = concatenate_chunks(vcfs_chunks_annotated.groupTuple(by: [0, 1]))

	concatenate_chromosomes(vcfs_annotated.map{ it -> [it[1], it[2], it[3]] }.groupTuple(by: 0))

	if (params.enable_summary == true) {
		vcfs_summaries = summarize(vcfs_annotated)
		notebooks(vcfs_summaries.groupTuple(by: 0))
	}
}
