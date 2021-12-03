def range = 1..300000000
def step = params.window_size


windows = Channel.from(range.by(step)).map { w -> [ w, w + step - 1 ] }


process vcf_by_chrom {
	cache "lenient"
	executor "local"
	cpus 1

	input:
	set file(vcf), file(vcf_index) from Channel.fromPath(params.vcfs).map{ vcf -> [ vcf, vcf + ".tbi" ] }

	output:
	set stdout, file(vcf), file(vcf_index) into vcfs

	"""
	tabix -l ${vcf} | tr "\n" ","
	"""
}


vcfs = vcfs.flatMap({ chroms, vcf, vcf_index -> chroms.split(',').collect { [it, vcf, vcf_index] }}).toSortedList({a, b -> a[0] <=> b[0]}).flatMap({it})
chunks =  vcfs.combine(windows)


process annotate_chunks {
	cache "lenient"
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return "retry" }
	maxRetries 3
	cpus 1

	containerOptions "-B ${params.vep_cache}:/opt/vep/.vep"

	input:
	set val(chrom), file(vcf), file(vcf_index), val(start), val(stop) from chunks

	output:
	tuple val(chrom), val("${vcf.getBaseName()}"), file("${vcf.getBaseName()}.${chrom}_${start}_${stop}.vep.vcf.gz") into annotated_chunks
	file "*.vep.log"

        publishDir "results/logs", pattern: "*.vep.log", mode: "copy"

	script:
	if (params.assembly == "GRCh38")
		"""
		export PERL5LIB=/opt/vep/.vep/Plugins/:$PERL5LIB
		loftee_args=human_ancestor_fa:/opt/vep/.vep/loftee_db_${params.assembly}/human_ancestor.fa.gz,gerp_bigwig:/opt/vep/.vep/loftee_db_${params.assembly}/gerp_conservation_scores.homo_sapiens.GRCh38.bw,conservation_file:/opt/vep/.vep/loftee_db_${params.assembly}/loftee.sql${params.loftee_flags}
		bcftools view -G ${vcf} ${chrom}:${start}-${stop} | vep --cache --offline --assembly ${params.assembly} --format vcf --vcf --compress_output bgzip --force_overwrite --no_stats --dir_cache /opt/vep/.vep/ --plugin LoF,loftee_path:/opt/vep/.vep/loftee_${params.assembly},\${loftee_args} --dir_plugins /opt/vep/.vep/loftee_${params.assembly} --plugin CADD,/opt/vep/.vep/CADD_${params.assembly}/whole_genome_SNVs.tsv.gz,/opt/vep/.vep/CADD_${params.assembly}/InDels.tsv.gz ${params.vep_flags} ${params.vep_flags} --warning_file STDERR --output_file STDOUT > ${vcf.getBaseName()}.${chrom}_${start}_${stop}.vep.vcf.gz 2> ${vcf.getBaseName()}.${chrom}_${start}_${stop}.vep.log
	 	"""
	else if (params.assembly == "GRCh37")
		"""
		export PERL5LIB=/opt/vep/.vep/Plugins/:$PERL5LIB
		loftee_args=human_ancestor_fa:/opt/vep/.vep/loftee_db_${params.assembly}/human_ancestor.fa.gz,conservation_file:/opt/vep/.vep/loftee_db_${params.assembly}/phylocsf_gerp.sql${params.loftee_flags}
		bcftools view -G ${vcf} ${chrom}:${start}-${stop} | vep --cache --offline --assembly ${params.assembly} --format vcf --vcf --compress_output bgzip --force_overwrite --no_stats --dir_cache /opt/vep/.vep/ --plugin LoF,loftee_path:/opt/vep/.vep/loftee_${params.assembly},\${loftee_args} --dir_plugins /opt/vep/.vep/loftee_${params.assembly} --plugin CADD,/opt/vep/.vep/CADD_${params.assembly}/whole_genome_SNVs.tsv.gz,/opt/vep/.vep/CADD_${params.assembly}/InDels.tsv.gz ${params.vep_flags} --warning_file STDERR --output_file STDOUT > ${vcf.getBaseName()}.${chrom}_${start}_${stop}.vep.vcf.gz 2> ${vcf.getBaseName()}.${chrom}_${start}_${stop}.vep.log

		"""
	else
		error "Invalid assembly name: ${params.assembly}"
}


process concatenate_chunks {
	cache "lenient"
	errorStrategy "retry"
	maxRetries 3
	cpus 1

	input:
	set val(chrom), val(base_name), file(annotated_vcfs) from annotated_chunks.groupTuple(by: [0, 1])

        output:
        set file("${base_name}.${chrom}.vep.vcf.gz"), file("${base_name}.${chrom}.vep.vcf.gz.csi") into concatenated

        publishDir "results", mode: "copy"
	
	"""
	for f in ${annotated_vcfs}; do bcftools index \${f}; done
	for f in ${annotated_vcfs}; do echo "\${f}"; done | sort -V > files.txt
	bcftools concat -a -f files.txt -Oz -o ${base_name}.${chrom}.vep.vcf.gz
        bcftools index ${base_name}.${chrom}.vep.vcf.gz
	"""
}
