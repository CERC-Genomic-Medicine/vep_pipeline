
process vcf_by_chrom {
	executor "local"
	cpus 1

	input:
	set file(vcf), file(vcf_index), val(base) from Channel.fromPath(params.vcfs).map{ vcf -> [ vcf, vcf + ".tbi", vcf.getParent() ] }

	output:
	set stdout, val(base), file(vcf), file(vcf_index) into vcfs

	"""
	tabix -l ${vcf} | tr "\n" ","
	"""
}


windows = Channel.from((1..300000000).by(params.window_size)).map { w -> [ w, w + params.window_size - 1 ] }
vcfs = vcfs.flatMap { chroms, base, vcf, vcf_index -> chroms.split(',').collect { [it, base, vcf, vcf_index] }}
chunks = vcfs.combine(windows)


process annotate_chunks {
	cpus 1

	containerOptions "-B ${params.vep_cache}:/opt/vep/.vep"

	input:
	set val(chrom), val(base), file(vcf), file(vcf_index), val(start), val(stop) from chunks

	output:
	tuple val(chrom), val("${vcf.getBaseName()}"), file("${vcf.getBaseName()}.${chrom}_${start}_${stop}.vep.vcf.gz") into annotated_chunks
	file "*.vep.log"

        publishDir "results/logs", pattern: "*.vep.log", mode: "move"

	script:
	if (params.assembly == "GRCh38")
		"""
		loftee_args=human_ancestor_fa:/opt/vep/.vep/loftee_db_${params.assembly}/human_ancestor.fa.gz,gerp_bigwig:/opt/vep/.vep/loftee_db_${params.assembly}/gerp_conservation_scores.homo_sapiens.GRCh38.bw,conservation_file:/opt/vep/.vep/loftee_db_${params.assembly}/loftee.sql${params.loftee_flags}
		bcftools view -G ${vcf} ${chrom}:${start}-${stop} | vep --cache --offline --assembly ${params.assembly} --format vcf --vcf --compress_output bgzip --force_overwrite --no_stats --dir_cache /opt/vep/.vep/ --plugin LoF,loftee_path:/opt/vep/.vep/loftee_${params.assembly},\${loftee_args} --dir_plugins /opt/vep/.vep/loftee_${params.assembly} ${params.vep_flags} --warning_file STDERR --output_file STDOUT > ${vcf.getBaseName()}.${chrom}_${start}_${stop}.vep.vcf.gz 2> ${vcf.getBaseName()}.${chrom}_${start}_${stop}.vep.log
	 	"""
	else if (params.assembly == "GRCh37")
		"""
		loftee_args=human_ancestor_fa:/opt/vep/.vep/loftee_db_${params.assembly}/human_ancestor.fa.gz,conservation_file:/opt/vep/.vep/loftee_db_${params.assembly}/phylocsf_gerp.sql${params.loftee_flags}
		bcftools view -G ${vcf} ${chrom}:${start}-${stop} | vep --cache --offline --assembly ${params.assembly} --format vcf --vcf --compress_output bgzip --force_overwrite --no_stats --dir_cache /opt/vep/.vep/ --plugin LoF,loftee_path:/opt/vep/.vep/loftee_${params.assembly},\${loftee_args} --dir_plugins /opt/vep/.vep/loftee_${params.assembly} ${params.vep_flags} --warning_file STDERR --output_file STDOUT > ${vcf.getBaseName()}.${chrom}_${start}_${stop}.vep.vcf.gz 2> ${vcf.getBaseName()}.${chrom}_${start}_${stop}.vep.log

		"""
	else
		error "Invalid assembly name: ${params.assembly}"
}


process concatenate_chunks {
	cpus 1

	input:
	set val(chrom), val(base_name), file(annotated_vcfs) from annotated_chunks.groupTuple(by: [0, 1])

        output:
        set file("${base_name}.${chrom}.vep.vcf.gz"), file("${base_name}.${chrom}.vep.vcf.gz.csi") into concatenated

        publishDir "results", mode: "move"
	
	"""
	for f in ${annotated_vcfs}; do bcftools index \${f}; done
	for f in ${annotated_vcfs}; do echo "\${f}"; done | sort -V > files.txt
	bcftools concat -f files.txt -Oz -o ${base_name}.${chrom}.vep.vcf.gz
        bcftools index ${base_name}.${chrom}.vep.vcf.gz
	"""
}
