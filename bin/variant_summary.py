#!/usr/bin/env python

import argparse
import pysam
from collections import Counter


argparser = argparse.ArgumentParser(description = 'Variant level summary from annotated VCF/BCF files.')
argparser.add_argument('-a', '--annotations', metavar = 'file', dest = 'in_VCF_CSQ', required = True, help = 'Input VCF/BCF file with annotations in CSQ INFO field (e.g. from Variant Effect Predictor). Requires tabix index. Genotypes are not required.')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'out_file', required = 'True', help = 'Name for the tab-delimited output file.')

# Consequences in order of severity from http://useast.ensembl.org/info/genome/variation/predicted_data.html#consequence_type_table
CONSEQUENCES = [
    'transcript_ablation',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'stop_gained',
    'frameshift_variant',
    'stop_lost',
    'start_lost',
    'transcript_amplification',
    'inframe_insertion',
    'inframe_deletion',
    'missense_variant',
    'protein_altering_variant',
    'splice_region_variant',
    'incomplete_terminal_codon_variant',
    'stop_retained_variant',
    'start_retained_variant',
    'synonymous_variant',
    'coding_sequence_variant',
    'mature_miRNA_variant',
    '5_prime_UTR_variant',
    '3_prime_UTR_variant',
    'non_coding_transcript_exon_variant',
    'intron_variant',
    'NMD_transcript_variant',
    'non_coding_transcript_variant',
    'upstream_gene_variant',
    'downstream_gene_variant',
    'TFBS_ablation',
    'TFBS_amplification',
    'TF_binding_site_variant',
    'regulatory_region_ablation',
    'regulatory_region_amplification',
    'feature_elongation',
    'regulatory_region_variant',
    'feature_truncation',
    'intergenic_variant',
    'ALL',
]
CONSEQUENCES_SEVERITY = {consequence:i for i, consequence in enumerate(CONSEQUENCES)}


INTERGENIC_CONSEQUENCES = {
    'downstream_gene_variant',
    'upstream_gene_variant',
    'intergenic_variant',
    'regulatory_region_variant',
    'regulatory_region_ablation',
    'regulatory_region_amplification',
    'TF_binding_site_variant',
    'TFBS_ablation',
    'TFBS_amplification',
    'intergenic_variant'
}


CONTEXT_TYPES = [
    'CODING',
    'NONCODING',
    'INTERGENIC',
    'ALL',
]
CONTEXT_TYPES_PRIORITY = {context_type:i for i, context_type in enumerate(CONTEXT_TYPES)}


VARIANT_FILTERS = ['PASS', 'FAIL']
VARIANT_TYPES = ['ALL', 'SNV', 'INDEL']


counts = {}
for variant_filter in VARIANT_FILTERS:
    counts[variant_filter] = {}
    for variant_type in VARIANT_TYPES:
        counts[variant_filter][variant_type] = {}
        for context_type in CONTEXT_TYPES:
            counts[variant_filter][variant_type][context_type] = {}
            for variant_category in ['ALL', 'ALL_KNOWN', 'ALL_NOVEL', 'SINGLETONS', 'SINGLETONS_KNOWN', 'SINGLETONS_NOVEL']:
                counts[variant_filter][variant_type][context_type][variant_category] = {}
                for consequence in CONSEQUENCES:
                    counts[variant_filter][variant_type][context_type][variant_category][consequence] = Counter()


def read_vcf(filename):
    with pysam.VariantFile(filename) as ifile:
        if not 'CSQ' in ifile.header.info:
            print(f'Missing CSQ INFO field meta-information.')
            sys.exit(1)
        vep_field_names = ifile.header.info['CSQ'].description.split(':', 1)[-1].strip().split('|')
        for record in ifile:
            allele_effects = {}
            for allele_effect in record.info['CSQ']:
                allele_effect = allele_effect.split('|')
                assert len(vep_field_names) == len(allele_effect), (vep_field_names, allele_effect)
                allele_effect = dict(zip(vep_field_names, allele_effect))
                allele_effects.setdefault(int(allele_effect['ALLELE_NUM']) - 1, []).append(allele_effect)
            variant_filter  = 'PASS' if 'PASS' in record.filter else 'FAIL'
            for i, alt_allele in enumerate(record.alts):
                if record.info['AN'] == 0: # skip if all genotypes are missing
                    continue
                if record.info['AC'][i] == 0: # skip monomorphic
                    continue
                variant_types = [ 'ALL' ] 
                variant_types.append('SNV' if len(alt_allele) == 1 and len(record.ref) == 1 else 'INDEL')
                variant_categories = [ 'ALL' ]
                if record.info['AC'][i] == 1:
                    variant_categories.append('SINGLETONS')                
                seq_5mer = set()
                rsIds = set()
                allele_consequences = set()
                if i not in allele_effects:
                    continue
                for allele_effect in  allele_effects[i]:
                    if allele_effect['SEQ_5MER']:
                        seq_5mer.add(allele_effect['SEQ_5MER'])
                    if allele_effect['Existing_variation']:
                        rsIds.add(allele_effect['Existing_variation'])
                    for x in allele_effect['Consequence'].split('&'):
                        if allele_effect['BIOTYPE'] == '' or x in INTERGENIC_CONSEQUENCES:
                            sequence_type = 'INTERGENIC'
                        elif allele_effect['BIOTYPE'] == 'protein_coding':
                            sequence_type = 'CODING'
                        else:
                            sequence_type = 'NONCODING'
                        allele_consequences.add((x, sequence_type))
                most_severe_consequence = min(allele_consequences, key = lambda x: (CONSEQUENCES_SEVERITY[x[0]], CONTEXT_TYPES_PRIORITY[x[1]]))
                variant_consequences = [ most_severe_consequence, (most_severe_consequence[0], 'ALL'), ('ALL', most_severe_consequence[1]), ('ALL', 'ALL') ]
                if rsIds:
                    variant_categories += [ f'{x}_KNOWN' for x in variant_categories]
                else:
                    variant_categories += [ f'{x}_NOVEL' for x in variant_categories]
                if 'SNV' in variant_types:
                    assert len(seq_5mer) == 1
                    seq_5mer = seq_5mer.pop()
                    isCpG = 'CG' in seq_5mer[1:4]
                else:
                    assert len(seq_5mer) == 0
                    isCpG = None
                for variant_type in variant_types:
                    for variant_consequence in variant_consequences:
                        for variant_category in variant_categories:
                            if variant_type == 'SNV':
                                count_types = [f'{record.ref}_{alt_allele}']
                                if isCpG:
                                    count_types.append('CpG')
                            else:
                                count_types = ['N']
                            for count_type in count_types:
                                counts[variant_filter][variant_type][variant_consequence[1]][variant_category][variant_consequence[0]][count_type] += 1


if __name__ == '__main__':
    args = argparser.parse_args()
    read_vcf(args.in_VCF_CSQ)
    with open(args.out_file, 'wt') as ofile:
        ofile.write('FILTER\tVARIANT_TYPE\tCONTEXT_TYPE\tCONSEQUENCE\tCOUNT\t{}\n'.format('\t'.join(['ALL', 'ALL_KNOWN', 'ALL_NOVEL', 'SINGLETONS', 'SINGLETONS_KNOWN', 'SINGLETONS_NOVEL'])))
        for variant_filter in VARIANT_FILTERS:
            for variant_type in VARIANT_TYPES:
                for context_type in CONTEXT_TYPES:
                    for consequence in CONSEQUENCES:
                        n = []
                        cpg = []
                        ts = []
                        for variant_category in  ['ALL', 'ALL_KNOWN', 'ALL_NOVEL', 'SINGLETONS', 'SINGLETONS_KNOWN', 'SINGLETONS_NOVEL']:
                            consequence_counts = counts[variant_filter][variant_type][context_type][variant_category][consequence]
                            n.append(sum(consequence_counts.values()) - consequence_counts['CpG'])
                            cpg.append(consequence_counts['CpG'])
                            ts.append(sum(consequence_counts[x] for x in ['A_G', 'G_A', 'C_T', 'T_C']))
                        ofile.write('{}\t{}\t{}\t{}\tN\t{}\n'.format(variant_filter, variant_type, context_type, consequence, '\t'.join(str(x) for x in n)))
                        if variant_type == 'SNV':
                            ofile.write('{}\t{}\t{}\t{}\tCpG\t{}\n'.format(variant_filter, variant_type, context_type, consequence, '\t'.join(str(x) for x in cpg)))
                            ofile.write('{}\t{}\t{}\t{}\tTS\t{}\n'.format(variant_filter, variant_type, context_type, consequence, '\t'.join(str(x) for x in ts)))
                        else:
                            assert sum(ts) == 0

               


