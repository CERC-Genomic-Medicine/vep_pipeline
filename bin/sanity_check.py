#!/usr/bin/env python

import argparse
import gzip
import re

argparser = argparse.ArgumentParser(description = 'Minimal check of VEP annotated VCFs for the presence of non-truncated CSQ INFO fields.')
argparser.add_argument('-a', '--annotations', metavar = 'file', dest = 'in_VCF_CSQ', required = True, help = 'Input VCF/BCF file with annotations in CSQ INFO field (e.g. from Variant Effect Predictor).')

if __name__ == '__main__':
    args = argparser.parse_args()

    with gzip.open(args.in_VCF_CSQ, 'rt') as ivcf:
        # find meta-information line
        n_csq_columns = 0
        for line in ivcf:
            if line.startswith('##'):
                if line.startswith('##INFO=<ID=CSQ'):
                    n_csq_columns = len(line.split('|'))
            else:
                break
                         
        assert n_csq_columns > 1 # CSQ meta-information line is present
        assert line.startswith('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO') # VCF header line is present and have minimal number of fields

        # scan variants
        for line in ivcf:
            fields = line.split('\t', 8)
            assert len(fields) >= 8 # Minimal number of VCF fields is present
            info_fields = fields[7]
            csq_field = None
            for info_field in info_fields.split(';'):
                if info_field.startswith('CSQ='):
                    csq_field = info_field.split('=', 1)[1]
                    break
            assert csq_field is not None # CSQ INFO annotation is present
            for consequence in csq_field.split(','):
                assert n_csq_columns == len(consequence.split('|')) # All CSQ INFO annoations have declared number of columns

