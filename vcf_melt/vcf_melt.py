#!/usr/bin/env python
""" Melt a VCF file into a tab delimited set of calls, one per line

VCF files have all the calls from different samples on one line.  This
script reads vcf on stdin and writes all calls to stdout in tab delimited
format with one call in one sample per line.  This makes it easy to find
a given sample's genotype with, say, grep.
"""

import argparse
import sys
import csv
import collections
import json
import vcf


def flatten(x):
    if type(x) == type([]):
        x = ','.join(map(str, x))
    return x


def parse_annotation(ann, ann_fields, ann_slash_fields):
    """
    """
    
    ann_split = ann.split('|')

    parsed_ann = collections.OrderedDict()
    for idx, field in enumerate(ann_fields):
        if field in ann_slash_fields.keys():
            [f1, f2] = ann_slash_fields[field]
            if ann_split[idx]:
                value_split = ann_split[idx].split('/')
            else:
                value_split = ["",""]
            parsed_ann[f1] = value_split[0]
            parsed_ann[f2] = value_split[1]
        else:
            parsed_ann[field] = ann_split[idx]

    return parsed_ann



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('vcf')
    args = parser.parse_args()

    out = csv.writer(sys.stdout, delimiter=',', quotechar='"')
    if args.vcf:
        inp = open(args.vcf)
    else:
        inp = sys.stdin
    reader = vcf.VCFReader(inp)

    formats = list(reader.formats.keys())

    infos = list(reader.infos.keys())

    has_info_ann = 'ANN' in infos
    has_info_lof = 'LOF' in infos
    has_info_nmd = 'NMD' in infos

    ann_fields = [
        'ANN_ALLELE',
        'ANN_ANNOTATION',
        'ANN_IMPACT',
        'ANN_GENE_NAME',
        'ANN_GENE_ID',
        'ANN_FEATURE_TYPE',
        'ANN_FEATURE_ID',
        'ANN_TRANSCRIPT_BIOTYPE',
        'ANN_RANK',
        'ANN_HGVS_C',
        'ANN_HGVS_P',
        'ANN_CDNA_POSITION_LENGTH',
        'ANN_CDS_POSITION_LENGTH',
        'ANN_AA_POSITION_LENGTH',
        'ANN_DISTANCE',
        'ANN_ERRORS_WARNINGS_INFO',
    ]

    maybe_multiple_alts_info_fields = [
        'AC',
        'AF',
        'AO',
        'PAO',
        'QA',
        'PQA',
        'SAF',
        'SAR',
        'SAP',
        'AB',
        'ABP',
        'RUN',
        'RPP',
        'RPL',
        'RPR',
        'EPP',
        'DPRA',
        'TYPE',
        'CIGAR',
        'MEANALT',
        'MQM',
        'PAIRED',
        'LEN',
        'VAF',
    ]

    maybe_multiple_alts_format_fields = [
        'AD',
        'AO',
        'QA',
    ]

    ref_plus_multiple_alts_format_fields = [
        'GL',
    ]

    ann_slash_fields = {
        'ANN_CDNA_POSITION_LENGTH': ['ANN_CDNA_POSITION', 'ANN_CDNA_LENGTH'],
        'ANN_CDS_POSITION_LENGTH': ['ANN_CDS_POSITION', 'ANN_CDS_LENGTH'],
        'ANN_AA_POSITION_LENGTH': ['ANN_AA_POSITION', 'ANN_AA_LENGTH'],
    }

    snpeff_lof_fields = [
        'lof_gene_name',
        'lof_gene_id',
        'lof_num_transcripts_in_gene',
        'lof_percent_transcripts_affected',
    ]

    snpeff_nmd_fields = [
        'nmd_gene_name',
        'nmd_gene_id',
        'nmd_num_transcripts_in_gene',
        'nmd_percent_transcripts_affected',
    ]

    excluded_formats = {'GQ', 'MIN_DP'}
    formats = list(filter(lambda x: x not in excluded_formats, formats))

    excluded_infos = {'END', 'MIN_DP'}
    infos = list(filter(lambda x: x not in excluded_infos, infos))

    if has_info_ann:
        infos = list(filter(lambda x: x != 'ANN', infos))

    if has_info_lof:
        infos = list(filter(lambda x: x != 'LOF', infos))

    if has_info_nmd:
        infos = list(filter(lambda x: x != 'NMD', infos))

    header = ["SAMPLE"] + ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER'] + ['INFO_' + i for i in infos]

    if has_info_ann:
        header += ann_fields[0:11]
        for slash_field in ann_slash_fields.keys():
            header += ann_slash_fields[slash_field]
        header += ann_fields[14:]

    for f in formats:
        if f in ref_plus_multiple_alts_format_fields:
            header.append('FORMAT_REF_' + f)
            header.append('FORMAT_ALT_' + f)
        else:
            header.append('FORMAT_' + f)

    out.writerow(header)

    for record in reader:
        
        annotations = []
        if has_info_ann:
            annotations = [parse_annotation(x, ann_fields, ann_slash_fields) for x in record.INFO.get('ANN', None)]

        for alt_idx, alt in enumerate(record.ALT):

            if not record.ID:
                record_id = '.'
            fixed = [record.CHROM, record.POS, record_id, record.REF, alt, record.QUAL]
            info = []
            for i in infos:
                if i in maybe_multiple_alts_info_fields:
                    info.append(record.INFO.get(i, None)[alt_idx])
                else:
                    info.append(flatten(record.INFO.get(i, None)))

            for sample in record.samples:
                formats_row = []
                for f in formats:
                    if f in maybe_multiple_alts_format_fields:
                        if type(getattr(sample.data, f, None)) == type([]):
                            formats_row.append(getattr(sample.data, f, None)[alt_idx])
                        else:
                            formats_row.append(getattr(sample.data, f, None))
                    elif f in ref_plus_multiple_alts_format_fields:
                        formats_row.append(getattr(sample.data, f, None)[0])
                        formats_row.append(getattr(sample.data, f, None)[alt_idx + 1])
                    else:
                        formats_row.append(flatten(getattr(sample.data, f, None)))

                if has_info_ann:
                    for ann in annotations:
                        row = [sample.sample]
                        # Format fields not present will simply end up "blank"
                        # in the output
                        row += fixed
                        row += [record.FILTER or '.']
                        row += info
                        row += ann.values()
                        row += formats_row

                        out.writerow(row)
                else:
                    row = [sample.sample]
                    # Format fields not present will simply end up "blank"
                    # in the output
                    row += fixed
                    row += [record.FILTER or '.']
                    row += info
                    row += formats_row

                    out.writerow(row)


if __name__ == '__main__':
    main()
