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


def parse_snpeff_annotation(ann, snpeff_ann_fields, snpeff_ann_slash_fields):
    """
    """
    
    ann_split = ann.split('|')

    parsed_ann = collections.OrderedDict()
    for idx, field in enumerate(snpeff_ann_fields):
        if field in snpeff_ann_slash_fields.keys():
            [f1, f2] = snpeff_ann_slash_fields[field]
            if ann_split[idx]:
                value_split = ann_split[idx].split('/')
            else:
                value_split = ["",""]
            parsed_ann[f1] = value_split[0]
            parsed_ann[f2] = value_split[1]
        else:
            parsed_ann[field] = ann_split[idx]

    return parsed_ann
        


def main(args):

    out = csv.writer(sys.stdout, delimiter=',', quotechar='"')
    if args.vcf:
        inp = open(args.vcf)
    else:
        inp = sys.stdin
    reader = vcf.VCFReader(inp)

    formats = list(reader.formats.keys())
    formats_lookup = {
        'GT':     'genotype',
        'GQ':     'genotype_quality',
        'GL':     'genotype_likelihood',
        'DP':     'read_depth',
        'AD':     'allele_depth',
        'RO':     'ref_allele_observation_count',
        'QR':     'sum_quality_of_ref_observations',
        'AO':     'alt_allele_observation_count',
        'QA':     'sum_quality_of_alt_observations',
        'MIN_DP': 'minimum_depth_gvcf',
    }

    infos = list(reader.infos.keys())
    infos_lookup = {
        'NS':      'num_samples',
        'DP':      'total_depth',
        'DPB':     'total_depth_per_bp',
        'AC':      'alt_alleles',
        'AN':      'num_alleles',
        'AF':      'estimated_allele_freq',
        'RO':      'num_ref_haplotype',
        'AO':      'num_alt_haplotype',
        'PRO':     'num_ref_haplotypes_with_partial',
        'PAO':     'num_alt_haplotypes_with_partial',
        'QR':      'sum_quality_of_ref_observations',
        'QA':      'sum_quality_of_alt_observations',
        'PQR':     'sum_quality_of_partial_ref_observations',
        'PQA':     'sum_quality_of_partial_alt_observations',
        'SRF':     'num_ref_observations_fwd_strand',
        'SRR':     'num_ref_observations_rev_strand',
        'SAF':     'num_alt_observations_fwd_strand',
        'SAR':     'num_alt_observations_rev_strand',
        'SRP':     'strand_balance_probability_ref_allele',
        'SAP':     'strand_balance_probability_alt_allele',
        'AB':      'allele_balance_het_sites',
        'ABP':     'allele_balance_probability_het_sites',
        'RUN':     'run_length',
        'RPP':     'read_placement_probability',
        'RPPR':    'read_placement_probability_ref',
        'RPL':     'reads_placed_left',
        'RPR':     'reads_placed_right',
        'EPP':     'end_placement_probability',
        'EPPR':    'end_placement_probability_ref',
        'DPRA':    'alt_allele_depth_ratio',
        'ODDS':    'log_odds_ratio',
        'GTI':     'genotyping_iterations',
        'TYPE':    'variant_type',
        'CIGAR':   'cigar_alt',
        'NUMALT':  'num_alt_alleles',
        'MEANALT': 'mean_num_alt_alleles',
        'LEN':     'allele_length',
        'MQM':     'mean_mapping_quality_alt',
        'MQMR':    'mean_mapping_quality_ref',
        'PAIRED':  'proportion_properly_paired_alt',
        'PAIREDR': 'proportion_properly_paired_ref',
        'MIN_DP':  'min_depth_gvcf',
        'END':     'end_gvcf',
        'VAF':     'variant_allele_fraction',
        'ANN':     'functional_annotation',
        'LOF':     'loss_of_function',
        'NMD':     'nonsense_mediated_decay',
    }

    has_info_ann = 'ANN' in infos
    has_info_lof = 'LOF' in infos
    has_info_nmd = 'NMD' in infos

    """
    Allele               | Annotation         | Annotation_Impact | Gene_Name                | Gene_ID | Feature_Type           | 
    Feature_ID           | Transcript_BioType | Rank              | HGVS.c                   | HGVS.p  | cDNA.pos / cDNA.length | 
    CDS.pos / CDS.length | AA.pos / AA.length | Distance          | ERRORS / WARNINGS / INFO |
    """
    snpeff_ann_fields = [
        'snpeff_ann_allele',
        'snpeff_ann_annotation',
        'snpeff_ann_annotation_impact',
        'snpeff_ann_gene_name',
        'snpeff_ann_gene_id',
        'snpeff_ann_feature_type',
        'snpeff_ann_feature_id',
        'snpeff_ann_transcript_biotype',
        'snpeff_ann_rank',
        'snpeff_ann_hgvs_c',
        'snpeff_ann_hgvs_p',
        'snpeff_ann_cdna_pos_length',
        'snpeff_ann_cds_pos_length',
        'snpeff_ann_aa_pos_length',
        'snpeff_ann_distance',
        'snpeff_ann_errors_warnings_info',
    ]

    snpeff_ann_slash_fields = {
        'snpeff_ann_cdna_pos_length': ['snpeff_ann_cdna_position', 'snpeff_ann_cdna_length'],
        'snpeff_ann_cds_pos_length': ['snpeff_ann_cds_position', 'snpeff_ann_cds_length'],
        'snpeff_ann_aa_pos_length': ['snpeff_ann_aa_position', 'snpeff_ann_aa_length'],
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

    header = ["library_id"] + ['ref_accession', 'position', 'variant_id', 'ref_allele', 'alt_allele', 'quality', 'filter'] + [infos_lookup[i] for i in infos]

    if has_info_ann:
        header += snpeff_ann_fields[0:11]
        for slash_field in snpeff_ann_slash_fields.keys():
            header += snpeff_ann_slash_fields[slash_field]
        header += snpeff_ann_fields[14:]

    header += [formats_lookup[f] for f in formats]

    out.writerow(header)

    for record in reader:
        info_row = [flatten(record.INFO.get(i, None)) for i in infos]
        annotations = []
        if has_info_ann:
            annotations = [parse_snpeff_annotation(x, snpeff_ann_fields, snpeff_ann_slash_fields) for x in record.INFO.get('ANN', None)]
        
        [alt] = record.ALT
        if not record.ID:
            record_id = '.'
        fixed = [record.CHROM, record.POS, record_id, record.REF, alt, record.QUAL]

        for sample in record.samples:
            if has_info_ann:
                for ann in annotations:
                    row = [sample.sample]
                    # Format fields not present will simply end up "blank"
                    # in the output
                    row += fixed
                    row += [record.FILTER or '.']
                    row += info_row
                    row += ann.values()
                    formats_row = [flatten(getattr(sample.data, f, None)) for f in formats]
                    row += formats_row

                    out.writerow(row)
            else:
                row = [sample.sample]
                # Format fields not present will simply end up "blank"
                # in the output
                row += fixed
                row += [record.FILTER or '.']
                row += info_row
                formats_row = [flatten(getattr(sample.data, f, None)) for f in formats]
                row += formats_row

                out.writerow(row)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('vcf')
    args = parser.parse_args()
    main(args)
