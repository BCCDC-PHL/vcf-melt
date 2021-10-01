# vcf-melt
This repo is based on the [`vcf_melt`](https://github.com/jamescasbon/PyVCF/blob/master/scripts/vcf_melt) script from [PyVCF](https://github.com/jamescasbon/PyVCF).
It has been expanded upon to support snpEff-annotated vcf files.

The 'melt' in vcf-melt refers to the process of 'melting' a wide dataset into a long dataset, as is commonly done when preparing [Tidy data](https://r4ds.had.co.nz/tidy-data.html).
A similar procedure is provided by the [`melt`](https://pandas.pydata.org/docs/reference/api/pandas.melt.html) function in the popular [pandas](https://pandas.pydata.org/docs/index.html) python library.

This tool does not completely 'melt' its input, but does resolve several issues with one-to-many relationships that arise when coverting vcf files to csv format:
  - Each variant for each sample is reported on a separate line.
  - Each annotation for each variant is reported on a separate line.
  - In cases where there are multiple alternate alleles for a variant, each alternate allele is reported on a separate line.

## Dependencies
The `vcf_melt.py` script depends on [PyVCF](https://github.com/jamescasbon/PyVCF). Dependencies are listed in the [`setup.py`](https://github.com/BCCDC-PHL/vcf-melt/blob/main/setup.py) file.

Note that PyVCF is not compatible with Python > 3.5. A future version of this script may switch to either manually parsing the vcf contents using only python standard library features, or switching to using [cyvcf2](https://github.com/brentp/cyvcf2) to parse the input.

## Usage
```
usage: vcf-melt [-h] vcf

positional arguments:
  vcf

optional arguments:
  -h, --help  show this help message and exit
```

The `vcf-melt` script takes a vcf file as input and prints the 'melted' csv-formatted output to stdandard output. Eg:

```
vcf-melt a_vcf_file.vcf > melted_vcf.csv
```
## Output

Output is in `.csv` format.

The first field is the sample ID:

| Header   |  Description      |
|:---------|:------------------|
| `SAMPLE` | Sample identifier |

The following 7 fields correspond to 7 of the 8 fixed mandatory columns as defined in the [VCF spec](https://samtools.github.io/hts-specs/):

| Header   |  Description                                                                                                             |
|:---------|:-------------------------------------------------------------------------------------------------------------------------|
| `CHROM`  | Chromosome (or contig) on which the variant occurs                                                                       |  
| `POS`    | Position of the variant, relative to the reference sequence                                                              |
| `ID`     | An identifier for the variant, or `.` if no identifier is available                                                      |
| `REF`    | Reference base(s)                                                                                                        |
| `ALT`    | Alternate base(s)                                                                                                        |
| `QUAL`   | [Phred-scaled](https://en.wikipedia.org/wiki/Phred_quality_score) quality score for the evidence of the alternate allele |
| `FILTER` | `PASS` if the variant passed all quality filters, or `.` if no filters have been applied                                 |

Refer to the [VCF spec](https://samtools.github.io/hts-specs/) for more detailed descriptions of these fields.

The eighth mandatory column (`INFO`) is split out into separate fields, according to the `INFO` fields present in the file. These fields may vary depending on the tool(s) used to produce the `.vcf` file. Some examples:

| Header  | Description                                                                           |
|:--------|:--------------------------------------------------------------------------------------|
| `NS`    | Number of samples with data                                                           |
| `DP`    | Total read depth at the locus                                                         |
| `DPB`   | Total read depth per bp at the locus; bases in reads overlapping / bases in haplotype |
| `AC`    | Total number of alternate alleles in called genotypes                                 |
| `AN`    | Total number of alleles in called genotypes                                           |
| `AF`    | Estimated allele frequency                                                            |

All `INFO` fields are prefixed with `INFO_` in the output `.csv` file.

If the `.vcf` file includes annotation information, the `ANN` field is split out into several separate fields according to the [VCF Annotation Format standard](https://pcingola.github.io/SnpEff/adds/VCFannotationformat_v1.0.pdf).

| Header                  | Description                                                                                                                                   |
|:------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------|
| `ALLELE`                | Allele that the annotation applies to                                                                                                         |
| `ANNOTATION`            | Effect or consequence of the variant, annotated using terms from the [Sequence Ontology](http://www.sequenceontology.org/)                    |
| `IMPACT`                | Estimation of the functional impact of the variant (one of: `HIGH`, `MODERATE`, `LOW`, `MODIFIER`)                                            |
| `GENE_NAME`             | Common gene name                                                                                                                              |
| `GENE_ID`               | Identifier for the gene that the variant occurs in                                                                                            |
| `FEATURE_TYPE`          | Type of feature that is referred to in the `FEATURE_ID` field, using terms from the [Sequence Ontology](http://www.sequenceontology.org/)     |
| `FEATURE_ID`            | Identifier for the sequence feature that the variant occurs in                                                                                |
| `TRANSCRIPT_BIOTYPE`    | Type of transcript (coding, non-coding, etc.), corresponding to [Ensembl biotypes](https://m.ensembl.org/info/genome/genebuild/biotypes.html) |
| `RANK`                  | Exon or Intron rank / total number of exons or introns                                                                                        |
| `HGVS_C`                | Variant using HGVS notation (DNA level)                                                                                                       |
| `HGVS_P`                | If variant is coding, the variant using HGVS notation (Protein level).                                                                        |
| `CDS_POSITION`          | Variant position in coding sequence (CDS)                                                                                                     |
| `CDS_LENGTH`            | Coding sequence (CDS) length                                                                                                                  |
| `CDNA_POSITION`         | Variant position in cDNA sequence                                                                                                             |
| `CDNA_LENGTH`           | cDNA length                                                                                                                                   |
| `AA_POSITION`           | Variant position in amino acid sequence                                                                                                       |
| `AA_LENGTH`             | Amino acid length                                                                                                                             |
| `DISTANCE`              | Distance to feature                                                                                                                           |
| `ERRORS_WARNINGS_INFO`  | Any errors, warnings or info related to variant annotation                                                                                    |

All `ANN` fields are prefixed with `ANN_` in the output `.csv` file.

Note that in the VCF Annotation Format standard, the fields `CDS_POSITION` and `CDS_LENGTH` are combined into a single slash-delimited field, as are the `CDNA` and `AA` position and length fields. They are parsed into separate fields in the `.csv` output of this tool.

As with the `INFO` field, the `FORMAT` field may include several sub-fields, which differ depending on the tool(s) used to produce the `.vcf` file. Some examples:

| Header                    | Description                                            |
|:--------------------------|--------------------------------------------------------|
| `GT`                      | Genotype                                               |
| `REF_GL`                  | Genotype likelihood (log-10 scaled) for the REF allele |
| `ALT_GL`                  | Genotype likelihood (log-10 scaled) for the ALT allele |
| `DP`                      | Read depth                                             |
| `AD`                      | Allele depth                                           |

All `FORMAT` fields are prefixed with `FORMAT_` in the output `.csv` file.

### Example

The following `.vcf` file:
```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##reference=ref.fa
##contig=<ID=MN908947.3,length=29903>
##phasing=none
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus">
##INFO=<ID=DPB,Number=1,Type=Float,Description="Total read depth per bp at the locus; bases in reads overlapping / bases in haplotype">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1]">
##INFO=<ID=RO,Number=1,Type=Integer,Description="Count of full observations of the reference haplotype.">
##INFO=<ID=AO,Number=A,Type=Integer,Description="Count of full observations of this alternate haplotype.">
##INFO=<ID=PRO,Number=1,Type=Float,Description="Reference allele observation count, with partial observations recorded fractionally">
##INFO=<ID=PAO,Number=A,Type=Float,Description="Alternate allele observations, with partial observations recorded fractionally">
##INFO=<ID=QR,Number=1,Type=Float,Description="Reference allele quality sum in phred">
##INFO=<ID=QA,Number=A,Type=Integer,Description="Alternate allele quality sum in phred">
##INFO=<ID=PQR,Number=1,Type=Float,Description="Reference allele quality sum in phred for partial observations">
##INFO=<ID=PQA,Number=A,Type=Float,Description="Alternate allele quality sum in phred for partial observations">
##INFO=<ID=SRF,Number=1,Type=Integer,Description="Number of reference observations on the forward strand">
##INFO=<ID=SRR,Number=1,Type=Integer,Description="Number of reference observations on the reverse strand">
##INFO=<ID=SAF,Number=A,Type=Integer,Description="Number of alternate observations on the forward strand">
##INFO=<ID=SAR,Number=A,Type=Integer,Description="Number of alternate observations on the reverse strand">
##INFO=<ID=SRP,Number=1,Type=Float,Description="Strand balance probability for the reference allele: Phred-scaled upper-bounds estimate of the probability of observing the deviation between SRF and SRR given E(SRF/SRR) ~ 0.5, derived using Hoeffding's inequality">
##INFO=<ID=SAP,Number=A,Type=Float,Description="Strand balance probability for the alternate allele: Phred-scaled upper-bounds estimate of the probability of observing the deviation between SAF and SAR given E(SAF/SAR) ~ 0.5, derived using Hoeffding's inequality">
##INFO=<ID=AB,Number=A,Type=Float,Description="Allele balance at heterozygous sites: a number between 0 and 1 representing the ratio of reads showing the reference allele to all reads, considering only reads from individuals called as heterozygous">
##INFO=<ID=ABP,Number=A,Type=Float,Description="Allele balance probability at heterozygous sites: Phred-scaled upper-bounds estimate of the probability of observing the deviation between ABR and ABA given E(ABR/ABA) ~ 0.5, derived using Hoeffding's inequality">
##INFO=<ID=RUN,Number=A,Type=Integer,Description="Run length: the number of consecutive repeats of the alternate allele in the reference genome">
##INFO=<ID=RPP,Number=A,Type=Float,Description="Read Placement Probability: Phred-scaled upper-bounds estimate of the probability of observing the deviation between RPL and RPR given E(RPL/RPR) ~ 0.5, derived using Hoeffding's inequality">
##INFO=<ID=RPPR,Number=1,Type=Float,Description="Read Placement Probability for reference observations: Phred-scaled upper-bounds estimate of the probability of observing the deviation between RPL and RPR given E(RPL/RPR) ~ 0.5, derived using Hoeffding's inequality">
##INFO=<ID=RPL,Number=A,Type=Float,Description="Reads Placed Left: number of reads supporting the alternate balanced to the left (5') of the alternate allele">
##INFO=<ID=RPR,Number=A,Type=Float,Description="Reads Placed Right: number of reads supporting the alternate balanced to the right (3') of the alternate allele">
##INFO=<ID=EPP,Number=A,Type=Float,Description="End Placement Probability: Phred-scaled upper-bounds estimate of the probability of observing the deviation between EL and ER given E(EL/ER) ~ 0.5, derived using Hoeffding's inequality">
##INFO=<ID=EPPR,Number=1,Type=Float,Description="End Placement Probability for reference observations: Phred-scaled upper-bounds estimate of the probability of observing the deviation between EL and ER given E(EL/ER) ~ 0.5, derived using Hoeffding's inequality">
##INFO=<ID=DPRA,Number=A,Type=Float,Description="Alternate allele depth ratio.  Ratio between depth in samples with each called alternate allele and those without.">
##INFO=<ID=ODDS,Number=1,Type=Float,Description="The log odds ratio of the best genotype combination to the second-best.">
##INFO=<ID=GTI,Number=1,Type=Integer,Description="Number of genotyping iterations required to reach convergence or bailout.">
##INFO=<ID=TYPE,Number=A,Type=String,Description="The type of allele, either snp, mnp, ins, del, or complex.">
##INFO=<ID=CIGAR,Number=A,Type=String,Description="The extended CIGAR representation of each alternaote allele, with the exception that '=' is replaced by 'M' to ease VCF parsing.  Note that INDEL alleles do not have the first matched base (which is provided by default, per the spec) referred to by the CIGAR.">
##INFO=<ID=NUMALT,Number=1,Type=Integer,Description="Number of unique non-reference alleles in called genotypes at this position.">
##INFO=<ID=MEANALT,Number=A,Type=Float,Description="Mean number of unique non-reference allele observations per sample with the corresponding alternate alleles.">
##INFO=<ID=LEN,Number=A,Type=Integer,Description="allele length">
##INFO=<ID=MQM,Number=A,Type=Float,Description="Mean mapping quality of observed alternate alleles">
##INFO=<ID=MQMR,Number=1,Type=Float,Description="Mean mapping quality of observed reference alleles">
##INFO=<ID=PAIRED,Number=A,Type=Float,Description="Proportion of observed alternate alleles which are supported by properly paired read fragments">
##INFO=<ID=PAIREDR,Number=1,Type=Float,Description="Proportion of observed reference alleles which are supported by properly paired read fragments">
##INFO=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum depth in gVCF output block.">
##INFO=<ID=END,Number=1,Type=Integer,Description="Last position (inclusive) in gVCF output record.">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype">
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Number of observation for each allele">
##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Reference allele observation count">
##FORMAT=<ID=QR,Number=1,Type=Float,Description="Sum of quality of the reference observations">
##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observation count">
##FORMAT=<ID=QA,Number=A,Type=Integer,Description="Sum of quality of the alternate observations">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum depth in gVCF output block.">
##INFO=<ID=VAF,Number=A,Type=Float,Description="Variant allele fraction, called from observed reference/alt reads">
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ">
##INFO=<ID=LOF,Number=.,Type=String,Description="Predicted loss of function effects for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected'">
##INFO=<ID=NMD,Number=.,Type=String,Description="Predicted nonsense mediated decay effects for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected'">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample-01
MN908947.3	241	.	C	T	20097.1	.	AB=0;ABP=0;AC=1;AF=1;AN=1;AO=610;CIGAR=1X;DP=610;DPB=610;DPRA=0;EPP=229.071;EPPR=0;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=0;NS=1;NUMALT=1;ODDS=4627.53;PAIRED=1;PAIREDR=0;PAO=0;PQA=0;PQR=0;PRO=0;QA=22681;QR=0;RO=0;RPL=266;RPP=24.6681;RPPR=0;RPR=344;RUN=1;SAF=283;SAP=9.90206;SAR=327;SRF=0;SRP=0;SRR=0;TYPE=snp;VAF=1;ANN=T|upstream_gene_variant|MODIFIER|orf1ab|Gene_265_21554|transcript|QHD43415.1|protein_coding||c.-25C>T|||||25|,T|intergenic_region|MODIFIER|CHR_START-orf1ab|CHR_START-Gene_265_21554|intergenic_region|CHR_START-Gene_265_21554|||n.241C>T||||||	GT:DP:AD:RO:QR:AO:QA:GL	1:610:0,610:0:0:610:22681:-2040.3,0
MN908947.3	3037	.	C	T	17421.0	.	AB=0;ABP=0;AC=1;AF=1;AN=1;AO=515;CIGAR=1X;DP=515;DPB=515;DPRA=0;EPP=163.341;EPPR=0;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=0;NS=1;NUMALT=1;ODDS=4011.34;PAIRED=1;PAIREDR=0;PAO=0;PQA=0;PQR=0;PRO=0;QA=19716;QR=0;RO=0;RPL=206;RPP=47.7426;RPPR=0;RPR=309;RUN=1;SAF=198;SAP=62.7195;SAR=317;SRF=0;SRP=0;SRR=0;TYPE=snp;VAF=1;ANN=T|synonymous_variant|LOW|orf1ab|Gene_265_21554|transcript|QHD43415.1|protein_coding|1/2|c.2772C>T|p.F924F|2772/21291|2772/21291|924/7096||	GT:DP:AD:RO:QR:AO:QA:GL	1:515:0,515:0:0:515:19716:-1773.38,0
MN908947.3	6248	.	A	C	8674.36	.	AB=0;ABP=0;AC=1;AF=1;AN=1;AO=265;CIGAR=1X;DP=265;DPB=265;DPRA=0;EPP=15.4737;EPPR=0;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=0;NS=1;NUMALT=1;ODDS=1997.35;PAIRED=1;PAIREDR=0;PAO=0;PQA=0;PQR=0;PRO=0;QA=10146;QR=0;RO=0;RPL=176;RPP=65.0325;RPPR=0;RPR=89;RUN=1;SAF=36;SAP=308.237;SAR=229;SRF=0;SRP=0;SRR=0;TYPE=snp;VAF=1;ANN=C|missense_variant|MODERATE|orf1ab|Gene_265_21554|transcript|QHD43415.1|protein_coding|1/2|c.5983A>C|p.N1995H|5983/21291|5983/21291|1995/7096||	GT:DP:AD:RO:QR:AO:QA:GL	1:265:0,265:0:0:265:10146:-912.785,0
MN908947.3	14408	.	C	T	11344.9	.	AB=0;ABP=0;AC=1;AF=1;AN=1;AO=333;CIGAR=1X;DP=333;DPB=333;DPRA=0;EPP=30.5613;EPPR=0;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=0;NS=1;NUMALT=1;ODDS=2612.26;PAIRED=1;PAIREDR=0;PAO=0;PQA=0;PQR=0;PRO=0;QA=12733;QR=0;RO=0;RPL=129;RPP=39.6906;RPPR=0;RPR=204;RUN=1;SAF=179;SAP=7.08589;SAR=154;SRF=0;SRP=0;SRR=0;TYPE=snp;VAF=1;ANN=T|missense_variant|MODERATE|orf1ab|Gene_265_21554|transcript|QHD43415.1|protein_coding|2/2|c.14144C>T|p.P4715L|14144/21291|14144/21291|4715/7096||	GT:DP:AD:RO:QR:AO:QA:GL	1:333:0,333:0:0:333:12733:-1145.43,0
MN908947.3	20268	.	A	G	7808.6	.	AB=0;ABP=0;AC=1;AF=1;AN=1;AO=230;CIGAR=1X;DP=230;DPB=230;DPRA=0;EPP=3.95442;EPPR=0;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=0;NS=1;NUMALT=1;ODDS=1798;PAIRED=1;PAIREDR=0;PAO=0;PQA=0;PQR=0;PRO=0;QA=8736;QR=0;RO=0;RPL=131;RPP=12.6781;RPPR=0;RPR=99;RUN=1;SAF=123;SAP=5.42724;SAR=107;SRF=0;SRP=0;SRR=0;TYPE=snp;VAF=1;ANN=G|synonymous_variant|LOW|orf1ab|Gene_265_21554|transcript|QHD43415.1|protein_coding|2/2|c.20004A>G|p.L6668L|20004/21291|20004/21291|6668/7096||,G|upstream_gene_variant|MODIFIER|S|Gene_21562_25383|transcript|QHD43416.1|protein_coding||c.-1295A>G|||||1295|	GT:DP:AD:RO:QR:AO:QA:GL	1:230:0,230:0:0:230:8736:-785.999,0
MN908947.3	23403	.	A	G	17799.4	.	AB=0;ABP=0;AC=1;AF=1;AN=1;AO=530;CIGAR=1X;DP=530;DPB=530;DPRA=0;EPP=124.219;EPPR=0;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=0;NS=1;NUMALT=1;ODDS=4098.47;PAIRED=1;PAIREDR=0;PAO=0;PQA=0;PQR=0;PRO=0;QA=20075;QR=0;RO=0;RPL=253;RPP=5.37024;RPPR=0;RPR=277;RUN=1;SAF=192;SAP=90.3445;SAR=338;SRF=0;SRP=0;SRR=0;TYPE=snp;VAF=1;ANN=G|missense_variant|MODERATE|S|Gene_21562_25383|transcript|QHD43416.1|protein_coding|1/1|c.1841A>G|p.D614G|1841/3822|1841/3822|614/1273||,G|upstream_gene_variant|MODIFIER|ORF3a|Gene_25392_26219|transcript|QHD43417.1|protein_coding||c.-1990A>G|||||1990|,G|upstream_gene_variant|MODIFIER|E|Gene_26244_26471|transcript|QHD43418.1|protein_coding||c.-2842A>G|||||2842|,G|upstream_gene_variant|MODIFIER|M|Gene_26522_27190|transcript|QHD43419.1|protein_coding||c.-3120A>G|||||3120|,G|upstream_gene_variant|MODIFIER|ORF6|Gene_27201_27386|transcript|QHD43420.1|protein_coding||c.-3799A>G|||||3799|,G|upstream_gene_variant|MODIFIER|ORF7a|Gene_27393_27758|transcript|QHD43421.1|protein_coding||c.-3991A>G|||||3991|,G|upstream_gene_variant|MODIFIER|ORF8|Gene_27893_28258|transcript|QHD43422.1|protein_coding||c.-4491A>G|||||4491|,G|upstream_gene_variant|MODIFIER|N|Gene_28273_29532|transcript|QHD43423.2|protein_coding||c.-4871A>G|||||4871|,G|downstream_gene_variant|MODIFIER|orf1ab|Gene_265_21554|transcript|QHD43415.1|protein_coding||c.*1848A>G|||||1848|	GT:DP:AD:RO:QR:AO:QA:GL	1:530:0,530:0:0:530:20075:-1805.71,0
```

..produces the following 'melted' `.csv`:

```csv
SAMPLE,CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO_NS,INFO_DP,INFO_DPB,INFO_AC,INFO_AN,INFO_AF,INFO_RO,INFO_AO,INFO_PRO,INFO_PAO,INFO_QR,INFO_QA,INFO_PQR,INFO_PQA,INFO_SRF,INFO_SRR,INFO_SAF,INFO_SAR,INFO_SRP,INFO_SAP,INFO_AB,INFO_ABP,INFO_RUN,INFO_RPP,INFO_RPPR,INFO_RPL,INFO_RPR,INFO_EPP,INFO_EPPR,INFO_DPRA,INFO_ODDS,INFO_GTI,INFO_TYPE,INFO_CIGAR,INFO_NUMALT,INFO_MEANALT,INFO_LEN,INFO_MQM,INFO_MQMR,INFO_PAIRED,INFO_PAIREDR,INFO_VAF,ANN_ALLELE,ANN_ANNOTATION,ANN_IMPACT,ANN_GENE_NAME,ANN_GENE_ID,ANN_FEATURE_TYPE,ANN_FEATURE_ID,ANN_TRANSCRIPT_BIOTYPE,ANN_RANK,ANN_HGVS_C,ANN_HGVS_P,ANN_CDS_POSITION,ANN_CDS_LENGTH,ANN_AA_POSITION,ANN_AA_LENGTH,ANN_CDNA_POSITION,ANN_CDNA_LENGTH,ANN_DISTANCE,ANN_ERRORS_WARNINGS_INFO,FORMAT_GT,FORMAT_REF_GL,FORMAT_ALT_GL,FORMAT_DP,FORMAT_AD,FORMAT_RO,FORMAT_QR,FORMAT_AO,FORMAT_QA
sample-01,MN908947.3,241,.,C,T,20097.1,.,1,610,610.0,1,1,1.0,0,610,0.0,0.0,0.0,22681,0.0,0.0,0,0,283,327,0.0,9.90206,0.0,0.0,1,24.6681,0.0,266.0,344.0,229.071,0.0,0.0,4627.53,0,snp,1X,1,1.0,1,60.0,0.0,1.0,0.0,1.0,T,upstream_gene_variant,MODIFIER,orf1ab,Gene_265_21554,transcript,QHD43415.1,protein_coding,,c.-25C>T,,,,,,,,25,,1,-2040.3,0.0,610,0,0,0.0,610,22681
sample-01,MN908947.3,241,.,C,T,20097.1,.,1,610,610.0,1,1,1.0,0,610,0.0,0.0,0.0,22681,0.0,0.0,0,0,283,327,0.0,9.90206,0.0,0.0,1,24.6681,0.0,266.0,344.0,229.071,0.0,0.0,4627.53,0,snp,1X,1,1.0,1,60.0,0.0,1.0,0.0,1.0,T,intergenic_region,MODIFIER,CHR_START-orf1ab,CHR_START-Gene_265_21554,intergenic_region,CHR_START-Gene_265_21554,,,n.241C>T,,,,,,,,,,1,-2040.3,0.0,610,0,0,0.0,610,22681
sample-01,MN908947.3,3037,.,C,T,17421.0,.,1,515,515.0,1,1,1.0,0,515,0.0,0.0,0.0,19716,0.0,0.0,0,0,198,317,0.0,62.7195,0.0,0.0,1,47.7426,0.0,206.0,309.0,163.341,0.0,0.0,4011.34,0,snp,1X,1,1.0,1,60.0,0.0,1.0,0.0,1.0,T,synonymous_variant,LOW,orf1ab,Gene_265_21554,transcript,QHD43415.1,protein_coding,1/2,c.2772C>T,p.F924F,2772,21291,2772,21291,924,7096,,,1,-1773.38,0.0,515,0,0,0.0,515,19716
sample-01,MN908947.3,6248,.,A,C,8674.36,.,1,265,265.0,1,1,1.0,0,265,0.0,0.0,0.0,10146,0.0,0.0,0,0,36,229,0.0,308.237,0.0,0.0,1,65.0325,0.0,176.0,89.0,15.4737,0.0,0.0,1997.35,0,snp,1X,1,1.0,1,60.0,0.0,1.0,0.0,1.0,C,missense_variant,MODERATE,orf1ab,Gene_265_21554,transcript,QHD43415.1,protein_coding,1/2,c.5983A>C,p.N1995H,5983,21291,5983,21291,1995,7096,,,1,-912.785,0.0,265,0,0,0.0,265,10146
sample-01,MN908947.3,14408,.,C,T,11344.9,.,1,333,333.0,1,1,1.0,0,333,0.0,0.0,0.0,12733,0.0,0.0,0,0,179,154,0.0,7.08589,0.0,0.0,1,39.6906,0.0,129.0,204.0,30.5613,0.0,0.0,2612.26,0,snp,1X,1,1.0,1,60.0,0.0,1.0,0.0,1.0,T,missense_variant,MODERATE,orf1ab,Gene_265_21554,transcript,QHD43415.1,protein_coding,2/2,c.14144C>T,p.P4715L,14144,21291,14144,21291,4715,7096,,,1,-1145.43,0.0,333,0,0,0.0,333,12733
sample-01,MN908947.3,20268,.,A,G,7808.6,.,1,230,230.0,1,1,1.0,0,230,0.0,0.0,0.0,8736,0.0,0.0,0,0,123,107,0.0,5.42724,0.0,0.0,1,12.6781,0.0,131.0,99.0,3.95442,0.0,0.0,1798.0,0,snp,1X,1,1.0,1,60.0,0.0,1.0,0.0,1.0,G,synonymous_variant,LOW,orf1ab,Gene_265_21554,transcript,QHD43415.1,protein_coding,2/2,c.20004A>G,p.L6668L,20004,21291,20004,21291,6668,7096,,,1,-785.999,0.0,230,0,0,0.0,230,8736
sample-01,MN908947.3,20268,.,A,G,7808.6,.,1,230,230.0,1,1,1.0,0,230,0.0,0.0,0.0,8736,0.0,0.0,0,0,123,107,0.0,5.42724,0.0,0.0,1,12.6781,0.0,131.0,99.0,3.95442,0.0,0.0,1798.0,0,snp,1X,1,1.0,1,60.0,0.0,1.0,0.0,1.0,G,upstream_gene_variant,MODIFIER,S,Gene_21562_25383,transcript,QHD43416.1,protein_coding,,c.-1295A>G,,,,,,,,1295,,1,-785.999,0.0,230,0,0,0.0,230,8736
sample-01,MN908947.3,23403,.,A,G,17799.4,.,1,530,530.0,1,1,1.0,0,530,0.0,0.0,0.0,20075,0.0,0.0,0,0,192,338,0.0,90.3445,0.0,0.0,1,5.37024,0.0,253.0,277.0,124.219,0.0,0.0,4098.47,0,snp,1X,1,1.0,1,60.0,0.0,1.0,0.0,1.0,G,missense_variant,MODERATE,S,Gene_21562_25383,transcript,QHD43416.1,protein_coding,1/1,c.1841A>G,p.D614G,1841,3822,1841,3822,614,1273,,,1,-1805.71,0.0,530,0,0,0.0,530,20075
sample-01,MN908947.3,23403,.,A,G,17799.4,.,1,530,530.0,1,1,1.0,0,530,0.0,0.0,0.0,20075,0.0,0.0,0,0,192,338,0.0,90.3445,0.0,0.0,1,5.37024,0.0,253.0,277.0,124.219,0.0,0.0,4098.47,0,snp,1X,1,1.0,1,60.0,0.0,1.0,0.0,1.0,G,upstream_gene_variant,MODIFIER,ORF3a,Gene_25392_26219,transcript,QHD43417.1,protein_coding,,c.-1990A>G,,,,,,,,1990,,1,-1805.71,0.0,530,0,0,0.0,530,20075
sample-01,MN908947.3,23403,.,A,G,17799.4,.,1,530,530.0,1,1,1.0,0,530,0.0,0.0,0.0,20075,0.0,0.0,0,0,192,338,0.0,90.3445,0.0,0.0,1,5.37024,0.0,253.0,277.0,124.219,0.0,0.0,4098.47,0,snp,1X,1,1.0,1,60.0,0.0,1.0,0.0,1.0,G,upstream_gene_variant,MODIFIER,E,Gene_26244_26471,transcript,QHD43418.1,protein_coding,,c.-2842A>G,,,,,,,,2842,,1,-1805.71,0.0,530,0,0,0.0,530,20075
sample-01,MN908947.3,23403,.,A,G,17799.4,.,1,530,530.0,1,1,1.0,0,530,0.0,0.0,0.0,20075,0.0,0.0,0,0,192,338,0.0,90.3445,0.0,0.0,1,5.37024,0.0,253.0,277.0,124.219,0.0,0.0,4098.47,0,snp,1X,1,1.0,1,60.0,0.0,1.0,0.0,1.0,G,upstream_gene_variant,MODIFIER,M,Gene_26522_27190,transcript,QHD43419.1,protein_coding,,c.-3120A>G,,,,,,,,3120,,1,-1805.71,0.0,530,0,0,0.0,530,20075
sample-01,MN908947.3,23403,.,A,G,17799.4,.,1,530,530.0,1,1,1.0,0,530,0.0,0.0,0.0,20075,0.0,0.0,0,0,192,338,0.0,90.3445,0.0,0.0,1,5.37024,0.0,253.0,277.0,124.219,0.0,0.0,4098.47,0,snp,1X,1,1.0,1,60.0,0.0,1.0,0.0,1.0,G,upstream_gene_variant,MODIFIER,ORF6,Gene_27201_27386,transcript,QHD43420.1,protein_coding,,c.-3799A>G,,,,,,,,3799,,1,-1805.71,0.0,530,0,0,0.0,530,20075
sample-01,MN908947.3,23403,.,A,G,17799.4,.,1,530,530.0,1,1,1.0,0,530,0.0,0.0,0.0,20075,0.0,0.0,0,0,192,338,0.0,90.3445,0.0,0.0,1,5.37024,0.0,253.0,277.0,124.219,0.0,0.0,4098.47,0,snp,1X,1,1.0,1,60.0,0.0,1.0,0.0,1.0,G,upstream_gene_variant,MODIFIER,ORF7a,Gene_27393_27758,transcript,QHD43421.1,protein_coding,,c.-3991A>G,,,,,,,,3991,,1,-1805.71,0.0,530,0,0,0.0,530,20075
sample-01,MN908947.3,23403,.,A,G,17799.4,.,1,530,530.0,1,1,1.0,0,530,0.0,0.0,0.0,20075,0.0,0.0,0,0,192,338,0.0,90.3445,0.0,0.0,1,5.37024,0.0,253.0,277.0,124.219,0.0,0.0,4098.47,0,snp,1X,1,1.0,1,60.0,0.0,1.0,0.0,1.0,G,upstream_gene_variant,MODIFIER,ORF8,Gene_27893_28258,transcript,QHD43422.1,protein_coding,,c.-4491A>G,,,,,,,,4491,,1,-1805.71,0.0,530,0,0,0.0,530,20075
sample-01,MN908947.3,23403,.,A,G,17799.4,.,1,530,530.0,1,1,1.0,0,530,0.0,0.0,0.0,20075,0.0,0.0,0,0,192,338,0.0,90.3445,0.0,0.0,1,5.37024,0.0,253.0,277.0,124.219,0.0,0.0,4098.47,0,snp,1X,1,1.0,1,60.0,0.0,1.0,0.0,1.0,G,upstream_gene_variant,MODIFIER,N,Gene_28273_29532,transcript,QHD43423.2,protein_coding,,c.-4871A>G,,,,,,,,4871,,1,-1805.71,0.0,530,0,0,0.0,530,20075
sample-01,MN908947.3,23403,.,A,G,17799.4,.,1,530,530.0,1,1,1.0,0,530,0.0,0.0,0.0,20075,0.0,0.0,0,0,192,338,0.0,90.3445,0.0,0.0,1,5.37024,0.0,253.0,277.0,124.219,0.0,0.0,4098.47,0,snp,1X,1,1.0,1,60.0,0.0,1.0,0.0,1.0,G,downstream_gene_variant,MODIFIER,orf1ab,Gene_265_21554,transcript,QHD43415.1,protein_coding,,c.*1848A>G,,,,,,,,1848,,1,-1805.71,0.0,530,0,0,0.0,530,20075
```