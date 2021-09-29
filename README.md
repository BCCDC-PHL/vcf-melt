# vcf-melt
This repo is based on the [`vcf_melt`](https://github.com/jamescasbon/PyVCF/blob/master/scripts/vcf_melt) script from [PyVCF](https://github.com/jamescasbon/PyVCF).
It has been expanded upon to support snpEff-annotated vcf files.

The 'melt' in vcf-melt refers to the process of 'melting' a wide dataset into a long dataset, as is commonly done when preparing [Tidy data](https://r4ds.had.co.nz/tidy-data.html).
A similar procedure is provided by the [`melt`](https://pandas.pydata.org/docs/reference/api/pandas.melt.html) function in the popular [pandas](https://pandas.pydata.org/docs/index.html) python library.

This tool does not completely 'melt' its input, but does resolve several issues with one-to-many relationships that arise when coverting vcf files to csv format:
  - Each variant for each sample is reported on a separate line.
  - Each annotation for each variant is reported on a separate line.
  - In cases where there are multiple alternate alleles for a variant, each alternate allele is reported on a separate line.

## Usage
```
usage: vcf_melt.py [-h] vcf

positional arguments:
  vcf

optional arguments:
  -h, --help  show this help message and exit
```

The `vcf_melt.py` script takes a vcf file as input and prints the 'melted' csv-formatted output to stdandard output. Eg:

```
vcf_melt.py a_vcf_file.vcf > melted_vcf.csv
```
