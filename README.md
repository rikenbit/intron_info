# intron_info_from_gtf.R

`intron_info_from_gtf.R` is the R script to extract and summarise information on intron length.

## Dependency

- `tidyverse`
- `rtracklayer`

## Example usage

```R
path_gtf <- "https://bioinformatics.riken.jp/ramdaq/ramdaq_annotation/mouse/gencode.vM27.primary_assembly.annotation.ERCC.gtf"
path_outdir <- "intron_info/gencode.vM27"

source("src/intron_info_from_gtf.R")
```

### Example of very short intron

- <https://asia.ensembl.org/Mus_musculus/Transcript/Exons?db=core;g=ENSMUSG00000050122;r=1:37192930-37226686;t=ENSMUST00000169057>


## Outputs

### Transcript support level (TSL)

> TSL Categories
> The following categories are assigned to each of the evaluated annotations:
> 
> tsl1 – all splice junctions of the transcript are supported by at least one non-suspect mRNA
> tsl2 – the best supporting mRNA is flagged as suspect or the support is from multiple ESTs
> tsl3 – the only support is from a single EST
> tsl4 – the best supporting EST is flagged as suspect
> tsl5 – no single transcript supports the model structure
> tslNA – the transcript was not analysed for one of the following reasons:
> pseudogene annotation, including transcribed pseudogenes
> human leukocyte antigen (HLA) transcript
> immunoglobin gene transcript
> T-cell receptor transcript
> single-exon transcript (will be included in a future version)

<http://asia.ensembl.org/info/genome/genebuild/transcript_quality_tags.html>

## Contact

Haruka Ozaki <harukao.cb_at_gmail.com>
