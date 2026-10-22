# mTRACE analysis pipeline

## Workflow overview

1. Split raw-Read QC, adapter trimming, UMI processing, alignment, deduplication and coverage-track generation.
2. Call peaks in WT samples to define endogenous m6A-enriched regions; restore strand information to MACS2 peaks and merge WT regions. 
3. Remove WT-m6A-overlapping regions from experimental BAM files.
4. Repeat strand-aware peak calling on WT-m6A-excluded BAM files.
5. Assess peak replicate reproducibility with IDR and apply a fold-enrichment threshold.
6. Identify sample-specific peaks by excluding peaks overlapping control peak sets.
7. Annotate retained peak summits to genomic features and genes; quantify exon-level reads with featureCounts.

Reference files expected by the scripts include the hg38 genome/index, a matching hg38 GTF, feature-annotation BED files, and a gene-body BED file. Ensure all reference files use a consistent genome build and chromosome naming convention.

hg38 genome: `Homo_sapiens.GRCh38.dna.primary_assembly.fa (https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz)`

hg38 gtf: `Homo_sapiens.GRCh38.109.gtf (https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz)`

feature-annotation BED files: `reference/hg38_features_IntergenicExcldTSS500bp.bed`

gene BED file: `reference/hg38_gene.bed`

## Scripts and steps
### `01-reads_trim_and_align.sh` — QC, trimming, alignment and deduplication

### `02.1-WT_stranded_switch.sh` — Split WT BAMs by strand and normalize flags

The standard MACS2 callpeak output does not reliably provide the transcript strand needed for this strand-specific analysis. Therefore, the paired-end BAMs are split into positive-/forward- and negative-/reverse- strand groups using the paired-read SAM flags before peak calling. The script selects flag combinations 99/147 and 83/163, merges the appropriate pairs into positive- and negative-strand BAMs, and then changes the second-in-pair flags to their first-in-pair equivalents before MACS2 is run.

This R2-to-R1 flag conversion is used to make reads from both mates available to the strand-specific peak-calling procedure rather than relying only on reads carrying an R1 flag. Specifically, 163 is changed to 99 and 147 is changed to 83. This is a flag-field reassignment for the downstream peak-calling workflow, not a change to the underlying sequence or genomic coordinates; use it only with the same strand-assignment logic and validate the resulting BAMs if adapting the workflow.

### `02.2-WT_peakcalling.sh` — Call WT peaks and build a WT peak mask

### `03.1-SampleBam_exclude_WTm6Aregion.sh` — Remove WT peak regions from sample BAMs

### `03.2-excldWT_stranded_switch.sh` — Strand splitting and flag normalization after WT exclusion

Repeats the strand-splitting and R2-to-R1 flag conversion procedure from `02.1`, but uses the WT-excluded BAMs as input. 

### `03.3-excldWTsample_peakcalling.sh` — Call peaks in WT-excluded samples

### `04-excldWTsample_peak_replicate_process.sh` — Replicate reproducibility with IDR

### `05-excldWTsample_getUniqPeak.sh` — Identify sample-specific peaks and count reads

#### Required downstream step: `limma` filtering of significant peaks

### `06-RBP_bindingSite_annotation.sh` — Annotate peak summits

feature-annotated peak tables, parsed feature labels, gene-annotation tables, and gene lists.

The Python script reads a tab-delimited BED-like annotation table, examines the comma-separated feature labels in column 5 (zero-based index 4), and assigns one label using this priority: `mRNA exon > ncRNA exon > intron > TSS_upstream500 > start_codon > stop_codon > Intergenic`.

If all three mRNA exon categories (`five_prime_utr`, `CDS`, and `three_prime_utr`) overlap, the site is labelled `mRNA_exons_mix`; if one or more of those categories overlaps together with exon, it is labelled `exons_mix`.

### `07-featureCounts.sh` — Gene-level exon counts
