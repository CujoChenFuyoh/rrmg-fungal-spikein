# Base-resolved gene evidence (0.3.4)

Step3/4 now create the inputs needed to investigate nucleotide differences in
contigs, metagenomic reads and synthetic target-only reads. They are observations,
not a validated FungAMR caller, resistance phenotype or automatic false-positive
decision. Existing coverage R figures are not changed by this update.

## Existing project: add evidence without rebuilding analysis

Activate the same analysis environment, set the absolute code/project paths,
then run on a compute node (or a local machine with the required files):

```bash
conda activate fungamr_meta
CODE=/path/to/rrmg-fungal-spikein
WD=/path/to/SpikeIN
bash "$CODE/rrmg_backfill_cyp51a_evidence.sh" -w "$WD" -r 44 -t 4
bash "$CODE/rrmg_backfill_cyp51a_evidence.sh" -w "$WD" -r all -t 8
column -t -s $'\t' "$WD/result/CYP51A_evidence_backfill_status.tsv" | head
bash "$CODE/rrmg_collect_figure_data.sh" --wd "$WD" \
  --gene "$WD/ref/CYP51A_wt.fasta" --require-base-evidence
```

`-r all` follows actual `run_index` values, including a non-contiguous plan.
Rows run serially; `-t` applies to contig minimap2 only. Missing files are logged
per source and row, with nonzero exit status if any row failed. Existing derived
files with matching fingerprints/checksums are skipped; `--force` rebuilds only
the derived evidence. Status lives in `result/`, outside managed `figure_data/`.

Inputs per run:

- `asm/contigs.fasta`
- `readmap/reads_to_GENE.bam`
- `readmap/target_only_reads_to_GENE.bam`

The reference defaults to `ref/CYP51A_wt.fasta`; override with `--reference`.
BAM indexing is not required for the sequential pileup. No FASTQ, fastp,
assembly or read remapping is performed. Original PAFs, BAMs, mapping summaries
and AFR files are not changed. For zero-target designs only, a missing target
BAM is replaced by explicit zero evidence, never for a positive-target design.

## New runs and resume

- Step3 uses one `minimap2 -x asm5 -c --cs=long` gene alignment for both new
  base-resolved evidence and the normal contig summary. The derived filename is
  `asm/contigs.fasta.to_GENE.cs.paf`; the existing PAF filename is preserved.
- A completed Step3 rerun only fills the derived cs PAF; it does not replace the
  old summary/PAF or trigger QC/assembly. Old unaudited assemblies still require
  the explicit backfill path, which labels generation as UNVERIFIED.
- Step4 adds raw and Q20 evidence for each read set, even when its mapping is
  already complete. Missing evidence alone never triggers FASTQ remapping.
  Use backfill if the old FASTQs have already been deleted.
- Every derived output has an `.evidence.json` sidecar with reference sequence
  SHA256, source path/size/mtime, generation, tool version, options and output
  checksum. This detects ordinary source changes, not malicious content edits
  that preserve size/mtime. A new fingerprint cannot prove which reference was
  originally used to create a legacy BAM; retain the original logs/checksums.

## Read-count definition

Each prefix has `.base_counts.raw.tsv.gz` and `.base_counts.q20.tsv.gz`.
Raw uses MAPQ/baseQ 0/0; Q20 uses 20/20. Both use `mpileup -A -B -x --ff 0xF04
-d 0 -aa`: count anomalous pairs, disable BAQ, count each overlapping mate,
exclude unmapped/secondary/supplementary/QC-fail/duplicate-flagged alignments,
and remove the depth cap. These are read-level, not molecule-level counts.
This is a diagnostic evidence policy, not a variant-calling filter validation.

Tables include all reference positions, strand-specific A/C/G/T/N counts,
REF/non-REF counts, non-REF fraction and anchored insertion/deletion alleles.
No coverage remains `ObservationStatus=NO_COVERAGE`, never a WT call. Indels
are separate from SNV counts; N is retained but is not treated as an ALT allele.

## Portable collector outputs

Under `figure_data/CYP51A_alignment/`:

- each `row_*/` retains the legacy plotting inputs and adds the complete cs PAF,
  four base-count tables and their checksum sidecars;
- `base_evidence_status.tsv`: presence and recorded generation for every source;
- `contig_recovery.tsv`: a separate recovery definition, not a replacement for
  legacy `HasHit`;
- `contig_differences.tsv`: all aligned contig substitutions/indels, including
  low-quality alignments explicitly flagged `PassRecoveryFilter=N`.

Default legacy collection is retained; inspect `base_evidence_complete` in both
COMPLETE.json files. Use `--require-base-evidence` to reject any missing/invalid
source before accepting a snapshot. It proves evidence completeness, not that
legacy generation is verified. Rerun collection only when no other process is
writing these inputs, and download the whole bundle.

## Recovery and review rules

Qualifying alignments: MAPQ >=20, alignment block >=150 bp, identity >=95%.
Merge target intervals separately for each contig, then across contigs:

| RecoveryClass | Definition |
| --- | --- |
| FULL_SINGLE | One qualifying contig, single-contig breadth >=90% |
| MULTI_CONTIG_WITH_FULL | Multiple qualifying contigs, at least one >=90% |
| FRAGMENTED_COMPLETE | No single contig >=90%, combined breadth >=90% |
| FRAGMENTED_PARTIAL | Qualifying contig(s), combined breadth <90% |
| NOT_DETECTED | No qualifying alignment |

These are explicit initial recovery filters, not validated species specificity
or SNV detection limits. A 90% recovered gene is not necessarily a complete ORF.
`DiscordanceFlag=Y` means combined contig breadth >=90% and read breadth <50%
or a >=50 percentage-point contig/read gap. This is a review flag, never an
automatic false-positive label. `GenerationStatus=MISMATCH/UNVERIFIED` prevents
claiming a paired assembly/read result. Current-read support, background
homology and provenance must be investigated before calling a false positive.

cs difference positions are 1-based relative to the gene reference, oriented
to that reference even for reverse-strand contigs. Insertions are anchored to
the preceding reference base (0 for a leading insertion). Deletions start at
the first deleted base. Coordinates are not amino-acid positions multiplied by
three. Known FungAMR annotation requires validated CDS/codon coordinates and
actual observed codons; covering a labelled site is not detecting its ALT allele.

Technical definitions: [minimap2 cs format](https://github.com/lh3/minimap2#the-cs-optional-tag),
[samtools mpileup](https://www.htslib.org/doc/samtools-mpileup.html).
