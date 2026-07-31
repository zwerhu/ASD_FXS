# FXS-ASD Transcriptome-to-Drug-Repurposing Pipeline

Integrative transcriptome-wide association and Connectivity Map analysis
nominating candidate repurposable drugs for Fragile X syndrome (FXS) and
autism spectrum disorder (ASD).

This repository contains the complete analysis pipeline used to: (1)
meta-analyze population case-control ASD GWAS; (2) test genetic overlap
between population ASD risk and monozygotic (MZ) twin symptom-discordance
GWAS; (3) impute tissue-specific gene expression via S-PrediXcan across 49
GTEx v8 tissues for 8 GWAS traits, including a plasma pQTL for the Fragile X
mental retardation protein (FMRP); and integrate association signal across
tissues and test Hallmark pathway enrichment; and (4) query the CMap LINCS
L1000 Level 5 dataset for small molecules predicted to computationally
reverse the disease-associated transcriptional signatures, cross-referenced
against mechanism-of-action annotations.

## Pipeline overview

<img width="361" height="509" alt="image" src="https://github.com/user-attachments/assets/307b5f63-3c2b-4ade-89ec-254cc2e1e49d" />


## Repository structure

```
.
├── 01_data_preparation/         GWAS-VCF -> TSV conversion
├── 02_meta_analysis/            Fixed-effects IVW meta-analysis
├── 03_genetic_correlation/      LD-aware genetic overlap testing
├── 04_transcriptome_imputation/ S-PrediXcan (GTEx v8 MASHR), incl. numpy2 patches
├── 05_integrative_analysis/     Cross-tissue ACAT integration + Hallmark GSEA
├── 06_drug_repurposing/         CMap L1000 connectivity mapping + MOA enrichment
└── requirements.txt
```

Scripts within each stage are numbered in the order they are meant to be
run. Every script takes explicit `--input`/`--output`-style command-line
arguments (no hardcoded paths), so they can be chained with any directory
layout you prefer.

## Data sources

| Trait (short name) | Source | Access |
|---|---|---|
| ASD_caseControl_iPSYCH2017 | Grove et al. 2019, iPSYCH-PGC | [OpenGWAS ieu-a-1185](https://gwas.mrcieu.ac.uk/datasets/ieu-a-1185/) |
| ASD_FinnGen_strict | FinnGen R5, KRA_PSY_AUTISM_EXMORE | [OpenGWAS](https://gwas.mrcieu.ac.uk) |
| MZtwin_ASD_all / MZtwin_ASD_children | Assary et al., Within Families Consortium | [OpenGWAS ieu-b-5151 / ieu-b-5153](https://gwas.mrcieu.ac.uk) |
| FMRP_protein | Sun et al. 2018, INTERVAL pQTL | [OpenGWAS prot-a-1156](https://gwas.mrcieu.ac.uk) |
| ASD_vs_Tourette / ASD_vs_OCD / ASD_vs_MDD | Peyrot & Robinson 2021, CC-GWAS | [GWAS Catalog](https://www.ebi.ac.uk/gwas/) (GCST90016603 / 90016602 / 90016614) |
| GTEx v8 MASHR models | GTEx Consortium 2020 | [predictdb.org](https://predictdb.org) |
| CMap LINCS L1000 Level 5 | Broad Institute | [GSE92742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742) |
| Drug MOA annotations | Broad Repurposing Hub (2020-03-24) | mirrored at [github.com/ncats/drug_rep](https://github.com/ncats/drug_rep) |

**Note on GWAS dataset verification:** one dataset in the GWAS Catalog
(GCST90016603) is superficially labelled "autism-related" but is in fact a
case-case (CC-GWAS) contrast between ASD and Tourette syndrome cases, not a
population case-control GWAS. This pipeline treats it accordingly (see
`03_genetic_correlation/` and the manuscript Discussion, Section 4.3, for
why this distinction matters and what happens if it is missed). Always
verify a GWAS dataset's actual case/control definition against its source
publication, not just its Catalog/OpenGWAS trait label.

## Installation

```bash
pip install -r requirements.txt

# S-PrediXcan is not distributed via PyPI:
git clone https://github.com/hakyimlab/MetaXcan.git
cd MetaXcan/software
patch -p2 < ../../04_transcriptome_imputation/patches/numpy2_compatibility.patch
```

See [`04_transcriptome_imputation/patches/README.md`](04_transcriptome_imputation/patches/README.md)
for details on why this patch is required under NumPy >= 2.0.

## Usage

### 1. Convert GWAS-VCF to TSV

```bash
python 01_data_preparation/batch_convert_vcf.py \
    --input_dir raw_vcf/ --output_dir tables/ \
    --converter 01_data_preparation/vcf_to_table.py
```

### 2. Meta-analyse the population case-control ASD GWAS

```bash
python 02_meta_analysis/meta_analysis_ivw.py \
    --study1 tables/ieu-a-1185.table.gz \
    --study2 tables/finn-b-KRA_PSY_AUTISM_EXMORE.table.gz \
    --output meta_groupA_ASD.tsv.gz
```

### 3. Test genetic overlap with MZ-twin symptom discordance

```bash
# Pruned genome-wide + naive (uncorrected) suggestive-loci test
python 03_genetic_correlation/ld_pruned_correlation.py \
    --reference meta_groupA_ASD.tsv.gz \
    --comparison tables/ieu-b-5151.table.gz \
    --label ieu-b-5151

# LD-corrected independent lead-loci test (use this number, not the naive one above)
python 03_genetic_correlation/independent_lead_loci_test.py \
    --reference meta_groupA_ASD.tsv.gz \
    --comparison tables/ieu-b-5151.table.gz \
    --label ieu-b-5151
```

### 4. Transcriptome imputation

```bash
# Build the rsID<->varID map once (reusable across all traits)
python 04_transcriptome_imputation/build_rsid_to_varid_map.py \
    --models_dir mashr_models/mashr --output rsid_to_varid_map.tsv.gz

# Recode each GWAS trait
python 04_transcriptome_imputation/recode_gwas_to_varid.py \
    --study ieu-a-1185 --mapping rsid_to_varid_map.tsv.gz

# Run S-PrediXcan across all 49 tissues
bash 04_transcriptome_imputation/run_spredixcan_batch.sh \
    MetaXcan/software mashr_models/mashr tables_recoded spredixcan_out \
    ieu-a-1185 finn-b-KRA_PSY_AUTISM_EXMORE ieu-b-5151 ieu-b-5153 \
    prot-a-1156 ebi-a-GCST90016602 ebi-a-GCST90016614 ebi-a-GCST90016603
```

### 5. Cross-tissue integration and pathway enrichment

```bash
python 05_integrative_analysis/cross_tissue_acat_integration.py \
    --input spredixcan_ALL_RESULTS.tsv.gz \
    --mhc_gene_list mhc_gene_list.txt \
    --output integrative_ACAT_results.tsv.gz

python 05_integrative_analysis/hallmark_gsea.py \
    --integrative_results integrative_ACAT_results.tsv.gz \
    --genesets_json genesets.json \
    --combined_output gsea_ALL_hallmark_results.tsv.gz
```

### 6. Drug repurposing via CMap connectivity mapping

```bash
python 06_drug_repurposing/00_prepare_disease_signatures.py \
    --integrative_results integrative_ACAT_results.tsv.gz \
    --cmap_gene_info GSE92742_Broad_LINCS_gene_info.txt.gz \
    --output disease_signatures_for_cmap.tsv.gz \
    --gene_list_output overlap_gene_list.txt

# Run on a machine with ~25-45GB free disk (see script docstring)
python 06_drug_repurposing/01_merge_filter_gctx.py \
    --sig_info GSE92742_Broad_LINCS_sig_info.txt.gz \
    --sig_metrics GSE92742_Broad_LINCS_sig_metrics.txt.gz \
    --gene_info GSE92742_Broad_LINCS_gene_info.txt.gz \
    --gene_list overlap_gene_list.txt \
    --output_matrix level5_filtered_matrix.tsv.gz

python 06_drug_repurposing/02_connectivity_analysis.py \
    --disease_signatures disease_signatures_for_cmap.tsv.gz \
    --cmap_matrix level5_filtered_matrix.parquet \
    --output cmap_drug_repurposing_results.tsv.gz

python 06_drug_repurposing/03_significance_and_moa_enrichment.py \
    --connectivity_results cmap_drug_repurposing_results.tsv.gz \
    --moa_table broad_moa.txt \
    --output_significant cmap_significant_drugs.tsv \
    --output_moa_enrichment cmap_moa_enrichment.tsv

python 06_drug_repurposing/04_drug_disease_gene_tissue_links.py \
    --significant_drugs cmap_significant_drugs.tsv \
    --moa_table broad_moa.txt \
    --disease_signatures disease_signatures_for_cmap.tsv.gz \
    --cmap_matrix level5_filtered_matrix.parquet \
    --integrative_results integrative_ACAT_results.tsv.gz \
    --output cmap_drug_disease_gene_tissue.tsv
```


## Key findings

- Corrected ASD case-control meta-analysis (iPSYCH-PGC 2017 + FinnGen,
  N=46,350): 30 genome-wide significant SNPs anchored at **chr20q13.13**.
- Cross-tissue integration nominates ***XRN2*** as the top ASD-associated
  gene and ***PPP3CA*** (calcineurin catalytic subunit) as the top gene
  associated with the FMRP plasma pQTL signal, independent of any ASD
  case-control comparison.
- No evidence of shared common-variant architecture between population ASD
  risk and MZ-twin symptom discordance once linkage-disequilibrium-driven
  inflation is corrected for (naive analysis: p<1e-12; LD-corrected
  independent-locus analysis: p>0.4).
- 2,788 candidate reversal compounds nominated across 8 traits; glutamate
  receptor antagonism, angiotensin receptor antagonism, and T-type calcium
  channel blockade show the strongest (nominal, not FDR-significant)
  mechanism-of-action enrichment.

See the manuscript (coming soon) for full methods, results,
and discussion, including a detailed account of two methodological
corrections made during this analysis (a mischaracterized CC-GWAS trait, and
LD-driven inflation of a naive genetic-overlap test) that we believe
generalize as cautions for similar integrative genomics pipelines.

## Known limitations

- The genetic-overlap analysis in `03_genetic_correlation/` is an LD-pruning
  proxy, **not** a substitute for LD Score Regression (LDSC).
- The CMap connectivity analysis aggregates signatures to one representative
  profile per compound across all cell lines by default; a neural-lineage
  restricted re-analysis (`--cell_lines NPC NEU` in
  `01_merge_filter_gctx.py`) is recommended as a follow-up.
- MOA annotation covers only ~9% of the tested compound library (limited by
  the Broad Repurposing Hub's coverage of named vs. unnamed probe
  compounds), limiting MOA enrichment power.
- All input GWAS are of predominantly European ancestry.

## Citation

If you use this pipeline, please cite the accompanying manuscript (to be posted) and the original data sources listed above.

## License

Code in this repository is released under the MIT License. GWAS summary
statistics, GTEx models, and CMap data retain their original data-use
licenses; consult each source's terms before redistribution.
