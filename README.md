# MOrDOr_project
Seed­‐carried Microbiota from Orobanche during the early parasitic (Orobanche) cycle

Folders:
- Genotype: SSR data and analysis
- Phenotype: data and analysis of RUC files
- Raw_fastq: sequencing output as raw fastq files
- Snakemake_workflow: fastq files assigned to 16S/ITS with jason files needed for Snakemake workflow (https://gitlab.univ-nantes.fr/bird_pipeline_registry/microSysMics)
- Qiime2_output: output files obtained from Snakemake_workflow fastq by Qiime2 workflow 
- Krona: Krona charts created from seed data

Scripts:
1. Import and subset data: Import_subset_data.R
2. Estimate sequencing depth: Seq_depth.R
3. Calculate prevalence and abundance of ASV: Prev_abund.R
4. Estimate alpha-diversity: Alpha_div.R
5. Estimate beta-diversity: Beta_div.R
6. Taxonomic composition overview: Taxo_compo.R
7. Create upset plot: Upset_plot.R
