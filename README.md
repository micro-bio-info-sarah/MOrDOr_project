# MOrDOr_project
Seed­‐carried Microbiota from Orobanche during the early parasitic (Orobanche) cycle

Folders:
- Raw_fastq: sequencing output as raw fastq files
- Snakemake_workflow: fastq files assigned to 16S/ITS with jason files needed for Snakemake workflow (https://gitlab.univ-nantes.fr/bird_pipeline_registry/microSysMics)
- Qiime2_output: output files obtained from Raw_data by Qiime2 workflow 
- Krona: Krona charts created from seed data

1. Import and subset data: Import_subset_data.R
2. Estimate sequencing depth: seq_depth
3. Calculate prevalence and abundance of ASV: prev_abund.R
4. Estimate alpha-diversity: alpha_div.R
5. Estimate beta-diversity: beta-div.R
6. Taxonomic composition overview: taxo_compo.R
