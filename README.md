This repository contains the custom-made scripts used for the work presented in "A pangenomic approach identifies new structural variants associated with gene expression in the CARTaGENE cohort" by Laflamme et al.

"update_ids.py" is the script used to modify the SV ids in a VCF.

"Nextflow_conditional_analysis" contains the nextflow pipeline used for the conditional eQTL analysis.

"permutations.1apr26.sh" is the script used to validate the enrichment of SV-eQTLs overlaping regulatory elements (by creating sets of randomly selected matched SVs from those included in the eQTL analysis and testing their overlap with regulatory elements).

"OpenTargets_query.py" is the script used to query small variants (those in high LD with the SV-eQTLs) in the Open Targets datasets.

"figures.R" is the script used to generate the figures (main and suppplementary). 
"chm13.lengths_chr.txt" and "hg38.chr_sizes.txt" are used with this script.
