# cagablea
Monitoring SARS-CoV-2 variants alterations in Nice neighborhoods by wastewater nanopore sequencing
https://www.medrxiv.org/content/10.1101/2021.07.09.21257475v1

Background: Wastewater surveillance was proposed as an epidemiological tool to define the prevalence and evolution of the SARS-CoV-2 epidemics. However, most implemented SARS-CoV-2 wastewater surveillance projects was based on qPCR measurement of virus titers and did not address the mutational spectrum of SARS-CoV-2 circulating in the population. 

Methods: We have implemented a nanopore RNA sequencing monitoring system in the city of Nice (France, 550,000 inhabitants). Between October 2020 and March 2021, we monthly analyzed the SARS-CoV-2 variants in 113 wastewater samples collected in the main wastewater treatment plant and 20 neighborhoods.

Findings: We initially detected the lineages predominant in Europe at the end of 2020 (B.1.160, B.1.177, B.1.367, B.1.474, and B.1.221). In January, a localized emergence of a variant (Spike:A522S) of the B.1.1.7 lineage occurred in one neighborhood. It rapidly spread and became dominant all over the city. Other variants of concern (B.1.351, P.1) were also detected in some neighborhoods, but at low frequency. Comparison with individual clinical samples collected during the same week showed that wastewater sequencing correctly identified the same lineages as those found in COVID-19 patients.

Interpretation: Wastewater sequencing allowed to document the diversity of SARS-CoV-2 sequences within the different neighborhoods of the city of Nice. Our results illustrate how sequencing of sewage samples can be used to track pathogen sequence diversity in the current pandemics and in future infectious disease outbreaks.


The project details the bioinformatic approaches used in this study. 

It contains the following folders:

data_base: contains the files needed to build a database.
----------------------------------------------------------
agg_data_lineage_All_2019-12-15-2020-12-31_top_2000.csv
 is a file that was downloaded from covidcg webpage. it contains the information about the alterations describing every lineage. the total of 200 lineages are described. 
 
 covidCG_lineages.R
 Containes the code that was used to format data from covidcg project and upload them into the database in the format where every line corresponds to an alteration and a name of a lineage it is connected to.
 
 description_to_SQL.R
 Demands a file in .csv format with the samples description, transfers information from it to the database
 
populate_SQL.R
Formats and inserts the ivar output to created before Ww_ivar table in the database.
It takes a path to a folder with the ivar output as an argument
The decoding table should be created before addind ivar data.

wastewater: contains the files to build a shiny application. 
-------------------------------------------------------------
server.R
defines the functions needed to vizualize plots and maps displayed by the application

ui.R
frontend part, defines the appearance of the application

ww_samples_analysis: contains the scripts used to analyse the data.
--------------------------------------------------------------------- 
contains folders: heatmaps, quality research and variant calling
heatmap_VOC_function_upd.R is the code used to build the heatmap shown in the paper
my_strain_col.csv is the file used by the function. It binds the mutations with their lineages.

ivar_bash.sh in variant_calling folder is the script to perform variant calling. The code takes as an argument the path to the folder containing .bam files of one run.
Sars_cov_2.ASM985889v3.101.gff3 contains the information about sars-cov-2 protein coding regions
sars-cov2.fa is the reference genome used in the studies

summary_lineage_composition.R was used to summarize the percentage of several choosen lineages in the sample.


 
