# Spatiotemporal evolution of SARS-CoV-2 in the Bangkok metropolitan region, Thailand, 2020—2022: implications for future outbreak preparedness

18/12/2023

Description
==============
R scripts used to generate figures and results in the manuscript entitled 
"Spatiotemporal evolution of SARS-CoV-2 in the Bangkok metropolitan region, Thailand, 2020—2022: implications for future outbreak preparedness"

The analyses were run using the following R packages:
==============
adephylo==1.1.13,
ape==5.7.1,
car==3.1.2,
cowplot==1.1.1,
dplyr==1.0.9,
emmeans==1.6.1,
forcats==1.0.0,
ggbreak==0.1.1,
ggforce==0.4.1,
ggimage==0.3.1,
ggnewscale==0.4.9,
ggplot2==3.4.0,
ggplotify==0.1.2,
ggpmisc==0.5.2,
ggtree==3.5.0.901,
lsmeans==2.30.0,
lubridate==1.9.0,
sf==1.0.14,
tibble==3.1.7,
tidytree==0.4.5,
treeio==1.21.0,
wesanderson==0.3.6,

in R 4.0.4

INSTALLATION INSTRUCTIONS
==============
Download the "code" folder, and put it somewhere appropriate. 
If necessary, install required packages in R.

INSTRUCTIONS FOR USE
==============
Upzip the file "tha_admbnda_adm1_rtsd_20220121.zip" in the "data/Th_map" directory, and leave the file "tha_admbnda_adm1_rtsd_20220121.shp" there. 
Run the scrtips in R. 
To execute the codes, set the variable "path_to_wd" to the path of the "code" directory. 
Results will be available to you in "code/out".
