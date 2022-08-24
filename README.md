Dissertation Project
   ---------------------
   Thesis, code and supporting documents.
   
   Author name and contact: 
   -----------------------
   Kate Griffin, kate.griffin21@imperial.ac.uk

   Brief description:
   ------------------
   This repository contains all the data and code used to produce this thesis. The results are also provided in .pdf form.

   Languages: 
   ---------
   R (v 3.6.4), LaTeX
   
   Dependencies:
   -------------
   tidyverse, dpplyr, ggplot2, phytools, cowplot & gridExtra


   Project structure and Usage: 
   ---------------------------
   Code: This folder contains (1) all the code used in this analysis and (2) source code for the LaTeX write.
   
   1. get_phylogeny_from_OTL_and_timetree.R: Used to get a time calibrated tree used for tree_script.R (run this before running tree_script.R)
   2. TPC_fitting_pipeline_all_in.R: Used to fit TPC to data (run before running tree_script.R)
   3. tree_script.R: main code
   4. compare_species.R: additional script for plots showing different E estimates for (i) juvenilles vs adults and (ii) AMR vs RMR
   5. Thesis folder: Contains all the code and supporting documents for the project write up 

   
   Data: This folder contains the metabolic rate data set that was produced for this thesis. 
   
   Results: Folder with thesis and supplementary information in .pdf form
