# Mallard multistate model 
Data and code to develop multistate model for mallards in the LMAV. Model illustrates functional connectivity of spatial sanctuaries by mallards within a large protected area network (~12,000 km2)

## Sponsorship
This project is in collaboration with and supported by the Tennessee Wildlife Resources Agency, U.S. Fish and Wildlife Service, Southeastern Region, and Tennessee Technological University. All rights reserved. 

### Folder, Files, and WorkFlow
Files outside of folders are the `.gitignore` file, `.Rhistory`, the `multistate.sh` script to run on the HPC for parallel processing, and the `refuge_connectivity.Rproj` file.

All necessary data and scripts in their associated folders:
1. The "data" folder contains necessary GPS and PAD-US spatial files in .rds, .csv, and .shp formats
2. The "scripts" folder contains 5 required scripts to reproduce analyses.
   - The `DataPrep_MultistateModel.R` prepares daily encounter histories and calculates distance matrices (km) and sanctuary sizes (km2). `all_dat_oct2019_march2023.rds` is the gps file too large to upload. However, the file produces the `datum.rds` and `datum_dist.rds` needed for following scripts.
  - `ModelSpecification_multistate.R` specifies multistate model in R and writes it to .txt file
  - `Functions_MultistateModel.R` prepares custom functions needed for the analyses
  - `Analyses_Mallard_MultistateModel.R` sources scripts, runs model (need `multistate.sh`), loads/processes results, and plots
  - `NumRefsUsedAnalysis.R` fits maximum liklihood truncated Poisson model to evaluate seasonal Refuge use and plots
3. "Results" contain fitted model results `FullMultStatMod_v500i.RData` needed for plotting in `Analyses_Mallard_MultistateModel.R`
4. "Figures" contain figures included and not included in the mansucript


