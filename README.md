# Mallard multistate model 
Data and code to develop multistate model for mallards in the LMAV 
Multistate model identifies functional connectivity of spatial sanctuaries
by mallards within a large protected area network (~12,000 ha)

## Sponsorship
This project is in collaboration with and  supported by the
Tennessee Wildlife Resources Agency, U.S. Fish and Wildlife Service, Southeastern Region, 
and Tennessee Technological University. All rights reserved. 

### Folder, Files, and WorkFlow
Files outside of folders are the `.gitignore` file, `.Rhistory`, JAGS model specification file (`fullconstant_MallardMSM`),
the `multistate.sh` bash script to run on the HPC for parallel processing, and the `refuge_connectivity.Rproj` file for 
reproducibility on any device.

All necessary data and scripts in their associated folders:
1. The "data" and "geo_data" folder contains necessary GPS and PAD-US spatial files in .RDS and .shp formats, respectively
2. The "scripts" folder contains 6 required scripts to reproduce analysis in totality
  - `DataPrep_MultistateModel.R` prepares daily encounter histories and calculates distance matrices (km) and sanctuary sizes (km2)
  - `ModelSpecification_multistate.R` specifies multistate model in R and writes it to .txt file
  - `Functions_MultistateModel.R` prepares custom functions needed for the Analyses
  - `Analyses_Mallard_MultistateModel.R` sources scripts, runs model (need `multistate.sh`), loads/processes results, and plots
  - `NumRefsUsedAnalysis.R` fits maximum liklihood Poisson model to evaluate seasonal Refuge use and plots
3. "Results" contain fitted model results `FullMultStatMod_v500i.RData` needed for plotting in `Analyses_Mallard_MultistateModel.R`
4. "Figures" contain Figures (included and not included in the mansucript)


