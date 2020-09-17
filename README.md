# wood_formation_phenology
Exploring the phenology of xylogenesis, using data from Huang et al. 2020 [doi:10.1073/pnas.2007058117](https://doi.org/10.1073/pnas.2007058117). Methodological details can be found in the methods.pdf document.

## Usage notes
- To run the cross validation you will first need to add the ERA_export_pnas_clim_08_20_2020.csv to your own personal /data directory; *it is too big for git* and will be ignored in the .gitignore
    - It can be downloaded from [here](https://drive.google.com/file/d/165Owv5gtYrBYTmPJXwPGJx9TSAGP8iLo/view?usp=sharing)
- paths in scripts relative to the top project directory unless otherwise noted


## Directory structure
- data Supp mat tables from Huang et al. 
- scripts
    - cv\_huang Analysis of data using cross-validation
    - gddsims_huang simulations
    - GEE_ERA_extraction.txt Google Earth Engine script to download ERA 5 reanalysis daily climate data for sites    
- plots 
    - figures from simulations and cross validations + some checks on calculations
- output & calculations generate some intermediary files, can speed up rerunning code when tweaking figures
