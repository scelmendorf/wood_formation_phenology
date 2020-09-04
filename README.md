# wood_formation_phenology
Exploring the phenology of xylogenesis, using data from Huang et al. 2020 [doi:10.1073/pnas.2007058117](https://doi.org/10.1073/pnas.2007058117)

## Usage notes
- add the ERA_export_pnas_clim_08_20_2020.csv to your own personal /data directory; *it is too big for git* and will be ignored in the .gitignore
- paths in scripts relative to the top project directory unless otherwise noted


## Directory structure
- data Supp mat tables from Huang et al. 
- scripts
    - cv\_huang Analysis of data using cross-validation
    - gddsims_huang simulations
    - GEE_ERA_extraction.txt Google Earth Engine script to download ERA 5 reanalysis daily climate data for sites    
- plots 
    - figures from simulations and cross validations + some checks on calculations
