# optimalcores: Analyse optimal ice core locations in a climate model simulation

## Overview

**optimalcores** is an R software project to analyse the temperature and isotope time series in the isotope-enabled ECHAM5/MPI-OM-wiso past1000 climate model simulation and, specifically, to determine optimal spatial sampling configurations for Antarctic ice cores which maximize the correlation with a target site temperature time series.

The project is subdivided into three main components:
- the `data/` folder provides the climate model data as an R data file,
- the `lib/` folder contains R library functions which provide the main functionality, and
- the `analysis/` folder contains R code for the actual analyses.

The **optimalcores** software is the basis for the results published in Münch, Werner and Laepple, How precipitation intermittency sets an optimal spatial sampling configuration for Antarctic ice cores, xxx. 

All code has been written by [Dr. Thomas Münch](https://www.awi.de/ueber-uns/organisation/mitarbeiter/thomas-muench.html) at the [Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research](https://www.awi.de/). For further information, code enhancements or potential bugs, please write an email or open an issue here. This work was supported by Helmholtz funding through the Polar Regions and Coasts in the Changing Earth System (PACES) programme of the Alfred Wegener Institute.

The original climate model data used here is archived under xxx.

## Getting started

You can start with the analyses in **optimalcores** after just a few steps:

- Download or clone the repository to your machine into a directory of your choice.
- Run the `dependencies.R` to install all packages required for working with **optimalcores**.
- Update the `setup.R`: 
  - Set the `SRCPATH` variable to the direcory into which you copied the **optimalcores** project;
  - Set the `SAVEPATH` variable to a directory where you want to save analysis plots.

## Starting an analysis

Each new analysis using **optimalcores** starts with running the code in `setup.R` and `init.R` by calling `source("setup.R")` followed by `source("init.R")` from within R. The latter step provides all the project's functionality by loading the relevant R packages and the **optimalcores**' function library in `lib/`.

## Extending optimalcores

To extend **optimalcores** by new climate model data you want to analyse, add the data as an R data file to the `data/` folder and update the library function `selectData()` so that the data can be easily loaded. Note that the data need to be [`pField` objects](https://github.com/EarthSystemDiagnostics/pfields) to comply with the **optimalcores** command syntax.
