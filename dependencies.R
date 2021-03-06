##
## aim:
## install all package dependencies for "optimalcores".
## relation:
## https://github.com/EarthSystemDiagnostics/optimalcores
##

install.packages("abind")
install.packages("arrangements")
install.packages("egg")
install.packages("magrittr")

# install.packages("remotes")
remotes::install_github("EarthSystemDiagnostics/geostools")
remotes::install_github("EarthSystemDiagnostics/grfxtools")
remotes::install_github("EarthSystemDiagnostics/pfields")
remotes::install_github("EarthSystemDiagnostics/stattools")

