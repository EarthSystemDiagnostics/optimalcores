##
## aim:
## install all package dependencies for "optimalcores".
## relation:
## https://github.com/EarthSystemDiagnostics/optimalcores
##

install.packages("abind")
install.packages("arrangements")
install.packages("magrittr")
install.packages("RColorBrewer")

devtools::install_github("EarthSystemDiagnostics/ecustools")
devtools::install_github("EarthSystemDiagnostics/pfields")

