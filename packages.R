# Verifica devtools e BiocManager
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

filter.installed.packages <- function(packageList)  {
    if("-f" %in% commandArgs(trailingOnly = TRUE)) {
        return(packageList)
    } else {
        return(setdiff(packageList, rownames(installed.packages())))
    }
}

remove.installed.packages <- function(pack) {
    for (path in .libPaths()) {
        tryCatch({
            remove.packages(pack, lib = path)
            message(paste("Removed package:", pack))
        }, error = function(e) {})
    }
}

reinstall.package.from.github <- function(package, url) {
    if (package %in% rownames(installed.packages())) {
        remove.installed.packages(package)
    }
    devtools::install_github(url, upgrade = "never")
}

# Bioconductor packages
bioc_pkgs <- filter.installed.packages(c("BiRewire", "graph", "Rgraphviz"))
if (length(bioc_pkgs) > 0) {
    BiocManager::install(bioc_pkgs, ask = FALSE)
}

# CRAN packages
cran_pkgs <- filter.installed.packages(c(
    "statnet", "ggplot2", "tm", "optparse", "igraph", "zoo", "xts", "lubridate", "xtable",
    "reshape", "wordnet", "stringr", "yaml", "plyr", "scales", "gridExtra", "RMySQL",
    "RCurl", "mgcv", "shiny", "dtw", "httpuv", "png", "rjson", "lsa", "RJSONIO",
    "corrgram", "logging", "testthat"
))
if (length(cran_pkgs) > 0) {
    install.packages(cran_pkgs, dependencies = TRUE)
}

# GitHub packages
reinstall.package.from.github("tm.plugin.mail", "wolfgangmauerer/tm-plugin-mail/pkg")
reinstall.package.from.github("snatm", "wolfgangmauerer/snatm/pkg")
reinstall.package.from.github("shinyGridster", "wch/shiny-gridster")
reinstall.package.from.github("shinybootstrap2", "rstudio/shinybootstrap2")
reinstall.package.from.github("Rgraphviz", "mitchell-joblin/Rgraphviz")