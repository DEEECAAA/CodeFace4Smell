# Verifica e installa strumenti base
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Bioconductor
BiocManager::install(c("BiRewire", "graph", "Rgraphviz"), ask = FALSE)

# CRAN
install.packages(c("statnet", "logging", "corrgram"), dependencies = TRUE)

# GitHub (pacchetti privati presumibilmente)
devtools::install_github("wolfgangmauerer/snatm/pkg")
devtools::install_github("wolfgangmauerer/tm-plugin-mail/pkg")

# Altro GitHub
devtools::install_github("wch/shiny-gridster")