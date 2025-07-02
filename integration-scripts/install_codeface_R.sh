#!/bin/bash
set -e

echo "▶ Updating package lists..."
sudo apt-get update -qq

echo "▶ Fixing held or broken packages..."
sudo apt-get install -y --allow-downgrades pkg-config || {
  echo "❌ Failed to install pkg-config (possibly due to pkgconf conflict)."
  exit 1
}

echo "▶ Installing R base and system dependencies for R packages..."
sudo apt-get install -y \
  r-base r-base-dev libx11-dev libssh2-1-dev \
  libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev \
  libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev \
  libglpk-dev libgsl-dev libudunits2-dev libgdal-dev libgeos-dev \
  libproj-dev graphviz libgraphviz-dev libmysqlclient-dev || {
    echo "❌ Failed installing one or more system dependencies."
    exit 1
  }

echo "✅ System dependencies installed."

echo "▶ Creating user R library and setting R_LIBS_USER"
R_LIBS_USER="/home/vagrant/R/x86_64-pc-linux-gnu-library/4.5"
mkdir -p "$R_LIBS_USER"
echo "R_LIBS_USER=$R_LIBS_USER" >> /home/vagrant/.Renviron
chown vagrant:vagrant /home/vagrant/.Renviron

echo "▶ Installing core R packages from CRAN..."
Rscript -e 'install.packages(c(
  "testthat", "devtools", "pbkrtest", "lme4", "nloptr",
  "animation", "magick", "pkgdown", "ragg", "textshaping", "svglite"
), lib="'"$R_LIBS_USER"'", repos="https://cloud.r-project.org")' || {
  echo "❌ Failed to install base CRAN packages"
  exit 1
}

echo "▶ Installing additional project dependencies from CRAN..."
Rscript -e 'install.packages(c(
  "igraph", "logging", "reshape", "optparse", "lubridate",
  "colorspace", "RMySQL", "xts", "dtw", "gridExtra"
), lib="'"$R_LIBS_USER"'", repos="https://cloud.r-project.org")' || {
  echo "❌ Failed to install additional CRAN packages"
  exit 1
}

echo "▶ Installing Bioconductor packages..."
Rscript -e '
Sys.setenv(R_LIBS_USER="'"$R_LIBS_USER"'")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib=Sys.getenv("R_LIBS_USER"))
BiocManager::install(c("graph", "Rgraphviz", "BiRewire"), lib=Sys.getenv("R_LIBS_USER"), ask=FALSE, update=FALSE)' || {
  echo "❌ Failed to install Bioconductor packages"
  exit 1
}

PACKAGE_SCRIPT="/vagrant/packages.R"
echo "▶ Checking for project-specific package file at: $PACKAGE_SCRIPT"

if [ -f "$PACKAGE_SCRIPT" ]; then
  echo "▶ Running packages.R..."
  Rscript "$PACKAGE_SCRIPT" || {
    echo "❌ Failed running packages.R"
    exit 1
  }
  echo "✅ Finished packages.R"
else
  echo "❌ packages.R not found at $PACKAGE_SCRIPT"
  ls -la /vagrant
  exit 1
fi

echo "✅ All R packages installed and configured successfully."
