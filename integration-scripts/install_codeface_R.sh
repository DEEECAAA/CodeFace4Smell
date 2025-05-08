#!/bin/sh
set -e

echo "Installing R base system and development libraries"

# Install base R and core system dependencies
sudo apt-get update -qq
sudo DEBIAN_FRONTEND=noninteractive apt-get -qqy install \
  r-base r-base-dev libx11-dev libssh2-1-dev

# Install the core testthat package
echo "Installing required R packages manually"
Rscript -e 'install.packages("testthat", repos="https://cloud.r-project.org")'

# Determine script and parent directory
SCRIPT_DIR="$( cd "$( dirname "$0" )" && pwd )"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Run additional R package installs from packages.R
if [ -f "$PROJECT_DIR/packages.R" ]; then
  echo "Running R package installer script at $PROJECT_DIR/packages.R"
  Rscript "$PROJECT_DIR/packages.R"
else
  echo "Warning: packages.R not found in $PROJECT_DIR"
fi
