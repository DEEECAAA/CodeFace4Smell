#!/bin/sh
set -e

echo "Installing R base system and development libraries"

# Install base R and core system dependencies
sudo apt-get update -qq
sudo DEBIAN_FRONTEND=noninteractive apt-get -qqy install \
  r-base r-base-dev libx11-dev libssh2-1-dev

# Tutti gli altri pacchetti R saranno installati da packages.R
if [ -f "packages.R" ]; then
  echo "Running R package installer script"
  sudo Rscript packages.R
else
  echo "Warning: packages.R not found"
fi
