#!/bin/sh
set -e

sudo apt-get update -qq
sudo DEBIAN_FRONTEND=noninteractive apt-get -qqy install software-properties-common curl gnupg lsb-release

# Aggiungi repository ufficiale CRAN
echo "Adding CRAN R repository"
sudo mkdir -p /etc/apt/keyrings
curl -fsSL https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | gpg --dearmor | sudo tee /etc/apt/keyrings/cran-archive-keyring.gpg >/dev/null
echo "deb [signed-by=/etc/apt/keyrings/cran-archive-keyring.gpg] https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" | \
  sudo tee /etc/apt/sources.list.d/cran.list > /dev/null

# Aggiungi Node.js da NodeSource
echo "Adding Node.js 18 repository"
curl -fsSL https://deb.nodesource.com/setup_18.x | sudo -E bash -
sudo apt-get install -y nodejs  # 👈 INSTALLAZIONE MANCANTE

# Aggiorna repository
sudo apt-get update -qq
