#!/bin/bash

# Aggiorna sistema
sudo apt update -y && sudo apt upgrade -y

# Installa dipendenze Python
pip install lxml || true

echo "Provisioning completato con successo."