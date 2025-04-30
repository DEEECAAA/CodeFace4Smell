#!/bin/sh
set -e

echo "▶ Setting up Node.js dependencies for id_service..."

cd /vagrant

# Verifica che la cartella id_service esista
if [ ! -d "id_service" ]; then
  echo "❌ Error: 'id_service' directory not found!"
  exit 1
fi

# Verifica che npm sia disponibile
if ! command -v npm >/dev/null 2>&1; then
  echo "❌ Error: npm is not installed or not in PATH"
  exit 1
fi

cd id_service

# Installa dipendenze solo se non esistono già
if [ ! -d "node_modules" ]; then
  echo "▶ Installing npm dependencies..."
  npm install --no-bin-links
else
  echo "ℹ️ npm dependencies already installed. Skipping."
fi
