#!/bin/sh
set -e

echo "Setting up Python environment for CodeFace"

# Evita /vagrant: crea venv altrove (es. /home/vagrant/codeface_venv)
VENV_PATH="/home/vagrant/codeface_venv"
cd /vagrant

# Ensure virtualenv is installed
if ! command -v virtualenv >/dev/null 2>&1; then
  echo "virtualenv not found, installing..."
  pip install virtualenv
fi

# Create and activate virtualenv
if [ ! -d "$VENV_PATH" ]; then
  virtualenv -p python3 "$VENV_PATH"
fi

. "$VENV_PATH/bin/activate"

# Install dependencies
pip install --upgrade setuptools mock

# Development install
if [ -f "setup.py" ]; then
  python setup.py develop
else
  echo "Error: setup.py not found"
  exit 1
fi
