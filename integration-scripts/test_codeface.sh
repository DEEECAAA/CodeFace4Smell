#!/bin/bash
set -e

echo "▶ Running CodeFace system test"

cd /vagrant

# Attiva l'ambiente virtuale
if [ -f /home/vagrant/codeface_venv/bin/activate ]; then
  source /home/vagrant/codeface_venv/bin/activate
else
  echo "❌ Virtualenv not found at /home/vagrant/codeface_venv"
  exit 1
fi

# Verifica presenza cartella id_service
if [ ! -d "id_service" ]; then
  echo "❌ Error: 'id_service' directory not found!"
  exit 1
fi

# Verifica dipendenze
if ! command -v node >/dev/null 2>&1; then
  echo "❌ Error: Node.js is not installed"
  exit 1
fi

if ! command -v codeface >/dev/null 2>&1; then
  echo "❌ Error: 'codeface' command not available in virtualenv"
  exit 1
fi

# Avvia id_service
cd id_service
node id_service.js ../codeface.conf &
node_pid=$!
cd ..

sleep 2

# Esegui test
echo "▶ Running codeface test..."
if codeface test -c codeface.conf; then
  echo "✅ Codeface test passed"
else
  echo "❌ Codeface test failed"
fi

# Termina node
echo "▶ Shutting down id_service..."
kill "$node_pid" || echo "⚠️ Node process not running"
