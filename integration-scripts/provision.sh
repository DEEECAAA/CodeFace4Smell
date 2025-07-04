#!/bin/bash

# Aggiorna sistema
sudo apt update -y && sudo apt upgrade -y

# Imposta permessi del virtualenv
sudo chown -R vagrant:vagrant /home/vagrant/codeface_venv

# Attiva virtualenv
source /home/vagrant/codeface_venv/bin/activate

# Installa dipendenze Python
pip install lxml || true

# Assicurati che la directory del socket esista
sudo mkdir -p /var/run/mysqld
sudo chown mysql:mysql /var/run/mysqld

# Abilita local_infile in MySQL
sudo sed -i '/^\[mysqld\]/a local_infile=1' /etc/mysql/mysql.conf.d/mysqld.cnf
sudo systemctl restart mysql

# Crea database, utente e assegna privilegi
mysql -uroot -prootpass123 <<EOF
CREATE DATABASE IF NOT EXISTS codeface;
CREATE DATABASE IF NOT EXISTS codeface_testing;
CREATE USER IF NOT EXISTS 'codeface'@'localhost' IDENTIFIED BY 'codeface';
GRANT ALL PRIVILEGES ON codeface.* TO 'codeface'@'localhost';
GRANT ALL PRIVILEGES ON codeface_testing.* TO 'codeface'@'localhost';
FLUSH PRIVILEGES;
EOF

# Carica schema se presente
if [ -f /vagrant/datamodel/codeface_schema.sql ]; then
  mysql -ucodeface -pcodeface codeface < /vagrant/datamodel/codeface_schema.sql
  sed 's/codeface/codeface_testing/g' /vagrant/datamodel/codeface_schema.sql | mysql -ucodeface -pcodeface codeface_testing
fi

# Verifica che local_infile sia attivo
sudo mysql -ucodeface -pcodeface -e "SHOW VARIABLES LIKE 'local_infile';"
sudo mysql -ucodeface -pcodeface -e "ALTER TABLE cluster ADD COLUMN label TEXT DEFAULT NULL;"

echo "Provisioning completato con successo."
