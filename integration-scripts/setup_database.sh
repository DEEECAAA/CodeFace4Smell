#!/bin/sh
set -e

echo "Setting up CodeFace MySQL database"

cd /vagrant  # 🔍 Assicura che i percorsi relativi funzionino

# Parametri
DB_USER="codeface"
DB_PASS="codeface"
DB_NAME="codeface"
DB_TEST="codeface_testing"
SCHEMA_FILE="/vagrant/datamodel/codeface_schema.sql"
MYSQL_ROOT_PASS="root"

# Crea utente solo se non esiste
if ! sudo mysql -u root -p"$MYSQL_ROOT_PASS" -e "SELECT user FROM mysql.user WHERE user = '$DB_USER';" | grep -q "$DB_USER"; then
    echo "Creating user $DB_USER..."
    sudo mysql -u root -p"$MYSQL_ROOT_PASS" -e "CREATE USER '$DB_USER'@'localhost' IDENTIFIED BY '$DB_PASS';"
else
    echo "User $DB_USER already exists, skipping."
fi

# Crea i database se non esistono
for DB in "$DB_NAME" "$DB_TEST"; do
  if ! sudo mysql -u root -p"$MYSQL_ROOT_PASS" -e "SHOW DATABASES LIKE '$DB';" | grep -q "$DB"; then
    echo "Creating database $DB..."
    sudo mysql -u root -p"$MYSQL_ROOT_PASS" -e "CREATE DATABASE $DB;"
  else
    echo "Database $DB already exists, skipping."
  fi
done

# Concedi privilegi
sudo mysql -u root -p"$MYSQL_ROOT_PASS" -e "GRANT ALL PRIVILEGES ON $DB_NAME.* TO '$DB_USER'@'localhost';"
sudo mysql -u root -p"$MYSQL_ROOT_PASS" -e "GRANT ALL PRIVILEGES ON $DB_TEST.* TO '$DB_USER'@'localhost';"
sudo mysql -u root -p"$MYSQL_ROOT_PASS" -e "FLUSH PRIVILEGES;"

# Carica schema solo se tabella principale non esiste
if [ -f "$SCHEMA_FILE" ]; then
  if ! mysql -u"$DB_USER" -p"$DB_PASS" "$DB_NAME" -e "SHOW TABLES LIKE 'commit';" | grep -q "commit"; then
    echo "Importing schema into $DB_NAME and $DB_TEST..."
    mysql -u"$DB_USER" -p"$DB_PASS" "$DB_NAME" < "$SCHEMA_FILE"
    sed "s/$DB_NAME/$DB_TEST/g" "$SCHEMA_FILE" | mysql -u"$DB_USER" -p"$DB_PASS" "$DB_TEST"
  else
    echo "Tables already exist, skipping schema import."
  fi
else
  echo "Schema file $SCHEMA_FILE not found!"
  exit 1
fi
