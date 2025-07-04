# Guida di installazione per Codeface4Smell

## Installazione di Codeface
L'installazione di codeface richiede l'utilizzo di vagrant. Come prima cosa si dovrà clonare la repository ed eseguire
il seguento comando:
    
    vagrant up

Al fine di poter avere una GUI Linux funzionante eseguire il seguente comando:

    sudo reboot

### Step Aggiuntivo per Guest Additions
Dopo aver riavviato la macchina installare le guest additions seguendo gli step seguenti:
1. Andare in impostazioni della macchina virtuale
2. Andare nella voce archiviazione ed aggiungere un disco ottico
3. Avviare la VM, cliccare la voce dispositivi e inserire le guest additions
4. Dopo aver installato le guest additions eseguire il "reboot" della macchina


## Setup dell'ambiente
Una volta avviata la vm, da shell eseguire il seguento comando:

    sudo chown -R vagrant:vagrant /home/vagrant/codeface_venv

Il comando consentirà l'attivazione dell'enviroment di codeface. Successivamente eseguire il seguente 
comando per attivare l'env:

    source /home/vagrant/codeface_venv/bin/activate

Una volta eseguito il comando, si dovrà navigare nella cartella contente il file
setup.py (//vagrant), ed eseguire il seguente comando: 

    python setup.py develop

Nel caso di mancato funzionamento del comando eseguire il comande seguente:

    sudo apt upgrade

Tramite le guest additions è possibile attivare le note condivise tra Host e Guest. Tramite questa funzionalità 
spostare la cartella CppStatsFilesUpdate all'interno della cartella 

    /home/vagrant/vendor/

Infine, per visualizzare e testare codeface.conf eseguire il seguente comando:

    codeface  test -c codeface.conf 


## Popolazione del Database
Eseguire il seguente comando dopo aver avviato la vm:

    vagrant ssh

Successivamente eseguire il seguente comando per visualizzare la corretta installazione di MySql:
    
    sudo systemctl status msql

### Se root non ha password
Come primo step fermare eventuali processi errati con i seguenti comandi:

    sudo killall -9 mysql
    sudo rm -f /var/run/mysqld/mysqld.sock
    sudo rm -f /var/run/mysqld/mysqld.pid

Creare la directory del socket (Se serve):

    sudo mkdir -p /var/run/mysqld
    sudo chown mysql:mysql /var/run/mysqld

Avviare MySql senza controlli:
    
    sudo mysqld_safe --skip-grant-tables &

Entrare in MySql ed effettuare il reset della password (aprendo un nuovo terminale):

    sudo mysql

Poi dentro MySQL:
    
    sql
    FLUSH PRIVILEGES;
    ALTER USER 'root'@'localhost' IDENTIFIED BY 'codeface';
    EXIT;

O in alternativa:
    
    ALTER USER 'codeface'@'localhost'
    IDENTIFIED WITH mysql_native_password
    BY 'codeface';

Dopodiché, riavviare MySQL normalmente:
    
    sudo systemctl restart mysql

### Creazione dell'utente codeface e dei due Database
Aprire una nuova bash, ed eseguire i seguenti comandi:

    mysql -u root -p
    # Inserire la password: rootpass123

Oppure, se il comando non funziona eseguire:
    
    sudo mysql -u root -p
    # Inserire la password: rootpass123

Dopodichè creare i due database con i seguenti comandi:
    
    sql
    CREATE DATABASE codeface;
    CREATE DATABASE codeface_testing;
    
    CREATE USER 'codeface'@'localhost' IDENTIFIED BY 'codeface';
    GRANT ALL PRIVILEGES ON codeface.* TO 'codeface'@'localhost';
    GRANT ALL PRIVILEGES ON codeface_testing.* TO 'codeface'@'localhost';
    FLUSH PRIVILEGES;
    EXIT;

Nel caso in cui il comando sql non funzioni, eseguire il comando sottostante e riscrivere
la query riportata precedentemente:

    sudo sql 

### Caricamento dello schema nel database principale
Asicurarsi che il file codeface_schema.sql sia in /vagrant/datamodel o nel path corretto.
Dopodichè aprire la shell ed eseguire i seguenti comandi:

    mysql -ucodeface -pcodeface codeface < /vagrant/datamodel/codeface_schema.sql
    sudo mysql -ucodeface -pcodeface -e "ALTER TABLE cluster ADD COLUMN label TEXT DEFAULT NULL;"

oppure 

    sudo mysql < /vagrant/datamodel/codeface_schema.sql


### Crazione schema nel database di test
Eseguire i seguenti comandi:

    cat /vagrant/datamodel/codeface_schema.sql | sed 's/codeface/codeface_testing/g' | mysql -ucodeface -pcodeface codeface_testing

oppure

    sudo bash -c "sed 's/codeface/codeface_testing/g' /vagrant/datamodel/codeface_schema.sql | mysql codeface_testing"

### Verifica che i database siano correttamente popolati
Aprire la shell ed eseguire i seguenti comandi

    mysql -ucodeface -pcodeface -e "SHOW TABLES;" codeface
    mysql -ucodeface -pcodeface -e "SHOW TABLES;" codeface_testing

## Step Aggiuntivo: cambiare metodo di autenticazione
Entra in MySQL come root. Aprire la shell ed eseguire i seguenti comandi:

    mysql -u root -p
    # Inserisci: rootpass123
oppure

    sudo mysql

### Cambia il metodo di autenticazione dell'utente codeface
Esegui i seguenti comandi:    
    
    sql
    ALTER USER 'codeface'@'localhost'
    IDENTIFIED WITH mysql_native_password
    BY 'codeface';
    FLUSH PRIVILEGES;
    EXIT;

Questo va ad impostare codeface per usare mysql_native_password, compatibile con il client
Python/MySQL usato da CodeFace.

### Comandi per la modifica dei permessi in LOAD DI R SU DB
Attivare nel server:
    
    sudo sed -i '/^\[mysqld\]/a local_infile=1' /etc/mysql/mysql.conf.d/mysqld.cnf
    sudo systemctl restart mysql

Per verificare:

    mysql -u codeface -p -e "SHOW VARIABLES LIKE 'local_infile';"

## Avvio del Server
Aprire una nuova shell, e navigare nella seguente cartella:

    cd /vagrant/id_service

Una volta navigato nella cartella eseguire il seguente comando:

    nodejs id_service.js ../codeface.conf    


### Installazione lxml
Installazione libreria per cppstats:
    
    pip install lxml

## Eseguire codeface
Aprire la shell ed eseguire il seguente comando:

    ps aux | grep node

Creare due directory. La prima directory conterrà i risultati di codeface, mentre la seconda conterrà 
la directory clonata da github. Dopodichè eseguire il seguente comando:
    
    install -d -m 0777 /tua/dir-risultato

Di seguito riportiamo un esempio del file .conf del progetto da analizzare:

    mailinglists: []
        project: Pixel-Arena
        proxyPort: 83
    rcs:
        -null
        -null
        -null
    repo: ../repos/Pixel-Arena
    revisions:
        -71a4adc
        -357bbc4
        -06c0f58
    sleepTime: 1000
    sloccount: false
    tagging: tag
    understand: false

Eseguire il comando:
    
    codeface run -p /vagrant/PixelArena.conf/home/result/home/repos