# -*- mode: ruby -*-
# vi: set ft=ruby :

Vagrant.configure("2") do |config|
  # Box base Ubuntu 22.04 LTS
  config.vm.box = "ubuntu/jammy64"

  # Impostazioni di VirtualBox
  config.vm.provider :virtualbox do |vbox|
    vbox.memory = 4096
    vbox.cpus = 4
    vbox.gui = true # Attiva la GUI della VM
  end

  # Porte inoltrate (facoltative)
  config.vm.network "forwarded_port", guest: 8081, host: 8081
  config.vm.network "forwarded_port", guest: 8100, host: 8100

  # Cartella sincronizzata
  config.vm.synced_folder ".", "/vagrant", type: "virtualbox"

  # Provisioning script
  config.vm.provision "shell", path: "integration-scripts/install_repositories.sh", run: "always"
  config.vm.provision "shell", path: "integration-scripts/install_common.sh", run: "always"
  config.vm.provision "shell", path: "integration-scripts/install_codeface_R.sh", run: "always"
  config.vm.provision "shell", path: "integration-scripts/install_codeface_node.sh", run: "always"
  config.vm.provision "shell", path: "integration-scripts/install_codeface_python.sh", run: "always"
  config.vm.provision "shell", path: "integration-scripts/install_cppstats.sh", run: "always"
  config.vm.provision "shell", path: "integration-scripts/setup_database.sh", run: "always"

  # Installazione GUI
  config.vm.provision "shell", path: "integration-scripts/install_gui.sh", run: "always"

  # Test finale
  config.vm.provision "shell", path: "integration-scripts/test_codeface.sh", run: "always"
end
