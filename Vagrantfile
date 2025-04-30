# -*- mode: ruby -*-
# vi: set ft=ruby :

Vagrant.configure("2") do |config|
  config.vm.box = "ubuntu/jammy64" # Ubuntu 22.04 LTS

  config.vm.provider :virtualbox do |vbox|
    vbox.memory = 4096
    vbox.cpus = 2
  end

  # Port forwarding
  config.vm.network "forwarded_port", guest: 8081, host: 8081
  config.vm.network "forwarded_port", guest: 8100, host: 8100

  # Sync scripts dir (opzionale)
  config.vm.synced_folder ".", "/vagrant", type: "virtualbox"

  # Setup provisioning
  config.vm.provision "shell", path: "integration-scripts/install_repositories.sh"
  config.vm.provision "shell", path: "integration-scripts/install_common.sh"
  config.vm.provision "shell", path: "integration-scripts/install_codeface_R.sh"
  config.vm.provision "shell", path: "integration-scripts/install_codeface_node.sh"
  config.vm.provision "shell", path: "integration-scripts/install_codeface_python.sh"
  config.vm.provision "shell", path: "integration-scripts/install_cppstats.sh"
  config.vm.provision "shell", path: "integration-scripts/setup_database.sh"

  # Test execution
  config.vm.provision "shell", path: "integration-scripts/test_codeface.sh"
end
