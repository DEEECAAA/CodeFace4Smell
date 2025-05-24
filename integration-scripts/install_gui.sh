#!/bin/bash
set -euo pipefail
trap 'echo "Errore alla riga $LINENO"; exit 1' ERR

# Abilita "universe" (solo se Ubuntu)
if grep -qi ubuntu /etc/os-release; then
  sudo add-apt-repository universe -y
fi

# Aggiorna pacchetti
sudo apt update

# Pacchetti GUI minimi + miglioramenti estetici leggeri
sudo apt install -y --no-install-recommends \
  xserver-xorg-core \
  xserver-xorg-input-all \
  xinit \
  openbox \
  lxsession \
  lightdm \
  lightdm-gtk-greeter \
  lxappearance \
  pcmanfm \
  fonts-dejavu \
  tint2 \
  lxterminal \
  build-essential dkms wget linux-headers-$(uname -r)

# Tenta di installare guest additions via pacchetti (se disponibili)
if apt-cache show virtualbox-guest-dkms &>/dev/null; then
  sudo apt install -y virtualbox-guest-dkms virtualbox-guest-x11
else
  echo " Pacchetti VirtualBox Guest Additions non trovati."

  # Ottieni versione VBox (solo se file presente)
  if [ -f /home/vagrant/.vbox_version ]; then
    VBOX_VERSION=$(cat /home/vagrant/.vbox_version)
    wget "https://download.virtualbox.org/virtualbox/${VBOX_VERSION}/VBoxGuestAdditions_${VBOX_VERSION}.iso" -O /tmp/VBoxGuestAdditions.iso
    sudo mkdir -p /mnt/vbox
    sudo mount -o loop /tmp/VBoxGuestAdditions.iso /mnt/vbox
    sudo sh /mnt/vbox/VBoxLinuxAdditions.run || echo " Fallita installazione Guest Additions"
    sudo umount /mnt/vbox
    rm /tmp/VBoxGuestAdditions.iso
  else
    echo " File /home/vagrant/.vbox_version non trovato, impossibile installare Guest Additions"
  fi
fi

# Imposta LightDM come display manager
echo "/usr/sbin/lightdm" | sudo tee /etc/X11/default-display-manager

# File di sessione LXDE
echo "exec startlxsession" | sudo tee /home/vagrant/.xinitrc
sudo chown vagrant:vagrant /home/vagrant/.xinitrc
sudo chmod +x /home/vagrant/.xinitrc

# Imposta Openbox come WM e configura layout semplice
mkdir -p /home/vagrant/.config/openbox
cat <<EOF | sudo tee /home/vagrant/.config/openbox/autostart
tint2 &
lxterminal &
EOF
sudo chown -R vagrant:vagrant /home/vagrant/.config

# Pulizia
sudo apt autoremove -y
