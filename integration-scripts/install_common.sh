#!/bin/sh
set -e

echo "Installing common system dependencies"

# Update e installazione pacchetti
sudo apt-get update -qq
sudo DEBIAN_FRONTEND=noninteractive apt-get install -y \
  mysql-server mysql-client default-jdk \
  texlive graphviz doxygen exuberant-ctags git subversion sloccount \
  libxml2-dev libcurl4-openssl-dev libmysqlclient-dev libcairo2-dev \
  libxt-dev astyle xsltproc build-essential libyaml-dev gfortran \
  python3 python3-dev python3-pip python3-setuptools python3-numpy python3-matplotlib \
  python3-lxml gcc libpoppler-dev libpoppler-glib-dev libx11-dev libglu1-mesa-dev \
  libgles2-mesa-dev xorg-dev screen

# Pulizia pacchetti duplicati o obsoleti rimossi

echo "Installing universal-ctags from source"
sudo apt-get install -y autoconf automake pkg-config libseccomp-dev libjansson-dev \
                        libyaml-dev libxml2-dev

echo "Installing dbus dependencies"
sudo apt-get install -y libdbus-1-dev libdbus-glib-1-dev
pip install dbus-python

git clone https://github.com/universal-ctags/ctags.git /tmp/ctags
cd /tmp/ctags
./autogen.sh
./configure
make
sudo make install
cd -
rm -rf /tmp/ctags
