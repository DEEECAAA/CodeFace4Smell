#!/bin/sh
set -e

CPPSTATS_VERSION=0.8.4
INSTALL_DIR="$PWD/vendor"

echo "Installing cppstats $CPPSTATS_VERSION"

mkdir -p "$INSTALL_DIR"
cd "$INSTALL_DIR"

# Download cppstats
wget --quiet https://codeload.github.com/clhunsen/cppstats/tar.gz/v$CPPSTATS_VERSION -O cppstats.tar.gz
tar -xf cppstats.tar.gz
CPPSTATS="$INSTALL_DIR/cppstats-$CPPSTATS_VERSION"

# Create wrapper script
cat <<EOF > "$CPPSTATS/cppstats"
#!/bin/bash
cd "$CPPSTATS"
PYTHONPATH="\$PYTHONPATH:$CPPSTATS/lib" ./cppstats.py "\$@"
EOF

chmod +x "$CPPSTATS/cppstats"
sudo ln -sf "$CPPSTATS/cppstats" /usr/local/bin/cppstats

# Download latest srcML manually
echo ">> You should download a recent srcML manually from: https://www.srcml.org/download/"
echo ">> And extract it into: $CPPSTATS/lib/srcml/linux/"
echo ">> This script does NOT fetch old Ubuntu 12.04 binaries anymore."

cd ..
