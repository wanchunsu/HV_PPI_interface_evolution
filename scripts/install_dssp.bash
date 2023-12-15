#!/bin/bash

## Installs DSSP and any dependencies (see https://github.com/PDB-REDO/dssp for more info)

# install libmcfp (required for mrc)
git clone https://github.com/mhekkel/libmcfp.git
cd libmcfp
mkdir build
cd build
cmake ..
cmake --build .
sudo cmake --install . #need to have sudo for permission

# install mrc (required for libcifpp)
git clone https://github.com/mhekkel/mrc.git 
cd mrc
mkdir build
cd build
cmake ..
cmake --build .
sudo cmake --install . #need to have sudo for permission

# First install libsifpp which is required for building DSSP
git clone https://github.com/PDB-REDO/libcifpp.git
cd libcifpp
mkdir build
cd build
cmake ..
cmake --build . --config Release
ctest -C Release
sudo cmake --install . #need to have sudo for permission

# Now install DSSP (this installs the mkdssp program in $HOME/.local/bin)
git clone https://github.com/PDB-REDO/dssp.git
cd dssp
mkdir build
cd build
cmake ..
cmake --build . --config Release
ctest -C Release
sudo cmake --install . #need to have sudo for permission