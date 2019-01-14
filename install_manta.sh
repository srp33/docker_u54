#! /bin/bash

wget https://github.com/Illumina/manta/releases/download/v${MANTA_VERSION}/\
manta-${MANTA_VERSION}.release_src.tar.bz2
tar -xjf manta-${MANTA_VERSION}.release_src.tar.bz2
mkdir build && cd build

../manta-${MANTA_VERSION}.release_src/configure --prefix=/path/to/install
make -j4 install