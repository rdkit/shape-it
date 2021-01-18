#!/bin/bash

cmake \
    -D PYTHON_INSTDIR=$SP_DIR \
    -D BUILD_RDKIT_SUPPORT=ON \
    -D BUILD_PYTHON_SUPPORT=ON \
    -D CMAKE_SYSTEM_PREFIX_PATH=$PREFIX \
    -D CMAKE_INSTALL_PREFIX=$PREFIX \
    -D RDKIT_INCLUDE_DIR=$PREFIX/include/rdkit \
    -D CMAKE_BUILD_TYPE=Release \
    .


make -j$CPU_COUNT install
ctest -j$CPU_COUNT --output-on-failure
