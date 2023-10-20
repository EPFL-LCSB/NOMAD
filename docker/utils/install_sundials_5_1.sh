#!/bin/bash

wget https://github.com/LLNL/sundials/releases/download/v5.1.0/sundials-5.1.0.tar.gz
tar -xzf sundials-5.1.0.tar.gz -C $HOME
cd $HOME/sundials-5.1.0

mkdir $HOME/sundials-5.1.0/builddir
cd $HOME/sundials-5.1.0/builddir
cmake -DLAPACK_ENABLE=ON \
      -DSUNDIALS_INDEX_SIZE=64 \
      -DCMAKE_INSTALL_PREFIX=$HOME/sundials-5.1.0 \
      $HOME/sundials-5.1.0

make
make install


ls $HOME/sundials-5.1.0/include
