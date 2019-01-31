# ==================================================================
#  tubex-lib - IBEX-lib installation
# ==================================================================

#!/bin/bash
echo 'Installing the IBEX-lib...';

  # CMake Ibex from Benoit Desrochers
  git clone https://github.com/benEnsta/ibex-lib
  cd ibex-lib
  git checkout with_cmake_update_21012019
  mkdir make ; cd make ; cmake .. ; make -j4 ; sudo make install

  # # Master Ibex with waf install
  # 
  # git clone https://github.com/ibex-team/ibex-lib
  # cd ibex-lib
  # ./waf configure --with-debug --interval-lib=filib
  # sudo ./waf install
  # cd ..


##set -x # debugging
#
#if [ ! -e "$HOME/ibex/lib/libibex.a" ]; then
#  echo 'Installing the IBEX-lib...';
#  git clone https://github.com/benEnsta/ibex-lib.git
#  cd ibex-lib
#  git checkout pyIbex_version_3
#  mkdir build
#  cd build
#  cmake -DBUILD_TESTS=OFF -DCMAKE_INSTALL_PREFIX=${HOME}/ibex ../
#  make -j4
#  make install
#fi