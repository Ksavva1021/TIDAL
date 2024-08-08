#!/bin/bash

export BASE_DIR=$PWD
export LIBRARY_PATH=$LIBRARY_PATH:$PWD/TauAnalysis/ClassicSVfit/lib
make -f TauAnalysis/ClassicSVfit/Makefile -j4
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/TauAnalysis/ClassicSVfit/lib
./TauAnalysis/ClassicSVfit/exec/testClassicSVfit

cd TauAnalysis/ClassicSVfit/wrapper/pybind11
mkdir build
cd build
cmake ..
make check -j4

cd $BASE_DIR

export LIBRARY_PATH=$LIBRARY_PATH:$PWD/TauAnalysis/ClassicSVfit/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/TauAnalysis/ClassicSVfit/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.30.04/x86_64-centosstream9-gcc113-opt/lib/
cmake -S TauAnalysis/ClassicSVfit/wrapper/ -B TauAnalysis/ClassicSVfit/wrapper/ -DBASE_DIR=$BASE_DIR
make -C TauAnalysis/ClassicSVfit/wrapper/
