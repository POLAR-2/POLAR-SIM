#!/bin/sh

geant4_base=/hxmt/work/POLAR/PSDC/software/Geant4/Geant4.10.02.p02

Geant4_DIR=$geant4_base/lib64/Geant4-10.2.2
CMAKE_MODULE_PATH=$Geant4_DIR/Modules

cmake $1 -DGeant4_DIR=$Geant4_DIR -DCMAKE_MODULE_PATH=$CMAKE_MODULE_PATH
