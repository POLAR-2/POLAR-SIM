# POLAR-SIM
## Introduction

This is the final merged simulation software package of POLAR project. It will be maintained by POLAR simulation group. 

## How to compile and run the G4 simulation

* cd POLARsim
* mkdir build
* cd build
* ../g4cmake.sh ..  # change the value geant4_base in this script to the right one your computer
* make
* ./POLARsim gps.mac run.mac -c config.mac -o output.root
* the simulation data file is then stored in directore G4out

## How to compile and run the analysis part

* cd Analysis
* make
* ./bin/sim_digitalize -c config.cfg sim_data.root -t IHEP -o digi_data.root  # for IHEP output data structure
* ./bin/sim_digitalize -c config.cfg sim_data.root -t DPNC -o digi_data.root  # for DPNC output data structure

## Note for using Geant4-10.2.0

Change line 66 in POLARSensitiveDetector.cc into:
    "G4double          energyVis         = fEmSaturation_->VisibleEnergyDeposition(aStep);"
