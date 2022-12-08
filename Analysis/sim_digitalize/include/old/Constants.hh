#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <iostream>
#include <stdint.h>

using namespace std;

const string SW_NAME = "sim_digitalize";

const string SW_VERSION = "v2.0.0";

const string RELEASE_DATE = "2017 Nov 25";

// ################################################################
// all the constants for digitization can be put here
// ################################################################

// parameters that are fixed, no need to tune
const double opt_photon_keV = 9.2;      // Taken from EJ-248 data sheet
const double QE = 0.24;                 // Taken from PMT data sheet
const double light_yield_eff = 0.32;    // Conclusion of Merlin's optical simulation
const double PMT_resolution = 0.6;      // Taken from data measured by CTA coll. using the same PMT

// parameters that can be tuned
const double intrinsic_res = 1.0;
const double xtalk_res = 1.0;

const int ped_frequency = 2000;

// ################################################################

#endif
