#ifndef Config_HH
#define Config_HH

#include <iostream>

using namespace std;

class Config {
public:
    // config parameters
    string ped_temp_filename;
    string noise_filename;
    string threshold_filename;
    string threshold_corr_filename;
    string crosstalk_filename;
    string gain_hv_filename;
    string nonlin_filename;

    int   pre_c;
    int   pre_s;
    float too_many_mean;
    float too_many_width;

    float temp[25];
    float hv[25];

public:
    Config();
    ~Config();

    void clear();
    bool parse(const char* config_fn);
    void print();

};

#endif
