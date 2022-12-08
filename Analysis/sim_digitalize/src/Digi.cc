#include "Digi.hh"
#include "Constants.hh"
#include <iomanip>
#include <cmath>

Digi::Digi() {
    set_seed(0);
    clear_event();

    pre_s_count_ = 0;
    pre_c_count_ = 0;

    for (int i = 0; i < 25; i++) {
        raw_rate_counter_[i] = 0;
    }
}

Digi::~Digi() {

}

void Digi::set_seed(int seed) {
    digi_random_.SetSeed(seed);
}

void Digi::clear_event() {
    general_event.event_id = 0;
    general_event.is_ped = false;
    general_event.type = 0;
    general_event.trigger_n = 0;
    general_event.pkt_count = 0;
    for (int i = 0; i < 25; i++) {
        general_event.trig_accepted[i] = false;
        general_event.multiplicity[i] = 0;
        general_event.compress[i] = -1;
        general_event.raw_rate[i] = -1;
        general_event.t_out_too_many[i] = false;
        general_event.t_out_1[i] = false;
        general_event.t_out_2[i] = false;
        general_event.fe_hv[i] = 0;
        general_event.fe_temp[i] = 100;
        for (int j = 0; j < 64; j++) {
            general_event.trigger_bit[i][j] = false;
            general_event.energy_adc[i][j] = -1;
            general_event.raw_energy[i][j] = 0;
        }
    }
}

void Digi::energy_response(const vector<EventReader::HitsCol_T>& hitscol_vec, const Calib& calib_obj, const Config& config_obj) {
    // energy -> ADC

    if (hitscol_vec.size() < 1) return;

    general_event.event_id = hitscol_vec[0].EventID;

// ==============================================================================
// ==============================================================================

    // (1) energy to optical photon
    double bar_optical_photon[25][64];
    for (int i = 0; i < 25; i++) {
        for (int j = 0; j < 64; j++) {
            general_event.raw_energy[i][j] = 0;
            bar_optical_photon[i][j] = 0;
        }
    }
    for (size_t i = 0; i < hitscol_vec.size(); i++) {
        int mod_id = hitscol_vec[i].ModID;
        int bar_id = hitscol_vec[i].BarID;
        double energy = hitscol_vec[i].EnergyVis;
        // double energy = hitscol_vec[i].EnergyDep;
        double locpos = hitscol_vec[i].LocPosZ;

        // according to the conclusion of Merlin's optical photon simulation
        double cur_eff = light_yield_eff;
        if (locpos < 4) {
            cur_eff = 0.6 - locpos * 0.28 / 4.0;
        } else if (locpos > 172) {
            cur_eff=0.3 + (locpos - 172) * 0.28 / 4.0;
        }

        general_event.raw_energy[mod_id][bar_id] += energy;

        double mean_photon = opt_photon_keV * energy;
        double rand_photon = digi_random_.Gaus(mean_photon, config_obj.intrinsic_res * TMath::Sqrt(mean_photon));
        if (rand_photon > 0) {
            bar_optical_photon[mod_id][bar_id] += round(rand_photon * cur_eff);
        }
    }

// ==============================================================================

    // (2) apply crosstalk
    double bar_optical_photon_ax[25][64];

    // # method 1: keep the total photons of one bar unchanged
    // for (int i = 0; i < 25; i++) {
    //     for (int j = 0; j < 64; j++) {
    //         bar_optical_photon_ax[i][j] = bar_optical_photon[i][j];
    //     }
    // }
    // for (int i = 0; i < 25; i++) {
    //     for (int jx = 0; jx < 64; jx++) {
    //         for (int jy = 0; jy < 64; jy++) {
    //             if (jy == jx) continue;
    //             if (bar_optical_photon[i][jy] > 0) {
    //                 double tmp_value = calib_obj.photon_xtalk[i](jx, jy) * bar_optical_photon[i][jy];
    //                 // double neigh_opt = digi_random_.PoissonD(tmp_value);
    //                 double neigh_opt = digi_random_.Gaus(tmp_value, TMath::Sqrt(tmp_value));
    //                 neigh_opt = round(neigh_opt);
    //                 if (neigh_opt < 0) neigh_opt = 0;
    //                 bar_optical_photon_ax[i][jx] += neigh_opt;
    //                 bar_optical_photon_ax[i][jy] -= neigh_opt;
    //             }
    //         }
    //     }
    // }

    // # method 2: randomize all
    for (int i = 0; i < 25; i++) {
        for (int j = 0; j < 64; j++) {
            bar_optical_photon_ax[i][j] = 0;
        }
    }
    for (int i = 0; i < 25; i++) {
        for (int jx = 0; jx < 64; jx++) {
            for (int jy = 0; jy < 64; jy++) {
                if (bar_optical_photon[i][jy] > 0) {
                    double tmp_value = calib_obj.photon_xtalk[i](jx, jy) * bar_optical_photon[i][jy];
                    // bar_optical_photon_ax[i][jx] += digi_random_.PoissonD(tmp_value);
                    double neigh_opt = digi_random_.Gaus(tmp_value, config_obj.xtalk_res * TMath::Sqrt(tmp_value));
                    neigh_opt = round(neigh_opt);
                    if (neigh_opt < 0) neigh_opt = 0;
                    bar_optical_photon_ax[i][jx] += neigh_opt;
                }
            }
        }
    }

// ==============================================================================

    // (3) opt_photon to NPe
    double bar_NPe[25][64];
    for (int i = 0; i < 25; i++) {
        for (int j = 0; j < 64; j++) {
            if (bar_optical_photon_ax[i][j] > 0) {
                double tmp_value = QE * bar_optical_photon_ax[i][j];
                bar_NPe[i][j] = digi_random_.PoissonD(tmp_value);
                // bar_NPe[i][j] = round(digi_random_.Gaus(tmp_value, TMath::Sqrt(tmp_value)));
                // if (bar_NPe[i][j] < 0) bar_NPe[i][j] = 0;
            } else {
                bar_NPe[i][j] = 0;
            }
        }
    }

// ==============================================================================

    // (4) NPe to ADC
    double bar_ADC[25][64];
    for (int i = 0; i < 25; i++) {
        for (int j = 0; j < 64; j++) {
            if (bar_NPe[i][j] > 0) {

                double rand_NPe = digi_random_.Gaus(bar_NPe[i][j], PMT_resolution * TMath::Sqrt(bar_NPe[i][j]));
                bar_ADC[i][j] = rand_NPe * calib_obj.gain_true[i](j) / (opt_photon_keV * QE * light_yield_eff);
                if (bar_ADC[i][j] < 0) bar_ADC[i][j] = 0;

            } else {
                bar_ADC[i][j] = 0;
            }
        }
    }

// ==============================================================================

    // (5) add intrinsic_noise, then give it to general_event for later trigger decision
    for (int i = 0; i < 25; i++) {
        for (int j = 0; j < 64; j++) {
            double intrinsic_noise = digi_random_.Gaus(0, calib_obj.intrinsic_noise[i](j));
            general_event.energy_adc[i][j] = bar_ADC[i][j] + intrinsic_noise;
        }
    }

// ==============================================================================

    // (6) decide the trigger_bit according to trigger efficiency
    for (int i = 0; i < 25; i++) {
        for (int j = 0; j < 64; j++) {
            double cur_vthr = digi_random_.Gaus(calib_obj.vthr_mean[i](j), calib_obj.vthr_width[i](j));
            // double cur_vthr = calib_obj.vthr_mean[i](j);
            if (general_event.energy_adc[i][j] > cur_vthr) {
                general_event.trigger_bit[i][j] = true;
                general_event.multiplicity[i]++;
                general_event.trigger_n++;
            }
        }
    }

// ==============================================================================
    // (7) apply the non-linearity based on the max ADC
    for (int i = 0; i < 25; i++) {
        // only when DAC1 triggered
        if (general_event.multiplicity[i] < 1) continue;
        // find the max ADC of this module
        double max_ADC = 0;
        for (int j = 0; j < 64; j++) {
            if (max_ADC < general_event.energy_adc[i][j]) {
                max_ADC = general_event.energy_adc[i][j];
            }
        }

        // calculate and apply the non-linearity
        for (int j = 0; j < 64; j++) {
            // nonlinearity parameters
            double nonlin_p0 = calib_obj.nonlin_p0[i](j);
            double nonlin_p1 = calib_obj.nonlin_p1[i](j);
            double nonlin_p2 = calib_obj.nonlin_p2[i](j);
            // convert the true max ADC to the measured one
            double maxADC_meas = get_measured_maxADC_(max_ADC, nonlin_p0, nonlin_p1, nonlin_p2);
            double ratio = nonlin_p0 * (1.0 + nonlin_p1 * maxADC_meas) * (1 + TMath::Erf(maxADC_meas / nonlin_p2)) / 2;
            // double ratio = 1.0;
            general_event.energy_adc[i][j] *= ratio;
        }
    }

// ==============================================================================

    // (8) add common noise and pedestal, then check overflow
    for (int i = 0; i < 25; i++) {
        double common_noise = digi_random_.Gaus(0, calib_obj.common_noise(i));
        for (int j = 0; j < 64; j++) {
            general_event.energy_adc[i][j] += (common_noise + calib_obj.pedestal_mean[i](j));
            general_event.energy_adc[i][j] = round(general_event.energy_adc[i][j]); // raw ADC should be an integer
            if (general_event.energy_adc[i][j] > 4095) general_event.energy_adc[i][j] = 4095;
        }
    }

// ==============================================================================
// ==============================================================================

    // (9) decide trig_accepted
    for (int i = 0; i < 25; i++) {
        if (general_event.multiplicity[i] > 0) {
            general_event.trig_accepted[i] = true;
            general_event.compress[i] = 0;
            general_event.pkt_count++;
        }
    }

// ==============================================================================

    // (10) this event is not a pedestal event
    general_event.is_ped = false;

// ==============================================================================

//    // debug start
//    if (general_event.trigger_n < 1) return;
//    cout << "----------------------------------------------------------------------" << endl;
//    cout << "EventID: " << general_event.event_id << endl;
//    for (int i = 0; i < 25; i++) {
//        if (general_event.trig_accepted[i]) {
//            cout << Form("CT_%02d", i + 1) << endl;
//            for (int j = 0; j < 64; j++) {
//                cout << Form("bar_%02d", j)
//                     << Form("%8.2f", bar_energy_vis[i][j]                                           )
//                     << Form("%8.2f", bar_optical_photon[i][j]                                       )
//                     << Form("%8.2f", bar_optical_photon_ax[i][j]                                    )
//                     << Form("%8.2f", bar_NPe[i][j]                                                  )
//                     << Form("%8.2f", bar_ADC[i][j]                                                  )
//                     << Form("%8.2f", bar_ADC_nonlinear[i][j]                                        )
//                     << Form("%8.2f", bar_ADC_noise[i][j]                                            )
//                     << Form("%8.2f", general_event.energy_adc[i][j]                                 )
//                     << Form("%3d", general_event.trigger_bit[i][j]                                  )
//                     << Form("%8.2f", calib_obj.vthr_mean[i][j] + calib_obj.vthr_mean_corr[i][j]     )
//                     << Form("%8.2f", calib_obj.vthr_width[i][j] + calib_obj.vthr_width_corr[i][j]   )
//                     << endl;
//            }
//        }
//    }
//    // debug stop

}

bool Digi::event_selection(const Config& config_obj) {
    // decide if this event is acceptable

    bool is_accepted = true;

    general_event.type = 0x00FF;

    if (general_event.trigger_n < 1) {
        is_accepted = false;
    } else if (general_event.trigger_n == 1) {
        // single event
        general_event.type = 0xF000;
        pre_s_count_++;
        if (pre_s_count_ % config_obj.pre_s > 0) {
            is_accepted = false;
        }
    }

    bool is_cosmic = false;
    for (int i = 0; i < 25; i++) {
        if (general_event.multiplicity[i] > 0) general_event.t_out_1[i] = true;
        if (general_event.multiplicity[i] > 1) general_event.t_out_2[i] = true;
        double too_many_cut = digi_random_.Gaus(config_obj.too_many_mean, config_obj.too_many_width);
        if (general_event.multiplicity[i] > too_many_cut) {
            general_event.t_out_too_many[i] = true;
            is_cosmic = true;
        }
    }

    if (is_cosmic) {
        general_event.type = 0xFF00;
        pre_c_count_++;
        if (pre_c_count_ % config_obj.pre_c > 0) {
            is_accepted = false;
        }
    }

    // copy temperature and HV
    for (int i = 0; i < 25; i++) {
        general_event.fe_temp[i] = config_obj.temp[i];
        general_event.fe_hv[i] = config_obj.hv[i];
    }

    if (is_accepted) {
        for (int i = 0; i < 25; i++) {
            if (general_event.trig_accepted[i]) {
                general_event.raw_rate[i] = raw_rate_counter_[i];
                raw_rate_counter_[i] = 0;
            }
        }
        return true; // record this event
    } else {
        for (int i = 0; i < 25; i++) {
            if (general_event.trig_accepted[i]) {
                raw_rate_counter_[i] += 1;
            }
        }
        return false; // throw out this event
    }
}

void Digi::generate_ped_event(const Calib& calib_obj, const Config& config_obj) {
    // generate pedestal event, actually this is not needed
    general_event.event_id = -1;
    general_event.is_ped = true;
    general_event.pkt_count = 25;
    general_event.type = 0x00F0;
    for (int i = 0; i < 25; i++) {
        general_event.trig_accepted[i] = true;
        general_event.compress[i] = 2;
        general_event.raw_rate[i] = raw_rate_counter_[i];
        raw_rate_counter_[i] = 0;
        double common_noise = digi_random_.Gaus(0, calib_obj.common_noise(i));
        for (int j = 0; j < 64; j++) {
            double intrinsic_noise = digi_random_.Gaus(0, calib_obj.intrinsic_noise[i](j));
            general_event.energy_adc[i][j] = round(calib_obj.pedestal_mean[i](j) + common_noise + intrinsic_noise);
        }
        general_event.fe_hv[i] = config_obj.hv[i];
        general_event.fe_temp[i] = config_obj.temp[i];
    }
}

double Digi::get_measured_maxADC_(const double maxADC, const double p0, const double p1, const double p2) {
    // solve the equation: maxADC_true = maxADC_meas / (p0 * (1 + p1 * maxADC_meas) * (1 + erf(maxADC_meas / p2)) / 2)
    // to get the maxADC_meas from maxADC_true

    const double meas_max = 5000;
    const double precision = 1.0E-3;
    const double step_l = 100;

    int step_n = meas_max / step_l;

    if (maxADC < 1) return 0;
    if (maxADC > meas_max) return meas_max;

    double head(0), tail(0), center(0);
    bool found = false;
    for (int i = 0; i < step_n; i++) {
        head = i * step_l;
        tail = i * step_l + step_l;
        double head_adc = head / (p0 * (1 + p1 * head) * (1 + TMath::Erf(head / p2)) / 2);
        double tail_adc = tail / (p0 * (1 + p1 * tail) * (1 + TMath::Erf(tail / p2)) / 2);
        if (abs(head_adc - maxADC) < precision) return head;
        if (abs(tail_adc - maxADC) < precision) return tail;
        if (maxADC > head_adc && maxADC < tail_adc) {
            found = true;
            break;
        }
    }
    if (!found) return meas_max;

    while(true) {
        center = (head + tail) / 2;
        double center_adc = center / (p0 * (1 + p1 * center) * (1 + TMath::Erf(center / p2)) / 2);
        if (abs(center_adc - maxADC) < precision) return center;
        if (center_adc < maxADC) {
            head = center;
        } else {
            tail = center;
        }
    }

}
