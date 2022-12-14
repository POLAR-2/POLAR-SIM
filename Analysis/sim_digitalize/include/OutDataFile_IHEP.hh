#ifndef OutDataFile_IHEP_HH
#define OutDataFile_IHEP_HH

#include "OutDataFile.hh"

class OutDataFile_IHEP : public OutDataFile {
public:
    struct SAPType {
        // from SCI
        Long64_t        event_id;
        Double_t        event_time;
        Int_t           type;
        Bool_t          is_ped;
        Int_t           packet_num;
        Int_t           ct_pre_is_bad;
        Int_t           fe_pre_is_bad[25];
        Double_t        ct_time_second;
        Double_t        ct_time_wait;
        Float_t         ct_dead_ratio;
        Double_t        fe_time_second[25];
        Double_t        fe_time_wait[25];
        Float_t         fe_dead_ratio[25];
        Bool_t          trig_accepted[25];
        Bool_t          time_aligned[25];
        Int_t           raw_rate[25];
        Int_t           pkt_count;
        Int_t           lost_count;
        Bool_t          trigger_bit[25][64];
        Int_t           trigger_n;
        Int_t           multiplicity[25];
        Float_t         energy_value[25][64];
        UShort_t        channel_status[25][64];
        Float_t         common_noise[25];
        Int_t           compress[25];
        Bool_t          dy12_too_high[25];
        Bool_t          t_out_too_many[25];
        Bool_t          t_out_2[25];
        Bool_t          t_out_1[25];
        // from AUX
        Float_t         aux_interval;
        Int_t           obox_mode;
        Float_t         fe_hv[25];
        Float_t         fe_thr[25];
        Float_t         fe_temp[25];
        // from PPD
        Float_t         ppd_interval;
        Double_t        wgs84_xyz[3];
        Double_t        det_z_radec[2];
        Double_t        det_x_radec[2];
        Double_t        earth_radec[2];
        Double_t        sun_radec[2];
    };

private:
    TFile* t_file_out_;
    TTree* t_pol_event_tree_;
    TTree* t_raw_energy_tree_;

private:
    SAPType t_pol_event_;
    Float_t t_raw_energy_[25][64];

private:
    long int total_events_;

private:
    void clear_data_();

public:
    OutDataFile_IHEP();
    virtual ~OutDataFile_IHEP();

    bool open(const char* filename);
    void write_event(const Digi::GeneralEvent_T& general_event);
    void write_tree();
    void close();
    long int get_total_events();

};

#endif
