#include <iostream>
#include "EventReader.hh"
#include "OptionsManager.hh"
#include "Config.hh"
#include "Calib.hh"
#include "Digi.hh"
#include "OutDataFile.hh"
#include "OutDataFile_IHEP.hh"
#include "OutDataFile_DPNC.hh"

using namespace std;

int main(int argc, char** argv) {
    // process command parameters
    OptionsManager options_mgr;
    if (!options_mgr.parse(argc, argv)) {
        if (options_mgr.get_version_flag()) {
            options_mgr.print_version();
        } else {
            options_mgr.print_help();
        }
        return 2;
    }

    // configuration
    Config config_obj;
    if (!config_obj.parse(options_mgr.config_filename.Data())) {
        cout << "Configuration file reading failed: " << options_mgr.config_filename.Data() << endl;
        return 1;
    }

    // calibration
    Calib calib_obj;
    if (!calib_obj.read_calib(config_obj)) {
        cout << "Calibration file reading failed: " << options_mgr.config_filename.Data() << endl;
        return 1;
    }
    calib_obj.generate_ped_gain(config_obj);
    calib_obj.uncouple_xtalk_gain();

    // simulation data reading
    EventReader eventReader;
    if (!eventReader.open(options_mgr.simdata_filename.Data())) {
        cout << "simulation data open failed: " << options_mgr.simdata_filename.Data() << endl;
        return 1;
    }
    eventReader.primary_deactive_all();
    eventReader.hitscol_deactive_all();
    eventReader.hitscol_active("EventID");
    eventReader.hitscol_active("BarID");
    eventReader.hitscol_active("ModID");
    eventReader.hitscol_active("LocPosZ");
    eventReader.hitscol_active("EnergyVis");
    eventReader.hitscol_active("EnergyDep");

    // output file opening
    OutDataFile* datafile = NULL;
    if (options_mgr.output_type == "IHEP") {
        datafile = new OutDataFile_IHEP();
    } else if (options_mgr.output_type == "DPNC") {
        datafile = new OutDataFile_DPNC();
    } else {
        cout << "output type must be IHEP or DPNC." << endl;
        datafile = NULL;
        return 1;
    }
    if (!datafile->open(options_mgr.digiout_filename.Data())) {
        cout << "output file open failed: " << options_mgr.digiout_filename.Data() << endl;
        return 1;
    }

    // start digitization
    Digi digi_obj;
    digi_obj.set_seed(options_mgr.rand_seed);

    int pre_percent = 0;
    int cur_percent = 0;
    cout << "digitalizing ..." << endl;
    cout << "[ " << flush;
    eventReader.event_start();
    while (eventReader.event_next()) {
        cur_percent = static_cast<int>(eventReader.hitscol_get_cur_entry() * 100.0 / eventReader.hitscol_get_tot_entries());
        if (cur_percent - pre_percent > 0 && cur_percent % 2 == 0) {
            pre_percent = cur_percent;
            cout << "#" << flush;
        }
        digi_obj.clear_event();
        digi_obj.energy_response(eventReader.HitsCol, calib_obj);
        if (digi_obj.event_selection(config_obj)) {
            datafile->write_event(digi_obj.general_event);
            if (datafile->get_total_events() % ped_frequency == 0) {
                digi_obj.clear_event();
                digi_obj.generate_ped_event(calib_obj, config_obj);
                datafile->write_event(digi_obj.general_event);
            }
        }
    }
    cout << " DONE ]" << endl;
    // finish digitization

    datafile->write_tree();
    datafile->close();
    delete datafile;
    datafile = NULL;

    return 0;
}
