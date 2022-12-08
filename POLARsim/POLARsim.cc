#include "G4RunManager.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#endif

#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4VisExecutive.hh"


#include "POLARDetectorConstruction.hh"
#include "POLARPhysicsList.hh"
#include "POLARActionInitialization.hh"
#include "POLARGlobalConfig.hh"

#include "OptionsManager.hh"

int main(int argc,char** argv) {
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

    POLARGlobalConfig* fPOLARGlobalConfig = POLARGlobalConfig::Instance();

    G4UImanager* uiManager = G4UImanager::GetUIpointer();
    G4String execute = "/control/execute ";
    if (options_mgr.config_file != "") {
        uiManager->ApplyCommand(execute + options_mgr.config_file);
        fPOLARGlobalConfig->print_config();
    }

    POLARDetectorConstruction* fPOLARDetectorConstruction = new POLARDetectorConstruction();
    POLARPhysicsList*          fPOLARPhysicsList          = new POLARPhysicsList(fPOLARGlobalConfig->phys_verbose);
    POLARActionInitialization* fPOLARActionInitialization = new POLARActionInitialization(options_mgr.gps_flag, options_mgr.output_file, options_mgr.fixed_name);

#ifdef G4MULTITHREADED
    G4MTRunManager* runManagerMT = NULL;
    G4RunManager*   runManagerSQ = NULL;
    if (options_mgr.mt_flag && options_mgr.gps_flag) {
        runManagerMT = new G4MTRunManager();
        G4cout << "#### using multi-threaded mode with " << options_mgr.num_of_thread << " threads ####" << G4endl;
        runManagerMT->SetNumberOfThreads(options_mgr.num_of_thread);
        // set mandatory initialization classes
        runManagerMT->SetUserInitialization(fPOLARDetectorConstruction);
        runManagerMT->SetUserInitialization(fPOLARPhysicsList);
        runManagerMT->SetUserInitialization(fPOLARActionInitialization);
        // initialize G4 kernel
        G4cout << "initialize G4 kernel" << G4endl;
        runManagerMT->Initialize();
    } else {
        runManagerSQ = new G4RunManager();
        G4cout << "#### using sequential mode ####" << G4endl;
        // set mandatory initialization classes
        runManagerSQ->SetUserInitialization(fPOLARDetectorConstruction);
        runManagerSQ->SetUserInitialization(fPOLARPhysicsList);
        runManagerSQ->SetUserInitialization(fPOLARActionInitialization);
        // initialize G4 kernel
        G4cout << "initialize G4 kernel" << G4endl;
        runManagerSQ->Initialize();
    }
#else
    G4RunManager* runManager = new G4RunManager();
    G4cout << "#### WARNING: only sequential mode can be used ####" << G4endl;
    G4cout << "#### using sequential mode ####" << G4endl;
    // set mandatory initialization classes
    runManager->SetUserInitialization(fPOLARDetectorConstruction);
    runManager->SetUserInitialization(fPOLARPhysicsList);
    runManager->SetUserInitialization(fPOLARActionInitialization);
    // initialize G4 kernel
    G4cout << "initialize G4 kernel" << G4endl;
    runManager->Initialize();
#endif

    G4cout << "#### G4 KERNEL INITIALIZED ####" << G4endl;

    G4VisManager* visManager = NULL;
//     if (options_mgr.gui_flag || options_mgr.ter_flag) {
//         visManager = new G4VisExecutive("warnings"); // "quiet", "errors", "warnings"
//         visManager->Initialize();
//         G4cout << "#### G4 VIS INITIALIZED ####" << G4endl;
//     }

    if (options_mgr.gui_flag) {
//#ifdef G4UI_USE
  G4UIExecutive* ui      = new G4UIExecutive(1, argv);  //define UI session
  G4UImanager* UImanager = G4UImanager::GetUIpointer(); //get the pointer to the UI manager

  G4VisManager* visManager = new G4VisExecutive("warnings");        //initialize visualization
  visManager->Initialize();
  
  //UImanager->ApplyCommand("/control/verbose 2");
  //UImanager->ApplyCommand("/control/saveHistory");
  //UImanager->ApplyCommand("/run/verbose 2");
  //UImanager->ApplyCommand("/run/initialize");

  G4String command = "/control/execute " + options_mgr.vis_mac_file;
  UImanager->ApplyCommand(command);

  ui->SessionStart();

  delete ui;
  delete visManager;
//#else
//        G4cout << " ******** No UI_USE ********" << endl;
//#endif
    } else if (options_mgr.ter_flag) {
        G4UIsession* session = NULL;
#ifdef G4UI_USE_TCSH
        session = new G4UIterminal(new G4UItcsh);
#else
        session = new G4UIterminal();
#endif
        uiManager->ApplyCommand(execute + options_mgr.vis_mac_file);
        session->SessionStart();
        delete session;
    } else {
        // start simulation
        G4cout<< "************** good morning **************" << G4endl;

        if (options_mgr.gps_flag) {
            G4cout << "==== using general particle source ====" << G4endl;
            uiManager->ApplyCommand(execute + options_mgr.primary_file);
        } else {
            G4cout << "==== using root file of particle data ====" << G4endl;
            if (!fPOLARActionInitialization->GetParticleDataFileR()->open(options_mgr.primary_file.c_str())) {
                G4cout << "particle data file open failed." << G4endl;
                return 1;
            }
        }

        uiManager->ApplyCommand(execute + options_mgr.run_mac_file);

        if (!options_mgr.gps_flag) {
            fPOLARActionInitialization->GetParticleDataFileR()->close();
        }

        G4cout<< "************** good night **************" << G4endl;
    }

    // job termination
    if (options_mgr.gui_flag || options_mgr.ter_flag) {
        delete visManager;
    }

#ifdef G4MULTITHREADED
    if (options_mgr.mt_flag) {
        delete runManagerMT;
    } else {
        delete runManagerSQ;
    }
#else
    delete runManager;
#endif
    return 0;
}

