#include "POLARActionInitialization.hh"
#include "POLARPrimaryGeneratorAction.hh"
#include "POLARRunAction.hh"
#include "POLAREventAction.hh"
#include "POLARSteppingAction.hh"
#include "g4root.hh"
#include "G4RunManager.hh"
#include "POLARGlobalConfig.hh"

POLARActionInitialization::POLARActionInitialization(G4bool the_gps_flag, G4String the_output_file, G4bool fixed_name) {
    fParticleDataFileR_  = new ParticleDataFileR();
    gps_flag_       = the_gps_flag;
    output_file_    = the_output_file;
    fixed_filename_ = fixed_name;
}

POLARActionInitialization::~POLARActionInitialization() {
    delete fParticleDataFileR_;
}

ParticleDataFileR* POLARActionInitialization::GetParticleDataFileR() {
    return fParticleDataFileR_;
}

void POLARActionInitialization::Build() const {
    SetUserAction(new POLARPrimaryGeneratorAction(gps_flag_, fParticleDataFileR_));
    SetUserAction(new POLARRunAction(output_file_, fixed_filename_));
    SetUserAction(new POLAREventAction());
    SetUserAction(new POLARSteppingAction());

    POLARGlobalConfig* fPOLARGlobalConfig = POLARGlobalConfig::Instance();

    // prepare file
    G4AnalysisManager* fAnalysisManager = G4AnalysisManager::Instance();
    if (G4RunManager::GetRunManager()->GetRunManagerType() == G4RunManager::RMType::sequentialRM) {
        fAnalysisManager->CreateNtuple("Keyword", "Keyword");
        fAnalysisManager->CreateNtupleIColumn("RunID");                     // 0
        fAnalysisManager->CreateNtupleIColumn("NumberOfEvents");            // 1
        fAnalysisManager->CreateNtupleIColumn("TotalHits");                 // 2
        fAnalysisManager->CreateNtupleIColumn("TotalEvents");               // 3
        fAnalysisManager->FinishNtuple();
    }

    fAnalysisManager->CreateNtuple("Primary", "Primary");
    // necessary
    fAnalysisManager->CreateNtupleIColumn("EventID");                       // 0
    fAnalysisManager->CreateNtupleIColumn("HitsCount");                     // 1
    // optional
    if (fPOLARGlobalConfig->simout_more) {
        fAnalysisManager->CreateNtupleSColumn("PName_Primary");             // 2
        fAnalysisManager->CreateNtupleFColumn("X_Primary");                 // 3
        fAnalysisManager->CreateNtupleFColumn("Y_Primary");                 // 4
        fAnalysisManager->CreateNtupleFColumn("Z_Primary");                 // 5
        fAnalysisManager->CreateNtupleDColumn("Time_Primary");              // 6
        fAnalysisManager->CreateNtupleFColumn("Energy_Primary");            // 7
        fAnalysisManager->CreateNtupleFColumn("POLX_Primary");              // 8
        fAnalysisManager->CreateNtupleFColumn("POLY_Primary");              // 9
        fAnalysisManager->CreateNtupleFColumn("POLZ_Primary");              // 10
        fAnalysisManager->CreateNtupleFColumn("Theta_Primary");             // 11
        fAnalysisManager->CreateNtupleFColumn("Phi_Primary");               // 12
    }
    fAnalysisManager->FinishNtuple();

    fAnalysisManager->CreateNtuple("HitsCol", "HitsCollection");
    // necessary
    fAnalysisManager->CreateNtupleIColumn("EventID");                       // 0
    fAnalysisManager->CreateNtupleIColumn("TrackID");                       // 1
    fAnalysisManager->CreateNtupleIColumn("ParticleCode");                  // 2
    fAnalysisManager->CreateNtupleIColumn("BarID");                         // 3
    fAnalysisManager->CreateNtupleIColumn("ModID");                         // 4
    fAnalysisManager->CreateNtupleFColumn("LocPosX");                       // 5
    fAnalysisManager->CreateNtupleFColumn("LocPosY");                       // 6
    fAnalysisManager->CreateNtupleFColumn("LocPosZ");                       // 7
    fAnalysisManager->CreateNtupleFColumn("EnergyVis");                     // 8
    fAnalysisManager->CreateNtupleFColumn("EnergyDep");                     // 9
    // optional
    if (fPOLARGlobalConfig->simout_more) {
        fAnalysisManager->CreateNtupleSColumn("ParticleName");              // 10
        fAnalysisManager->CreateNtupleDColumn("PreGlobalTime");             // 11
        fAnalysisManager->CreateNtupleIColumn("StepStatus");                // 12
        fAnalysisManager->CreateNtupleFColumn("PreMomTheta");               // 13
        fAnalysisManager->CreateNtupleFColumn("PreMomPhi");                 // 14
        fAnalysisManager->CreateNtupleFColumn("PostMomTheta");              // 15
        fAnalysisManager->CreateNtupleFColumn("PostMomPhi");                // 16
        fAnalysisManager->CreateNtupleFColumn("PreKinEnergy");              // 17
        fAnalysisManager->CreateNtupleFColumn("PostKinEnergy");             // 18
        fAnalysisManager->CreateNtupleSColumn("ProcessName");               // 19
        fAnalysisManager->CreateNtupleDColumn("DeltaTime");                 // 20
    }
    fAnalysisManager->FinishNtuple();

}

void POLARActionInitialization::BuildForMaster() const {
    SetUserAction(new POLARRunAction(output_file_, fixed_filename_));
    // prepare file
    G4AnalysisManager* fAnalysisManager = G4AnalysisManager::Instance();
    fAnalysisManager->CreateNtuple("Keyword", "Keyword");
    fAnalysisManager->CreateNtupleIColumn("RunID");
    fAnalysisManager->CreateNtupleIColumn("NumberOfEvents");
    fAnalysisManager->CreateNtupleIColumn("TotalHits");
    fAnalysisManager->CreateNtupleIColumn("TotalEvents");
    fAnalysisManager->FinishNtuple();

}
