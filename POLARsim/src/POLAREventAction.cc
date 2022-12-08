#include "POLAREventAction.hh"
#include "POLARGlobalConfig.hh"

POLAREventAction::POLAREventAction() {

}

POLAREventAction::~POLAREventAction() {

}

void POLAREventAction::BeginOfEventAction(const G4Event* anEvent) {
    POLARGlobalConfig* fPOLARGlobalConfig = POLARGlobalConfig::Instance();
    if (fPOLARGlobalConfig->event_verbose < 1) return;
    G4int nEvent = anEvent -> GetEventID();
    if((nEvent % 10000000 == 0) && (nEvent != 0)) {
        G4cout << "INFORMATION: " << nEvent << " event in progress..." << G4endl;
    }
}

void POLAREventAction::EndOfEventAction(const G4Event* anEvent) {
    POLARGlobalConfig* fPOLARGlobalConfig = POLARGlobalConfig::Instance();
    if (fPOLARGlobalConfig->event_verbose < 1) return;
    G4int nmbEvents = (anEvent -> GetEventID()) + 1;
    if(nmbEvents % 10000 == 0)
    {
        G4cout << "INFORMATION: " << nmbEvents << " events processed." << G4endl;
    }
}
