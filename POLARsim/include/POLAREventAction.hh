#ifndef POLAREventAction_hh
#define POLAREventAction_hh 1

#include "G4UserEventAction.hh"
#include "G4Event.hh"

class POLAREventAction: public G4UserEventAction {
public:
    POLAREventAction();
    ~POLAREventAction();

    void BeginOfEventAction(const G4Event* anEvent);
    void EndOfEventAction(const G4Event* anEvent);

};


#endif
