#ifndef POLARSteppingAction_h
#define POLARSteppingAction_h 1

#include "G4UserSteppingAction.hh"

class POLARSteppingAction: public G4UserSteppingAction {
public:
    POLARSteppingAction();
    ~POLARSteppingAction();

    void UserSteppingAction(const G4Step* aStep);

};

#endif
