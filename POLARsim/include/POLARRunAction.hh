#ifndef POLARRunAction_h
#define POLARRunAction_h

#include "globals.hh"
#include "G4UserRunAction.hh"

class POLARRunAction: public G4UserRunAction {
private:
    G4String output_file_;
    G4bool   fixed_output_;
    static G4int randomSeed_;
    static G4int runId_;

public:
    POLARRunAction(G4String the_output_file_, G4bool fixed_name = false);
    ~POLARRunAction();

    G4Run* GenerateRun();
    void BeginOfRunAction(const G4Run* aRun);
    void EndOfRunAction(const G4Run* aRun);

};

#endif /* POLARRunAction_h */
