#ifndef POLARSensitiveDetector_hh
#define POLARSensitiveDetector_hh 1

#include "G4VSensitiveDetector.hh"
#include "G4EmSaturation.hh"
#include "POLARHit.hh"

class POLARSensitiveDetector: public G4VSensitiveDetector {
public:
    POLARSensitiveDetector(G4String SDName, G4String HCName);
    ~POLARSensitiveDetector();

    void Initialize(G4HCofThisEvent* hitCollection);
    G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
    void EndOfEvent(G4HCofThisEvent* hitCollection);

private:
    POLARHitsCollection * fPOLARHitsCollection_;
    G4EmSaturation* fEmSaturation_;

};

#endif
