#ifndef POLARPrimaryGeneratorAction_h
#define POLARPrimaryGeneratorAction_h 1

#include "globals.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleTable.hh"
#include "G4ThreeVector.hh"

class G4GeneralParticleSource;
class G4ParticleGun;
class ParticleDataFileR;

class POLARPrimaryGeneratorAction: public G4VUserPrimaryGeneratorAction {
private:
    G4bool                   gpsFlag_;
    ParticleDataFileR*       fParticleDataFileR_;
    G4ParticleTable*         fParticleTable_;

    G4ThreeVector            PolarizationVec_;
    G4ThreeVector            EmitPositionVec_;
    G4ThreeVector            MomDirectionVec_;

private:
    G4GeneralParticleSource* fParticleSource_;
    G4ParticleGun*           fParticleGun_;

public:
    POLARPrimaryGeneratorAction(
            G4bool             the_gps_flag = true,
            ParticleDataFileR*  theParticleDataFileR = NULL
            );
    ~POLARPrimaryGeneratorAction();

public:
    void GeneratePrimaries(G4Event* anEvent);

};

#endif /* POLARPrimaryGeneratorAction_h */
