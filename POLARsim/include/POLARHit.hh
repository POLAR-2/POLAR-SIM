#ifndef POLARHit_hh
#define POLARHit_hh 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "tls.hh"
#include "G4ThreeVector.hh"

class POLARHit: public G4VHit {
public:
    POLARHit();
    ~POLARHit();

    inline void* operator new(size_t);
    inline void  operator delete(void* aHit);

public:
    // hit information
    G4int             TrackID;
    G4String          ParticleName;
    G4int             ParticleCode;
    G4double          GlobalTime;
    G4int             BarID;
    G4int             ModID;
    G4ThreeVector     LocalPosition;
    G4bool            IsEntering;
    G4bool            IsLeaving;
    G4double          PreMomTheta;
    G4double          PreMomPhi;
    G4double          PostMomTheta;
    G4double          PostMomPhi;
    G4double          PreKinEnergy;
    G4double          PostKinEnergy;
    G4String          ProcessName;
    G4double          EnergyDep;
    G4double          EnergyVis;
    G4double          DeltaTime;

private:
    static G4ThreadLocal G4Allocator<POLARHit>* POLARHitAllocator_;

};

typedef G4THitsCollection<POLARHit> POLARHitsCollection;

inline void* POLARHit::operator new(size_t) {
    if (POLARHitAllocator_ == NULL)
        POLARHitAllocator_ = new G4Allocator<POLARHit>();
    return static_cast<void*>(POLARHitAllocator_->MallocSingle());
}

inline void POLARHit::operator delete(void* aHit) {
    POLARHitAllocator_->FreeSingle(static_cast<POLARHit*>(aHit));
}

#endif
