#ifndef POLARPhysicsList_h
#define POLARPhysicsList_h 1

#include "G4VModularPhysicsList.hh"

class POLARPhysicsList: public G4VModularPhysicsList {
public:
    POLARPhysicsList(G4int ver = 1);
    ~POLARPhysicsList();

    void SetCuts();

};

#endif /* POLARPhysicsList_h */
