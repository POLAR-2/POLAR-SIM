#ifndef POLARDetectorConstruction_h
#define POLARDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

class G4Material;

class POLARDetectorConstruction : public G4VUserDetectorConstruction {
public:
    POLARDetectorConstruction();
    ~POLARDetectorConstruction();

public:
    G4VPhysicalVolume* Construct();
    void ConstructSDandField();

private: // Materials
    bool materials_defined_;
    G4Material* Vacuum_;
    G4Material* EJ_248_;
    G4Material* CarbonFiber_;
    G4Material* Steel_;
    G4Material* BorosilicateGlass_;
    G4Material* Rubber_;
    G4Material* PEEK_Plastic_;
    G4Material* G4_Al_;
    G4Material* G4_Cu_;
    G4Material* PMTdynodeMat_;
    G4Material* VikuitiESR_;
    G4Material* Polythene_Al_;
    G4Material* Resin_;
    G4Material* FR4_;
    G4Material* Sylgard_184_;
    G4Material* AlloyAl_7075_;
    G4Material* AlloyAl_2219_;

private: // Volumes
    G4LogicalVolume* WorldLog_;
    G4LogicalVolume* ScintillatorLog_;
    G4LogicalVolume* ModuleLog_;
    G4LogicalVolume* DetectorLog_;
    G4LogicalVolume* SpaceLabLog_;

private:
    void DefineMaterials_();
    void ConstructModule_();
    void ConstructDetector_();
    void ConstructSpaceLab_();

};

#endif /* POLARDetectorConstruction_h */
