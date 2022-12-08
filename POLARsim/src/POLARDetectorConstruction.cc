#include "POLARDetectorConstruction.hh"

//Unit
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

//material
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4NistManager.hh"

//geometry
#include "G4Orb.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Trap.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4Sphere.hh"
#include "G4Para.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4MultiUnion.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include <cstdio>
#include <cmath>

//transformation
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"

//sensitive detector
#include "G4SDManager.hh"
#include "POLARSensitiveDetector.hh"

//visual attribute
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

//global function
#include "globals.hh"

// configure
#include "POLARGlobalConfig.hh"

POLARDetectorConstruction::POLARDetectorConstruction() {
    WorldLog_          = NULL;
    ScintillatorLog_   = NULL;
    ModuleLog_         = NULL;
    DetectorLog_       = NULL;
    SpaceLabLog_       = NULL;
    materials_defined_ = false;
}

POLARDetectorConstruction::~POLARDetectorConstruction() {

}

G4VPhysicalVolume* POLARDetectorConstruction::Construct() {
    //-----------------------------------------------------------------------------------------------------------
    //--------Define the detector geometry
    //-----------------------------------------------------------------------------------------------------------

    POLARGlobalConfig* fPOLARGlobalConfig = POLARGlobalConfig::Instance();

    DefineMaterials_();

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    // Construct WorldPhys
    G4double world_hx = (fPOLARGlobalConfig->spacelab ? 7.5*m : 500*mm);
    G4double world_hy = (fPOLARGlobalConfig->spacelab ? 7.5*m : 500*mm);
    G4double world_hz = (fPOLARGlobalConfig->spacelab ? 7.5*m : 500*mm);
    G4Box* WorldBox = new G4Box("WorldBox", world_hx, world_hy, world_hz);
    G4LogicalVolume* WorldLog = new G4LogicalVolume(WorldBox, Vacuum_, "WorldLog");
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    
    // Rotated mother volume
    G4RotationMatrix* RotWorld;
    RotWorld = new G4RotationMatrix();
// RotWorld->rotateY(18.7*deg);         //170206A
// RotWorld->rotateZ(147.3*deg);   
/*    RotWorld->rotateY(84.8*deg);
    RotWorld->rotateZ(67.1 *deg);  */ 
//      RotWorld->rotateY(80.3*deg);
//      RotWorld->rotateZ(-106.9*deg);
//      RotWorld->rotateY(5.6*deg);
//      RotWorld->rotateZ(68.7*deg); 

//         RotWorld->rotateY((54.31)*deg);
//         RotWorld->rotateZ((247.5)*deg);
//  RotWorld->rotateY((27.1)*deg);	//170114 
//  RotWorld->rotateZ(6.1*deg); 
//    RotWorld->rotateY((26.3+3.7)*deg);	//170114 (alt)
//    RotWorld->rotateZ(5.44*deg); 
   
//   RotWorld->rotateY(26.3*deg);	//170114 (alt2)
//   RotWorld->rotateZ((5.44+3.7)*deg); 

//      RotWorld->rotateY(41.8* deg);     //170127C
//      RotWorld->rotateZ(157.6*deg); 
//     RotWorld->rotateY(31.4*deg);     //170305A
//     RotWorld->rotateZ(239.1*deg); 
//      RotWorld->rotateY(70.6*deg);     //170207A
//      RotWorld->rotateZ(-2.2*deg); 

//        RotWorld->rotateY(24.3*deg);     //161218A
//        RotWorld->rotateZ(356.6*deg); 
       
//        RotWorld->rotateY(77.7*deg);     //161218B
//        RotWorld->rotateZ(252.4*deg); 


    G4LogicalVolume* l_Rotate;
    const G4String RotateName = "Rotated World";
    G4Box* s_Rotate = new G4Box(RotateName, world_hx, world_hy, world_hz);
    l_Rotate = new G4LogicalVolume(s_Rotate, Vacuum_, RotateName);					// Using air for on-ground measurements
    G4VPhysicalVolume* p_Rotate = new G4PVPlacement(RotWorld,	G4ThreeVector(0.*mm,0.0*mm,0*mm),  l_Rotate,"Rotated World", WorldLog, false, 0);  

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    if (!fPOLARGlobalConfig->primary_only) {
        // Construct Detector
        ConstructModule_();
        ConstructDetector_();
        new G4PVPlacement(NULL, G4ThreeVector(), DetectorLog_, "Detector", l_Rotate, false, 0);
        // Construct SpaceLab
        if (fPOLARGlobalConfig->spacelab) {
            G4cout << "Building SpaceLab TG02 ..." << G4endl;
            ConstructSpaceLab_();
            new G4PVPlacement(NULL, G4ThreeVector(), SpaceLabLog_, "SpaceLab", l_Rotate, false, 0);
        }
        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    }

    // Return WorldPhys
    G4VPhysicalVolume* WorldPhys = new G4PVPlacement(NULL, G4ThreeVector(), WorldLog, "World", 0, false, 0);
    return WorldPhys;

}

void POLARDetectorConstruction::ConstructSDandField() {
    POLARGlobalConfig* fPOLARGlobalConfig = POLARGlobalConfig::Instance();
    if (!fPOLARGlobalConfig->primary_only) {
        POLARSensitiveDetector* POLAR_SD = new POLARSensitiveDetector("POLARSD", "POLARHitsCollection");
        G4SDManager::GetSDMpointer()->AddNewDetector(POLAR_SD);
        SetSensitiveDetector(ScintillatorLog_, POLAR_SD);
    }
}

void POLARDetectorConstruction::DefineMaterials_() {
    //-----------------------------------------------------------------------------------------------------------
    //--------Define Material
    //-----------------------------------------------------------------------------------------------------------

    POLARGlobalConfig* fPOLARGlobalConfig = POLARGlobalConfig::Instance();

    G4NistManager* man = G4NistManager::Instance();

    G4double density;                         // matter density
    G4double fractionmass;                    // compound proportion
    G4String name;                            // material name
    G4int ncomponents;                        // component number
    G4int natoms;                             // atom proportion

    //Vacuum_
    Vacuum_ = man->FindOrBuildMaterial("G4_Galactic");  // for space
    // Vacuum_ = man->FindOrBuildMaterial("G4_AIR");       // for ground

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    G4Element* elC  = man->FindOrBuildElement("C");
    G4Element* elH  = man->FindOrBuildElement("H");
    G4Element* elO  = man->FindOrBuildElement("O");
    G4Element* elN  = man->FindOrBuildElement("N");
    G4Element* elAl = man->FindOrBuildElement("Al");
    G4Element* elSi = man->FindOrBuildElement("Si");
    G4Element* elFe = man->FindOrBuildElement("Fe");
    G4Element* elCr = man->FindOrBuildElement("Cr");
    G4Element* elMn = man->FindOrBuildElement("Mn");
    G4Element* elMg = man->FindOrBuildElement("Mg");
    G4Element* elSb = man->FindOrBuildElement("Sb");
    G4Element* elCu = man->FindOrBuildElement("Cu");
    G4Element* elZn = man->FindOrBuildElement("Zn");
    G4Element* elTi = man->FindOrBuildElement("Ti");
    G4Element* elZr = man->FindOrBuildElement("Zr");
    G4Element* elV  = man->FindOrBuildElement("V");

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    // EJ_248_
    EJ_248_ = new G4Material(name = "EJ_248", density = 1.023*g/cm3, ncomponents = 2);
    EJ_248_->AddElement(elH, fractionmass = 0.08483);  // 5.18E22 atoms per cm3
    EJ_248_->AddElement(elC, fractionmass = 0.91517);  // 4.69E22 atoms per cm3
    EJ_248_->GetIonisation()->SetBirksConstant(fPOLARGlobalConfig->birks_constant*mm/MeV);
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    // define other materials that are not sensitive

    // CarbonFiber
    CarbonFiber_ = new G4Material(name = "CarbonFiber", density = 1.7*g/cm3, ncomponents = 3);
    CarbonFiber_->AddElement(elC, natoms = 3);
    CarbonFiber_->AddElement(elH, natoms = 1);
    CarbonFiber_->AddElement(elN, natoms = 1);

    // Steel
    Steel_ = man->FindOrBuildMaterial("G4_STAINLESS-STEEL");

    // BorosilicateGlass
    BorosilicateGlass_ = man->FindOrBuildMaterial("G4_Pyrex_Glass");

    // Rubber
    Rubber_ = man->FindOrBuildMaterial("G4_RUBBER_NEOPRENE");

    // PEEK_Plastic
    PEEK_Plastic_ = new G4Material(name = "PEEK_Plastic", density = 1.32*g/cm3, ncomponents = 3);
    PEEK_Plastic_->AddElement(elC, natoms = 20);
    PEEK_Plastic_->AddElement(elH, natoms = 12);
    PEEK_Plastic_->AddElement(elO, natoms = 3);

    // G4_Al
    G4_Al_ = man->FindOrBuildMaterial("G4_Al");

    // G4_Cu
    G4_Cu_ = man->FindOrBuildMaterial("G4_Cu");

    // PMTdynodeMat
    PMTdynodeMat_ = new G4Material(name = "PMTdynodeMat", density = 7.6*g/cm3, ncomponents = 2);
    PMTdynodeMat_->AddElement(elFe, 99.0*perCent);
    PMTdynodeMat_->AddElement(elSb, 1.0*perCent);

    // VikuitiESR
    VikuitiESR_ = new G4Material(name = "VikuitiESR", density = 1.29*g/cm3, ncomponents = 2);
    VikuitiESR_->AddElement(elC, natoms = 1);
    VikuitiESR_->AddElement(elH, natoms = 2);

    // Polythene
    Polythene_Al_ = new G4Material(name = "Polythene", density = 0.7*g/cm3, ncomponents = 3);
    Polythene_Al_->AddElement(elC, fractionmass = 85.68*perCent);
    Polythene_Al_->AddElement(elH, fractionmass = 14.28*perCent);
    Polythene_Al_->AddElement(elAl, fractionmass = 0.04*perCent);

    // Resin
    Resin_ = new G4Material(name = "Resin", density = 1.2*g/cm3, ncomponents = 3);
    Resin_->AddElement(elC, natoms = 2);
    Resin_->AddElement(elH, natoms = 3);
    Resin_->AddElement(elO, natoms = 1);

    //FR4 for PCB
    //Silicon Oxide
    G4Material* SiliconOxide = new G4Material("SiliconOxide", density = 2.65*g/cm3, ncomponents = 2);
    SiliconOxide->AddElement(elSi, natoms=1);
    SiliconOxide->AddElement(elO,  natoms=2);
    //Diglycidyl Ether of Bisphenol A (First compound of epoxy resin Epotek 301-1)
    G4Material* Epoxy_1 = new G4Material("Epoxy_1", density = 1.16*g/cm3, ncomponents = 3);
    Epoxy_1->AddElement(elC, natoms=19);
    Epoxy_1->AddElement(elH, natoms=20);
    Epoxy_1->AddElement(elO, natoms=4);
    //1,4-Butanediol Diglycidyl Ether (Second compound of epoxy resin Epotek 301-1)
    G4Material* Epoxy_2 = new G4Material("Epoxy_2", density = 1.10*g/cm3, ncomponents = 3);
    Epoxy_2->AddElement(elC, natoms=10);
    Epoxy_2->AddElement(elH, natoms=18);
    Epoxy_2->AddElement(elO, natoms=4);
    //1,6-Hexanediamine 2,2,4-trimetyl (Third compound of epoxy resin Epotek 301-1)
    G4Material* Epoxy_3 = new G4Material("Epoxy_3", density = 1.16*g/cm3, ncomponents = 3);
    Epoxy_3->AddElement(elC, natoms=9);
    Epoxy_3->AddElement(elH, natoms=22);
    Epoxy_3->AddElement(elN, natoms=2);
    //Epoxy resin Epotek 301-1
    G4Material* Epoxy_Resin = new G4Material("Epoxy_Resin", density = 1.19*g/cm3, ncomponents = 3);
    Epoxy_Resin->AddMaterial(Epoxy_1, fractionmass=56*perCent);
    Epoxy_Resin->AddMaterial(Epoxy_2, fractionmass=24*perCent);
    Epoxy_Resin->AddMaterial(Epoxy_3, fractionmass=20*perCent);
    //FR4 PCB material
    FR4_ = new G4Material("FR4", density = 1.8*g/cm3, ncomponents=2);
    FR4_->AddMaterial(SiliconOxide, fractionmass=60*perCent);
    FR4_->AddMaterial(Epoxy_Resin,  fractionmass=40*perCent);

    // Sylgard 184
    Sylgard_184_ = new G4Material(name = "Sylgard_184", density = 0.965*g/cm3, ncomponents = 4);
    Sylgard_184_->AddElement(elC,  natoms = 2);
    Sylgard_184_->AddElement(elH,  natoms = 6);
    Sylgard_184_->AddElement(elO,  natoms = 1);
    Sylgard_184_->AddElement(elSi, natoms = 1);

    // AlloyAl_7075
    AlloyAl_7075_ = new G4Material(name = "AlloyAl_7075", density = 2.81*g/cm3, ncomponents = 9);
    AlloyAl_7075_->AddElement(elSi, fractionmass = 0.4*perCent);
    AlloyAl_7075_->AddElement(elFe, fractionmass = 0.5*perCent);
    AlloyAl_7075_->AddElement(elCu, fractionmass = 1.6*perCent);
    AlloyAl_7075_->AddElement(elMn, fractionmass = 0.3*perCent);
    AlloyAl_7075_->AddElement(elMg, fractionmass = 2.5*perCent);
    AlloyAl_7075_->AddElement(elCr, fractionmass = 0.23*perCent);
    AlloyAl_7075_->AddElement(elZn, fractionmass = 5.6*perCent);
    AlloyAl_7075_->AddElement(elTi, fractionmass = 0.2*perCent);
    AlloyAl_7075_->AddElement(elAl, fractionmass = 88.67*perCent);

    // AlloyAl_2219
    AlloyAl_2219_ = new G4Material("AlloyAl_2219", 2.85*g/cm3, 7);
    AlloyAl_2219_->AddElement(elV,  fractionmass = 0.1*perCent);
    AlloyAl_2219_->AddElement(elSi, fractionmass = 0.1*perCent);
    AlloyAl_2219_->AddElement(elFe, fractionmass = 0.15*perCent);
    AlloyAl_2219_->AddElement(elZr, fractionmass = 0.15*perCent);
    AlloyAl_2219_->AddElement(elMn, fractionmass = 0.3*perCent);
    AlloyAl_2219_->AddElement(elCu, fractionmass = 6.3*perCent);
    AlloyAl_2219_->AddElement(elAl, fractionmass = 92.90*perCent);

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    materials_defined_ = true;
}

void POLARDetectorConstruction::ConstructModule_() {
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    // build sensitive detector: ScintillatorLog_
    if (!materials_defined_) return;
    // define ModuleBox
    G4double MBzPlane[6] = {-47.895*mm, -46.0*mm,
                            -46.0*mm, 3.642*mm,
                            3.642*mm, (3.642 + 178.0)*mm};
    G4double MBrOuter[6] = {59.6*mm / 2.0, 59.6*mm / 2,
                            55.5*mm / 2.0, 55.5*mm / 2,
                            52.5*mm / 2.0, 52.5*mm / 2};
    G4double MBrInner[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    G4Polyhedra* ModuleBox = new G4Polyhedra("ModuleBox", 45.0*deg, 360.0*deg, 4, 6, MBzPlane, MBrInner, MBrOuter);
    G4Tubs* ModuleTopHole = new G4Tubs("ModuleTopHole", 0.0*mm, 1.25*mm, 5.0*mm, 0.0*deg, 360.0*deg);
    G4VSolid* ModuleSolid = new G4SubtractionSolid("ModuleSolid",
            ModuleBox, ModuleTopHole, NULL, G4ThreeVector(0.0*mm, 0.0*mm, (3.642 + 178.0)*mm));
    ModuleLog_ = new G4LogicalVolume(ModuleSolid, Vacuum_, "ModuleLog");
    G4VisAttributes* VisAtt_ModuleLog = new G4VisAttributes(true, G4Colour(1.0, 1.0, 0, 1.0));
    VisAtt_ModuleLog->SetForceWireframe(true);
    // VisAtt_ModuleLog->SetForceSolid(true);
    ModuleLog_->SetVisAttributes(VisAtt_ModuleLog);

    // define Scintillator
    G4double SBzPlane[4] = {0.0*mm, 5.0*mm, (5.0 + 166.0)*mm, (5.0 + 166.0 + 5.0)*mm};
    G4double SBrOuter[4] = {5.0*mm / 2.0, 5.8*mm / 2.0, 5.8*mm / 2.0, 5.0*mm / 2.0};
    G4double SBrInner[4] = {0.0, 0.0, 0.0, 0.0};
    G4Polyhedra* ScintillatorBox = new G4Polyhedra("ScintillatorBox", 45.0*deg, 360.0*deg, 4, 4, SBzPlane, SBrInner, SBrOuter);
    ScintillatorLog_ = new G4LogicalVolume(ScintillatorBox, EJ_248_, "ScintillatorLog");
    G4VisAttributes* VisAtt_ScintillatorLog = new G4VisAttributes(true, G4Color(0.0, 0.66, 1.0, 0.2));
    VisAtt_ScintillatorLog->SetForceSolid(true);
    ScintillatorLog_->SetVisAttributes(VisAtt_ScintillatorLog);

    // place 64 Scintillators
    G4double BarD = 6.08*mm;
    const G4double x_start = -(BarD * 3 + BarD / 2);
    const G4double y_start = -(BarD * 3 + BarD / 2);
    for (int j = 0; j < 64; j++) {
        G4double cur_x = x_start + (7 - j % 8) * BarD;
        G4double cur_y = y_start + (j / 8) * BarD;
        new G4PVPlacement(NULL, G4ThreeVector(cur_x, cur_y, 0), ScintillatorLog_, "Scintillator", ModuleLog_, false, j);
    }
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    // build nonsensitive components

    // (1) Top Grid
    G4VisAttributes* VisAtt_TopGrid = new G4VisAttributes(true, G4Color(0.8, 0.8, 0.8, 1.0));
    VisAtt_TopGrid->SetForceSolid(true);
    G4double TGIBzPlane[2] = {0.0*mm, 3.0*mm};
    G4double TGIBrOuter[2] = {6.08*mm / 2, 6.08*mm / 2};
    G4double TGIBrInner[2] = {5.68*mm / 2, 5.20*mm / 2};
    G4Polyhedra* TopGridInnerBox = new G4Polyhedra("TopGridInnerBox", 45.0*deg, 360.0*deg, 4, 2, TGIBzPlane, TGIBrInner, TGIBrOuter);
    G4LogicalVolume* TopGridInnerLog = new G4LogicalVolume(TopGridInnerBox, PEEK_Plastic_, "TopGridInnerLog");
    TopGridInnerLog->SetVisAttributes(VisAtt_TopGrid);
    G4double TGOBzPlane[2] = {0.0*mm, 3.0*mm};
    G4double TGOBrOuter[2] = {50.0*mm / 2, 50.0*mm / 2};
    G4double TGOBrInner[2] = {50.0*mm / 2 - (1.12 - 0.88 / 2)*mm, 50.0*mm / 2 - (1.12 - 0.88 / 2)*mm};
    G4Polyhedra* TopGridOuterBox = new G4Polyhedra("TopGridOuterBox", 45.0*deg, 360.0*deg, 4, 2, TGOBzPlane, TGOBrInner, TGOBrOuter);
    G4LogicalVolume* TopGridOuterLog = new G4LogicalVolume(TopGridOuterBox, PEEK_Plastic_, "TopGridOuterLog");
    TopGridOuterLog->SetVisAttributes(VisAtt_TopGrid);
    new G4PVPlacement(NULL, G4ThreeVector(0, 0, (176.0 - 3.0)*mm), TopGridOuterLog, "TopGridOuter", ModuleLog_, false, 0);
    for (int j = 0; j < 64; j++) {
        G4double cur_x = x_start + (7 - j % 8) * BarD;
        G4double cur_y = y_start + (j / 8) * BarD;
        new G4PVPlacement(NULL, G4ThreeVector(cur_x, cur_y, (176.0 - 3.0)*mm), TopGridInnerLog, "TopGridInner", ModuleLog_, false, j);
    }

    // (2) Top Seal
    G4VisAttributes* VisAtt_TopSeal = new G4VisAttributes(true, G4Color(0.4, 0.4, 0.4, 1.0));
    VisAtt_TopSeal->SetForceSolid(true);
    G4RotationMatrix* zRot45  = new G4RotationMatrix();
    zRot45->rotateZ(45.0*deg);
    G4RotationMatrix* zRot90  = new G4RotationMatrix();
    zRot90->rotateZ(90.0*deg);
    G4RotationMatrix* zRot180 = new G4RotationMatrix();
    zRot180->rotateZ(180.0*deg);
    G4RotationMatrix* zRot270 = new G4RotationMatrix();
    zRot270->rotateZ(270.0*deg);
    G4double TSBSzPlane[2] = {0.0*mm, (3.0 - 0.065)*mm};
    G4double TSBSrOuter[2] = {11.0*mm / 2, 11.0*mm / 2};
    G4double TSBSrInner[2] = {0.0*mm, 0.0*mm};
    G4Polyhedra* TopSealSBox = new G4Polyhedra("TopSealSBox", 45.0*deg, 360.0*deg, 4, 2, TSBSzPlane, TSBSrInner, TSBSrOuter);
    G4LogicalVolume* TopSealSLog = new G4LogicalVolume(TopSealSBox, Rubber_, "TopSealSBox");
    TopSealSLog->SetVisAttributes(VisAtt_TopSeal);
    new G4PVPlacement(NULL, G4ThreeVector(-18.0*mm, 18.0*mm, (176.0 + 0.065)*mm), TopSealSLog, "TopSealS", ModuleLog_, false, 0);
    new G4PVPlacement(NULL, G4ThreeVector(-6.0*mm, 18.0*mm, (176.0 + 0.065)*mm), TopSealSLog, "TopSealS", ModuleLog_, false, 1);
    new G4PVPlacement(NULL, G4ThreeVector(6.0*mm, 18.0*mm, (176.0 + 0.065)*mm), TopSealSLog, "TopSealS", ModuleLog_, false, 2);
    new G4PVPlacement(NULL, G4ThreeVector(18.0*mm, 18.0*mm, (176.0 + 0.065)*mm), TopSealSLog, "TopSealS", ModuleLog_, false, 3);
    new G4PVPlacement(NULL, G4ThreeVector(-18.0*mm, 6.0*mm, (176.0 + 0.065)*mm), TopSealSLog, "TopSealS", ModuleLog_, false, 4);
    new G4PVPlacement(NULL, G4ThreeVector(18.0*mm, 6.0*mm, (176.0 + 0.065)*mm), TopSealSLog, "TopSealS", ModuleLog_, false, 5);
    new G4PVPlacement(NULL, G4ThreeVector(-18.0*mm, -6.0*mm, (176.0 + 0.065)*mm), TopSealSLog, "TopSealS", ModuleLog_, false, 6);
    new G4PVPlacement(NULL, G4ThreeVector(18.0*mm, -6.0*mm, (176.0 + 0.065)*mm), TopSealSLog, "TopSealS", ModuleLog_, false, 7);
    new G4PVPlacement(NULL, G4ThreeVector(-18.0*mm, -18.0*mm, (176.0 + 0.065)*mm), TopSealSLog, "TopSealS", ModuleLog_, false, 8);
    new G4PVPlacement(NULL, G4ThreeVector(-6.0*mm, -18.0*mm, (176.0 + 0.065)*mm), TopSealSLog, "TopSealS", ModuleLog_, false, 9);
    new G4PVPlacement(NULL, G4ThreeVector(6.0*mm, -18.0*mm, (176.0 + 0.065)*mm), TopSealSLog, "TopSealS", ModuleLog_, false, 10);
    new G4PVPlacement(NULL, G4ThreeVector(18.0*mm, -18.0*mm, (176.0 + 0.065)*mm), TopSealSLog, "TopSealS", ModuleLog_, false, 11);
    G4double TSubBzPlane[2] = {-1.0*mm, (4.0 - 0.065)*mm};
    G4double TSubBrOuter[2] = {11.0*mm / 2, 11.0*mm / 2};
    G4double TSubBrInner[2] = {0.0*mm, 0.0*mm};
    G4Polyhedra* TopSubBox = new G4Polyhedra("TopSubBox", 45.0*deg, 360.0*deg, 4, 2, TSubBzPlane, TSubBrInner, TSubBrOuter);
    G4VSolid*    TopSealPBox = new G4SubtractionSolid("TopSealPBox", TopSealSBox, TopSubBox, zRot45, G4ThreeVector(6.89*mm, 6.89*mm, 0.0*mm));
    G4LogicalVolume* TopSealPLog = new G4LogicalVolume(TopSealPBox, Rubber_, "TopSealPBox");
    TopSealPLog->SetVisAttributes(VisAtt_TopSeal);
    new G4PVPlacement(NULL, G4ThreeVector(-6.0*mm, -6.0*mm, (176.0 + 0.065)*mm), TopSealPLog, "TopSealP", ModuleLog_, false, 0);
    new G4PVPlacement(zRot270, G4ThreeVector(6.0*mm, -6.0*mm, (176.0 + 0.065)*mm), TopSealPLog, "TopSealP", ModuleLog_, false, 1);
    new G4PVPlacement(zRot180, G4ThreeVector(6.0*mm, 6.0*mm, (176.0 + 0.065)*mm), TopSealPLog, "TopSealP", ModuleLog_, false, 2);
    new G4PVPlacement(zRot90, G4ThreeVector(-6.0*mm, 6.0*mm, (176.0 + 0.065)*mm), TopSealPLog, "TopSealP", ModuleLog_, false, 3);

    // (3) Top Plate
    G4VisAttributes* VisAtt_TopPlate = new G4VisAttributes(true, G4Color(0.6, 0.6, 0.6, 1.0));
    VisAtt_TopPlate->SetForceSolid(true);
    G4double TPBzPlane[2] = {0.0*mm, 1.642*mm};
    G4double TPBrOuter[2] = {50.5*mm / 2, 50.5*mm / 2};
    G4double TPBrInner[2] = {0.0*mm, 0.0*mm};
    G4Polyhedra* TopPlateBox = new G4Polyhedra("TopPlateBox", 45.0*deg, 360.0*deg, 4, 2, TPBzPlane, TPBrInner, TPBrOuter);
    G4Tubs* TopPlateTubs = new G4Tubs("TopPlateTubs", 0.0*mm, 3.0*mm, 1.0*mm, 0.0*deg, 360.0*deg);
    G4VSolid* TopPlateSolid1 = new G4UnionSolid("TopPlateSolid1", TopPlateBox, TopPlateTubs, NULL, G4ThreeVector(0.0*mm, 0.0*mm, -1*mm));
    G4Tubs* TopPlateHole = new G4Tubs("TopPlateHole", 0.0*mm, 1.25*mm, 3*mm, 0*deg, 360.0*deg);
    G4VSolid* TopPlateSolid2 = new G4SubtractionSolid("TopPlateSolid2", TopPlateSolid1, TopPlateHole, NULL, G4ThreeVector());
    G4LogicalVolume* TopPlateLog = new G4LogicalVolume(TopPlateSolid2, PEEK_Plastic_, "TopPlateLog");
    TopPlateLog->SetVisAttributes(VisAtt_TopPlate);
    new G4PVPlacement(NULL, G4ThreeVector(0.0*mm, 0.0*mm, 179.0*mm), TopPlateLog, "TopPlate", ModuleLog_, false, 0);

    // (4) Bottom Grid
    G4VisAttributes* VisAtt_BottomGrid = new G4VisAttributes(true, G4Color(0.7, 0.7, 0.7, 1.0));
    VisAtt_BottomGrid->SetForceSolid(true);
    G4double BGIBzPlane[2] = {0.0*mm, 3.0*mm};
    G4double BGIBrOuter[2] = {6.08*mm / 2, 6.08*mm / 2};
    G4double BGIBrInner[2] = {5.20*mm / 2, 5.68*mm / 2};
    G4Polyhedra* BottomGridInnerBox = new G4Polyhedra("BottomGridInnerBox", 45.0*deg, 360.0*deg, 4, 2, BGIBzPlane, BGIBrInner, BGIBrOuter);
    G4LogicalVolume* BottomGridInnerLog = new G4LogicalVolume(BottomGridInnerBox, PEEK_Plastic_, "BottomGridInnerLog");
    BottomGridInnerLog->SetVisAttributes(VisAtt_BottomGrid);
    for (int j = 0; j < 64; j++) {
        G4double cur_x = x_start + (7 - j % 8) * BarD;
        G4double cur_y = y_start + (j / 8) * BarD;
        new G4PVPlacement(NULL, G4ThreeVector(cur_x, cur_y, 0.0*mm), BottomGridInnerLog, "BottomGridInner", ModuleLog_, false, j);
    }
    G4double BGOBzPlane[6] = {-4.0*mm, 0.0*mm, 0.0*mm, 1.0*mm, 1.0*mm, 3.0*mm};
    G4double BGOBrOuter[6] = {53.2*mm / 2, 53.2*mm / 2, 53.20*mm / 2, 53.20*mm / 2, 50.00*mm / 2, 50.00*mm / 2};
    G4double BGOBrInner[6] = {52.2*mm / 2, 52.2*mm / 2, 48.64*mm / 2, 48.64*mm / 2, 48.64*mm / 2, 48.64*mm / 2};
    G4Polyhedra* BottomGridOuterBox = new G4Polyhedra("BottomGridOuterBox", 45.0*deg, 360.0*deg, 4, 6, BGOBzPlane, BGOBrInner, BGOBrOuter);
    G4LogicalVolume* BottomGridOuterLog = new G4LogicalVolume(BottomGridOuterBox, PEEK_Plastic_, "BottomGridOuterLog");
    BottomGridOuterLog->SetVisAttributes(VisAtt_BottomGrid);
    new G4PVPlacement(NULL, G4ThreeVector(0.0*mm, 0.0*mm, 0.0*mm), BottomGridOuterLog, "BottomGridOuter", ModuleLog_, false, 0);

    // (5) BottomSealOpt
    G4VisAttributes* VisAtt_BottomSealOpt = new G4VisAttributes(true, G4Color(0.7, 0.7, 0.0, 0.3));
    VisAtt_BottomSealOpt->SetForceSolid(true);
    G4double BSOBzPlane[2] = {-0.5*mm, 0.0*mm};
    G4double BSOBrOuter[2] = {49.0*mm / 2, 49.0*mm / 2};
    G4double BSOBrInner[2] = {0.0*mm, 0.0*mm};
    G4Polyhedra* BottomSealOptBox = new G4Polyhedra("BottomSealOptBox", 45.0*deg, 360.0*deg, 4, 2, BSOBzPlane, BSOBrInner, BSOBrOuter);
    G4LogicalVolume* BottomSealOptLog = new G4LogicalVolume(BottomSealOptBox, Sylgard_184_, "BottomSealOptLog");
    BottomSealOptLog->SetVisAttributes(VisAtt_BottomSealOpt);
    new G4PVPlacement(NULL, G4ThreeVector(0.0*mm, 0.0*mm, 0.0*mm), BottomSealOptLog, "BottomSealOpt", ModuleLog_, false, 0);

    // (6) PMT glass
    G4VisAttributes* VisAtt_PMTGlass = new G4VisAttributes(true, G4Color(0.7, 0.7, 0.0, 0.3));
    VisAtt_PMTGlass->SetForceSolid(true);
    G4double PGBzPlane[2] = {0.0*mm, 1.0*mm};
    G4double PGBrOuter[2] = {51.2*mm / 2, 51.2*mm / 2};
    G4double PGBrInner[2] = {0.0*mm, 0.0*mm};
    G4Polyhedra* PMTGlassBox = new G4Polyhedra("PMTGlassBox", 45.0*deg, 360.0*deg, 4, 2, PGBzPlane, PGBrInner, PGBrOuter);
    G4LogicalVolume* PMTGlassLog = new G4LogicalVolume(PMTGlassBox, BorosilicateGlass_, "PMTGlassLog");
    PMTGlassLog->SetVisAttributes(VisAtt_PMTGlass);
    new G4PVPlacement(NULL, G4ThreeVector(0.0*mm, 0.0*mm, -1.5*mm), PMTGlassLog, "PMTGlass", ModuleLog_, false, 0);

    // (7) PMT inner box
    G4VisAttributes* VisAtt_PMTInner = new G4VisAttributes(true, G4Color(0.65, 0.64, 0.97, 1.0));
    VisAtt_PMTInner->SetForceSolid(true);
    G4double PIBzPlane[2] = {-15.6*mm, 0.0*mm};
    G4double PIBrOuter[2] = {49.2*mm / 2, 49.2*mm / 2};
    G4double PIBrInner[2] = {0.0*mm, 0.0*mm};
    G4Polyhedra* PMTInnerBox = new G4Polyhedra("PMTInnerBox", 45.0*deg, 360.0*deg, 4, 2, PIBzPlane, PIBrInner, PIBrOuter);
    G4LogicalVolume* PMTInnerLog = new G4LogicalVolume(PMTInnerBox, Vacuum_, "PMTInnerLog");
    PMTInnerLog->SetVisAttributes(VisAtt_PMTInner);
    new G4PVPlacement(NULL, G4ThreeVector(0.0*mm, 0.0*mm, -1.5*mm), PMTInnerLog, "PMTInner", ModuleLog_, false, 0);

    // (8) PMT dynode
    G4VisAttributes* VisAtt_PMTdynode = new G4VisAttributes(true, G4Color(0.5, 0.5, 0.5, 1.0));
    VisAtt_PMTdynode->SetForceSolid(true);
    G4Box* PMTdynodeBox = new G4Box("PMTdynodeBox", 48.0*mm / 2. , 48.0*mm / 2. , 0.3*mm / 2.);
    G4LogicalVolume* PMTdynodeLog = new G4LogicalVolume(PMTdynodeBox, PMTdynodeMat_, "PMTdynodeLog");
    PMTdynodeLog->SetVisAttributes(VisAtt_PMTdynode);
    for (int i = 0; i < 12; i++) {
        new G4PVPlacement(NULL, G4ThreeVector(0.0*mm, 0.0*mm, (-0.3 - i * 0.6)*mm), PMTdynodeLog, "PMTdynode", PMTInnerLog, false, i);
    }

    // (9) PMT outer box
    G4VisAttributes* VisAtt_PMTOuter = new G4VisAttributes(true, G4Color(0.54, 0.54, 0.82, 1.0));
    VisAtt_PMTOuter->SetForceSolid(true);
    G4double POBzPlane[2] = {-14.1*mm, 0.0*mm};
    G4double POBrOuter[2] = {51.2*mm / 2, 51.2*mm / 2};
    G4double POBrInner[2] = {49.2*mm / 2, 49.2*mm / 2};
    G4Polyhedra* PMTOuterBox = new G4Polyhedra("PMTOuterBox", 45.0*deg, 360.0*deg, 4, 2, POBzPlane, POBrInner, POBrOuter);
    G4LogicalVolume* PMTOuterLog = new G4LogicalVolume(PMTOuterBox, Steel_, "PMTOuterLog");
    PMTOuterLog->SetVisAttributes(VisAtt_PMTOuter);
    new G4PVPlacement(NULL, G4ThreeVector(0.0*mm, 0.0*mm, -1.5*mm), PMTOuterLog, "PMTOuter", ModuleLog_, false, 0);

    // (10) FEESidePlate
    G4VisAttributes* VisAtt_FEESidePlate = new G4VisAttributes(true, G4Color(0.70, 0.54, 0.20, 1.0));
    VisAtt_FEESidePlate->SetForceSolid(true);
    G4Box* FEESidePlateBoxUp = new G4Box("FEESidePlateBoxUp", 2.5*mm / 2, 30.0*mm / 2, 16.292*mm / 2);
    G4Trap* FEESidePlateTrapDown = new G4Trap("FEESidePlateTrapDown", 30*mm, 3*mm, 4.505*mm, 2.773*mm);
    G4RotationMatrix* xRot = new G4RotationMatrix();
    xRot->rotateX(-90.0*deg);
    G4VSolid* FEESidePlateSolid1 = new G4UnionSolid("FEESidePlateSolid1",
            FEESidePlateBoxUp, FEESidePlateTrapDown, xRot, G4ThreeVector(0.5695*mm, 0.0*mm, -9.646*mm));
    G4Box* FEESidePlateBoxDown1 = new G4Box("FEESidePlateBoxDown1", 4.505*mm / 2, 30.0*mm / 2, 2.008*mm / 2);
    G4VSolid* FEESidePlateSolid2 = new G4UnionSolid("FEESidePlateSolid2",
            FEESidePlateSolid1, FEESidePlateBoxDown1, NULL, G4ThreeVector(1.0025*mm, 0.0*mm, -12.15*mm));
    G4Box* FEESidePlateBoxDown2 = new G4Box("FEESidePlateBoxDown2", 4.505*mm / 2, 11.8*mm / 2, 1.492*mm / 2);
    G4VSolid* FEESidePlateSolid3 = new G4UnionSolid("FEESidePlateSolid3",
            FEESidePlateSolid2, FEESidePlateBoxDown2, NULL, G4ThreeVector(1.0025*mm, 0.0*mm, -13.9*mm));
    G4Box* FEESidePlateSlotUp = new G4Box("FEESidePlateSlotUp", 1.118*mm, 31.0*mm / 2, 1.7*mm / 2);
    G4VSolid* FEESidePlateSolid4 = new G4SubtractionSolid("FEESidePlateSolid4",
            FEESidePlateSolid3, FEESidePlateSlotUp, NULL, G4ThreeVector(2.5*mm / 2, 0.0*mm, 1.496*mm));
    G4Box* FEESidePlateSlotDown = new G4Box("FEESidePlateSlotDown", 1.118*mm, 31.0*mm / 2, 1.2*mm / 2);
    G4VSolid* FEESidePlateSolid5 = new G4SubtractionSolid("FEESidePlateSolid5",
            FEESidePlateSolid4, FEESidePlateSlotDown, NULL, G4ThreeVector(2.5*mm / 2, 0.0*mm, -6.554*mm));
    G4LogicalVolume* FEESidePlateLog = new G4LogicalVolume(FEESidePlateSolid5, AlloyAl_7075_, "FEESidePlateLog");
    FEESidePlateLog->SetVisAttributes(VisAtt_FEESidePlate);
    new G4PVPlacement(NULL, G4ThreeVector(-24.75*mm, 0.0*mm, -25.354*mm), FEESidePlateLog, "FEESidePlate", ModuleLog_, false, 0);
    G4RotationMatrix* zRot = new G4RotationMatrix();
    zRot->rotateZ(180.0*deg);
    new G4PVPlacement(zRot, G4ThreeVector(24.75*mm, 0.0*mm, -25.354*mm), FEESidePlateLog, "FEESidePlate", ModuleLog_, false, 1);

    // (11) FEERodSept
    G4VisAttributes* VisAtt_FEERodSept = new G4VisAttributes(true, G4Color(0.6, 0.6, 0.6, 1.0));
    VisAtt_FEERodSept->SetForceSolid(true);
    G4double FEERSSzPlane[6] = {-23.5*mm, -17.5*mm, -16.75*mm, 16.75*mm, 17.5*mm, 23.5*mm};
    G4double FEERSSrOuter[6] = {2.0*mm, 2.0*mm, 1.25*mm, 1.25*mm, 2.0*mm, 2.0*mm};
    G4double FEERSSrInner[6] = {0.0*mm, 0.0*mm, 0.0*mm, 0.0*mm, 0.0*mm, 0.0*mm};
    G4Polycone* FEERodSeptSolid = new G4Polycone("FEERodSeptSolid", 0.0*deg, 360.0*deg, 6, FEERSSzPlane, FEERSSrInner, FEERSSrOuter);
    G4LogicalVolume* FEERodSeptLog = new G4LogicalVolume(FEERodSeptSolid, AlloyAl_7075_, "FEERodSeptLog");
    FEERodSeptLog->SetVisAttributes(VisAtt_FEERodSept);
    G4RotationMatrix* yRot = new G4RotationMatrix();
    yRot->rotateY(90.0*deg);
    new G4PVPlacement(yRot, G4ThreeVector(0.0*mm, 0.0*mm, -(25.354 + 2.154)*mm), FEERodSeptLog, "FEERodSept", ModuleLog_, false, 0);

    // (12) PMTStitch
    G4VisAttributes* VisAtt_PMTStitch = new G4VisAttributes(true, G4Color(1.0, 1.0, 0.0, 1.0));
    VisAtt_PMTStitch->SetForceSolid(true);
    G4Tubs* PMTStitchTubs = new G4Tubs("PMTStitchTubs", 0.0*mm, 0.7*mm / 2, 5.908*mm / 2, 0.0*deg, 360.0*deg);
    G4LogicalVolume* PMTStitchLog = new G4LogicalVolume(PMTStitchTubs, G4_Cu_, "PMTStitchLog");
    PMTStitchLog->SetVisAttributes(VisAtt_PMTStitch);
    for (int j = 0; j < 64; j++) {
        G4double cur_x = x_start + (7 - j % 8) * BarD;
        G4double cur_y = y_start + (j / 8) * BarD;
        new G4PVPlacement(NULL, G4ThreeVector(cur_x, cur_y, -20.054*mm), PMTStitchLog, "PMTStitch", ModuleLog_, false, j);
    }

    // (13) PCBHVRC
    G4VisAttributes* VisAtt_PCBHVRC = new G4VisAttributes(true, G4Color(0.0, 1.0, 0.0, 1.0));
    VisAtt_PCBHVRC->SetForceSolid(true);
    G4Box* PCBHVRCBox = new G4Box("PCBHVRCBox", 51.0*mm / 2, 51.0*mm / 2, 1.56*mm / 2);
    G4Box* PCBHVRCSideSlot = new G4Box("PCBHVRCSideSlot", 1.0*mm, 30.0*mm / 2, 1.0*mm);
    G4VSolid* PCBHVRCSolid1 = new G4SubtractionSolid("PCBHVRCSolid1",
            PCBHVRCBox, PCBHVRCSideSlot, NULL, G4ThreeVector(51.0*mm / 2, 0.0*mm, 0.0*mm));
    G4VSolid* PCBHVRCSolid2 = new G4SubtractionSolid("PCBHVRCSolid2",
            PCBHVRCSolid1, PCBHVRCSideSlot, NULL, G4ThreeVector(-51.0*mm / 2, 0.0*mm, 0.0*mm));
    G4Tubs* PCBHVRCCenterHole = new G4Tubs("PCBHVRCCenterHole", 0.0*mm, 3.5*mm, 1.0*mm, 0.0*deg, 360.0*deg);
    G4VSolid* PCBHVRCSolid3 = new G4SubtractionSolid("PCBHVRCSolid3",
            PCBHVRCSolid2, PCBHVRCCenterHole, NULL, G4ThreeVector(0.0*mm, 0.0*mm, 0.0*mm));
    G4LogicalVolume* PCBHVRCLog = new G4LogicalVolume(PCBHVRCSolid3, FR4_, "PCBHVRCLog");
    PCBHVRCLog->SetVisAttributes(VisAtt_PCBHVRC);
    new G4PVPlacement(NULL, G4ThreeVector(0.0*mm, 0.0*mm, -23.788*mm), PCBHVRCLog, "PCBHVRC", ModuleLog_, false, 0);

    // (14) PCBHVRCSlot
    G4VisAttributes* VisAtt_PCBHVRCSlot = new G4VisAttributes(true, G4Color(0.6, 0.6, 0.6, 1.0));
    VisAtt_PCBHVRCSlot->SetForceSolid(true);
    G4Box* PCBHVRCSlotBox = new G4Box("PCBHVRCSlotBox", 28.575*mm / 2, 5.715*mm / 2, 6.74*mm / 2);
    G4LogicalVolume* PCBHVRCSlotLog = new G4LogicalVolume(PCBHVRCSlotBox, Resin_, "PCBHVRCSlotLog");
    PCBHVRCSlotLog->SetVisAttributes(VisAtt_PCBHVRCSlot);
    new G4PVPlacement(NULL, G4ThreeVector(0.0*mm, 12.016*mm, -27.938*mm), PCBHVRCSlotLog, "PCBHVRCSlot", ModuleLog_, false, 0);
    new G4PVPlacement(NULL, G4ThreeVector(0.0*mm, -12.016*mm, -27.938*mm), PCBHVRCSlotLog, "PCBHVRCSlot", ModuleLog_, false, 1);

    // (15) PCBASIC
    G4VisAttributes* VisAtt_PCBASIC = new G4VisAttributes(true, G4Color(0.0, 1.0, 0.0, 1.0));
    VisAtt_PCBASIC->SetForceSolid(true);
    G4Box* PCBASICBox = new G4Box("PCBASICBox", 52.0*mm / 2, 52.0*mm / 2, 1.08*mm / 2);
    G4Box* PCBASICSideSlot = new G4Box("PCBASICSideSlot", 1.5*mm, 33.4*mm / 2, 1.0*mm);
    G4VSolid* PCBASICSolid1 = new G4SubtractionSolid("PCBASICSolid1",
            PCBASICBox, PCBASICSideSlot, NULL, G4ThreeVector(52.0*mm / 2, 0.0*mm, 0.0*mm));
    G4VSolid* PCBASICSolid2 = new G4SubtractionSolid("PCBASICSolid2",
            PCBASICSolid1, PCBASICSideSlot, NULL, G4ThreeVector(-52.0*mm / 2, 0.0*mm, 0.0*mm));
    G4LogicalVolume* PCBASICLog = new G4LogicalVolume(PCBASICSolid2, FR4_, "PCBASICLog");
    PCBASICLog->SetVisAttributes(VisAtt_PCBASIC);
    new G4PVPlacement(NULL, G4ThreeVector(0.0*mm, 0.0*mm, -31.848*mm), PCBASICLog, "PCBASIC", ModuleLog_, false, 0);

    // (16) PCBIF
    G4VisAttributes* VisAtt_PCBIF = new G4VisAttributes(true, G4Color(0.0, 1.0, 0.0, 1.0));
    VisAtt_PCBIF->SetForceSolid(true);
    G4Box* PCBIFBox = new G4Box("PCBIFBox", 52.0*mm / 2, 52.0*mm / 2, 1.492*mm / 2);
    G4Box* PCBIFSideSlot = new G4Box("PCBIFSideSlot", 4.537*mm, 12.0*mm / 2, 1.0*mm);
    G4VSolid* PCBIFSolid1 = new G4SubtractionSolid("PCBIFSolid1",
            PCBIFBox, PCBIFSideSlot, NULL, G4ThreeVector(52.0*mm / 2, 0.0*mm, 0.0*mm));
    G4VSolid* PCBIFSolid2 = new G4SubtractionSolid("PCBIFSolid2",
            PCBIFSolid1, PCBIFSideSlot, NULL, G4ThreeVector(-52.0*mm / 2, 0.0*mm, 0.0*mm));
    G4LogicalVolume* PCBIFLog = new G4LogicalVolume(PCBIFSolid2, FR4_, "PCBIFLog");
    PCBIFLog->SetVisAttributes(VisAtt_PCBIF);
    new G4PVPlacement(NULL, G4ThreeVector(0.0*mm, 0.0*mm, -39.254*mm), PCBIFLog, "PCBIF", ModuleLog_, false, 0);

    // (17) PCBIFSlot
    G4VisAttributes* VisAtt_PCBIFSlot = new G4VisAttributes(true, G4Color(0.6, 0.6, 0.6, 1.0));
    VisAtt_PCBIFSlot->SetForceSolid(true);
    G4Box* PCBIFSlotBox1 = new G4Box("PCBIFSlotBox1", 30.012*mm / 2, 10.16*mm / 2, 7.8*mm / 2);
    G4Box* PCBIFSlotBox2 = new G4Box("PCBIFSlotBox2", 18.0*mm / 2, 4.73*mm / 2, 4.7*mm / 2);
    G4VSolid* PCBIFSlotSolid = new G4UnionSolid("PCBIFSlotSolid",
            PCBIFSlotBox1, PCBIFSlotBox2, NULL, G4ThreeVector(0.0*mm, 7.445*mm, 0.0*mm));
    G4LogicalVolume* PCBIFSlotLog = new G4LogicalVolume(PCBIFSlotSolid, G4_Cu_, "PCBIFSlotLog");
    PCBIFSlotLog->SetVisAttributes(VisAtt_PCBIFSlot);
    new G4PVPlacement(NULL, G4ThreeVector(0.0*mm, -12.733*mm, -43.9*mm), PCBIFSlotLog, "PCBIFSlot", ModuleLog_, false, 0);

    // (18) PMTSeal
    G4VisAttributes* VisAtt_PMTSeal = new G4VisAttributes(true, G4Color(0.6, 0.6, 0.6, 1.0));
    VisAtt_PMTSeal->SetForceSolid(true);
    G4Box* PMTSealBox1 = new G4Box("PMTSealBox1", 53.0*mm / 2, 53.0*mm / 2, 2.0*mm / 2);
    G4Box* PMTSealBox2 = new G4Box("PMTSealBox2", 44.0*mm / 2, 44.0*mm / 2, 2.0*mm);
    G4VSolid* PMTSealSolid = new G4SubtractionSolid("PMTSealSolid",
            PMTSealBox1, PMTSealBox2, NULL, G4ThreeVector(0.0*mm, 0.0*mm, 0.0*mm));
    G4LogicalVolume* PMTSealLog = new G4LogicalVolume(PMTSealSolid, Rubber_, "PMTSealLog");
    PMTSealLog->SetVisAttributes(VisAtt_PMTSeal);
    new G4PVPlacement(NULL, G4ThreeVector(0.0*mm, 0.0*mm, -41.0*mm), PMTSealLog, "PMTSeal", ModuleLog_, false, 0);

    // (19) EndCap
    G4VisAttributes* VisAtt_EndCap = new G4VisAttributes(true, G4Color(0.5, 0.5, 0.5, 1.0));
    VisAtt_EndCap->SetForceSolid(true);
    G4double ECBzPlane[4] = {0.0*mm, 1.5*mm, 1.5*mm, 5.5*mm};
    G4double ECBrOuter[4] = {59.6*mm / 2, 59.6*mm / 2, 53.0*mm / 2, 53.0*mm / 2};
    G4double ECBrInner[4] = {47.0*mm / 2, 43.0*mm / 2, 43.0*mm / 2, 43.0*mm / 2};
    G4Polyhedra* EndCapBox = new G4Polyhedra("EndCapBox", 45.0*deg, 360.0*deg, 4, 4, ECBzPlane, ECBrInner, ECBrOuter);
    G4LogicalVolume* EndCapLog = new G4LogicalVolume(EndCapBox, AlloyAl_7075_, "EndCapLog");
    EndCapLog->SetVisAttributes(VisAtt_EndCap);
    new G4PVPlacement(NULL, G4ThreeVector(0.0*mm, 0.0*mm, -47.5*mm), EndCapLog, "EndCap", ModuleLog_, false, 0);

    // (20) CarbonShell
    G4VisAttributes* VisAtt_CarbonShell = new G4VisAttributes(true, G4Color(0.2, 0.2, 0.2, 1.0));
    VisAtt_CarbonShell->SetForceSolid(true);
    G4double CSBzPlane[8] = {-46.0*mm, 2.642*mm, 2.642*mm, 3.642*mm, 3.642*mm, 180.642*mm, 180.642*mm, 181.642*mm};
    G4double CSBrOuter[8] = {55.5*mm / 2, 55.5*mm / 2, 55.5*mm / 2, 55.5*mm / 2, 52.5*mm / 2, 52.5*mm / 2, 52.5*mm / 2, 52.5*mm / 2};
    G4double CSBrInner[8] = {53.5*mm / 2, 53.5*mm / 2, 50.5*mm / 2, 50.5*mm / 2, 50.5*mm / 2, 50.5*mm / 2, 0.0*mm, 0.0*mm};
    G4Polyhedra* CarbonShellBox = new G4Polyhedra("CarbonShellBox", 45.0*deg, 360.0*deg, 4, 8, CSBzPlane, CSBrInner, CSBrOuter);
    G4Tubs* CarbonShellHole = new G4Tubs("CarbonShellHole", 0.0*mm, 1.6*mm, 3*mm, 0*deg, 360.0*deg);
    G4VSolid* CarbonShellSolid = new G4SubtractionSolid("CarbonShellSolid",
            CarbonShellBox, CarbonShellHole, NULL, G4ThreeVector(0.0*mm, 0.0*mm, 181.142*mm));
    G4LogicalVolume* CarbonShellLog = new G4LogicalVolume(CarbonShellSolid, CarbonFiber_, "CarbonShellLog");
    CarbonShellLog->SetVisAttributes(VisAtt_CarbonShell);
    new G4PVPlacement(NULL, G4ThreeVector(0.0*mm, 0.0*mm, 0.0*mm), CarbonShellLog, "CarbonShell", ModuleLog_, false, 0);

    // (21) ESR
    G4VisAttributes* VisAtt_ESR = new G4VisAttributes(true, G4Color(1.0, 1.0, 1.0, 0.1));
    VisAtt_ESR->SetForceSolid(true);
    G4Box* ESRInnerBox = new G4Box("ESRInnerBox", 0.065*mm / 2, 6.0*mm / 2, 166.0*mm / 2);
    G4LogicalVolume* ESRInnerLog = new G4LogicalVolume(ESRInnerBox, VikuitiESR_, "ESRInnerLog");
    ESRInnerLog->SetVisAttributes(VisAtt_ESR);
    G4RotationMatrix* ESRzRot = new G4RotationMatrix();
    ESRzRot->rotateZ(90.0*deg);
    for (int i = 0; i < 7; i++) {
        for (int j = 0; j < 8; j++) {
            G4double cur_x = -3.0 * BarD + i * BarD;
            G4double cur_y = -3.5 * BarD + j * BarD;
            new G4PVPlacement(NULL, G4ThreeVector(cur_x, cur_y, 88.0*mm), ESRInnerLog, "ESRInner", ModuleLog_, false, 2 * (i * 8 + j));
            new G4PVPlacement(ESRzRot, G4ThreeVector(cur_y, cur_x, 88.0*mm), ESRInnerLog, "ESRInner", ModuleLog_, false, 2 * (i * 8 + j) + 1);
        }
    }
    G4Box* ESROuterBox = new G4Box("ESROuterBox",  0.065*mm / 2, 48.57*mm / 2, 166.0*mm / 2);
    G4LogicalVolume* ESROuterLog = new G4LogicalVolume(ESROuterBox, VikuitiESR_, "ESROuterLog");
    ESROuterLog->SetVisAttributes(VisAtt_ESR);
    new G4PVPlacement(NULL, G4ThreeVector(-4.0 * BarD, 0.0*mm, 88.0*mm), ESROuterLog, "ESROuter", ModuleLog_, false, 0);
    new G4PVPlacement(NULL, G4ThreeVector( 4.0 * BarD, 0.0*mm, 88.0*mm), ESROuterLog, "ESROuter", ModuleLog_, false, 1);
    new G4PVPlacement(ESRzRot, G4ThreeVector(0.0*mm, -4.0 * BarD, 88.0*mm), ESROuterLog, "ESROuter", ModuleLog_, false, 2);
    new G4PVPlacement(ESRzRot, G4ThreeVector(0.0*mm,  4.0 * BarD, 88.0*mm), ESROuterLog, "ESROuter", ModuleLog_, false, 3);
    G4Box* ESRTopBox = new G4Box("ESRTopBox", 48.6*mm / 2, 48.6*mm / 2, 0.065*mm / 2);
    G4LogicalVolume* ESRTopLog = new G4LogicalVolume(ESRTopBox, VikuitiESR_, "ESRTopLog");
    ESRTopLog->SetVisAttributes(VisAtt_ESR);
    new G4PVPlacement(NULL, G4ThreeVector(0.0*mm,  0.0*mm, (176.0 + 0.065 / 2)*mm), ESRTopLog, "ESRTop", ModuleLog_, false, 0);

}

void POLARDetectorConstruction::ConstructDetector_() {
    if (ModuleLog_ == NULL) return;
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    // define DetectorBox
    G4double DBzPlane[2] = {-74.0*mm, 190.0*mm};
    G4double DBrOuter[2] = {460.0*mm / 2.0, 460.0*mm / 2.0};
    G4double DBrInner[2] = {0.0, 0.0};
    G4Polyhedra* DetectorBox = new G4Polyhedra("DetectorBox", 45.0*deg, 360.0*deg, 4, 2, DBzPlane, DBrInner, DBrOuter);
    DetectorLog_ = new G4LogicalVolume(DetectorBox, Vacuum_, "DetectorLog");
    G4VisAttributes* VisAtt_DetectorLog = new G4VisAttributes(true, G4Colour(1.0, 0.0, 0.0, 1.0));
    VisAtt_DetectorLog->SetForceWireframe(false);
    VisAtt_DetectorLog->SetVisibility(false);
    DetectorLog_->SetVisAttributes(VisAtt_DetectorLog);

    
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    // build other components in the detector but outside modules

    // (1) OBOXGrid
    G4double ModD = 60.0*mm;
    const G4double x_start = -2 * ModD;
    const G4double y_start = -2 * ModD;
    
    G4VisAttributes* VisAtt_OBOXGrid = new G4VisAttributes(true, G4Colour(0.97, 0.97, 0.97, 1.0));
    VisAtt_OBOXGrid->SetForceSolid(true);
    G4double OGIBzPlane[2] = {-46.0*mm, 6.0*mm};
    G4double OGIBrOuter[2] = {60.0*mm / 2, 60.0*mm / 2};
    G4double OGIBrInner[2] = {56.0*mm / 2, 56.0*mm / 2};
    G4Polyhedra* OBOXGridInnerBox = new G4Polyhedra("OBOXGridInnerBox", 45.0*deg, 360.0*deg, 4, 2, OGIBzPlane, OGIBrInner, OGIBrOuter);
    G4LogicalVolume* OBOXGridInnerLog = new G4LogicalVolume(OBOXGridInnerBox, AlloyAl_7075_, "OBOXGridInnerLog");
    OBOXGridInnerLog->SetVisAttributes(VisAtt_OBOXGrid);
    for (int i = 0; i < 25; i++) {
        G4double cur_x = x_start + (i / 5) * ModD;
        G4double cur_y = y_start + (4 - i % 5) * ModD;
        new G4PVPlacement(NULL, G4ThreeVector(cur_x, cur_y, 0), OBOXGridInnerLog, "OBOXGridInner", DetectorLog_, false, i);
    }
    G4double OGOBzPlane[6] = {-64.0*mm, -56.55*mm, -56.55*mm, 1.0*mm, 1.0*mm, 6.0*mm};
    G4double OGOBrOuter[6] = {458.5*mm / 2, 458.5*mm / 2, 415.0*mm / 2, 415.0*mm / 2, 415.0*mm / 2, 415.0*mm / 2};
    G4double OGOBrInner[6] = {407.0*mm / 2, 407.0*mm / 2, 407.0*mm / 2, 407.0*mm / 2, 300.0*mm / 2, 300.0*mm / 2};
    G4Polyhedra* OBOXGridOuterBox = new G4Polyhedra("OBOXGridOuterBox", 45.0*deg, 360.0*deg, 4, 6, OGOBzPlane, OGOBrInner, OGOBrOuter);
    G4Box* OBOXSlot1 = new G4Box("OBOXSlot1", 8.0*mm, 22.25*mm / 2, 13.25*mm / 2);
    G4VSolid* OBOXSolid1 = new G4SubtractionSolid("OBOXSolid1",
            OBOXGridOuterBox, OBOXSlot1, NULL, G4ThreeVector(-412.0*mm / 2, 74.875*mm, -25.125*mm));
    G4Box* OBOXSlot2 = new G4Box("OBOXSlot2", 8.0*mm, 32.5*mm / 2, 13.5*mm / 2);
    G4VSolid* OBOXSolid2 = new G4SubtractionSolid("OBOXSolid2",
            OBOXSolid1, OBOXSlot2, NULL, G4ThreeVector(-412.0*mm / 2, -72.0*mm, -25.0*mm));
    G4Box* OBOXSlot3 = new G4Box("OBOXSlot3", 8.0*mm, 26.9*mm / 2, 7.52*mm / 2);
    G4VSolid* OBOXSolid3 = new G4SubtractionSolid("OBOXSolid3",
            OBOXSolid2, OBOXSlot3, NULL, G4ThreeVector(-412.0*mm / 2, -150.0*mm, -24.86*mm));
    G4LogicalVolume* OBOXGridOuterLog = new G4LogicalVolume(OBOXSolid3, AlloyAl_7075_, "OBOXGridOuterLog");
    OBOXGridOuterLog->SetVisAttributes(VisAtt_OBOXGrid);
    new G4PVPlacement(NULL, G4ThreeVector(0.0*mm, 0.0*mm, 0.0*mm), OBOXGridOuterLog, "OBOXGridOuter", DetectorLog_, false, 0);

    // place 25 Modules
    for (int i = 0; i < 25; i++) {
        G4double cur_x = x_start + (i / 5) * ModD;
        G4double cur_y = y_start + (4 - i % 5) * ModD;
        new G4PVPlacement(NULL, G4ThreeVector(cur_x, cur_y, 0), ModuleLog_, "Module", DetectorLog_, false, i);
    }
    
    
    // (2) OBOXSpacer
    G4VisAttributes* VisAtt_OBOXSpacer = new G4VisAttributes(true, G4Colour(0.7, 0.7, 0.7, 1.0));
    VisAtt_OBOXSpacer->SetForceSolid(true);
    G4double OSBzPlane[2] = {0.0*mm, 1.0*mm};
    G4double OSBrOuter[2] = {322.0*mm / 2, 322.0*mm / 2};
    G4double OSBrInner[2] = {296.0*mm / 2, 296.0*mm / 2};
    G4Polyhedra* OBOXSpacerBox = new G4Polyhedra("OBOXSpacerBox", 45.0*deg, 360.0*deg, 4, 2, OSBzPlane, OSBrInner, OSBrOuter);
    G4LogicalVolume* OBOXSpacerLog = new G4LogicalVolume(OBOXSpacerBox, AlloyAl_7075_, "OBOXSpacerLog");
    OBOXSpacerLog->SetVisAttributes(VisAtt_OBOXSpacer);
    new G4PVPlacement(NULL, G4ThreeVector(0.0*mm, 0.0*mm, 6.0*mm), OBOXSpacerLog, "OBOXSpacer", DetectorLog_, false, 0);

    // (3) OBOXCarbonShell
    G4VisAttributes* VisAtt_OBOXCarbonShell = new G4VisAttributes(true, G4Color(0.2, 0.2, 0.2, 1.0));
    VisAtt_OBOXCarbonShell->SetForceSolid(true);
    G4double OCSBSzPlane[4] = {0.0*mm, 3.0*mm, 3.0*mm, 175.0*mm};
    G4double OCSBSrOuter[4] = {323.0*mm / 2, 323.0*mm / 2, 300.0*mm / 2, 300.0*mm / 2};
    G4double OCSBSrInner[4] = {294.0*mm / 2, 294.0*mm / 2, 294.0*mm / 2, 294.0*mm / 2};
    G4Polyhedra* OBOXCarbonShellSurBox = new G4Polyhedra("OBOXCarbonShellSurBox", 45.0*deg, 360.0*deg, 4, 4, OCSBSzPlane, OCSBSrInner, OCSBSrOuter);
    G4LogicalVolume* OBOXCarbonShellSurLog = new G4LogicalVolume(OBOXCarbonShellSurBox, CarbonFiber_, "OBOXCarbonShellSurLog");
    OBOXCarbonShellSurLog->SetVisAttributes(VisAtt_OBOXCarbonShell);
    new G4PVPlacement(NULL, G4ThreeVector(0.0*mm, 0.0*mm, 7.0*mm), OBOXCarbonShellSurLog, "OBOXCarbonShellSur", DetectorLog_, false, 0);
    G4Box* OBOXCarbonShellTopBox = new G4Box("OBOXCarbonShellTopBox", 60.0*mm / 2, 60.0*mm / 2, 3.0*mm / 2);
    G4double OTSHzPlane[3] = {-5.0*mm, 1.4*mm, 3.4*mm};
    G4double OTSHrOuter[3] = {2.0*mm, 2.0*mm, 4.1875*mm};
    G4double OTSHrInner[3] = {0.0*mm, 0.0*mm, 0.0*mm};
    G4Polycone* OBOXTopScrewHole = new G4Polycone("OBOXTopScrewSolid", 0.0*deg, 360.0*deg, 3, OTSHzPlane, OTSHrInner, OTSHrOuter);
    G4VSolid* OBOXCarbonShellTopSolid = new G4SubtractionSolid("OBOXCarbonShellTopSolid",
            OBOXCarbonShellTopBox, OBOXTopScrewHole, NULL, G4ThreeVector(0.0*mm, 0.0*mm, -3.0*mm / 2));
    G4LogicalVolume* OBOXCarbonShellTopLog = new G4LogicalVolume(OBOXCarbonShellTopSolid, CarbonFiber_, "OBOXCarbonShellTopLog");
    OBOXCarbonShellTopLog->SetVisAttributes(VisAtt_OBOXCarbonShell);
    for (int i = 0; i < 25; i++) {
        G4double cur_x = x_start + (i / 5) * ModD;
        G4double cur_y = y_start + (4 - i % 5) * ModD;
        new G4PVPlacement(NULL, G4ThreeVector(cur_x, cur_y, 183.5*mm), OBOXCarbonShellTopLog, "OBOXCarbonShellTop", DetectorLog_, false, i);
    }

    // (4) OBOXTopScrew
    G4VisAttributes* VisAtt_OBOXTopScrew = new G4VisAttributes(true, G4Color(0.7, 0.7, 0.7, 1.0));
    VisAtt_OBOXTopScrew->SetForceSolid(true);
    G4double OTSSzPlane[4] = {-5.0*mm, -3.5*mm, 1.4*mm, 3.0*mm};
    G4double OTSSrOuter[4] = {0.2*mm, 1.25*mm, 1.25*mm, 7.2*mm / 2};
    G4double OTSSrInner[4] = {0.0*mm, 0.0*mm, 0.0*mm, 0.0*mm};
    G4Polycone* OBOXTopScrewSolid = new G4Polycone("OBOXTopScrewSolid", 0.0*deg, 360.0*deg, 4, OTSSzPlane, OTSSrInner, OTSSrOuter);
    G4LogicalVolume* OBOXTopScrewLog = new G4LogicalVolume(OBOXTopScrewSolid, Steel_, "OBOXTopScrewLog");
    OBOXTopScrewLog->SetVisAttributes(VisAtt_OBOXTopScrew);
    for (int i = 0; i < 25; i++) {
        G4double cur_x = x_start + (i / 5) * ModD;
        G4double cur_y = y_start + (4 - i % 5) * ModD;
        new G4PVPlacement(NULL, G4ThreeVector(cur_x, cur_y, 182.0*mm), OBOXTopScrewLog, "OBOXTopScrew", DetectorLog_, false, i);
    }

    // (5) AbsorberInner
    G4VisAttributes* VisAtt_AbsorberInner = new G4VisAttributes(true, G4Color(0.7, 0.7, 0.7, 1.0));
    VisAtt_AbsorberInner->SetForceSolid(true);
    G4Tubs* AbsorberInnerTubs = new G4Tubs("AbsorberInnerTubs", 0.0*mm, 5.2*mm / 2, 10.0*mm / 2, 0.0*deg, 360.0*deg);
    G4LogicalVolume* AbsorberInnerLog = new G4LogicalVolume(AbsorberInnerTubs, Steel_, "AbsorberInnerLog");
    AbsorberInnerLog->SetVisAttributes(VisAtt_AbsorberInner);
    for (int i = 0; i < 5; i++) {
        new G4PVPlacement(NULL, G4ThreeVector((-219.0 + i * 88.0)*mm, -219.0*mm, -51.55*mm), AbsorberInnerLog, "AbsorberInner", DetectorLog_, false, 8 * i + 0);
        new G4PVPlacement(NULL, G4ThreeVector((-219.0 + i * 88.0)*mm, -219.0*mm, -69.00*mm), AbsorberInnerLog, "AbsorberInner", DetectorLog_, false, 8 * i + 1);
        new G4PVPlacement(NULL, G4ThreeVector( 219.0*mm, (-219.0 + i * 88.0)*mm, -51.55*mm), AbsorberInnerLog, "AbsorberInner", DetectorLog_, false, 8 * i + 2);
        new G4PVPlacement(NULL, G4ThreeVector( 219.0*mm, (-219.0 + i * 88.0)*mm, -69.00*mm), AbsorberInnerLog, "AbsorberInner", DetectorLog_, false, 8 * i + 3);
        new G4PVPlacement(NULL, G4ThreeVector(( 219.0 - i * 88.0)*mm,  219.0*mm, -51.55*mm), AbsorberInnerLog, "AbsorberInner", DetectorLog_, false, 8 * i + 4);
        new G4PVPlacement(NULL, G4ThreeVector(( 219.0 - i * 88.0)*mm,  219.0*mm, -69.00*mm), AbsorberInnerLog, "AbsorberInner", DetectorLog_, false, 8 * i + 5);
        new G4PVPlacement(NULL, G4ThreeVector(-219.0*mm, ( 219.0 - i * 88.0)*mm, -51.55*mm), AbsorberInnerLog, "AbsorberInner", DetectorLog_, false, 8 * i + 6);
        new G4PVPlacement(NULL, G4ThreeVector(-219.0*mm, ( 219.0 - i * 88.0)*mm, -69.00*mm), AbsorberInnerLog, "AbsorberInner", DetectorLog_, false, 8 * i + 7);
    }

    // (6) AbsorberOuter
    G4VisAttributes* VisAtt_AbsorberOuter = new G4VisAttributes(true, G4Color(0.8, 0.0, 0.0, 1.0));
    VisAtt_AbsorberOuter->SetForceSolid(true);
    G4Tubs* AbsorberOuterTubs = new G4Tubs("AbsorberOuterTubs", 5.2*mm / 2, 21.5*mm / 2, 10.0*mm / 2, 0.0*deg, 360.0*deg);
    G4LogicalVolume* AbsorberOuterLog = new G4LogicalVolume(AbsorberOuterTubs, Rubber_, "AbsorberOuterLog");
    AbsorberOuterLog->SetVisAttributes(VisAtt_AbsorberOuter);
    for (int i = 0; i < 5; i++) {
        new G4PVPlacement(NULL, G4ThreeVector((-219.0 + i * 88.0)*mm, -219.0*mm, -51.55*mm), AbsorberOuterLog, "AbsorberOuter", DetectorLog_, false, 8 * i + 0);
        new G4PVPlacement(NULL, G4ThreeVector((-219.0 + i * 88.0)*mm, -219.0*mm, -69.00*mm), AbsorberOuterLog, "AbsorberOuter", DetectorLog_, false, 8 * i + 1);
        new G4PVPlacement(NULL, G4ThreeVector( 219.0*mm, (-219.0 + i * 88.0)*mm, -51.55*mm), AbsorberOuterLog, "AbsorberOuter", DetectorLog_, false, 8 * i + 2);
        new G4PVPlacement(NULL, G4ThreeVector( 219.0*mm, (-219.0 + i * 88.0)*mm, -69.00*mm), AbsorberOuterLog, "AbsorberOuter", DetectorLog_, false, 8 * i + 3);
        new G4PVPlacement(NULL, G4ThreeVector(( 219.0 - i * 88.0)*mm,  219.0*mm, -51.55*mm), AbsorberOuterLog, "AbsorberOuter", DetectorLog_, false, 8 * i + 4);
        new G4PVPlacement(NULL, G4ThreeVector(( 219.0 - i * 88.0)*mm,  219.0*mm, -69.00*mm), AbsorberOuterLog, "AbsorberOuter", DetectorLog_, false, 8 * i + 5);
        new G4PVPlacement(NULL, G4ThreeVector(-219.0*mm, ( 219.0 - i * 88.0)*mm, -51.55*mm), AbsorberOuterLog, "AbsorberOuter", DetectorLog_, false, 8 * i + 6);
        new G4PVPlacement(NULL, G4ThreeVector(-219.0*mm, ( 219.0 - i * 88.0)*mm, -69.00*mm), AbsorberOuterLog, "AbsorberOuter", DetectorLog_, false, 8 * i + 7);
    }

    // (7) BackPlate
    G4VisAttributes* VisAtt_BackPlate = new G4VisAttributes(true, G4Color(0.7, 0.7, 0.7, 1.0));
    VisAtt_BackPlate->SetForceSolid(true);
    G4Box * BackPlateBox = new G4Box("BackPlateBox", 407.0*mm / 2, 407.0*mm / 2, 3.0*mm / 2);
    G4LogicalVolume* BackPlateLog = new G4LogicalVolume(BackPlateBox, AlloyAl_7075_, "BackPlateLog");
    BackPlateLog->SetVisAttributes(VisAtt_BackPlate);
    new G4PVPlacement(NULL, G4ThreeVector(0.0*mm, 0.0*mm, -62.5*mm), BackPlateLog, "BackPlate", DetectorLog_, false, 0);

    // (8) ThermalCoating
    G4VisAttributes* VisAtt_ThermalCoating = new G4VisAttributes(true, G4Color(0.0, 0.0, 0.5, 0.1));
    VisAtt_ThermalCoating->SetForceSolid(true);
    G4double TCBzPlane[4] = {0.0*mm, 175.0*mm, 175.0*mm, 178.0*mm};
    G4double TCBrInner[4] = {300.0*mm / 2, 300.0*mm / 2, 0.0*mm, 0.0*mm};
    G4double TCBrOuter[4] = {306.0*mm / 2, 306.0*mm / 2, 306.0*mm / 2, 306.0*mm / 2};
    G4Polyhedra* ThermalCoatingBox = new G4Polyhedra("ThermalCoatingBox", 45.0*deg, 360.0*deg, 4, 4, TCBzPlane, TCBrInner, TCBrOuter);
    G4LogicalVolume* ThermalCoatingLog = new G4LogicalVolume(ThermalCoatingBox, Polythene_Al_, "ThermalCoatingLog");
    ThermalCoatingLog->SetVisAttributes(VisAtt_ThermalCoating);
    new G4PVPlacement(NULL, G4ThreeVector(0.0*mm, 0.0*mm, 10.0*mm), ThermalCoatingLog, "ThermalCoating", DetectorLog_, false, 0);

    // (9) Na22 Position
    G4VisAttributes* VisAtt_Na22Pos = new G4VisAttributes(true, G4Color(1.0, 0.0, 0.0, 1.0));
    VisAtt_Na22Pos->SetForceSolid(true);
    G4Box* Na22PosBox = new G4Box("Na22PosBox", 0.5*mm, 0.5*mm, 0.5*mm);
    G4LogicalVolume* Na22PosLog = new G4LogicalVolume(Na22PosBox, Vacuum_, "Na22PosLog");
    Na22PosLog->SetVisAttributes(VisAtt_Na22Pos);
    G4ThreeVector Na22Pos_ID_1( 86.9*mm,  85.4*mm, 88.0*mm);
    G4ThreeVector Na22Pos_ID_2(-86.8*mm,  86.6*mm, 88.0*mm);
    G4ThreeVector Na22Pos_ID_3(-86.9*mm, -82.5*mm, 88.0*mm);
    G4ThreeVector Na22Pos_ID_4( 85.1*mm, -86.9*mm, 88.0*mm);
    new G4PVPlacement(NULL, Na22Pos_ID_1, Na22PosLog, "Na22Pos", DetectorLog_, false, 1);
    new G4PVPlacement(NULL, Na22Pos_ID_2, Na22PosLog, "Na22Pos", DetectorLog_, false, 2);
    new G4PVPlacement(NULL, Na22Pos_ID_3, Na22PosLog, "Na22Pos", DetectorLog_, false, 3);
    new G4PVPlacement(NULL, Na22Pos_ID_4, Na22PosLog, "Na22Pos", DetectorLog_, false, 4);

}

void POLARDetectorConstruction::ConstructSpaceLab_() {
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    // define SpaceLab

    G4Box* SpaceLabBox = new G4Box("SpaceLabBox", 6.5*m, 6.5*m, 6.5*m);
    G4double DSBzPlane[2] = {-74.0*mm, 190.0*mm};
    G4double DSBrOuter[2] = {460.0*mm / 2.0, 460.0*mm / 2.0};
    G4double DSBrInner[2] = {0.0, 0.0};
    G4Polyhedra* DetectorSpaceBox = new G4Polyhedra("DetectorSpaceBox", 45.0*deg, 360.0*deg, 4, 2, DSBzPlane, DSBrInner, DSBrOuter);
    G4VSolid* SpaceLabWorld = new G4SubtractionSolid("SpaceLabSolid", SpaceLabBox, DetectorSpaceBox, NULL, G4ThreeVector());
    SpaceLabLog_ = new G4LogicalVolume(SpaceLabWorld, Vacuum_, "SpaceLabLog");
    G4VisAttributes* VisAtt_SpaceLabWord = new G4VisAttributes(true, G4Colour(0.0, 0.0, 1.0, 1.0));
    VisAtt_SpaceLabWord->SetForceWireframe(true);
    SpaceLabLog_->SetVisAttributes(VisAtt_SpaceLabWord);

    // (1) SpaceBody
    G4VisAttributes* VisAtt_SpaceLabBody = new G4VisAttributes(true, G4Color(0.8, 0.8, 0.8, 1.0));
    VisAtt_SpaceLabBody->SetForceSolid(true);
    G4double SLBzPlane[8] = {-3930.0*mm, -3747.0*mm, -3747.0*mm, -730.0*mm, 773.0*mm, 5117.0*mm, 5617.0*mm, 6317.0*mm};
    G4double SLBrOuter[8] = {2800.0*mm / 2, 2800.0*mm / 2, 2800.0*mm / 2, 2800.0*mm / 2, 3350.0*mm / 2, 3350.0*mm / 2, 1500.0*mm / 2, 1500.0*mm / 2};
    G4double SLBrInner[8] = {0.0*mm / 2,    0.0*mm / 2,    2434.0*mm / 2, 2434.0*mm / 2, 2984.0*mm / 2, 2984.0*mm / 2, 1134.0*mm / 2, 1134.0*mm / 2};
    G4Polycone* SpaceLabBody = new G4Polycone("SpaceLabBody", 0.0*deg, 360.0*deg, 8, SLBzPlane, SLBrInner, SLBrOuter);
    G4LogicalVolume* SpaceLabBodyLog = new G4LogicalVolume(SpaceLabBody, AlloyAl_2219_, "SpaceLabBodyLog");
    SpaceLabBodyLog->SetVisAttributes(VisAtt_SpaceLabBody);
    G4RotationMatrix* xRot2 = new G4RotationMatrix();
    xRot2->rotateX(-90.0*deg);
    new G4PVPlacement(xRot2, G4ThreeVector(0.0*mm, 0.0*mm, -1681.0*mm), SpaceLabBodyLog, "SpaceLabBody", SpaceLabLog_, false, 0);

    // (2) Detector Base
    G4double DBPzPlane[2] = {-7.0*mm, 0.0*mm};
    G4double DBPrOuter[2] = {227.0*mm, 227.0*mm};
    G4double DBPrInner[2] = {0.0*mm,   0.0*mm  };
    G4Polyhedra* DetectorBasePlatform = new G4Polyhedra("DetectorBasePlatform", 45.0*deg, 360.0*deg, 4, 2, DBPzPlane, DBPrInner, DBPrOuter);
    G4LogicalVolume* DetectorBasePlatformLog = new G4LogicalVolume(DetectorBasePlatform, AlloyAl_2219_, "DetectorBasePlatformLog");
    DetectorBasePlatformLog->SetVisAttributes(VisAtt_SpaceLabBody);
    new G4PVPlacement(NULL, G4ThreeVector(0.0*mm, 0.0*mm, -74.0*mm), DetectorBasePlatformLog, "DetectorBasePlatform", SpaceLabLog_, false, 0);
    G4double DBLLzPlane[2] = {-125.5*mm, -7.0*mm};
    G4double DBLLrOuter[2] = {10.0*mm, 10.0*mm};
    G4double DBLLrInner[2] = {0.0*mm,   0.0*mm  };
    G4Polyhedra* DetectorBaseLegLong = new G4Polyhedra("DetectorBaseLegLong", 45.0*deg, 360.0*deg, 4, 2, DBLLzPlane, DBLLrInner, DBLLrOuter);
    G4LogicalVolume* DetectorBaseLegLongLog = new G4LogicalVolume(DetectorBaseLegLong, AlloyAl_2219_, "DetectorBaseLegLongLog");
    DetectorBaseLegLongLog->SetVisAttributes(VisAtt_SpaceLabBody);
    new G4PVPlacement(NULL, G4ThreeVector(217.0*mm, 217.0*mm, -74.0*mm), DetectorBaseLegLongLog, "DetectorBaseLegLong", SpaceLabLog_, false, 0);
    new G4PVPlacement(NULL, G4ThreeVector(-217.0*mm, 217.0*mm, -74.0*mm), DetectorBaseLegLongLog, "DetectorBaseLegLong", SpaceLabLog_, false, 1);
    G4double DBLSzPlane[2] = {-45.5*mm, -7.0*mm};
    G4double DBLSrOuter[2] = {10.0*mm, 10.0*mm};
    G4double DBLSrInner[2] = {0.0*mm,   0.0*mm  };
    G4Polyhedra* DetectorBaseLegShort = new G4Polyhedra("DetectorBaseLegShort", 45.0*deg, 360.0*deg, 4, 2, DBLSzPlane, DBLSrInner, DBLSrOuter);
    G4LogicalVolume* DetectorBaseLegShortLog = new G4LogicalVolume(DetectorBaseLegShort, AlloyAl_2219_, "DetectorBaseLegShortLog");
    DetectorBaseLegShortLog->SetVisAttributes(VisAtt_SpaceLabBody);
    new G4PVPlacement(NULL, G4ThreeVector(-217.0*mm, -217.0*mm, -74.0*mm), DetectorBaseLegShortLog, "DetectorBaseLegShort", SpaceLabLog_, false, 0);
    new G4PVPlacement(NULL, G4ThreeVector(217.0*mm, -217.0*mm, -74.0*mm), DetectorBaseLegShortLog, "DetectorBaseLegShort", SpaceLabLog_, false, 1);

    // (3) StarTracker
    G4VisAttributes* VisAtt_StarTracker = new G4VisAttributes(true, G4Color(0.9, 0.9, 0.9, 1.0));
    VisAtt_StarTracker->SetForceSolid(true);
    G4RotationMatrix* yRot1 = new G4RotationMatrix();
    yRot1->rotateY(11.0*deg);
    // Star Tracker Base
    G4double first_z = 2.0*mm;
    G4double slope = (3350.0 - 2800.0) / 2.0 / (773.0 + 730.0);
    G4double radius = 1487.0*mm;
    G4Box* StarTrackerBaseSolid1 = new G4Box("StarTrackerBaseSolid1", 72.0*mm, 2.5*mm, (first_z * 2 + 0.0*mm   * slope) / 2.0);
    G4LogicalVolume* StarTrackerBaseLog1 = new G4LogicalVolume(StarTrackerBaseSolid1, G4_Al_, "StarTrackerBaseLog1");
    StarTrackerBaseLog1->SetVisAttributes(VisAtt_StarTracker);
    new G4PVPlacement(yRot1,
            G4ThreeVector(-1 * (radius + first_z - 0) * sin(11.0*deg),
                257.0*mm,      (radius + first_z - 0) * cos(11.0*deg) - 1681.0*mm),
            StarTrackerBaseLog1, "StarTrackerBase1", SpaceLabLog_, false, 0);
    G4Box* StarTrackerBaseSolid2 = new G4Box("StarTrackerBaseSolid2", 72.0*mm, 2.5*mm, (first_z * 2 + 145.0*mm * slope) / 2.0);
    G4LogicalVolume* StarTrackerBaseLog2 = new G4LogicalVolume(StarTrackerBaseSolid2, G4_Al_, "StarTrackerBaseLog2");
    StarTrackerBaseLog2->SetVisAttributes(VisAtt_StarTracker);
    new G4PVPlacement(yRot1,
            G4ThreeVector(-1 * (radius + first_z - 145.0*mm * slope / 2.0) * sin(11.0*deg),
                (257.0 + 145.0)*mm, (radius + first_z - 145.0*mm * slope / 2.0) * cos(11.0*deg) - 1681.0*mm),
            StarTrackerBaseLog2, "StarTrackerBase2", SpaceLabLog_, false, 0);
    G4Box* StarTrackerBaseSolid3 = new G4Box("StarTrackerBaseSolid3", 72.0*mm, 2.5*mm, (first_z * 2 + 180.0*mm * slope) / 2.0);
    G4LogicalVolume* StarTrackerBaseLog3 = new G4LogicalVolume(StarTrackerBaseSolid3, G4_Al_, "StarTrackerBaseLog3");
    StarTrackerBaseLog3->SetVisAttributes(VisAtt_StarTracker);
    new G4PVPlacement(yRot1,
            G4ThreeVector(-1 * (radius + first_z - 180.0*mm * slope / 2.0) * sin(11.0*deg),
                (257.0 + 180.0)*mm, (radius + first_z - 180.0*mm * slope / 2.0) * cos(11.0*deg) - 1681.0*mm),
            StarTrackerBaseLog3, "StarTrackerBase3", SpaceLabLog_, false, 0);
    G4Box* StarTrackerBaseSolid4 = new G4Box("StarTrackerBaseSolid4", 72.0*mm, 2.5*mm, (first_z * 2 + 325.0*mm * slope) / 2.0);
    G4LogicalVolume* StarTrackerBaseLog4 = new G4LogicalVolume(StarTrackerBaseSolid4, G4_Al_, "StarTrackerBaseLog4");
    StarTrackerBaseLog4->SetVisAttributes(VisAtt_StarTracker);
    new G4PVPlacement(yRot1,
            G4ThreeVector(-1 * (radius + first_z - 325.0*mm * slope / 2.0) * sin(11.0*deg),
                (257.0 + 325.0)*mm, (radius + first_z - 325.0*mm * slope / 2.0) * cos(11.0*deg) - 1681.0*mm),
            StarTrackerBaseLog4, "StarTrackerBase4", SpaceLabLog_, false, 0);
    // Star Tracker
    G4double STFSzPlane[8] = {0.0*mm, 17.0*mm, 17.0*mm, 54.0*mm, 54.0*mm, 157.0*mm, 157.0*mm, 260.0*mm};
    G4double STFSrOuter[8] = {57.0*mm, 57.0*mm, 57.0*mm, 57.0*mm, 90.5*mm, 90.5*mm, 111.0*mm, 111.0*mm};
    G4double STFSrInner[8] = {0.0*mm,  0.0*mm,  53.5*mm, 53.5*mm, 53.5*mm, 53.5*mm, 53.5*mm,  53.5*mm};
    G4Polycone* StarTrackerFrontSolid = new G4Polycone("StarTrackerFrontSolid", 0.0*deg, 360.0*deg, 8, STFSzPlane, STFSrInner, STFSrOuter);
    G4Box* StarTrackerBackSolid = new G4Box("StarTrackerBackSolid", 85.0*mm, 75.0*mm, 72.5*mm);
    G4VSolid* StarTrackerUpperSolid = new G4UnionSolid("StarTrackerUpperSolid",
            StarTrackerBackSolid, StarTrackerFrontSolid, NULL, G4ThreeVector(0.0*cm, 0.0*cm, 72.5*mm));
    G4Box* StarTrackerHolderSolid1 = new G4Box("StarTrackerHolderSolid1", 72.0*mm, 2.5*mm, 60.0*mm);
    G4Box* StarTrackerHolderSolid2 = new G4Box("StarTrackerHolderSolid2", 72.0*mm, 2.5*mm, 60.0*mm);
    G4RotationMatrix* STyRot1 = new G4RotationMatrix();
    STyRot1->rotateY(66.12*deg);
    G4VSolid* StarTrackerSolid1 = new G4UnionSolid("StarTrackerSolid1",
            StarTrackerHolderSolid1, StarTrackerUpperSolid, STyRot1, G4ThreeVector(30.0*mm, 72.5*mm, 65.0*mm));
    G4VSolid* StarTrackerSolid2 = new G4UnionSolid("StarTrackerSolid2",
            StarTrackerSolid1, StarTrackerHolderSolid2, NULL, G4ThreeVector(0.0*mm, 145.0*mm, 0.0*mm));
    G4LogicalVolume* StarTrackerLog1 = new G4LogicalVolume(StarTrackerSolid2, G4_Al_, "StarTrackerLog1");
    StarTrackerLog1->SetVisAttributes(VisAtt_StarTracker);
    new G4PVPlacement(yRot1,
            G4ThreeVector(-1551.0*mm * sin(11.0*deg), 257.0*mm, 1551.0*mm * cos(11.0*deg) - 1681.0*mm),
            StarTrackerLog1, "StarTracker1", SpaceLabLog_, false, 0);
    G4Box* StarTrackerHolderSolid3 = new G4Box("StarTrackerHolderSolid3", 72.0*mm, 2.5*mm, 60.0*mm);
    G4Box* StarTrackerHolderSolid4 = new G4Box("StarTrackerHolderSolid4", 72.0*mm, 2.5*mm, 60.0*mm);
    G4RotationMatrix* STyRot2 = new G4RotationMatrix();
    STyRot2->rotateY(-66.12*deg);
    G4VSolid* StarTrackerSolid3 = new G4UnionSolid("StarTrackerSolid3",
            StarTrackerHolderSolid3, StarTrackerUpperSolid, STyRot2, G4ThreeVector(-30.0*mm, 72.5*mm, 65.0*mm));
    G4VSolid* StarTrackerSolid4 = new G4UnionSolid("StarTrackerSolid4",
            StarTrackerSolid3, StarTrackerHolderSolid4, NULL, G4ThreeVector(0.0*mm, 145.0*mm, 0.0*mm));
    G4LogicalVolume* StarTrackerLog2 = new G4LogicalVolume(StarTrackerSolid4, G4_Al_, "StarTrackerLog2");
    StarTrackerLog2->SetVisAttributes(VisAtt_StarTracker);
    new G4PVPlacement(yRot1, G4ThreeVector(-1551.0*mm * sin(11.0*deg), 437.0*mm, 1551.0*mm * cos(11.0*deg) - 1681.0*mm),
            StarTrackerLog2, "StarTracker2", SpaceLabLog_, false, 0);

    // (4) Antenna
    G4VisAttributes* VisAtt_Antenna = new G4VisAttributes(true, G4Color(0.9, 0.9, 0.9, 1.0));
    VisAtt_Antenna->SetForceSolid(true);
    // Antenna Base
    G4Tubs* AntennaBaseSolid1 = new G4Tubs("AntennaBaseSolid1", 0.0*cm, 100.0*mm, 9.0*mm, 0.0*deg, 360.0*deg);
    G4Box* AntennaPillarBottomSolid = new G4Box("AntennaPillarBottomSolid", 39.5*mm, 15.0*mm, 167.5*mm);
    G4VSolid* AntennaBaseSolid2 = new G4UnionSolid("AntennaBaseSolid2",
            AntennaBaseSolid1, AntennaPillarBottomSolid, NULL, G4ThreeVector(0.0*mm, 0.0*mm, 176.5*mm));
    G4Para* AntennaPillarMiddleSolid = new G4Para("AntennaPillarMiddleSolid", 15.0*mm, 39.5*mm, 221.5*mm, 0.0*deg, 10.0*deg, 0.0*deg);
    G4RotationMatrix* zRot1 = new G4RotationMatrix();
    zRot1->rotateZ(-90.0*deg);
    G4VSolid* AntennaBaseSolid3 = new G4UnionSolid("AntennaBaseSolid3",
            AntennaBaseSolid2, AntennaPillarMiddleSolid, zRot1, G4ThreeVector(0.0*mm, 221.5*mm * tan(10.0*deg), 565.5*mm));
    G4Box* AntennaPillarTopSolid = new G4Box("AntennaPillarTopSolid", 39.5*mm, 15.0*mm, 88.5*mm);
    G4VSolid* AntennaBaseSolid4 = new G4UnionSolid("AntennaBaseSolid4",
            AntennaBaseSolid3, AntennaPillarTopSolid, NULL, G4ThreeVector(0.0*mm, 443.0*mm * tan(10.0*deg), 875.5*mm));
    G4double AHSzPlane[4] = {0.0*mm, 10.0*mm, 10.0*mm, 70.0*mm};
    G4double AHSrOuter[4] = {60.0*mm, 60.0*mm, 50.0*mm, 50.0*mm};
    G4double AHSrInner[4] = {0.0*mm,  0.0*mm,  0.0*mm, 0.0*mm};
    G4Polycone* AntennaHolderSolid = new G4Polycone("AntennaHolderSolid", 0.0*deg, 360.0*deg, 4, AHSzPlane, AHSrInner,AHSrOuter);
    G4VSolid* AntennaBaseSolid5 = new G4UnionSolid("AntennaBaseSolid5",
            AntennaBaseSolid4, AntennaHolderSolid, NULL, G4ThreeVector(0.0*mm, 443.0*mm * tan(10.0*deg), 964.0*mm));
    G4LogicalVolume* AntennaBaseLog = new G4LogicalVolume(AntennaBaseSolid5, AlloyAl_2219_, "AntennaBaseLog");
    AntennaBaseLog->SetVisAttributes(VisAtt_Antenna);
    new G4PVPlacement(yRot1, G4ThreeVector(-1409.0*mm * sin(11.0*deg), 1414.5*mm, 1409.0*mm * cos(11.0*deg) - 1681.0*mm),
            AntennaBaseLog, "AntennaBase", SpaceLabLog_, false, 0);
    // Antenna Head
    G4Orb* AntennaCenterSolid = new G4Orb("AntennaCenterSolid", 80.0*mm);
    G4double antenna_thickness = 19.0*mm;
    G4Sphere* AntennaOuterSolid = new G4Sphere("AntennaOuterSolid", 632.156*mm - antenna_thickness, 632.156*mm, 0.0*deg, 360.0*deg, 126.0*deg, 180.0*deg);
    G4Tubs* AntennaInnerSolid = new G4Tubs("AntennaInnerSolid", 0.0*cm, 25.0*mm, 191.0*mm, 0.0*deg, 360.0*deg);
    G4VSolid* AntennaHeadSolid1 = new G4UnionSolid("AntennaHeadSolid1",
            AntennaCenterSolid, AntennaOuterSolid, NULL, G4ThreeVector(0.0*mm, 0.0*mm, 680*mm));
    G4VSolid* AntennaHeadSolid2 = new G4UnionSolid("AntennaHeadSolid2",
            AntennaHeadSolid1, AntennaInnerSolid, NULL, G4ThreeVector(0.0*mm, 0.0*mm, (191.0 + 80.0)*mm));
    G4LogicalVolume* AntennaHeadLog = new G4LogicalVolume(AntennaHeadSolid2, CarbonFiber_, "AntennaHeadLog");
    AntennaHeadLog->SetVisAttributes(VisAtt_Antenna);
    G4RotationMatrix* AntennaRot = new G4RotationMatrix();
    POLARGlobalConfig* fPOLARGlobalConfig = POLARGlobalConfig::Instance();
    AntennaRot->rotateY(11.0*deg);
    AntennaRot->rotateZ(-1.0 * fPOLARGlobalConfig->antenna_angle_lr);
    AntennaRot->rotateX(-1.0 * fPOLARGlobalConfig->antenna_angle_ud);
    G4double antenna_height = (2443.0 + 80.0) * mm;
    new G4PVPlacement(AntennaRot,
            G4ThreeVector(-antenna_height * sin(11.0*deg),
                1414.5*mm + 443.0*mm * tan(10.0*deg),
                antenna_height * cos(11.0*deg) - 1681.0*mm),
            AntennaHeadLog, "AntennaHead", SpaceLabLog_, false, 0);

}

