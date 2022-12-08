#ifndef POLARConfigMessenger_HH
#define POLARConfigMessenger_HH

#include "globals.hh"
#include "G4UImessenger.hh"
#include "POLARGlobalConfig.hh"

class POLARGlobalConfig;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;

class POLARConfigMessenger: public G4UImessenger {
private:
    POLARGlobalConfig* fGlobalConfig_;

public:
    POLARConfigMessenger(POLARGlobalConfig* theGlobalConfig);
    ~POLARConfigMessenger();

    void SetNewValue(G4UIcommand* command, G4String newValue);

private:
    // command list
    G4UIdirectory*             fGlobalConfigDir_;
    G4UIcmdWithADoubleAndUnit* fGlobalConfigHitThrCmd_;
    G4UIcmdWithADoubleAndUnit* fGlobalConfigBarThrCmd_;
    G4UIcmdWithAString*        fGlobalConfigOutDirCmd_;
    G4UIcmdWithADoubleAndUnit* fGlobalConfigPolAglCmd_;
    G4UIcmdWithADoubleAndUnit* fGlobalConfigIncThetaCmd_;
    G4UIcmdWithADoubleAndUnit* fGlobalConfigIncPhiCmd_;
    G4UIcmdWithAnInteger*      fGlobalConfigEventVerCmd_;
    G4UIcmdWithABool*          fGlobalConfigPriOnlyCmd_;
    G4UIcmdWithADouble*        fGlobalConfigBirksCmd_;
    G4UIcmdWithABool*          fGlobalConfigSpaceLabCmd_;
    G4UIcmdWithAnInteger*      fGlobalConfigPhysVerCmd_;
    G4UIcmdWithABool*          fGlobalConfigSimOutMoreCmd_;
    G4UIcmdWithABool*          fGlobalConfigPhysMoreCmd_;
    G4UIcmdWithADoubleAndUnit* fGlobalConfigAntAglUDCmd_;
    G4UIcmdWithADoubleAndUnit* fGlobalConfigAntAglLRCmd_;

};

#endif // POLARConfigMessenger_HH
