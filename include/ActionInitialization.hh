#ifndef ACTION_INITIALIZATION_HH
#define ACTION_INITIALIZATION_HH

#include "G4VUserActionInitialization.hh"
#include "globals.hh"

class ActionInitialization : public G4VUserActionInitialization
{
public:
    ActionInitialization(G4double energyGeV,
                         G4double thetaDeg,
                         G4double phiDeg,
                         G4double sourceRadius,
                         G4double coneHalfAngleDeg);
    ~ActionInitialization() override = default;

    void BuildForMaster() const override;
    void Build() const override;

private:
    G4double fEnergyGeV;
    G4double fThetaDeg;
    G4double fPhiDeg;
    G4double fSourceRadius;
    G4double fConeHalfAngleDeg;
};

#endif
