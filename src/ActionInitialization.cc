// ActionInitialization.cc
#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"

ActionInitialization::ActionInitialization(G4double energyGeV,
                                           G4double thetaDeg,
                                           G4double phiDeg,
                                           G4double sourceRadius,
                                           G4double coneHalfAngleDeg,
                                           G4double planeXY_m)
: G4VUserActionInitialization(),
  fEnergyGeV(energyGeV),
  fThetaDeg(thetaDeg),
  fPhiDeg(phiDeg),
  fSourceRadius(sourceRadius),
  fConeHalfAngleDeg(coneHalfAngleDeg),
  fPlaneXY_m(planeXY_m)
{}

void ActionInitialization::Build() const
{
    auto* prim = new PrimaryGeneratorAction(fEnergyGeV,
                                            fThetaDeg,
                                            fPhiDeg,
                                            fSourceRadius,
                                            fConeHalfAngleDeg,
                                            fPlaneXY_m);
    SetUserAction(prim);
}
