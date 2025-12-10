// ActionInitialization.cc
#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"

ActionInitialization::ActionInitialization(G4double energyGeV,
                                           G4double thetaDeg,
                                           G4double phiDeg,
                                           G4double sourceRadius,
                                           G4double coneHalfAngleDeg,
                                           G4double planeXY_m,
                                           const G4String& outputFileName)
: G4VUserActionInitialization(),
  fEnergyGeV(energyGeV),
  fThetaDeg(thetaDeg),
  fPhiDeg(phiDeg),
  fSourceRadius(sourceRadius),
  fConeHalfAngleDeg(coneHalfAngleDeg),
  fPlaneXY_m(planeXY_m),
  fOutputFileName(outputFileName)
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

    auto* runAction = new RunAction(fOutputFileName);
    SetUserAction(runAction);

    auto* steppingAction = new SteppingAction();
    SetUserAction(steppingAction);
}

void ActionInitialization::BuildForMaster() const
{
    auto* runAction = new RunAction(fOutputFileName);
    SetUserAction(runAction);
}
