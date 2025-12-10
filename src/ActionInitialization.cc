#include "ActionInitialization.hh"

#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"

ActionInitialization::ActionInitialization(G4double energyGeV,
                                           G4double thetaDeg,
                                           G4double phiDeg,
                                           G4double sourceRadius,
                                           G4double coneHalfAngleDeg)
: G4VUserActionInitialization(),
  fEnergyGeV(energyGeV),
  fThetaDeg(thetaDeg),
  fPhiDeg(phiDeg),
  fSourceRadius(sourceRadius),
  fConeHalfAngleDeg(coneHalfAngleDeg)
{}

void ActionInitialization::BuildForMaster() const
{
    SetUserAction(new RunAction());
}

void ActionInitialization::Build() const
{
    // Primary generator with cone
    SetUserAction(new PrimaryGeneratorAction(
        fEnergyGeV,
        fThetaDeg,
        fPhiDeg,
        fSourceRadius,
        fConeHalfAngleDeg
    ));

    // RunAction for analysis
    SetUserAction(new RunAction());
}
