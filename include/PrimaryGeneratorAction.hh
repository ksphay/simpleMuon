#ifndef PRIMARY_GENERATOR_ACTION_HH
#define PRIMARY_GENERATOR_ACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"   // <-- include instead of forward-decl
#include "globals.hh"

class G4ParticleGun;
class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction(G4double energyGeV,
                           G4double thetaDeg,
                           G4double phiDeg,
                           G4double sourceRadius,
                           G4double coneHalfAngleDeg);
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event* event) override;

private:
    G4ParticleGun* fGun;

    G4double fEnergyGeV;
    G4double fThetaDeg;        // zenith in WORLD frame (deg)
    G4double fPhiDeg;          // azimuth in WORLD frame (deg)
    G4double fSourceRadius;    // distance from origin along ray
    G4double fConeHalfAngleDeg;// cone half-angle (deg), 0 => pencil beam

    G4ThreeVector SampleDirectionInCone(const G4ThreeVector& axisWorld) const;
};

#endif
