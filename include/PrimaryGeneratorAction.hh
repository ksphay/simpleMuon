// PrimaryGeneratorAction.hh
#ifndef PRIMARY_GENERATOR_ACTION_HH
#define PRIMARY_GENERATOR_ACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
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
                           G4double coneHalfAngleDeg,
                           G4double planeXY_m);
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event* event) override;

private:
    // Sample direction in cone of half-angle fConeHalfAngleDeg around axisWorld
    G4ThreeVector SampleDirectionInCone(const G4ThreeVector& axisWorld) const;

    G4ParticleGun* fGun;

    G4double fEnergyGeV;
    G4double fThetaDeg;
    G4double fPhiDeg;
    G4double fSourceRadius;       // in world units (mm)
    G4double fConeHalfAngleDeg;   // cone half-angle in degrees
    G4double fPlaneHalfXY;        // half-length of offset square on plane (world units)
};

#endif
