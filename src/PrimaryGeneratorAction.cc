#include "PrimaryGeneratorAction.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include "CLHEP/Units/SystemOfUnits.h"

#include <cmath>

PrimaryGeneratorAction::PrimaryGeneratorAction(G4double energyGeV,
                                               G4double thetaDeg,
                                               G4double phiDeg,
                                               G4double sourceRadius,
                                               G4double coneHalfAngleDeg,
                                               G4double hitPlaneHalfSize)
: G4VUserPrimaryGeneratorAction(),
  fGun(nullptr),
  fEnergyGeV(energyGeV),
  fThetaDeg(thetaDeg),
  fPhiDeg(phiDeg),
  fSourceRadius(sourceRadius),
  fConeHalfAngleDeg(coneHalfAngleDeg),
  fHitPlaneHalfSize(hitPlaneHalfSize)
{
    // One-particle gun
    fGun = new G4ParticleGun(1);

    // Default particle: mu-
    auto* particleTable = G4ParticleTable::GetParticleTable();
    auto* muMinus       = particleTable->FindParticle("mu-");
    fGun->SetParticleDefinition(muMinus);

    G4cout << "PrimaryGeneratorAction: generating mu- with E = "
           << fEnergyGeV << " GeV, "
           << "theta = " << fThetaDeg << " deg, "
           << "phi = "   << fPhiDeg   << " deg, "
           << "source radius = " << fSourceRadius / m << " m, "
           << "cone half-angle = " << fConeHalfAngleDeg << " deg, "
           << "hit plane half-size = " << fHitPlaneHalfSize / m << " m."
           << G4endl;
    // Set default energy; direction/position per-event
    fGun->SetParticleEnergy(fEnergyGeV * GeV);

    // Dummy defaults (overwritten in GeneratePrimaries)
    fGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));
    fGun->SetParticlePosition(G4ThreeVector(0., 0., 1. * m));
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fGun;
}

// Sample a direction in a cone centered at ORIGIN,
// axis given in WORLD coordinates (axisWorld),
// with half-angle fConeHalfAngleDeg, uniform in solid angle.
G4ThreeVector PrimaryGeneratorAction::SampleDirectionInCone(const G4ThreeVector& axisWorld) const
{
    G4ThreeVector axis = axisWorld;
    if (axis.mag2() == 0.) {
        axis = G4ThreeVector(0., 0., -1.);   // fallback
    }
    axis = axis.unit();

    // Pencil beam: just return axis
    if (fConeHalfAngleDeg <= 0.0) {
        return axis;
    }

    G4double alpha  = fConeHalfAngleDeg * deg;
    G4double cosMax = std::cos(alpha);

    // Uniform in cos(theta') between cosMax and 1:
    //   cosθ' = 1 - u*(1 - cosMax), u ~ U[0,1)
    G4double u      = G4UniformRand();
    G4double cosThP = 1.0 - u * (1.0 - cosMax);
    if (cosThP >  1.0)    cosThP = 1.0;
    if (cosThP < cosMax)  cosThP = cosMax;

    G4double sinThP = std::sqrt(std::max(0.0, 1.0 - cosThP * cosThP));

    // Uniform phi' in [0, 2π)
    G4double phiP = 2.0 * CLHEP::pi * G4UniformRand();
    G4double cphP = std::cos(phiP);
    G4double sphP = std::sin(phiP);

    // Build an orthonormal basis (u1, u2, axis)
    G4ThreeVector tmp(0., 0., 1.);
    if (std::fabs(axis.dot(tmp)) > 0.9) {
        tmp = G4ThreeVector(1., 0., 0.);
    }

    G4ThreeVector u1 = tmp.cross(axis);
    if (u1.mag2() == 0.) {
        u1 = G4ThreeVector(1., 0., 0.);
    } else {
        u1 = u1.unit();
    }

    G4ThreeVector u2 = axis.cross(u1);
    u2 = u2.unit();

    // Direction in world coords:
    //   d = sinθ' cosφ' u1 + sinθ' sinφ' u2 + cosθ' axis
    G4ThreeVector dir =
        sinThP * cphP * u1 +
        sinThP * sphP * u2 +
        cosThP * axis;

    return dir.unit();
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
    // 1) Central cone axis in WORLD coordinates from (fThetaDeg, fPhiDeg)
    //
    // Standard spherical convention:
    //   theta: zenith from +Z, phi: azimuth from +X toward +Y.
    //
    G4double theta = fThetaDeg * deg;
    G4double phi   = fPhiDeg   * deg;

    G4double sTh = std::sin(theta);
    G4double cTh = std::cos(theta);
    G4double cPh = std::cos(phi);
    G4double sPh = std::sin(phi);

    // Axis direction at the origin (WORLD frame)
    G4ThreeVector axisWorld(sTh * cPh,   // x
                            sTh * sPh,   // y
                            cTh);        // z

    if (axisWorld.mag2() == 0.) {
        axisWorld = G4ThreeVector(0., 0., -1.);
    }
    axisWorld = axisWorld.unit();

    // 2) Sample direction inside cone around this axis
    G4ThreeVector direction = SampleDirectionInCone(axisWorld);

    // 3) Pick a hit point on the Z=0 plane within +/-fHitPlaneHalfSize
    //    in X and Y. A zero or negative half-size defaults to the origin.
    G4double dx = 0.0;
    G4double dy = 0.0;
    if (fHitPlaneHalfSize > 0.0) {
        auto symUniform = [this]() {
            return 2.0 * G4UniformRand() - 1.0; // in [-1, 1]
        };
        dx = symUniform() * fHitPlaneHalfSize;
        dy = symUniform() * fHitPlaneHalfSize;
    }

    G4ThreeVector hitPoint(dx, dy, 0.0);

    // 4) Place source at distance fSourceRadius BACK along this direction
    //    so that the ray passes exactly through the hit point:
    //
    //    ray: x(t) = x0 + t * direction
    //    require x(t0) = hitPoint => x0 = hitPoint - t0 * direction
    //    choose t0 = fSourceRadius
    //
    G4ThreeVector position = hitPoint - fSourceRadius * direction;

    fGun->SetParticlePosition(position);
    fGun->SetParticleMomentumDirection(direction);
    fGun->SetParticleEnergy(fEnergyGeV * GeV);

    fGun->GeneratePrimaryVertex(event);
}
