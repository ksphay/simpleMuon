// PrimaryGeneratorAction.cc
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4RunManager.hh"
#include "Randomize.hh"
#include "CLHEP/Units/SystemOfUnits.h"

#include <cmath>

PrimaryGeneratorAction::PrimaryGeneratorAction(G4double energyGeV,
                                               G4double thetaDeg,
                                               G4double phiDeg,
                                               G4double sourceRadius,
                                               G4double coneHalfAngleDeg,
                                               G4double planeXY_m)
: G4VUserPrimaryGeneratorAction(),
  fGun(nullptr),
  fEnergyGeV(energyGeV),
  fThetaDeg(thetaDeg),
  fPhiDeg(phiDeg),
  fSourceRadius(sourceRadius),
  fConeHalfAngleDeg(coneHalfAngleDeg),
  fPlaneHalfXY(0.0)
{
    // One-particle gun
    fGun = new G4ParticleGun(1);

    auto* particleTable = G4ParticleTable::GetParticleTable();
    auto* muMinus       = particleTable->FindParticle("mu-");
    fGun->SetParticleDefinition(muMinus);

    // Store plane half-length in internal units
    if (planeXY_m > 0.0) {
        fPlaneHalfXY = 0.5 * planeXY_m * m;
    } else {
        fPlaneHalfXY = 0.0;
    }

    G4cout << "PrimaryGeneratorAction: generating mu- with E = "
           << fEnergyGeV << " GeV, "
           << "theta = " << fThetaDeg << " deg, "
           << "phi = "   << fPhiDeg   << " deg, "
           << "source radius = " << fSourceRadius / m << " m, "
           << "cone half-angle = " << fConeHalfAngleDeg << " deg, "
           << "planeXY full = " << planeXY_m << " m (half = "
           << fPlaneHalfXY / m << " m)"
           << G4endl;

    // Default energy; direction/position per-event
    fGun->SetParticleEnergy(fEnergyGeV * GeV);

    // Dummy defaults (overwritten in GeneratePrimaries)
    fGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));
    fGun->SetParticlePosition(G4ThreeVector(0., 0., 1. * m));
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fGun;
}

// -----------------------------------------------------------------------------
// SampleDirectionInCone
// -----------------------------------------------------------------------------
G4ThreeVector PrimaryGeneratorAction::SampleDirectionInCone(const G4ThreeVector& axisWorld) const
{
    G4ThreeVector axis = axisWorld;
    if (axis.mag2() == 0.) {
        axis = G4ThreeVector(0., 0., -1.);   // fallback
    }
    axis = axis.unit();

    // Pencil beam
    if (fConeHalfAngleDeg <= 0.0) {
        return axis;
    }

    G4double alpha  = fConeHalfAngleDeg * deg;
    G4double cosMax = std::cos(alpha);

    // Uniform in cos(theta') between cosMax and 1
    G4double u      = G4UniformRand();
    G4double cosThP = 1.0 - u * (1.0 - cosMax);
    if (cosThP >  1.0)    cosThP = 1.0;
    if (cosThP < cosMax)  cosThP = cosMax;

    G4double sinThP = std::sqrt(std::max(0.0, 1.0 - cosThP * cosThP));

    // Uniform phi' in [0, 2π)
    G4double phiP = 2.0 * CLHEP::pi * G4UniformRand();
    G4double cphP = std::cos(phiP);
    G4double sphP = std::sin(phiP);

    // Build orthonormal basis (u1, u2, axis)
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

// -----------------------------------------------------------------------------
// GeneratePrimaries
// -----------------------------------------------------------------------------
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
    // 1) Get plane basis from DetectorConstruction
    G4ThreeVector ex_plane(1., 0., 0.);
    G4ThreeVector ey_plane(0., 1., 0.);
    G4ThreeVector ez_plane(0., 0., 1.);

    auto* runManager = G4RunManager::GetRunManager();
    const auto* det =
        dynamic_cast<const DetectorConstruction*>(runManager->GetUserDetectorConstruction());

    if (det) {
        det->GetPlaneBasis(ex_plane, ey_plane, ez_plane);
    }

    // Beam axis in WORLD coordinates from (theta, phi)
    G4double theta     = fThetaDeg * deg;
    G4double phiOffset = (fPhiDeg + 90.0) * deg;

    G4double sTh = std::sin(theta);
    G4double cTh = std::cos(theta);
    G4double cPh = std::cos(phiOffset);
    G4double sPh = std::sin(phiOffset);

    // Standard spherical coordinates around world +Z
    G4ThreeVector axisWorld(sTh * cPh, sTh * sPh, cTh);
    if (axisWorld.mag2() == 0.) {
        axisWorld = ez_plane.unit(); // fallback to plane normal
    } else {
        axisWorld = axisWorld.unit();
    }

    // 2) Sample direction inside cone around this axis
    G4ThreeVector direction = SampleDirectionInCone(axisWorld);

    // 3) Sample a random offset in the plane:
    //    disX, disY uniform in [-fPlaneHalfXY, +fPlaneHalfXY]
    G4double disX = 0.0;
    G4double disY = 0.0;

    if (fPlaneHalfXY > 0.0) {
        disX = (2.0 * G4UniformRand() - 1.0) * fPlaneHalfXY;
        disY = (2.0 * G4UniformRand() - 1.0) * fPlaneHalfXY;
    }

    // Offset vector in WORLD coordinates
    G4ThreeVector offset = disX * ex_plane + disY * ey_plane;

    // 4) Place source at distance fSourceRadius BACK along 'direction',
    //    then add the offset in the plane:
    //
    //    ray: x(t) = offset - fSourceRadius*direction + t*direction
    //    at t = fSourceRadius, x = offset  (point in the plane)
    //
    G4ThreeVector position = offset - fSourceRadius * direction;

    fGun->SetParticlePosition(position);
    fGun->SetParticleMomentumDirection(direction);
    fGun->SetParticleEnergy(fEnergyGeV * GeV);

    fGun->GeneratePrimaryVertex(event);
}
