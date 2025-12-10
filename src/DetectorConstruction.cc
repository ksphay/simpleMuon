#include "DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SDManager.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

#include "PlaneSD.hh"

#include <cmath>

DetectorConstruction::DetectorConstruction(G4double halfX,
                                           G4double halfY,
                                           G4double halfZ,
                                           G4double thetaDeg,
                                           G4double phiDeg)
: G4VUserDetectorConstruction(),
  fHalfX(halfX),
  fHalfY(halfY),
  fHalfZ(halfZ),
  fThetaDeg(thetaDeg),
  fPhiDeg(phiDeg),
  fPlaneLogical(nullptr)
{}

DetectorConstruction::~DetectorConstruction() = default;

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    auto* nist = G4NistManager::Instance();
    auto* air  = nist->FindOrBuildMaterial("G4_AIR");

    // -------------------
    // World: box centered at origin with given half-lengths
    // -------------------
    auto* solidWorld = new G4Box("World", fHalfX, fHalfY, fHalfZ);
    auto* logicWorld = new G4LogicalVolume(solidWorld, air, "World");

    auto* physWorld = new G4PVPlacement(
        nullptr,            // no rotation
        {},                 // at (0,0,0)
        logicWorld,
        "World",
        nullptr,
        false,
        0,
        true
    );

    auto* worldVis = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8));
    worldVis->SetForceWireframe(true);
    logicWorld->SetVisAttributes(worldVis);

    // -------------------
    // Sensitive plane
    // -------------------
    G4double minHalfXY   = (fHalfX < fHalfY) ? fHalfX : fHalfY;
    G4double halfXY      = 0.2 * minHalfXY;   // 20% of min transverse size
    G4double halfZPlane  = 0.5 * mm;
    if (halfZPlane > 0.25 * fHalfZ) {
        halfZPlane = 0.25 * fHalfZ;
    }

    auto* solidPlane = new G4Box("Plane", halfXY, halfXY, halfZPlane);
    fPlaneLogical    = new G4LogicalVolume(solidPlane, air, "Plane");

    G4cout << "[DetectorConstruction] Placing sensitive plane with halfX/Y = "
           << halfXY / m << " m, halfZ = " << halfZPlane / m << " m" << G4endl;
    G4cout << "[DetectorConstruction] Plane angles (theta,phi) = ("
           << fThetaDeg << " deg, " << fPhiDeg << " deg)" << G4endl;

    // =========================================================
    // Build rotation with order:
    //
    //   R = Rx(theta) * Rz'(phi)
    //
    // Active view: first rotate around global X by the zenith
    // angle θ, then rotate around the intermediate z' axis by
    // the azimuth φ. This matches the requested convention of
    // tilting the plane around X and then spinning it around its
    // new normal.
    //
    // theta, phi in radians
    // =========================================================
    G4double theta = fThetaDeg * deg;
    G4double phi   = fPhiDeg   * deg;

    G4double cth = std::cos(theta);
    G4double sth = std::sin(theta);
    G4double cph = std::cos(phi);
    G4double sph = std::sin(phi);

    // Columns of R = Rx(theta) * Rz'(phi): world coordinates of the
    // local basis vectors after first tilting by θ around X and then
    // spinning by φ around the new z' axis.
    //
    // From symbolic multiplication:
    //   ex' = ( cosφ,          cosθ sinφ,   sinθ sinφ )
    //   ey' = ( -sinφ,         cosθ cosφ,   sinθ cosφ )
    //   ez' = ( 0,             -sinθ,       cosθ      )
    //
    G4ThreeVector ex(cph,        cth * sph,  sth * sph);
    G4ThreeVector ey(-sph,       cth * cph,  sth * cph);
    G4ThreeVector ez(0.0,        -sth,       cth);

    auto* rot = new G4RotationMatrix(ex, ey, ez);

    // Plane centered at origin with this rotation
    new G4PVPlacement(
        rot,                // rotation (owned by geometry)
        G4ThreeVector(),    // at (0,0,0)
        fPlaneLogical,
        "Plane",
        logicWorld,
        false,
        0,
        true
    );

    auto* planeVis = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));
    planeVis->SetForceSolid(true);
    planeVis->SetVisibility(true);
    planeVis->SetForceAuxEdgeVisible(true);
    fPlaneLogical->SetVisAttributes(planeVis);

    return physWorld;
}

void DetectorConstruction::ConstructSDandField()
{
    if (!fPlaneLogical) return;

    auto* sdManager = G4SDManager::GetSDMpointer();

    auto* planeSD = new PlaneSD("PlaneSD");
    sdManager->AddNewDetector(planeSD);
    fPlaneLogical->SetSensitiveDetector(planeSD);
}
