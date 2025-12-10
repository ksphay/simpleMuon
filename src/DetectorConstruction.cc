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
#include "G4Material.hh"

#include "PlaneSD.hh"

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

DetectorConstruction::DetectorConstruction(G4double halfX,
                                           G4double halfY,
                                           G4double halfZ,
                                           G4double thetaDeg,
                                           G4double detXY_m,
                                           const G4String& demPath,
                                           G4double colDX_m,
                                           G4double colDY_m,
                                           G4double colDZ_m,
                                           G4double terrainPhiDeg,
                                           G4double rockDensity_g_cm3)
: G4VUserDetectorConstruction(),
  fStandardRock(nullptr),
  fHalfX(halfX),
  fHalfY(halfY),
  fHalfZ(halfZ),
  fThetaDeg(thetaDeg),
  fPlaneEx(1., 0., 0.),
  fPlaneEy(0., 1., 0.),
  fPlaneEz(0., 0., 1.),
  fPlaneLogical(nullptr),
  fDetXY_m(detXY_m),
  fDEMPath(demPath),
  fColDX_m(colDX_m),
  fColDY_m(colDY_m),
  fColDZ_m(colDZ_m),
  fTerrainPhiDeg(terrainPhiDeg),
  fRockDensity_g_cm3(rockDensity_g_cm3)
{}

DetectorConstruction::~DetectorConstruction() = default;

// -----------------------------------------------------------------------------
// DefineMaterials
// -----------------------------------------------------------------------------
void DetectorConstruction::DefineMaterials()
{
    auto* nist = G4NistManager::Instance();

    // Make sure AIR exists for the world
    nist->FindOrBuildMaterial("G4_AIR");

    // -----------------------------------------------------------------------
    // Standard rock composition
    // Default density: 2.65 g/cm^3 unless user overrides via CLI
    // -----------------------------------------------------------------------
    G4double density_g_cm3 = fRockDensity_g_cm3;
    if (density_g_cm3 <= 0.) {
        density_g_cm3 = 2.65; // fallback
    }
    G4double density = density_g_cm3 * g/cm3;

    G4Element* elO  = nist->FindOrBuildElement("O");
    G4Element* elCa = nist->FindOrBuildElement("Ca");
    G4Element* elC  = nist->FindOrBuildElement("C");
    G4Element* elMg = nist->FindOrBuildElement("Mg");

    // Try to reuse if already defined
    fStandardRock = G4Material::GetMaterial("StandardRock", false);
    if (!fStandardRock) {
        fStandardRock = new G4Material("StandardRock", density, 4, kStateSolid);
        fStandardRock->AddElement(elO,  0.52);
        fStandardRock->AddElement(elCa, 0.27);
        fStandardRock->AddElement(elC,  0.12);
        fStandardRock->AddElement(elMg, 0.09);
    } else {
        if (std::fabs(fStandardRock->GetDensity()/(g/cm3) - density_g_cm3) > 1e-3) {
            G4cout << "[DetectorConstruction] WARNING: StandardRock already exists "
                   << "with density = " << fStandardRock->GetDensity()/(g/cm3)
                   << " g/cm^3; requested " << density_g_cm3
                   << " g/cm^3 will be ignored for this run."
                   << G4endl;
        }
    }

    G4cout << "[DetectorConstruction] StandardRock density (g/cm^3): "
           << fStandardRock->GetDensity()/(g/cm3) << G4endl;
}

// -----------------------------------------------------------------------------
// Construct
// -----------------------------------------------------------------------------
G4VPhysicalVolume* DetectorConstruction::Construct()
{
    DefineMaterials();

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
    worldVis->SetVisibility(true);
    logicWorld->SetVisAttributes(worldVis);

    // -------------------
    // Sensitive plane
    // -------------------
    // Decide half-size from user detXY or fallback to 20% of world transverse
    G4double minHalfXY   = (fHalfX < fHalfY) ? fHalfX : fHalfY;
    G4double halfXY      = 0.0;

    if (fDetXY_m > 0.0) {
        halfXY = 0.5 * fDetXY_m * m;
        if (halfXY > minHalfXY) {
            G4cout << "[DetectorConstruction] WARNING: requested detXY = "
                   << fDetXY_m << " m is larger than world transverse half-size ("
                   << minHalfXY / m << " m). Clamping." << G4endl;
            halfXY = 0.99 * minHalfXY;
        }
    } else {
        // Old behavior: 20% of min world half-size
        halfXY = 0.2 * minHalfXY;
    }

    G4double halfZPlane  = 0.5 * mm;
    if (halfZPlane > 0.25 * fHalfZ) {
        halfZPlane = 0.25 * fHalfZ;
    }

    auto* solidPlane = new G4Box("Plane", halfXY, halfXY, halfZPlane);
    fPlaneLogical    = new G4LogicalVolume(solidPlane, air, "Plane");

    G4cout << "[DetectorConstruction] Placing sensitive plane with halfX/Y = "
           << halfXY / m << " m, halfZ = " << halfZPlane / m << " m" << G4endl;
    G4cout << "[DetectorConstruction] Plane zenith angle theta = "
           << fThetaDeg << " deg; plane phi is fixed at 90 deg." << G4endl;

    // =========================================================
    // Plane normal in WORLD coordinates with phi = 90 deg:
    //
    //   n(theta, phi=90) = (0, sinθ, cosθ)
    //
    // Then build an orthonormal basis:
    //   ez = n
    //   ex = z × n     (fallback to (1,0,0) if parallel)
    //   ey = ez × ex
    // =========================================================
    G4double theta = fThetaDeg * deg;

    G4double sTh = std::sin(theta);
    G4double cTh = std::cos(theta);

    // Plane normal (WORLD), phi fixed at 90 deg
    G4ThreeVector n(0.0, sTh, cTh);
    if (n.mag2() == 0.) {
        n = G4ThreeVector(0., 0., 1.);
    }
    n = n.unit();

    // ex from z × n
    G4ThreeVector zaxis(0., 0., 1.);
    G4ThreeVector ex = zaxis.cross(n);
    if (ex.mag2() < 1e-12) {
        ex = G4ThreeVector(1., 0., 0.);
    }
    ex = ex.unit();

    // ey = n × ex
    G4ThreeVector ey = n.cross(ex);
    if (ey.mag2() < 1e-12) {
        ey = G4ThreeVector(0., 1., 0.);
    }
    ey = ey.unit();

    // Store basis
    fPlaneEx = ex;
    fPlaneEy = ey;
    fPlaneEz = n;  // local +Z' is plane normal

    // Rotation matrix: columns = ex, ey, ez (world coords of local axes)
    auto* rot = new G4RotationMatrix(fPlaneEx, fPlaneEy, fPlaneEz);

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

    // -------------------
    // DEM TERRAIN (columns of StandardRock)
    // -------------------
    if (!fDEMPath.empty() && fColDX_m > 0.0 && fColDY_m > 0.0) {
        G4cout << "[DetectorConstruction] Loading DEM from: " << fDEMPath << G4endl;
        G4cout << "[DetectorConstruction] Column size (m): "
               << "DX=" << fColDX_m << " DY=" << fColDY_m
               << " below z=0 by " << fColDZ_m << " m, "
               << "rotation phi=" << fTerrainPhiDeg << " deg."
               << G4endl;

        std::ifstream in(fDEMPath);
        if (!in) {
            G4cerr << "[DetectorConstruction] ERROR: cannot open DEM file: "
                   << fDEMPath << G4endl;
        } else {
            G4double colHalfX = 0.5 * fColDX_m * m;
            G4double colHalfY = 0.5 * fColDY_m * m;
            G4double bottom_m = -fColDZ_m; // bottom at z = -cdz

            G4double phiRad = fTerrainPhiDeg * deg;
            G4double cph    = std::cos(phiRad);
            G4double sph    = std::sin(phiRad);

            G4String solidName  = "DEMColumnSolid";
            G4String logicName  = "DEMColumnLV";
            G4String physName   = "DEMColumnPV";

            // Shared rock vis attributes (brown-ish solid)
            static G4VisAttributes* rockVis = nullptr;
            if (!rockVis) {
                rockVis = new G4VisAttributes(G4Colour(0.6, 0.4, 0.2));
                rockVis->SetForceSolid(true);
                rockVis->SetVisibility(true);
                rockVis->SetForceAuxEdgeVisible(true);
            }

            std::string line;
            G4int lineNo = 0;
            G4int placed = 0;

            while (std::getline(in, line)) {
                ++lineNo;
                if (line.empty()) continue;
                if (line[0] == '#') continue;

                std::istringstream iss(line);
                G4double x_m, y_m, z_m;
                if (!(iss >> x_m >> y_m >> z_m)) {
                    G4cerr << "[DetectorConstruction] WARNING: malformed DEM line "
                           << lineNo << " : \"" << line << "\"" << G4endl;
                    continue;
                }

                // Top and bottom in meters
                G4double top_m    = z_m;
                G4double bottom_m_local = bottom_m;

                if (top_m <= bottom_m_local) {
                    // Non-physical or zero-height column; skip
                    continue;
                }

                G4double centerZ_m = 0.5 * (top_m + bottom_m_local);
                G4double halfZ_m   = 0.5 * (top_m - bottom_m_local);

                G4double centerZ = centerZ_m * m;
                G4double halfZ   = halfZ_m   * m;

                // Rotated XY in WORLD coordinates (around origin, about +Z)
                G4double xr_m =  cph * x_m - sph * y_m;
                G4double yr_m =  sph * x_m + cph * y_m;

                G4double xr = xr_m * m;
                G4double yr = yr_m * m;

                // Create a logical volume for this column height
                auto* logicCol = new G4LogicalVolume(
                    new G4Box(solidName + std::to_string(lineNo),
                              colHalfX,
                              colHalfY,
                              halfZ),
                    fStandardRock,
                    logicName + std::to_string(lineNo)
                );

                logicCol->SetVisAttributes(rockVis);

                new G4PVPlacement(
                    nullptr,
                    G4ThreeVector(xr, yr, centerZ),
                    logicCol,
                    physName,
                    logicWorld,
                    false,
                    lineNo,
                    true
                );

                ++placed;
            }

            G4cout << "[DetectorConstruction] DEM: placed " << placed
                   << " rock columns." << G4endl;
        }
    } else {
        G4cout << "[DetectorConstruction] No DEM columns: "
               << "DEMpath=\"" << fDEMPath << "\", cdx=" << fColDX_m
               << ", cdy=" << fColDY_m << G4endl;
    }

    return physWorld;
}

// -----------------------------------------------------------------------------
// Sensitive detector
// -----------------------------------------------------------------------------
void DetectorConstruction::ConstructSDandField()
{
    if (!fPlaneLogical) return;

    auto* sdManager = G4SDManager::GetSDMpointer();

    auto* planeSD = new PlaneSD("PlaneSD");
    sdManager->AddNewDetector(planeSD);
    fPlaneLogical->SetSensitiveDetector(planeSD);
}
