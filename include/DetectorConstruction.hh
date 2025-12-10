#ifndef DETECTOR_CONSTRUCTION_HH
#define DETECTOR_CONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction(G4double halfX,
                         G4double halfY,
                         G4double halfZ,
                         G4double thetaDeg,
                         G4double detXY_m,
                         const G4String& demPath,
                         G4double colDX_m,
                         G4double colDY_m,
                         G4double colDZ_m,
                         G4double terrainPhiDeg,
                         G4double rockDensity_g_cm3);
    ~DetectorConstruction() override;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    // Plane basis in WORLD coordinates (local ex, ey, ez as seen in world)
    inline void GetPlaneBasis(G4ThreeVector& ex,
                              G4ThreeVector& ey,
                              G4ThreeVector& ez) const
    {
        ex = fPlaneEx;
        ey = fPlaneEy;
        ez = fPlaneEz;
    }

private:
    // Materials
    void DefineMaterials();
    G4Material* fStandardRock;

    // World half-sizes (internal units)
    G4double fHalfX;
    G4double fHalfY;
    G4double fHalfZ;

    // Plane orientation (zenith angle, degrees)
    G4double fThetaDeg;

    // Plane basis vectors in WORLD coordinates:
    G4ThreeVector fPlaneEx;
    G4ThreeVector fPlaneEy;
    G4ThreeVector fPlaneEz;

    // Sensitive plane logical volume
    G4LogicalVolume* fPlaneLogical;

    // Plane size control (full side length in meters; 0 => auto)
    G4double fDetXY_m;

    // DEM / terrain parameters
    G4String fDEMPath;
    G4double fColDX_m;        // full size in X (meters)
    G4double fColDY_m;        // full size in Y (meters)
    G4double fColDZ_m;        // extension below z=0 (meters)
    G4double fTerrainPhiDeg;  // rotation of DEM XY around world Z (degrees)
    G4double fRockDensity_g_cm3; // density for StandardRock (g/cm3)
};

#endif
