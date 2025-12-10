#ifndef DETECTOR_CONSTRUCTION_HH
#define DETECTOR_CONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    // halfX/Y/Z in internal units (e.g. m)
    // thetaDeg, phiDeg in degrees (zenith, azimuth)
    DetectorConstruction(G4double halfX,
                         G4double halfY,
                         G4double halfZ,
                         G4double thetaDeg,
                         G4double phiDeg);
    ~DetectorConstruction() override;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

private:
    G4double         fHalfX;
    G4double         fHalfY;
    G4double         fHalfZ;
    G4double         fThetaDeg;
    G4double         fPhiDeg;
    G4LogicalVolume* fPlaneLogical;
};

#endif
