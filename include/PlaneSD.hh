#ifndef PLANE_SD_HH
#define PLANE_SD_HH

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class PlaneSD : public G4VSensitiveDetector
{
public:
    explicit PlaneSD(const G4String& name);
    ~PlaneSD() override;

    G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;
};

#endif
