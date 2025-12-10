#include "PlaneSD.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4TouchableHistory.hh"
#include "G4AffineTransform.hh"
#include "G4AnalysisManager.hh"

#include <cmath>

PlaneSD::PlaneSD(const G4String& name)
: G4VSensitiveDetector(name)
{}

PlaneSD::~PlaneSD() = default;

G4bool PlaneSD::ProcessHits(G4Step* step, G4TouchableHistory* /*history*/)
{
    G4StepPoint* pre  = step->GetPreStepPoint();
    G4StepPoint* post = step->GetPostStepPoint();

    // Only first step when entering this volume
    if (pre->GetStepStatus() != fGeomBoundary)
        return false;

    G4Track* track = step->GetTrack();
    if (!track) return false;

    G4ParticleDefinition* particle = track->GetDefinition();
    if (!particle) return false;

    // Only muons (mu- or mu+)
    G4int pdg = particle->GetPDGEncoding();
    if (pdg != 13 && pdg != -13)
        return false;

    // World coordinates at plane entry/exit
    G4ThreeVector posInWorld  = pre->GetPosition();
    G4ThreeVector posOutWorld = post->GetPosition();
    G4ThreeVector dirInWorld  = pre->GetMomentumDirection();
    G4ThreeVector dirOutWorld = post->GetMomentumDirection();

    // World coordinates at SOURCE (vertex)
    G4ThreeVector vertexWorld   = track->GetVertexPosition();
    G4ThreeVector vertexDirWorld = track->GetVertexMomentumDirection();

    // Plane local frame via touchable
    G4TouchableHandle touchable = pre->GetTouchableHandle();
    const G4AffineTransform& topTransform = touchable->GetHistory()->GetTopTransform();

    // Transform positions/directions to plane-local coordinates
    G4ThreeVector posInLocal    = topTransform.TransformPoint(posInWorld);
    G4ThreeVector posOutLocal   = topTransform.TransformPoint(posOutWorld);
    G4ThreeVector vertexLocal   = topTransform.TransformPoint(vertexWorld);
    G4ThreeVector dirInLocal    = topTransform.TransformAxis(dirInWorld);
    G4ThreeVector dirOutLocal   = topTransform.TransformAxis(dirOutWorld);
    G4ThreeVector dirInitLocal  = topTransform.TransformAxis(vertexDirWorld);

    // ---- Compute intersection with plane z' = 0 in local coords ----
    //
    // segment: P(t) = P_in + t (P_out - P_in), t in [0,1]
    // want P(t_plane) with z' = 0
    //   z_in + t (z_out - z_in) = 0  =>  t = -z_in / (z_out - z_in)
    //
    G4double z1 = posInLocal.z();
    G4double z2 = posOutLocal.z();

    G4ThreeVector posPlaneLocal;
    G4double      tPlane = 0.0;

    if (std::fabs(z2 - z1) < 1e-12) {
        // Degenerate: almost no thickness or step parallel to plane.
        // Just use entry position as best guess.
        posPlaneLocal = posInLocal;
    } else {
        tPlane = -z1 / (z2 - z1);
        // Clamp t to [0,1] just in case
        if (tPlane < 0.0) tPlane = 0.0;
        if (tPlane > 1.0) tPlane = 1.0;
        posPlaneLocal = posInLocal + tPlane * (posOutLocal - posInLocal);
    }

    // Direction at plane center (linear blend between entry/exit)
    G4ThreeVector dirPlaneLocal = dirInLocal + tPlane * (dirOutLocal - dirInLocal);
    if (dirPlaneLocal.mag2() == 0.) {
        dirPlaneLocal = dirOutLocal;
    }
    dirPlaneLocal = dirPlaneLocal.unit();

    // Energies in GeV
    G4double E_in_GeV    = pre->GetKineticEnergy()  / GeV;
    G4double E_out_GeV   = post->GetKineticEnergy() / GeV;
    G4double E_init_GeV  = track->GetVertexKineticEnergy() / GeV;
    G4double E_plane_GeV = E_in_GeV + tPlane * (E_out_GeV - E_in_GeV);

    // Event / track IDs
    G4int eventID = -1;
    if (auto* runManager = G4RunManager::GetRunManager()) {
        if (auto* evt = runManager->GetCurrentEvent()) {
            eventID = evt->GetEventID();
        }
    }
    G4int trackID = track->GetTrackID();

    // Zenith angle at firing (world frame, from +Z)
    G4ThreeVector vdir0 = track->GetVertexMomentumDirection();  // unit vector
    G4double cosZen = vdir0.z();
    if (cosZen >  1.0) cosZen =  1.0;
    if (cosZen < -1.0) cosZen = -1.0;
    G4double zenithRad = std::acos(cosZen);
    G4double zenithDeg = zenithRad / deg;

    // ---- Fill ntuple ----
    auto* analysis = G4AnalysisManager::Instance();
    if (!analysis) return false;

    // 0,1: IDs
    analysis->FillNtupleIColumn(0, eventID);
    analysis->FillNtupleIColumn(1, trackID);

    // 2-4: SOURCE point in plane-local coordinates (entry = source)
    analysis->FillNtupleDColumn(2, vertexLocal.x());
    analysis->FillNtupleDColumn(3, vertexLocal.y());
    analysis->FillNtupleDColumn(4, vertexLocal.z());

    // 5-7: intersection with plane z' = 0 (local)
    analysis->FillNtupleDColumn(5, posPlaneLocal.x());
    analysis->FillNtupleDColumn(6, posPlaneLocal.y());
    analysis->FillNtupleDColumn(7, posPlaneLocal.z());

    // 8-10: initial direction at source (local)
    analysis->FillNtupleDColumn(8,  dirInitLocal.x());
    analysis->FillNtupleDColumn(9,  dirInitLocal.y());
    analysis->FillNtupleDColumn(10, dirInitLocal.z());

    // 11-13: direction at plane center (local)
    analysis->FillNtupleDColumn(11, dirPlaneLocal.x());
    analysis->FillNtupleDColumn(12, dirPlaneLocal.y());
    analysis->FillNtupleDColumn(13, dirPlaneLocal.z());

    // 14-15: energies
    analysis->FillNtupleDColumn(14, E_init_GeV);
    analysis->FillNtupleDColumn(15, E_plane_GeV);

    // 16: zenith of the *initial* direction (world)
    analysis->FillNtupleDColumn(16, zenithDeg);

    analysis->AddNtupleRow();

    return true;
}
