// SteppingAction.cc
#include "SteppingAction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    auto* track = step->GetTrack();
    if (track->GetParentID() > 0 && track->GetCurrentStepNumber() == 1) {
        track->SetTrackStatus(fStopAndKill);
    }
}
