#include "RunAction.hh"

#include "G4Run.hh"
#include "G4AnalysisManager.hh"
#include "G4String.hh"

RunAction::RunAction(const G4String& outputFileName)
: G4UserRunAction()
{
    auto* analysisManager = G4AnalysisManager::Instance();

    analysisManager->SetDefaultFileType("root");
    analysisManager->SetFileName(outputFileName);
    analysisManager->SetVerboseLevel(1);
    analysisManager->SetNtupleMerging(true);  // MT friendly

    // Define ntuple structure
    analysisManager->CreateNtuple("Plane", "Plane crossings");

    analysisManager->CreateNtupleIColumn("eventID"); // 0
    analysisManager->CreateNtupleIColumn("trackID"); // 1

    // Positions (local to plane)
    analysisManager->CreateNtupleDColumn("x_init");     // 2
    analysisManager->CreateNtupleDColumn("y_init");     // 3
    analysisManager->CreateNtupleDColumn("z_init");     // 4

    analysisManager->CreateNtupleDColumn("x_plane");    // 5
    analysisManager->CreateNtupleDColumn("y_plane");    // 6
    analysisManager->CreateNtupleDColumn("z_plane");    // 7

    // Directions (local to plane)
    analysisManager->CreateNtupleDColumn("dx_init");    // 8
    analysisManager->CreateNtupleDColumn("dy_init");    // 9
    analysisManager->CreateNtupleDColumn("dz_init");    // 10

    analysisManager->CreateNtupleDColumn("dx_plane");   // 11
    analysisManager->CreateNtupleDColumn("dy_plane");   // 12
    analysisManager->CreateNtupleDColumn("dz_plane");   // 13

    // Energies (GeV)
    analysisManager->CreateNtupleDColumn("E_init_GeV");   // 14
    analysisManager->CreateNtupleDColumn("E_plane_GeV");  // 15

    // NEW: zenith angle at firing (world frame, degrees)
    analysisManager->CreateNtupleDColumn("zenith_world_deg"); // 16

    analysisManager->FinishNtuple();
}

RunAction::~RunAction()
{
    delete G4AnalysisManager::Instance();
}

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
    auto* analysisManager = G4AnalysisManager::Instance();
    analysisManager->OpenFile();
}

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
    auto* analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();
}
