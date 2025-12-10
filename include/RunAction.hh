#ifndef RUN_ACTION_HH
#define RUN_ACTION_HH

#include "G4UserRunAction.hh"
#include "G4String.hh"

class G4Run;

class RunAction : public G4UserRunAction
{
public:
    explicit RunAction(const G4String& outputFileName);
    ~RunAction() override;

    void BeginOfRunAction(const G4Run* run) override;
    void EndOfRunAction(const G4Run* run) override;
};

#endif
