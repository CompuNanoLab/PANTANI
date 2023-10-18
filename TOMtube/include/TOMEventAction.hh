#ifndef TOMEventAction_h
#define TOMEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class TOMEventAction : public G4UserEventAction
{
  public:
    TOMEventAction();
    virtual ~TOMEventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);
};
#endif


