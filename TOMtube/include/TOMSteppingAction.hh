#ifndef TOMSteppingAction_h
#define TOMSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include <variant>

class G4LogicalVolume;
class TOMEventAction;

class TOMSteppingAction : public G4UserSteppingAction
{
  public:
    TOMSteppingAction(TOMEventAction* eventAction);
    virtual ~TOMSteppingAction();
    virtual void UserSteppingAction(const G4Step* step);
};
#endif
