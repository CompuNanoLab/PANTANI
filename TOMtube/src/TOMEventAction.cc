#include "TOMEventAction.hh"
#include "TOMRunAction.hh"
#include "G4MTRunManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "Randomize.hh"
#include <iomanip>
#include "G4Event.hh"
#include "G4RunManager.hh"

TOMEventAction::TOMEventAction()
 : G4UserEventAction()
{}

TOMEventAction::~TOMEventAction()
{}

void TOMEventAction::BeginOfEventAction(const G4Event*)
{}

void TOMEventAction::EndOfEventAction(const G4Event*)
{}

