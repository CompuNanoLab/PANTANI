#ifndef TOMDetectorConstruction_h
#define TOMDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

class TOMDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    TOMDetectorConstruction();
    ~TOMDetectorConstruction() override;

    G4VPhysicalVolume* Construct() override;

    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

  protected:
    G4LogicalVolume* fScoringVolume = nullptr;
};
#endif
