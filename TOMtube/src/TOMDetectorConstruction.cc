#include "TOMDetectorConstruction.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"
#include "G4AssemblyVolume.hh"

TOMDetectorConstruction::TOMDetectorConstruction()
 : G4VUserDetectorConstruction()
{}

TOMDetectorConstruction::~TOMDetectorConstruction()
{}

G4VPhysicalVolume* TOMDetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  //Define materials
  //
  nist->FindOrBuildMaterial("G4_W");
  nist->FindOrBuildMaterial("G4_Al");
  nist->FindOrBuildMaterial("G4_Mo");
  nist->FindOrBuildMaterial("G4_Re");
  nist->FindOrBuildMaterial("G4_Ti");
  nist->FindOrBuildMaterial("G4_Zr");
  nist->FindOrBuildMaterial("G4_C");      
  //Define vacuum
  // 
  G4double a,z,density, fractionmass;
  G4int ncomponents;
  G4UnitDefinition::BuildUnitsTable();
  
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);  

  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = true;

  // World
  G4double world_sizeXYZ = 5*m;
  G4Material* world_mat = G4Material::GetMaterial("Galactic");

  G4Box* solidWorld =
    new G4Box("World",                       //its name
       0.5*world_sizeXYZ, 0.5*world_sizeXYZ, 0.5*world_sizeXYZ);     //its size

  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name

  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

  // Anode
  G4AssemblyVolume* assemblyAnode = new G4AssemblyVolume();
  G4ThreeVector Ta;
  G4Transform3D Tr;
  
  density = 19.3*g/cm3;
  G4Material* anode1_mat = new G4Material("Anode1", density, ncomponents=2);
  anode1_mat->AddMaterial(G4Material::GetMaterial("G4_Re"),fractionmass=5*perCent);
  anode1_mat->AddMaterial(G4Material::GetMaterial("G4_W"),fractionmass=95*perCent);
  G4double anode1_sizeX  = 1*mm;
  G4double anode1_sizeYZ = 14*mm;
  G4RotationMatrix* anode1_pRot = new G4RotationMatrix();
  anode1_pRot->rotateY(0.*deg);
  Ta.setX( 0. ); Ta.setY( 0. ); Ta.setZ( 0. );
  Tr = G4Transform3D(*anode1_pRot,Ta);
  G4VSolid* solidAnode1 
    = new G4Box("Anode1",        
                 0.5*anode1_sizeX, 0.5*anode1_sizeYZ, 0.5*anode1_sizeYZ);

  G4LogicalVolume* logicAnode1 =
    new G4LogicalVolume(solidAnode1,         //its solid
                        anode1_mat,          //its material
                        "Anode1");           //its name
                        
  assemblyAnode->AddPlacedVolume( logicAnode1, Tr );
  density = 10.28*g/cm3;                                     
  G4Material* anode2_mat = new G4Material("Anode2", density, ncomponents=4);
  anode2_mat->AddMaterial(G4Material::GetMaterial("G4_Mo"),fractionmass=99.4*perCent);
  anode2_mat->AddMaterial(G4Material::GetMaterial("G4_Ti"),fractionmass=0.5*perCent);
  anode2_mat->AddMaterial(G4Material::GetMaterial("G4_Zr"),fractionmass=0.08*perCent);
  anode2_mat->AddMaterial(G4Material::GetMaterial("G4_C"),fractionmass=0.02*perCent);    
  G4double anode2_sizeX  = 1*cm;
  G4double anode2_sizeYZ = 13*mm;
  G4RotationMatrix* anode2_pRot = new G4RotationMatrix();
  anode2_pRot->rotateY(0.*deg);
  Ta.setX( 0.5*(anode2_sizeX+anode1_sizeX) ); Ta.setY( 0. ); Ta.setZ( 0. );
  Tr = G4Transform3D(*anode2_pRot,Ta);
  G4VSolid* solidAnode2 
    = new G4Box("Anode2",        
                 0.5*anode2_sizeX, 0.5*anode2_sizeYZ, 0.5*anode2_sizeYZ);

  G4LogicalVolume* logicAnode2 =
    new G4LogicalVolume(solidAnode2,         //its solid
                        anode2_mat,          //its material
                        "Anode2");           //its name

  assemblyAnode->AddPlacedVolume( logicAnode2, Tr );
  G4RotationMatrix* anode_pRot = new G4RotationMatrix();
  anode_pRot->rotateY(-13.*deg);
  G4ThreeVector Tm(0,0,0);
  Tr = G4Transform3D(*anode_pRot,Tm);
  assemblyAnode->MakeImprint( logicWorld, Tr );    
                                  
  // Filter
  G4Material* filter_mat = G4Material::GetMaterial("G4_Al");
  G4double filter_sizeZ  = 2.4*mm;
  G4double filter_sizeXY = 50*cm;
  G4RotationMatrix* filter_pRot = new G4RotationMatrix();
  filter_pRot->rotateY(0.*deg);  

  G4VSolid* solidFilter 
    = new G4Box("Filter",        
                 0.5*filter_sizeXY, 0.5*filter_sizeXY, 0.5*filter_sizeZ);

  G4LogicalVolume* logicFilter =
    new G4LogicalVolume(solidFilter,         //its solid
                        filter_mat,          //its material
                        "Filter");           //its name

  new G4PVPlacement(filter_pRot,                       //no rotation
                    G4ThreeVector(0,0,-27*mm),    //at position
                    logicFilter,               //its logical volume
                    "Filter",                  //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
                   
  /*// Collimator
  G4Material* collim1_mat = G4Material::GetMaterial("G4_Al");
  G4double collim1_sizeZ  = 1*cm;
  G4double collim1_sizeXY = 100*cm;

  G4VSolid* solidCollim1 
    = new G4Box("Filter",        
                 0.5*collim1_sizeXY, 0.5*collim1_sizeXY, 0.5*collim1_sizeZ);

  G4LogicalVolume* logicCollim1 =
    new G4LogicalVolume(solidCollim1,         //its solid
                        collim1_mat,          //its material
                        "Collimator1");           //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(-0.5*(collim1_sizeXY+1*mm),0,-10*cm),    //at position
                    logicCollim1,               //its logical volume
                    "Collimator1",                  //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  
  G4Material* collim2_mat = G4Material::GetMaterial("G4_Al");
  G4double collim2_sizeZ  = 1*cm;
  G4double collim2_sizeXY = 100*cm;

  G4VSolid* solidCollim2 
    = new G4Box("Filter",        
                 0.5*collim2_sizeXY, 0.5*(collim2_sizeXY+1*mm), 0.5*collim2_sizeZ);

  G4LogicalVolume* logicCollim2 =
    new G4LogicalVolume(solidCollim2,         //its solid
                        collim2_mat,          //its material
                        "Collimator2");           //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0.5*collim2_sizeXY,0,-10*cm),    //at position
                    logicCollim2,               //its logical volume
                    "Collimator2",                  //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking*/
                    
  return physWorld;
}

