#include "TOMDetectorConstruction.hh"
#include "TOMActionInitialization.hh"
#include "G4RunManagerFactory.hh"
#include "G4SteppingVerbose.hh"
#include "G4UImanager.hh"
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "TOMPhysicsList.hh"
#include "Randomize.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"

int main(int argc,char** argv)
{
  G4String macro;
  #ifdef G4MULTITHREADED
    G4int nThreads = 0;
  #endif
  for (G4int i = 1; i < argc; i = i+2 ) 
  {
    if ( G4String(argv[i]) == "-m" ) macro = argv[i+1];
    #ifdef G4MULTITHREADED
      else if ( G4String(argv[i]) == "-t" ) 
      {
        nThreads = G4UIcommand::ConvertToInt(argv[i+1]);
      }
    #endif
    else 
    {
      G4cerr << " Usage: " << G4endl;
      G4cerr << " tom [-m macro ] [-u UIsession] [-t nThreads]" << G4endl;
      G4cerr << "   note: -t option is available only for multi-threaded mode."<< G4endl;
      return 1;
    }
  }
  // Detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = nullptr;
  if ( argc == 1 ) { ui = new G4UIExecutive(argc, argv); }

  G4int precision = 4;
  G4SteppingVerbose::UseBestUnit(precision);

  // Construct the default run manager
  #ifdef G4MULTITHREADED
    G4MTRunManager * runManager = new G4MTRunManager;
    if ( nThreads > 0 ) 
    { 
      runManager->SetNumberOfThreads(nThreads);
    }  
  #else
    G4RunManager * runManager = new G4RunManager;
  #endif

  // Set mandatory initialization classes
  runManager->SetUserInitialization(new TOMDetectorConstruction());
  G4VModularPhysicsList* physicsList = new TOMPhysicsList;
  physicsList->SetVerboseLevel(0);
  runManager->SetUserInitialization(physicsList);
  runManager->SetUserInitialization(new TOMActionInitialization());

  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  if ( ! ui ) 
  {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = macro;
    UImanager->ApplyCommand(command+fileName);
  }
  else 
  {
    // interactive mode
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  delete visManager;
  delete runManager;
  return EXIT_SUCCESS;
}
