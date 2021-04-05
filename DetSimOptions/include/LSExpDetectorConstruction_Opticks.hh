#pragma once

class G4Opticks ; 
class G4VPhysicalVolume ; 
class G4VSensitiveDetector ;

struct LSExpDetectorConstruction_Opticks
{
    static G4Opticks* Setup(const G4VPhysicalVolume* world, const G4VSensitiveDetector* sd_, int opticksMode ); 
}; 

