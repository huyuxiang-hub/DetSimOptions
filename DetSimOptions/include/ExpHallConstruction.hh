#ifndef ExpHallConstruction_hh
#define ExpHallConstruction_hh

#include "SniperKernel/ToolBase.h"
#include "DetSimAlg/IDetElement.h"
#include "globals.hh"

class G4Box;
class G4Sphere;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;

class ExpHallConstruction : public IDetElement,
                            public ToolBase{
public:
    G4LogicalVolume* getLV();
    bool inject(std::string motherName, IDetElement* other, IDetElementPos* pos);

    ExpHallConstruction(const std::string& name);
    ~ExpHallConstruction();
    double m_poolcoverHeight;
    double m_pooR;
private:

    void initMaterials();
    void initVariables();

    void makeExpHallLogical();
    void makePoolCoverLogical();
    void makePoolCoverPhysical();
private:

    double m_expHallX;
    double m_expHallY;
    double m_expHallZ;

    double m_TopTrackerOffsetZ;

    G4Box*           solidExpHall;
    G4LogicalVolume* logicExpHall;
    G4Tubs*           solidPoolCover;
    G4LogicalVolume* logicPoolCover;

 
    G4Material*      Air;
    G4Material*      Steel;
};

#endif
