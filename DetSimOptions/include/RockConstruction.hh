#ifndef RockConstruction_hh
#define RockConstruction_hh

/*  
 * In current implmentation, the Rock is divided into two parts.
 *
 * * Top Rock (Box)
 *   + Exp. Hall
 *
 * * Bottom Rock (Tube)
 *   + Water Pool Lining 
 *     
 *     + Central Detector
 */

#include "SniperKernel/ToolBase.h"
#include "DetSimAlg/IDetElement.h"
#include "globals.hh"

class G4Box;
class G4Sphere;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4VSolid;

class TopRockConstruction : public IDetElement,
                            public ToolBase {
public:
    G4LogicalVolume* getLV();
    bool inject(std::string motherName, IDetElement* other, IDetElementPos* pos);

    double geom_info(const std::string& param);
    double geom_info(const std::string& , int) { return 0.0; }

    TopRockConstruction(const std::string& name);
    ~TopRockConstruction();

private:

    void initMaterials();
    void initVariables();

    void makeTopRockLogical();

    G4VSolid*        solidTopRock;
    G4LogicalVolume* logicTopRock;

    G4Material* Rock;

    double m_topRockX;
    double m_topRockY;
    double m_topRockZ;

    double m_offset;
    double m_offset_in_world;
};

class BottomRockConstruction : public IDetElement,
                               public ToolBase {
public:
    G4LogicalVolume* getLV();
    bool inject(std::string motherName, IDetElement* other, IDetElementPos* pos);

    BottomRockConstruction(const std::string& name);
    ~BottomRockConstruction();

private:

    void initMaterials();
    void initVariables();

    void makeBottomRockLogical();
    void makePoolLiningLogical();
    void makePoolLiningPhysical();

    G4VSolid*        solidBottomRock;
    G4LogicalVolume* logicBottomRock;

    G4VSolid*        solidPoolLining;
    G4LogicalVolume* logicPoolLining;
    G4VPhysicalVolume* physiPoolLining;

    G4Material* Rock;
    G4Material* Tyvek;

    double m_btmRockR;
    double m_btmRockH;
    double m_btmRockZtop;
    double m_btmRockZbottom;
    double m_poolLiningR;
    double m_poolLiningH;
    double m_poolLiningZtop;
    double m_poolLiningZbottom;
    double m_waterPoolRadius;
};

#endif
