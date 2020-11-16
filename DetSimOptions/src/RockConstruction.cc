#include "RockConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Polycone.hh"
#include "G4VisAttributes.hh"

#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/ToolFactory.h"
#include "SniperKernel/SniperLog.h"

#include "DetSimAlg/DetSimAlg.h"
#include "DetSimAlg/IDetElement.h"

DECLARE_TOOL(TopRockConstruction);
DECLARE_TOOL(BottomRockConstruction);

// = Top Rock =
TopRockConstruction::TopRockConstruction(const std::string& name)
    : ToolBase(name)
{
    logicTopRock = 0;
}

TopRockConstruction::~TopRockConstruction() {

}

G4LogicalVolume*
TopRockConstruction::getLV() {
    if (logicTopRock) {
        return logicTopRock;
    }
    initMaterials();
    initVariables();
    makeTopRockLogical();
    return logicTopRock;
}

bool
TopRockConstruction::inject(std::string motherName, IDetElement* other, IDetElementPos* pos) {
    G4LogicalVolume* mothervol = 0;
    G4LogicalVolume* daughtervol = 0;

    if (motherName == "lTopRock") {
        mothervol = logicTopRock;
    } else {
        LogWarn << "Can't find < " 
                << motherName
                << "> in TopRock."
                << std::endl;
        return false;
    }

    daughtervol = other->getLV();
    LogDebug << "Exp Hall. LV: " << daughtervol << std::endl;
    if (not daughtervol) {
        LogWarn << "Can't find any injected LV in TopRock" << std::endl;
        return false;
    }

    // Place the Exp.Hall into the top rock
    if (pos == 0 and daughtervol->GetName() == "lExpHall") {
        LogDebug << "Create pExpHall Begin:" << std::endl;
        LogDebug << "pExpHall offset : " << m_offset << std::endl;
        LogDebug << "lExpHall: " << daughtervol << std::endl;
        LogDebug << "lTopRock: " << mothervol << std::endl;
        new G4PVPlacement(0, // no rotation
                G4ThreeVector(0, 0, m_offset), // offset
                daughtervol,
                "pExpHall",
                mothervol,
                false,
                0
                );
        LogDebug << "Create pExpHall End." << std::endl;
    }

    return true;
}

double
TopRockConstruction::geom_info(const std::string& param) {
    // make sure the geom is setup
    assert (getLV());
    if (param == "TopRockOffset.Z") {
        return m_offset_in_world;
    }
    throw SniperException(std::string("Can't find param: ") + param);
    return 0;

}

void
TopRockConstruction::initMaterials() {
    Rock = G4Material::GetMaterial("Rock");
}

void
TopRockConstruction::initVariables() {
    SniperPtr<DetSimAlg> detsimalg(*getParent(), "DetSimAlg");
    if (detsimalg.invalid()) {
      std::cout << "Can't Load DetSimAlg" << std::endl;
      assert(0);
    }
    ToolBase* t = detsimalg->findTool("GlobalGeomInfo");
    assert(t);
    IDetElement* globalinfo = dynamic_cast<IDetElement*>(t);

    double m_expHallX = globalinfo->geom_info("ExpHall.+X");//Half length
    double m_expHallY = globalinfo->geom_info("ExpHall.+Y");//Half length
    double m_expHallZ = globalinfo->geom_info("ExpHall.+Z");//Half length
    m_topRockX = m_expHallX + 3*m;//Half length
    m_topRockY = m_expHallY + 3*m;
    m_topRockZ = m_expHallZ + 1.5*m;
    
    m_offset = m_expHallZ - m_topRockZ;
    m_offset_in_world = m_topRockZ + globalinfo->geom_info("WaterPool.R");
}

void
TopRockConstruction::makeTopRockLogical() {
    solidTopRock = new G4Box("sTopRock", m_topRockX, m_topRockY, m_topRockZ);
    logicTopRock = new G4LogicalVolume(solidTopRock,
        Rock,
        "lTopRock",
        0,0,0);

}

// = Bottom Rock =
BottomRockConstruction::BottomRockConstruction(const std::string& name)
    : ToolBase(name)
{
    logicBottomRock = 0;
}

BottomRockConstruction::~BottomRockConstruction() {

}

G4LogicalVolume*
BottomRockConstruction::getLV() {
    if (logicBottomRock) {
        return logicBottomRock;
    }
    initMaterials();
    initVariables();
    makeBottomRockLogical();

    makePoolLiningLogical();
    makePoolLiningPhysical();
    return logicBottomRock;
}

bool
BottomRockConstruction::inject(std::string motherName, IDetElement* other, IDetElementPos* pos) {
    G4LogicalVolume* mothervol = 0;
    G4LogicalVolume* daughtervol = 0;

    if (motherName == "lPoolLining") {
        mothervol = logicPoolLining;
    } else {
        return false;
    }

    if (not other) {
        return false;
    }
    daughtervol = other->getLV();
    if (not daughtervol) {
        return false;
    }

    // Place the Water Pool into the bottom rock
    if (pos == 0) {
        new G4PVPlacement(0, // no rotation
                G4ThreeVector(0, 0, 0), // offset
                daughtervol,
                "pOuterWaterPool",
                mothervol,
                false,
                0
                );
    }

    return true;
}

void
BottomRockConstruction::initMaterials() {
    Rock = G4Material::GetMaterial("Rock");
    Tyvek = G4Material::GetMaterial("Tyvek");
}

void
BottomRockConstruction::initVariables() {
    SniperPtr<DetSimAlg> detsimalg(*getParent(), "DetSimAlg");
    if (detsimalg.invalid()) {
      std::cout << "Can't Load DetSimAlg" << std::endl;
      assert(0);
    }
    ToolBase* t = detsimalg->findTool("GlobalGeomInfo");
    assert(t);
    IDetElement* globalinfo = dynamic_cast<IDetElement*>(t);

    // TODO: load these data from geom svc
    m_waterPoolRadius = globalinfo->geom_info("WaterPool.R");
    m_btmRockH = globalinfo->geom_info("WaterPool.H") + 3*m;
    m_btmRockR = m_waterPoolRadius + 3*m;
    m_btmRockZtop = m_waterPoolRadius;
    m_btmRockZbottom = m_waterPoolRadius - m_btmRockH;

    m_poolLiningH = globalinfo->geom_info("WaterPool.H") + 3*mm;
    m_poolLiningZtop = m_waterPoolRadius;
    m_poolLiningZbottom = m_waterPoolRadius - m_poolLiningH;
    m_poolLiningR = m_waterPoolRadius + 3*mm;//3*mm thickness steel tube out around water pool
}

void
BottomRockConstruction::makeBottomRockLogical() {
    // solidBottomRock = new G4Tubs("sBottomRock", 0, m_btmRockR, m_btmRockZ, 0, 360*deg);
    // please note, increase order
    static G4double  zPlane[] = {m_btmRockZbottom, m_btmRockZtop};
    static G4double  rInner[] = {0., 0.};
    static G4double  rOuter[] = {m_btmRockR, m_btmRockR};
    solidBottomRock = new G4Polycone("sBottomRock", //const G4String& pName,
                                     0, //G4double  phiStart,
                                     360*deg, //G4double  phiTotal,
                                     2, //G4int         numZPlanes,
                                     zPlane, //const G4double  zPlane[],
                                     rInner, //const G4double  rInner[],
                                     rOuter //const G4double  rOuter[])
                                    );
    logicBottomRock = new G4LogicalVolume(solidBottomRock,
        Rock,
        "lBtmRock",
        0, 0, 0);

}

void
BottomRockConstruction::makePoolLiningLogical()
{
    // solidPoolLining = new G4Tubs("sPoolLining", 0, m_poolLiningR, m_poolLiningZ, 0, 360*deg);
    static G4double  zPlane[] = {m_poolLiningZbottom, m_poolLiningZtop};
    static G4double  rInner[] = {0., 0.};
    static G4double  rOuter[] = {m_poolLiningR, m_poolLiningR};
    solidPoolLining = new G4Polycone("sPoolLining", //const G4String& pName,
                                     0, //G4double  phiStart,
                                     360*deg, //G4double  phiTotal,
                                     2, //G4int         numZPlanes,
                                     zPlane, //const G4double  zPlane[],
                                     rInner, //const G4double  rInner[],
                                     rOuter //const G4double  rOuter[])
                                    );
    logicPoolLining = new G4LogicalVolume(solidPoolLining,
        Tyvek,
        "lPoolLining",
        0, 0, 0);
}

void
BottomRockConstruction::makePoolLiningPhysical()
{
    physiPoolLining = new G4PVPlacement(0,              // no rotation
        G4ThreeVector(0,0,0), // at (x,y,z)
        logicPoolLining,    // its logical volume
        "pPoolLining",       // its name
        logicBottomRock,      // its mother  volume
        false,           // no boolean operations
        0);              // no particular field
}
