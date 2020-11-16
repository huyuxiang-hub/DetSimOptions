#include "ExpHallConstruction.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/ToolFactory.h"

#include "DetSimAlg/DetSimAlg.h"
#include "DetSimAlg/IDetElement.h"

DECLARE_TOOL(ExpHallConstruction);

ExpHallConstruction::ExpHallConstruction(const std::string& name)
    : ToolBase(name)
{
    logicExpHall = 0;
}

ExpHallConstruction::~ExpHallConstruction() {

}

G4LogicalVolume*
ExpHallConstruction::getLV() {
    if (logicExpHall) {
        return logicExpHall;
    }
    initVariables();
    initMaterials();

    makeExpHallLogical();

    return logicExpHall;
}

bool
ExpHallConstruction::inject(std::string motherName, IDetElement* other, IDetElementPos* pos) {
    G4LogicalVolume* mothervol = 0;
    G4LogicalVolume* daughtervol = 0;

    if (motherName == "lExpHall") {
        mothervol = logicExpHall;
    } else {
        return false;
    }

    if (not other) {
        return false;
    }
    daughtervol = other->getLV();
int copyno=0;

    // Place the Top Tracker into the Exp. Hall 
    if (pos == 0) {
      //G4RotationMatrix MRot; 
      //MRot.rotateZ(M_PI/2.*rad);
      //new G4PVPlacement(G4Transform3D(MRot,G4ThreeVector(0,0, m_TopTrackerOffsetZ)),
      new G4PVPlacement(0,G4ThreeVector(0,0, m_TopTrackerOffsetZ),
                daughtervol,
                "pTopTracker",
                mothervol,
                false,
                0
                );
    }
    else {
    
       while (pos->hasNext()) {
        new G4PVPlacement(
            pos->next(),
            daughtervol,
            daughtervol->GetName()+"_phys",
            mothervol,
            false,
            copyno++
    //        true
                );
    } 
    
    }
    return true;
}

void
ExpHallConstruction::makeExpHallLogical() {

    solidExpHall = new G4Box("sExpHall", m_expHallX, m_expHallY, m_expHallZ);

    logicExpHall = new G4LogicalVolume(solidExpHall,
            Air,
            "lExpHall",
            0,0,0
            );
}
void
ExpHallConstruction::makePoolCoverLogical() {
 // pool cover
   SniperPtr<DetSimAlg> detsimalg(*getParent(), "DetSimAlg");
    if (detsimalg.invalid()) {
      std::cout << "Can't Load DetSimAlg" << std::endl;
      assert(0);
    }
    ToolBase* t = detsimalg->findTool("GlobalGeomInfo");
    assert(t);
    //IDetElement* globalinfo = dynamic_cast<IDetElement*>(t);

   G4double fR   = m_pooR +0.05*m;
   G4double CHradius = 0.75*m;
   G4VSolid*PoolCover   =  new G4Tubs("coversolid",0.*m,fR ,m_poolcoverHeight ,0.*deg, 360.*deg);
   G4VSolid*CentralHole =  new G4Tubs("CentralHole1",0,CHradius,m_poolcoverHeight+2*mm,0*deg,360*deg);
   PoolCover =  new G4SubtractionSolid("PoolCoversub",PoolCover,CentralHole,0, G4ThreeVector(0*m,0*m,0*m));
   logicPoolCover = new G4LogicalVolume(PoolCover, Steel, "lPoolCover",0,0,0);
 }

void
ExpHallConstruction::makePoolCoverPhysical() {
 // G4VPhysicalVolume* physPoolCover = 
      new G4PVPlacement(0, G4ThreeVector(0, 0, -m_expHallZ+m_poolcoverHeight/2),
          logicPoolCover,"pPoolCover",logicExpHall , false, 0);
}

void
ExpHallConstruction::initMaterials() {
    Air   = G4Material::GetMaterial("Air");
    Steel = G4Material::GetMaterial("Steel");
}

void 
ExpHallConstruction::initVariables() {
    SniperPtr<DetSimAlg> detsimalg(*getParent(), "DetSimAlg");
    if (detsimalg.invalid()) {
      std::cout << "Can't Load DetSimAlg" << std::endl;
      assert(0);
    }
    ToolBase* t = detsimalg->findTool("GlobalGeomInfo");
    assert(t);
    IDetElement* globalinfo = dynamic_cast<IDetElement*>(t);

    // the original dim is 24m. 
    m_expHallX = globalinfo->geom_info("ExpHall.+X");// 24*m;//Half length
    m_expHallY = globalinfo->geom_info("ExpHall.+Y");// 24*m;
    m_expHallZ = globalinfo->geom_info("ExpHall.+Z");// 9.3*m;
    m_pooR     = globalinfo->geom_info("WaterPool.R");
    double m_TTVerticalSize = 8.40*m;
    m_poolcoverHeight = 2*mm;
    m_TopTrackerOffsetZ = -(m_expHallZ - m_TTVerticalSize/2.)+m_poolcoverHeight;

}
