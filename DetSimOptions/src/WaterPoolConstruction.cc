
#include "WaterPoolConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "RoundBottomFlaskSolidMaker.hh"
#include "G4Transform3D.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/ToolFactory.h"
#include "SniperKernel/SniperLog.h"

#include "DetSimAlg/DetSimAlg.h"
#include "DetSimAlg/IDetElement.h"

DECLARE_TOOL(WaterPoolConstruction);

WaterPoolConstruction::WaterPoolConstruction(const std::string& name)
    : ToolBase(name)
{
    logicWaterPool = 0;

    declProp("enabledLatticedShell", m_enabledLatticedShell=false);
}

WaterPoolConstruction::~WaterPoolConstruction() {

}

G4LogicalVolume* 
WaterPoolConstruction::getLV() {
    if (logicWaterPool) {
        return logicWaterPool;
    }
    initVariables();
    initMaterials();

    makeWPLogical();

    if (m_enabledLatticedShell) { makeLatticedShell(); }

    return logicWaterPool;
}

bool
WaterPoolConstruction::inject(std::string motherName, IDetElement* other, IDetElementPos* pos) {
    // = Outer Water Pool =
    // == CD ==
    // == Veto PMTs ==
    G4LogicalVolume* mothervol = 0;
    G4int copyno = 0;

    if (motherName == "lWaterPool") {
        // place CD
        mothervol = logicWaterPool;
    } else if (motherName == "OuterWaterVeto") {
        // place Veto PMTs
        mothervol = logicWaterPool;
        copyno = 30000;
    }
    if (not mothervol) {
        // don't find the volume.
        LogError << "Can't find the mother volume named " << motherName << std::endl;
        return false;
    }

    // retrieve the daughter's LV
    G4LogicalVolume* daughtervol = other->getLV();

    if (not daughtervol) {
        LogError << "Can't find the daughter volume in " << motherName << std::endl;
        return false;
    }

    // Place the Central Detector (with water) into the Water Pool
    if (pos == 0) {
        new G4PVPlacement(0, // no rotation
                G4ThreeVector(0, 0, 0), // offset
                daughtervol,
                "pCentralDetector",
                mothervol,
                false,
                0
                );
    } else {
        while (pos->hasNext()) {
            new G4PVPlacement(
                pos->next(),
                daughtervol,
                daughtervol->GetName()+"_phys",
                mothervol,
                false,
                copyno++
                    );
        }
    }

    if(motherName == "OuterWaterVeto") {
        LogInfo << "Veto PMT Number = " << copyno-30000 << std::endl;
    }

    return true;
}

void
WaterPoolConstruction::initVariables() {
    SniperPtr<DetSimAlg> detsimalg(*getParent(), "DetSimAlg");
    if (detsimalg.invalid()) {
      std::cout << "Can't Load DetSimAlg" << std::endl;
      assert(0);
    }
    ToolBase* t = detsimalg->findTool("GlobalGeomInfo");
    assert(t);
    IDetElement* globalinfo = dynamic_cast<IDetElement*>(t);

    m_radWP = globalinfo->geom_info("WaterPool.R"); // This value will be changed in the future
    m_heightWP = globalinfo->geom_info("WaterPool.H");

}

void 
WaterPoolConstruction::initMaterials() {
    vetoWater = G4Material::GetMaterial("vetoWater");
}

void
WaterPoolConstruction::makeWPLogical() {
    solidWaterPool = new G4Tubs("sOuterWaterPool",
                                0*m,
                                m_radWP,
                                m_heightWP/2,
                                0.*deg,
                                360*deg);

    logicWaterPool = new G4LogicalVolume(solidWaterPool,
                                         vetoWater,
                                         "lOuterWaterPool",
                                         0,
                                         0,
                                         0);
}

struct ShellData {
    std::string name; // model.label
    double theta;     // theta
};

ShellData all_data_layer[] = {
    {"GLb3.up11", 90.-8.4-8.4-8.4-8.4-8.2-8.0-7.8-7.6-7.0-6.2-6.0},
    {"GLb4.up10", 90.-8.4-8.4-8.4-8.4-8.2-8.0-7.8-7.6-7.0-6.2},
    {"GLb3.up09", 90.-8.4-8.4-8.4-8.4-8.2-8.0-7.8-7.6-7.0},
    {"GLb2.up08", 90.-8.4-8.4-8.4-8.4-8.2-8.0-7.8-7.6},
    {"GLb2.up07", 90.-8.4-8.4-8.4-8.4-8.2-8.0-7.8},
    {"GLb2.up06", 90.-8.4-8.4-8.4-8.4-8.2-8.0},
    {"GLb1.up05", 90.-8.4-8.4-8.4-8.4-8.2},
    {"GLb1.up04", 90.-8.4-8.4-8.4-8.4},
    {"GLb1.up03", 90.-8.4-8.4-8.4},
    {"GLb1.up02", 90.-8.4-8.4},
    {"GLb1.up01", 90.-8.4},
    {"GLb2.equ", 90}, // equator
    {"GLb2.bt01", 90.+8.4}, 
    {"GLb1.bt02", 90.+8.4+8.4}, 
    {"GLb2.bt03", 90.+8.4+8.4+8.4}, 
    {"GLb2.bt04", 90.+8.4+8.4+8.4+8.4}, // FIXME
    {"GLb1.bt05", 90.+8.4+8.4+8.4+8.4+8.2}, 
    {"GLb1.bt06", 90.+8.4+8.4+8.4+8.4+8.2+8.0}, 
    {"GLb1.bt07", 90.+8.4+8.4+8.4+8.4+8.2+8.0+7.8}, 
    {"GLb1.bt08", 90.+8.4+8.4+8.4+8.4+8.2+8.0+7.8+7.6}, 
    {"GLb3.bt09", 90.+8.4+8.4+8.4+8.4+8.2+8.0+7.8+7.6+7.0}, 
    {"GLb3.bt10", 90.+8.4+8.4+8.4+8.4+8.2+8.0+7.8+7.6+7.0+6.2}, 
    {"GLb3.bt11", 90.+8.4+8.4+8.4+8.4+8.2+8.0+7.8+7.6+7.0+6.2+6.0}, 
    {"END", 180}
};

ShellData all_data_col[] = {
    {"GLw1.up10_up11", 90.-8.4-8.4-8.4-8.4-8.2-8.0-7.8-7.6-7.0-6.2-6.0/2},
    {"GLw1.up09_up10", 90.-8.4-8.4-8.4-8.4-8.2-8.0-7.8-7.6-7.0-6.2/2},
    {"GLw1.up08_up09", 90.-8.4-8.4-8.4-8.4-8.2-8.0-7.8-7.6-7.0/2},
    {"GLw1.up07_up08", 90.-8.4-8.4-8.4-8.4-8.2-8.0-7.8-7.6/2},
    {"GLw1.up06_up07", 90.-8.4-8.4-8.4-8.4-8.2-8.0-7.8/2},
    {"GLw1.up05_up06", 90.-8.4-8.4-8.4-8.4-8.2-8.0/2},
    {"GLw1.up04_up05", 90.-8.4-8.4-8.4-8.4-8.2/2},
    {"GLw1.up03_up04", 90.-8.4-8.4-8.4-8.4/2},
    {"GLw1.up02_up03", 90.-8.4-8.4-8.4/2},
    {"GLw1.up01_up02", 90.-8.4-8.4/2},
    {"GLw2.equ_up01",  90.-8.4/2}, 
    {"GLw2.equ_bt01",  90.+8.4/2}, 
    {"GLw3.bt01_bt02", 90.+8.4+8.4/2}, 
    {"GLw3.bt02_bt03", 90.+8.4+8.4+8.4/2}, 
    {"GLw2.bt03_bt04", 90.+8.4+8.4+8.4+8.4/2}, // FIXME
    {"GLw2.bt04_bt05", 90.+8.4+8.4+8.4+8.4+8.2/2}, 
    {"GLw1.bt05_bt06", 90.+8.4+8.4+8.4+8.4+8.2+8.0/2}, 
    {"GLw1.bt06_bt07", 90.+8.4+8.4+8.4+8.4+8.2+8.0+7.8/2}, 
    {"GLw1.bt07_bt08", 90.+8.4+8.4+8.4+8.4+8.2+8.0+7.8+7.6/2}, 
    {"GLw1.bt08_bt09", 90.+8.4+8.4+8.4+8.4+8.2+8.0+7.8+7.6+7.0/2}, 
    {"GLw1.bt09_bt10", 90.+8.4+8.4+8.4+8.4+8.2+8.0+7.8+7.6+7.0+6.2/2}, 
    {"GLw1.bt10_bt11", 90.+8.4+8.4+8.4+8.4+8.2+8.0+7.8+7.6+7.0+6.2+6.0/2}, 
    {"END", 180}
};

struct PillarData {
    std::string name; // model.label
};

PillarData all_pillar_data[] = {


    {"GZ1.A01_02"},
    {"GZ1.A02_03"},

    {"GZ1.A03_04"},
    {"GZ1.A04_05"},
    {"GZ1.A05_06"},
    {"GZ1.A06_07"},

    {"GZ1.B01_02"},
    {"GZ1.B02_03"},
    {"GZ1.B03_04"},
    {"GZ1.B04_05"},
    {"GZ1.B05_06"},
    {"GZ1.B06_07"},

    {"ZC2.A02_B02"},

    {"ZC2.A03_B03"},
    {"ZC2.A04_B04"},
    {"ZC2.A05_B05"},
    {"ZC2.A06_B06"},

    {"ZC2.A02_B03"},

    {"ZC2.A03_B04"},
    {"ZC2.A04_B05"},
    {"ZC2.A05_B06"},
    {"ZC2.A06_B07"},

    {"ZC2.B01_B01"},
    {"ZC2.B03_B03"},
    {"ZC2.B05_B05"},

    {"ZC2.A03_A03"},
    {"ZC2.A05_A05"},

    {"END"}
};

void
WaterPoolConstruction::makeLatticedShell() {
    SniperPtr<DetSimAlg> detsimalg(*getParent(), "DetSimAlg");
    if (detsimalg.invalid()) {
      std::cout << "Can't Load DetSimAlg" << std::endl;
      assert(0);
    }

    // retrieve detector elements
    //
    for (ShellData* p = all_data_col; p->name!="END"; ++p) {
        std::string name = p->name;

        IDetElement* de = dynamic_cast<IDetElement*>(detsimalg->findTool(name));
        double theta = de->geom_info("theta");

        G4LogicalVolume* daughtervol = de->getLV();
        assert(daughtervol);

        int N = 30;
        double step = 12.*deg;
        double init_phi = 0.0;
        if (name == "GLw1.up10_up11" || name =="GLw1.bt10_bt11") {
            N = 10;
            step = 36.*deg;
            init_phi = 12.*deg;
        }

        double R = 20050.*mm;
        R += de->geom_info("H")/2; // half thickness of H beam

        for (int copyno = 0; copyno < N; ++copyno) { // phi
            double phi = copyno*step+6.0*deg + init_phi;
            double x = R * sin(theta) * cos(phi);
            double y = R * sin(theta) * sin(phi);
            double z = R * cos(theta);

            G4ThreeVector pos(x, y, z);
            G4RotationMatrix rot;
            rot.rotateX(90.*deg);
            rot.rotateZ(90.*deg);
            rot.rotateY(theta+90.*deg);
            rot.rotateZ(phi);
            G4Transform3D trans3d(rot, pos);
            new G4PVPlacement( trans3d,
                               daughtervol,
                               daughtervol->GetName()+"_phys",
                               logicWaterPool,
                               false,
                               copyno
                             );
        }

    }

    for (ShellData* p = all_data_layer; p->name!="END"; ++p) {
        std::string name = p->name;

        IDetElement* de = dynamic_cast<IDetElement*>(detsimalg->findTool(name));
        double theta = de->geom_info("theta");

        G4LogicalVolume* daughtervol = de->getLV();
        assert(daughtervol);

        double R = 20050.*mm;
        R += de->geom_info("H")/2; // half thickness of H beam

        for (int copyno = 0; copyno < 30; ++copyno) { // phi
            double phi = copyno*12.*deg;
            double x = R * sin(theta) * cos(phi);
            double y = R * sin(theta) * sin(phi);
            double z = R * cos(theta);

            G4ThreeVector pos(x, y, z);
            G4RotationMatrix rot;
            rot.rotateY(theta);
            rot.rotateZ(phi);
            G4Transform3D trans3d(rot, pos);
            new G4PVPlacement( trans3d,
                               daughtervol,
                               daughtervol->GetName()+"_phys",
                               logicWaterPool,
                               false,
                               copyno
                             );
        }

    }

    // =======================================================================
    // Pillars
    // =======================================================================
    for (PillarData* p = all_pillar_data; p->name!="END"; ++p) {
        IDetElement* de = dynamic_cast<IDetElement*>(detsimalg->findTool(p->name));
        assert(de);

        G4LogicalVolume* daughtervol = de->getLV();
        assert(daughtervol);

        double R = de->geom_info("R");
        double z = de->geom_info("z");
        double theta = de->geom_info("theta");
        double direct_theta = de->geom_info("direct_theta");
        // if direct_theta is not zero, the beam is connect different column
        double start_phi = 6.*deg;
        if (direct_theta!=0.0) {
            start_phi = 0.;
        }

        for (int copyno = 0; copyno < 30; ++copyno) { // phi
            double phi = start_phi + copyno*12.*deg;

            double x = R * cos(phi);
            double y = R * sin(phi);

            G4ThreeVector pos(x, y, z);
            G4RotationMatrix rot;
            if (direct_theta==0.0) {
                rot.rotateX(90.*deg);
                rot.rotateZ(90.*deg);
                rot.rotateY(theta);
            } else {
                rot.rotateY(direct_theta);
            }
            rot.rotateZ(phi);
            G4Transform3D trans3d(rot, pos);
            new G4PVPlacement( trans3d,
                               daughtervol,
                               daughtervol->GetName()+"_phys",
                               logicWaterPool,
                               false,
                               copyno
                             );

        }

    }



}
