//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

// **********************************************************************

#include <boost/python.hpp>
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "LSExpParticleGun.hh"
#include "Randomize.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "LSExpPrimaryGeneratorAction.hh"
#include "G4StateManager.hh"
#include "G4OpticalPhoton.hh"
#include "G4ProcessManager.hh"
#include "G4RadioactiveDecay.hh"
#include "G4DecayTable.hh"
#include "G4GenericIon.hh"
#include "G4NuclideTable.hh"

// FIXME: This is a temporary solution to get data.
#include "SniperKernel/SniperDataPtr.h"
#include "SniperKernel/SniperLog.h"
#include "EvtNavigator/NavBuffer.h"
#include "Event/GenHeader.h"
#include "Event/GenEvent.h"

#include <fstream>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

LSExpPrimaryGeneratorAction::LSExpPrimaryGeneratorAction()
{
  particleGun = new LSExpParticleGun (1);

  particleName = "e+";
  xpos = 0;
  ypos = 0;
  zpos = 0;

  m_scope = 0;


}

LSExpPrimaryGeneratorAction::~LSExpPrimaryGeneratorAction()
{
}

void LSExpPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    // special case: no task in the MT mode
    if (m_isMT && !m_scope) {
        G4ParticleTable* particletbl = G4ParticleTable::GetParticleTable();
        G4ParticleDefinition* particle_def = particletbl->FindParticle(particleName);
        particleGun->SetParticleDefinition(particle_def);
        particleGun->GeneratePrimaryVertex(anEvent);
        return;
    }

    // normal case: load data from event data buffer

    HepMC::GenEvent* gep = 0;
    gep = load_gen_event();
    if (not gep) {
        // TODO raise an Error
        assert(gep);
        return;
    }
    if (SniperLog::logLevel() <= 2) {
        gep->print();
    }

    // set the event id
    anEvent->SetEventID( gep->event_number() );

    // Refer to G4DataHelpers in Dayabay
    // Loop over vertex first
    //     Loop over particles in vertex
      
    // Loop over vertices in the event
    HepMC::GenEvent::vertex_const_iterator
        iVtx = (*gep).vertices_begin(),
        doneVtx = (*gep).vertices_end();
    for (/*nop*/; doneVtx != iVtx; ++iVtx) {
        const HepMC::FourVector& v = (*iVtx)->position();
        G4PrimaryVertex* g4vtx = new G4PrimaryVertex(v.x(), v.y(), v.z(), v.t());

        // Loop over particles in the vertex
        HepMC::GenVertex::particles_out_const_iterator 
            iPart = (*iVtx)->particles_out_const_begin(),
            donePart = (*iVtx)->particles_out_const_end();
        for (/*nop*/; donePart != iPart; ++iPart) {

            // Only keep particles that are important for tracking
            // Use status to pass messages.
            int istatus = (*iPart)->status();
            if (istatus == 0x1000) {
                // NEW: the normal particle, need to use G4 to do radioactivity decay simulation

                // Update the threshold of half life in G4NuclideTable
                G4NuclideTable* nuclide_table = G4NuclideTable::GetInstance();
                nuclide_table->SetThresholdOfHalfLife(0.0);
                

            } else if (istatus != 1) {
                continue;
            }

            G4int pdgcode= (*iPart)-> pdg_id();
            // check the pdgid
            G4ParticleTable* particletbl = G4ParticleTable::GetParticleTable();
            G4ParticleDefinition* particle_def = particletbl->FindParticle(pdgcode);

            // try to look up in ion tables
            //
            // In principle, we need to let nucleus do simulation.
            // However, GenDecay already helps us separate decay chain, 
            // we don't need Geant4 to simulate decay again.
            //
            // Only if the kinetic energy of this nucleus is important,
            // we need to enable following part.
            //
            // -- Tao Lin <lintao@ihep.ac.cn>, 29th Dec 2016
            //
            // if (!particle_def and pdgcode!=20022) {
            //     G4IonTable *theIonTable = particletbl->GetIonTable();
            //     particle_def = theIonTable->GetIon(pdgcode);
            // }

            if (particle_def == 0 and pdgcode != 20022) {

                // if this is a ion, but we don't want it decay,
                // here we just set null ponter for decay table

                // G4cout << "----------- debug ------------" << G4endl;

                G4IonTable *theIonTable = particletbl->GetIonTable();
                particle_def = theIonTable->GetIon(pdgcode);

                // // TODO: force Geant4 simulation stops at Pa234m
                // // The current Geant4 10.4 can't stop simulation from excited to ground state.
                // // So the solution is enabled RadioAnaMgr to kill the secondaries.
                // // -- Tao Lin <lintao@ihep.ac.cn>, 20 May 2020
                //
                // if (pdgcode == 1000902340) { // only if parent is Th234
                //     // Find Pa234m first
                //     G4int Z = 91;
                //     G4int A = 234;
                //     G4double E = 73.92*CLHEP::keV;
                //     G4ParticleDefinition* particle_Pa234m = theIonTable->GetIon(Z, A, E, G4Ions::FloatLevelBase('X'));

                //     G4RadioactiveDecay* g4radiodecay = 0;

                //     // get the radioactivity decay
                //     G4ProcessManager* pmanager = particle_Pa234m->GetProcessManager();
                //     G4int MAXofAtRestLoops  = pmanager->GetAtRestProcessVector()->entries();
                //     G4ProcessVector* fAtRestDoItVector = pmanager->GetAtRestProcessVector();
                //     for (G4int i=0; i<MAXofAtRestLoops; ++i) {
                //         G4VProcess* fCurrentProcess = (*fAtRestDoItVector)[i];
                //         G4RadioactiveDecay* grd = dynamic_cast<G4RadioactiveDecay*>(fCurrentProcess);
                //         if (grd) {
                //             G4cout << "Found G4RadioactiveDecay." << G4endl;
                //             g4radiodecay = grd;
                //             break;
                //         }
                //     }

                //     // // It turns out that grdm will always add a channel from excited state to ground state.
                //     // g4radiodecay->AddUserDecayDataFile(Z, A, "offline/Simulation/DetSimV2/DetSimOptions/data/norad");

                //     // force load
                //     // G4DecayTable* decaytbl = g4radiodecay->LoadDecayTable(*particle_Pa234m);
                //     G4DecayTable* decaytbl = g4radiodecay->GetDecayTable(particle_Pa234m);
                //     if (decaytbl) {
                //         G4cout << "After load decay table, got " << decaytbl << G4endl;
                //         decaytbl->DumpInfo();

                //         // now, create an empty decay table and set it to Pa234m

                //         G4DecayTable* empty_decay_table = new G4DecayTable();
                //         particle_Pa234m->SetDecayTable(empty_decay_table);
                //     }

                //     G4DecayTable* decaytbl_Pa234m = particle_Pa234m->GetDecayTable();
                //     G4cout << "Decay table for Pa234m " << decaytbl_Pa234m << G4endl;
                //     particle_Pa234m->SetPDGStable(true);
                //     particle_Pa234m->SetPDGLifeTime(-1.0);
                // }


                if (!particle_def) {
                    // maybe it is an Isomer > 0
                    theIonTable->CreateAllIon();
                    theIonTable->CreateAllIsomer();

                    // =======================================================
                    // find the radioactive decay process
                    // =======================================================
                    G4RadioactiveDecay* g4radiodecay = 0;

                    // get the radioactivity decay
                    G4ProcessManager* pmanager = G4GenericIon::GenericIon()->GetProcessManager();
                    G4int MAXofAtRestLoops  = pmanager->GetAtRestProcessVector()->entries();
                    G4ProcessVector* fAtRestDoItVector = pmanager->GetAtRestProcessVector();
                    for (G4int i=0; i<MAXofAtRestLoops; ++i) {
                        G4VProcess* fCurrentProcess = (*fAtRestDoItVector)[i];
                        G4RadioactiveDecay* grd = dynamic_cast<G4RadioactiveDecay*>(fCurrentProcess);
                        if (grd) {
                            G4cout << "Found G4RadioactiveDecay." << G4endl;
                            g4radiodecay = grd;
                            break;
                        }
                    }


                    G4int Z=0, A=0, L=0, lvl=0;
                    G4double E=0;
                    theIonTable->GetNucleusByEncoding(pdgcode, Z, A, E, lvl);
                    G4cout << " [x] GetNucleusByEncoding pdg " << pdgcode
                           << " Z: " << Z
                           << " A: " << A
                           << " E: " << E
                           << " lvl: " << lvl
                           << G4endl;

                    // Pa234m
                    if (Z == 91 && A == 234 && lvl == 1) {
                        E = 73.92*CLHEP::keV;
                        particle_def = theIonTable->GetIon(Z, A, E, G4Ions::FloatLevelBase('X'));
                    } else if (Z == 36 && A == 83 && (lvl == 1 || lvl == 2)) { // Kr83m
                        E = 41.5569*CLHEP::keV;
                        particle_def = theIonTable->GetIon(Z, A, E, G4Ions::FloatLevelBase('\0'));

                        if (g4radiodecay) {
                            G4DecayTable* decaytbl = g4radiodecay->GetDecayTable(particle_def);
                            particle_def->SetDecayTable(decaytbl);
                            G4cout << __FILE__ << ":" << __LINE__ << ": "
                                   << "Add decay table " << decaytbl << " to Kr83m " << G4endl;
                            decaytbl->DumpInfo();
                        }

                        // debug
                        G4cout << __FILE__ << ":" << __LINE__ << ": "
                               << " Ion Life Time: " << particle_def->GetIonLifeTime()
                               << G4endl;
                        G4DecayTable* decaytbl = particle_def->GetDecayTable();
                        G4cout << __FILE__ << ":" << __LINE__ << ": "
                               << " Kr83m decay table: " << decaytbl
                               << G4endl;

                        E = 9.4053*CLHEP::keV;
                        G4ParticleDefinition* Kr83_9405eV = theIonTable->GetIon(Z, A, E, G4Ions::FloatLevelBase('\0'));
                        if (g4radiodecay) {
                            G4DecayTable* decaytbl = g4radiodecay->GetDecayTable(Kr83_9405eV);
                            Kr83_9405eV->SetDecayTable(decaytbl);
                            G4cout << __FILE__ << ":" << __LINE__ << ": "
                                   << "Add decay table " << decaytbl << " to Kr83[9.405keV] " << G4endl;
                            decaytbl->DumpInfo();
                        }

                        // In G4, there is threshold_of_half_life(1000.0*ns) in G4NuclideTable. so the 154ns may be ignored?
                        Kr83_9405eV->SetPDGLifeTime(154.4/log(2)*CLHEP::ns);

                        G4cout << __FILE__ << ":" << __LINE__ << ": "
                               << " Life Time of Kr83[9.405keV]: " << Kr83_9405eV->GetPDGLifeTime()
                               << " Ion Life Time of Kr83[9.405keV]: " << Kr83_9405eV->GetIonLifeTime()
                               << G4endl;
                    } else if (Z == 54 && A == 129 && (lvl == 1 || lvl == 2)) { // Xe129m
                        E = 236.140*CLHEP::keV;
                        particle_def = theIonTable->GetIon(Z, A, E, G4Ions::FloatLevelBase('\0'));

                        if (g4radiodecay) {
                            G4DecayTable* decaytbl = g4radiodecay->GetDecayTable(particle_def);
                            particle_def->SetDecayTable(decaytbl);
                            G4cout << __FILE__ << ":" << __LINE__ << ": "
                                   << "Add decay table " << decaytbl << " to Xe129m " << G4endl;
                            decaytbl->DumpInfo();
                        }

                        // debug
                        G4cout << __FILE__ << ":" << __LINE__ << ": "
                               << " Ion Life Time: " << particle_def->GetIonLifeTime()
                               << G4endl;
                        G4DecayTable* decaytbl = particle_def->GetDecayTable();
                        G4cout << __FILE__ << ":" << __LINE__ << ": "
                               << " Xe129m decay table: " << decaytbl
                               << G4endl;
                        
                        E = 39.578*CLHEP::keV;
                        G4ParticleDefinition* Xe129_39578eV = theIonTable->GetIon(Z, A, E, G4Ions::FloatLevelBase('\0'));
                        if (g4radiodecay) {
                            G4DecayTable* decaytbl = g4radiodecay->GetDecayTable(Xe129_39578eV);
                            Xe129_39578eV->SetDecayTable(decaytbl);
                            G4cout << __FILE__ << ":" << __LINE__ << ": "
                                   << "Add decay table " << decaytbl << " to Xe129[39.578keV] " << G4endl;
                            decaytbl->DumpInfo();
                        }

                        // In G4, there is threshold_of_half_life(1000.0*ns) in G4NuclideTable. so the 154ns may be ignored?
                        Xe129_39578eV->SetPDGLifeTime(0.97/log(2)*CLHEP::ns);

                        G4cout << __FILE__ << ":" << __LINE__ << ": "
                               << " Life Time of Xe129[39.578keV]: " << Xe129_39578eV->GetPDGLifeTime()
                               << " Ion Life Time of Xe129[39.578keV]: " << Xe129_39578eV->GetIonLifeTime()
                               << G4endl;
                    } else if (Z == 54 && A == 131 && lvl == 1) { // Xe131m
                        E = 163.930*CLHEP::keV;
                        particle_def = theIonTable->GetIon(Z, A, E, G4Ions::FloatLevelBase('\0'));
                        
                        if (g4radiodecay) {
                            G4DecayTable* decaytbl = g4radiodecay->GetDecayTable(particle_def);
                            particle_def->SetDecayTable(decaytbl);
                            G4cout << __FILE__ << ":" << __LINE__ << ": "
                                   << "Add decay table " << decaytbl << " to Xe131m " << G4endl;
                            decaytbl->DumpInfo();
                        }

                        // debug
                        G4cout << __FILE__ << ":" << __LINE__ << ": "
                               << " Ion Life Time: " << particle_def->GetIonLifeTime()
                               << G4endl;
                        G4DecayTable* decaytbl = particle_def->GetDecayTable();
                        G4cout << __FILE__ << ":" << __LINE__ << ": "
                               << " Xe131m decay table: " << decaytbl
                               << G4endl;
                    } else if (Z == 54 && A == 133 && lvl == 1) { // Xe133m
                        E = 233.221*CLHEP::keV;
                        particle_def = theIonTable->GetIon(Z, A, E, G4Ions::FloatLevelBase('\0'));

                        if (g4radiodecay) {
                            G4DecayTable* decaytbl = g4radiodecay->GetDecayTable(particle_def);
                            particle_def->SetDecayTable(decaytbl);
                            G4cout << __FILE__ << ":" << __LINE__ << ": "
                                   << "Add decay table " << decaytbl << " to Xe133m " << G4endl;
                            decaytbl->DumpInfo();
                        }

                        // debug
                        G4cout << __FILE__ << ":" << __LINE__ << ": "
                               << " Ion Life Time: " << particle_def->GetIonLifeTime()
                               << G4endl;
                        G4DecayTable* decaytbl = particle_def->GetDecayTable();
                        G4cout << __FILE__ << ":" << __LINE__ << ": "
                               << " Xe133m decay table: " << decaytbl
                               << G4endl;
						
                        G4ParticleDefinition* Xe133 = theIonTable->GetIon(Z, A);
                        if (g4radiodecay) {
                            G4DecayTable* decaytbl = g4radiodecay->GetDecayTable(Xe133);
                            Xe133->SetDecayTable(decaytbl);
                            G4cout << __FILE__ << ":" << __LINE__ << ": "
                                   << "Add decay table " << decaytbl << " to Xe133 " << G4endl;
                            decaytbl->DumpInfo();
                        }

                        // In G4, there is threshold_of_half_life(1000.0*ns) in G4NuclideTable. so the 154ns may be ignored?
                        Xe133->SetPDGStable(true);
                        Xe133->SetPDGLifeTime(-1.0);

                        G4cout << __FILE__ << ":" << __LINE__ << ": "
                               << " Life Time of Xe133: " << Xe133->GetPDGLifeTime()
                               << " Ion Life Time of Xe133: " << Xe133->GetIonLifeTime()
                               << G4endl;
                    } else {
                        G4cout << " WARNING: Ion pdgcode " << pdgcode << ", which is isomer or excited," 
                               << " is not supported in the current simulation. " << G4endl;
                        G4cout << " Please contact with Tao Lin <lintao@ihep.ac.cn>, " << G4endl;
                        G4cout << " or you can register the ion in file "
                               << __FILE__ << ":" << __LINE__ << G4endl;
                        particle_def = theIonTable->GetIon(Z, A);
                    }

                    if (particle_def) {
                        G4cout << "Found " << pdgcode
                               << " named " << particle_def->GetParticleName()
                               << G4endl;
                    }

                }

                if (particle_def) {

                    if (istatus == 0x1000) {
                        // Using G4 to do decay.

                        // Force the decay time to zero for the parent nuclide
                        particle_def->SetPDGLifeTime(0.0);
                    } else {

                        // force stable (we don't want to G4 decay the paticle from GenDecay)

                        // G4cout << "theLife: " << particle_def->GetPDGLifeTime() << G4endl;
                        // G4cout << "stable: " << particle_def->GetPDGStable() << G4endl;

                        // check is loaded or not

                        G4DecayTable* decaytbl = particle_def->GetDecayTable();
                        if (decaytbl) {
                            G4cout << "=== Ion pdgcode: [" << pdgcode 
                                   << "] set decay table empty "<< G4endl;

                            particle_def->SetDecayTable(0);
                        }

                        if (!particle_def->GetPDGStable()) {
                            // force stable
                            G4cout << "=== Ion pdgcode: [" << pdgcode 
                                   << "] force stable to avoid decay "<< G4endl;
                            particle_def->SetPDGStable(true);
                            particle_def->SetPDGLifeTime(-1.0);

                        }

                    }
                } else {
                    G4cout << "=== Unknown pdgcode: [" << pdgcode 
                           << "] skip tracking"<< G4endl;
                    // skip this particle
                    continue;

                }



            } else if (pdgcode == 20022) {
                particle_def = G4OpticalPhoton::Definition();
            }
            //
            const HepMC::FourVector& p = (*iPart)->momentum();
            // TODO: What's the unit!
            G4PrimaryParticle* g4prim=new G4PrimaryParticle(particle_def, p.px(), p.py(), p.pz());

            HepMC::ThreeVector pol = (*iPart)->polarization().normal3d();
            g4prim->SetPolarization(pol.x(),pol.y(),pol.z());

            g4vtx->SetPrimary(g4prim);
        }

        if (SniperLog::logLevel() <= 2) {
            g4vtx->Print();
        }

        anEvent->AddPrimaryVertex(g4vtx);

    }
}

HepMC::GenEvent*
LSExpPrimaryGeneratorAction::load_gen_event() {
    // FIXME: Don't know the scope
    SniperDataPtr<JM::NavBuffer>  navBuf(*m_scope, "/Event");
    if (navBuf.invalid()) {
        return 0;
    }
    JM::EvtNavigator* evt_nav = navBuf->curEvt();
    if (not evt_nav) {
        return 0;
    }
    JM::GenHeader* gen_header = dynamic_cast<JM::GenHeader*>(evt_nav->getHeader("/Event/Gen"));
    if (not gen_header) {
        return 0;
    }
    JM::GenEvent* gen_event = dynamic_cast<JM::GenEvent*>(gen_header->event());
    if (not gen_event) {
        return 0;
    }
    return gen_event->getEvent();
}

void LSExpPrimaryGeneratorAction::setMTmode(bool flag) {
    m_isMT = flag;
}
