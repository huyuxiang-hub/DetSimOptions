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


//
// **********************************************************************

#ifndef LSExpPrimaryGeneratorAction_h
#define LSExpPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

#include <vector>
#include <map>
#include <string>

namespace HepMC {
    class GenEvent;
}

class G4ParticleGun;
class LSExpParticleGun;
class G4Event;

class Task;

class LSExpPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  LSExpPrimaryGeneratorAction();
  ~LSExpPrimaryGeneratorAction();

public:
  void GeneratePrimaries(G4Event* anEvent);

public:
  void setScope(Task* scope) {m_scope = scope;}
  Task* getScope() {return m_scope;}

public:
  HepMC::GenEvent* load_gen_event();
public:
    // MT related
    void setMTmode(bool flag);

private:
  LSExpParticleGun* particleGun;

  G4String particleName;
  G4double xpos;
  G4double ypos;
  G4double zpos;

  // Verbosity
  G4int m_verbosity;
private:
  Task* m_scope;

private:
    bool m_isMT; // disabled by default

};

#endif



