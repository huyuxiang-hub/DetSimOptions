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
// CHANGE HISTORY
// --------------

#include "G4ios.hh"
#include "LSExpEventAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4UImanager.hh"

#include "G4ios.hh"


#include "DetSimAlg/MgrOfAnaElem.h"

using namespace std;

extern G4bool drawEvent;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

LSExpEventAction::LSExpEventAction()
  : drawFlag("all")
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

LSExpEventAction::~LSExpEventAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void LSExpEventAction::BeginOfEventAction(const G4Event* evt)
{
  MgrOfAnaElem::instance().BeginOfEventAction(evt);
  //if(evt->GetEventID()%1000 == 0) G4cout << "Begin of Event --> " << evt->GetEventID() << G4endl;
  G4cout << "Begin of Event --> " << evt->GetEventID() << G4endl;

  if(m_eventID[0]!=-1)
   { 
      std::vector<int>::iterator ret;
      ret = std::find(m_eventID.begin(), m_eventID.end(), evt->GetEventID());
     
      if (ret!=m_eventID.end())
        {
          G4UImanager * UImanager = G4UImanager::GetUIpointer();
          UImanager->ApplyCommand("/tracking/verbose 2");      
        }
      else
        {
          G4UImanager * UImanager = G4UImanager::GetUIpointer();
          UImanager->ApplyCommand("/tracking/verbose 0");
        }
   }

  //G4UImanager* uiMgr = G4UImanager::GetUIpointer();
  //if(evt->GetEventID() == 1537)
  // uiMgr->ApplyCommand("/tracking/verbose 1");
  //else
  //  uiMgr->ApplyCommand("/tracking/verbose 0");


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void LSExpEventAction::EndOfEventAction(const G4Event* evt)
{
  MgrOfAnaElem::instance().EndOfEventAction(evt);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

