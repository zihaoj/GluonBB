#ifndef SetGluonBBEventCuts_H
#define SetGluonBBEventCuts_H

#include <EventLoop/StatusCode.h>
#include <EventLoop/Algorithm.h>
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"
#include <xAODTruth/TruthParticleContainer.h>

// Particle 
#include "xAODParticleEvent/Particle.h"
#include "xAODParticleEvent/ParticleContainer.h"
#include "xAODParticleEvent/ParticleAuxContainer.h"

#include <TH1.h>

#include <AsgTools/AnaToolHandle.h>
#include <JetInterface/IJetExecuteTool.h>


// EDM include(s):
#ifndef __CINT__
#include "xAODJet/JetContainer.h"
#endif

class SetGluonBBEventCuts : public EL::Algorithm
{

 public:

  std::string m_name;

  xAOD::TEvent *m_event; //!
  xAOD::TStore *m_store; //!
  bool        m_debug;
  bool        m_dothirdJet;
  bool        m_dothirdTrkJet;

  //std::string m_outDiJetName; 
  std::string m_inJetName; 
  std::string m_inAlgoName; 
  std::string m_trackJet;
  std::string m_inputAlgo;
  std::string m_outJetName;
  std::string m_outputAlgo;
  std::string m_outLeadingAkt4Jet;
  std::string m_inTruthParticleName;
  std::string m_inTruthJetName;

  asg::AnaToolHandle<IJetExecuteTool> m_jetReclusteringTool;

  float m_leadingJetPtCut;
  float m_subleadingJetPtCut;
  float m_dEtaCut;
  float m_trackJetPtCut;
  float m_trackJetEtaCut;
  float m_trackJetBtagCut;
  bool  m_isMC;
  private:

  //
  //  No //! for these guys as they are configuration
  //

  EL::StatusCode GluonBBSelection (std::string systName);
  EL::StatusCode GhostMatchTruthTrkJetToLargeRJet( xAOD::ParticleContainer* TruthLargeRJets,  xAOD::ParticleContainer* TruthTrackJets);
  bool IsWeaklyDecayingBHadron(int mc_pdgId);

  //
  // sorting function
  //

  static bool SortPt(const xAOD::IParticle* j1, const xAOD::IParticle* j2) { return j1->pt() > j2->pt(); }

  //
  // track-jet selection
  //

  bool SelectTrackJet(const xAOD::Jet* TrackJet);

  public:

  // this is a standard constructor
  SetGluonBBEventCuts ();

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();

  // these are the functions not inherited from Algorithm
  virtual EL::StatusCode configure ();

  // this is needed to distribute the algorithm to the workers
  ClassDef(SetGluonBBEventCuts, 1);
};

#endif
