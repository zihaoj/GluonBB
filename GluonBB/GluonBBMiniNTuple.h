#ifndef GluonBBMiniNTuple_GluonBBMiniNTuple_H
#define GluonBBMiniNTuple_GluonBBMiniNTuple_H

#include <EventLoop/StatusCode.h>
#include <EventLoop/Algorithm.h>

//algorithm wrapper
#include "xAODAnaHelpers/Algorithm.h"

#include "JetSubStructureUtils/BoostedXbbTag.h"
#include "JetSubStructureUtils/BosonTag.h"
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>

// jet trimming
#include <fastjet/tools/Filter.hh>
#include <GluonBB/Helpers.h>

#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"
#include "TTree.h"
#include "TH1D.h"

#include <xAODTruth/TruthParticleContainer.h>
#include "xAODParticleEvent/ParticleContainer.h"
#include "xAODParticleEvent/ParticleAuxContainer.h"
#include <InDetTrackSystematicsTools/InDetTrackSmearingTool.h>

class HelpTreeBase;

class GluonBBMiniNTuple : public xAH::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:

  std::string m_name;                  
  TH1D* m_cutflowHist;    //!
  TH1D* m_cutflowHistW;   //!
  int m_cutflowFirst;     //!
  int m_iCutflow;         //!
  bool m_isMC;            //!

  xAOD::TEvent *m_event;               //!
  xAOD::TStore *m_store;               //!

  std::string m_resolvedJetsName;
  std::string m_resolvedHcandName;
  std::string m_lepTopCandName;
  std::string m_boostedGluoncandName;

  std::string m_evtDetailStr;
  std::string m_trigDetailStr;
  std::string m_resolvedJetDetailStr;
  std::string m_boostedJetDetailStr;
  std::string m_boostedFatJetDetailStr;
  std::string m_lepTopJetDetailStr;
  std::string m_truthDetailStr;
  std::string m_inCaloJetName;
  std::string m_inTrackJetName;
  std::string m_inTruthParticleName;
  std::string m_inTruthJetName;
  //std::string m_inOtherResolvedJetsName;
  std::string m_muonContainerName;
  std::string m_muonDetailStr;
  bool m_doResolutionStudy;
  bool m_debug;
  bool m_doGluonBBTagging;
  bool m_doFatJetMassCut;
  double m_FatJetPtSkimCut;
  double m_FatJetPtTruthCut;
  int m_FatJetTruthNTrkJetMatched;
  bool m_storeLeptonVeto;
  bool m_storeMETVeto;
  bool m_doResolved;
  bool m_doLeptop;
  bool m_doBoosted;
  bool m_ignoreRecoCuts;
  std::string m_eventCuts;
  std::string m_resolvedSysName;      
  std::string m_boostedSysName;     

  InDet::InDetTrackSmearingTool *m_InDetTrackSmearingTool;

private:

  std::map< std::string, HelpTreeBase* > m_helpTree; //!
  JetSubStructureUtils::BoostedXbbTag* m_higgsTaggerLoose = nullptr;  //! 
  JetSubStructureUtils::BoostedXbbTag* m_higgsTaggerMedium = nullptr; //!
  JetSubStructureUtils::BoostedXbbTag* m_higgsTaggerTight = nullptr;  //!
  JetSubStructureUtils::BosonTag*      m_WbosonTaggerMedium = nullptr; //!
  JetSubStructureUtils::BosonTag*      m_ZbosonTaggerMedium = nullptr; //!
  JetSubStructureUtils::BosonTag*      m_WbosonTaggerTight = nullptr; //!
  JetSubStructureUtils::BosonTag*      m_ZbosonTaggerTight = nullptr; //!

private: 

  // Truth variables
  double b_truth_mtt;


  // Event selection flags
  bool b_EventPass_GluonBBL1;
  bool b_EventPass_GluonBBHLT;
  bool b_EventPass_GluonBBTrig;
  float Total_Number_of_Events;


  std::vector< std::vector<int> >   b_jet_ak4emtopo_asso_idx_in_resolvedJets;



  // Leptonic Tops
  std::vector<float> b_leptop_pt;
  std::vector<float> b_leptop_eta;
  std::vector<float> b_leptop_phi;
  std::vector<float> b_leptop_m;
  std::vector<float> b_leptop_dRmj;

  std::vector<int>   b_leptop_jet_idx_in_lepTopJets;
  
  std::vector<float> b_leptop_muon_pt;
  std::vector<float> b_leptop_muon_eta;
  std::vector<float> b_leptop_muon_phi;
  std::vector<float> b_leptop_muon_m;
  std::vector<float> b_leptop_muon_ptcone20;

  std::vector<float> b_leptop_met      ;
  std::vector<float> b_leptop_metphi   ;
  std::vector<float> b_leptop_dPhimmet ;
  std::vector<float> b_leptop_Mt       ;

  // Higgs candidates from boosted selection
  int   b_gluon_boosted_n;
  float b_gluon_boosted_pt;
  float b_gluon_boosted_eta;
  float b_gluon_boosted_phi;
  float b_gluon_boosted_m;
  float b_gluon_boosted_dRjj;
  float b_gluon_boosted_D2;
  float b_gluon_boosted_Tau21;
  float b_gluon_boosted_Tau21WTA;
  int   b_gluon_boosted_nTrack;
  int   b_gluon_boosted_nTruthB;
  int   b_gluon_boosted_nTruthC;

  int b_gluon_boosted_truth_n;
  float b_gluon_boosted_truth_pt;
  float b_gluon_boosted_truth_eta;
  float b_gluon_boosted_truth_phi;
  float b_gluon_boosted_truth_m;
  float b_gluon_boosted_truth_dRjj;

  std::vector<int>                  b_jet_ak2track_asso_n;
  std::vector<int>                  b_jet_ak2track_asso_n_addl;
  std::vector< std::vector<float> > b_jet_ak2track_asso_pt;
  std::vector< std::vector<float> > b_jet_ak2track_asso_eta;
  std::vector< std::vector<float> > b_jet_ak2track_asso_phi;
  std::vector< std::vector<float> > b_jet_ak2track_asso_m;
  std::vector< std::vector<float> > b_jet_ak2track_asso_MV2c00;
  std::vector< std::vector<float> > b_jet_ak2track_asso_MV2c10;
  std::vector< std::vector<float> > b_jet_ak2track_asso_MV2c20;
  std::vector< std::vector<float> > b_jet_ak2track_asso_MV2c100;
  std::vector< std::vector<int> > b_jet_ak2track_asso_TruthB;
  std::vector< std::vector<int> > b_jet_ak2track_asso_TruthC;
  std::vector< std::vector<int> > b_jet_ak2track_asso_TruthL;
  std::vector< std::vector<int> > b_jet_ak2track_asso_TruthFlav;
  std::vector< std::vector< std::vector<float> > >  b_jet_ak2track_asso_sys;
  std::vector<std::string> b_jet_ak2track_asso_sysname;
  std::vector<float>  b_boosted_bevent_sys;

  float b_jet_leading_ak2track_asso_pt;
  float b_jet_leading_ak2track_asso_eta;
  float b_jet_leading_ak2track_asso_phi;
  float b_jet_leading_ak2track_asso_m;
  float b_jet_leading_ak2track_asso_Sd0;
  float b_jet_leading_ak2track_asso_Sd0_Sub;
  float b_jet_leading_ak2track_asso_IP3D;
  float b_jet_leading_ak2track_asso_IP2D;
  int   b_jet_leading_ak2track_asso_nMuon;
  float b_jet_leading_ak2track_asso_MV2c10;
  float b_jet_leading_ak2track_asso_MV2c20;
  int b_jet_leading_ak2track_asso_TruthFlav;
  int b_jet_leading_ak2track_asso_TruthB;
  int b_jet_leading_ak2track_asso_TruthC;
  int b_jet_leading_ak2track_asso_TruthL;
  std::vector< std::vector<float> > b_jet_leading_ak2track_asso_sys;
  std::vector<std::string> b_jet_leading_ak2track_asso_sysname;

  std::vector<int> b_jet_leading_ak2track_trk_IBLHits;
  std::vector<int> b_jet_leading_ak2track_trk_nBLHits;
  std::vector<int> b_jet_leading_ak2track_trk_expectBLHits;
  std::vector<int> b_jet_leading_ak2track_trk_nPixelHits;
  std::vector<int> b_jet_leading_ak2track_trk_nPixelSharedHits;
  std::vector<int> b_jet_leading_ak2track_trk_nPixelSplitHits;
  std::vector<int> b_jet_leading_ak2track_trk_nSCTHits;
  std::vector<float> b_jet_leading_ak2track_trk_pT;
  std::vector<float> b_jet_leading_ak2track_trk_eta;
  int b_jet_leading_ak2track_ntrk;

  std::vector<int> b_jet_subleading_ak2track_trk_IBLHits;
  std::vector<int> b_jet_subleading_ak2track_trk_nBLHits;
  std::vector<int> b_jet_subleading_ak2track_trk_expectBLHits;
  std::vector<int> b_jet_subleading_ak2track_trk_nPixelHits;
  std::vector<int> b_jet_subleading_ak2track_trk_nPixelSharedHits;
  std::vector<int> b_jet_subleading_ak2track_trk_nPixelSplitHits;
  std::vector<int> b_jet_subleading_ak2track_trk_nSCTHits;


  float b_jet_subleading_ak2track_asso_pt;
  float b_jet_subleading_ak2track_asso_eta;
  float b_jet_subleading_ak2track_asso_phi;
  float b_jet_subleading_ak2track_asso_m;
  float b_jet_subleading_ak2track_asso_Sd0;
  float b_jet_subleading_ak2track_asso_Sd0_Sub;
  float b_jet_subleading_ak2track_asso_IP3D;
  float b_jet_subleading_ak2track_asso_IP2D;
  int   b_jet_subleading_ak2track_asso_nMuon;
  float b_jet_subleading_ak2track_asso_MV2c10;
  float b_jet_subleading_ak2track_asso_MV2c20;
  int b_jet_subleading_ak2track_asso_TruthFlav;
  int b_jet_subleading_ak2track_asso_TruthB;
  int b_jet_subleading_ak2track_asso_TruthC;
  int b_jet_subleading_ak2track_asso_TruthL;
  std::vector< std::vector<float> > b_jet_subleading_ak2track_asso_sys;
  std::vector<std::string> b_jet_subleading_ak2track_asso_sysname;

  // truth B-hadron tagged track jets 
  float b_jet_leading_ak2track_asso_truth_pt;
  float b_jet_leading_ak2track_asso_truth_eta;
  float b_jet_leading_ak2track_asso_truth_phi;
  float b_jet_leading_ak2track_asso_truth_m;
  int b_jet_leading_ak2track_asso_truth_TruthFlav;
  int b_jet_leading_ak2track_asso_truth_TruthB;
  int b_jet_leading_ak2track_asso_truth_TruthC;
  int b_jet_leading_ak2track_asso_truth_TruthL;

  float b_jet_subleading_ak2track_asso_truth_pt;
  float b_jet_subleading_ak2track_asso_truth_eta;
  float b_jet_subleading_ak2track_asso_truth_phi;
  float b_jet_subleading_ak2track_asso_truth_m;
  int b_jet_subleading_ak2track_asso_truth_TruthFlav;
  int b_jet_subleading_ak2track_asso_truth_TruthB;
  int b_jet_subleading_ak2track_asso_truth_TruthC;
  int b_jet_subleading_ak2track_asso_truth_TruthL;



  float  b_Leading_akt4_pt;
  float  b_Leading_akt4_eta;
  float  b_Leading_akt4_phi;
  float  b_Leading_akt4_m;


  // Loose lepton for veto purpose
  int                 b_n_muons_veto;
  int                 b_n_electrons_veto;

  // MET for veto purpose
  double              b_METsum;
  double              b_METphi;

  float m_mcEventWeight; 
  float m_weight;
  float m_weight_xs;



  /* //Boosted Branches */
  /* std::vector<int> b_n_asso_with_sv; */
  /* int b_n_pass_htagger_loose; */
  /* int b_n_pass_htagger_medium; */
  /* int b_n_pass_htagger_tight; */

  //
  // For cutflow counting
  //
  void passCut();

  //
  // For mtt slices
  //
  bool getEventMtt(const xAOD::TruthParticleContainer* Truths, double& mtt);

  //
  // KillPileupFluctuation (remove high weight event in JZXW). Inherited from ProofAna
  //
  int killPileupFluctuation();    // return 1 means event is good; 0 means event is not good

  //
  // Substructure Computation 
  //

  double getD2(const xAOD::Jet* inputJet);
  double getTau21(const xAOD::Jet* inputJet, bool doWTA);

  std::vector<int> JetMuondRMatch(std::vector<const xAOD::Jet*> inputjets, const xAOD::MuonContainer* muons, float radius);
  //float BTAGIP(const xAOD::Jet* inputjet, std::string flavor);

public:

  // this is a standard constructor
  GluonBBMiniNTuple ();      

  bool IsWeaklyDecayingBHadron(int mc_pdgId);              //!
  bool GhostMatchTruthTrkJetToLargeRJet(  std::vector<TLorentzVector> TruthTrackJets4Vector_Preselected,
					  std::vector<int> TruthTrackJets4Vector_BMatched);
					  

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);           //!
  virtual EL::StatusCode fileExecute ();                    //!
  virtual EL::StatusCode histInitialize ();                 //!
  virtual EL::StatusCode changeInput (bool firstFile);      //!
  virtual EL::StatusCode initialize ();                     //!
  virtual EL::StatusCode execute ();                        //!
  virtual EL::StatusCode postExecute ();                    //!
  virtual EL::StatusCode finalize ();                       //!
  virtual EL::StatusCode histFinalize ();                   //!

  virtual EL::StatusCode AddTree (std::string);      //!


  // these are the functions not inherited from Algorithm
  EL::StatusCode executeSingle(std::string resolvedSys="", 
                               std::string boostedSys ="", 
			       bool countEvents = false);


  // this is needed to distribute the algorithm to the workers
  ClassDef(GluonBBMiniNTuple, 1);                                 //!
};

#endif
