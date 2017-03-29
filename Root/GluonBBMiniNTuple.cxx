#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <EventLoop/OutputStream.h>

#include <xAODJet/JetContainer.h>
#include <xAODTracking/VertexContainer.h>
#include <xAODTracking/TrackParticlexAODHelpers.h>
#include <xAODTracking/TrackParticle.h>

#include <xAODEventInfo/EventInfo.h>
#include <AthContainers/ConstDataVector.h>
#include <SampleHandler/MetaFields.h>

#include "xAODParticleEvent/ParticleContainer.h"

#include <xAODAnaHelpers/HelpTreeBase.h>
#include <GluonBB/GluonBBMiniNTuple.h>

#include <xAODAnaHelpers/HelperFunctions.h>
#include <xAODAnaHelpers/HelperClasses.h>
#include <xAODAnaHelpers/tools/ReturnCheck.h>
#include <xAODAnaHelpers/tools/ReturnCheckConfig.h>

#include "xAODParticleEvent/ParticleContainer.h"
#include "xAODParticleEvent/ParticleAuxContainer.h"


#include "TEnv.h"
#include "TSystem.h"

//#define DEBUG std::cerr << __FILE__ << "::" << __LINE__ << std::endl
 
// this is needed to distribute the algorithm to the workers
ClassImp(GluonBBMiniNTuple)

using std::cout;  using std::endl;
using std::string; using std::vector;

//class InDetTrackSmearingTool;

GluonBBMiniNTuple :: GluonBBMiniNTuple () :
  m_name(""),
  m_cutflowHist(0),
  m_cutflowHistW(0),
  m_cutflowFirst(0),
  m_iCutflow(0),
  m_isMC(false),
  m_resolvedJetsName(""),
  m_resolvedHcandName(""),
  m_lepTopCandName(""),
  m_boostedGluoncandName(""),
  m_evtDetailStr(""),
  m_trigDetailStr(""),
  m_resolvedJetDetailStr(""),
  m_boostedJetDetailStr(""),
  m_lepTopJetDetailStr(""),
  m_truthDetailStr(""),
  m_inCaloJetName(""),
  m_inTrackJetName(""),
  m_inTruthParticleName(""),
  m_muonContainerName(""),
  m_muonDetailStr(""),
  m_inTruthJetName(""),
  m_doResolutionStudy(false),
  m_debug(false), 
  m_doGluonBBTagging(true), //defualt turned on!
  m_doFatJetMassCut(true),
  m_FatJetPtSkimCut(450.),  // in GeV
  m_FatJetPtTruthCut(450.),  // in GeV
  m_FatJetTruthNTrkJetMatched(2),  
  m_storeLeptonVeto(false),
  m_storeMETVeto(false),
  m_doResolved(true),
  m_doLeptop(true),
  m_doBoosted(true),
  m_eventCuts(""),
  m_resolvedSysName(""),
  m_boostedSysName(""),
  m_weight(1.0),
  m_weight_xs(1.0),
  m_InDetTrackSmearingTool(nullptr),
  m_ignoreRecoCuts(false)
{
  this->SetName("GluonBBMiniNTuple"); // needed if you want to retrieve this algo with wk()->getAlg(ALG_NAME) downstream
}

EL::StatusCode GluonBBMiniNTuple :: setupJob (EL::Job& job)
{
  Info("setupJob()", "Calling setupJob \n");
  job.useXAOD();
  xAOD::Init("GluonBBMiniNTuple").ignore();

  EL::OutputStream outForTree(m_name.c_str());
  job.outputAdd (outForTree);
  Info("setupJob()", "Ready\n");

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode GluonBBMiniNTuple :: initialize ()
{
  Info("initialize()", m_name.c_str());
  m_event = wk()->xaodEvent();
  m_store = wk()->xaodStore();

  if(m_debug) Info("initialize()", "after add store");
  //if (m_resolvedSysName.empty() && m_boostedSysName.empty()) this->AddTree("");

  //string recs_file = "$ROOTCOREBIN/data/JetSubStructureUtils/config_13TeV_20150710_Htagging.dat"; //this one is old
  string recs_file = "$ROOTCOREBIN/data/JetSubStructureUtils/config_13TeV_Htagging_MC15_Prerecommendations_20150812.dat";

  if(m_debug) Info("initialize()", "after rec file");

  //
  // Set isMC flag
  //
  const xAOD::EventInfo* eventInfo(nullptr);
  RETURN_CHECK("ResonanceAlgorithm::initialize()", HelperFunctions::retrieve(eventInfo, "EventInfo", m_event, m_store, m_debug), "");
  m_isMC = ( eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ) ? true : false;

  //
  // Interlock on m_doResolutionStudy
  //
  if(!m_isMC) m_doResolutionStudy = false;
 
  //
  // Cut Flow 
  //
  TFile *file = wk()->getOutputFile("cutflow");
  m_cutflowHist  = (TH1D*)file->Get("cutflow");
  m_cutflowHistW = (TH1D*)file->Get("cutflow_weighted");

  m_cutflowFirst = m_cutflowHist->GetXaxis()->FindBin("GluonBBMiniTreeAll");
  m_cutflowHistW->GetXaxis()->FindBin("GluonBBMiniTreeAll");

  m_cutflowHist ->GetXaxis()->FindBin("GluonBBMiniTreePreSel");
  m_cutflowHistW->GetXaxis()->FindBin("GluonBBMiniTreePreSel");


  // set up indet_smearing_tool
  std::string indet_smearing_tool_name = std::string("indet_smearing_tool");
  m_InDetTrackSmearingTool = new InDet::InDetTrackSmearingTool( indet_smearing_tool_name);
  RETURN_CHECK( "InDetTrackSmearingTool::initialize()", m_InDetTrackSmearingTool->initialize(), "Failed to properly initialize the InDetTrackSmearingTool");
  CP::SystematicSet systSet = {InDet::TrackSystematicMap[InDet::TRK_RES_D0_MEAS]};
			       //InDet::TrackSystematicMap[InDet::TRK_RES_Z0_MEAS]};
  auto systCode = m_InDetTrackSmearingTool->applySystematicVariation(systSet);


  if(m_debug) Info("initialize()", "left");
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode GluonBBMiniNTuple::AddTree(string syst = "")
{
  Info("AddTree()", "%s", m_name.c_str() );
  // needed here and not in initalize since this is called first

  string treeName = "GluonBBMiniNTuple";
  if (!syst.empty()) treeName += syst;
  TTree * outTree = new TTree(treeName.c_str(), treeName.c_str());
  if( !outTree ) {
    Error("AddTree()","Failed to instantiate output tree!");
    return EL::StatusCode::FAILURE;
  }

  // get the file we created already
  TFile* treeFile = wk()->getOutputFile( m_name.c_str() );
  outTree->SetDirectory( treeFile );


  m_helpTree[syst] = new HelpTreeBase( m_event, outTree, treeFile );

  // if want to add to same file as ouput histograms
  // wk()->addOutput( outTree );

  m_helpTree[syst]->AddEvent(m_evtDetailStr);
  m_helpTree[syst]->AddTrigger(m_trigDetailStr);
  m_helpTree[syst]->AddJets(m_boostedJetDetailStr,     "boostedJets");

  //if(!m_muonContainerName.empty())
  //m_helpTree[syst]->AddMuons(m_muonDetailStr);

  // Event selection flags
  outTree->Branch("EventPass_GluonBBL1", &b_EventPass_GluonBBL1);
  outTree->Branch("EventPass_GluonBBHLT", &b_EventPass_GluonBBHLT);
  outTree->Branch("EventPass_GluonBBTrig", &b_EventPass_GluonBBTrig);
  outTree->Branch("Total_Number_of_Events", &Total_Number_of_Events);

  // Higgs candidates from boosted selection
  outTree->Branch("gluon_boosted_n", &b_gluon_boosted_n, "gluon_boosted_n/I");
  outTree->Branch("gluon_boosted_pt", &b_gluon_boosted_pt, "gluon_boosted_pt/F");
  outTree->Branch("gluon_boosted_eta", &b_gluon_boosted_eta, "gluon_boosted_eta/F");
  outTree->Branch("gluon_boosted_phi", &b_gluon_boosted_phi, "gluon_boosted_phi/F");
  outTree->Branch("gluon_boosted_m", &b_gluon_boosted_m, "gluon_boosted_m/F");
  outTree->Branch("gluon_boosted_D2", &b_gluon_boosted_D2, "gluon_boosted_D2/F");
  outTree->Branch("gluon_boosted_Tau21", &b_gluon_boosted_Tau21, "gluon_boosted_Tau21/F");
  outTree->Branch("gluon_boosted_Tau21WTA", &b_gluon_boosted_Tau21WTA, "gluon_boosted_Tau21WTA/F");
  outTree->Branch("gluon_boosted_nTrack", &b_gluon_boosted_nTrack, "gluon_boosted_nTrack/I");
  outTree->Branch("gluon_boosted_nTruthB", &b_gluon_boosted_nTruthB, "gluon_boosted_nTruthB/I");
  outTree->Branch("gluon_boosted_nTruthC", &b_gluon_boosted_nTruthC, "gluon_boosted_nTruthC/I");

  outTree->Branch("gluon_boosted_truth_n", &b_gluon_boosted_truth_n, "gluon_boosted_truth_n/I");
  outTree->Branch("gluon_boosted_truth_pt", &b_gluon_boosted_truth_pt, "gluon_boosted_truth_pt/F");
  outTree->Branch("gluon_boosted_truth_eta", &b_gluon_boosted_truth_eta, "gluon_boosted_truth_eta/F");
  outTree->Branch("gluon_boosted_truth_phi", &b_gluon_boosted_truth_phi, "gluon_boosted_truth_phi/F");
  outTree->Branch("gluon_boosted_truth_m", &b_gluon_boosted_truth_m, "gluon_boosted_truth_m/F");
  outTree->Branch("gluon_boosted_truth_dRjj", &b_gluon_boosted_truth_dRjj, "gluon_boosted_truth_dRjj/F");

  outTree->Branch("jet_ak2track_asso_n", &b_jet_ak2track_asso_n);
  outTree->Branch("jet_ak2track_asso_n_addl", &b_jet_ak2track_asso_n_addl);
  outTree->Branch("jet_ak2track_asso_pt", &b_jet_ak2track_asso_pt);
  outTree->Branch("jet_ak2track_asso_eta", &b_jet_ak2track_asso_eta);
  outTree->Branch("jet_ak2track_asso_phi", &b_jet_ak2track_asso_phi);
  outTree->Branch("jet_ak2track_asso_m", &b_jet_ak2track_asso_m);
  outTree->Branch("jet_ak2track_asso_MV2c00", &b_jet_ak2track_asso_MV2c00);
  outTree->Branch("jet_ak2track_asso_MV2c10", &b_jet_ak2track_asso_MV2c10);
  outTree->Branch("jet_ak2track_asso_MV2c20", &b_jet_ak2track_asso_MV2c20);
  outTree->Branch("jet_ak2track_asso_MV2c100", &b_jet_ak2track_asso_MV2c100);
  outTree->Branch("jet_ak2track_asso_sys", &b_jet_ak2track_asso_sys);
  outTree->Branch("jet_ak2track_asso_sysname", &b_jet_ak2track_asso_sysname);
  outTree->Branch("jet_ak2track_asso_TruthFlav", &b_jet_ak2track_asso_TruthFlav);
  outTree->Branch("jet_ak2track_asso_TruthB", &b_jet_ak2track_asso_TruthB);
  outTree->Branch("jet_ak2track_asso_TruthC", &b_jet_ak2track_asso_TruthC);
  outTree->Branch("jet_ak2track_asso_TruthL", &b_jet_ak2track_asso_TruthL);

  outTree->Branch("jet_leading_ak2track_asso_pt", &b_jet_leading_ak2track_asso_pt, "jet_leading_ak2track_asso_pt/F");
  outTree->Branch("jet_leading_ak2track_asso_eta", &b_jet_leading_ak2track_asso_eta, "jet_leading_ak2track_asso_eta/F");
  outTree->Branch("jet_leading_ak2track_asso_phi", &b_jet_leading_ak2track_asso_phi, "jet_leading_ak2track_asso_phi/F");
  outTree->Branch("jet_leading_ak2track_asso_m", &b_jet_leading_ak2track_asso_m, "jet_leading_ak2track_asso_m/F");
  outTree->Branch("jet_leading_ak2track_asso_Sd0", &b_jet_leading_ak2track_asso_Sd0, "jet_leading_ak2track_asso_Sd0/F");
  outTree->Branch("jet_leading_ak2track_asso_Sd0_Sub", &b_jet_leading_ak2track_asso_Sd0_Sub, "jet_leading_ak2track_asso_Sd0_Sub/F");
  outTree->Branch("jet_leading_ak2track_asso_MV2c10", &b_jet_leading_ak2track_asso_MV2c10, "jet_leading_ak2track_asso_MV2c10/F");
  outTree->Branch("jet_leading_ak2track_asso_MV2c20", &b_jet_leading_ak2track_asso_MV2c20, "jet_leading_ak2track_asso_MV2c20/F");
  outTree->Branch("jet_leading_ak2track_asso_IP3D", &b_jet_leading_ak2track_asso_IP3D, "jet_leading_ak2track_asso_IP3D/F");
  outTree->Branch("jet_leading_ak2track_asso_IP2D", &b_jet_leading_ak2track_asso_IP2D, "jet_leading_ak2track_asso_IP2D/F");
  outTree->Branch("jet_leading_ak2track_asso_nMuon", &b_jet_leading_ak2track_asso_nMuon, "jet_leading_ak2track_asso_nMuon/I");
  outTree->Branch("jet_leading_ak2track_asso_sys", &b_jet_leading_ak2track_asso_sys);
  outTree->Branch("jet_leading_ak2track_asso_sysname", &b_jet_leading_ak2track_asso_sysname);
  outTree->Branch("jet_leading_ak2track_asso_TruthFlav", &b_jet_leading_ak2track_asso_TruthFlav, "jet_leading_ak2track_asso_TruthFlav/I");
  outTree->Branch("jet_leading_ak2track_asso_TruthB", &b_jet_leading_ak2track_asso_TruthB, "jet_leading_ak2track_asso_TruthB/I");
  outTree->Branch("jet_leading_ak2track_asso_TruthC", &b_jet_leading_ak2track_asso_TruthC, "jet_leading_ak2track_asso_TruthC/I");

  outTree->Branch("jet_leading_ak2track_trk_IBLHits", &b_jet_leading_ak2track_trk_IBLHits);
  outTree->Branch("jet_leading_ak2track_trk_nBLHits", &b_jet_leading_ak2track_trk_nBLHits);
  outTree->Branch("jet_leading_ak2track_trk_expectBLHits", &b_jet_leading_ak2track_trk_expectBLHits);
  outTree->Branch("jet_leading_ak2track_trk_nPixelHits", &b_jet_leading_ak2track_trk_nPixelHits);
  outTree->Branch("jet_leading_ak2track_trk_nPixelSharedHits", &b_jet_leading_ak2track_trk_nPixelSharedHits);
  outTree->Branch("jet_leading_ak2track_trk_nPixelSplitHits", &b_jet_leading_ak2track_trk_nPixelSplitHits);
  outTree->Branch("jet_leading_ak2track_trk_nSCTHits", &b_jet_leading_ak2track_trk_nSCTHits);
  outTree->Branch("jet_leading_ak2track_trk_pT", &b_jet_leading_ak2track_trk_pT);
  outTree->Branch("jet_leading_ak2track_trk_eta", &b_jet_leading_ak2track_trk_eta);
  outTree->Branch("jet_leading_ak2track_ntrk",  &b_jet_leading_ak2track_ntrk);

  outTree->Branch("jet_subleading_ak2track_asso_pt", &b_jet_subleading_ak2track_asso_pt, "jet_subleading_ak2track_asso_pt/F");
  outTree->Branch("jet_subleading_ak2track_asso_eta", &b_jet_subleading_ak2track_asso_eta, "jet_subleading_ak2track_asso_eta/F");
  outTree->Branch("jet_subleading_ak2track_asso_phi", &b_jet_subleading_ak2track_asso_phi, "jet_subleading_ak2track_asso_phi/F");
  outTree->Branch("jet_subleading_ak2track_asso_m", &b_jet_subleading_ak2track_asso_m, "jet_subleading_ak2track_asso_m/F");
  outTree->Branch("jet_subleading_ak2track_asso_Sd0", &b_jet_subleading_ak2track_asso_Sd0, "jet_subleading_ak2track_asso_Sd0/F");
  outTree->Branch("jet_subleading_ak2track_asso_Sd0_Sub", &b_jet_subleading_ak2track_asso_Sd0_Sub, "jet_subleading_ak2track_asso_Sd0_Sub/F");
  outTree->Branch("jet_subleading_ak2track_asso_MV2c10", &b_jet_subleading_ak2track_asso_MV2c10, "jet_subleading_ak2track_asso_MV2c10/F");
  outTree->Branch("jet_subleading_ak2track_asso_MV2c20", &b_jet_subleading_ak2track_asso_MV2c20, "jet_subleading_ak2track_asso_MV2c20/F");
  outTree->Branch("jet_subleading_ak2track_asso_IP3D", &b_jet_subleading_ak2track_asso_IP3D, "jet_subleading_ak2track_asso_IP3D/F");
  outTree->Branch("jet_subleading_ak2track_asso_IP2D", &b_jet_subleading_ak2track_asso_IP2D, "jet_subleading_ak2track_asso_IP2D/F");
  outTree->Branch("jet_subleading_ak2track_asso_nMuon", &b_jet_subleading_ak2track_asso_nMuon, "jet_subleading_ak2track_asso_nMuon/I");
  outTree->Branch("jet_subleading_ak2track_asso_sys", &b_jet_subleading_ak2track_asso_sys);
  outTree->Branch("jet_subleading_ak2track_asso_sysname", &b_jet_subleading_ak2track_asso_sysname);
  outTree->Branch("jet_subleading_ak2track_asso_TruthFlav", &b_jet_subleading_ak2track_asso_TruthFlav, "jet_subleading_ak2track_asso_TruthFlav/I");
  outTree->Branch("jet_subleading_ak2track_asso_TruthB", &b_jet_subleading_ak2track_asso_TruthB, "jet_subleading_ak2track_asso_TruthB/I");
  outTree->Branch("jet_subleading_ak2track_asso_TruthC", &b_jet_subleading_ak2track_asso_TruthC, "jet_subleading_ak2track_asso_TruthC/I");
  outTree->Branch("jet_subleading_ak2track_asso_TruthL", &b_jet_subleading_ak2track_asso_TruthL, "jet_subleading_ak2track_asso_TruthL/I");

  outTree->Branch("jet_subleading_ak2track_trk_IBLHits", &b_jet_subleading_ak2track_trk_IBLHits);
  outTree->Branch("jet_subleading_ak2track_trk_nBLHits", &b_jet_subleading_ak2track_trk_nBLHits);
  outTree->Branch("jet_subleading_ak2track_trk_expectBLHits", &b_jet_subleading_ak2track_trk_expectBLHits);
  outTree->Branch("jet_subleading_ak2track_trk_nPixelHits", &b_jet_subleading_ak2track_trk_nPixelHits);
  outTree->Branch("jet_subleading_ak2track_trk_nPixelSharedHits", &b_jet_subleading_ak2track_trk_nPixelSharedHits);
  outTree->Branch("jet_subleading_ak2track_trk_nPixelSplitHits", &b_jet_subleading_ak2track_trk_nPixelSplitHits);
  outTree->Branch("jet_subleading_ak2track_trk_nSCTHits", &b_jet_subleading_ak2track_trk_nSCTHits);

  outTree->Branch("jet_leading_ak2track_asso_truth_pt", &b_jet_leading_ak2track_asso_truth_pt, "jet_leading_ak2track_asso_truth_pt/F");
  outTree->Branch("jet_leading_ak2track_asso_truth_eta", &b_jet_leading_ak2track_asso_truth_eta, "jet_leading_ak2track_asso_truth_eta/F");
  outTree->Branch("jet_leading_ak2track_asso_truth_phi", &b_jet_leading_ak2track_asso_truth_phi, "jet_leading_ak2track_asso_truth_phi/F");
  outTree->Branch("jet_leading_ak2track_asso_truth_m", &b_jet_leading_ak2track_asso_truth_m, "jet_leading_ak2track_asso_truth_m/F");
  outTree->Branch("jet_leading_ak2track_asso_truth_TruthB", &b_jet_leading_ak2track_asso_truth_TruthB, "jet_leading_ak2track_asso_truth_TruthB/I");
  outTree->Branch("jet_leading_ak2track_asso_truth_TruthC", &b_jet_leading_ak2track_asso_truth_TruthC, "jet_leading_ak2track_asso_truth_TruthC/I");
  outTree->Branch("jet_leading_ak2track_asso_truth_TruthL", &b_jet_leading_ak2track_asso_truth_TruthL, "jet_leading_ak2track_asso_truth_TruthL/I");

  outTree->Branch("jet_subleading_ak2track_asso_truth_pt", &b_jet_subleading_ak2track_asso_truth_pt, "jet_subleading_ak2track_asso_truth_pt/F");
  outTree->Branch("jet_subleading_ak2track_asso_truth_eta", &b_jet_subleading_ak2track_asso_truth_eta, "jet_subleading_ak2track_asso_truth_eta/F");
  outTree->Branch("jet_subleading_ak2track_asso_truth_phi", &b_jet_subleading_ak2track_asso_truth_phi, "jet_subleading_ak2track_asso_truth_phi/F");
  outTree->Branch("jet_subleading_ak2track_asso_truth_m", &b_jet_subleading_ak2track_asso_truth_m, "jet_subleading_ak2track_asso_truth_m/F");
  outTree->Branch("jet_subleading_ak2track_asso_truth_TruthB", &b_jet_subleading_ak2track_asso_truth_TruthB, "jet_subleading_ak2track_asso_truth_TruthB/I");
  outTree->Branch("jet_subleading_ak2track_asso_truth_TruthC", &b_jet_subleading_ak2track_asso_truth_TruthC, "jet_subleading_ak2track_asso_truth_TruthC/I");
  outTree->Branch("jet_subleading_ak2track_asso_truth_TruthL", &b_jet_subleading_ak2track_asso_truth_TruthL, "jet_subleading_ak2track_asso_truth_TruthL/I");


  outTree->Branch("boosted_bevent_sys", &b_boosted_bevent_sys);

  outTree->Branch("jet_leading_akt4_pt", &b_Leading_akt4_pt);
  outTree->Branch("jet_leading_akt4_eta", &b_Leading_akt4_eta);
  outTree->Branch("jet_leading_akt4_phi", &b_Leading_akt4_phi);
  outTree->Branch("jet_leading_akt4_m", &b_Leading_akt4_m);


  outTree->Branch("weight",    &m_weight,    "weight/F");
  outTree->Branch("weight_xs", &m_weight_xs, "weight_xs/F");
    

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode GluonBBMiniNTuple :: histInitialize () { return EL::StatusCode::SUCCESS; }
EL::StatusCode GluonBBMiniNTuple :: fileExecute    () { return EL::StatusCode::SUCCESS; }
EL::StatusCode GluonBBMiniNTuple :: changeInput    (bool /*firstFile*/) { return EL::StatusCode::SUCCESS; }

EL::StatusCode GluonBBMiniNTuple::execute ()
{

  const xAOD::EventInfo* eventInfo(0);
  RETURN_CHECK("GluonBBMiniNTuple::execute()", HelperFunctions::retrieve(eventInfo, "EventInfo", m_event, m_store), "");
  //
  // JZXW Pile-up Fluctuation Killer
  //
  if(m_isMC && (killPileupFluctuation() != 1)) return EL::StatusCode::SUCCESS;

  // 
  // For the cut flow
  //
  m_iCutflow = m_cutflowFirst;

  //
  // Set this first, as it's needed by passCut()
  //
  if(m_isMC){
    m_mcEventWeight = eventInfo->mcEventWeight();
  }  else{
    m_mcEventWeight = 1;
  }

  if (m_resolvedSysName.empty() && m_boostedSysName.empty()) {

    executeSingle("", "", true);

  //
  // Do Systematics
  //
  } else{

    //
    // Do Nominal
    //   (do the cut flow calculation here)
    //
    executeSingle("", "", true);

    //

    //
    // Do Boosted Systematics
    //
    if(!m_boostedSysName.empty()){

      // get vector of string giving the names
      vector<string>* boostedSystNames(nullptr);
      RETURN_CHECK("SelectDiJets::execute()", HelperFunctions::retrieve(boostedSystNames, m_boostedSysName, 0, m_store, m_debug) ,"");
    
      // Loop over resolved systematics 
      for ( string& boostSystName : *boostedSystNames ) {
      
		if(boostSystName.empty()) continue;

		if(m_debug) Info("execute",  "systName %s", boostSystName.c_str());

		executeSingle("", boostSystName, false);
	
      }

    }

  }
  
  return EL::StatusCode::SUCCESS;

}

EL::StatusCode GluonBBMiniNTuple::executeSingle(string resolvedSys, string boostedSys, bool countEvents) {
  if(m_debug) cout << " In executeSingle"  << resolvedSys << " " << boostedSys << endl;

  //
  //  Count All Events
  //
  if(countEvents) passCut(); 
  
  //
  // Start with stuff which is common to boosted and resolved
  //
  string syst = "";
  if(!boostedSys.empty())  syst = "Boosted_" +boostedSys;

  if( m_helpTree.find( syst ) == m_helpTree.end() ) { AddTree( syst ); }

  const xAOD::EventInfo* eventInfo(0);
  RETURN_CHECK("GluonBBMiniNTuple::execute()", HelperFunctions::retrieve(eventInfo, "EventInfo", m_event, m_store), "");
  
  const xAOD::VertexContainer* vertices(0);
  RETURN_CHECK("GluonBBMiniNTuple::execute()", HelperFunctions::retrieve(vertices, "PrimaryVertices", m_event, m_store), "");

  const xAOD::Vertex *pv = 0;
  pv = vertices->at( HelperFunctions::getPrimaryVertexLocation( vertices ) );

  //if(m_debug) cout << " Getting Muons"  << endl;
  //const xAOD::MuonContainer* muons(0);
  //RETURN_CHECK("GluonBBMiniNTuple::execute()", HelperFunctions::retrieve(muons, "Muons", m_event, m_store), "");

  const xAOD::MuonContainer* muons(nullptr);
  RETURN_CHECK("GluonBBMiniNTuple::execute()", HelperFunctions::retrieve(muons, m_muonContainerName, m_event, m_store), (m_muonContainerName).c_str());

  const xAOD::ParticleContainer* GluonBBCands(0);
  if(m_debug) cout << " Getting boosted H candidates: "  << m_boostedGluoncandName+boostedSys << endl;  
  RETURN_CHECK("GluonBBMiniNTuple::execute()", HelperFunctions::retrieve(GluonBBCands, m_boostedGluoncandName+boostedSys, m_event, m_store), "");  

  // Retrieve the container of resolved jets
  const xAOD::JetContainer* resolvedJets(0);
  if(m_debug) cout << " Getting resolved Jets: "  << m_resolvedJetsName+resolvedSys << endl;  
  RETURN_CHECK("GluonBBMiniNTuple::execute()", HelperFunctions::retrieve(resolvedJets, m_resolvedJetsName+resolvedSys, m_event, m_store), "");  

  const xAOD::TrackParticleContainer* InDetTrackParticles(0);
  RETURN_CHECK("GluonBBMiniNTuple::execute()", HelperFunctions::retrieve(InDetTrackParticles, "InDetTrackParticles", m_event, m_store), "");  

  //
  // Count+  the boosted jet kinematics
  //
  unsigned int nBoostedJetPassPtCut = 0;
  for(auto gluon : *GluonBBCands) {
    const xAOD::Jet* fatJet = gluon->auxdecor< const xAOD::Jet* >("caloJet");
    if(fatJet) {
      if(fatJet->pt() > m_FatJetPtSkimCut*1000) nBoostedJetPassPtCut++;
    }
  }
  //  only fill ntup if:
  //    - 2 resolved predijets and 2 bjets
  //    - 2 fatjets, 1 jet with pt > 300 GeV
  //    - 1 leptop, 1 predijets and 2 bjets
  //
  bool PassBoostedNGluonCands     = (GluonBBCands->size()     > 0);
  bool PassBoostedGluonCandPt     = (nBoostedJetPassPtCut      > 0);

  bool PassBoostedPreSel     = (PassBoostedNGluonCands     && PassBoostedGluonCandPt && m_doBoosted) || (m_ignoreRecoCuts);


  if(m_debug){ 
    cout << "Run/Event " << eventInfo->runNumber() << " / " << eventInfo->eventNumber() << endl;
  }

  if(m_debug){ 
    cout << " PassBoostedPreSel: "  << PassBoostedPreSel << endl;  
  }

  if (!PassBoostedPreSel) return EL::StatusCode::SUCCESS;
  if(m_debug) cout << " Pass GluonBBPresection: "  << endl;  
  
  //
  //  Count All Events
  //
  if(countEvents) passCut(); 

  //
  // Fill Event info
  //
  if(m_debug) cout << " Filling event " << endl;  
  m_helpTree[syst]->FillEvent( eventInfo );

  double xs            = wk()->metaData()->castDouble(SH::MetaFields::crossSection    ,1); 
  double filtEff       = wk()->metaData()->castDouble(SH::MetaFields::filterEfficiency,1); 

  m_weight_xs = xs * filtEff;
  m_weight    = m_mcEventWeight * xs * filtEff;

  //
  // Fill Trigger Info
  //
  if(m_debug) cout << " Filling trigger " << endl;  
  m_helpTree[syst]->FillTrigger( eventInfo );


  //
  // Fill Muon
  //
//  if(!m_muonContainerName.empty()){
//    if(m_debug) cout << " Filling muons " << endl;  
//    const xAOD::MuonContainer* muons(nullptr);
//    RETURN_CHECK("GluonBBMiniNTuple::execute()", HelperFunctions::retrieve(muons, m_muonContainerName, m_event, m_store), (m_muonContainerName).c_str());
//    m_helpTree[syst]->FillMuons(  muons, HelperFunctions::getPrimaryVertex( vertices ) );
//  }


  if(m_debug) cout << " Finished Filling Truth Particles"  << endl;  
   

  //
  // Fill Event Selection Flags
  //

  // Event selection flags
  if(m_debug) cout << " Filling Event Selection Flags"  << endl;  
  b_EventPass_GluonBBL1   = eventInfo->auxdecor<bool>("EventPass_GluonBBL1");
  b_EventPass_GluonBBHLT  = eventInfo->auxdecor<bool>("EventPass_GluonBBHLT");
  b_EventPass_GluonBBTrig = eventInfo->auxdecor<bool>("EventPass_GluonBBTrig");
  if (m_isMC){
    int mcChannelNumber = eventInfo->mcChannelNumber();

    // pythia
    if (mcChannelNumber == 361023){
      Total_Number_of_Events = 7349799;
    }
    if (mcChannelNumber == 361024){
      Total_Number_of_Events = 7975217;
    }
    if (mcChannelNumber == 361025){
      Total_Number_of_Events = 7977567;
    }
    if (mcChannelNumber == 361026){
      Total_Number_of_Events = 1893389;
    }
    if (mcChannelNumber == 361027){
      Total_Number_of_Events = 1770193;
    }
    if (mcChannelNumber == 361028){
      Total_Number_of_Events = 1743200;
    }

    // herwig
    if (mcChannelNumber == 426043){
      Total_Number_of_Events = 1784347;
    }
    if (mcChannelNumber == 426044){
      Total_Number_of_Events = 1925560;
    }
    if (mcChannelNumber == 426045){
      Total_Number_of_Events = 1926497;
    }
    if (mcChannelNumber == 426046){
      Total_Number_of_Events = 1930489;
    }
    if (mcChannelNumber == 426047){
      Total_Number_of_Events = 1927485;
    }
    if (mcChannelNumber == 426048){
      Total_Number_of_Events = 1928489;
    }

    // sherpa
    if (mcChannelNumber == 426133){
      Total_Number_of_Events = 1557533;
    }
    if (mcChannelNumber == 426134){
      Total_Number_of_Events = 1978198;
    }
    if (mcChannelNumber == 426135){
      Total_Number_of_Events = 1961488;
    }
    if (mcChannelNumber == 426136){
      Total_Number_of_Events = 1963794;
    }
    if (mcChannelNumber == 426137){
      Total_Number_of_Events = 954993;
    }
    if (mcChannelNumber == 426138){
      Total_Number_of_Events = 977392;
    }

    // ggF
    if (mcChannelNumber == 345342){
      Total_Number_of_Events = 55449982.3789;
    }
    
  }

  
  if(m_debug) cout << " Finished Filling Event Selection Flags"  << endl;  


  //
  // Fill Boosted Gluon Candidates
  //
  if(m_debug) cout << " Filling Boosted Gluon Candidates"  << endl;  
  b_gluon_boosted_n = 0;
  b_gluon_boosted_pt=-999;
  b_gluon_boosted_eta=-999;
  b_gluon_boosted_phi=-999;
  b_gluon_boosted_m=-999;
  b_gluon_boosted_D2=-999;
  b_gluon_boosted_Tau21=-999;
  b_gluon_boosted_Tau21WTA=-999;
  b_gluon_boosted_nTrack=-999;
  b_gluon_boosted_nTruthB=-999;
  b_gluon_boosted_nTruthC=-999;

  b_gluon_boosted_truth_n=0;
  b_gluon_boosted_truth_pt=-999;
  b_gluon_boosted_truth_eta=-999;
  b_gluon_boosted_truth_phi=-999;
  b_gluon_boosted_truth_m=-999;


  b_jet_ak2track_asso_n.clear();
  b_jet_ak2track_asso_n_addl.clear();
  b_jet_ak2track_asso_pt.clear();
  b_jet_ak2track_asso_eta.clear();
  b_jet_ak2track_asso_phi.clear();
  b_jet_ak2track_asso_m.clear();
  b_jet_ak2track_asso_MV2c00.clear();
  b_jet_ak2track_asso_MV2c10.clear();
  b_jet_ak2track_asso_MV2c20.clear();
  b_jet_ak2track_asso_MV2c100.clear();
  b_jet_ak2track_asso_sys.clear();
  b_jet_ak2track_asso_sysname.clear();
  b_jet_ak2track_asso_TruthFlav.clear();
  b_jet_ak2track_asso_TruthB.clear();
  b_jet_ak2track_asso_TruthC.clear();
  b_jet_ak2track_asso_TruthL.clear();

  b_jet_leading_ak2track_asso_pt=-999;
  b_jet_leading_ak2track_asso_eta=-999;
  b_jet_leading_ak2track_asso_phi=-999;
  b_jet_leading_ak2track_asso_m=-999;
  b_jet_leading_ak2track_asso_Sd0=-999;
  b_jet_leading_ak2track_asso_Sd0_Sub=-999;
  b_jet_leading_ak2track_asso_MV2c10=-999;
  b_jet_leading_ak2track_asso_MV2c20=-999;
  b_jet_leading_ak2track_asso_IP2D=-999;
  b_jet_leading_ak2track_asso_IP3D=-999;
  b_jet_leading_ak2track_asso_sys.clear();
  b_jet_leading_ak2track_asso_sysname.clear();
  b_jet_leading_ak2track_asso_TruthFlav=-999;
  b_jet_leading_ak2track_asso_TruthB=-999;
  b_jet_leading_ak2track_asso_TruthC=-999;
  b_jet_leading_ak2track_asso_TruthL=-999;
  b_jet_leading_ak2track_asso_nMuon=-999;
  b_jet_leading_ak2track_ntrk = -999;

  b_jet_leading_ak2track_trk_IBLHits.clear();
  b_jet_leading_ak2track_trk_nBLHits.clear();
  b_jet_leading_ak2track_trk_expectBLHits.clear();
  b_jet_leading_ak2track_trk_nPixelHits.clear();
  b_jet_leading_ak2track_trk_nPixelSharedHits.clear();
  b_jet_leading_ak2track_trk_nPixelSplitHits.clear();
  b_jet_leading_ak2track_trk_nSCTHits.clear();
  b_jet_leading_ak2track_trk_pT.clear();
  b_jet_leading_ak2track_trk_eta.clear();


  b_jet_subleading_ak2track_asso_pt=-999;
  b_jet_subleading_ak2track_asso_eta=-999;
  b_jet_subleading_ak2track_asso_phi=-999;
  b_jet_subleading_ak2track_asso_m=-999;
  b_jet_subleading_ak2track_asso_Sd0=-999;
  b_jet_subleading_ak2track_asso_Sd0_Sub=-999;
  b_jet_subleading_ak2track_asso_MV2c10=-999;
  b_jet_subleading_ak2track_asso_MV2c20=-999;
  b_jet_subleading_ak2track_asso_IP2D=-999;
  b_jet_subleading_ak2track_asso_IP3D=-999;
  b_jet_subleading_ak2track_asso_sys.clear();
  b_jet_subleading_ak2track_asso_sysname.clear();
  b_jet_subleading_ak2track_asso_TruthFlav=-999;
  b_jet_subleading_ak2track_asso_TruthB=-999;
  b_jet_subleading_ak2track_asso_TruthC=-999;
  b_jet_subleading_ak2track_asso_TruthL=-999;
  b_jet_subleading_ak2track_asso_nMuon=-999;

  b_jet_subleading_ak2track_trk_IBLHits.clear();
  b_jet_subleading_ak2track_trk_nBLHits.clear();
  b_jet_subleading_ak2track_trk_expectBLHits.clear();
  b_jet_subleading_ak2track_trk_nPixelHits.clear();
  b_jet_subleading_ak2track_trk_nPixelSharedHits.clear();
  b_jet_subleading_ak2track_trk_nPixelSplitHits.clear();
  b_jet_subleading_ak2track_trk_nSCTHits.clear();

  b_jet_leading_ak2track_asso_truth_pt=-999;
  b_jet_leading_ak2track_asso_truth_eta=-999;
  b_jet_leading_ak2track_asso_truth_phi=-999;
  b_jet_leading_ak2track_asso_truth_m=-999;
  b_jet_leading_ak2track_asso_truth_TruthFlav=-999;
  b_jet_leading_ak2track_asso_truth_TruthL=-999;
  b_jet_leading_ak2track_asso_truth_TruthB=-999;
  b_jet_leading_ak2track_asso_truth_TruthC=-999;

  b_jet_subleading_ak2track_asso_truth_pt=-999;
  b_jet_subleading_ak2track_asso_truth_eta=-999;
  b_jet_subleading_ak2track_asso_truth_phi=-999;
  b_jet_subleading_ak2track_asso_truth_m=-999;
  b_jet_subleading_ak2track_asso_truth_TruthFlav=-999;
  b_jet_subleading_ak2track_asso_truth_TruthL=-999;
  b_jet_subleading_ak2track_asso_truth_TruthB=-999;
  b_jet_subleading_ak2track_asso_truth_TruthC=-999;


  b_boosted_bevent_sys.clear();

  b_Leading_akt4_pt=-99;
  b_Leading_akt4_eta=-99;
  b_Leading_akt4_phi=-99;
  b_Leading_akt4_m=-99;

  m_helpTree[syst]->ClearJets("boostedJets");

  // loop over boosted Higgs candidates
  for(auto gluon : *GluonBBCands) {  

    const xAOD::Jet* fatJet = gluon->auxdecor< const xAOD::Jet* >("caloJet");
    std::vector<const xAOD::Jet*> assotrkjets_fatJet = gluon->auxdecor< std::vector<const xAOD::Jet*> >("allTrkJet");

    if(fatJet) {
      //fill the jet information first    
      m_helpTree[syst]->FillJet(fatJet, pv, HelperFunctions::getPrimaryVertexLocation( vertices ), "boostedJets");
      //fill the boosted information

      b_gluon_boosted_pt=( fatJet->pt() );
      b_gluon_boosted_eta=( fatJet->eta() );
      b_gluon_boosted_phi=(fatJet->phi() );
      b_gluon_boosted_m=( fatJet->m() );

      b_gluon_boosted_D2=( getD2(fatJet) );
      b_gluon_boosted_Tau21=( getTau21(fatJet, false) );
      b_gluon_boosted_Tau21WTA=( getTau21(fatJet, true) );
      b_gluon_boosted_nTrack=( fatJet->auxdata<int>("GhostTrackCount") );

      if (m_isMC){
	std::vector<const xAOD::TruthParticle*> BHadrons = fatJet->getAssociatedObjects<xAOD::TruthParticle>("GhostBHadronsFinal");
	std::vector<const xAOD::TruthParticle*> CHadrons = fatJet->getAssociatedObjects<xAOD::TruthParticle>("GhostCHadronsFinal");
	b_gluon_boosted_nTruthB=( BHadrons.size());
	b_gluon_boosted_nTruthC=( CHadrons.size());
      }

      const xAOD::Jet* fatJetParentJet = 0;
      try{
        auto el = fatJet->auxdata<ElementLink<xAOD::JetContainer> >("Parent");
        if(!el.isValid()){
          Warning("executeSingle()", "Invalid link to \"Parent\" from fat-jet.");
        }
        else{
          fatJetParentJet = (*el);
        }
      }catch(...){
        Warning("executeSingle()", "Unable to get parent jet of fat-jet for truth labeling. Trimmed jet area would be used!");
        fatJetParentJet = fatJet;
      }


      int ak2track_asso_n=0;
      std::vector<float> ak2track_asso_pt;
      std::vector<float> ak2track_asso_eta;
      std::vector<float> ak2track_asso_phi;
      std::vector<float> ak2track_asso_m;
      std::vector<float> ak2track_asso_MV2c00;
      std::vector<float> ak2track_asso_MV2c10;
      std::vector<float> ak2track_asso_MV2c20;
      std::vector<float> ak2track_asso_MV2c100;
      std::vector<float> ak2track_asso_IP3D;
      std::vector<float> ak2track_asso_IP2D;
      std::vector<vector<float> > ak2track_asso_sys;
      std::vector<int>  ak2track_asso_TruthFlav;
      std::vector<int>  ak2track_asso_TruthB;
      std::vector<int>  ak2track_asso_TruthC;
      std::vector<int>  ak2track_asso_TruthL;
      b_boosted_bevent_sys = eventInfo->auxdecor<std::vector<float>>("BTag_SF_FixedCutBEff_77_GLOBAL");

      //retrive the boosted b-tag systematic names
      vector<std::string>* boostedbtagSystNames(nullptr);
      RETURN_CHECK("SelectDiJets::execute()", HelperFunctions::retrieve(boostedbtagSystNames, "FTSys_FixedCutBEff_77", 0, m_store, m_debug) ,"");
      b_jet_ak2track_asso_sysname = *boostedbtagSystNames;
      b_jet_leading_ak2track_asso_sysname = *boostedbtagSystNames;
      b_jet_subleading_ak2track_asso_sysname = *boostedbtagSystNames;

      std::vector<float> dummy_sys_vec;
      
      int ntrkjet = assotrkjets_fatJet.size();
      b_gluon_boosted_n = ntrkjet;

      if (ntrkjet ==0){
	b_jet_leading_ak2track_asso_pt=(-999);
	b_jet_leading_ak2track_asso_eta=(-999);
	b_jet_leading_ak2track_asso_phi=(-999);
	b_jet_leading_ak2track_asso_m=(-999);
	b_jet_leading_ak2track_asso_MV2c10=(-999);
	b_jet_leading_ak2track_asso_MV2c20=(-999);
	b_jet_leading_ak2track_asso_Sd0=(-999);
	b_jet_leading_ak2track_asso_Sd0_Sub=(-999);
	b_jet_leading_ak2track_asso_IP3D=(-999);
	b_jet_leading_ak2track_asso_IP2D=(-999);
	b_jet_leading_ak2track_asso_TruthFlav=(-999);
	b_jet_leading_ak2track_asso_TruthB=(-999);
	b_jet_leading_ak2track_asso_TruthC=(-999);
	b_jet_leading_ak2track_asso_TruthL=(-999);
	b_jet_leading_ak2track_asso_nMuon=(-999);
	b_jet_leading_ak2track_asso_sys.push_back(dummy_sys_vec);

	b_jet_subleading_ak2track_asso_pt=(-999);
	b_jet_subleading_ak2track_asso_eta=(-999);
	b_jet_subleading_ak2track_asso_phi=(-999);
	b_jet_subleading_ak2track_asso_m=(-999);
	b_jet_subleading_ak2track_asso_Sd0=(-999);
	b_jet_subleading_ak2track_asso_Sd0_Sub=(-999);
	b_jet_subleading_ak2track_asso_MV2c10=(-999);
	b_jet_subleading_ak2track_asso_MV2c20=(-999);
	b_jet_subleading_ak2track_asso_IP3D=(-999);
	b_jet_subleading_ak2track_asso_IP2D=(-999);

	b_jet_subleading_ak2track_asso_TruthFlav=(-999);
	b_jet_subleading_ak2track_asso_TruthB=(-999);
	b_jet_subleading_ak2track_asso_TruthC=(-999);
	b_jet_subleading_ak2track_asso_TruthL=(-999);
	b_jet_subleading_ak2track_asso_nMuon=(-999);
	b_jet_subleading_ak2track_asso_sys.push_back(dummy_sys_vec);
      }

      if (ntrkjet ==1){
	b_jet_subleading_ak2track_asso_pt=(-999);
	b_jet_subleading_ak2track_asso_eta=(-999);
	b_jet_subleading_ak2track_asso_phi=(-999);
	b_jet_subleading_ak2track_asso_m=(-999);
	b_jet_subleading_ak2track_asso_Sd0=(-999);
	b_jet_subleading_ak2track_asso_Sd0_Sub=(-999);
	b_jet_subleading_ak2track_asso_MV2c10=(-999);
	b_jet_subleading_ak2track_asso_MV2c20=(-999);
	b_jet_subleading_ak2track_asso_IP2D=(-999);
	b_jet_subleading_ak2track_asso_IP3D=(-999);

	b_jet_subleading_ak2track_asso_TruthFlav=(-999);
	b_jet_subleading_ak2track_asso_TruthB=(-999);
	b_jet_subleading_ak2track_asso_TruthC=(-999);
	b_jet_subleading_ak2track_asso_TruthL=(-999);
	b_jet_subleading_ak2track_asso_sys.push_back(dummy_sys_vec);
	b_jet_subleading_ak2track_asso_nMuon=(-999);
      }

      int itrkjet = 0;

      std::vector<int> JetMuonMatchVec = GluonBBMiniNTuple::JetMuondRMatch(assotrkjets_fatJet, muons, 0.2);
      
      for(auto fatJetTrackJet : assotrkjets_fatJet){
	itrkjet++;

	const xAOD::BTagging* BTag = fatJetTrackJet->btagging();

	float maxSd0=0;
	float secmaxSd0=-1;
	float LeadSd0=-99;
	float SubSd0=-99;

	xAOD::JetConstituentVector consVec = fatJetTrackJet->getConstituents();
	b_jet_leading_ak2track_ntrk = consVec.size();

	try{
	  for (auto cons : consVec){
	    const xAOD::TrackParticle* track = dynamic_cast <const xAOD::TrackParticle*> (cons->rawConstituent());

	    xAOD::TrackParticle track_smear = *track;
	    m_InDetTrackSmearingTool->applyCorrection( track_smear);

	    double d0sig = xAOD::TrackingHelpers::d0significance( track, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), eventInfo->beamPosSigmaXY());
	    double sign_rphi = sin( fatJetTrackJet->p4().Phi() - track->p4().Phi() ) * track->d0();
	    int sign = (sign_rphi>0 ? 1 : -1);
	    d0sig = sign* fabs(d0sig);

	    double d0sig_corr = xAOD::TrackingHelpers::d0significance( &track_smear, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), eventInfo->beamPosSigmaXY());
	    double sign_rphi_smear = sin( fatJetTrackJet->p4().Phi() - track_smear.p4().Phi() ) * track_smear.d0();
	    int sign_smear = (sign_rphi_smear>0 ? 1 : -1);
	    d0sig_corr = sign_smear* fabs(d0sig_corr);

	    if (not m_isMC){
	      d0sig_corr = d0sig;
	    }

	    if ( fabs(d0sig_corr)>maxSd0){
	      LeadSd0 = d0sig_corr;
	      maxSd0 = fabs(d0sig_corr);
	    }
	    if ( fabs(d0sig_corr)<maxSd0 and fabs(d0sig_corr)>secmaxSd0){
	      SubSd0 = d0sig_corr;
	      secmaxSd0 = fabs(d0sig_corr);
	    }

	    uint8_t numberOfIBLHits;
	    track->summaryValue(numberOfIBLHits,xAOD::numberOfInnermostPixelLayerHits);

	    uint8_t expectBLayerHits;
	    track->summaryValue(expectBLayerHits,xAOD::expectInnermostPixelLayerHit);

	    uint8_t numberOfBLHits;
	    track->summaryValue(numberOfBLHits, xAOD::numberOfNextToInnermostPixelLayerHits);

	    uint8_t numberOfPixelHits;
	    track->summaryValue(numberOfPixelHits,xAOD::numberOfPixelHits);
	    
	    uint8_t numberOfPixelSharedHits;
	    track->summaryValue(numberOfPixelSharedHits,xAOD::numberOfPixelSharedHits);

	    uint8_t numberOfPixelSplitHits;
	    track->summaryValue(numberOfPixelSplitHits,xAOD::numberOfPixelSharedHits);

	    uint8_t numberOfSCTHits;
	    track->summaryValue(numberOfSCTHits,xAOD::numberOfSCTHits);
	    
	    if (itrkjet ==1){
	      b_jet_leading_ak2track_trk_IBLHits.push_back( int(numberOfIBLHits) );
	      b_jet_leading_ak2track_trk_nBLHits.push_back( int(numberOfBLHits) );
	      b_jet_leading_ak2track_trk_expectBLHits.push_back( int(expectBLayerHits) );
	      b_jet_leading_ak2track_trk_nPixelHits.push_back( int(numberOfPixelHits));
	      b_jet_leading_ak2track_trk_nPixelSharedHits.push_back( int(numberOfPixelSharedHits));
	      b_jet_leading_ak2track_trk_nPixelSplitHits.push_back( int(numberOfPixelSplitHits));
	      b_jet_leading_ak2track_trk_nSCTHits.push_back( int(numberOfSCTHits));
	      b_jet_leading_ak2track_trk_pT.push_back(track->p4().Pt());
	      b_jet_leading_ak2track_trk_eta.push_back(track->p4().Eta());
	    }
	    if (itrkjet ==2){
	      b_jet_subleading_ak2track_trk_IBLHits.push_back( int(numberOfIBLHits) );
	      b_jet_subleading_ak2track_trk_nBLHits.push_back( int(numberOfBLHits) );
	      b_jet_subleading_ak2track_trk_expectBLHits.push_back( int(expectBLayerHits) );
	      b_jet_subleading_ak2track_trk_nPixelHits.push_back( int(numberOfPixelHits));
	      b_jet_subleading_ak2track_trk_nPixelSharedHits.push_back( int(numberOfPixelSharedHits));
	      b_jet_subleading_ak2track_trk_nPixelSplitHits.push_back( int(numberOfPixelSplitHits));
	      b_jet_subleading_ak2track_trk_nSCTHits.push_back( int(numberOfSCTHits));
	    }
	  }
	}
	catch(...){
	  Warning("executeSingle()", "element link problem for track converstion from iparticle");
	}

        ak2track_asso_n++; 
        ak2track_asso_pt.push_back( fatJetTrackJet->pt() );
        ak2track_asso_eta.push_back( fatJetTrackJet->eta() );
        ak2track_asso_phi.push_back( fatJetTrackJet->phi() );
        ak2track_asso_m.push_back( fatJetTrackJet->m() );
        ak2track_asso_MV2c00.push_back( MV2(fatJetTrackJet, "MV2c00") );
        ak2track_asso_MV2c10.push_back( MV2(fatJetTrackJet, "MV2c10") );
        ak2track_asso_MV2c20.push_back( MV2(fatJetTrackJet, "MV2c20") );
        ak2track_asso_MV2c100.push_back( MV2(fatJetTrackJet, "MV2c100") );
        ak2track_asso_sys.push_back(fatJetTrackJet->auxdecor<std::vector<float> >("BTag_SF_FixedCutBEff_77"));

	if (itrkjet==1){
	  b_jet_leading_ak2track_asso_pt=( fatJetTrackJet->pt() );
	  b_jet_leading_ak2track_asso_Sd0=( LeadSd0);
	  b_jet_leading_ak2track_asso_Sd0_Sub=( SubSd0);
	  b_jet_leading_ak2track_asso_eta=( fatJetTrackJet->eta() );
	  b_jet_leading_ak2track_asso_phi=( fatJetTrackJet->phi() );
	  b_jet_leading_ak2track_asso_m=( fatJetTrackJet->m() );
	  b_jet_leading_ak2track_asso_MV2c10=( MV2(fatJetTrackJet, "MV2c10") );
	  b_jet_leading_ak2track_asso_MV2c20=( MV2(fatJetTrackJet, "MV2c20") );
	  b_jet_leading_ak2track_asso_IP3D=( BTAGIP(fatJetTrackJet, "IP3D") );
	  b_jet_leading_ak2track_asso_IP2D=( BTAGIP(fatJetTrackJet, "IP2D") );
	  b_jet_leading_ak2track_asso_nMuon=( JetMuonMatchVec[0]);
	  b_jet_leading_ak2track_asso_sys.push_back(fatJetTrackJet->auxdecor<std::vector<float> >("BTag_SF_FixedCutBEff_77"));
	}

	if (itrkjet==2){
	  b_jet_subleading_ak2track_asso_pt=( fatJetTrackJet->pt() );
	  b_jet_subleading_ak2track_asso_Sd0=( LeadSd0);
	  b_jet_subleading_ak2track_asso_Sd0_Sub=( SubSd0);
	  b_jet_subleading_ak2track_asso_eta=( fatJetTrackJet->eta() );
	  b_jet_subleading_ak2track_asso_phi=( fatJetTrackJet->phi() );
	  b_jet_subleading_ak2track_asso_m=( fatJetTrackJet->m() );
	  b_jet_subleading_ak2track_asso_MV2c10=( MV2(fatJetTrackJet, "MV2c10") );
	  b_jet_subleading_ak2track_asso_MV2c20=( MV2(fatJetTrackJet, "MV2c20") );
	  b_jet_subleading_ak2track_asso_IP3D=( BTAGIP(fatJetTrackJet, "IP3D") );
	  b_jet_subleading_ak2track_asso_IP2D=( BTAGIP(fatJetTrackJet, "IP2D") );
	  b_jet_subleading_ak2track_asso_nMuon=( JetMuonMatchVec[1]);
	  b_jet_subleading_ak2track_asso_sys.push_back(fatJetTrackJet->auxdecor<std::vector<float> >("BTag_SF_FixedCutBEff_77"));
	}

	if (m_isMC){
	  std::vector<const xAOD::TruthParticle*> BHadrons = fatJetTrackJet->getAssociatedObjects<xAOD::TruthParticle>("ConeExclBHadronsFinal");
	  std::vector<const xAOD::TruthParticle*> CHadrons = fatJetTrackJet->getAssociatedObjects<xAOD::TruthParticle>("ConeExclCHadronsFinal");
	  if ( BHadrons.size()>0){
	    ak2track_asso_TruthFlav.push_back(5);
	    ak2track_asso_TruthB.push_back(1);
	    ak2track_asso_TruthC.push_back(0);
	    ak2track_asso_TruthL.push_back(0);
	    if (itrkjet ==1){
	      b_jet_leading_ak2track_asso_TruthFlav=(5);
	      b_jet_leading_ak2track_asso_TruthB=(1);
	      b_jet_leading_ak2track_asso_TruthC=(0);
	      b_jet_leading_ak2track_asso_TruthL=(0);
	    }
	    if (itrkjet ==2){
	      b_jet_subleading_ak2track_asso_TruthFlav=(5);
	      b_jet_subleading_ak2track_asso_TruthB=(1);
	      b_jet_subleading_ak2track_asso_TruthC=(0);
	      b_jet_subleading_ak2track_asso_TruthL=(0);
	    }
	  }

	  else if ( CHadrons.size()>0){
	    ak2track_asso_TruthFlav.push_back(4);
	    ak2track_asso_TruthB.push_back(0);
	    ak2track_asso_TruthC.push_back(1);
	    ak2track_asso_TruthL.push_back(0);
	    if (itrkjet ==1){
	      b_jet_leading_ak2track_asso_TruthFlav=(4);
	      b_jet_leading_ak2track_asso_TruthB=(0);
	      b_jet_leading_ak2track_asso_TruthC=(1);
	      b_jet_leading_ak2track_asso_TruthL=(0);
	    }
	    if (itrkjet ==2){
	      b_jet_subleading_ak2track_asso_TruthFlav=(4);
	      b_jet_subleading_ak2track_asso_TruthB=(0);
	      b_jet_subleading_ak2track_asso_TruthC=(1);
	      b_jet_subleading_ak2track_asso_TruthL=(0);
	    }
	  }
	  else{
	    ak2track_asso_TruthFlav.push_back(1);
	    ak2track_asso_TruthB.push_back(0);
	    ak2track_asso_TruthC.push_back(0);
	    ak2track_asso_TruthL.push_back(1);
	    if (itrkjet ==1){
	      b_jet_leading_ak2track_asso_TruthFlav=(1);
	      b_jet_leading_ak2track_asso_TruthB=(0);
	      b_jet_leading_ak2track_asso_TruthC=(0);
	      b_jet_leading_ak2track_asso_TruthL=(1);
	    }
	    if (itrkjet ==2){
	      b_jet_subleading_ak2track_asso_TruthFlav=(1);
	      b_jet_subleading_ak2track_asso_TruthB=(0);
	      b_jet_subleading_ak2track_asso_TruthC=(0);
	      b_jet_subleading_ak2track_asso_TruthL=(1);
	    }
	  }
	}
      }
      //save all the vectors 
      b_jet_ak2track_asso_n.push_back( ak2track_asso_n );
      b_jet_ak2track_asso_n_addl.push_back( gluon->auxdecor< int > ("ak2track_asso_n") - ak2track_asso_n );
      b_jet_ak2track_asso_pt.push_back( ak2track_asso_pt );
      b_jet_ak2track_asso_eta.push_back( ak2track_asso_eta );
      b_jet_ak2track_asso_phi.push_back( ak2track_asso_phi );
      b_jet_ak2track_asso_m.push_back( ak2track_asso_m );
      b_jet_ak2track_asso_MV2c00.push_back( ak2track_asso_MV2c00 );
      b_jet_ak2track_asso_MV2c10.push_back( ak2track_asso_MV2c10 );
      b_jet_ak2track_asso_MV2c20.push_back( ak2track_asso_MV2c20 );
      b_jet_ak2track_asso_MV2c100.push_back( ak2track_asso_MV2c100 );
      b_jet_ak2track_asso_sys.push_back( ak2track_asso_sys );
      b_jet_ak2track_asso_TruthFlav.push_back(ak2track_asso_TruthFlav);
      b_jet_ak2track_asso_TruthB.push_back(ak2track_asso_TruthB);
      b_jet_ak2track_asso_TruthC.push_back(ak2track_asso_TruthC);
      b_jet_ak2track_asso_TruthL.push_back(ak2track_asso_TruthL);
    }
    else{
      cout<<"WARNING : NO FAT-JET ASSOCIATED TO A BOOSTED HIGGS CANDIDATE!"<<endl;
    }
  }
  if(m_debug) cout << " Finished Filling Boosted Gluon Candidates"  << endl;  


  // Fill resolved jets/ leading anti-kt 4 jets for pT slices combining only
  if(!m_resolvedJetsName.empty()){
    if(m_debug) cout << " Filling resolved jets " << endl;  
    for(auto jet_itr: *resolvedJets){
      b_Leading_akt4_pt = jet_itr->pt();
      b_Leading_akt4_eta = jet_itr->eta();
      b_Leading_akt4_phi = jet_itr->phi();
      b_Leading_akt4_m = jet_itr->m();
      break;
    }
  }

  //// get charged truth particles
  if (m_isMC){

    std::vector<xAOD::TruthParticle> ChargedTruthParticle;

    const xAOD::TruthParticleContainer* truth_particles_unsort(0);
    xAOD::TruthParticleContainer* truth_particles(0);

    if (!m_inTruthParticleName.empty() && m_isMC){
      if(m_debug) cout << " Getting Truth Particles"  << endl;  
      RETURN_CHECK("GluonBBMiniNTuple::execute()", HelperFunctions::retrieve(truth_particles_unsort, m_inTruthParticleName, m_event, m_store), m_inTruthParticleName.c_str());
      truth_particles = new xAOD::TruthParticleContainer(*truth_particles_unsort);//memory leak

      for(unsigned int i=0; i<truth_particles->size(); i++){
	const xAOD::TruthParticle* curr_truth = truth_particles->at(i);
	if(!curr_truth) continue;
	int status = curr_truth->status();
	double charge = curr_truth->charge();
	double pt = curr_truth->pt();
	double eta = curr_truth->eta();
	int barcode = curr_truth->barcode();
	
	// apply selection according to https://twiki.cern.ch/twiki/bin/view/AtlasProtected/TrackingCPICHEP2016#Truth_definitions
	// we apply phase space and primary particles
	if (status != 1 ) continue;
	if (charge == 0 ) continue;
	if (barcode <0 or barcode>2e5) continue;
	if (pt <500 or fabs(eta)>2.5) continue;
	
	ChargedTruthParticle.push_back( *curr_truth);
      }
    }

    // cluster truth particles to track jets
    std::vector<TLorentzVector> TruthTrackJets4Vector = HelperFunctions::TruthTrackJetReclustering(ChargedTruthParticle, 0.2, 0, fastjet::antikt_algorithm);
    // pre-select to track jets
    std::vector<TLorentzVector> TruthTrackJets4Vector_Preselected;
    for(auto trkjet : TruthTrackJets4Vector){
      if (trkjet.Pt()<10) continue;
      if (fabs(trkjet.Eta())>2.5) continue;
      TruthTrackJets4Vector_Preselected.push_back(trkjet);
    }


    std::vector<int> TruthTrackJets4Vector_BMatched;
    int nB_matched_truth_trackjet=0;

    // loop over track jets to dR match truth b-hadron
    for(auto trkjet : TruthTrackJets4Vector_Preselected){
      bool matched = false;
      for (auto truthpart : *truth_particles_unsort){
	double dR = sqrt( pow(truthpart->eta()-trkjet.Eta(), 2) +pow(truthpart->phi()-trkjet.Phi(), 2));

	if ( IsWeaklyDecayingBHadron(truthpart->pdgId()) and dR<0.3){
	  matched = true;
	  nB_matched_truth_trackjet++;
	}
      }
      if (matched) TruthTrackJets4Vector_BMatched.push_back(1);
      else TruthTrackJets4Vector_BMatched.push_back(0);
    }
  
    bool PassTruthGbbSel = false;
    PassTruthGbbSel = GhostMatchTruthTrkJetToLargeRJet( TruthTrackJets4Vector_Preselected, TruthTrackJets4Vector_BMatched);
  
    if (!PassTruthGbbSel) return EL::StatusCode::SUCCESS;
    if(m_debug) cout << " Pass TruthGluonBBPresection: "  << endl;  
  }

  // fill the tree
  if(m_debug) cout << " Filling Tree"  << endl;
  m_helpTree[syst]->Fill();

  if(m_debug) cout << " Filled Tree"  << endl;

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode GluonBBMiniNTuple :: postExecute () { return EL::StatusCode::SUCCESS; }

EL::StatusCode GluonBBMiniNTuple :: finalize () {

  Info("finalize()", m_name.c_str());

  if (!m_helpTree.empty()){
    for( auto tree : m_helpTree) {
      if (tree.second) delete tree.second;
    }
  }

  return EL::StatusCode::SUCCESS;
}

 

EL::StatusCode GluonBBMiniNTuple :: histFinalize ()
{
  TFile * treeFile = wk()->getOutputFile( m_name.c_str() );
  
  //
  // Write out the cut flow histograms
  //
  TH1F* thisCutflowHist = (TH1F*) m_cutflowHist->Clone();
  std::string thisName = thisCutflowHist->GetName();
  thisCutflowHist->SetName( (thisName+"_GluonBBMiniNTuple").c_str() );
  thisCutflowHist->SetDirectory( treeFile );
  
  TH1F* thisCutflowHistW = (TH1F*) m_cutflowHistW->Clone();
  thisName = thisCutflowHistW->GetName();
  thisCutflowHistW->SetName( (thisName+"_GluonBBMiniNTuple").c_str() );
  thisCutflowHistW->SetDirectory( treeFile );
  
  //
  // Write out the Metadata
  //
  TFile *fileMD = wk()->getOutputFile("metadata");
  TH1D* m_histEventCount   = (TH1D*)fileMD->Get("MetaData_EventCount");
  TH1F* thisHistEventCount = (TH1F*) m_histEventCount->Clone();
  thisName = thisHistEventCount->GetName();
  thisHistEventCount->SetName( (thisName+"_GluonBBMiniNTuple").c_str() );
  thisHistEventCount->SetDirectory( treeFile );

  //
  // Write out the jet multiplicities
  //
//  TH1* m_nJet  = wk()->getOutputHist("GluonBBCutFlow/nJet");
//  TH1F* nJetClone = (TH1F*) m_nJet->Clone();
//  thisName = nJetClone->GetName();
//  nJetClone->SetName( (thisName).c_str() );
//  nJetClone->SetDirectory( treeFile );
//
//  TH1* m_nBJet  = wk()->getOutputHist("GluonBBCutFlow/nBJet");
//  TH1F* nBJetClone = (TH1F*) m_nBJet->Clone();
//  thisName = nBJetClone->GetName();
//  nBJetClone->SetName( (thisName ).c_str() );
//  nBJetClone->SetDirectory( treeFile );
  
  return EL::StatusCode::SUCCESS;
}

//
// Easy method for automatically filling cutflow and incrementing counter
//
void GluonBBMiniNTuple::passCut(){
  m_cutflowHist ->Fill(m_iCutflow, 1);
  m_cutflowHistW->Fill(m_iCutflow, m_mcEventWeight);
  m_iCutflow++;
}

int GluonBBMiniNTuple::killPileupFluctuation(){
  const xAOD::EventInfo* eventInfo(0);
  RETURN_CHECK("GluonBBMiniNTuple::killPileupFluctuation()", HelperFunctions::retrieve(eventInfo, "EventInfo", m_event, m_store), "");

  const xAOD::JetContainer* AntiKt4TruthJets(0);
  RETURN_CHECK("GluonBBMiniNTuple::killPileupFluctuation()", HelperFunctions::retrieve(AntiKt4TruthJets, "AntiKt4TruthJets", m_event, m_store), "");

  const xAOD::JetContainer* AntiKt4EMTopoJets(0);
  RETURN_CHECK("GluonBBMiniNTuple::killPileupFluctuation()", HelperFunctions::retrieve(AntiKt4EMTopoJets, "AntiKt4EMTopoJets_Calib", m_event, m_store), "");

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // all jet collections are pre-assumed to be sorted by pT

  int mcChannelNumber = eventInfo->mcChannelNumber();

  if ((mcChannelNumber >= 361020) && (mcChannelNumber <= 361032)){
    if((AntiKt4TruthJets->size() >= 1) && (AntiKt4EMTopoJets->size() >= 2)){
      float ptave = (AntiKt4EMTopoJets->at(0)->pt() + AntiKt4EMTopoJets->at(1)->pt())/2.;
      float pttruth = AntiKt4TruthJets->at(0)->pt();
      float ratio = ptave/pttruth;

      if(ratio > 1.4) return 0;
    }
  }

  return 1;
}

double GluonBBMiniNTuple::getD2(const xAOD::Jet* jet){
  if(jet->getAttribute<double>("ECF2") == 0.){
    Warning("getD2()", "ECF2 returns 0. -1 will be returned");
    return -1.;
  }

  return (jet->getAttribute<double>("ECF3") * TMath::Power(jet->getAttribute<double>("ECF1"), 3) )/(TMath::Power(jet->getAttribute<double>("ECF2"), 3));
}

double GluonBBMiniNTuple::getTau21(const xAOD::Jet* jet, bool doWTA){
  std::string WTAappendix = (doWTA ? "_wta" : "");

  if(jet->getAttribute<double>("Tau1"+WTAappendix) == 0.){
    Warning("getTau21", "Tau1%s returns 0. -1 will be returned", WTAappendix.data());
  }

  return jet->getAttribute<double>("Tau2"+WTAappendix)/jet->getAttribute<double>("Tau1"+WTAappendix);
}


///=========================================
/// truth weakly decaying B-hadron pdg id's
///=========================================
bool GluonBBMiniNTuple::IsWeaklyDecayingBHadron(int mc_pdgId)
{
    
  int abs_mc_pdgId = fabs(mc_pdgId);
  if( abs_mc_pdgId==511  || abs_mc_pdgId==521  || abs_mc_pdgId==531  || abs_mc_pdgId==541  || 
      abs_mc_pdgId==5112 || abs_mc_pdgId==5122 || abs_mc_pdgId==5132 || abs_mc_pdgId==5142 ||
      abs_mc_pdgId==5212 || abs_mc_pdgId==5222 || abs_mc_pdgId==5232 || abs_mc_pdgId==5242 ||
      abs_mc_pdgId==5332)
    return true;
  
  return false;
}


///=========================================
/// deltaR match for n muons
///=========================================

std::vector<int> GluonBBMiniNTuple::JetMuondRMatch(std::vector<const xAOD::Jet*> inputjets, const xAOD::MuonContainer* muons, float radius)
{

  std::vector<int> nMatched(30, 0);

  for (auto muon : *muons){
    int matchedind = -1;
    int ijet = -1;
    int mindR = 1000;

    for (auto inputjet: inputjets){
      ijet ++;
      double dR = inputjet->p4().DeltaR(muon->p4());
      if (dR<mindR and dR<radius){
	matchedind = ijet;
      }
    }
    
    if (matchedind>=0){
      nMatched[matchedind] ++;
    }
  }
  return nMatched;
    
}

bool GluonBBMiniNTuple::GhostMatchTruthTrkJetToLargeRJet(  std::vector<TLorentzVector> TruthTrackJets4Vector_Preselected,
							   std::vector<int> TruthTrackJets4Vector_BMatched)
{

  double Pi = 3.14159265358979323846;
  const xAOD::JetContainer* TruthJets = 0;
  RETURN_CHECK("SetGluonBBEventCuts::GluonBBSelection()", HelperFunctions::retrieve(TruthJets, m_inTruthJetName, m_event, m_store, m_debug), "");

  const xAOD::Jet* leadingJet = *(TruthJets->begin());

  if (leadingJet->pt() <  m_FatJetPtTruthCut*1000.0){
    return false;
  }

  TLorentzVector Leading_TruthTrackJets4Vector;
  TLorentzVector SubLeading_TruthTrackJets4Vector;

  xAOD::JetConstituentVector consVec = leadingJet->getConstituents();  
  std::vector<fastjet::PseudoJet> ConstsAndGhosts;

  // get large R jet constituents
  for (auto cons : consVec){
    TLorentzVector cons_4m;
    cons_4m.SetPtEtaPhiE(cons->pt(), cons->eta(), cons->phi(), cons->e());
    fastjet::PseudoJet cons_cj( cons_4m.Px(), cons_4m.Py(), cons_4m.Pz(), cons_4m.E());
    ConstsAndGhosts.push_back( cons_cj);
  }

  // get track jets constituents
  for (auto trkjet : TruthTrackJets4Vector_Preselected){
    TLorentzVector p;
    p.SetPtEtaPhiM(1e-12, trkjet.Eta(), trkjet.Phi(), 0.0);
    fastjet::PseudoJet ghost(p);
    ConstsAndGhosts.push_back( ghost);
  }

  // cluster
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 1.0);
  fastjet::ClusterSequence cs(ConstsAndGhosts, jet_def);
  vector<fastjet::PseudoJet> csjets = fastjet::sorted_by_pt(cs.inclusive_jets());

  int nMatchedTrkJet = 0;
  int nMatchedTrkJet_B = 0;

  b_gluon_boosted_truth_pt = leadingJet->pt();
  b_gluon_boosted_truth_eta = leadingJet->eta();
  b_gluon_boosted_truth_phi = leadingJet->phi();
  b_gluon_boosted_truth_m = leadingJet->m();

  if (csjets.size()==0){
    return false;
  }
  else{
    for(auto csjet : csjets){
      double cs_phi = csjet.phi();
      if (cs_phi > Pi) cs_phi = cs_phi-2*Pi;

      double dR = sqrt( pow( (leadingJet->eta() - csjet.eta()), 2) + pow( (leadingJet->phi() - cs_phi  ), 2) );

      vector<fastjet::PseudoJet> matched_constituents = csjet.constituents();
      //std::cout<<"cs jet pt "<< csjet.pt() << " eta " << csjet.eta() << " phi "<<csjet.phi() <<std::endl;
      //std::cout<<"dR cs jet reco jet "<< dR<<std::endl;

      if (dR>0.1) continue;

      int itrkjet =-1;
      for(auto trkjet: TruthTrackJets4Vector_Preselected){
	itrkjet ++;
	for(auto matched_cons : matched_constituents){

	  if( fabs(trkjet.Eta()-matched_cons.eta())<0.0001){
	    if (nMatchedTrkJet==0){
	      b_jet_leading_ak2track_asso_truth_pt= trkjet.Pt()*1000.;
	      b_jet_leading_ak2track_asso_truth_eta= trkjet.Eta();
	      b_jet_leading_ak2track_asso_truth_phi= trkjet.Phi();
	      b_jet_leading_ak2track_asso_truth_m= trkjet.M();
	      b_jet_leading_ak2track_asso_truth_TruthB= TruthTrackJets4Vector_BMatched[itrkjet];
	      b_jet_leading_ak2track_asso_truth_TruthC= 1-TruthTrackJets4Vector_BMatched[itrkjet];
	      b_jet_leading_ak2track_asso_truth_TruthL= 1-TruthTrackJets4Vector_BMatched[itrkjet];
	      if (b_jet_leading_ak2track_asso_truth_TruthB ==1) nMatchedTrkJet_B++;
	    }
	    if (nMatchedTrkJet==1){
	      b_jet_subleading_ak2track_asso_truth_pt= trkjet.Pt()*1000.;
	      b_jet_subleading_ak2track_asso_truth_eta= trkjet.Eta();
	      b_jet_subleading_ak2track_asso_truth_phi= trkjet.Phi();
	      b_jet_subleading_ak2track_asso_truth_m= trkjet.M();
	      b_jet_subleading_ak2track_asso_truth_TruthB= TruthTrackJets4Vector_BMatched[itrkjet];
	      b_jet_subleading_ak2track_asso_truth_TruthC= 1-TruthTrackJets4Vector_BMatched[itrkjet];
	      b_jet_subleading_ak2track_asso_truth_TruthL= 1-TruthTrackJets4Vector_BMatched[itrkjet];
	      if (b_jet_subleading_ak2track_asso_truth_TruthB ==1) nMatchedTrkJet_B++;
	    }
	    nMatchedTrkJet ++;
	  }
	}
      }
    }
  }
  
  b_gluon_boosted_truth_n = nMatchedTrkJet;
  std::cout<<"nMatchedTrkJet "<<nMatchedTrkJet<<std::endl;
  if (nMatchedTrkJet < m_FatJetTruthNTrkJetMatched){
    return false;
  }

  return true;
}

