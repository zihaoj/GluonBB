#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <EventLoop/OutputStream.h>

#include <xAODJet/JetContainer.h>
#include <xAODTracking/VertexContainer.h>
#include <xAODTracking/TrackParticlexAODHelpers.h>
#include <xAODEventInfo/EventInfo.h>
#include <AthContainers/ConstDataVector.h>
#include <SampleHandler/MetaFields.h>

#include "xAODParticleEvent/ParticleContainer.h"

#include <xAODAnaHelpers/HelpTreeBase.h>
#include <GluonBB/GluonBBTemplateHisto.h>


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
ClassImp(GluonBBTemplateHisto)

using std::cout;  using std::endl;
using std::string; using std::vector;

GluonBBTemplateHisto :: GluonBBTemplateHisto () :
  m_name(""),
  m_isMC(false),
  m_inTrackJetName(""),
  m_boostedSysName(""),
  m_muonContainerName(""),
  m_InDetTrackSmearingTool(nullptr),
  m_debug(false)
{
  this->SetName("GluonBBTemplateHisto"); // needed if you want to retrieve this algo with wk()->getAlg(ALG_NAME) downstream
}

EL::StatusCode GluonBBTemplateHisto :: setupJob (EL::Job& job)
{
  Info("setupJob()", "Calling setupJob \n");
  job.useXAOD();
  xAOD::Init("GluonBBTemplateHisto").ignore();

  EL::OutputStream outForTree(m_name.c_str());
  job.outputAdd (outForTree);
  Info("setupJob()", "Ready\n");

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode GluonBBTemplateHisto :: initialize ()
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

  std::string indet_smearing_tool_name = std::string("indet_smearing_tool_template");
  m_InDetTrackSmearingTool = new InDet::InDetTrackSmearingTool( indet_smearing_tool_name);
  RETURN_CHECK( "InDetTrackSmearingTool::initialize()", m_InDetTrackSmearingTool->initialize(), "Failed to properly initialize the InDetTrackSmearingTool");
  CP::SystematicSet systSet = {InDet::TrackSystematicMap[InDet::TRK_RES_D0_MEAS]};
			       //InDet::TrackSystematicMap[InDet::TRK_RES_Z0_MEAS]};
  auto systCode = m_InDetTrackSmearingTool->applySystematicVariation(systSet);


  if(m_debug) Info("initialize()", "left");
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode GluonBBTemplateHisto::AddTree(string syst = "")
{
  Info("AddTree()", "%s", m_name.c_str() );
  // needed here and not in initalize since this is called first

  string treeName = "GluonBBTemplateHisto";
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


  // Event selection flags
  outTree->Branch("jet_ak2track_asso_n", &b_jet_ak2track_asso_n);
  outTree->Branch("jet_ak2track_asso_pt", &b_jet_ak2track_asso_pt);
  outTree->Branch("jet_ak2track_asso_eta", &b_jet_ak2track_asso_eta);
  outTree->Branch("jet_ak2track_asso_phi", &b_jet_ak2track_asso_phi);
  outTree->Branch("jet_ak2track_asso_m", &b_jet_ak2track_asso_m);
  outTree->Branch("jet_ak2track_asso_Sd0", &b_jet_ak2track_asso_Sd0);
  outTree->Branch("jet_ak2track_asso_Sd0_Sub", &b_jet_ak2track_asso_Sd0_Sub);
  outTree->Branch("jet_ak2track_asso_MV2c00", &b_jet_ak2track_asso_MV2c00);
  outTree->Branch("jet_ak2track_asso_MV2c10", &b_jet_ak2track_asso_MV2c10);
  outTree->Branch("jet_ak2track_asso_MV2c20", &b_jet_ak2track_asso_MV2c20);
  outTree->Branch("jet_ak2track_asso_IP3D", &b_jet_ak2track_asso_IP3D);
  outTree->Branch("jet_ak2track_asso_IP2D", &b_jet_ak2track_asso_IP2D);
  outTree->Branch("jet_ak2track_asso_sys", &b_jet_ak2track_asso_sys);
  outTree->Branch("jet_ak2track_asso_sysname", &b_jet_ak2track_asso_sysname);
  outTree->Branch("jet_ak2track_asso_TruthFlav", &b_jet_ak2track_asso_TruthFlav);
  outTree->Branch("jet_ak2track_asso_TruthB", &b_jet_ak2track_asso_TruthB);
  outTree->Branch("jet_ak2track_asso_TruthC", &b_jet_ak2track_asso_TruthC);
  outTree->Branch("jet_ak2track_asso_TruthL", &b_jet_ak2track_asso_TruthL);
  outTree->Branch("jet_ak2track_asso_nMuon", &b_jet_ak2track_asso_nMuon);

  outTree->Branch("weight",    &m_weight,    "weight/F");
  outTree->Branch("weight_xs", &m_weight_xs, "weight_xs/F");
    

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode GluonBBTemplateHisto :: histInitialize () { return EL::StatusCode::SUCCESS; }
EL::StatusCode GluonBBTemplateHisto :: fileExecute    () { return EL::StatusCode::SUCCESS; }
EL::StatusCode GluonBBTemplateHisto :: changeInput    (bool /*firstFile*/) { return EL::StatusCode::SUCCESS; }

EL::StatusCode GluonBBTemplateHisto::execute ()
{

  const xAOD::EventInfo* eventInfo(0);
  RETURN_CHECK("GluonBBTemplateHisto::execute()", HelperFunctions::retrieve(eventInfo, "EventInfo", m_event, m_store), "");
  //
  // JZXW Pile-up Fluctuation Killer
  //
  //if(m_isMC && (killPileupFluctuation() != 1)) return EL::StatusCode::SUCCESS;

  // 
  // For the cut flow
  //
  m_iCutflow = m_cutflowFirst;

  //
  //
  if(m_isMC){
    m_mcEventWeight = eventInfo->mcEventWeight();
  }  else{
    m_mcEventWeight = 1;
  }

  if (m_boostedSysName.empty()) {

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

EL::StatusCode GluonBBTemplateHisto::executeSingle(string resolvedSys, string boostedSys, bool countEvents) {
  if(m_debug) cout << " In executeSingle"  << resolvedSys << " " << boostedSys << endl;


  // Start with stuff which is common to boosted and resolved
  //
  string syst = "";
  if(!boostedSys.empty())  syst = "Boosted_" +boostedSys;

  if( m_helpTree.find( syst ) == m_helpTree.end() ) { AddTree( syst ); }

  const xAOD::EventInfo* eventInfo(0);
  RETURN_CHECK("GluonBBTemplateHisto::execute()", HelperFunctions::retrieve(eventInfo, "EventInfo", m_event, m_store), "");
  
  const xAOD::VertexContainer* vertices(0);
  RETURN_CHECK("GluonBBTemplateHisto::execute()", HelperFunctions::retrieve(vertices, "PrimaryVertices", m_event, m_store), "");

  const xAOD::Vertex *pv = 0;
  pv = vertices->at( HelperFunctions::getPrimaryVertexLocation( vertices ) );

  const xAOD::MuonContainer* muons(nullptr);
  RETURN_CHECK("GluonBBMiniNTuple::execute()", HelperFunctions::retrieve(muons, m_muonContainerName, m_event, m_store), (m_muonContainerName).c_str());


  // Retrieve the container of track jets jets
  const xAOD::JetContainer* TrackJets(0);
  if(m_debug) cout << " Getting track Jets: "  << m_inTrackJetName<<std::endl;
  RETURN_CHECK("GluonBBTemplateHisto::execute()", HelperFunctions::retrieve(TrackJets, m_inTrackJetName, m_event, m_store), "");  

  if(m_debug){ 
    cout << "Run/Event " << eventInfo->runNumber() << " / " << eventInfo->eventNumber() << endl;
  }

  b_jet_ak2track_asso_n = 0;
  b_jet_ak2track_asso_pt.clear();
  b_jet_ak2track_asso_eta.clear();
  b_jet_ak2track_asso_phi.clear();
  b_jet_ak2track_asso_m.clear();
  b_jet_ak2track_asso_Sd0.clear();
  b_jet_ak2track_asso_Sd0_Sub.clear();
  b_jet_ak2track_asso_MV2c00.clear();
  b_jet_ak2track_asso_MV2c10.clear();
  b_jet_ak2track_asso_MV2c20.clear();
  b_jet_ak2track_asso_IP3D.clear();
  b_jet_ak2track_asso_IP2D.clear();
  b_jet_ak2track_asso_sys.clear();
  b_jet_ak2track_asso_sysname.clear();
  b_jet_ak2track_asso_TruthFlav.clear();
  b_jet_ak2track_asso_TruthB.clear();
  b_jet_ak2track_asso_TruthC.clear();
  b_jet_ak2track_asso_TruthL.clear();
  b_jet_ak2track_asso_nMuon.clear();

  b_jet_ak2track_asso_n = TrackJets->size();

  
  std::vector<int> JetMuonMatchVec = GluonBBTemplateHisto::JetMuondRMatch(TrackJets, muons, 0.2);

  for(unsigned int i=0; i<TrackJets->size(); i++){
    const xAOD::Jet* trkjet = TrackJets->at(i);

    const xAOD::BTagging* BTag = trkjet->btagging();

    float maxSd0=0;
    float secmaxSd0=-1;
    float LeadSd0=-99;
    float SubSd0=-99;

    xAOD::JetConstituentVector consVec = trkjet->getConstituents();

    try{
      for (auto cons : consVec){
	const xAOD::TrackParticle* track = dynamic_cast <const xAOD::TrackParticle*> (cons->rawConstituent());

	xAOD::TrackParticle track_smear = *track;
	m_InDetTrackSmearingTool->applyCorrection( track_smear);

	double d0sig = xAOD::TrackingHelpers::d0significance( track, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), eventInfo->beamPosSigmaXY());
	double sign_rphi = sin( trkjet->p4().Phi() - track->p4().Phi() ) * track->d0();
	int sign = (sign_rphi>0 ? 1 : -1);
	d0sig = sign* fabs(d0sig);

	double d0sig_corr = xAOD::TrackingHelpers::d0significance( &track_smear, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), eventInfo->beamPosSigmaXY());
	double sign_rphi_smear = sin( trkjet->p4().Phi() - track_smear.p4().Phi() ) * track_smear.d0();
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
      }
    }
    catch(...){
      Warning("executeSingle()", "element link problem for track converstion from iparticle");
    }


    if (m_isMC){
      std::vector<const xAOD::TruthParticle*> BHadrons = trkjet->getAssociatedObjects<xAOD::TruthParticle>("ConeExclBHadronsFinal");
      std::vector<const xAOD::TruthParticle*> CHadrons = trkjet->getAssociatedObjects<xAOD::TruthParticle>("ConeExclCHadronsFinal");
      if ( BHadrons.size()>0){
	b_jet_ak2track_asso_TruthFlav.push_back(5);
	b_jet_ak2track_asso_TruthB.push_back(1);
	b_jet_ak2track_asso_TruthC.push_back(0);
	b_jet_ak2track_asso_TruthL.push_back(0);
      }
      else if ( CHadrons.size()>0){
	b_jet_ak2track_asso_TruthFlav.push_back(4);
	b_jet_ak2track_asso_TruthB.push_back(0);
	b_jet_ak2track_asso_TruthC.push_back(1);
	b_jet_ak2track_asso_TruthL.push_back(0);
      }
      else{
	b_jet_ak2track_asso_TruthFlav.push_back(1);
	b_jet_ak2track_asso_TruthB.push_back(0);
	b_jet_ak2track_asso_TruthC.push_back(0);
	b_jet_ak2track_asso_TruthL.push_back(1);
      }
      
      b_jet_ak2track_asso_pt.push_back(   trkjet->pt());
      b_jet_ak2track_asso_eta.push_back(  trkjet->eta());
      b_jet_ak2track_asso_phi.push_back(  trkjet->phi());
      b_jet_ak2track_asso_m.push_back(    trkjet->m());
      b_jet_ak2track_asso_Sd0.push_back(LeadSd0);
      b_jet_ak2track_asso_Sd0_Sub.push_back(SubSd0);
      b_jet_ak2track_asso_MV2c00.push_back( MV2(trkjet, "MV2c00") );
      b_jet_ak2track_asso_MV2c10.push_back( MV2(trkjet, "MV2c10") );
      b_jet_ak2track_asso_MV2c20.push_back( MV2(trkjet, "MV2c20") );
      b_jet_ak2track_asso_IP3D.push_back( BTAGIP(trkjet, "IP3D") );
      b_jet_ak2track_asso_IP2D.push_back( BTAGIP(trkjet, "IP2D") );
      b_jet_ak2track_asso_nMuon.push_back( JetMuonMatchVec[i]);
      b_jet_ak2track_asso_sys.push_back(trkjet->auxdecor<std::vector<float> >("BTag_SF_FixedCutBEff_77"));

    }
  }



  // fill the tree
  if(m_debug) cout << " Filling Tree"  << endl;
  m_helpTree[syst]->Fill();

  if(m_debug) cout << " Filled Tree"  << endl;

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode GluonBBTemplateHisto :: postExecute () { return EL::StatusCode::SUCCESS; }

EL::StatusCode GluonBBTemplateHisto :: finalize () {

  Info("finalize()", m_name.c_str());

  if (!m_helpTree.empty()){
    for( auto tree : m_helpTree) {
      if (tree.second) delete tree.second;
    }
  }

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode GluonBBTemplateHisto :: histFinalize ()
{
  TFile * treeFile = wk()->getOutputFile( m_name.c_str() );
  
  return EL::StatusCode::SUCCESS;
}



std::vector<int> GluonBBTemplateHisto::JetMuondRMatch(const xAOD::JetContainer* inputjets, const xAOD::MuonContainer* muons, float radius)
{

  std::vector<int> nMatched(50, 0);

  for (auto muon : *muons){
    int matchedind = -1;
    int ijet = -1;
    int mindR = 1000;

    for (auto inputjet: *inputjets){
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
