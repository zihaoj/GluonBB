#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <GluonBB/SetGluonBBEventCuts.h>

// EDM includes:
#include "AthContainers/ConstDataVector.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODCore/ShallowCopy.h"
#include "xAODRootAccess/TStore.h"
#include <xAODTruth/TruthParticleContainer.h>
#include <xAODTruth/TruthParticle.h>
#include <vector>
#include <JetReclustering/JetReclusteringTool.h>
//using std::vector;
#include "xAODAnaHelpers/HelperFunctions.h"
#include <xAODAnaHelpers/tools/ReturnCheck.h>

// Jet xAOD EDM container
#include "xAODJet/JetContainer.h"
#include "xAODBTagging/BTagging.h"

#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>

// jet trimming
#include <fastjet/tools/Filter.hh>


// Particle 
#include "xAODParticleEvent/Particle.h"
#include "xAODParticleEvent/ParticleContainer.h"
#include "xAODParticleEvent/ParticleAuxContainer.h"

#include <TSystem.h> // used to define JERTool calibration path
#include "TEnv.h"

#include "GluonBB/errorcheck.h"

// this is needed to distribute the algorithm to the workers
ClassImp(SetGluonBBEventCuts)


using namespace std;


SetGluonBBEventCuts :: SetGluonBBEventCuts ():
  m_name(""),
  m_debug(false),
  m_dothirdJet(false),
  m_dothirdTrkJet(false),
  m_isMC(false),
  m_inJetName(""),
  m_inTruthJetName(""),
  m_trackJet("GhostAntiKt2TrackJet"),
  m_inputAlgo(""),
  m_outJetName(""),
  m_outLeadingAkt4Jet(""),
  m_outputAlgo(""),
  m_leadingJetPtCut(350e3),
  m_subleadingJetPtCut(250e3),
  m_trackJetPtCut(10e3),
  m_trackJetEtaCut(2.5),
  m_inTruthParticleName(""),
  m_trackJetBtagCut(0.3706)
{
}


EL::StatusCode SetGluonBBEventCuts :: setupJob (EL::Job& job)
{
  job.useXAOD();
  xAOD::Init("SetGluonBBEventCuts").ignore(); // call before opening first file
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode SetGluonBBEventCuts :: histInitialize ()
{
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode SetGluonBBEventCuts :: fileExecute ()
{
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode SetGluonBBEventCuts :: changeInput (bool )
{
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode SetGluonBBEventCuts :: configure () 
{
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode SetGluonBBEventCuts :: initialize ()
{

  if ( this->configure() == EL::StatusCode::FAILURE ) 
  {
    Error("initialize()", "Failed to properly configure. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  

  //
  // Set isMC flag
  //
  m_event = wk()->xaodEvent();
  m_store = wk()->xaodStore();

  const xAOD::EventInfo* eventInfo(nullptr);
  RETURN_CHECK("ResonanceAlgorithm::initialize()", HelperFunctions::retrieve(eventInfo, "EventInfo", m_event, m_store, m_debug), "");
  m_isMC = ( eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ) ? true : false;

//
//  const std::string JetRecToolName ="JetRecTool";
//  if ( asg::ToolStore::contains<IJetExecuteTool>(JetRecToolName.c_str()) ) {
//  m_jetReclusteringTool(


  return EL::StatusCode::SUCCESS;
}

EL::StatusCode SetGluonBBEventCuts :: execute ()
{
  if(m_debug) cout << "In Select Jets" << endl;
  //
  // No systematics
  //
  if ( m_inputAlgo.empty() ) {

    GluonBBSelection("");
  //
  // Do Systemaitcs
  //
  }else{

    // get vector of string giving the names
    vector<string>* systNames(nullptr);
    RETURN_CHECK("GluonBBSelection::execute()", HelperFunctions::retrieve(systNames, m_inputAlgo, 0, m_store, m_debug) ,"");

    // loop over systematics
    vector< string >* vecOutContainerNames = new vector< string >;
    for ( string& systName : *systNames ) {
      
      if(m_debug) Info("execute",  "systName %s", systName.c_str());
      
      GluonBBSelection(systName);
      
      vecOutContainerNames->push_back( systName );
    }

    // save list of systs that should be considered down stream
    RETURN_CHECK( "MakeDiJets::execute()", m_store->record( vecOutContainerNames, m_outputAlgo), "Failed to record vector of output container names.");
  }


  return EL::StatusCode::SUCCESS;
}


EL::StatusCode SetGluonBBEventCuts :: GluonBBSelection (std::string systName)
{
  if(m_debug) cout << "In SetGluonBBEventCuts (Boosted) " << endl;

  const xAOD::EventInfo* eventInfo = 0;

  if ( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ) {
    this->Error("execute()", "Failed to retrieve event info collection. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  // Record a container for final objects (Gluon candidates)
  xAOD::ParticleContainer*     GluonJetsAll    = new xAOD::ParticleContainer();
  xAOD::ParticleAuxContainer*  GluonJetsAllAux = new xAOD::ParticleAuxContainer();
  GluonJetsAll->setStore( GluonJetsAllAux ); //< Connect the two 
  if( ! m_store->record( GluonJetsAll,    m_outJetName+systName       ) ) { return EL::StatusCode::FAILURE; }
  if( ! m_store->record( GluonJetsAllAux, m_outJetName+systName+"Aux." ) ) { return EL::StatusCode::FAILURE; } 

  // Initialize the cutflow
  std::vector<std::string> cuts = {
    "GluonBBPass_LeadJetPt",
    "GluonBBPass_LeadJet_AssoLeadTrackJet"
    "GluonBBPass_LeadJet_AssoSubLeadTrackJet",
    "GluonBBPass_0Btag",
    "GluonBBPass_1Btag",
    "GluonBBPass_2Btag",
  };

  eventInfo->auxdecor< std::vector<std::string> >("BoostedCuts") = cuts;  
  for(auto c : cuts) eventInfo->auxdecor< bool >(c) = false;

  // N jets
  //
  const xAOD::JetContainer* jets = 0;
  RETURN_CHECK("SetGluonBBEventCuts::GluonBBSelection()", HelperFunctions::retrieve(jets, m_inJetName+systName, m_event, m_store, m_debug), "");

  unsigned int nJets = jets->size();

  if(nJets < 1) return EL::StatusCode::SUCCESS;
  bool subleadingJet_exist = nJets > 1;//use this to pick out single jet events

  //
  // leading jet pt cut                                                                                         
  //
  const xAOD::Jet* leadingJet = *(jets->begin());
  bool passLeadJetPt = (leadingJet->pt() > m_leadingJetPtCut);
  eventInfo->auxdecor< bool > ("GluonBBPass_LeadJetPt") = passLeadJetPt;

  // check associated track-jets
  //
  // get un-trimmed jet area -- leading calo-jet
  const xAOD::Jet* jet_parent_leadingJet = 0;

  try{
    auto el = leadingJet->auxdata<ElementLink<xAOD::JetContainer> >("Parent");
    if(!el.isValid()){
      Warning("execute()", "Invalid link to \"Parent\" from leading calo-jet");
    }
    else{
      jet_parent_leadingJet = (*el);
    }
  }
  catch (...){
    Warning("execute()", "Unable to fetch \"Parent\" link from leading calo-jet");
  }

  if(jet_parent_leadingJet == 0){
    Warning("execute()", "Trimmed jet area will be used for leading calo-jet");
    jet_parent_leadingJet = leadingJet;
  }

  // initiate to appropriate pointer
  const xAOD::Jet* leadingTrackJet_leadingJet       = 0;
  const xAOD::Jet* subleadingTrackJet_leadingJet    = 0;
  std::vector<const xAOD::Jet*> sel_assotrkjets_leadingJet;
  std::vector<const xAOD::Jet*> sel_assotrkjets_subleadingJet;


  // get associated track-jets; leading jet
  std::vector<const xAOD::Jet*> assotrkjets_leadingJet;
  try{
    assotrkjets_leadingJet = jet_parent_leadingJet->getAssociatedObjects<xAOD::Jet>(m_trackJet);
  }
  catch (...){
    Warning("execute()", "Unable to fetch \"%s\" link from leading calo-jet", m_trackJet.data());
  }
  std::sort(assotrkjets_leadingJet.begin(), assotrkjets_leadingJet.end(), SetGluonBBEventCuts::SortPt);

  // apply selection on associated track-jets
  for(auto TrackJet : assotrkjets_leadingJet){
    if(!SelectTrackJet(TrackJet)) continue;
    sel_assotrkjets_leadingJet.push_back(TrackJet);
  }
  if(sel_assotrkjets_leadingJet.size() >= 1){
    leadingTrackJet_leadingJet = sel_assotrkjets_leadingJet[0];
    if(sel_assotrkjets_leadingJet.size() >= 2){
      subleadingTrackJet_leadingJet = sel_assotrkjets_leadingJet[1];
    }
  }

  eventInfo->auxdecor< bool > ("GluonBBPass_LeadJet_AssoLeadTrackJet") = leadingTrackJet_leadingJet;
  eventInfo->auxdecor< bool > ("GluonBBPass_LeadJet_AssoSubLeadTrackJet") = subleadingTrackJet_leadingJet;

  // leading track-jet associated to leading calo-jet
  double leadingTrackJet_leadingJet_BTag       = -1000;
  if( leadingTrackJet_leadingJet ) 
    leadingTrackJet_leadingJet->btagging()->MVx_discriminant("MV2c10", leadingTrackJet_leadingJet_BTag);
  // subleading track-jet associated to leading calo-jet
  double subleadingTrackJet_leadingJet_BTag    = -1000;
  if( subleadingTrackJet_leadingJet ) 
    subleadingTrackJet_leadingJet->btagging()->MVx_discriminant("MV2c10", subleadingTrackJet_leadingJet_BTag);

  // if (mvx > 0.473) // MV2c10 @ 60%
  // if (mvx > -0.046)  // MV2c10 @ 70%
  // if (mvx > -0.819) // MV2c10 @ 85%


  ////////////////////////////////////////////////
  // Now we define signal and control regions
  // Regions 0,1,2,3,4 b-tagged events
  /////////////////////////////////////////////////     

  int numOfBTags=0;
  if(leadingTrackJet_leadingJet_BTag       > m_trackJetBtagCut) numOfBTags++;
  if(subleadingTrackJet_leadingJet_BTag    > m_trackJetBtagCut) numOfBTags++;

  eventInfo->auxdecor< bool >("GluonBBPass_0Btag") = bool(numOfBTags==0);   
  eventInfo->auxdecor< bool >("GluonBBPass_1Btag") = bool(numOfBTags==1);   
  eventInfo->auxdecor< bool >("GluonBBPass_2Btag") = bool(numOfBTags==2);   

  // Finally create the container and its auxiliary store to hold the jets
  //
  xAOD::Particle* leadingJetParticle = new xAOD::Particle();

  // store leading gluon jet
  leadingJetParticle->makePrivateStore();
  leadingJetParticle->setPxPyPzE(leadingJet->p4().Px(),leadingJet->p4().Py(),leadingJet->p4().Pz(),leadingJet->p4().E());
  leadingJetParticle->auxdecor< const xAOD::Jet* >("caloJet") = leadingJet;
  leadingJetParticle->auxdecor< const xAOD::Jet* >("leadJet") = leadingTrackJet_leadingJet;          
  leadingJetParticle->auxdecor< const xAOD::Jet* >("sublJet") = subleadingTrackJet_leadingJet;       
  leadingJetParticle->auxdecor< std::vector<const xAOD::Jet*> >("allTrkJet") = sel_assotrkjets_leadingJet;
  float dRjj=-1;
  if(leadingTrackJet_leadingJet && subleadingTrackJet_leadingJet)
    dRjj = leadingTrackJet_leadingJet->p4().DeltaR(subleadingTrackJet_leadingJet->p4());
  leadingJetParticle->auxdecor< float > ("dRjj") = dRjj;
  leadingJetParticle->auxdecor< int > ("ak2track_asso_n") = sel_assotrkjets_leadingJet.size();
  GluonJetsAll->push_back(leadingJetParticle);


  //// get charged truth particles

  //std::cout<<" looping truth particles "<<m_isMC<<std::endl;

//  std::vector<xAOD::TruthParticle> ChargedTruthParticle;
//
//  const xAOD::TruthParticleContainer* truth_particles_unsort(0);
//  xAOD::TruthParticleContainer* truth_particles(0);
//  if (!m_inTruthParticleName.empty() && m_isMC){
//    if(m_debug) cout << " Getting Truth Particles"  << endl;  
//    RETURN_CHECK("GluonBBMiniNTuple::execute()", HelperFunctions::retrieve(truth_particles_unsort, m_inTruthParticleName, m_event, m_store), m_inTruthParticleName.c_str());
//    truth_particles = new xAOD::TruthParticleContainer(*truth_particles_unsort);//memory leak
//
//    //truth_particles->sort(sort_pt<xAOD::TruthParticle>());
//
//    for(unsigned int i=0; i<truth_particles->size(); i++){
//      const xAOD::TruthParticle* curr_truth = truth_particles->at(i);
//      if(!curr_truth) continue;
//      int status = curr_truth->status();
//      double charge = curr_truth->charge();
//      double pt = curr_truth->pt();
//      double eta = curr_truth->eta();
//      int barcode = curr_truth->barcode();
//      
//      // apply selection according to https://twiki.cern.ch/twiki/bin/view/AtlasProtected/TrackingCPICHEP2016#Truth_definitions
//      // we apply phase space and primary particles
//      if (status != 1 ) continue;
//      if (charge == 0 ) continue;
//      if (barcode <0 or barcode>2e5) continue;
//      if (pt <500 and fabs(eta)>2.5) continue;
//
//      ChargedTruthParticle.push_back( *curr_truth);
//    }
//  }
//
//
//  // Record a container for truth track jets
//  xAOD::ParticleContainer*     TruthTrackJets    = new xAOD::ParticleContainer();
//  xAOD::ParticleAuxContainer*  TruthTrackJetsAux = new xAOD::ParticleAuxContainer();
//  TruthTrackJets->setStore( TruthTrackJetsAux ); //< Connect the two 
//  std::vector<TLorentzVector> TruthTrackJets4Vector = HelperFunctions::TruthTrackJetReclustering(ChargedTruthParticle, 0.2, 0, fastjet::antikt_algorithm);
//
//  // loop over track jets to dR match truth b-hadron
//  for(auto trkjet : TruthTrackJets4Vector){
//    xAOD::Particle* trkJetParticle = new xAOD::Particle();    
//    trkJetParticle->makePrivateStore();
//    trkJetParticle->setPxPyPzE( trkjet.Px()*1000. ,trkjet.Py()*1000., trkjet.Pz()*1000.,trkjet.E()*1000.);
//    
//    std::cout<<"truth track jet eta "<<trkjet.Eta()<<" phi "<<trkjet.Phi()<<" pt "<<trkjet.Pt()<<std::endl;
//
//    for (auto truthpart : *truth_particles_unsort){
//      if (truthpart->pt()<5000) continue;
//
//      double dR = sqrt( pow(truthpart->eta()-trkjet.Eta(), 2) +pow(truthpart->phi()-trkjet.Phi(), 2));
//      //std::cout<<"truth particle eta "<<truthpart->eta()<<" truth particle phi "<<truthpart->phi()<<" truth particle pt "<<truthpart->pt()<<std::endl;
//      
//      //if (IsWeaklyDecayingBHadron(truthpart->pdgId() ) and dR<0.2){
//
//      TruthTrackJets->push_back(trkJetParticle);
//      //std::cout<<"found truth b-tagged track jets"<<std::endl;
//      //}
//    }
//  }
//  
////
////  // Record a container for truth large R jets
////  xAOD::ParticleContainer*     TruthLargeRJets    = new xAOD::ParticleContainer();
////  xAOD::ParticleAuxContainer*  TruthLargeRJetsAux = new xAOD::ParticleAuxContainer();
////  TruthLargeRJets->setStore( TruthLargeRJetsAux);
////  if( ! m_store->record( TruthLargeRJets,    "TruthBMatchedLargeRJets"+systName       ) ) { return EL::StatusCode::FAILURE; }
////  if( ! m_store->record( TruthLargeRJetsAux, "TruthBMatchedLargeRJets"+systName+"Aux." ) ) { return EL::StatusCode::FAILURE; } 
//
//  
//  //GhostMatchTruthTrkJetToLargeRJet( TruthLargeRJets, TruthTrackJets);
//
//  if( ! m_store->record( TruthTrackJets,    "TruthBMatchedTrackJets"+systName       ) ) { return EL::StatusCode::FAILURE; }
//  if( ! m_store->record( TruthTrackJetsAux, "TruthBMatchedTrackJets"+systName+"Aux." ) ) { return EL::StatusCode::FAILURE; } 

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode SetGluonBBEventCuts :: postExecute ()
{
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode SetGluonBBEventCuts :: finalize ()
{
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode SetGluonBBEventCuts :: histFinalize ()
{
  return EL::StatusCode::SUCCESS;
}


bool SetGluonBBEventCuts::SelectTrackJet(const xAOD::Jet* TrackJet)
{
  if( TrackJet->pt() < m_trackJetPtCut )            return false;
  if( fabs(TrackJet->eta()) > m_trackJetEtaCut )    return false;
  if( TrackJet->numConstituents() < 2 )             return false;

  return true;
}


