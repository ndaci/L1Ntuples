// -*- C++ -*-
//
// Package:    L1TriggerDPG/L1Ntuples
// Class:      L1GenTreeProducer
//
/**\class L1GenTreeProducer L1GenTreeProducer.cc L1TriggerDPG/L1Ntuples/src/L1GenTreeProducer.cc

 Description: Produces tree containing generated quantities
 Copied from L1RecoTreeProducer.cc + adapted to generated info

*/
//
// Original Author:  Nadir Daci
//         Created:
// $Id: L1GenTreeProducer.cc,v 1.0 2014/06/18 14:43:31 econte Exp $
//
//

// system include files
#include <memory>
#include <vector>

// framework
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"

// data formats
#include "DataFormats/JetReco/interface/GenJetCollection.h"
// #include "DataFormats/VertexReco/interface/Vertex.h"
// #include "DataFormats/VertexReco/interface/VertexFwd.h"
// #include "DataFormats/TrackReco/interface/Track.h"
// #include "DataFormats/TrackReco/interface/TrackFwd.h"
// #include "DataFormats/Common/interface/ValueMap.h"

// ROOT output stuff
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TTree.h"
#include "TF1.h"

//local  data formats
#include "L1TriggerDPG/L1Ntuples/interface/L1AnalysisGenJet.h"
#include "L1TriggerDPG/L1Ntuples/interface/L1AnalysisGenMet.h"
// #include "L1TriggerDPG/L1Ntuples/interface/L1AnalysisRecoVertex.h"
// #include "L1TriggerDPG/L1Ntuples/interface/L1AnalysisRecoTrack.h"

//
// class declaration
//

class L1GenTreeProducer : public edm::EDAnalyzer {
public:
  explicit L1GenTreeProducer(const edm::ParameterSet&);
  ~L1GenTreeProducer();


private:
  virtual void beginJob(void) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

public:
  L1Analysis::L1AnalysisGenJet*        jet;
  L1Analysis::L1AnalysisGenMet*        metCalo;
  L1Analysis::L1AnalysisGenMet*        metTrue;
  // L1Analysis::L1AnalysisRecoVertex*     vertices;
  // L1Analysis::L1AnalysisRecoTrack*      tracks;

  L1Analysis::L1AnalysisGenJetDataFormat*     jet_data;
  L1Analysis::L1AnalysisGenMetDataFormat*     metCalo_data;
  L1Analysis::L1AnalysisGenMetDataFormat*     metTrue_data;
  // L1Analysis::L1AnalysisRecoVertexDataFormat*           vertices_data;
  // L1Analysis::L1AnalysisRecoTrackDataFormat*            tracks_data;

private:

  // output file
  edm::Service<TFileService> fs_;

  // tree
  TTree * tree_;

  // EDM input tags
  edm::InputTag jetTag_;
  edm::InputTag metTagCalo_;
  edm::InputTag metTagTrue_;
  // edm::InputTag verticesTag_;
  // edm::InputTag tracksTag_;

  // debug stuff
  bool jetsMissing_;
  double jetptThreshold_;
  unsigned int maxJet_;
  std::vector<unsigned int> excludePdgId_;
  // unsigned int maxVtx_;
  // unsigned int maxTrk_;
};



L1GenTreeProducer::L1GenTreeProducer(const edm::ParameterSet& iConfig):
  jetTag_(    iConfig.getUntrackedParameter("jetTag",edm::InputTag("ak5GenJets"))),
  metTagCalo_(iConfig.getUntrackedParameter("metTagCalo",edm::InputTag("genMetCalo"))),
  metTagTrue_(iConfig.getUntrackedParameter("metTagTrue",edm::InputTag("genMetTrue"))),
  // verticesTag_(iConfig.getUntrackedParameter("verticesTag",edm::InputTag("offlinePrimaryVertices"))),
  // tracksTag_(iConfig.getUntrackedParameter("tracksTag",edm::InputTag("generalTracks"))),
  jetsMissing_(false)
{

  jetptThreshold_ = iConfig.getParameter<double>      ("jetptThreshold");
  maxJet_         = iConfig.getParameter<unsigned int>("maxJet");
  excludePdgId_   = iConfig.getParameter< std::vector<unsigned int> > ("excludePdgId");
  // maxVtx_         = iConfig.getParameter<unsigned int>("maxVtx");
  // maxTrk_         = iConfig.getParameter<unsigned int>("maxTrk");

  jet           = new L1Analysis::L1AnalysisGenJet();
  metCalo       = new L1Analysis::L1AnalysisGenMet();
  metTrue       = new L1Analysis::L1AnalysisGenMet();
  // vertices      = new L1Analysis::L1AnalysisRecoVertex();
  // tracks        = new L1Analysis::L1AnalysisRecoTrack();

  jet_data           = jet->getData();
  metCalo_data       = metCalo->getData();
  metTrue_data       = metTrue->getData();
  // vertices_data      = vertices->getData();
  // tracks_data        = tracks->getData();

  // set up output
  tree_=fs_->make<TTree>("GenTree", "GenTree");
  tree_->Branch("Jet",           "L1Analysis::L1AnalysisGenJetDataFormat",         &jet_data,                32000, 3);
  tree_->Branch("MetTrue",       "L1Analysis::L1AnalysisGenMetDataFormat",         &metTrue_data,            32000, 3);
  tree_->Branch("MetCalo",       "L1Analysis::L1AnalysisGenMetDataFormat",         &metCalo_data,            32000, 3);
  // tree_->Branch("Vertices",      "L1Analysis::L1AnalysisRecoVertexDataFormat",      &vertices_data,           32000, 3);
  // tree_->Branch("Tracks",        "L1Analysis::L1AnalysisRecoTrackDataFormat",       &tracks_data,             32000, 3);

}


L1GenTreeProducer::~L1GenTreeProducer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void L1GenTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Reset L1Analysis containers
  jet->Reset();
  metCalo->Reset();
  metTrue->Reset();
  // vertices->Reset();
  // tracks->Reset();

  // Handles
  edm::Handle<reco::GenJetCollection> genJets;
  edm::Handle<reco::GenMETCollection> genMetCalo; // missing transverse energy
  edm::Handle<reco::GenMETCollection> genMetTrue; // missing transverse energy
  // edm::Handle<reco::VertexCollection> recoVertices;
  // edm::Handle<reco::TrackCollection>  recoTrkRef;

  // Get collections using InputTags
  iEvent.getByLabel(jetTag_, genJets);
  iEvent.getByLabel(metTagCalo_, genMetCalo);
  iEvent.getByLabel(metTagTrue_, genMetTrue);
  // iEvent.getByLabel(verticesTag_, recoVertices);
  // iEvent.getByLabel(tracksTag_, recoTrkRef);

  // Store info in L1Analysis containers
  if (genJets.isValid()) {
    jet->SetGenJet(iEvent, iSetup, genJets, maxJet_);
    metCalo->SetHtMht(genJets, jetptThreshold_, excludePdgId_);
    metTrue->SetHtMht(genJets, jetptThreshold_, excludePdgId_);
  }
  else {
    if (!jetsMissing_) {edm::LogWarning("MissingProduct") << "GenJets not found.  Branch will not be filled" << std::endl;}
    jetsMissing_ = true;
  }
  
  if (genMetCalo.isValid()) {
    metCalo->SetMet(genMetCalo);
  }
  else {
    edm::LogWarning("MissingProduct") << "metCalo not found.  Branch will not be filled"<<std::endl;
  }

  if (genMetTrue.isValid()) {
    metTrue->SetMet(genMetTrue);
  }
  else {
    edm::LogWarning("MissingProduct") << "metTrue not found.  Branch will not be filled"<<std::endl;
  }
  
  // if (recoTrkRef.isValid()) {
  //   const reco::TrackCollection recoTracks(*recoTrkRef.product());
  //   tracks->SetTracks(recoTracks, maxTrk_);
  //  }
  //  else {
  //    edm::LogWarning("MissingProduct") << "tracks not found.  Branch will not be filled"<<std::endl;
  //   }

   // if (recoVertices.isValid()) {
   //  vertices->SetVertices(recoVertices, maxVtx_);
   // }
   // else {
   //   edm::LogWarning("MissingProduct") << "vertices not found.  Branch will not be filled"<<std::endl;
   //  }

  tree_->Fill();

}

// ------------ method called once each job just before starting event loop  ------------
void
L1GenTreeProducer::beginJob(void)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
L1GenTreeProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1GenTreeProducer);
