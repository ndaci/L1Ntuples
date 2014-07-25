#ifndef __L1Analysis_L1AnalysisGenJet_H__
#define __L1Analysis_L1AnalysisGenJet_H__

//-------------------------------------------------------------------------------
// Created 18/06/2014 - Nadir Daci
// 
//
// Copied from L1AnalysisRecoJet.h
//-------------------------------------------------------------------------------

#include "DataFormats/JetReco/interface/GenJetCollection.h"
/// #include "RecoJets/JetProducers/interface/JetMatchingTools.h"

#include "L1AnalysisGenJetDataFormat.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

//class JetCorrector;

namespace L1Analysis
{
  class L1AnalysisGenJet
  {
  public:
    L1AnalysisGenJet();
    ~L1AnalysisGenJet();
    
    //void Print(std::ostream &os = std::cout) const;
    void SetGenJet(const edm::Event& event,
		   const edm::EventSetup& setup,
		   const edm::Handle<reco::GenJetCollection> genJets, 
		   unsigned maxJet);

    L1AnalysisGenJetDataFormat * getData() {return &genJet_;}
    void Reset() {genJet_.Reset();}

  private :
    L1AnalysisGenJetDataFormat genJet_;
  }; 
}
#endif


