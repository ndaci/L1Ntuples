#ifndef __L1Analysis_L1AnalysisGenMet_H__
#define __L1Analysis_L1AnalysisGenMet_H__

//-------------------------------------------------------------------------------
// Created 03/03/2010 - A.C. Le Bihan
// 
//
// Addition of met reco information
//-------------------------------------------------------------------------------

// standard
#include <vector>
// CMSSW
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "L1AnalysisGenMetDataFormat.h"
// Root
#include <TVector2.h>

namespace L1Analysis
{
  class L1AnalysisGenMet 
  {
  public:
    L1AnalysisGenMet();
    ~L1AnalysisGenMet();
    
    void SetMet(const edm::Handle<reco::GenMETCollection> genMet);
    void SetHtMht(const edm::Handle<reco::GenJetCollection> genJets, float jetptThreshold, std::vector<unsigned int> excludePdgId);

    L1AnalysisGenMetDataFormat * getData() {return &genMet_;}
    void Reset() {genMet_.Reset();}

  private :
    L1AnalysisGenMetDataFormat genMet_;
  }; 
}
#endif


