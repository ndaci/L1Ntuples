#include "L1TriggerDPG/L1Ntuples/interface/L1AnalysisGenJet.h"

//#include "JetMETCorrections/Objects/interface/JetCorrector.h"

L1Analysis::L1AnalysisGenJet::L1AnalysisGenJet()
{
}


L1Analysis::L1AnalysisGenJet::~L1AnalysisGenJet()
{
}


void L1Analysis::L1AnalysisGenJet::SetGenJet(const edm::Event& event,
                 const edm::EventSetup& setup,
                 edm::Handle<reco::GenJetCollection> genJets,
                 unsigned maxJet)
{

  genJet_.nJets=0;

  for(reco::GenJetCollection::const_iterator it=genJets->begin();
      it!=genJets->end() && genJet_.nJets < maxJet;
      ++it) {

    genJet_.et.push_back(it->et());
    genJet_.pt.push_back(it->pt());
    genJet_.eta.push_back(it->eta());
    genJet_.phi.push_back(it->phi());
    genJet_.e.push_back(it->energy());

    genJet_.eEM.push_back(it->emEnergy());
    genJet_.eHad.push_back(it->hadEnergy());
    genJet_.eInv.push_back(it->invisibleEnergy());
    genJet_.eAux.push_back(it->auxiliaryEnergy());
    //genJet_.etaDet.push_back(it->detectorEta());

    genJet_.nJets++;

  }
}

