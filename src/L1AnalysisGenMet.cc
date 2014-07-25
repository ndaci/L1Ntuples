#include "L1TriggerDPG/L1Ntuples/interface/L1AnalysisGenMet.h"

L1Analysis::L1AnalysisGenMet::L1AnalysisGenMet() 
{
}

L1Analysis::L1AnalysisGenMet::~L1AnalysisGenMet()
{
}

void L1Analysis::L1AnalysisGenMet::SetMet(const edm::Handle<reco::GenMETCollection> genMet)
{ 
  const reco::GenMETCollection *metCol = genMet.product();
  const reco::GenMET theMet = metCol->front();

  genMet_.met    = theMet.et();
  genMet_.metPhi = theMet.phi();
  genMet_.sumEt  = theMet.sumEt();

  genMet_.neutralEM_Et_fraction  = theMet.NeutralEMEtFraction();
  genMet_.chargedEM_Et_fraction  = theMet.ChargedEMEtFraction();
  genMet_.neutralHad_Et_fraction = theMet.NeutralHadEtFraction();
  genMet_.chargedHad_Et_fraction = theMet.ChargedHadEtFraction();
  genMet_.muon_Et_fraction       = theMet.MuonEtFraction();
  genMet_.inv_Et_fraction        = theMet.InvisibleEtFraction();

}

void L1Analysis::L1AnalysisGenMet::SetHtMht(const edm::Handle<reco::GenJetCollection> genJets, float jetptThreshold, std::vector<unsigned int> excludePdgId)
{  
  float mHx = 0.;
  float mHy = 0.;

  genMet_.Ht     = 0;
  genMet_.mHt    = -999;
  genMet_.mHtPhi = -999;

  // Check if genJet must be excluded
  bool skipJet=false;
  std::vector<const reco::GenParticle*> constituents;

  for (reco::GenJetCollection::const_iterator genjet = genJets->begin(); genjet!=genJets->end(); ++genjet)
  {
    // Check if genJet must be excluded
    skipJet=false;
    //std::vector<const reco::GenParticle*> constituents = genjet->getGenConstituents();
    if(excludePdgId.size()!=0) 
      constituents = genjet->getGenConstituents();
    //
    for(unsigned int iP=0 ; iP<excludePdgId.size() ; iP++) {
      for(unsigned int iC=0 ; iC<constituents.size() ; iC++) {
	if(std::abs(constituents[iC]->pdgId())==excludePdgId[iP]) {
	  skipJet=true;
	  break; // stop looping over this jet's constituents
	}
      }
      if(skipJet) break; // stop looping over PdgIds to exclude
    }
    constituents.clear(); // clear for next iteration
    if(skipJet) continue; // go to next genJet
    ///////////////////////////////////

    // Add genJet passing threshold to HT and MHT
    if (genjet->pt()>jetptThreshold){
      mHx += -1.*genjet->px();
      mHy += -1.*genjet->py();
      genMet_.Ht  += genjet->pt();
    }
  }

  TVector2 *tv2 = new TVector2(mHx,mHy);
  
  genMet_.mHt	= tv2->Mod();
  genMet_.mHtPhi= tv2->Phi();

}

