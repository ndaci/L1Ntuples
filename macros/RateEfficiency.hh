#ifndef RATEEFFICIENCY
#define RATEEFFICIENCY

#include "L1Ntuple.h"
#include "L1AlgoFactory.h"
#include <algorithm>
#include <map>

#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TVector2.h"
#include "TLorentzVector.h"

typedef map<TString,TH1F*> MAPHISTO1D;
typedef map<TString,TH2F*> MAPHISTO2D;
typedef map<TString,TH3F*> MAPHISTO3D;

using namespace std;

// CLASS //
class RateEfficiency : public L1Ntuple {
  
public:
  RateEfficiency(string filename) : L1Ntuple(filename) 
  {
  }

  RateEfficiency()  {}
  ~RateEfficiency() {}

  Int_t run(bool runOnData, std::string resultTag, Int_t minLs, Int_t maxLs, 
	    Float_t crossSec, Float_t avPU, Int_t nBunches, Int_t isCrossSec, Int_t nEvents = 0, bool doRate=true,
	    UInt_t iPart=0, UInt_t nPart=1);
  
private:

  // Trigger methods
  Int_t GetExtraInfo();
  Int_t OrderJets();
  Int_t ScanMaxExtra();
  Int_t EvalAlgos();
  Int_t RescaleAlgos(Float_t scaleFactor);
  Int_t WriteValues(ofstream &outstream);

  // Geometry
  Float_t computeDeltaPhi(Float_t phi1, Float_t phi2);
  Float_t computeDeltaR(  Float_t phi1, Float_t phi2, Float_t eta1, Float_t eta2);

  // Initialize members
  Int_t InitVar();

  // Initialize and rescale histograms
  Int_t   InitHistos1D();
  Int_t   InitHistos2D();
  Int_t   InitHistos3D();
  void  RescaleHistos1D(Float_t scaleFactor);
  void  RescaleHistos2D(Float_t scaleFactor);
  void  RescaleHistos3D(Float_t scaleFactor);
  Int_t   FillHistos();

  // Rescaling functions
  Float_t ScaleFactor(Float_t nZeroBias, Float_t nBunches);
  void  setRateError(TH1F* histo);
  Float_t computeAvgLumi(Float_t xSec, Float_t avPU, Int_t nBunches) { return 11246. * avPU * nBunches / (1E7 * xSec); };

  // Histogram maps and iterators
  MAPHISTO1D hTH1F;
  MAPHISTO2D hTH2F;
  MAPHISTO3D hTH3F;
  MAPHISTO1D::iterator it_hTH1F;
  MAPHISTO2D::iterator it_hTH2F;
  MAPHISTO3D::iterator it_hTH3F;

  // Trigger event maxima
  Int_t maxJetPt, maxJetPtC, maxETM, maxHTM;
  Int_t maxJetPt_DPhi, maxJetPtC_DPhi;
  Double_t maxDPhi_Jet_ETM, maxDPhi_Jet_HTM, maxDPhi_JetC_ETM, maxDPhi_JetC_HTM;
  Double_t dphi, mindphi;

  vector<Double_t> minDPhi_Jet_ETM, minDPhi_Jet_HTM, minDPhi_JetC_ETM, minDPhi_JetC_HTM;
  vector<UInt_t> iMiniFirst, iMiniLast;

  // Counters for algos
  Float_t n_ETM60_NoQCD2_OR_ETM70, n_ETM58_NoQCD4_OR_ETM70, 
    n_ETM60_NoQCD2_OR_ETM80, n_ETM60_NoQCD4_OR_ETM80, 
    n_ETM65_NoQCD2_OR_ETM75, n_ETM65_NoQCD4_OR_ETM75;

  // L1Extra informations
  UInt_t          nAllJets;
  vector<Double_t>  allJetEt;
  vector<Double_t>  allJetEta;
  vector<Double_t>  allJetPhi;
  vector<Int_t>     allJetBx;

  UInt_t          nIsoEm;
  vector<Double_t>  isoEmEt;
  vector<Double_t>  isoEmEta;
  vector<Double_t>  isoEmPhi;
  vector<Int_t>     isoEmBx;

  UInt_t          nNonIsoEm;
  vector<Double_t>  nonIsoEmEt;
  vector<Double_t>  nonIsoEmEta;
  vector<Double_t>  nonIsoEmPhi;
  vector<Int_t>     nonIsoEmBx;

  UInt_t          nCenJets;
  vector<Double_t>  cenJetEt;
  vector<Double_t>  cenJetEta;
  vector<Double_t>  cenJetPhi;
  vector<Int_t>     cenJetBx;

  UInt_t          nFwdJets;
  vector<Double_t>  fwdJetEt;
  vector<Double_t>  fwdJetEta;
  vector<Double_t>  fwdJetPhi;
  vector<Int_t>     fwdJetBx;

  UInt_t          nTauJets;
  vector<Double_t>  tauJetEt;
  vector<Double_t>  tauJetEta;
  vector<Double_t>  tauJetPhi;
  vector<Int_t>     tauJetBx;

  UInt_t          nMuons;
  vector<Double_t>  muonEt;
  vector<Double_t>  muonEta;
  vector<Double_t>  muonPhi;
  vector<Int_t>     muonChg;
  vector<UInt_t> muonIso;
  vector<UInt_t> muonFwd;
  vector<UInt_t> muonMip;
  vector<UInt_t> muonRPC;
  vector<Int_t>     muonBx;
  vector<Int_t>     muonQuality;

  vector<Double_t>  hfEtSum;
  vector<UInt_t> hfBitCnt;
  vector<Int_t>     hfBx;

  UInt_t    nMet;
  vector<Double_t>  et;
  vector<Double_t>  met;
  vector<Double_t>  metPhi;
  vector<Double_t>  metBx;

  UInt_t          nMht;
  vector<Double_t>  ht;
  vector<Double_t>  mht;
  vector<Double_t>  mhtPhi;
  vector<Double_t>  mhtBx;

};

struct HighestPt{
  bool operator()( TLorentzVector v1, TLorentzVector v2 ) const{
    return v1.Pt() > v2.Pt();
  }
};

#endif
