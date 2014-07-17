#ifndef RATEEFFICIENCY
#define RATEEFFICIENCY

#include "L1Ntuple.h"
#include "L1AlgoFactory.h"
#include <algorithm>
#include <map>

typedef map<TString,TH1F*> MAPHISTO1D;
typedef map<TString,TH2F*> MAPHISTO2D;

using namespace std;

class RateEfficiency : public L1Ntuple {
  
public:
  RateEfficiency(string filename) : L1Ntuple(filename) 
  {
    InitHistos1D();
    InitHistos2D();
  }

  RateEfficiency()  {}
  ~RateEfficiency() {}

  int run(bool runOnData, std::string resultTag, int minLs, int maxLs, 
	  float crossSec, float avPU, int nBunches, int isCrossSec, int nEvents = 0, bool doRate=true);
  
private:

  // Trigger methods
  int GetExtraInfo();
  int ScanMaxExtra();

  // Initialize and rescale histograms
  int   InitHistos1D();
  int   InitHistos2D();
  void  RescaleHistos1D(float scaleFactor);
  void  RescaleHistos2D(float scaleFactor);
  int   FillHistos();

  // Rescaling functions
  float ScaleFactor(float nZeroBias, float nBunches);
  void  setRateError(TH1F* histo);
  float computeAvgLumi(float xSec, float avPU, int nBunches) { return 11246. * avPU * nBunches / (1E7 * xSec); };

  // Histogram maps and iterators
  MAPHISTO1D hTH1F;
  MAPHISTO2D hTH2F;
  MAPHISTO1D::iterator it_hTH1F;
  MAPHISTO2D::iterator it_hTH2F;

  // Trigger event maxima
  int maxJetPt, maxJetPtC, maxETM, maxHTM;
  int maxJetPt_DPhi, maxJetPtC_DPhi;
  double maxDPhi_Jet_ETM, maxDPhi_Jet_HTM, maxDPhi_JetC_ETM, maxDPhi_JetC_HTM;
  double dphi;

  // L1Extra informations
  UInt_t          nIsoEm;
  vector<double>  isoEmEt;
  vector<double>  isoEmEta;
  vector<double>  isoEmPhi;
  vector<int>     isoEmBx;
  UInt_t          nNonIsoEm;
  vector<double>  nonIsoEmEt;
  vector<double>  nonIsoEmEta;
  vector<double>  nonIsoEmPhi;
  vector<int>     nonIsoEmBx;
  UInt_t          nCenJets;
  vector<double>  cenJetEt;
  vector<double>  cenJetEta;
  vector<double>  cenJetPhi;
  vector<int>     cenJetBx;
  UInt_t          nFwdJets;
  vector<double>  fwdJetEt;
  vector<double>  fwdJetEta;
  vector<double>  fwdJetPhi;
  vector<int>     fwdJetBx;
  UInt_t          nTauJets;
  vector<double>  tauJetEt;
  vector<double>  tauJetEta;
  vector<double>  tauJetPhi;
  vector<int>     tauJetBx;
  UInt_t          nMuons;
  vector<double>  muonEt;
  vector<double>  muonEta;
  vector<double>  muonPhi;
  vector<int>     muonChg;
  vector<unsigned int> muonIso;
  vector<unsigned int> muonFwd;
  vector<unsigned int> muonMip;
  vector<unsigned int> muonRPC;
  vector<int>     muonBx;
  vector<int>     muonQuality;
  vector<double>  hfEtSum;
  vector<unsigned int> hfBitCnt;
  vector<int>     hfBx;
  UInt_t          nMet;
  vector<double>  et;
  vector<double>  met;
  vector<double>  metPhi;
  vector<double>  metBx;
  UInt_t          nMht;
  vector<double>  ht;
  vector<double>  mht;
  vector<double>  mhtPhi;
  vector<double>  mhtBx;

};

#endif
