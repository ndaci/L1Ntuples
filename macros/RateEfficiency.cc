#include "RateEfficiency.hh"

int verbose=1;
int CUTETM=70; // L1_ETM70 seed from L1 group
int CUTPHI=2;
int FREQ=100;
int MAXETM=100;
int MAXJET=200;

int RateEfficiency::run(bool runOnData, string resultTag, int minLs, int maxLs, float crossSec, float avPU, int nBunches, int isCrossSec, int nEvents, bool doRate) {

  // Prepare outputs
  if(verbose>0) cout << "- Prepare outputs" << endl;
  system("mkdir -p results");
  string resultName = "results_" + resultTag + (isCrossSec ? "_XSEC" : "_RATE") + ".root";
  if(verbose>0) cout << "-- try to create file : " << "results/" + resultName << "..." << endl;
  TFile *outFile = new TFile(("results/" + resultName).c_str(),"recreate");
  if(verbose>0) cout << "-- outFile->cd()... " ;
  outFile->cd();
  if(verbose>0) cout << "done !" << endl;

  // Prepare histograms
  if(verbose>0) cout << "- Initialize Histograms : do 1D, " ;
  InitHistos1D();
  if(verbose>0) cout << " done. Do 2D, " ;
  InitHistos2D();
  if(verbose>0) cout << "done." << endl;

  // Process the inputs
  float nZeroBias = 0;
  int nevents = nEvents == 0 ? GetEntries() : nEvents;
  if(verbose>0) cout << "Running on " << nevents << " events." << endl;
  if(runOnData) if(verbose>0) cout << "Run on data" << endl;

  // Loop over input tree
  for (Long64_t event=0; event<nevents; ++event) { 

    // Get current entry
    Long64_t eventEntry = LoadTree(event); 
    if (eventEntry < 0) break;
    GetEntry(event);

    // Warn regularly
    if (event%FREQ == 0 && verbose>0) cout << "Processed " << event << " events." << endl;

    // Limits in terms of LS
    if ( event_->lumi < minLs || event_->lumi > maxLs ) continue;

    // Pile-up reweighting
    double weight = event_->puWeight > -0.001 ? event_->puWeight : 1; 
    nZeroBias += weight;

    // Compute triggers and fill histograms
    GetExtraInfo();
    ScanMaxExtra();
    FillHistos();

  } // end loop over entries

  // Rescale histograms
  if(verbose>0) cout << "# of zero bias events (weighted) used for rate computation : " << nZeroBias << endl;
  //
  float scaleFactor = ScaleFactor(nZeroBias,nBunches);
  //
  if (isCrossSec) 
    scaleFactor /= (computeAvgLumi(crossSec,avPU,nBunches)*10000) ; // CB lumi is in 1E34 units
  //
  if(!doRate) scaleFactor = nZeroBias!=0 ? 1/float(nZeroBias) : 1;
  //
  RescaleHistos1D(scaleFactor);
  RescaleHistos2D(scaleFactor);

  // Write and close output file
  outFile->Write();
  outFile->Close();
  delete outFile;

  return 0;
}

int RateEfficiency::ScanMaxExtra()
{
  
  maxJetPt = maxJetPtC = maxJetPt_DPhi = maxJetPtC_DPhi = 0;
  
  // Scan central jets
  for(unsigned int iJ=0 ; iJ<nCenJets ; iJ++) {
    if(cenJetEt[iJ]>maxJetPtC)
      maxJetPtC = cenJetEt[iJ];
    if(cenJetEt[iJ]>maxJetPtC_DPhi && TMath::Abs(cenJetPhi[iJ]-metPhi[0])>CUTPHI )
      maxJetPtC_DPhi = cenJetEt[iJ];
  }
  
  // Save central jets contribution to {all jets} set
  maxJetPt = maxJetPtC;
  maxJetPt_DPhi = maxJetPtC_DPhi;

  // Scan forward jets
  for(unsigned int iJ=0 ; iJ<nFwdJets ; iJ++) {
    if(fwdJetEt[iJ]>maxJetPt)
      maxJetPt = fwdJetEt[iJ];
    if(fwdJetEt[iJ]>maxJetPt_DPhi && TMath::Abs(fwdJetPhi[iJ]-metPhi[0])>CUTPHI ) 
      maxJetPt_DPhi = fwdJetEt[iJ];
  }
  
  // Determine max met/mht in event
  maxETM=0;
  for(unsigned int iETM=0 ; iETM<nMet ; iETM++)
    if(met[iETM]>maxETM) maxETM = met[iETM];
  
  maxHTM=0;
  for(unsigned int iHTM=0 ; iHTM<nMht ; iHTM++)
    if(mht[iHTM]>maxHTM) maxHTM = mht[iHTM];
  
  if(verbose>1) 
    cout << "-- maxima :"
	 << " maxJetPt="  << maxJetPt
	 << " maxJetPtC=" << maxJetPtC
	 << " maxETM="    << maxETM
	 << " maxHTM="    << maxHTM
	 << endl;
  
  // DeltaPhi Study //
  
  // Initialize
  maxDPhi_Jet_ETM = maxDPhi_Jet_HTM = maxDPhi_JetC_ETM = maxDPhi_JetC_HTM = 0;
  dphi = 0;
  
  // Scan all jet eT cut values
  for(int iJetCut=0 ; iJetCut<=MAXJET ; iJetCut++) {
    
    // Look central jets
    if(verbose>1) cout << "--- Look central jets : iJetCut=" << iJetCut << endl;
    for(uint iJ=0 ; iJ<nCenJets ; iJ++) {
      if( cenJetEt[iJ]>=iJetCut ) {
	dphi = TMath::Abs( cenJetPhi[iJ] - metPhi[0] );
	if(dphi > maxDPhi_Jet_ETM ) maxDPhi_Jet_ETM  = dphi;
	if(dphi > maxDPhi_JetC_ETM) maxDPhi_JetC_ETM = dphi;
	dphi = TMath::Abs( cenJetPhi[iJ] - mhtPhi[0] );
	if(dphi > maxDPhi_Jet_HTM ) maxDPhi_Jet_HTM  = dphi;
	if(dphi > maxDPhi_JetC_HTM) maxDPhi_JetC_HTM = dphi;
      }
    }
    
    // Look forward jets
    if(verbose>1) cout << "--- Look forward jets : iJetCut=" << iJetCut << endl;
    for(uint iJ=0 ; iJ<nFwdJets ; iJ++) {
      if(verbose>2) cout << "---- iJ=" << iJ << endl;
      if( fwdJetEt[iJ]>=iJetCut ) {
	dphi = TMath::Abs( fwdJetPhi[iJ] - metPhi[0] );
	if(dphi > maxDPhi_Jet_ETM ) maxDPhi_Jet_ETM  = dphi;
	dphi = TMath::Abs( fwdJetPhi[iJ] - mhtPhi[0] );
	if(dphi > maxDPhi_Jet_HTM ) maxDPhi_Jet_HTM  = dphi;
      }
    }

    // Fill DeltaPhi histograms
    if(verbose>1) cout << "--- Fill DeltaPhi iJetCut=" << iJetCut << endl;

    hTH2F["hMaxDPhi_Jet_ETM"]  -> Fill( iJetCut, maxDPhi_Jet_ETM ); 
    hTH2F["hMaxDPhi_Jet_HTM"]  -> Fill( iJetCut, maxDPhi_Jet_HTM ); 
    hTH2F["hMaxDPhi_JetC_ETM"] -> Fill( iJetCut, maxDPhi_JetC_ETM ); 
    hTH2F["hMaxDPhi_JetC_HTM"] -> Fill( iJetCut, maxDPhi_JetC_HTM ); 
  }

  return 0;
}

void RateEfficiency::RescaleHistos1D(float scaleFactor)
{
  for(it_hTH1F=hTH1F.begin() ; it_hTH1F!=hTH1F.end() ; it_hTH1F++) {
    setRateError(it_hTH1F->second);
    it_hTH1F->second->Scale(scaleFactor);
  }
}

void RateEfficiency::RescaleHistos2D(float scaleFactor)
{
  for(it_hTH2F=hTH2F.begin() ; it_hTH2F!=hTH2F.end() ; it_hTH2F++) {
    it_hTH2F->second->Scale(scaleFactor);
  }
}

// scale factor computed w.r.t. ZeroBias rate fratcion and # bunches 
float RateEfficiency::ScaleFactor(float nZeroBias, float nBunches) {

  float scal = 11246.; // ZB per bunch in Hz
  scal /= nZeroBias;
  scal *= nBunches;

  return scal;
}

void RateEfficiency::setRateError(TH1F* histo) {

  int nBins = histo->GetNbinsX();

  for (int iBin=1; iBin<=nBins; ++iBin) {
    float value = histo->GetBinContent(iBin);
    float error = sqrt(value);

    histo->SetBinError(iBin,error);
  }

}

int RateEfficiency::InitHistos1D()
{

  const uint nH=4;
  const uint nA=2;

  TString name_histos[ nH] = {"hEff_Jet", "hEff_JetC", "hEff_ETM", "hEff_HTM"};
  TString title_histos[nH] = {"L1_Jet efficiency", "L1_JetC efficiency", 
			      "L1_ETM efficiency", "L1_HTM efficiency"};

  TString xytitle[nH][nA] = { {"Jet eT cut [GeV]", "Efficiency"},
			      {"Central Jet eT cut [GeV]", "Efficiency"},
			      {"ETM cut [GeV]", "Efficiency"},
			      {"HTM cut [GeV]", "Efficiency"} };

  Int_t nbins[nH] = {MAXJET,MAXJET,MAXETM,MAXETM};
  Int_t xmax[ nH] = {MAXJET,MAXJET,MAXETM,MAXETM};
  Int_t xmin[ nH] = {0,0,0,0};

  for(uint iH=0 ; iH<nH ; iH++) {
    hTH1F[name_histos[iH]] = new TH1F(name_histos[iH],title_histos[iH],
				      nbins[iH],xmin[iH],xmax[iH]);
    hTH1F[name_histos[iH]]->SetXTitle(xytitle[iH][0]);
    hTH1F[name_histos[iH]]->SetYTitle(xytitle[iH][1]);
  }

  return 0;
}

int RateEfficiency::InitHistos2D()
{

  const int nA=2;
  const int nH=16;
  TString name_histos[nH] = {"hEff_Jet_ETM", "hEff_Jet_HTM", "hEff_JetC_ETM", "hEff_JetC_HTM",
			     "hEff_Jet_ETM_DPhi", "hEff_Jet_HTM_DPhi", "hEff_JetC_ETM_DPhi", "hEff_JetC_HTM_DPhi",
			     "hEff_Jet_ETM_OR70", "hEff_JetC_ETM_OR70", "hEff_Jet_ETM_DPhi_OR70", "hEff_JetC_ETM_DPhi_OR70",
			     "hMaxDPhi_Jet_ETM","hMaxDPhi_Jet_HTM","hMaxDPhi_JetC_ETM","hMaxDPhi_JetC_HTM"};

  TString title_histos[nH] = {"L1_Jet_ETM efficiency" ,     "L1_Jet_HTM efficiency", 
			      "L1_JetC_ETM efficiency",     "L1_JetC_HTM efficiency",
			      "L1_Jet_ETM_DPhi efficiency", "L1_Jet_HTM_DPhi efficiency",
			      "L1_JetC_ETM_DPhi efficiency","L1_JetC_HTM_DPhi efficiency",
			      "(L1_Jet_ETM || L1_ETM70) efficiency",
			      "(L1_JetC_ETM || L1_ETM70) efficiency",
			      "(L1_Jet_ETM_DPhi || L1_ETM70) efficiency",
			      "(L1_JetC_ETM_DPhi || L1_ETM70) efficiency",
			      "Max | #Delta#Phi(Jet, ETM) |", "Max | #Delta#Phi(Jet, HTM) |",
			      "Max | #Delta#Phi(JetC, ETM) |","Max | #Delta#Phi(JetC, HTM) |"};

  TString xytitle[nH][nA] = { {"Jet eT cut [GeV]", "ETM cut [GeV]"},
			      {"Jet eT cut [GeV]", "HTM cut [GeV]"},
			      {"Central Jet eT cut [GeV]", "ETM cut [GeV]"},
			      {"Central Jet eT cut [GeV]", "HTM cut [GeV]"},
			      {"Jet eT cut [GeV]", "ETM cut [GeV]"},
			      {"Jet eT cut [GeV]", "HTM cut [GeV]"},
			      {"Central Jet eT cut [GeV]", "ETM cut [GeV]"},
			      {"Central Jet eT cut [GeV]", "HTM cut [GeV]"},
			      {"Jet eT cut [GeV]", "ETM cut [GeV]"},
			      {"Central Jet eT cut [GeV]", "ETM cut [GeV]"},
			      {"Jet eT cut [GeV]", "ETM cut [GeV]"},
			      {"Central Jet eT cut [GeV]", "ETM cut [GeV]"},
			      {"Jet eT cut [GeV]", "Max | #Delta#Phi(Jet, ETM) |"},
			      {"Jet eT cut [GeV]", "Max | #Delta#Phi(Jet, HTM) |"},
			      {"Central Jet eT cut [GeV]", "Max | #Delta#Phi(JetC, ETM) |"},
			      {"Central Jet eT cut [GeV]", "Max | #Delta#Phi(JetC, HTM) |"} };

  Int_t nbins[2]={MAXJET, MAXETM};
  Int_t xmax[2] ={MAXJET, MAXETM};
  Int_t xmin =0;

  for(uint iH=0 ; iH<nH ; iH++) {
    hTH2F[name_histos[iH]] = new TH2F(name_histos[iH],title_histos[iH], 
				      nbins[0],xmin,xmax[0],
				      nbins[1],xmin,xmax[1]);
    hTH2F[name_histos[iH]]->SetXTitle(xytitle[iH][0]);
    hTH2F[name_histos[iH]]->SetYTitle(xytitle[iH][1]);
  }

  return 0;
}

int RateEfficiency::GetExtraInfo()
{

  nIsoEm        = l1extra_-> nIsoEm;       
  isoEmEt	= l1extra_-> isoEmEt;	
  isoEmEta	= l1extra_-> isoEmEta;	
  isoEmPhi	= l1extra_-> isoEmPhi;	
  isoEmBx	= l1extra_-> isoEmBx;	
  nNonIsoEm	= l1extra_-> nNonIsoEm;	
  nonIsoEmEt	= l1extra_-> nonIsoEmEt;	
  nonIsoEmEta	= l1extra_-> nonIsoEmEta;	
  nonIsoEmPhi	= l1extra_-> nonIsoEmPhi;	
  nonIsoEmBx	= l1extra_-> nonIsoEmBx;	
  nCenJets	= l1extra_-> nCenJets;	
  cenJetEt	= l1extra_-> cenJetEt;	
  cenJetEta	= l1extra_-> cenJetEta;	
  cenJetPhi	= l1extra_-> cenJetPhi;	
  cenJetBx	= l1extra_-> cenJetBx;	
  nFwdJets	= l1extra_-> nFwdJets;	
  fwdJetEt	= l1extra_-> fwdJetEt;	
  fwdJetEta	= l1extra_-> fwdJetEta;	
  fwdJetPhi	= l1extra_-> fwdJetPhi;	
  fwdJetBx	= l1extra_-> fwdJetBx;	
  nTauJets	= l1extra_-> nTauJets;	
  tauJetEt	= l1extra_-> tauJetEt;	
  tauJetEta	= l1extra_-> tauJetEta;	
  tauJetPhi	= l1extra_-> tauJetPhi;	
  tauJetBx	= l1extra_-> tauJetBx;	
  nMuons	= l1extra_-> nMuons;	
  muonEt	= l1extra_-> muonEt;	
  muonEta	= l1extra_-> muonEta;	
  muonPhi	= l1extra_-> muonPhi;	
  muonChg	= l1extra_-> muonChg;	
  muonIso	= l1extra_-> muonIso;	
  muonFwd	= l1extra_-> muonFwd;	
  muonMip	= l1extra_-> muonMip;	
  muonRPC	= l1extra_-> muonRPC;	
  muonBx	= l1extra_-> muonBx;	
  muonQuality	= l1extra_-> muonQuality;	
  hfEtSum	= l1extra_-> hfEtSum;	
  hfBitCnt	= l1extra_-> hfBitCnt;	
  hfBx		= l1extra_-> hfBx;		
  nMet		= l1extra_-> nMet;		
  et		= l1extra_-> et;		
  met		= l1extra_-> met;		
  metPhi	= l1extra_-> metPhi;	
  metBx	        = l1extra_-> metBx;	
  nMht		= l1extra_-> nMht;		
  ht		= l1extra_-> ht;		
  mht		= l1extra_-> mht;		
  mhtPhi	= l1extra_-> mhtPhi;	
  mhtBx  	= l1extra_-> mhtBx;  	
  
  return 0;
}

int RateEfficiency::FillHistos()
{
  
  // Fill histograms
  for(int iJetPt=0 ; iJetPt<=MAXJET ; iJetPt++) {
    if(iJetPt<=maxJetPt)  hTH1F["hEff_Jet"]  -> Fill(iJetPt);
    if(iJetPt<=maxJetPtC) hTH1F["hEff_JetC"] -> Fill(iJetPt);
    
    for(int iETM=0 ; iETM<=MAXETM ; iETM++) {

      // Jet(C)+ETM (+DPhi)
      if(iJetPt<=maxJetPt        && iETM<=maxETM) hTH2F["hEff_Jet_ETM"]       -> Fill(iJetPt, iETM);
      if(iJetPt<=maxJetPtC       && iETM<=maxETM) hTH2F["hEff_JetC_ETM"]      -> Fill(iJetPt, iETM);
      if(iJetPt<=maxJetPt_DPhi   && iETM<=maxETM) hTH2F["hEff_Jet_ETM_DPhi"]  -> Fill(iJetPt, iETM);
      if(iJetPt<=maxJetPtC_DPhi  && iETM<=maxETM) hTH2F["hEff_JetC_ETM_DPhi"] -> Fill(iJetPt, iETM);

      // Jet(C)+ETM (+DPhi) OR ETM70
      if((iJetPt<=maxJetPt       && iETM<=maxETM) || maxETM>=CUTETM) hTH2F["hEff_Jet_ETM_OR70"] -> Fill(iJetPt, iETM);
      if((iJetPt<=maxJetPtC      && iETM<=maxETM) || maxETM>=CUTETM) hTH2F["hEff_JetC_ETM_OR70"]-> Fill(iJetPt, iETM);
      if((iJetPt<=maxJetPt_DPhi  && iETM<=maxETM) || maxETM>=CUTETM) hTH2F["hEff_Jet_ETM_DPhi_OR70"]-> Fill(iJetPt, iETM);
      if((iJetPt<=maxJetPtC_DPhi && iETM<=maxETM) || maxETM>=CUTETM) hTH2F["hEff_JetC_ETM_DPhi_OR70"]-> Fill(iJetPt, iETM);

      // Jet(C)+HTM (+DPhi)
      if(iJetPt<=maxJetPt       && iETM<=maxHTM) hTH2F["hEff_Jet_HTM"] -> Fill(iJetPt, iETM);
      if(iJetPt<=maxJetPtC      && iETM<=maxHTM) hTH2F["hEff_JetC_HTM"]-> Fill(iJetPt, iETM);
      if(iJetPt<=maxJetPt_DPhi  && iETM<=maxHTM) hTH2F["hEff_Jet_HTM_DPhi"] -> Fill(iJetPt, iETM);
      if(iJetPt<=maxJetPtC_DPhi && iETM<=maxHTM) hTH2F["hEff_JetC_HTM_DPhi"]-> Fill(iJetPt, iETM);
    }
  }
  
  for(int iETM=0 ; iETM<=maxHTM || iETM<=maxETM ; iETM++) {
    if(iETM<=maxETM) hTH1F["hEff_ETM"] -> Fill(iETM);
    if(iETM<=maxHTM) hTH1F["hEff_HTM"] -> Fill(iETM);
  }

  return 0;
}

