#include "RateEfficiency.hh"

int verbose=-1;
int CUTETM=70; // L1_ETM70 seed from L1 group
int CUTPHI=2.0;
int CUTPHI_H=1;
int FREQ=100;
int MAXETM=100;
int MAXJET=200;

int RateEfficiency::run(bool runOnData, string resultTag, int minLs, int maxLs, float crossSec, float avPU, int nBunches, int isCrossSec, int nEvents, bool doRate) {

  cout << "∫=== THE BEGINNING ===∫" << endl;

  // Prepare outputs
  if(verbose>0) cout << "- Prepare outputs" << endl;
  //system("mkdir -p results");
  TString resultName = "results_" + resultTag + (isCrossSec ? "_XSEC" : "_RATE");
  TFile *outFile = new TFile("results/"+resultName+".root","recreate");
  outFile->cd();
  ofstream outstream("results/"+resultName+".txt", ios::out);

  // Prepare histograms
  if(verbose>0) cout << "- Initialize Histograms : do 1D, " ;
  InitHistos1D();
  if(verbose>0) cout << " done. Do 2D, " ;
  InitHistos2D();
  if(verbose>0) cout << " done. Do 3D, " ;
  InitHistos3D();
  if(verbose>0) cout << "done." << endl;

  // Initialize algo counters
  n_ETM60_NoQCD2_OR_ETM70 = n_ETM58_NoQCD4_OR_ETM70 = 
    n_ETM60_NoQCD2_OR_ETM80 = n_ETM60_NoQCD4_OR_ETM80 = 
    n_ETM65_NoQCD2_OR_ETM75 = n_ETM65_NoQCD4_OR_ETM75 = 0;

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
    InitVar();
    GetExtraInfo();
    OrderJets();
    EvalAlgos();
    ScanMaxExtra();
    FillHistos();
    InitVar();

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
  if(verbose>1) cout << "- scaleFactor=" << scaleFactor << endl;
  //
  if(verbose>1) cout << "- Rescale histos" << endl;
  RescaleHistos1D(scaleFactor);
  RescaleHistos2D(scaleFactor);
  RescaleHistos3D(scaleFactor);
  //
  if(verbose>1) cout << "- Rescale algos" << endl;
  RescaleAlgos(scaleFactor);
  if(verbose>1) cout << "- Write values" << endl;
  WriteValues(outstream);
  //

  // Write and close output file
  if(verbose>1) cout << "- Write outFile" << endl;
  outFile->Write();

  if(verbose>1) cout << "- Close outFile" << endl;
  outFile->Close();

  if(verbose>1) cout << "- Delete outFile" << endl;
  delete outFile;

  if(verbose>1) cout << "∫=== THE END ===∫" << endl;

  return 0;
}

int RateEfficiency::WriteValues(ofstream &outstream)
{
  outstream << "n_ETM60_NoQCD2_OR_ETM70=" << n_ETM60_NoQCD2_OR_ETM70  << endl;
  outstream << "n_ETM58_NoQCD4_OR_ETM70=" << n_ETM58_NoQCD4_OR_ETM70  << endl;
  outstream << "n_ETM60_NoQCD2_OR_ETM80=" << n_ETM60_NoQCD2_OR_ETM80  << endl;
  outstream << "n_ETM60_NoQCD4_OR_ETM80=" << n_ETM60_NoQCD4_OR_ETM80  << endl;
  outstream << "n_ETM65_NoQCD2_OR_ETM75=" << n_ETM65_NoQCD2_OR_ETM75  << endl;
  outstream << "n_ETM65_NoQCD4_OR_ETM75=" << n_ETM65_NoQCD4_OR_ETM75  << endl;

  return 0;
}

int RateEfficiency::RescaleAlgos(float scaleFactor)
{
  
  n_ETM60_NoQCD2_OR_ETM70  = n_ETM60_NoQCD2_OR_ETM70 * scaleFactor ;
  n_ETM58_NoQCD4_OR_ETM70  = n_ETM58_NoQCD4_OR_ETM70 * scaleFactor ;
  n_ETM60_NoQCD2_OR_ETM80  = n_ETM60_NoQCD2_OR_ETM80 * scaleFactor ;
  n_ETM60_NoQCD4_OR_ETM80  = n_ETM60_NoQCD4_OR_ETM80 * scaleFactor ;
  n_ETM65_NoQCD2_OR_ETM75  = n_ETM65_NoQCD2_OR_ETM75 * scaleFactor ;
  n_ETM65_NoQCD4_OR_ETM75  = n_ETM65_NoQCD4_OR_ETM75 * scaleFactor ;

  return 0;
}

int RateEfficiency::EvalAlgos()
{

  bool veto2=false;
  bool veto4=false;

  for(UInt_t iJ=0 ; iJ<nCenJets ; iJ++) {
    dphi = computeDeltaPhi(metPhi[0] , cenJetPhi[iJ]);
    if(dphi<CUTPHI_H) {
      if(iJ<2) veto2=true;
      veto4=true;
    }
  }

  for(UInt_t iJ=0 ; iJ<nFwdJets ; iJ++) {
    dphi = computeDeltaPhi(metPhi[0] , fwdJetPhi[iJ]);
    if(dphi<CUTPHI_H) {
      if(iJ<2) veto2=true;
      veto4=true;
    }
  }

  // Test L1 algorithms
  //
  if( (!veto2 && met[0]>=60) || met[0]>=70) n_ETM60_NoQCD2_OR_ETM70 ++ ;
  if( (!veto4 && met[0]>=58) || met[0]>=70) n_ETM58_NoQCD4_OR_ETM70 ++ ;
  //
  if( (!veto2 && met[0]>=60) || met[0]>=80) n_ETM60_NoQCD2_OR_ETM80 ++ ;
  if( (!veto4 && met[0]>=60) || met[0]>=80) n_ETM60_NoQCD4_OR_ETM80 ++ ;
  //
  if( (!veto2 && met[0]>=65) || met[0]>=75) n_ETM65_NoQCD2_OR_ETM75 ++ ;
  if( (!veto4 && met[0]>=65) || met[0]>=75) n_ETM65_NoQCD4_OR_ETM75 ++ ;

  return 0;
}

int RateEfficiency::ScanMaxExtra()
{
  
  maxJetPt = maxJetPtC = maxJetPt_DPhi = maxJetPtC_DPhi = 0;
  
  // Scan central jets
  for(unsigned int iJ=0 ; iJ<nCenJets ; iJ++) {
    if(cenJetEt[iJ]>maxJetPtC)
      maxJetPtC = cenJetEt[iJ];
    if(cenJetEt[iJ]>maxJetPtC_DPhi && computeDeltaPhi(cenJetPhi[iJ], metPhi[0])>CUTPHI )
      maxJetPtC_DPhi = cenJetEt[iJ];
  }
  
  // Save central jets contribution to {all jets} set
  maxJetPt = maxJetPtC;
  maxJetPt_DPhi = maxJetPtC_DPhi;

  // Scan forward jets
  for(unsigned int iJ=0 ; iJ<nFwdJets ; iJ++) {
    if(fwdJetEt[iJ]>maxJetPt)
      maxJetPt = fwdJetEt[iJ];
    if(fwdJetEt[iJ]>maxJetPt_DPhi && computeDeltaPhi(fwdJetPhi[iJ], metPhi[0])>CUTPHI ) 
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
  dphi = 0;
  maxDPhi_Jet_ETM = maxDPhi_Jet_HTM = maxDPhi_JetC_ETM = maxDPhi_JetC_HTM = -1;

  mindphi = 1000000000;
  for(UInt_t iL=0 ; iL<4 ; iL++)
    minDPhi_Jet_ETM[iL] = minDPhi_Jet_HTM[iL] = minDPhi_JetC_ETM[iL] = minDPhi_JetC_HTM[iL] = 1000000000;
  
  // Scan all jet eT cut values
  for(int iJetCut=0 ; iJetCut<=MAXJET ; iJetCut++) {

    dphi = 0;
    maxDPhi_Jet_ETM = maxDPhi_Jet_HTM = maxDPhi_JetC_ETM = maxDPhi_JetC_HTM = -1;

    mindphi = 1000000000;
    for(UInt_t iL=0 ; iL<4 ; iL++)
      minDPhi_Jet_ETM[iL] = minDPhi_Jet_HTM[iL] = minDPhi_JetC_ETM[iL] = minDPhi_JetC_HTM[iL] = 1000000000;

    // Look all jets
    for(UInt_t iJ=0 ; iJ<nAllJets ; iJ++) {
      if( allJetEt[iJ]>=iJetCut ) {
	mindphi = computeDeltaPhi(allJetPhi[iJ], metPhi[0]);
	for(UInt_t iL=0 ; iL<4 ; iL++)
	  if(iJ<=iL && mindphi < minDPhi_Jet_ETM[iL]) minDPhi_Jet_ETM[iL] = mindphi;


	mindphi = computeDeltaPhi(allJetPhi[iJ], mhtPhi[0]);
	for(UInt_t iL=0 ; iL<4 ; iL++)
	  if(iJ<=iL && mindphi < minDPhi_Jet_HTM[iL]) minDPhi_Jet_HTM[iL] = mindphi;	
      }
    }
    
    // Look central jets
    if(verbose>1) cout << "--- Look central jets : iJetCut=" << iJetCut << endl;
    for(uint iJ=0 ; iJ<nCenJets ; iJ++) {

      if( cenJetEt[iJ]>=iJetCut ) {

	dphi = computeDeltaPhi(cenJetPhi[iJ], metPhi[0]);
	if(dphi > maxDPhi_Jet_ETM ) maxDPhi_Jet_ETM  = dphi;
	if(dphi > maxDPhi_JetC_ETM) maxDPhi_JetC_ETM = dphi;
	//
	mindphi = dphi;
	for(UInt_t iL=0 ; iL<4 ; iL++) {
	  if(iJ<=iL && mindphi < minDPhi_JetC_ETM[iL]) minDPhi_JetC_ETM[iL] = mindphi;
	}

	dphi = computeDeltaPhi(cenJetPhi[iJ], mhtPhi[0]);
	if(dphi > maxDPhi_Jet_HTM ) maxDPhi_Jet_HTM  = dphi;
	if(dphi > maxDPhi_JetC_HTM) maxDPhi_JetC_HTM = dphi;
	//
	mindphi = dphi;
	for(UInt_t iL=0 ; iL<4 ; iL++) {
	  if(iJ<=iL && mindphi < minDPhi_JetC_HTM[iL]) minDPhi_JetC_HTM[iL] = mindphi;
	}
	
      }
    }
    
    // Look forward jets
    if(verbose>1) cout << "--- Look forward jets : iJetCut=" << iJetCut << endl;
    for(uint iJ=0 ; iJ<nFwdJets ; iJ++) {
      if(verbose>2) cout << "---- iJ=" << iJ << endl;
      if( fwdJetEt[iJ]>=iJetCut ) {
	dphi = computeDeltaPhi(fwdJetPhi[iJ], metPhi[0]);
	if(dphi > maxDPhi_Jet_ETM ) maxDPhi_Jet_ETM  = dphi;
	//
	dphi = computeDeltaPhi(fwdJetPhi[iJ], mhtPhi[0]);
	if(dphi > maxDPhi_Jet_HTM ) maxDPhi_Jet_HTM  = dphi;
      }
    }

    // Printout values
    if(verbose>0) cout << "--- Values : min dphi for iJetCut=" << iJetCut << endl;
    for(UInt_t iL=0 ; iL<4 ; iL++) {
      if(verbose>0) cout << "minDPhi_Jet_ETM["  << iL << "]=" << minDPhi_Jet_ETM[iL]  << " ; "
			 << "minDPhi_JetC_ETM[" << iL << "]=" << minDPhi_JetC_ETM[iL] << " ; "
			 << "minDPhi_Jet_HTM["  << iL << "]=" << minDPhi_Jet_HTM[iL]  << " ; "
			 << "minDPhi_JetC_HTM[" << iL << "]=" << minDPhi_JetC_HTM[iL] << " ; "
			 << endl;
    }

    // Fill DeltaPhi histograms
    if(verbose>1) cout << "--- Fill DeltaPhi iJetCut=" << iJetCut << endl;
    //for(int iETM=0 ; iETM<MAXETM ; iETM++) {
      //if(met[0]>=iETM) {
    if(maxDPhi_Jet_ETM >=0) hTH3F["hMaxDPhi_Jet_ETM"]  -> Fill( iJetCut, met[0], maxDPhi_Jet_ETM ); 
    if(maxDPhi_JetC_ETM>=0) hTH3F["hMaxDPhi_JetC_ETM"] -> Fill( iJetCut, met[0], maxDPhi_JetC_ETM ); 

    // Scan all cases (iL leading jets)
    for(UInt_t iL=0 ; iL<4 ; iL++) {
      TString siL = TString(Form("%d",iL));
      if(minDPhi_Jet_ETM[iL] >=0) hTH3F["hMinDPhi_Jet_ETM"+siL]  -> Fill( iJetCut, met[0], minDPhi_Jet_ETM[iL] ); 
      if(minDPhi_JetC_ETM[iL]>=0) hTH3F["hMinDPhi_JetC_ETM"+siL] -> Fill( iJetCut, met[0], minDPhi_JetC_ETM[iL] ); 
    }
    //}
    //if(mht[0]>=iETM) {
    if(maxDPhi_Jet_HTM>=0)  hTH3F["hMaxDPhi_Jet_HTM"]  -> Fill( iJetCut, mht[0], maxDPhi_Jet_HTM ); 
    if(maxDPhi_JetC_HTM>=0) hTH3F["hMaxDPhi_JetC_HTM"] -> Fill( iJetCut, mht[0], maxDPhi_JetC_HTM ); 

    for(UInt_t iL=0 ; iL<4 ; iL++) {
      TString siL = TString(Form("%d",iL));
      if(minDPhi_Jet_HTM[iL] >=0) hTH3F["hMinDPhi_Jet_HTM"+siL] -> Fill( iJetCut, mht[0], minDPhi_Jet_HTM[iL] ); 
      if(minDPhi_JetC_HTM[iL]>=0) hTH3F["hMinDPhi_JetC_HTM"+siL]-> Fill( iJetCut, mht[0], minDPhi_JetC_HTM[iL]); 
    }
    //}
    //}
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

void RateEfficiency::RescaleHistos3D(float scaleFactor)
{
  for(it_hTH3F=hTH3F.begin() ; it_hTH3F!=hTH3F.end() ; it_hTH3F++) {
    it_hTH3F->second->Scale(scaleFactor);
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

int RateEfficiency::InitHistos3D()
{
  const int nA=3;
  const int nH=8;
  TString name_histos[nH] = {"hMaxDPhi_Jet_ETM","hMaxDPhi_Jet_HTM","hMaxDPhi_JetC_ETM","hMaxDPhi_JetC_HTM",
			     "hMinDPhi_Jet_ETM","hMinDPhi_Jet_HTM","hMinDPhi_JetC_ETM","hMinDPhi_JetC_HTM"};

  TString title_histos[nH] = {"Max | #Delta#Phi(Jet, ETM) |", "Max | #Delta#Phi(Jet, HTM) |",
			      "Max | #Delta#Phi(JetC, ETM) |","Max | #Delta#Phi(JetC, HTM) |",
			      "Min | #Delta#Phi(Jet, ETM) |", "Min | #Delta#Phi(Jet, HTM) |",
			      "Min | #Delta#Phi(JetC, ETM) |","Min | #Delta#Phi(JetC, HTM) |"};


  TString xytitle[nH][nA] = { {"Jet eT cut [GeV]", "ETM cut", "Max | #Delta#Phi(Jet, ETM) |"},
			      {"Jet eT cut [GeV]", "HTM cut", "Max | #Delta#Phi(Jet, HTM) |"},
			      {"Central Jet eT cut [GeV]", "ETM cut", "Max | #Delta#Phi(JetC, ETM) |"},
			      {"Central Jet eT cut [GeV]", "HTM cut", "Max | #Delta#Phi(JetC, HTM) |"},
			      {"Jet eT cut [GeV]", "ETM cut", "Min | #Delta#Phi(Jet, ETM) |"},
			      {"Jet eT cut [GeV]", "HTM cut", "Min | #Delta#Phi(Jet, HTM) |"},
			      {"Central Jet eT cut [GeV]", "ETM cut", "Min | #Delta#Phi(JetC, ETM) |"},
			      {"Central Jet eT cut [GeV]", "HTM cut", "Min | #Delta#Phi(JetC, HTM) |"} };

  Int_t nbins[2]={MAXJET, MAXETM};
  Int_t xmax[2] ={MAXJET, MAXETM};
  Int_t xmin =0;

  TString locname="";
  TString siL="";

  for(uint iH=0 ; iH<nH ; iH++) {
    
    if(iH>=4) {
      for(UInt_t iL=0 ; iL<4 ; iL++) {
	siL = TString(Form("%d",iL));
	hTH3F[name_histos[iH]+siL] = new TH3F(name_histos[iH]+siL,
					      title_histos[iH]+"("+siL+" jets)", 
					      nbins[0],xmin,xmax[0],
					      nbins[1],xmin,xmax[1],
					      50,0,10);
	
	hTH3F[name_histos[iH]+siL]->SetXTitle(xytitle[iH][0]);
	hTH3F[name_histos[iH]+siL]->SetYTitle(xytitle[iH][1]);
	hTH3F[name_histos[iH]+siL]->SetZTitle(xytitle[iH][2]);
      }
    }

    else {
      hTH3F[name_histos[iH]] = new TH3F(name_histos[iH],title_histos[iH], 
					nbins[0],xmin,xmax[0],
					nbins[1],xmin,xmax[1],
					50,0,10);
      
      hTH3F[name_histos[iH]]->SetXTitle(xytitle[iH][0]);
      hTH3F[name_histos[iH]]->SetYTitle(xytitle[iH][1]);
      hTH3F[name_histos[iH]]->SetZTitle(xytitle[iH][2]);
    }
  }

  return 0;
}

int RateEfficiency::InitHistos2D()
{

  const int nA=2;
  const int nH=12;
  TString name_histos[nH] = {"hEff_Jet_ETM", "hEff_Jet_HTM", "hEff_JetC_ETM", "hEff_JetC_HTM",
			     "hEff_Jet_ETM_DPhi", "hEff_Jet_HTM_DPhi", "hEff_JetC_ETM_DPhi", "hEff_JetC_HTM_DPhi",
			     "hEff_Jet_ETM_OR70", "hEff_JetC_ETM_OR70", "hEff_Jet_ETM_DPhi_OR70", "hEff_JetC_ETM_DPhi_OR70"};

  TString title_histos[nH] = {"L1_Jet_ETM efficiency" ,     "L1_Jet_HTM efficiency", 
			      "L1_JetC_ETM efficiency",     "L1_JetC_HTM efficiency",
			      "L1_Jet_ETM_DPhi efficiency", "L1_Jet_HTM_DPhi efficiency",
			      "L1_JetC_ETM_DPhi efficiency","L1_JetC_HTM_DPhi efficiency",
			      "(L1_Jet_ETM || L1_ETM70) efficiency",
			      "(L1_JetC_ETM || L1_ETM70) efficiency",
			      "(L1_Jet_ETM_DPhi || L1_ETM70) efficiency",
			      "(L1_JetC_ETM_DPhi || L1_ETM70) efficiency"};

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
			      {"Central Jet eT cut [GeV]", "ETM cut [GeV]"} };

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

int RateEfficiency::OrderJets()
{

  nAllJets=0;
  allJetEt.clear();
  allJetEta.clear();
  allJetPhi.clear();
  allJetBx.clear();

  vector<TLorentzVector> vecJets;
  TLorentzVector theJet;

  // Fill the vector with central jets
  for(UInt_t iJ=0 ; iJ<nCenJets ; iJ++) {
    theJet.SetPtEtaPhiM(cenJetEt[iJ], cenJetEta[iJ], cenJetPhi[iJ], 0);
    vecJets.push_back(theJet);
  }

  // Fill the vector with forward jets
  for(UInt_t iJ=0 ; iJ<nFwdJets ; iJ++) {
    theJet.SetPtEtaPhiM(fwdJetEt[iJ], fwdJetEta[iJ], fwdJetPhi[iJ], 0);
    vecJets.push_back(theJet);
  }

  // Sort the jets by pT
  sort(vecJets.begin(), vecJets.end(), HighestPt());
  
  for(UInt_t iJ=0 ; iJ<vecJets.size() ; iJ++) {
    allJetEt.push_back(vecJets[iJ].Pt());
    allJetEta.push_back(vecJets[iJ].Eta());
    allJetPhi.push_back(vecJets[iJ].Phi());
    nAllJets++ ;
  }

  if(verbose>0 && nFwdJets!=0) {
    cout << "--- nAllJets=" << nAllJets << endl;
    for(UInt_t iJ=0 ; iJ<nAllJets ; iJ++) {
      cout << "allJetEt["  << iJ << "]=" << allJetEt[iJ]  << " ; "
	   << "allJetPhi[" << iJ << "]=" << allJetPhi[iJ] << " ; "
	   << endl;
    }
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
  met		= l1extra_-> met;		
  metPhi	= l1extra_-> metPhi;	
  metBx         = l1extra_-> metBx;	

  if(verbose>2) cout << "%%%%%%%%%%%%%%%" << endl
		     << nMet     << " " << (l1extra_->met)[0] << " "
		     << met[0]   << " " << metPhi[0] << " " 
		     << metBx[0] << endl 
		     << "%%%%%%%%%%%%%%%" << endl;

  et		= l1extra_-> et;		
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

float RateEfficiency::computeDeltaPhi(float phi1, float phi2)
{
  // Return value in [0;pi] 
  /*
  float dphi0 = TMath::Abs( phi1 - phi2 );
  if(dphi0 > TMath::Pi()) return  TMath::TwoPi() - dphi0;
  else                    return dphi0;
  */
  return TMath::Abs( TVector2::Phi_mpi_pi(phi1 - phi2) );

}

float RateEfficiency::computeDeltaR(float phi1, float phi2, float eta1, float eta2)
{
  // Return value in [0;pi] 
  float dphi0 = TMath::Abs( phi1 - phi2 );
  if(dphi0 > TMath::Pi()) dphi0 = TMath::TwoPi() - dphi0;

  float deta = TMath::Abs(eta1 - eta2);

  return TMath::Sqrt( dphi0*dphi0 + deta*deta );
}

int RateEfficiency::InitVar()
{

  nAllJets=0;
  allJetEt.clear();
  allJetEta.clear();
  allJetPhi.clear();
  allJetBx.clear();

  nIsoEm=0;
  isoEmEt.clear();
  isoEmEta.clear();
  isoEmPhi.clear();
  isoEmBx.clear();

  nNonIsoEm=0;
  nonIsoEmEt.clear();
  nonIsoEmEta.clear();
  nonIsoEmPhi.clear();
  nonIsoEmBx.clear();

  nCenJets=0;
  cenJetEt.clear();
  cenJetEta.clear();
  cenJetPhi.clear();
  cenJetBx.clear();

  nFwdJets=0;
  fwdJetEt.clear();
  fwdJetEta.clear();
  fwdJetPhi.clear();
  fwdJetBx.clear();

  nTauJets=0;
  tauJetEt.clear();
  tauJetEta.clear();
  tauJetPhi.clear();
  tauJetBx.clear();

  nMuons=0;
  muonEt.clear();
  muonEta.clear();
  muonPhi.clear();
  muonChg.clear();
  muonIso.clear();
  muonFwd.clear();
  muonMip.clear();
  muonRPC.clear();
  muonBx.clear();
  muonQuality.clear();

  hfEtSum.clear();
  hfBitCnt.clear();
  hfBx.clear();

  nMet=0;
  et.clear();
  met.clear();
  metPhi.clear();
  metBx.clear();

  nMht=0;
  ht.clear();
  mht.clear();
  mhtPhi.clear();
  mhtBx.clear();

  return 0;
}
