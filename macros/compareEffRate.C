#include "/user/ndaci/WorkArea/Common/mymacros.h"

typedef map< float , vector< float > > MAPDPHI;
typedef map< float , MAPDPHI >         MAPJET;
typedef map< float , MAPJET  >         MAPETM;
typedef map< TString , MAPETM >        MAPALGO;

int doPalette(TPaletteAxis* palette);

int compareEffRate(TString path     = "results/", 
		   TString tagM1    = "MonojetAV_M1_40PU_25bx_RATE", 
		   TString tagM1000 = "MonojetAV_M1000_40PU_25bx_RATE", 
		   TString tagRate  = "crab_40PU_25bx_v2_RATE")
{
  ofstream outfile(path+"/scan/scan.txt", ios::out);

  // Get Inputs //
  MAPALGO mapalgo;

  const int nF=3;
  const int nH=20;

  TH3F *hMaxDPhi[nH][nF];
  TH2F *hTemp;
  TH1F *hTemp1;

  TString nameH[nH]={"hEff_ETM", "hEff_HTM",   "hEff_Jet", "hEff_JetC",
		     "hEff_Jet_ETM",           "hEff_JetC_ETM", 
		     "hEff_Jet_HTM",           "hEff_JetC_HTM",
		     "hEff_Jet_ETM_DPhi",      "hEff_JetC_ETM_DPhi", 
		     "hEff_Jet_HTM_DPhi",      "hEff_JetC_HTM_DPhi",
		     "hEff_Jet_ETM_OR70",      "hEff_JetC_ETM_OR70", 
		     "hEff_Jet_ETM_DPhi_OR70", "hEff_JetC_ETM_DPhi_OR70",
		     "hMaxDPhi_Jet_ETM",       "hMaxDPhi_Jet_HTM",
		     "hMaxDPhi_JetC_ETM",      "hMaxDPhi_JetC_HTM"};

  TString nameF[nF]={tagM1, tagM1000, tagRate};

  TString snameF[nF]={"MJ1", "MJ1000","NeuPU"};

  TFile *f[nF];
  TPaletteAxis *palette;
  TCanvas c("c","c",60,0,600,600);

  ofstream outlog(path+"/scan/values.txt", ios::out);
  
  for(int iF=0 ; iF<nF ; iF++) {
    f[iF] = new TFile(path+"/results_"+nameF[iF]+".root", "READ");

    if(f[iF]->IsZombie()) {
      cout << "ERROR OPENING FILE : " << nameF[iF] << endl;
      return 2;
    }

    outlog << endl << f[iF]->GetName() << endl << endl;

    for(int iH=0 ; iH<nH ; iH++) {

      if(iH<4) {
	f[iF]->GetObject(nameH[iH],hTemp1);
	if(!hTemp1) continue;

	hTemp1->Draw();
	gStyle->SetOptStat(0);

	c.Print(path+"/scan/"+nameH[iH]+"_"+snameF[iF]+".png");
	c.Print(path+"/scan/"+nameH[iH]+"_"+snameF[iF]+".pdf");
	continue;
      }

      else if(iH<16) {
	
	f[iF]->GetObject(nameH[iH],hTemp);
	if(!hTemp) continue;

	outlog << endl
	       << hTemp->GetName() << " (130,60) => " << hTemp->GetBinContent(131,61) << endl
	       << hTemp->GetName() << " (120,60) => " << hTemp->GetBinContent(121,61) << endl
	       << hTemp->GetName() << " (130,50) => " << hTemp->GetBinContent(131,51) << endl
	       << hTemp->GetName() << " (120,50) => " << hTemp->GetBinContent(121,51) << endl
	       << endl;

	if(iF==2) {
	  if(iH>=12 && iH<=15) {
	    hTemp->SetMinimum(5000);
	    hTemp->SetMaximum(13000);
	  }
	  else {
	    hTemp->SetMinimum(3000);
	    hTemp->SetMaximum(7000);
	  }
	}

	hTemp->Draw("colz");
	gStyle->SetOptStat(0);
	palette = (TPaletteAxis*) hTemp->GetListOfFunctions()->FindObject("palette");
	doPalette(palette);
	gPad->Update();

	c.Print(path+"/scan/"+nameH[iH]+"_"+snameF[iF]+".png");
	c.Print(path+"/scan/"+nameH[iH]+"_"+snameF[iF]+".pdf");
	
	continue;
      }
      
      f[iF]->GetObject(nameH[iH],hMaxDPhi[iH][iF]);
      if(!hMaxDPhi[iH][iF]) continue;

      // Scan ETM, Jet eT, DPhi cuts
      int ETM0=41;
      int JET0=101;
      int ETMD=5;
      int JETD=5;
      float jetcut=0;
      float etmcut=0;
      float dphicut=0;
      float integral=0;

      for(int iETM=0 ; iETM<7 ; iETM+=1) {

	// Project 3D=>2D at fixed ETM bin
	etmcut = ETM0+iETM*ETMD;
	hMaxDPhi[iH][iF]->GetYaxis()->SetRange(etmcut,etmcut);
	hTemp = (TH2F*) (hMaxDPhi[iH][iF]->Project3D("zx")) ;
	hTemp->GetYaxis()->SetRangeUser(0,4);
	
	// Plot the 2D projection
	//TCanvas c("c","c",60,0,600,600);
	hTemp->Draw("colz");
	gStyle->SetOptStat(0);
	palette = (TPaletteAxis*) hTemp->GetListOfFunctions()->FindObject("palette");
	doPalette(palette);
	gPad->Update();

	TString nameETM = Form("%d", int(etmcut-1) );
	c.Print(path+"/scan/"+nameH[iH]+"_"+snameF[iF]+"_ETM"+nameETM+".png");
	c.Print(path+"/scan/"+nameH[iH]+"_"+snameF[iF]+"_ETM"+nameETM+".pdf");

	// Scan jet cut and dphi cut values => integrate => efficiency/rate
	for(int iJetCut=0 ; iJetCut<21 ; iJetCut++) {
	  for(int iDPhi=0 ; iDPhi<16 ; iDPhi++) {
	    jetcut  = JET0+iJetCut*JETD;
	    dphicut = iDPhi*0.2;
	    //integral = hTemp->Integral(jetcut,jetcut,iDPhi+1,-1);
	    integral = hTemp->Integral(jetcut,jetcut,iDPhi+1,16);
	    mapalgo[nameH[iH]][etmcut-1][jetcut-1][dphicut].push_back(integral);
	  }
	}
      }
      ////////////////////////////////

    } // end loop over histos
    
  } // end loop over files

  // Print out values
  MAPETM  mapetm;
  MAPJET  mapjet;
  MAPDPHI mapdphi;
  vector<float> values;

  MAPALGO::iterator iterAlgo;
  MAPETM::iterator  iterEtm;
  MAPJET::iterator  iterJet;
  MAPDPHI::iterator iterDphi;

  for(int iH=0 ; iH<nH ; iH++) {
    mapetm = mapalgo[nameH[iH]];
    for(iterEtm=mapetm.begin() ; iterEtm!=mapetm.end() ; iterEtm++) {
      mapjet = iterEtm->second;
      for(iterJet=mapjet.begin() ; iterJet!=mapjet.end() ; iterJet++) {
	mapdphi = iterJet->second;
	for(iterDphi=mapdphi.begin() ; iterDphi!=mapdphi.end() ; iterDphi++) {
	  values = iterDphi->second;

	  outfile << nameH[iH]       << "   " 
		  << iterEtm ->first << "   "
		  << iterJet ->first << "   "
		  << iterDphi->first << "   " ;

	  for(uint iV=0 ; iV<values.size() ; iV++) {
	    outfile << values[iV] << "   " ;
	  }
	  if( (values[0]>=0.73 && values[1]>=0.78) && values[2]<=10000 )
	    outfile << "   WORTH" ;
	  outfile << endl;
	}
	outfile << endl;
      }
      outfile << endl;
    }
    outfile << endl;
  }

  return 0;
}

int doPalette(TPaletteAxis* palette)
{
  palette->SetLabelSize(0.03);
  palette->SetX1NDC(0.904362); 
  palette->SetY1NDC(0.101399); 
  palette->SetX2NDC(0.937919); 
  palette->SetY2NDC(0.90035);
  return 0;
}
