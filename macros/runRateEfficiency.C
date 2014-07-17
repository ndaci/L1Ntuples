#include "RateEfficiency.hh"

int runRateEfficiency(string fileType, int isCrossSec = false, int nEvents = 0) 
{
  int ctrl=-1;
  
  int nBunches50ns = 1368;
  int nBunches25ns = 2508; //2508 is what agreed with TSG for # bunches

  float xSec13TeV = isCrossSec ? 78.26 : 80.; // Using McM for cross section comparison and 80 (agreed with TSG) for rates
  float xSec8TeV  = 72.7; 

  // Launch the program on the relevant inputs //

  if (fileType == "crab_40PU_25bx_v2")
    {
      RateEfficiency myRateEfficiency("dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ndaci/Neutrino_Pt-2to20_gun/Rate_13TeV_40PU_25ns_62X_ReEmul2015_v2_10June2014/140610_144741/0000/L1Tree_14.root"); 
      ctrl = myRateEfficiency.run(false,"crab_40PU_25bx_v2",0,500000000,xSec13TeV,40,nBunches25ns,isCrossSec,nEvents,true);
    }

  else if (fileType == "crab_40PU_25bx_v2_list")
    {
      RateEfficiency myRateEfficiency;
      myRateEfficiency.OpenWithList("list_uct_25ns_v2.txt");
      ctrl = myRateEfficiency.run(false,"crab_40PU_25bx_v2",0,500000000,xSec13TeV,40,nBunches25ns,isCrossSec,nEvents,true);
    }

  else if (fileType == "crab_40PU_50bx_v2")
    {
      RateEfficiency myRateEfficiency("dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ndaci/Neutrino_Pt-2to20_gun/Rate_13TeV_40PU_25ns_62X_ReEmul2015_v2_10June2014/140610_153319/0000/L1Tree_79.root"); 
      ctrl = myRateEfficiency.run(false,"crab_40PU_50bx_v2",0,500000000,xSec13TeV,40,nBunches50ns,isCrossSec,nEvents,true);
    }

  else if (fileType == "crab_40PU_25bx_v4_Emul")
    {
      RateEfficiency myRateEfficiency("dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ndaci/Neutrino_Pt-2to20_gun/Rate_13TeV_40PU_25ns_62X_ReEmul2015_v4_12June2014_Emul/140611_161547/0000/L1Tree_38.root"); 
      ctrl = myRateEfficiency.run(false,fileType,0,500000000,xSec13TeV,40,nBunches25ns,isCrossSec,nEvents,true);
    }

  else if (fileType == "MonojetAV_M1_40PU_25bx")
    {
      RateEfficiency myRateEfficiency;
      myRateEfficiency.OpenWithList("list_MonojetM1.txt");
      ctrl = myRateEfficiency.run(false,fileType,0,500000000,xSec13TeV,40,nBunches25ns,isCrossSec,nEvents,false);
    }

  else if (fileType == "MonojetAV_M1000_40PU_25bx")
    {
      RateEfficiency myRateEfficiency;
      myRateEfficiency.OpenWithList("list_MonojetM1000.txt");
      ctrl = myRateEfficiency.run(false,fileType,0,500000000,xSec13TeV,40,nBunches25ns,isCrossSec,nEvents,false);
    }

  else if(fileType == "8TeV") {
      RateEfficiency myRateEfficiency("");
      ctrl = myRateEfficiency.run(false,fileType,0,500000000,xSec8TeV,40,nBunches25ns,isCrossSec,nEvents,true);
  }

  else 
    {
      cout << "Config param " << fileType << " invalid! \n"
		<< "Valid fileType values are : DATA, 8TEV_TF_DATA, 8TEV_TF_2012_RE-EMUL, "
		<< "8TEV_25PU_ORIG_RE-EMUL, 8TEV_25PU_2012_RE-EMUL, 8TEV_25PU_2012GCT10GEV_RE-EMUL, 8TEV_25PU_2015_RE-EMUL, "
		<< "13TEV_25PU_ORIG_RE-EMUL, 13TEV_25PU_2012_RE-EMUL, 13TEV_25PU_2012GCT10GEV_RE-EMUL, 13TEV_25PU_2015_RE-EMUL\n";
    }

  return ctrl;
}
