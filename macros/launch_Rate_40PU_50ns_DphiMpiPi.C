{

  gROOT->ProcessLine(".x ../../L1Ntuples/macros/initL1Analysis.C+");
  gSystem->Load("RateEfficiency_cc.so");
  gSystem->Load("runRateEfficiency_C.so");

  runRateEfficiency("crab_40PU_50bx_v2_PhiMpiPi", false, 0);

  return 0;
}
