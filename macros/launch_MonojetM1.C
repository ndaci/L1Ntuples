{

  gROOT->ProcessLine(".x ../../L1Ntuples/macros/initL1Analysis.C+");
  gSystem->Load("RateEfficiency_cc.so");
  gSystem->Load("runRateEfficiency_C.so");

  runRateEfficiency("MonojetAV_M1_40PU_25bx", false, 0);

  return 0;
}
