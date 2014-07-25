#ifndef __L1Analysis_L1AnalysisGenMetDataFormat_H__
#define __L1Analysis_L1AnalysisGenMetDataFormat_H__

//-------------------------------------------------------------------------------
// Created 19/06/2014 - Nadir Daci
// 
//
// Addition of met gen information
//-------------------------------------------------------------------------------

#include <vector>

namespace L1Analysis
{
  struct L1AnalysisGenMetDataFormat
  {
    L1AnalysisGenMetDataFormat(){Reset();};
    ~L1AnalysisGenMetDataFormat(){Reset();};
    
    void Reset()
    {
     met    = -999;
     metPhi = -999;
     Ht     = -999;
     mHt    = -999;
     mHtPhi = -999;
     sumEt  = -999;
     neutralEM_Et_fraction  = -999;
     chargedEM_Et_fraction  = -999;
     neutralHad_Et_fraction = -999;
     chargedHad_Et_fraction = -999;
     muon_Et_fraction       = -999;
     inv_Et_fraction        = -999;
    }
    
    double met;
    double metPhi;
    double Ht;
    double mHt;
    double mHtPhi;
    double sumEt;

    double neutralEM_Et_fraction;
    double chargedEM_Et_fraction;
    double neutralHad_Et_fraction;
    double chargedHad_Et_fraction;
    double muon_Et_fraction;
    double inv_Et_fraction;
  }; 
}
#endif
