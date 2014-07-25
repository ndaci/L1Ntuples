#ifndef __L1Analysis_L1AnalysisGenJetDataFormat_H__
#define __L1Analysis_L1AnalysisGenJetDataFormat_H__

//-------------------------------------------------------------------------------
// Created 18/06/2014 - Nadir Daci
//
//
// Copied from L1AnalysisRecoJetDataFormat.h
//-------------------------------------------------------------------------------

#include <vector>

namespace L1Analysis
{
  struct L1AnalysisGenJetDataFormat
  {
    L1AnalysisGenJetDataFormat(){Reset();};
    ~L1AnalysisGenJetDataFormat(){Reset();};

    void Reset()
    {
    nJets=0;

    e.clear();
    et.clear();
    pt.clear();
    eta.clear();
    phi.clear();

    eEM.clear();
    eHad.clear();
    eInv.clear();
    eAux.clear();
    //etaDet.clear();

    }

    unsigned nJets;
    std::vector<double> e;
    std::vector<double> et;
    std::vector<double> pt;
    std::vector<double> eta;
    std::vector<double> phi;

    std::vector<double> eEM;
    std::vector<double> eHad;
    std::vector<double> eInv;
    std::vector<double> eAux;
    //std::vector<double> etaDet;

  };
}
#endif


