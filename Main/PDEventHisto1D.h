/**
 * \file PDEventHisto1D.h
 *
 * \ingroup DataTypes
 * 
 * \brief Class def header for a class PDEventHisto1D
 *
 * @author jiangl
 */

/** \addtogroup DataTypes

    @{*/
#ifndef __MAIN_PDEVENTHISTO1D_H__
#define __MAIN_PDEVENTHISTO1D_H__
#include <iostream>
#include <sstream>
#include <string>
#include <ctime>
#include <vector>
#include <map>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <fstream>
#include <math.h>
#include <algorithm>

#include <TChain.h>
#include "THStack.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TFile.h>
namespace Main {

  /**
     \class PDEventHisto1D
     User defined class PDEventHisto1D ... these comments are used to generate
     doxygen documentation!
  */
  class PDEventHisto1D{
    
  public:
    
    /// Default constructor
    PDEventHisto1D(){}
    
    /// Default destructor
    ~PDEventHisto1D(){}
    ///  
    void InitializeBootstraps();
    TH1D *h_chi2_phypo_proton;
    TH1D *h_chi2_phypo_pionpm;
    TH1D *h_chi2_phypo_muon;
    TH1D *h_chi2_phypo_electron;
    TH1D *h_chi2_phypo_kaon;
    TH1D *h_chi2_phypo_other;

    TH1D *h_mom_gentruep;
    TH1D *h_mom_recotruep;


    TH1D *h_tmdqdx;
    TH1D *h_tmdqdx_proton;
    TH1D *h_tmdqdx_pionpm;

    TH2D *h_tmdqdxvsrange_proton;
    TH2D *h_tmdqdxvsrange_pionpm;


    TH1D  *h_trackscore_proton;
    TH1D  *h_trackscore_electron;
    TH1D  *h_trackscore_pionpm;
    TH1D  *h_trackscore_photon;
    TH1D  *h_trackscore_muon;
    TH1D  *h_trackscore_pion0;
    TH1D  *h_trackscore_other;

    TH1D  *h_trktovtx;
    TH1D  *h_trktovtx_true;

  protected:

  };
}

#endif
/** @} */ // end of doxygen group 

