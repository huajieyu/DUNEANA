/**
 * \file PDEventDataHisto1D.h
 *
 * \ingroup Main
 * 
 * \brief Class def header for a class PDEventDataHisto1D
 *
 * @author jiangl
 */

/** \addtogroup Main

    @{*/
#ifndef __MAIN_PDEVENTDATAHISTO1D_H__
#define __MAIN_PDEVENTDATAHISTO1D_H__
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
     \class PDEventDataHisto1D
     User defined class PDEventDataHisto1D ... these comments are used to generate
     doxygen documentation!
  */
  class PDEventDataHisto1D{
    
  public:
    
    /// Default constructor
    PDEventDataHisto1D(){}
    
    /// Default destructor
    ~PDEventDataHisto1D(){}
    void InitializeBootstraps();
    void InitializeHistograms();

    TH1D *h_trkstartx;
    TH1D *h_trkstarty;
    TH1D *h_trkstartz;

    TH1D *h_trkendx;
    TH1D *h_trkendy;
    TH1D *h_trkendz;

    TH1D *h_beamstartx;
    TH1D *h_beamstarty;
    TH1D *h_beamstartz;

    TH1D *h_beamendx;
    TH1D *h_beamendy;
    TH1D *h_beamendz;

    TH1D *h_trunmean_dqdx;
    TH1D *h_chi2_proton;
   
    TH1D *h_trk_trackscore; 

    TH1D *h_sel_pcandmom;
    TH1D *h_sel_pcandcostheta;
    TH1D *h_sel_pcandtheta;
    TH1D *h_sel_pcandphi;

    TH1D  *h_PiAbs_sel_energeticproton_mom;
    TH1D  *h_PiAbs_sel_energeticproton_costheta;
    TH1D  *h_PiAbs_sel_energeticproton_phi;
    TH1D  *h_PiAbs_sel_energeticproton_costhetax;
    
    TH1D  *h_beamend_trkstart_dist;
    TH1D  *h_beamend_trkstart_xdist;
    TH1D  *h_beamend_trkstart_ydist;
    TH1D  *h_beamend_trkstart_zdist;
 
    TH1D  *h_beamend_trkstart_angle;

    TH2D  *h_costheta_vs_phi;


    TH1D  *h_reco_daughter_deltax;
    TH1D  *h_reco_daughter_deltay;
    TH1D  *h_reco_daughter_deltaz;

    TH2D  *h_reco_daughter_distvsangle; 
    TH1D  *h_sel_pcosthetax;

    TH2D  *h_tof_momz;


    TH1D *h_chi2_phypo_data;

    TH1D  *h_PiAbs_sel_Ptmissing; 
    TH1D  *h_PiAbs_sel_Plongit;
    TH1D  *h_PiAbs_sel_pmult;
  };
}

#endif
/** @} */ // end of doxygen group 

