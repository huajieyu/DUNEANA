/**
 * \file Hanalysis.h
 *
 * \ingroup Main
 * 
 * \brief Class def header for a class Hanalysis
 *
 * @author jiangl
 */

/** \addtogroup Main

    @{*/
#ifndef __MAIN_HANALYSIS_H__
#define __MAIN_HANALYSIS_H__
#include "PDEventPro.h"
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

#include <TRandom1.h>
#include <TSystem.h>
#include <TBrowser.h>
#include <TApplication.h>
#include <TChain.h>
#include "TThread.h"
#include "THStack.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include <TString.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TCanvas.h>
#include "TMath.h"
#include <TVector3.h>
#include "PDEventDataHisto1D.h"


namespace Main {

  /**
     \class Hanalysis
     User defined class Hanalysis ... these comments are used to generate
     doxygen documentation!
  */
  class Hanalysis{
    
  public:
    
    /// Default constructor
    Hanalysis(){}
    
    /// Default destructor
    ~Hanalysis(){}

    void SetInputFile(std::string);
    void SetInputFile_add(std::string);
    void SetOutputFile(std::string);
    void SetInitialEntry(int);
    void SetEntries(int);
    void SetFVCut(bool); 
    void MakeFile();

  private:
    void DrawProgressBar(double progress, double barWidth);
    bool inFV(double x, double y, double z);
    double thetax(double theta, double phi);
    bool FVcuton = false; 

    std::string filen = "";
    std::string filen_add = "";
    std::string fileoutn = "";
    int _initial_entry = 0;
    int maxEntries = -1;

    PDEventDataHisto1D *_event_histo_1d;

    bool isProton = false;
    double Ecalcmiss(double Esum, double PTmiss, int np); 

    const double NeutronMass = 0.93956542; 
    const double ProtonMass = 0.938272;


    double cutAPA3_Z = 226.;
    double cut_trackScore = 0.35;
    //daughter Distance cut
    double cut_daughter_track_distance = 10.;
    double cut_daughter_shower_distance_low = 2.;
    double cut_daughter_shower_distance_high = 100.;
    double cut_primary_chi2 = 140.;
    double cut_secondary_chi2 = 50.;
    int cut_nHits_shower_low = 12;
    int cut_nHits_shower_high = 1000;
    //
    //For MC from Owen Goodwins studies
    double xlow = -3.,  xhigh = 7.,  ylow = -8.,  yhigh = 7.;
    double zlow = 27.5,  zhigh = 32.5,  coslow = 0.93;
    //
    //For Data from Owen Goodwin
    double data_xlow = 0., data_xhigh = 10., data_ylow= -5.;
    double data_yhigh= 10., data_zlow=30., data_zhigh=35., data_coslow=.93;

    bool data_beam_PID(const std::vector<int> *pidCandidates);
    bool manual_beamPos_data (int event,            double data_startX,
                              double data_startY,   double data_startZ,
                              double data_dirX,     double data_dirY,
                              double data_dirZ,     double data_BI_X,
                              double data_BI_Y,     double data_BI_dirX,
                              double data_BI_dirY,  double data_BI_dirZ,
                              int data_BI_nMomenta, int data_BI_nTracks); 
     bool endAPA3(double reco_beam_endZ); 
     bool has_shower_nHits_distance(const std::vector<double> &track_score,
                                    const std::vector<int> &nHits,
                                    const std::vector<double> &distance); 
 

  };
}

#endif
/** @} */ // end of doxygen group 

