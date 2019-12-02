/**
 * \file Maker.h
 *
 * \ingroup Main
 * 
 * \brief Class def header for a class Maker
 *
 * @author jiangl
 */

/** \addtogroup Main

    @{*/
#ifndef __MAIN_MAKER_H__
#define __MAIN_MAKER_H__
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

#include "PDEventHisto1D.h"

namespace Main {

  /**
     \class Maker
     User defined class Maker ... these comments are used to generate
     doxygen documentation!
  */
  class Maker{
    
  public:
    
    /// Default constructor
    Maker(){}
    
    /// Default destructor
    ~Maker(){}

    void SetInputFile(std::string);
    void SetInputFile_add(std::string);
    void SetOutputFile(std::string);
    void SetInitialEntry(int);
    void SetEntries(int);
    
    void MakeFile();

  private:
    void DrawProgressBar(double progress, double barWidth);
    std::string filen = "";
    std::string filen_add = "";
    std::string fileoutn = "";
    int _initial_entry = 0;
    int maxEntries = -1;

    PDEventHisto1D *_event_histo_1d;

    
  };
}

#endif
/** @} */ // end of doxygen group 

