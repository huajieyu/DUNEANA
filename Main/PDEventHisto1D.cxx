#ifndef __MAIN_PDEVENTHISTO1D_CXX__
#define __MAIN_PDEVENTHISTO1D_CXX__

#include "PDEventHisto1D.h"
using namespace Main;
	void PDEventHisto1D::InitializeBootstraps()
	{   std::cout<<"initializing histograms"<<std::endl;
	    h_chi2_phypo_proton= new TH1D("h_chi2_phypo_proton", "h_chi2_phypo_proton", 100, 0, 500);
	    h_chi2_phypo_pionpm= new TH1D("h_chi2_phypo_pionpm", "h_chi2_phypo_pionpm", 100, 0, 500);
	    h_chi2_phypo_muon= new TH1D("h_chi2_phypo_muon", "h_chi2_phypo_muon", 100, 0, 500);
	    h_chi2_phypo_electron= new TH1D("h_chi2_phypo_electron", "h_chi2_phypo_electron", 100, 0, 500);
	    h_chi2_phypo_kaon= new TH1D("h_chi2_phypo_kaon", "h_chi2_phypo_kaon", 100, 0, 500);
	    h_chi2_phypo_other= new TH1D("h_chi2_phypo_other", "h_chi2_phypo_other", 100, 0, 500);	    
	    h_mom_truep = new TH1D("h_mom_truep", "h_mom_truep", 50, 0.0, 3.0);  		
	    h_tmdqdx=new TH1D("h_tmdqdx", "h_tmdqdx", 100, 0.0, 1200);
	    h_tmdqdx_proton=new TH1D("h_tmdqdx_proton", "h_tmdqdx_proton", 100, 0.0, 1200);
	    h_tmdqdx_pionpm=new TH1D("h_tmdqdx_pionpm", "h_tmdqdx_pionpm", 100, 0.0, 1200);
	}
#endif
