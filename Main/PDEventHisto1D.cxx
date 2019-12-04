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

	    h_mom_gentruep = new TH1D("h_mom_gentruep", "h_mom_gentruep", 50, 0.0, 3.0);  		
	    h_mom_recotruep = new TH1D("h_mom_recotruep", "h_mom_recotruep", 50, 0.0, 3.0);  		

	    h_tmdqdx=new TH1D("h_tmdqdx", "h_tmdqdx", 100, 0.0, 1200);
	    h_tmdqdx_proton=new TH1D("h_tmdqdx_proton", "h_tmdqdx_proton", 100, 0.0, 1200);
	    h_tmdqdx_pionpm=new TH1D("h_tmdqdx_pionpm", "h_tmdqdx_pionpm", 100, 0.0, 1200);
	   

            h_tmdqdxvsrange_proton=new TH2D("h_tmdqdxvsrange_proton", "h_tmdqdxvsrange_proton", 100, 0.0, 1200, 100, 0.0, 300.0);
            h_tmdqdxvsrange_pionpm=new TH2D("h_tmdqdxvsrange_pionpm", "h_tmdqdxvsrange_pionpm", 100, 0.0, 1200, 100, 0.0, 300.0);
 
	    h_trackscore_proton=new TH1D("h_trackscore_proton", "h_trackscore_proton", 100, -0.0, 1.0);
	    h_trackscore_pionpm=new TH1D("h_trackscore_pionpm", "h_trackscore_pionpm", 100, -0.0, 1.0);
	    h_trackscore_electron=new TH1D("h_trackscore_electron", "h_trackscore_electron", 100, -0.0, 1.0);
	    h_trackscore_pion0=new TH1D("h_trackscore_pion0", "h_trackscore_pion0", 100, -0.0, 1.0);
	    h_trackscore_muon=new TH1D("h_trackscore_muon", "h_trackscore_muon", 100, -0.0, 1.0);
	    h_trackscore_other=new TH1D("h_trackscore_other", "h_trackscore_other", 100, -0.0, 1.0);
	    h_trackscore_photon=new TH1D("h_trackscore_photon", "h_trackscore_photon", 100, -0.0, 1.0);
            h_trktovtx = new TH1D("h_trktovtx", "h_trktovtx", 100, 0.0, 50.0);
            h_trktovtx_true = new TH1D("h_trktovtx_true", "h_trktovtx_true", 100, 0.0, 50.0);
								
	}
#endif
