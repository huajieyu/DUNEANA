#ifndef __MAIN_PDEVENTHISTO1D_CXX__
#define __MAIN_PDEVENTHISTO1D_CXX__

#include "PDEventHisto1D.h"
using namespace Main;
	void PDEventHisto1D::InitializeBootstraps()
	{   std::cout<<"initializing histograms"<<std::endl;
    	    h_orimom_proton = new TH1D("h_orimom_proton", "h_orimom_proton", 100, 0.0001, 0.8);
	    h_orimom_neutron = new TH1D("h_orimom_neutron", "h_orimom_neutron", 100, 0.0001, 0.8);
	    h_orimom_pionpm = new TH1D("h_orimom_pionpm", "h_orimom_pionpm", 100, 0.0001, 0.8);
	    h_orimom_pion0 = new TH1D("h_orimom_pion0", "h_orimom_pion0", 100, 0.0001, 0.8);
	    h_orimom_kaon = new TH1D("h_orimom_kaon", "h_orimom_kaon", 100, 0.0001, 0.8);
	    h_orimom_electron = new TH1D("h_orimom_electron", "h_orimom_electron", 100, 0.0001, 0.8);
	    h_orimom_muon = new TH1D("h_orimom_muon", "h_orimom_muon", 100, 0.0001, 0.8);
	    h_orimom_photon = new TH1D("h_orimom_photon", "h_orimom_photon", 100, 0.0001, 0.8);
	    h_orimom_other = new TH1D("h_orimom_other", "h_orimom_other", 100, 0.0001, 0.8);



	    h_chi2_phypo_proton= new TH1D("h_chi2_phypo_proton", "h_chi2_phypo_proton", 100, 0, 500);
	    h_chi2_phypo_pionpm= new TH1D("h_chi2_phypo_pionpm", "h_chi2_phypo_pionpm", 100, 0, 500);
	    h_chi2_phypo_muon= new TH1D("h_chi2_phypo_muon", "h_chi2_phypo_muon", 100, 0, 500);
	    h_chi2_phypo_electron= new TH1D("h_chi2_phypo_electron", "h_chi2_phypo_electron", 100, 0, 500);
	    h_chi2_phypo_kaon= new TH1D("h_chi2_phypo_kaon", "h_chi2_phypo_kaon", 100, 0, 500);
	    h_chi2_phypo_other= new TH1D("h_chi2_phypo_other", "h_chi2_phypo_other", 100, 0, 500);	    
	    h_chi2_phypo_mc= new TH1D("h_chi2_phypo_mc", "h_chi2_phypo_mc", 100, 0, 500);	    



	    h_recomom_selected_proton= new TH1D("h_recomom_selected_proton", "h_recomom_selected_proton", 100, 0, 1.5);
	    h_recomom_selected_pionpm= new TH1D("h_recomom_selected_pionpm", "h_recomom_selected_pionpm", 100, 0, 1.5);
	    h_recomom_selected_muon= new TH1D("h_recomom_selected_muon", "h_recomom_selected_muon", 100, 0, 1.5);
	    h_recomom_selected_electron= new TH1D("h_recomom_selected_electron", "h_recomom_selected_electron", 100, 0, 1.5);
	    h_recomom_selected_kaon= new TH1D("h_recomom_selected_kaon", "h_recomom_selected_kaon", 100, 0, 1.5);
	    h_recomom_selected_other= new TH1D("h_recomom_selected_other", "h_recomom_selected_other", 100, 0, 1.5);	    

	    h_recopimom_selected_proton= new TH1D("h_recopimom_selected_proton", "h_recopimom_selected_proton", 100, 0, 100);
	    h_recopimom_selected_pionpm= new TH1D("h_recopimom_selected_pionpm", "h_recopimom_selected_pionpm", 100, 0, 100);
	    h_recopimom_selected_muon= new TH1D("h_recopimom_selected_muon", "h_recopimom_selected_muon", 100, 0, 100);
	    h_recopimom_selected_electron= new TH1D("h_recopimom_selected_electron", "h_recopimom_selected_electron", 100, 0, 100);
	    h_recopimom_selected_kaon= new TH1D("h_recopimom_selected_kaon", "h_recopimom_selected_kaon", 100, 0, 100);
	    h_recopimom_selected_other= new TH1D("h_recopimom_selected_other", "h_recopimom_selected_other", 100, 0, 100);	    


	    h_mom_gentruep = new TH1D("h_mom_gentruep", "h_mom_gentruep", 50, 0.0,2.0);  		
	    h_thetax_gentruep = new TH1D("h_thetax_gentruep", "h_thetax_gentruep", 50, -3.14,3.14);  		


	    h_mom_recotruep = new TH1D("h_mom_recotruep", "h_mom_recotruep", 50, 0.0, 2.0);  		
	    h_mom_recotruep_test = new TH1D("h_mom_recotruep_test", "h_mom_recotruep_test", 50, 0.0, 2.0);  		
	    h_thetax_recotruep_test = new TH1D("h_thetax_recotruep_test", "h_thetax_recotruep_test", 50, -3.14,3.14);  		
            h_mom_nonrecop = new TH1D("h_mom_nonrecop", "h_mom_nonrecop", 50, 0.0, 2.0);       								

            h_mom_selectedtruep = new TH1D("h_mom_selectedtruep", "h_mom_selectedtruep", 50, 0.0, 2.0);
            h_mom_trkscoretruep = new TH1D("h_mom_trkscoretruep", "h_mom_trkscoretruep", 50, 0.0, 2.0);
            h_mom_tmdqdxtruep = new TH1D("h_mom_tmdqdxtruep", "h_mom_tmdqdxtruep", 50, 0.0, 2.0);
            h_mom_chi2truep = new TH1D("h_mom_chi2truep", "h_mom_chi2truep", 50, 0.0, 2.0);

            h_mom_selectedtruepipm = new TH1D("h_mom_selectedtruepipm", "h_mom_selectedtruepipm", 50, 0.0, 1.2);
            h_mom_trkscoretruepipm = new TH1D("h_mom_trkscoretruepipm", "h_mom_trkscoretruepipm", 50, 0.0, 1.2);
            h_mom_tmdqdxtruepipm = new TH1D("h_mom_tmdqdxtruepipm", "h_mom_tmdqdxtruepipm", 50, 0.0, 1.2);
            h_mom_chi2truepipm = new TH1D("h_mom_chi2truepipm", "h_mom_chi2truepipm", 50, 0.0, 1.2);





	    h_mom_gentruepionpm = new TH1D("h_mom_gentruepionpm", "h_mom_gentruepionpm", 50, 0.0,1.2);  		
	    h_mom_recotruepionpm = new TH1D("h_mom_recotruepionpm", "h_mom_recotruepionpm", 50, 0.0, 1.2);  		
	 

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

            h_nonrecop_momvsnhits = new TH2D("h_nonrecop_momvsnhits", "h_nonrecop_momvsnhits", 100, 0.0, 1.5, 100, 0.0, 100);
            h_nonrecop_momvsnhits -> GetXaxis()->SetTitle("True Momentum (proton)[GeV]");
            h_nonrecop_momvsnhits -> GetYaxis()->SetTitle("True number of hits(proton)");
            h_nonrecop_momvsnhits -> SetTitle("");

        }
	void PDEventHisto1D::InitializeHistograms()
        { 
            h_lownhitsp_thetax = new TH1D("h_lownhitsp_thetax", "h_lownhitsp_thetax", 100, -1.0, 1.);       								
            h_lownhitsp_thetay = new TH1D("h_lownhitsp_thetay", "h_lownhitsp_thetay", 100, -1.0, 1.);       								
            h_lownhitsp_thetaz = new TH1D("h_lownhitsp_thetaz", "h_lownhitsp_thetaz", 100, -1.0, 1.);       								

            h_lownhitsp_startX = new TH1D("h_lownhitsp_startX", "h_lownhitsp_startX", 100, -400.0, 400.);
            h_lownhitsp_startY = new TH1D("h_lownhitsp_startY", "h_lownhitsp_startY", 100, -.0, 700.); 
            h_lownhitsp_startZ = new TH1D("h_lownhitsp_startZ", "h_lownhitsp_startZ", 100, -.0, 800.);       							
            h_lownhitsp_endX = new TH1D("h_lownhitsp_endX", "h_lownhitsp_endX", 100, -400.0, 400.); 
            h_lownhitsp_endY = new TH1D("h_lownhitsp_endY", "h_lownhitsp_endY", 100, -.0, 700.);   
            h_lownhitsp_endZ = new TH1D("h_lownhitsp_endZ", "h_lownhitsp_endZ", 100, -.0, 800.);       								

            h_chi2cutp_thetax = new TH1D("h_chi2cutp_thetax", "h_chi2cutp_thetax", 100, -1.0, 1.);  
            h_chi2cutp_thetay = new TH1D("h_chi2cutp_thetay", "h_chi2cutp_thetay", 100, -1.0, 1.);     
            h_chi2cutp_thetaz = new TH1D("h_chi2cutp_thetaz", "h_chi2cutp_thetaz", 100, -1.0, 1.);       								
            

            h_mom_gentruepion0 = new TH1D("h_mom_gentruepion0", "h_mom_gentruepion0", 50, 0.0, 1.2);
            h_mom_recotruepion0 = new TH1D("h_mom_recotruepion0", "h_mom_recotruepion0", 50, 0.0, 1.2);
            h_mom_recotruenonpi0 = new TH1D("h_mom_recotruenonpi0", "h_mom_recotruenonpi0", 50, 0.0, 1.2);



	    h_mom_gentruepionpm = new TH1D("h_mom_gentruepionpm", "h_mom_gentruepionpm", 50, 0.0, 1.2);
            h_mom_recotruepionpm = new TH1D("h_mom_recotruepionpm", "h_mom_recotruepionpm", 50, 0.0, 1.2);

            h_PiAbs_sel_pmom = new TH1D ("h_PiAbs_sel_pmom", "h_PiAbs_sel_pmom", 50, 0.0, 1.2);
            h_PiAbs_sig_pmom = new TH1D ("h_PiAbs_sig_pmom", "h_PiAbs_sig_pmom", 50, 0.0, 1.2);
            h_PiAbs_chxbac_pmom = new TH1D ("h_PiAbs_chxbac_pmom", "h_PiAbs_chxbac_pmom", 50, 0.0, 1.2);
            h_PiAbs_reabac_pmom = new TH1D ("h_PiAbs_reabac_pmom", "h_PiAbs_reabac_pmom", 50, 0.0, 1.2);

            h_PiAbs_sel_pcostheta = new TH1D ("h_PiAbs_sel_pcostheta", "h_PiAbs_sel_pcostheta", 50, -1., 1.0);
            h_PiAbs_sig_pcostheta = new TH1D ("h_PiAbs_sig_pcostheta", "h_PiAbs_sig_pcostheta", 50, -1., 1.0);
            h_PiAbs_chxbac_pcostheta = new TH1D ("h_PiAbs_chxbac_pcostheta", "h_PiAbs_chxbac_pcostheta", 50, -1., 1.0);
            h_PiAbs_reabac_pcostheta = new TH1D ("h_PiAbs_reabac_pcostheta", "h_PiAbs_reabac_pcostheta", 50, -1., 1.0);

            h_PiAbs_sel_ptheta = new TH1D ("h_PiAbs_sel_ptheta", "h_PiAbs_sel_ptheta", 50, -0., 3.14);
            h_PiAbs_sig_ptheta = new TH1D ("h_PiAbs_sig_ptheta", "h_PiAbs_sig_ptheta", 50, -0., 3.14);
            h_PiAbs_chxbac_ptheta = new TH1D ("h_PiAbs_chxbac_ptheta", "h_PiAbs_chxbac_ptheta", 50, -0., 3.14);
            h_PiAbs_reabac_ptheta = new TH1D ("h_PiAbs_reabac_ptheta", "h_PiAbs_reabac_ptheta", 50, -0., 3.14);


            h_PiAbs_sel_pphi = new TH1D ("h_PiAbs_sel_pphi", "h_PiAbs_sel_pphi", 50, -3.14, 3.14);
            h_PiAbs_sig_pphi = new TH1D ("h_PiAbs_sig_pphi", "h_PiAbs_sig_pphi", 50, -3.14, 3.14);
            h_PiAbs_chxbac_pphi = new TH1D ("h_PiAbs_chxbac_pphi", "h_PiAbs_chxbac_pphi", 50, -3.14, 3.14);
            h_PiAbs_reabac_pphi = new TH1D ("h_PiAbs_reabac_pphi", "h_PiAbs_reabac_pphi", 50, -3.14, 3.14);

            h_PiAbs_sel_pcosthetax= new TH1D("h_PiAbs_sel_pcosthetax", "h_PiAbs_sel_pcosthetax", 20, -1., 1.);

  

            h_PiAbs_sel_pmult = new TH1D ("h_PiAbs_sel_pmult", "h_PiAbs_sel_pmult", 5, -0.5, 4.5);
            h_PiAbs_sig_pmult = new TH1D ("h_PiAbs_sig_pmult", "h_PiAbs_sig_pmult", 5, -0.5, 4.5);
            h_PiAbs_chxbac_pmult = new TH1D ("h_PiAbs_chxbac_pmult", "h_PiAbs_chxbac_pmult", 5, -0.5, 4.5);
            h_PiAbs_reabac_pmult = new TH1D ("h_PiAbs_reabac_pmult", "h_PiAbs_reabac_pmult", 5, -0.5, 4.5);

            h_PiAbs_sig_energeticproton_reco_mom = new TH1D("h_PiAbs_sig_energeticproton_reco_mom", "h_PiAbs_sig_energeticproton_reco_mom", 15, 0.0, 1.2);
            h_PiAbs_sig_energeticproton_reco_costheta = new TH1D("h_PiAbs_sig_energeticproton_reco_costheta", "h_PiAbs_sig_energeticproton_reco_costheta", 15, -1.0, 1.0);
            h_PiAbs_sig_energeticproton_reco_phi = new TH1D("h_PiAbs_sig_energeticproton_reco_phi", "h_PiAbs_sig_energeticproton_reco_phi", 15, -3.14, 3.14);

            h_PiAbs_sig_energeticproton_truevsreco_mom = new TH2D("h_PiAbs_sig_energeticproton_truevsreco_mom", "h_PiAbs_sig_energeticproton_truevsreco_mom", 30, 0.0, 1.2, 30, 0.0, 1.2);
            h_PiAbs_sig_energeticproton_truevsreco_costheta = new TH2D("h_PiAbs_sig_energeticproton_truevsreco_costheta", "h_PiAbs_sig_energeticproton_truevsreco_costheta", 30, -1.0, 1.0, 30, -1.0, 1.0);
            h_PiAbs_sig_energeticproton_truevsreco_phi = new TH2D("h_PiAbs_sig_energeticproton_truevsreco_phi", "h_PiAbs_sig_energeticproton_truevsreco_phi", 30, -3.14, 3.14, 30, -3.14, 3.14);

            h_PiAbs_sig_energeticproton_mom = new TH1D("h_PiAbs_sig_energeticproton_mom", "h_PiAbs_sig_energeticproton_mom", 30, 0.0, 1.2);
            h_PiAbs_sig_energeticproton_costheta = new TH1D("h_PiAbs_sig_energeticproton_costheta", "h_PiAbs_sig_energeticproton_costheta", 30, -1.0, 1.0);
            h_PiAbs_sig_energeticproton_phi = new TH1D("h_PiAbs_sig_energeticproton_phi", "h_PiAbs_sig_energeticproton_phi", 30, -3.14, 3.14);

            h_PiAbs_gen_sig_energeticproton_mom = new TH1D("h_PiAbs_gen_sig_energeticproton_mom", "h_PiAbs_gen_sig_energeticproton_mom", 30, 0.0, 1.2);
            h_PiAbs_gen_sig_energeticproton_costheta = new TH1D("h_PiAbs_gen_sig_energeticproton_costheta", "h_PiAbs_gen_sig_energeticproton_costheta", 30, -1.0, 1.0);
            h_PiAbs_gen_sig_energeticproton_phi = new TH1D("h_PiAbs_gen_sig_energeticproton_phi", "h_PiAbs_gen_sig_energeticproton_phi", 30, -3.14, 3.14);

            h_PiAbs_sel_energeticproton_mom = new TH1D("h_PiAbs_sel_energeticproton_mom", "h_PiAbs_sel_energeticproton_mom", 15, 0.0, 1.2);
            h_PiAbs_sel_energeticproton_costheta = new TH1D("h_PiAbs_sel_energeticproton_costheta", "h_PiAbs_sel_energeticproton_costheta", 30, -1.0, 1.0);
            h_PiAbs_sel_energeticproton_phi = new TH1D("h_PiAbs_sel_energeticproton_phi", "h_PiAbs_sel_energeticproton_phi", 30, -3.14, 3.14);
            h_PiAbs_sel_energeticproton_costhetax = new TH1D("h_PiAbs_sel_energeticproton_costhetax", "h_PiAbs_sel_energeticproton_costhetax", 30, -1.0, 1.0);

            h_PiAbs_chxbac_energeticproton_mom = new TH1D("h_PiAbs_chxbac_energeticproton_mom", "h_PiAbs_chxbac_energeticproton_mom", 15, 0.0, 1.2);
            h_PiAbs_chxbac_energeticproton_costheta = new TH1D("h_PiAbs_chxbac_energeticproton_costheta", "h_PiAbs_chxbac_energeticproton_costheta", 15, -1.0, 1.0);
            h_PiAbs_chxbac_energeticproton_phi = new TH1D("h_PiAbs_chxbac_energeticproton_phi", "h_PiAbs_chxbac_energeticproton_phi", 15, -3.14, 3.14);

            h_PiAbs_reabac_energeticproton_mom = new TH1D("h_PiAbs_reabac_energeticproton_mom", "h_PiAbs_reabac_energeticproton_mom", 15, 0.0, 1.2);
            h_PiAbs_reabac_energeticproton_costheta = new TH1D("h_PiAbs_reabac_energeticproton_costheta", "h_PiAbs_reabac_energeticproton_costheta", 15, -1.0, 1.0);
            h_PiAbs_reabac_energeticproton_phi = new TH1D("h_PiAbs_reabac_energeticproton_phi", "h_PiAbs_reabac_energeticproton_phi", 15, -3.14, 3.14);




            h_PiAbs_sig_nmult = new TH1D ("h_PiAbs_sig_nmult", "h_PiAbs_sig_nmult", 15, -0.5, 14.5);
            h_PiAbs_sig_nmult ->GetXaxis()->SetTitle("Neutron Multiplicity[True]");
            h_PiAbs_sig_nmult ->GetYaxis()->SetTitle("No. of Events");

            h_PiAbs_sig_truevsreco_pmult = new TH2D ("h_PiAbs_sig_truevsreco_pmult", "h_PiAbs_sig_truevsreco_pmult", 10, -0.5, 9.5, 10, -0.5, 9.5);
            h_PiAbs_sig_truevsreco_pmult ->GetXaxis()->SetTitle("Proton Multiplicity[True]");
            h_PiAbs_sig_truevsreco_pmult ->GetYaxis()->SetTitle("Proton_Multiplicity[Reco]");

            h_sig_pvsnmult = new TH2D("h_sig_pvsnmult", "h_sig_pvsnmult", 15, -0.5, 14.5, 15, -0.5, 14.5);
            h_sig_pvsnmult->GetXaxis()->SetTitle("No. of Protons");
            h_sig_pvsnmult->GetYaxis()->SetTitle("No. of Neutrons");

            h_PiAbs_sig_true_totalKE_protonvsneutron = new TH2D("h_PiAbs_sig_true_totalKE_protonvsneutron", "h_PiAbs_sig_true_totalKE_protonvsneutron", 50, 0., 1., 50, 0., 1.);
            h_PiAbs_sig_true_totalKE_protonvsneutron->GetXaxis()->SetTitle("Total KE:Protons[GeV]");
            h_PiAbs_sig_true_totalKE_protonvsneutron->GetYaxis()->SetTitle("Total KE:Neutrons[GeV]");


            h_PiAbs_sig_proton_mom = new TH1D("h_PiAbs_sig_proton_mom", "h_PiAbs_sig_proton_mom", 30, 0.0, 1.0);

            h_PiAbs_sig_proton_mom->GetXaxis()->SetTitle("True Momentum of Proton [GeV]");
            h_PiAbs_sig_proton_mom->GetYaxis()->SetTitle("No. of Tracks");

          

            h_PiAbs_sig_neutron_mom = new TH1D("h_PiAbs_sig_neutron_mom", "h_PiAbs_sig_neutron_mom", 30, 0.0, 1.0);

            h_PiAbs_sig_neutron_mom->GetXaxis()->SetTitle("True Momentum of Neutron [GeV]");
            h_PiAbs_sig_neutron_mom->GetYaxis()->SetTitle("No. of Tracks");

            h_PiAbs_sig_proton_costheta = new TH1D("h_PiAbs_sig_proton_costheta", "h_PiAbs_sig_proton_costheta", 30, -1., 1.0);

            h_PiAbs_sig_proton_costheta->GetXaxis()->SetTitle("True cos#theta of Proton [GeV]");
            h_PiAbs_sig_proton_costheta->GetYaxis()->SetTitle("No. of Tracks");

          

            h_PiAbs_sig_neutron_costheta = new TH1D("h_PiAbs_sig_neutron_costheta", "h_PiAbs_sig_neutron_costheta", 30, -1., 1.0);

            h_PiAbs_sig_neutron_costheta->GetXaxis()->SetTitle("True cos#theta of Neutron [GeV]");
            h_PiAbs_sig_neutron_costheta->GetYaxis()->SetTitle("No. of Tracks");




            h_PiAbs_sig_true_Ptmissing = new TH1D("h_PiAbs_sig_true_Ptmissing","h_PiAbs_sig_true_Ptmissing", 30, 0.0, 0.5);
            h_PiAbs_sig_true_Ptmissing->GetXaxis()->SetTitle("True missing-P_{T}[GeV]");
            h_PiAbs_sig_true_Ptmissing->GetYaxis()->SetTitle("No. of Tracks");

            h_PiAbs_sig_true_Ptmissing_onlyproton = new TH1D("h_PiAbs_sig_true_Ptmissing_onlyproton","h_PiAbs_sig_true_Ptmissing_onlyproton", 30, 0.0, 1.0);
            h_PiAbs_sig_true_Ptmissing_onlyproton->GetXaxis()->SetTitle("True missing_onlyproton-P_{T}[GeV]");
            h_PiAbs_sig_true_Ptmissing_onlyproton->GetYaxis()->SetTitle("No. of Tracks");


            h_sel_Ptmissing = new TH1D ("h_sel_Ptmissing", "h_sel_Ptmissing", 15, 0.0, 1.2);
            h_sel_Pmissing = new TH1D ("h_sel_Pmissing", "h_sel_Pmissing", 15, 0.0, 1.2);
            h_sel_Emissing = new TH1D ("h_sel_Emissing", "h_sel_Emissing", 15, 0.0, 1.2);
            h_sel_Plongit = new TH1D ("h_sel_Plongit", "h_sel_Plongit", 15, 0.0, 1.2);
 

            h_sig_Ptmissing = new TH1D ("h_sig_Ptmissing", "h_sig_Ptmissing", 15, 0.0, 1.2);
            h_sig_Pmissing = new TH1D ("h_sig_Pmissing", "h_sig_Pmissing", 15, 0.0, 1.2);
            h_sig_Emissing = new TH1D ("h_sig_Emissing", "h_sig_Emissing", 15, 0.0, 1.2);
            h_sig_Plongit = new TH1D ("h_sig_Plongit", "h_sig_Plongit", 15, 0.0, 1.2);
            
            h_chxbac_Ptmissing = new TH1D ("h_chxbac_Ptmissing", "h_chxbac_Ptmissing", 15, 0.0, 1.2);
            h_chxbac_Pmissing = new TH1D ("h_chxbac_Pmissing", "h_chxbac_Pmissing", 15, 0.0, 1.2);
            h_chxbac_Emissing = new TH1D ("h_chxbac_Emissing", "h_chxbac_Emissing", 15, 0.0, 1.2);
            h_chxbac_Plongit = new TH1D ("h_chxbac_Plongit", "h_chxbac_Plongit", 15, 0.0, 1.2);
            
            h_reabac_Ptmissing = new TH1D ("h_reabac_Ptmissing", "h_reabac_Ptmissing", 15, 0.0, 1.2);
            h_reabac_Pmissing = new TH1D ("h_reabac_Pmissing", "h_reabac_Pmissing", 15, 0.0, 1.2);
            h_reabac_Emissing = new TH1D ("h_reabac_Emissing", "h_reabac_Emissing", 15, 0.0, 1.2);
            h_reabac_Plongit = new TH1D ("h_reabac_Plongit", "h_reabac_Plongit", 15, 0.0, 1.2);
            
            h_reco_pionpm_rea = new TH1D("h_reco_pionpm_rea", "h_reco_pionpm_rea ", 10, -0.5, 9.5);
            h_reco_pionpm_rea -> GetXaxis() ->SetTitle("No. of #pi^{+}/#pi^{-}");
            h_reco_pionpm_rea -> GetYaxis() ->SetTitle("No. of Events");

            h_reco_photon_rea = new TH1D("h_reco_photon_rea", "h_reco_photon_rea ", 10, -0.5, 9.5);
            h_reco_photon_rea -> GetXaxis() ->SetTitle("No. of #pi^{0}");
            h_reco_photon_rea -> GetYaxis() ->SetTitle("No. of Events");
 
            h_grand_daughter_beam_dist = new TH1D("h_grand_daughter_beam_dist", "h_grand_daughter_beam_dist", 30, 0.0, 30.0); 
            h_daughter_beam_dist = new TH1D("h_daughter_beam_dist", "h_daughter_beam_dist", 30, 0.0, 30.0); 

            h_reco_grand_daughter_beam_dist = new TH1D("h_reco_grand_daughter_beam_dist", "h_reco_grand_daughter_beam_dist", 30, 0.0, 30.0); 
            h_reco_daughter_beam_dist = new TH1D("h_reco_daughter_beam_dist", "h_reco_daughter_beam_dist", 30, 0.0, 30.0); 

            h_reco_grand_daughter_angle = new TH1D("h_reco_grand_daughter_angle", "h_reco_grand_daughter_angle", 30, 0.0, 3.14); 
            h_reco_daughter_angle = new TH1D("h_reco_daughter_angle", "h_reco_daughter_angle", 30, 0.0, 3.14); 

            h_reco_daughter_deltax = new TH1D("h_reco_daughter_deltax", "h_reco_daughter_deltax", 30, -20., 40); 
            h_reco_daughter_deltay = new TH1D("h_reco_daughter_deltay", "h_reco_daughter_deltay", 30, -20., 40); 
            h_reco_daughter_deltaz = new TH1D("h_reco_daughter_deltaz", "h_reco_daughter_deltaz", 30, -20., 40); 






            h_reco_daughter_distvsangle=new TH2D("h_reco_daughter_distvsangle", "h_reco_daughter_distvsangle", 30, 0.0, 30.0, 30, 0.0, 3.14);


            h_gdfromproton= new TH1D("h_gdfromproton", "h_gdfromproton", 30, 0.0, 0.8);
            h_gdfromother= new TH1D("h_gdfromother", "h_gdfromother", 30, 0.0, 0.8);

            h_selproton_momreso = new TH2D("h_selproton_momreso", "h_selproton_momreso", 50, 0.0, 1.0, 50, 0.0, 1.0);
            h_selproton_costhetareso = new TH2D("h_selproton_costhetareso", "h_selproton_costhetareso", 50, -1.0, 1.0, 50, -1.0, 1.0);
            h_ngamma_frompi0=new TH1D("h_ngamma_frompi0", "h_ngamma_frompi0", 2.0, 0.5, 2.5);
            h_pgamma_frompi0=new TH1D("h_pgamma_frompi0", "h_pgamma_frompi0", 20, 0.0, 1.5);

            //selected_signal_percut = new TH1D("selected_signal_percut", "selected_percut", 4, 0, 4);
            //selected_percut = new TH1D("selected_percut", "selected_percut", 4, 0, 4);




	}
#endif
