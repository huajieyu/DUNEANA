#ifndef __MAIN_PDEVENTDATAHISTO1D_CXX__
#define __MAIN_PDEVENTDATAHISTO1D_CXX__

#include "PDEventDataHisto1D.h"
using namespace Main;

	void PDEventDataHisto1D::InitializeBootstraps()
	{ 
              h_trkstartx = new TH1D("h_trkstartx", "h_trkstartx", 50,-350, 350);   
              h_trkstarty = new TH1D("h_trkstarty", "h_trkstarty", 50, 0.0, 700);   
              h_trkstartz = new TH1D("h_trkstartz", "h_trkstartz", 50, 0.0, 600);   
 
              h_trkendx = new TH1D("h_trkendx", "h_trkendx", 50,-350, 350);   
              h_trkendy = new TH1D("h_trkendy", "h_trkendy", 50, 0.0, 700);   
              h_trkendz = new TH1D("h_trkendz", "h_trkendz", 50, 0.0, 600);   
 
              h_beamstartx = new TH1D("h_beamstartx", "h_beamstartx", 50,-350, 350);   
              h_beamstarty = new TH1D("h_beamstarty", "h_beamstarty", 50, 0.0, 700);   
              h_beamstartz = new TH1D("h_beamstartz", "h_beamstartz", 50, 0.0, 600);   
 
              h_beamendx = new TH1D("h_beamendx", "h_beamendx", 50,-350, 350);   
              h_beamendy = new TH1D("h_beamendy", "h_beamendy", 50, 0.0, 700);   
              h_beamendz = new TH1D("h_beamendz", "h_beamendz", 50, 0.0, 600);   

              h_trunmean_dqdx = new TH1D("h_trunmean_dqdx", "h_trunmean_dqdx", 100, 0.0, 1200);
              h_trk_trackscore = new TH1D("h_trk_trackscore", "h_trk_trackscore", 100, 0.0, 1.0);

         } 
	void PDEventDataHisto1D::InitializeHistograms()
	{ 
              h_sel_pcandmom = new TH1D("h_sel_pcandmom", "h_sel_pcandmom",50, 0.0, 1.2 );
              h_sel_pcandcostheta = new TH1D("h_sel_pcandcostheta", "h_sel_pcandcostheta",50, -1.0, 1.0 );
              h_sel_pcandtheta = new TH1D("h_sel_pcandtheta", "h_sel_pcandtheta", 50, -0., 3.14);
              h_sel_pcandphi = new TH1D("h_sel_pcandphi", "h_sel_pcandphi",50, -3.14, 3.14);

              h_PiAbs_sel_energeticproton_mom = new TH1D("h_PiAbs_sel_energeticproton_mom", "h_PiAbs_sel_energeticproton_mom", 15, 0.0, 1.2);
              h_PiAbs_sel_energeticproton_costheta = new TH1D("h_PiAbs_sel_energeticproton_costheta", "h_PiAbs_sel_energeticproton_costheta", 15, -1.0, 1.0);
              h_PiAbs_sel_energeticproton_phi = new TH1D("h_PiAbs_sel_energeticproton_phi", "h_PiAbs_sel_energeticproton_phi", 15, -3.14, 3.14);
              h_PiAbs_sel_energeticproton_costhetax = new TH1D("h_PiAbs_sel_energeticproton_costhetax", "h_PiAbs_sel_energeticproton_costhetax", 10, -1.0, 1.0);

              h_beamend_trkstart_dist = new TH1D("h_beamend_trkstart_dist","h_beamend_trkstart_dist", 30, 0.0, 30.0);
              h_beamend_trkstart_xdist = new TH1D("h_beamend_trkstart_xdist","h_beamend_trkstart_xdist", 30, -15., 15.0);
              h_beamend_trkstart_ydist = new TH1D("h_beamend_trkstart_ydist","h_beamend_trkstart_ydist", 30, -15., 15.0);
              h_beamend_trkstart_zdist = new TH1D("h_beamend_trkstart_zdist","h_beamend_trkstart_zdist", 30, -15, 15.0);

              h_beamend_trkstart_angle = new TH1D("h_beamend_trkstart_angle","h_beamend_trkstart_angle", 30, -1.0, 1.0);

              h_costheta_vs_phi= new TH2D("h_costheta_vs_phi", "h_costheta_vs_phi", 50, -1.0, 1.0, 50, 0.0, 3.14); 

              h_reco_daughter_deltax = new TH1D("h_reco_daughter_deltax", "h_reco_daughter_deltax", 30, -20., 40); 
              h_reco_daughter_deltay = new TH1D("h_reco_daughter_deltay", "h_reco_daughter_deltay", 30, -20., 40); 
              h_reco_daughter_deltaz = new TH1D("h_reco_daughter_deltaz", "h_reco_daughter_deltaz", 30, -20., 40); 

              h_reco_daughter_distvsangle=new TH2D("h_reco_daughter_distvsangle", "h_reco_daughter_distvsangle", 30, 0.0, 30.0, 30, 0.0, 3.14);

              h_sel_pcosthetax = new TH1D("h_sel_pcosthetax", "h_sel_pcosthetax",20, -1.0, 1.0 );

              h_tof_momz = new TH2D("h_tof_momz","h_tof_momz", 30, 0.0, 180., 30, 0.0, 1.5);

	      h_chi2_phypo_data= new TH1D("h_chi2_phypo_data", "h_chi2_phypo_data", 100, 0, 500);
              h_PiAbs_sel_Ptmissing = new TH1D("h_PiAbs_sel_Ptmissing", "h_PiAbs_sel_Ptmissing", 15, 0.0, 1.2);
              h_PiAbs_sel_Plongit = new TH1D("h_PiAbs_sel_Plongit", "h_PiAbs_sel_Plongit", 15, 0.0, 1.2);

              h_PiAbs_sel_pmult = new TH1D ("h_PiAbs_sel_pmult", "h_PiAbs_sel_pmult", 5, -0.5, 4.5);



        } 
#endif
