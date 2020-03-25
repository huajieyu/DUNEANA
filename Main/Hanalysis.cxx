#ifndef __MAIN_HANALYSIS_CXX__
#define __MAIN_HANALYSIS_CXX__

#include "Hanalysis.h"
void Main::Hanalysis::SetInputFile(std::string in)
{
  filen = in;
}
void Main::Hanalysis::SetInputFile_add(std::string in)
{
  filen_add = in;
}


void Main::Hanalysis::SetOutputFile(std::string in)
{
  fileoutn = in;
}

void Main::Hanalysis::SetEntries(int e)
{
  maxEntries = e;
}

void Main::Hanalysis::SetInitialEntry(int e)
{
  _initial_entry = e;
}
void Main::Hanalysis::DrawProgressBar(double progress, double barWidth) {
  
  std::cout << "[";
  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos) std::cout << "=";
    else if (i == pos) std::cout << ">";
    else std::cout << " ";
  }
  std::cout << "] " << int(progress * 100.0) << " %\r";
  std::cout.flush();
}

bool Main::Hanalysis::inFV(double x, double y, double z) {
  if (x>-330. && x<330. && y>50. && y<550. && z>50. && z<645.) {return true; }
  else {return false;}
}
void Main::Hanalysis::SetFVCut(bool v){
    FVcuton = v;
}

double Main::Hanalysis::Ecalcmiss(double Esum, double PTmiss, int np) {
   Esum *= 1000; //convert to MeV
   PTmiss *= 1000; //convert to MeV
   double Eexcit = 30.4; //in MeV
   double Mass = 0; // in MeV
   if(np == 0) Mass = 37.2050e3; //Ar40
   else if(np == 1) Mass = 36.2758e3; //Ar39
   else if(np == 2) Mass = 35.3669e3; //Cl38
   else if(np == 3) Mass = 34.4201e3; //S37
   else if(np == 4) Mass = 33.4957e3; //P36
   else if(np == 5) Mass = 32.5706e3; //Si35
   else if(np == 6) Mass = 31.6539e3; //Al34
   else if(np == 7) Mass = 30.7279e3; //Mg33
   else if(np == 8) Mass = 29.8111e3; //Na32
   else if(np == 9) Mass = 28.8918e3; //Ne31
   else if(np >= 10) Mass = 27.9789e3; //F30

   float Ekinrecoil = sqrt(PTmiss*PTmiss + Mass*Mass) - Mass;
   return Esum + Eexcit + Ekinrecoil; // return result in MeV
}

double Main::Hanalysis::thetax(double theta,double phi){
  TVector3 v;
  v.SetMagThetaPhi(1,theta,phi);
  TVector3 x_axis(1,0,0);
  double theta_x = v.Angle(x_axis);
  return theta_x;
}
bool Main::Hanalysis::data_beam_PID(const std::vector<int> *pidCandidates){
  auto pid_search = std::find(pidCandidates->begin(), pidCandidates->end(), 211);
  return (pid_search !=pidCandidates->end());
}

bool Main::Hanalysis:: manual_beamPos_data (int event,            double data_startX,
                              double data_startY,   double data_startZ,
                              double data_dirX,     double data_dirY,
                              double data_dirZ,     double data_BI_X,
                              double data_BI_Y,     double data_BI_dirX,
                              double data_BI_dirY,  double data_BI_dirZ,
                              int data_BI_nMomenta, int data_BI_nTracks) {

  double deltaX = data_startX - data_BI_X;
  double deltaY = data_startY - data_BI_Y;
  double cos = data_BI_dirX*data_dirX + data_BI_dirY*data_dirY +
               data_BI_dirZ*data_dirZ;

  if(data_BI_nMomenta != 1 || data_BI_nTracks != 1)
    return false;

  if( (deltaX < data_xlow) || (deltaX > data_xhigh) )
    return false;

  if ( (deltaY < data_ylow) || (deltaY > data_yhigh) )
    return false;

  if ( (data_startZ < data_zlow) || (data_startZ > data_zhigh) )
    return false;

  if (cos < data_coslow)
    return false;

  return true;

};


bool Main::Hanalysis::endAPA3(double reco_beam_endZ){
  return(reco_beam_endZ < cutAPA3_Z);

}

bool Main::Hanalysis::has_shower_nHits_distance(const std::vector<double> &track_score,
                                    const std::vector<int> &nHits,
                                    const std::vector<double> &distance) {

  if(track_score.empty() || nHits.empty())
    return false;

  for(size_t i = 0; i < track_score.size(); ++i){
     if ((track_score[i] < cut_trackScore) &&
         (nHits[i] > cut_nHits_shower_low) &&
         (nHits[i] < cut_nHits_shower_high) && (track_score[i] != -999.) &&
         (distance[i] < cut_daughter_shower_distance_high) &&
         (distance[i] > cut_daughter_shower_distance_low)) {
       return true;
     }
  }

  return false;


}



void Main::Hanalysis::MakeFile(){
  clock_t begin = clock();
  system("mkdir -p output/");

  gROOT->SetBatch(kTRUE);
  gROOT->ProcessLine("gErrorIgnoreLevel = 2001;"); // 1001: INFO, 2001: WARNINGS, 3001: ERRORS
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  std::cout << "Opening output file with name " << fileoutn << std::endl;
  TFile *file_out = new TFile(fileoutn.c_str(),"RECREATE");
  if ( file_out->IsOpen() ) {
    std::cout << "File opened successfully." << std::endl;
  } else {
    std::cout << "File not opened (maybe not found?). File: " << fileoutn << std::endl;
    exit(0);
  }
  
  string pattern = filen;
  string pattern_add = filen_add;
  TChain *chain_pdxsec;
  chain_pdxsec = new TChain("pionana/beamana");
  chain_pdxsec -> Add(pattern.c_str());
  //chain_pdxsec -> Add(pattern_add.c_str());

  PDEventPro * t=new PDEventPro(chain_pdxsec);

  _event_histo_1d = new PDEventDataHisto1D();
  _event_histo_1d->InitializeBootstraps();
  _event_histo_1d->InitializeHistograms();

  int evts = chain_pdxsec -> GetEntries();
  std::cout<<"total number of events is "<<evts<<std::endl;
  //loop over all the events
  int barWidth = 70;
  int NDatapion=0;
  int NDatamatched=0;
  int NDataSel=0;
 
  for(int i= _initial_entry; i<evts; i++){
	if(i !=0) DrawProgressBar((double)i/(double)evts, barWidth);

	chain_pdxsec->GetEntry(i);
	//-----------------------------------------------

        //_event_histo_1d->h_tof_momz->Fill(t->reco_beam_tofs->at(0), t->reco_beam_momenta->at(0));



        //if(t->reco_beam_type != 13) continue; //choose a track like beam
        //std::cout<<"vector size of tofs is "<<t->reco_beam_tofs->size()<<std::endl;
        if(!data_beam_PID(t->data_BI_PDG_candidates)) continue;
        if(!manual_beamPos_data(t->event, t->reco_beam_startX, t->reco_beam_startY, t->reco_beam_startZ,
                                t->reco_beam_trackDirX, t->reco_beam_trackDirY, t->reco_beam_trackDirZ,
                                t->data_BI_X, t->data_BI_Y, t->data_BI_dirX, t->data_BI_dirY,
                                t->data_BI_dirZ, t->data_BI_nMomenta, t->data_BI_nTracks)) continue;
        if(!endAPA3(t->reco_beam_endZ) )continue;
                         
        NDatapion++;

        


        //if(t->reco_beam_tofs->at(0)==-1) continue;


        //if(!inFV(t->reco_beam_endX, t->reco_beam_endY, t->reco_beam_endZ)) continue;
        //_event_histo_1d->h_tof_momz->Fill(t->reco_beam_tofs->at(0), t->reco_beam_momenta->at(0)*t->reco_beam_trackEndDirZ);
        //if(t->reco_beam_endX>226) continue;

        NDatamatched++;


        _event_histo_1d->h_beamstartx->Fill(t->reco_beam_startX);
        _event_histo_1d->h_beamstarty->Fill(t->reco_beam_startY);
        _event_histo_1d->h_beamstartz->Fill(t->reco_beam_startZ);

        _event_histo_1d->h_beamendx->Fill(t->reco_beam_endX);
        _event_histo_1d->h_beamendy->Fill(t->reco_beam_endY);
        _event_histo_1d->h_beamendz->Fill(t->reco_beam_endZ);

        //==============================================================================
	//loop over all the daughter particles and get the chi2
        for(unsigned int ipfp=0; ipfp<t->reco_daughter_allTrack_Chi2_proton->size(); ipfp++){ 
            if(t->reco_daughter_PFP_trackScore_collection->at(ipfp)<0.3) continue;
            if(t->reco_daughter_allTrack_Chi2_ndof->at(ipfp)<5) continue;
            _event_histo_1d->h_chi2_phypo_data->Fill(t->reco_daughter_allTrack_Chi2_proton->at(ipfp)/t->reco_daughter_allTrack_Chi2_ndof->at(ipfp));
        }//end of loop over all the daughter PDF particles  
        //--------------------------------------------------------------------------------
	for(unsigned int ishw=0; ishw<t->reco_daughter_allTrack_Chi2_proton->size(); ishw++){       
            if(t->reco_daughter_PFP_trackScore_collection->at(ishw)<0.3) continue;
            if(t->reco_daughter_allTrack_Chi2_ndof->at(ishw)<5) continue;
            _event_histo_1d->h_chi2_phypo_data->Fill(t->reco_daughter_allTrack_Chi2_proton->at(ishw)/t->reco_daughter_allTrack_Chi2_ndof->at(ishw));
        }  
	        
        
        
        
        
        
        //===============================================================================
        bool isPCand=false;

        int recopcand = 0.0;
        for(unsigned int ntrk=0; ntrk<t->reco_daughter_allTrack_startX->size(); ntrk++){


            _event_histo_1d->h_trkstartx->Fill(t->reco_daughter_allTrack_startX->at(ntrk));
            _event_histo_1d->h_trkstarty->Fill(t->reco_daughter_allTrack_startY->at(ntrk));
            _event_histo_1d->h_trkstartz->Fill(t->reco_daughter_allTrack_startZ->at(ntrk));

            _event_histo_1d->h_trkendx->Fill(t->reco_daughter_allTrack_endX->at(ntrk));
            _event_histo_1d->h_trkendy->Fill(t->reco_daughter_allTrack_endY->at(ntrk));
            _event_histo_1d->h_trkendz->Fill(t->reco_daughter_allTrack_endZ->at(ntrk));

            //Calculate truncated mean dqdx

            _event_histo_1d->h_trk_trackscore->Fill(t->reco_daughter_PFP_trackScore_collection->at(ntrk));

  	    std::vector<double> trunmeandqdx_test;
	    std::vector<double> poop_par;


            poop_par.clear();
	    if(t->reco_daughter_allTrack_Chi2_proton->at(ntrk)<-998.) continue;
            //get the dQdx of the all Tracks object and Fill in a vector
	    for(unsigned int kk=0; kk<t->reco_daughter_allTrack_dQdX->at(ntrk).size(); kk++){
		poop_par.push_back(t->reco_daughter_allTrack_dQdX->at(ntrk).at(kk));
	    }
	    double median_par=0;		
	    double RMS_par=TMath::RMS(poop_par.begin(), poop_par.end());
	    int N_par=poop_par.size(); //poop_par number of hits
	    std::sort(poop_par.begin(), poop_par.end());
	    if(N_par<=0) continue;
	    if(N_par % 2 == 0) {median_par=(poop_par[((N_par/2)-1)]+poop_par[N_par/2])/2; }
            else if(N_par == 1) {median_par=poop_par[N_par];}
            else {median_par = poop_par[(N_par+1)/2];}

	    std::vector<double> TLMean_par; 
            TLMean_par.clear();
            for(int ii=0; ii<int(poop_par.size()); ii++){
                 if(poop_par[ii]<median_par+2*RMS_par && poop_par[ii]>median_par-2*RMS_par)
                 {TLMean_par.push_back(poop_par[ii]);}
            }
            if(t->reco_daughter_PFP_trackScore_collection->at(ntrk)>0.3){ 
                   _event_histo_1d->h_trunmean_dqdx->Fill(TMath::Mean(TLMean_par.begin(), TLMean_par.end()));
            }
            isPCand = false;         
            if(t->reco_daughter_PFP_trackScore_collection->at(ntrk)<0.3) continue;
            double temp_tmdqdx = TMath::Mean(TLMean_par.begin(), TLMean_par.end());

            if(temp_tmdqdx>600) isPCand = true;
            if(temp_tmdqdx<200 || (temp_tmdqdx >400 && temp_tmdqdx<600)) {
                if(t->reco_daughter_allTrack_Chi2_proton ->at(ntrk)/t->reco_daughter_allTrack_Chi2_ndof->at(ntrk) < 40 &&
                   t->reco_daughter_allTrack_Chi2_proton ->at(ntrk)/t->reco_daughter_allTrack_Chi2_ndof->at(ntrk) > 0 ) {
                      isPCand = true;
                }
            }   



            if(isPCand == true) {
               recopcand ++;
            }

            //-----------------------------------------------------------------------

        }
        if((unsigned)recopcand != t->reco_daughter_allTrack_len->size())  continue;


          NDataSel++;     
         //--------------Fill histogram of the proton momentum and angle


         TVector3 reco_beam_endv3; TVector3 recostartmom;
         reco_beam_endv3.SetXYZ(t->reco_beam_endX, t->reco_beam_endY, t->reco_beam_endZ);

         TVector3 reco_beam_endmomv3;

         reco_beam_endmomv3.SetXYZ(t->reco_beam_trackEndDirX, t->reco_beam_trackEndDirY, t->reco_beam_trackEndDirZ);

         TVector3 reco_tpctrack_startv3;

         double temp_pmom=-999.0; double temp_pcostheta=-999.0; double temp_pphi=-999.0; double temp_pcosthetax=-999.;
         double temp_plen=-999.0;
         //for(unsigned int npfp=0; npfp<t->reco_daughter_allTrack_momByRange_proton->size(); npfp++){
         double total_pPx = 0.0; double total_pPy=0.0; double total_pPz=0.0;
         for(unsigned int npfp=0; npfp<t->reco_daughter_allTrack_len->size(); npfp++){

           //remove cosmic tracks


           //get the most energetic proton momentum and angle
           double temp_KE = 0.0; double temp_Momentum = 0.0; double trkrange = 0.0;
           if(t->reco_daughter_allTrack_len->at(npfp) > 0 && t->reco_daughter_allTrack_len->at(npfp) < 80){
                     trkrange = t->reco_daughter_allTrack_len->at(npfp); 
                     temp_KE = 29.9317 * std::pow(trkrange, 0.586304);
           } else if (t->reco_daughter_allTrack_len->at(npfp) > 80 && t->reco_daughter_allTrack_len->at(npfp) < 3.022E3){

                     trkrange = t->reco_daughter_allTrack_len->at(npfp);
                   temp_KE = 149.904 + (3.34146 * trkrange) + (-0.00318856 * trkrange * trkrange) +
                   (4.34587E-6 * trkrange * trkrange * trkrange) +
                   (-3.18146E-9 * trkrange * trkrange * trkrange * trkrange) +
                   (1.17854E-12 * trkrange * trkrange * trkrange * trkrange * trkrange) +
                   (-1.71763E-16 * trkrange * trkrange * trkrange * trkrange * trkrange *trkrange);
           }
                
           double M = ProtonMass*1000;
           temp_Momentum = TMath::Sqrt((temp_KE*temp_KE) + 2*M*temp_KE); 
           temp_Momentum = temp_Momentum/1000.0;
 



            double calc_reco_daughter_allTrack_startPx =temp_Momentum*TMath::Sin(t->reco_daughter_allTrack_Theta->at(npfp))*TMath::Cos(t->reco_daughter_allTrack_Phi->at(npfp));
            double calc_reco_daughter_allTrack_startPy =temp_Momentum*TMath::Sin(t->reco_daughter_allTrack_Theta->at(npfp))*TMath::Sin(t->reco_daughter_allTrack_Phi->at(npfp));
            double calc_reco_daughter_allTrack_startPz =temp_Momentum*TMath::Cos(t->reco_daughter_allTrack_Theta->at(npfp));
            total_pPx += calc_reco_daughter_allTrack_startPx;
            total_pPy += calc_reco_daughter_allTrack_startPy;
            total_pPz += calc_reco_daughter_allTrack_startPz;


            if(!inFV(t->reco_daughter_allTrack_startX->at(npfp),t->reco_daughter_allTrack_startY->at(npfp), t->reco_daughter_allTrack_startZ->at(npfp)) &&
            !inFV(t->reco_daughter_allTrack_endX->at(npfp),t->reco_daughter_allTrack_endY->at(npfp), t->reco_daughter_allTrack_endZ->at(npfp))) continue; 
 
             reco_tpctrack_startv3.SetXYZ(t->reco_daughter_allTrack_startX->at(npfp),t->reco_daughter_allTrack_startY->at(npfp), t->reco_daughter_allTrack_startZ->at(npfp));



             _event_histo_1d->h_beamend_trkstart_dist->Fill((reco_beam_endv3-reco_tpctrack_startv3).Mag());
             _event_histo_1d->h_beamend_trkstart_xdist->Fill(t->reco_daughter_allTrack_startX->at(npfp) - t->reco_beam_endX); 
             _event_histo_1d->h_beamend_trkstart_ydist->Fill(t->reco_daughter_allTrack_startY->at(npfp) - t->reco_beam_endY); 
             _event_histo_1d->h_beamend_trkstart_zdist->Fill(t->reco_daughter_allTrack_startZ->at(npfp) - t->reco_beam_endZ); 
             //recostartmom.SetXYZ(calc_reco_daughter_allTrack_startPx,calc_reco_daughter_allTrack_startPy,calc_reco_daughter_allTrack_startPz);
             //_event_histo_1d->h_beamend_trkstart_angle->Fill(recostartmom.Angle(-reco_beam_endv3+reco_tpctrack_startv3));
             _event_histo_1d->h_beamend_trkstart_angle->Fill(TMath::Cos(recostartmom.Angle(reco_beam_endmomv3)));
             _event_histo_1d->h_reco_daughter_distvsangle->Fill((reco_beam_endv3-reco_tpctrack_startv3).Mag(), recostartmom.Angle(-reco_beam_endv3+reco_tpctrack_startv3));


             _event_histo_1d->h_sel_pcosthetax->Fill(TMath::Cos(thetax(t->reco_daughter_allTrack_Theta->at(npfp), t->reco_daughter_allTrack_Phi->at(npfp))));
 

             _event_histo_1d->h_costheta_vs_phi->Fill(TMath::Cos(t->reco_daughter_allTrack_Theta->at(npfp)), t->reco_daughter_allTrack_Phi->at(npfp));

             _event_histo_1d->h_sel_pcandmom->Fill(t->reco_daughter_allTrack_len->at(npfp));
             _event_histo_1d->h_sel_pcandcostheta->Fill(TMath::Cos(t->reco_daughter_allTrack_Theta->at(npfp)));
             _event_histo_1d->h_sel_pcandtheta->Fill(t->reco_daughter_allTrack_Theta->at(npfp));
             _event_histo_1d->h_sel_pcandphi->Fill(t->reco_daughter_allTrack_Phi->at(npfp));

             _event_histo_1d->h_reco_daughter_deltax->Fill(t->reco_daughter_allTrack_startX->at(npfp) - t->reco_beam_endX);
             _event_histo_1d->h_reco_daughter_deltay->Fill(t->reco_daughter_allTrack_startY->at(npfp) - t->reco_beam_endY);
             _event_histo_1d->h_reco_daughter_deltaz->Fill(t->reco_daughter_allTrack_startZ->at(npfp) - t->reco_beam_endZ);


             if(t->reco_daughter_allTrack_len->at(npfp) > temp_plen){
                                temp_plen=t->reco_daughter_allTrack_len->at(npfp);
                                temp_pmom=temp_Momentum;
                                temp_pcostheta=TMath::Cos(t->reco_daughter_allTrack_Theta->at(npfp));
                                temp_pphi=t->reco_daughter_allTrack_Phi->at(npfp);
                                temp_pcosthetax=TMath::Cos(thetax(t->reco_daughter_allTrack_Theta->at(npfp), t->reco_daughter_allTrack_Phi->at(npfp)));
             }



          }//end of loop over all the final state particles    
          _event_histo_1d->h_PiAbs_sel_energeticproton_mom->Fill(temp_pmom);
          _event_histo_1d->h_PiAbs_sel_energeticproton_costheta->Fill(temp_pcostheta);
          _event_histo_1d->h_PiAbs_sel_energeticproton_phi->Fill(temp_pphi);
          _event_histo_1d->h_PiAbs_sel_energeticproton_costhetax->Fill(temp_pcosthetax);
          _event_histo_1d->h_PiAbs_sel_pmult->Fill(t->reco_daughter_allTrack_ID->size());
          if(t->reco_daughter_allTrack_ID->size()){
             _event_histo_1d->h_PiAbs_sel_Ptmissing->Fill(TMath::Sqrt(total_pPx*total_pPx+total_pPy*total_pPy));
             _event_histo_1d->h_PiAbs_sel_Plongit->Fill(total_pPz);
          }
  }

  file_out->cd();

  //file_out->WriteObject(_event_histo_1d, "PDEventHisto1D");

  file_out->Write();
 
  file_out->Close();

 
  std::cout<<"Total Number of Pion-like event is "<<NDatapion<<std::endl;
  std::cout<<"Total Number of matched event is "<<NDatamatched<<std::endl;
  std::cout<<"Total Number of selected event is "<<NDataSel<<std::endl;

  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  std::cout << std::endl << std::endl;
  std::cout << "Computing time: " << elapsed_secs << " seconds = " << elapsed_secs/60. << " minutes." << endl << endl;
 
  return;

} //end of the memeber function MakeFile
#endif
