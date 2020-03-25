#ifndef __MAIN_MAKER_CXX__
#define __MAIN_MAKER_CXX__

#include "Maker.h"

void Main::Maker::SetInputFile(std::string in)
{
  filen = in;
}
void Main::Maker::SetInputFile_add(std::string in)
{
  filen_add = in;
}


void Main::Maker::SetOutputFile(std::string in)
{
  fileoutn = in;
}

void Main::Maker::SetEntries(int e)
{
  maxEntries = e;
}

void Main::Maker::SetInitialEntry(int e)
{
  _initial_entry = e;
}
void Main::Maker::DrawProgressBar(double progress, double barWidth) {
  
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
double Calc_reco_beam_energy(TVector3 momentum, double mass){
  return TMath::Sqrt(mass*mass+momentum.Mag()*momentum.Mag());
}


TVector3 Calc_reco_beam_momentum(vector<double> *dedx, TVector3 *dir3){
      //loop over all the dEdx and get the energy, calculate momentum
      double PionMass=0.13957;
      double totalE=0.0;
      for(unsigned int ihit=0; ihit<dedx->size(); ihit++){
          totalE=dedx->at(ihit);
      }


      double total_mom=TMath::Sqrt(totalE*totalE-PionMass*PionMass);
      TVector3 momentumv3;
      momentumv3.SetXYZ(total_mom*dir3->Px(), total_mom*dir3->Py(), total_mom*dir3->Pz());
      return momentumv3;
}


bool Main::Maker::inFV(double x, double y, double z) {
  if (x>-330. && x<330. && y>50. && y<550. && z>50. && z<645.) {return true; }
  else {return false;}
}
void Main::Maker::SetFVCut(bool v){
    FVcuton = v;
}

double Main::Maker::Ecalcmiss(double Esum, double PTmiss, int np) {
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

double Main::Maker::thetax(double theta,double phi){
  TVector3 v;
  v.SetMagThetaPhi(1,theta,phi);
  TVector3 x_axis(1,0,0);
  double theta_x = v.Angle(x_axis);
  return theta_x;
}


bool Main::Maker::isBeamType(int reco_beam_type){
  return (reco_beam_type == 13);
};

bool Main::Maker::manual_beamPos_mc(double beam_startX, double beam_startY,
                            double beam_startZ, double beam_dirX,
                            double beam_dirY,   double beam_dirZ, 
                            double true_dirX,   double true_dirY,
                            double true_dirZ,   double true_startX,
                            double true_startY, double true_startZ) {
  double projectX = (true_startX + -1*true_startZ*(true_dirX/true_dirZ) );
  double projectY = (true_startY + -1*true_startZ*(true_dirY/true_dirZ) );
  double cos = true_dirX*beam_dirX + true_dirY*beam_dirY + true_dirZ*beam_dirZ;

  if ( (beam_startX - projectX) < xlow )
    return false;
  
  if ( (beam_startX - projectX) > xhigh )
    return false;

  if ( (beam_startY - projectY) < ylow )
    return false;

  if ( (beam_startY - projectY) > yhigh )
    return false;
  
  if (beam_startZ < zlow || zhigh < beam_startZ)
    return false;
  
  if ( cos < coslow)
    return false;

  return true;

};
bool Main::Maker::endAPA3(double reco_beam_endZ){
  return(reco_beam_endZ < cutAPA3_Z);

}

bool Main::Maker::has_shower_nHits_distance(const std::vector<double> &track_score,
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



void Main::Maker::MakeFile(){
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

  _event_histo_1d = new PDEventHisto1D();
  _event_histo_1d->InitializeBootstraps();
  _event_histo_1d->InitializeHistograms();

  int evts = chain_pdxsec -> GetEntries();
  std::cout<<"total number of events is "<<evts<<std::endl;
  //loop over all the events
  int barWidth = 70;

  int Noriabs = 0;
  int Norichx = 0;
  int Norirea = 0;

  int Nshwcutabs = 0;
  int Nshwcutchx = 0;
  int Nshwcutrea = 0;

  int Ntmcutabs = 0;
  int Ntmcutchx = 0;
  int Ntmcutrea = 0;



  int Nabs = 0;
  int Nchx = 0;
  int Nrea = 0;
  int TestGenSig=0;
  //for(int i= _initial_entry; i<50; i++){

  for(int i= _initial_entry; i<evts; i++){
	if(i !=0) DrawProgressBar((double)i/(double)evts, barWidth);

	chain_pdxsec->GetEntry(i);
	//-----------------------------------------------
        if(abs(t->true_beam_PDG) != 211) continue; 	
        if(!isBeamType(t->reco_beam_type)) continue;
        if(!manual_beamPos_mc(t->reco_beam_startX, t->reco_beam_startY, t->reco_beam_startZ, 
                              t->reco_beam_trackDirX, t->reco_beam_trackDirY, t->reco_beam_trackDirZ,
                              t->true_beam_startDirX, t->true_beam_startDirY, t->true_beam_startDirZ,
                              t->true_beam_startX, t->true_beam_startY, t->true_beam_startZ)) continue;
        if(!endAPA3(t->reco_beam_endZ) )continue;

        if(FVcuton && t->reco_beam_true_byHits_matched == 0) continue;

       
        //std::cout<<"Start to check the generated and reconstructed particles "<<t->true_beam_endProcess_libo<<std::endl;
        
        //Fill true proton momentum
        unsigned int ngen_proton=0;
        unsigned int ngen_par=0;
        
        bool allgenpinFV = false;
        //bool nolowmomcand = true;
	for(unsigned int ind=0; ind<t->true_beam_daughter_PDG->size(); ind++){
           //if(t->true_beam_daughter_startP->at(ind)<0.1) {nolowmomcand = false;}
           if(abs(t->true_beam_daughter_PDG->at(ind))==2212) 
           {_event_histo_1d->h_orimom_proton->Fill(t->true_beam_daughter_startP->at(ind));}
           else if(abs(t->true_beam_daughter_PDG->at(ind))==2112) 
           {_event_histo_1d->h_orimom_neutron->Fill(t->true_beam_daughter_startP->at(ind));}
           else if(abs(t->true_beam_daughter_PDG->at(ind))==211) 
           {_event_histo_1d->h_orimom_pionpm->Fill(t->true_beam_daughter_startP->at(ind));}
           else if(abs(t->true_beam_daughter_PDG->at(ind))==111) 
           {_event_histo_1d->h_orimom_pion0->Fill(t->true_beam_daughter_startP->at(ind));}
           else if(abs(t->true_beam_daughter_PDG->at(ind))==321) 
           {_event_histo_1d->h_orimom_kaon->Fill(t->true_beam_daughter_startP->at(ind));}
           else if(abs(t->true_beam_daughter_PDG->at(ind))==11) 
           {_event_histo_1d->h_orimom_electron->Fill(t->true_beam_daughter_startP->at(ind));}
           else if(abs(t->true_beam_daughter_PDG->at(ind))==13) 
           {_event_histo_1d->h_orimom_muon->Fill(t->true_beam_daughter_startP->at(ind));
           }
           else if(abs(t->true_beam_daughter_PDG->at(ind))==22) 
           {_event_histo_1d->h_orimom_photon->Fill(t->true_beam_daughter_startP->at(ind));}
           else {_event_histo_1d->h_orimom_other->Fill(t->true_beam_daughter_startP->at(ind));
           }
           ngen_par++;
           if(abs(t->true_beam_daughter_PDG->at(ind))!=2212) continue;
 		     ngen_proton++;
        }


        std::cout<<t->reco_beam_true_byHits_matched<<std::endl;
        if(FVcuton && t->reco_beam_true_byHits_matched == 0) continue;

        std::cout<<"Find a matched events with true pions"<<std::endl; 
        bool isSignal = false;
        bool isChxBKG = false;
        bool isReaBKG = false;
        if(t->true_daughter_nPi0 == 0 && t->true_daughter_nPiPlus == 0 
           && t->true_daughter_nPiMinus ==0 /*&& t->true_daughter_endProcess == "pi+Inelastic"*/ ) {
		Noriabs++;
		isSignal = true;
        }
        else if(t->true_daughter_nPi0 > 0 ) {Norichx++;
                isChxBKG = true;
        }
        if(t->true_daughter_nPi0 == 0 && (t->true_daughter_nPiPlus > 0 
           || t->true_daughter_nPiMinus >0)  ) {Norirea++;
                isReaBKG = true;
        }

        //================================================================


        //get the denomator of the energetic protons
        if(t->true_daughter_nPi0 == 0 && t->true_daughter_nPiPlus == 0 && t->true_daughter_nPiMinus ==0 )      {
          int temp_genpindex=-999; double temp_genpmom=-999.0;
          for(unsigned int ntrk=0; ntrk<t->true_beam_daughter_startP->size(); ntrk++){
            if(abs(t->true_beam_daughter_PDG->at(ntrk)) !=2212) continue;
            if(t->true_beam_daughter_startP->at(ntrk)>temp_genpmom){
               temp_genpindex = ntrk;
               temp_genpmom=t->true_beam_daughter_startP->at(ntrk);
            }           
          }
          if(temp_genpindex<0) continue;
            _event_histo_1d->h_PiAbs_gen_sig_energeticproton_mom->Fill(t->true_beam_daughter_startP->at(temp_genpindex));
            //_event_histo_1d->h_PiAbs_gen_sig_energeticproton_costheta->Fill(TMath::Cos(t->true_beam_daughter_Theta->at(temp_genpindex)));
            _event_histo_1d->h_PiAbs_gen_sig_energeticproton_costheta->Fill(t->true_beam_daughter_startPz->at(temp_genpindex)/t->true_beam_daughter_startP->at(temp_genpindex));
            //_event_histo_1d->h_PiAbs_gen_sig_energeticproton_phi->Fill(t->true_beam_daughter_Phi->at(temp_genpindex));
            _event_histo_1d->h_PiAbs_gen_sig_energeticproton_phi->Fill(TMath::ATan2(t->true_beam_daughter_startPy->at(temp_genpindex),t->true_beam_daughter_startPx->at(temp_genpindex)));

           TestGenSig++;


        }//end of selected signal before event selection

        //==================================================================



        TVector3 truebeamendv3;

        truebeamendv3.SetXYZ(t->reco_beam_endX, t->reco_beam_endY, t->reco_beam_endZ);


        TVector3 truedaughterstartv3;
        TVector3 truedaughterendv3;
        TVector3 truegranddaughterstartv3;
        TVector3 truegranddaughterendv3;


        vector<int> pargdid;
        pargdid.clear();
        for(unsigned ingd=0; ingd<t->true_beam_grand_daughter_ID->size(); ingd++){
            //truegranddaughterstartv3.SetXYZ(t->true_beam_grand_daughter_startX->at(ingd), t->true_beam_grand_daughter_startY->at(ingd),t->true_beam_grand_daughter_startZ->at(ingd));
            //truegranddaughterendv3.SetXYZ(t->true_beam_grand_daughter_endX->at(ingd), t->true_beam_grand_daughter_endY->at(ingd),t->true_beam_grand_daughter_endZ->at(ingd));
            //if((truebeamendv3-truegranddaughterstartv3).Mag()<(truebeamendv3-truegranddaughterendv3).Mag()) {
            //   _event_histo_1d->h_grand_daughter_beam_dist->Fill((truebeamendv3-truegranddaughterstartv3).Mag());
            //}
            //else {_event_histo_1d->h_grand_daughter_beam_dist->Fill((truebeamendv3-truegranddaughterendv3).Mag());}

            //if(abs(t->true_beam_grand_daughter_PDG->at(ingd))==111 || abs(t->true_beam_grand_daughter_PDG->at(ingd))==211) {
            //     pargdid.push_back(t->true_beam_grand_daughter_parID->at(ingd));
            //}
        } 
        //====================================================================
        //----------------------------------------------------------------------------  
        //check how many nucleons 
        for(unsigned indp=0; indp<t->true_beam_daughter_ID->size(); indp++){


             truedaughterstartv3.SetXYZ(t->true_beam_daughter_startX->at(indp), t->true_beam_daughter_startY->at(indp),t->true_beam_daughter_startZ->at(indp));
             truedaughterendv3.SetXYZ(t->true_beam_daughter_endX->at(indp), t->true_beam_daughter_endY->at(indp),t->true_beam_daughter_endZ->at(indp));
 
             if((truebeamendv3-truedaughterstartv3).Mag()<(truebeamendv3-truedaughterendv3).Mag()) {
               _event_histo_1d->h_daughter_beam_dist->Fill((truebeamendv3-truedaughterstartv3).Mag());
             }
             else {_event_histo_1d->h_daughter_beam_dist->Fill((truebeamendv3-truedaughterendv3).Mag());}




             std::vector<int>::iterator tempit;

             tempit = std::find(pargdid.begin(), pargdid.end(), t->true_beam_daughter_ID->at(indp));
             if(tempit !=pargdid.end()){
                          if(abs(t->true_beam_daughter_PDG->at(indp))==2212 || abs(t->true_beam_daughter_PDG->at(indp))==2112) {
                           _event_histo_1d->h_gdfromproton->Fill(t->true_beam_daughter_startP->at(indp));
                          }
                          else /*if(abs(t->true_beam_daughter_PDG->at(indp))==)*/ {
                                 _event_histo_1d->h_gdfromother->Fill(t->true_beam_daughter_startP->at(indp));
                                 //std::cout<<"GRANDDAUGHTER PDG CODE is "<<t->true_beam_daughter_PDG->at(indp)<<std::endl;
                          }
      
             } else {continue;}
                          
        }         
        //----------------------------------------------------------
	//loop over all the daughter particles and get the chi2
        for(unsigned int ipfp=0; ipfp<t->reco_daughter_allTrack_ID->size(); ipfp++){ 
           if(t->reco_daughter_PFP_trackScore_collection->at(ipfp)<0.3) continue;
           if(t->reco_daughter_allTrack_Chi2_ndof->at(ipfp)<5) continue;
             _event_histo_1d->h_chi2_phypo_mc->Fill(t->reco_daughter_allTrack_Chi2_proton->at(ipfp)/t->reco_daughter_allTrack_Chi2_ndof->at(ipfp));
 
           if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(ipfp))==2212) {
             _event_histo_1d->h_chi2_phypo_proton->Fill(t->reco_daughter_allTrack_Chi2_proton->at(ipfp)/t->reco_daughter_allTrack_Chi2_ndof->at(ipfp));

           }
           else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(ipfp))==211) {
              _event_histo_1d->h_chi2_phypo_pionpm->Fill(t->reco_daughter_allTrack_Chi2_proton->at(ipfp)/t->reco_daughter_allTrack_Chi2_ndof->at(ipfp));}
           else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(ipfp))==13) {
              _event_histo_1d->h_chi2_phypo_muon->Fill(t->reco_daughter_allTrack_Chi2_proton->at(ipfp)/t->reco_daughter_allTrack_Chi2_ndof->at(ipfp));
           }
           else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(ipfp))==11) {
              _event_histo_1d->h_chi2_phypo_electron->Fill(t->reco_daughter_allTrack_Chi2_proton->at(ipfp)/t->reco_daughter_allTrack_Chi2_ndof->at(ipfp));
           }
           else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(ipfp))==321) {
              _event_histo_1d->h_chi2_phypo_kaon->Fill(t->reco_daughter_allTrack_Chi2_proton->at(ipfp)/t->reco_daughter_allTrack_Chi2_ndof->at(ipfp));}
           else{_event_histo_1d->h_chi2_phypo_other->Fill(t->reco_daughter_allTrack_Chi2_proton->at(ipfp)/t->reco_daughter_allTrack_Chi2_ndof->at(ipfp));}
        }//end of loop over all the daughter PDF particles  
        //--------------------------------------------------------------------------------


        //std::cout<<"Calculate Truncated Mean dQdx"<<std::endl;
	std::vector<double> trunmeandqdx_test;
	std::vector<double> poop_par;

        //======================================================================
        //loop over all the true beam daughters and fill the histogram for the start momentum
        //of pi0, pipm and protons
     
        for(unsigned int indt=0; indt<t->true_beam_daughter_PDG->size(); indt++){

           if(abs(t->true_beam_daughter_PDG->at(indt)) == 111) { 
              _event_histo_1d->h_mom_gentruepion0->Fill(t->true_beam_daughter_startP->at(indt));
           }
           if(abs(t->true_beam_daughter_PDG->at(indt)) == 211) { 
              _event_histo_1d->h_mom_gentruepionpm->Fill(t->true_beam_daughter_startP->at(indt));
              auto piCSearch = std::find (t->reco_daughter_PFP_true_byHits_ID->begin(), t->reco_daughter_PFP_true_byHits_ID->end(), t->true_beam_daughter_ID->at(indt) );
              if(piCSearch !=t->reco_daughter_PFP_true_byHits_ID->end() ){
                  _event_histo_1d->h_mom_recotruepionpm->Fill(t->true_beam_daughter_startP->at(indt));
              }
           }


           if(abs(t->true_beam_daughter_PDG->at(indt))!= 2212) continue;
	   _event_histo_1d->h_mom_gentruep->Fill(t->true_beam_daughter_startP->at(indt));
           _event_histo_1d->h_thetax_gentruep->Fill(TMath::ACos(t->true_beam_daughter_startPx->at(indt)/t->true_beam_daughter_startP->at(indt)));
           //check if this proton is reaconstructed, then get the initial reconstruction efficiency
           int foundreco=0;
           auto itSearch = std::find( t->reco_daughter_PFP_true_byHits_ID->begin(), t->reco_daughter_PFP_true_byHits_ID->end(), t->true_beam_daughter_ID->at(indt) );
           if( itSearch != t->reco_daughter_PFP_true_byHits_ID->end() ){
                 foundreco++;
                 _event_histo_1d->h_mom_recotruep->Fill(t->true_beam_daughter_startP->at(indt));
           } else { 
            //study the non-reco protons            
            //if(t->true_beam_daughter_startP->at(indt)> 0.3 && t->true_beam_daughter_nHits->at(indt)  > 0 )
            {
                  _event_histo_1d->h_lownhitsp_thetax->Fill(t->true_beam_daughter_startPx->at(indt)/t->true_beam_daughter_startP->at(indt));
                  _event_histo_1d->h_lownhitsp_thetay->Fill(t->true_beam_daughter_startPy->at(indt)/t->true_beam_daughter_startP->at(indt));
                  _event_histo_1d->h_lownhitsp_thetaz->Fill(t->true_beam_daughter_startPz->at(indt)/t->true_beam_daughter_startP->at(indt));

                  _event_histo_1d->h_lownhitsp_startX->Fill(t->true_beam_daughter_startX->at(indt)); 
                  _event_histo_1d->h_lownhitsp_startY->Fill(t->true_beam_daughter_startY->at(indt)); 
                  _event_histo_1d->h_lownhitsp_startZ->Fill(t->true_beam_daughter_startZ->at(indt)); 

                  _event_histo_1d->h_lownhitsp_endX->Fill(t->true_beam_daughter_endX->at(indt)); 
                  _event_histo_1d->h_lownhitsp_endY->Fill(t->true_beam_daughter_endY->at(indt)); 
                  _event_histo_1d->h_lownhitsp_endZ->Fill(t->true_beam_daughter_endZ->at(indt)); 


                  //std::cout<<"libotest "<<t->run<<"   "<<t->subrun<<"   "<<t->event<<"  "
                  //std::cout<<t->true_beam_daughter_endProcess->at(indt)<<"  "
                  //std::cout<<t->true_beam_daughter_startP->at(indt)<<"  "
                  //std::cout<<t->true_beam_daughter_nHits->at(indt)<<std::endl;    
                  _event_histo_1d->h_nonrecop_momvsnhits->Fill(t->true_beam_daughter_startP->at(indt), t->true_beam_daughter_nHits->at(indt));
                  _event_histo_1d->h_mom_nonrecop->Fill(t->true_beam_daughter_startP->at(indt));
            }
           }
        }
        //======================================================================
        /*
        */
        //=========================================================================

        vector<int> vec_parid; //find out the parID of the photons but no loop
        vec_parid.clear();

        for(unsigned int mm=0; mm<t->reco_daughter_PFP_true_byHits_PDG->size(); mm++){
             if(t->reco_daughter_PFP_true_byHits_PDG->at(mm) !=22) continue;
              
             std::vector<int>::iterator it;

             it = std::find(vec_parid.begin(), vec_parid.end(), t->reco_daughter_PFP_true_byHits_parID->at(mm));
             if (it != vec_parid.end()) 
             {
               
             } else {
                   vec_parid.push_back(t->reco_daughter_PFP_true_byHits_parID->at(mm));
             }
        }     
        
        /*
        */

        vector<int> vec_pionparid;
        vec_pionparid.clear();
        //----------------------------------------------------------------------- 
        double Pgamma = 0.;
        double Egamma = 0.;
        //loop over all the parID from reco photons 
        if(vec_parid.size()>0){
        for(unsigned int ind_parid=0; ind_parid<vec_parid.size(); ind_parid++){
          // loop over all the true beam daughter ID
          for(unsigned int ind_cd=0; ind_cd<t->true_beam_daughter_ID->size(); ind_cd++){
               if(t->true_beam_daughter_ID->at(ind_cd) == vec_parid.at(ind_parid) && t->true_beam_daughter_PDG->at(ind_cd)==111) {
                        vec_pionparid.push_back(vec_parid.at(ind_parid));
                        _event_histo_1d->h_mom_recotruepion0->Fill(t->true_beam_daughter_startP->at(ind_cd));
	       } else if(t->true_beam_daughter_ID->at(ind_cd) == vec_parid.at(ind_parid)){
                        _event_histo_1d->h_mom_recotruenonpi0->Fill(t->true_beam_daughter_startP->at(ind_cd));
               }
          }
        }
        }
        // vec pion0 parID produce photons  
        if(vec_pionparid.size()>0){
        for(unsigned int ind_parid=0; ind_parid<vec_pionparid.size(); ind_parid++){
          Pgamma =0.;
          double Pxgamma = 0.0;
          double Pygamma = 0.0; 
          double Pzgamma = 0.0;
          int Ngamma=0;
          for(unsigned int hh=0; hh<t->reco_daughter_PFP_true_byHits_PDG->size(); hh++){
            if(t->reco_daughter_PFP_true_byHits_PDG->at(hh) !=22) continue;
            if(t->reco_daughter_PFP_true_byHits_parID->at(hh) != vec_pionparid.at(ind_parid)) continue;
            Pgamma += t->reco_daughter_PFP_true_byHits_startP->at(hh);
            Egamma += t->reco_daughter_PFP_true_byHits_startE->at(hh);
            Pxgamma += t->reco_daughter_PFP_true_byHits_startPx->at(hh);
            Pygamma += t->reco_daughter_PFP_true_byHits_startPy->at(hh);
            Pzgamma += t->reco_daughter_PFP_true_byHits_startPz->at(hh);
            Ngamma++;
          } // end of all the reco pf objects
          _event_histo_1d->h_ngamma_frompi0->Fill(Ngamma);
          if(Ngamma >=2 && Pgamma>0){
          _event_histo_1d->h_pgamma_frompi0->Fill(TMath::Sqrt(Egamma*Egamma-(Pxgamma*Pxgamma+Pygamma*Pygamma+Pzgamma*Pzgamma)));
                
          }
        } //loop over all the vector of parid
        }
        //======================================================================== 
        //Fill the distance between PFP daughters and end point of beam particles
	


        int nonprotoncand = 0;
        int nshwcand = 0;
	int ntmdqdxcand = 0;
        TVector3 recobeamendv3, recodauendv3, recodaustartv3, recogranddaustartv3, recogranddauendv3;
        recobeamendv3.SetXYZ(t->reco_beam_endX, t->reco_beam_endY, t->reco_beam_endZ);
	for(unsigned int jj=0; jj<t->reco_daughter_PFP_true_byHits_PDG->size(); jj++){
            if(FVcuton) {
            // if( !allgenpinFV ) continue;
            }


	    if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==2212){
	            _event_histo_1d->h_trackscore_proton->Fill(t->reco_daughter_PFP_trackScore_collection->at(jj));       
                    if(t->reco_daughter_PFP_true_byHits_parID->at(jj) == t->true_beam_ID){
                    _event_histo_1d->h_mom_recotruep_test->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj));
                    _event_histo_1d->h_thetax_recotruep_test->Fill(TMath::ACos(t->reco_daughter_PFP_true_byHits_startPx->at(jj)/t->reco_daughter_PFP_true_byHits_startP->at(jj)));
                    }
            } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==211){
        	    _event_histo_1d->h_trackscore_pionpm->Fill(t->reco_daughter_PFP_trackScore_collection->at(jj));  
            } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==11){
	            _event_histo_1d->h_trackscore_electron->Fill(t->reco_daughter_PFP_trackScore_collection->at(jj));       
            } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==22){
		    _event_histo_1d->h_trackscore_photon->Fill(t->reco_daughter_PFP_trackScore_collection->at(jj));       
            } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==13){
        	    _event_histo_1d->h_trackscore_muon->Fill(t->reco_daughter_PFP_trackScore_collection->at(jj));       
            } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==111){
	            _event_histo_1d->h_trackscore_pion0->Fill(t->reco_daughter_PFP_trackScore_collection->at(jj));       
            } else {
        	    _event_histo_1d->h_trackscore_other->Fill(t->reco_daughter_PFP_trackScore_collection->at(jj));       
            }


            //---------------------------------------------------------------------

            if(t->reco_daughter_PFP_trackScore_collection->at(jj)<0.3 /*&& t->reco_daughter_allTrack_to_vertex->at(jj)*/ ) {nshwcand ++; } 
            TVector3 recostartmom;
            //double calc_reco_daughter_allTrack_startPx =t->reco_daughter_allTrack_momByRange_proton->at(jj)*TMath::Sin(t->reco_daughter_allTrack_Theta->at(jj))*TMath::Cos(t->reco_daughter_allTrack_Phi->at(jj));
            //double calc_reco_daughter_allTrack_startPy =t->reco_daughter_allTrack_momByRange_proton->at(jj)*TMath::Sin(t->reco_daughter_allTrack_Theta->at(jj))*TMath::Sin(t->reco_daughter_allTrack_Phi->at(jj));
            //double calc_reco_daughter_allTrack_startPz =t->reco_daughter_allTrack_momByRange_proton->at(jj)*TMath::Cos(t->reco_daughter_allTrack_Theta->at(jj));



            /*
            
            */
            if(t->reco_daughter_PFP_trackScore_collection->at(jj)<0.3) continue;
            if(t->reco_daughter_PFP_true_byHits_parID->at(jj) == t->true_beam_ID){
                if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==2212) {
                      _event_histo_1d->h_mom_trkscoretruep->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj));
                }
                if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==211) {
                      _event_histo_1d->h_mom_trkscoretruepipm->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj));
                }
                   
            }
            poop_par.clear();
	    if(t->reco_daughter_allTrack_Chi2_proton->at(jj)<-998.) continue;
            //get the dQdx of the all Tracks object and Fill in a vector
	    for(unsigned int kk=0; kk<t->reco_daughter_allTrack_dQdX->at(jj).size(); kk++){
		poop_par.push_back(t->reco_daughter_allTrack_dQdX->at(jj).at(kk));
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


             

            _event_histo_1d->h_tmdqdx->Fill(TMath::Mean(TLMean_par.begin(), TLMean_par.end()));
	    if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==2212){
	    _event_histo_1d->h_tmdqdx_proton->Fill(TMath::Mean(TLMean_par.begin(), TLMean_par.end()));        
            //_event_histo_1d->h_tmdqdxvsrange_proton->Fill(TMath::Mean(TLMean_par.begin(), TLMean_par.end()), t->reco_daughter_allTrack_len->at(jj));

	    }else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==211){
	    _event_histo_1d->h_tmdqdx_pionpm->Fill(TMath::Mean(TLMean_par.begin(), TLMean_par.end()));        
            //_event_histo_1d->h_tmdqdxvsrange_pionpm->Fill(TMath::Mean(TLMean_par.begin(), TLMean_par.end()), t->reco_daughter_allTrack_len->at(jj));
	    }
            //------------------------------------------------------------------

            bool isPCand = false;
            bool isPiCand = false;
            double temp_tmdqdx=TMath::Mean(TLMean_par.begin(), TLMean_par.end());
            //for temp_dqdx>600 proton cand
            if(temp_tmdqdx>600) {
               isPCand = true;
               if(t->reco_daughter_PFP_true_byHits_parID->at(jj) == t->true_beam_ID){
                 if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==2212) {
                   _event_histo_1d->h_mom_selectedtruep->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj));
                   _event_histo_1d->h_mom_tmdqdxtruep->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj));
                  }
               }
            } 

            if(temp_tmdqdx>250 && temp_tmdqdx<400) {//pion region
               isPiCand = true;
               ntmdqdxcand++;
               if(t->reco_daughter_PFP_true_byHits_parID->at(jj) == t->true_beam_ID){
                  if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==211){
                     _event_histo_1d->h_mom_selectedtruepipm->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj));
                     _event_histo_1d->h_mom_tmdqdxtruepipm->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj));
                  }
               }
            }
            if((temp_tmdqdx<600 && temp_tmdqdx>400) || temp_tmdqdx<250) {  //transition region
                if(t->reco_daughter_PFP_true_byHits_parID->at(jj) == t->true_beam_ID){
                  if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==2212) {
                     _event_histo_1d->h_mom_tmdqdxtruep->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj));
                  }
                  if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==211){
                     _event_histo_1d->h_mom_tmdqdxtruepipm->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj));
                  }
                }//selected the primary daughters
                if(t->reco_daughter_allTrack_Chi2_proton->at(jj)/t->reco_daughter_allTrack_Chi2_ndof->at(jj) < 40
                   && t->reco_daughter_allTrack_Chi2_proton->at(jj)/t->reco_daughter_allTrack_Chi2_ndof->at(jj) > 0 ) 
                {
                  isPCand = true;
                  if(t->reco_daughter_PFP_true_byHits_parID->at(jj) == t->true_beam_ID){
                     if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==2212) {              
                        _event_histo_1d->h_mom_selectedtruep->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj));
                     }//end of if the tracks with chi2>40
                  }//end of selecting the true daughter particles 
                } //end of apply chi2 cut
                else if(t->reco_daughter_allTrack_Chi2_proton->at(jj)/t->reco_daughter_allTrack_Chi2_ndof->at(jj) > 40
                   && t->reco_daughter_allTrack_Chi2_proton->at(jj)/t->reco_daughter_allTrack_Chi2_ndof->at(jj) > 0 ) 
                { 
                 isPiCand = true;
                  if(t->reco_daughter_PFP_true_byHits_parID->at(jj) ==t->true_beam_ID) {
                     if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==211){
                         _event_histo_1d->h_mom_selectedtruepipm->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj));
                     }   
                  }
                }
                else {
                    isPiCand = true;
                    if(t->reco_daughter_PFP_true_byHits_parID->at(jj) == t->true_beam_ID){
                    if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==2212) {
                    if(t->reco_daughter_PFP_true_byHits_startP->at(jj) > 1.0) {
                        _event_histo_1d->h_chi2cutp_thetax->Fill(t->reco_daughter_PFP_true_byHits_startPx->at(jj)/t->reco_daughter_PFP_true_byHits_startP->at(jj));
                        _event_histo_1d->h_chi2cutp_thetay->Fill(t->reco_daughter_PFP_true_byHits_startPy->at(jj)/t->reco_daughter_PFP_true_byHits_startP->at(jj));
                        _event_histo_1d->h_chi2cutp_thetaz->Fill(t->reco_daughter_PFP_true_byHits_startPz->at(jj)/t->reco_daughter_PFP_true_byHits_startP->at(jj));
                        
                    }}} 
                }

            } //end of performing tmdqdx cut
            //if(isPiCand == true){
            //  ntmdqdxcand++;
            //}
            if(isPCand == true){
                //if((PFP_daughter_start_v3 - beam_end_v3).Mag()<10)
                {
                /*if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==2212){
                    _event_histo_1d->h_recomom_selected_proton->Fill(t->reco_daughter_allTrack_momByRange_proton->at(jj));
                } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==211){
                    _event_histo_1d->h_recomom_selected_pionpm->Fill(t->reco_daughter_allTrack_momByRange_proton->at(jj));
                } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==11){
                    _event_histo_1d->h_recomom_selected_electron->Fill(t->reco_daughter_allTrack_momByRange_proton->at(jj));
                } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==13){
                    _event_histo_1d->h_recomom_selected_muon->Fill(t->reco_daughter_allTrack_momByRange_proton->at(jj));
                } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==321){
                    _event_histo_1d->h_recomom_selected_kaon->Fill(t->reco_daughter_allTrack_momByRange_proton->at(jj));
                } else{
                    _event_histo_1d->h_recomom_selected_other->Fill(t->reco_daughter_allTrack_momByRange_proton->at(jj));
                }*/
  
                if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==2212){
                    _event_histo_1d->h_recomom_selected_proton->Fill(t->reco_daughter_allTrack_len->at(jj));
                } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==211){
                    _event_histo_1d->h_recomom_selected_pionpm->Fill(t->reco_daughter_allTrack_len->at(jj));
                } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==11){
                    _event_histo_1d->h_recomom_selected_electron->Fill(t->reco_daughter_allTrack_len->at(jj));
                } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==13){
                    _event_histo_1d->h_recomom_selected_muon->Fill(t->reco_daughter_allTrack_len->at(jj));
                } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==321){
                    _event_histo_1d->h_recomom_selected_kaon->Fill(t->reco_daughter_allTrack_len->at(jj));
                } else{
                    _event_histo_1d->h_recomom_selected_other->Fill(t->reco_daughter_allTrack_len->at(jj));
                }
                //----------------------------------------------------------------------------
                std::vector<int>::iterator temprecoit;

                temprecoit=std::find(t->true_beam_daughter_ID->begin(), t->true_beam_daughter_ID->end(), t->reco_daughter_PFP_true_byHits_parID->at(jj));
                if(temprecoit != t->true_beam_daughter_ID->end()){
                  recodaustartv3.SetXYZ(t->reco_daughter_allTrack_startX->at(jj),t->reco_daughter_allTrack_startY->at(jj), t->reco_daughter_allTrack_startZ->at(jj));
                  recodauendv3.SetXYZ(t->reco_daughter_allTrack_endX->at(jj),t->reco_daughter_allTrack_endY->at(jj), t->reco_daughter_allTrack_endZ->at(jj));
                  //recostartmom.SetXYZ(calc_reco_daughter_allTrack_startPx,calc_reco_daughter_allTrack_startPy,calc_reco_daughter_allTrack_startPz);


                  if((recobeamendv3-recodaustartv3).Mag()<(recobeamendv3-recodauendv3).Mag()){
                  _event_histo_1d->h_reco_grand_daughter_beam_dist->Fill((recobeamendv3-recodaustartv3).Mag());
                  //_event_histo_1d->h_reco_grand_daughter_angle->Fill(recostartmom.Angle(-recobeamendv3+recodaustartv3));
                  }
                  else {_event_histo_1d->h_reco_grand_daughter_beam_dist->Fill((recobeamendv3-recodauendv3).Mag());
                        //_event_histo_1d->h_reco_grand_daughter_angle->Fill(recostartmom.Angle(-recobeamendv3+recodauendv3));
                  }

                } else {
                  recogranddaustartv3.SetXYZ(t->reco_daughter_allTrack_startX->at(jj),t->reco_daughter_allTrack_startY->at(jj), t->reco_daughter_allTrack_startZ->at(jj));
                  recogranddauendv3.SetXYZ(t->reco_daughter_allTrack_endX->at(jj),t->reco_daughter_allTrack_endY->at(jj), t->reco_daughter_allTrack_endZ->at(jj));
                  //recostartmom.SetXYZ(calc_reco_daughter_allTrack_startPx,calc_reco_daughter_allTrack_startPy,calc_reco_daughter_allTrack_startPz);


                  if((recobeamendv3-recogranddaustartv3).Mag()<(recobeamendv3-recogranddauendv3).Mag()){
                   _event_histo_1d->h_reco_daughter_beam_dist->Fill((recobeamendv3-recogranddaustartv3).Mag());
                   //_event_histo_1d->h_reco_daughter_angle->Fill(recostartmom.Angle(-recobeamendv3+recodaustartv3));

                   _event_histo_1d->h_reco_daughter_deltax->Fill(t->reco_daughter_allTrack_startX->at(jj) - t->reco_beam_endX);
                   _event_histo_1d->h_reco_daughter_deltay->Fill(t->reco_daughter_allTrack_startY->at(jj) - t->reco_beam_endY);
                   _event_histo_1d->h_reco_daughter_deltaz->Fill(t->reco_daughter_allTrack_startZ->at(jj) - t->reco_beam_endZ);




                   //_event_histo_1d->h_reco_daughter_distvsangle->Fill((recobeamendv3-recogranddaustartv3).Mag(),recostartmom.Angle(-recobeamendv3+recodaustartv3));
                  }
                  else { _event_histo_1d->h_reco_daughter_beam_dist->Fill((recobeamendv3-recogranddauendv3).Mag());
                         //_event_histo_1d->h_reco_daughter_angle->Fill(recostartmom.Angle(-recobeamendv3+recodauendv3));

                  }
                }     
                //check the resolution of the selected protons

                double temp_KE = 0.0; double temp_Momentum = 0.0; double trkrange = 0.0;
                if(t->reco_daughter_allTrack_len->at(jj) > 0 && t->reco_daughter_allTrack_len->at(jj) < 80){
                     trkrange = t->reco_daughter_allTrack_len->at(jj); 
                     temp_KE = 29.9317 * std::pow(trkrange, 0.586304);
                } else if (t->reco_daughter_allTrack_len->at(jj) > 80 && t->reco_daughter_allTrack_len->at(jj) < 3.022E3){

                     trkrange = t->reco_daughter_allTrack_len->at(jj);
                temp_KE = 149.904 + (3.34146 * trkrange) + (-0.00318856 * trkrange * trkrange) +
                (4.34587E-6 * trkrange * trkrange * trkrange) +
                (-3.18146E-9 * trkrange * trkrange * trkrange * trkrange) +
                (1.17854E-12 * trkrange * trkrange * trkrange * trkrange * trkrange) +
                (-1.71763E-16 * trkrange * trkrange * trkrange * trkrange * trkrange *trkrange);
                }
                
                double M = ProtonMass*1000;
                temp_Momentum = TMath::Sqrt((temp_KE*temp_KE) + 2*M*temp_KE); 
                temp_Momentum = temp_Momentum/1000.0;
                _event_histo_1d->h_selproton_momreso->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj),temp_Momentum);

                _event_histo_1d->h_selproton_costhetareso->Fill(t->reco_daughter_PFP_true_byHits_startPz->at(jj)/t->reco_daughter_PFP_true_byHits_startP->at(jj), TMath::Cos(t->reco_daughter_allTrack_Theta->at(jj)));

                //-----------------------------------------------------------------------------
                }      
            } //end of if pcand is true
            else if(isPCand == false ) {nonprotoncand++;}
            if(!isPCand && isPiCand){
                if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==2212){
                    _event_histo_1d->h_recopimom_selected_proton->Fill(t->reco_daughter_allTrack_len->at(jj));
                } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==211){
                    _event_histo_1d->h_recopimom_selected_pionpm->Fill(t->reco_daughter_allTrack_len->at(jj));
                } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==11){
                    _event_histo_1d->h_recopimom_selected_electron->Fill(t->reco_daughter_allTrack_len->at(jj));
                } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==13){
                    _event_histo_1d->h_recopimom_selected_muon->Fill(t->reco_daughter_allTrack_len->at(jj));
                } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==321){
                    _event_histo_1d->h_recopimom_selected_kaon->Fill(t->reco_daughter_allTrack_len->at(jj));
                } else{
                    _event_histo_1d->h_recopimom_selected_other->Fill(t->reco_daughter_allTrack_len->at(jj));
                }

            }


            //--------------------------------------------------------------------
	}// end of loop over all the PFP_true_
        
        /*
        */

        if(nshwcand > 0) continue;
        if(isSignal) {Nshwcutabs++; }
        if(isChxBKG) {Nshwcutchx++; } 
        if(isReaBKG) {Nshwcutrea++; }

        if(ntmdqdxcand>0) continue;
        if(isSignal) {Ntmcutabs++; }
        if(isChxBKG) {Ntmcutchx++; } 
        if(isReaBKG) {Ntmcutrea++; }

        

        if(nonprotoncand > 0 || nshwcand > 0) continue;
        //calculate Emissing Pmissing here

        double true_Ptmissing = 0; double Emissing_Qsubtracted=0;
        double reco_Ptmissing = 0; double reco_Pz = 0;
        double total_truepPx=0; double total_truepPy=0; 
        double total_pPx = 0; double total_pPy = 0; double total_pPz=0;
        double total_truepKE=0.0; double total_pKE=0;
        //const double ProtonMass = 0.938272;
        const double PionMass = 0.13957;




        double temp_pmom=-999.0; double temp_pcostheta=-999.0; double temp_pphi=-999.0; double temp_pcosthetax=-999.;
        double temp_plen = -999.0;

        for(unsigned int ntrk=0; ntrk<t->reco_daughter_allTrack_len->size(); ntrk++){

             _event_histo_1d->h_PiAbs_sel_pmom->Fill(t->reco_daughter_allTrack_len->at(ntrk));
             //_event_histo_1d->h_PiAbs_sel_pmom->Fill(t->reco_daughter_allTrack_momByRange_proton->at(ntrk));
             _event_histo_1d->h_PiAbs_sel_pcostheta->Fill(TMath::Cos(t->reco_daughter_allTrack_Theta->at(ntrk)));
             _event_histo_1d->h_PiAbs_sel_ptheta->Fill(t->reco_daughter_allTrack_Theta->at(ntrk));

             _event_histo_1d->h_PiAbs_sel_pcosthetax->Fill(TMath::Cos(thetax(t->reco_daughter_allTrack_Theta->at(ntrk), t->reco_daughter_allTrack_Phi->at(ntrk))));
             _event_histo_1d->h_PiAbs_sel_pphi->Fill(t->reco_daughter_allTrack_Phi->at(ntrk));

             total_truepPx += t->reco_daughter_PFP_true_byHits_startPx->at(ntrk);
             total_truepPy += t->reco_daughter_PFP_true_byHits_startPy->at(ntrk);
             total_truepKE += t->reco_daughter_PFP_true_byHits_startE->at(ntrk) - ProtonMass;

             //calculate totaldedx here and sum over all the dEdx to get the total energy loss

             double temp_KE = 0.0; double temp_Momentum = 0.0; double trkrange = 0.0;
             if(t->reco_daughter_allTrack_len->at(ntrk) > 0 && t->reco_daughter_allTrack_len->at(ntrk) < 80){
                     trkrange = t->reco_daughter_allTrack_len->at(ntrk); 
                     temp_KE = 29.9317 * std::pow(trkrange, 0.586304);
             } else if (t->reco_daughter_allTrack_len->at(ntrk) > 80 && t->reco_daughter_allTrack_len->at(ntrk) < 3.022E3){

                     trkrange = t->reco_daughter_allTrack_len->at(ntrk);
                temp_KE = 149.904 + (3.34146 * trkrange) + (-0.00318856 * trkrange * trkrange) +
                (4.34587E-6 * trkrange * trkrange * trkrange) +
                (-3.18146E-9 * trkrange * trkrange * trkrange * trkrange) +
                (1.17854E-12 * trkrange * trkrange * trkrange * trkrange * trkrange) +
                (-1.71763E-16 * trkrange * trkrange * trkrange * trkrange * trkrange *trkrange);
             }
                
             double M = ProtonMass*1000;
             temp_Momentum = TMath::Sqrt((temp_KE*temp_KE) + 2*M*temp_KE); 
             temp_Momentum = temp_Momentum/1000.0;
 


             total_pPx += temp_Momentum*TMath::Sin(t->reco_daughter_allTrack_Theta->at(ntrk))*TMath::Cos(t->reco_daughter_allTrack_Phi->at(ntrk));
             total_pPy += temp_Momentum*TMath::Sin(t->reco_daughter_allTrack_Theta->at(ntrk))*TMath::Sin(t->reco_daughter_allTrack_Phi->at(ntrk));
             total_pPz += temp_Momentum*TMath::Cos(t->reco_daughter_allTrack_Theta->at(ntrk));



            

             total_pKE = temp_KE/1000.0;

             if(t->reco_daughter_allTrack_len->at(ntrk) > temp_plen){
                                temp_plen=t->reco_daughter_allTrack_len->at(ntrk);
                                temp_pmom = temp_Momentum;
                                temp_pcostheta=TMath::Cos(t->reco_daughter_allTrack_Theta->at(ntrk));
                                temp_pphi=t->reco_daughter_allTrack_Phi->at(ntrk);
                                temp_pcosthetax=TMath::Cos(thetax(t->reco_daughter_allTrack_Theta->at(ntrk), t->reco_daughter_allTrack_Phi->at(ntrk)));
             }

        }
        _event_histo_1d->h_PiAbs_sel_energeticproton_mom->Fill(temp_pmom);
        _event_histo_1d->h_PiAbs_sel_energeticproton_costheta->Fill(temp_pcostheta);
        _event_histo_1d->h_PiAbs_sel_energeticproton_phi->Fill(temp_pphi);
        _event_histo_1d->h_PiAbs_sel_energeticproton_costhetax->Fill(temp_pcosthetax);



        true_Ptmissing = TMath::Sqrt(total_truepPx*total_truepPx + total_truepPy*total_truepPy);
        reco_Ptmissing = TMath::Sqrt(total_pPx*total_pPx + total_pPy*total_pPy);
        reco_Pz = total_pPz; 
        Emissing_Qsubtracted = TMath::Sqrt(t->true_beam_startP*t->true_beam_startP + PionMass*PionMass) - total_truepKE;
        
        


        if(t->true_daughter_nPi0 == 0 && t->true_daughter_nPiPlus == 0 && t->true_daughter_nPiMinus ==0 ) {Nabs++;

                 _event_histo_1d->h_PiAbs_sig_pmult->Fill(t->reco_daughter_allTrack_ID->size());
                 _event_histo_1d->h_PiAbs_sig_nmult->Fill(t->true_daughter_nNeutron);
                 _event_histo_1d->h_PiAbs_sig_truevsreco_pmult->Fill(t->true_daughter_nProton, t->reco_daughter_allTrack_ID->size());

                 if( t->reco_daughter_allTrack_ID->size() == 3){
                       std::cout<<"checking events with 3 reco protons "<<t->true_daughter_nProton<<"     "<<t->true_daughter_nNeutron<<std::endl;
                 }

                 //true to get the true particle multiplicity
                 _event_histo_1d->h_sig_Emissing->Fill(Emissing_Qsubtracted);
                 if(t->reco_daughter_allTrack_ID->size()){
                       _event_histo_1d->h_sig_Ptmissing->Fill(reco_Ptmissing);
                       _event_histo_1d->h_sig_Plongit->Fill(reco_Pz);
                 }
                 _event_histo_1d->h_sig_pvsnmult->Fill(t->true_daughter_nProton, t->true_daughter_nNeutron);

                 //get the total kinetic energy of protons and neutrons



                 temp_plen = -999.0; temp_pmom=-999.0;  temp_pcostheta=-999.0;  temp_pphi=-999.0;

                 for(unsigned int tmk=0; tmk<t->reco_daughter_allTrack_ID->size(); tmk++){
                      //_event_histo_1d->h_PiAbs_sig_pmom->Fill(t->reco_daughter_allTrack_momByRange_proton->at(tmk)); 
                      _event_histo_1d->h_PiAbs_sig_pmom->Fill(t->reco_daughter_allTrack_len->at(tmk)); 
                      _event_histo_1d->h_PiAbs_sig_pcostheta->Fill(TMath::Cos(t->reco_daughter_allTrack_Theta->at(tmk)));
                      _event_histo_1d->h_PiAbs_sig_pphi->Fill(t->reco_daughter_allTrack_Phi->at(tmk));
                      _event_histo_1d->h_PiAbs_sig_ptheta->Fill(t->reco_daughter_allTrack_Theta->at(tmk));

                      double temp_KE = 0.0; double temp_Momentum = 0.0; double trkrange = 0.0;
                      if(t->reco_daughter_allTrack_len->at(tmk) > 0 && t->reco_daughter_allTrack_len->at(tmk) < 80){
                       trkrange = t->reco_daughter_allTrack_len->at(tmk); 
                       temp_KE = 29.9317 * std::pow(trkrange, 0.586304);
                      } else if (t->reco_daughter_allTrack_len->at(tmk) > 80 && t->reco_daughter_allTrack_len->at(tmk) < 3.022E3){

                       trkrange = t->reco_daughter_allTrack_len->at(tmk);
                       temp_KE = 149.904 + (3.34146 * trkrange) + (-0.00318856 * trkrange * trkrange) +
                        (4.34587E-6 * trkrange * trkrange * trkrange) +
                        (-3.18146E-9 * trkrange * trkrange * trkrange * trkrange) +
                        (1.17854E-12 * trkrange * trkrange * trkrange * trkrange * trkrange) +
                        (-1.71763E-16 * trkrange * trkrange * trkrange * trkrange * trkrange *trkrange);
                      }
                
                      double M = ProtonMass*1000;
                      temp_Momentum = TMath::Sqrt((temp_KE*temp_KE) + 2*M*temp_KE); 
                      temp_Momentum = temp_Momentum/1000.0;
 


                      if(t->reco_daughter_allTrack_len->at(tmk) > temp_plen){
                                temp_plen=t->reco_daughter_allTrack_len->at(tmk);
                                temp_pmom=temp_Momentum;
                                temp_pcostheta=TMath::Cos(t->reco_daughter_allTrack_Theta->at(tmk));
                                temp_pphi=t->reco_daughter_allTrack_Phi->at(tmk);
                       }

                 }//loop over all the proton candidates

                 _event_histo_1d->h_PiAbs_sig_energeticproton_reco_mom->Fill(temp_pmom);
                 _event_histo_1d->h_PiAbs_sig_energeticproton_reco_costheta->Fill(temp_pcostheta);
                 _event_histo_1d->h_PiAbs_sig_energeticproton_reco_phi->Fill(temp_pphi);



                 double temp_pmom2=-999.0; 
                 int temp_selpindex=-999;
                 
                 double totalPxtest=0.0;
                 double totalPytest=0.0;
                 double totalPxtest_onlyproton=0.0;
                 double totalPytest_onlyproton=0.0;
              
                 double totalKE_proton=0.;
                 double totalKE_neutron=0.;
                 for(unsigned int tmks=0; tmks<t->true_beam_daughter_startP->size(); tmks++){
                      totalPxtest += t->true_beam_daughter_startPx->at(tmks);
                      totalPytest += t->true_beam_daughter_startPy->at(tmks);
                      if(abs(t->true_beam_daughter_PDG->at(tmks) ) ==2212) {

                          totalPxtest_onlyproton += t->true_beam_daughter_startPx->at(tmks);
                          totalPytest_onlyproton += t->true_beam_daughter_startPy->at(tmks);
 
                          totalKE_proton += TMath::Sqrt(t->true_beam_daughter_startP->at(tmks)*t->true_beam_daughter_startP->at(tmks)+ProtonMass*ProtonMass) - ProtonMass;
                          _event_histo_1d->h_PiAbs_sig_proton_mom->Fill(t->true_beam_daughter_startP->at(tmks));
                          _event_histo_1d->h_PiAbs_sig_proton_costheta->Fill(t->true_beam_daughter_startPz->at(tmks)/t->true_beam_daughter_startP->at(tmks));
                          
                      }

                      if(abs(t->true_beam_daughter_PDG->at(tmks) ) ==2112) {
                          totalKE_neutron +=TMath::Sqrt(t->true_beam_daughter_startP->at(tmks)*t->true_beam_daughter_startP->at(tmks)+NeutronMass*NeutronMass) - NeutronMass;
                          _event_histo_1d->h_PiAbs_sig_neutron_mom->Fill(t->true_beam_daughter_startP->at(tmks));
                          _event_histo_1d->h_PiAbs_sig_neutron_costheta->Fill(t->true_beam_daughter_startPz->at(tmks)/t->true_beam_daughter_startP->at(tmks));
                       }

                      if(abs(t->true_beam_daughter_PDG->at(tmks)) !=2212) continue;
                      if(t->true_beam_daughter_startP->at(tmks)>temp_pmom2){
                          temp_selpindex = tmks;
                          temp_pmom2=t->true_beam_daughter_startP->at(tmks);
                      }
                 }


                 _event_histo_1d->h_PiAbs_sig_true_totalKE_protonvsneutron->Fill(totalKE_proton, totalKE_neutron);

                 _event_histo_1d->h_PiAbs_sig_true_Ptmissing->Fill(TMath::Sqrt(totalPxtest*totalPxtest+totalPytest*totalPytest));
                 _event_histo_1d->h_PiAbs_sig_true_Ptmissing_onlyproton->Fill(TMath::Sqrt(totalPxtest_onlyproton*totalPxtest_onlyproton+totalPytest_onlyproton*totalPytest_onlyproton));


                 if(temp_selpindex<0) continue;
                 _event_histo_1d->h_PiAbs_sig_energeticproton_mom->Fill(t->true_beam_daughter_startP->at(temp_selpindex));
                 //_event_histo_1d->h_PiAbs_sig_energeticproton_costheta->Fill(TMath::Cos(t->true_beam_daughter_Theta->at(temp_selpindex)));
                 _event_histo_1d->h_PiAbs_sig_energeticproton_costheta->Fill(t->true_beam_daughter_startPz->at(temp_selpindex)/t->true_beam_daughter_startP->at(temp_selpindex));
                 _event_histo_1d->h_PiAbs_sig_energeticproton_phi->Fill(TMath::ATan2(t->true_beam_daughter_startPy->at(temp_selpindex), t->true_beam_daughter_startPx->at(temp_selpindex)));

                 _event_histo_1d->h_PiAbs_sig_energeticproton_truevsreco_mom->Fill(t->true_beam_daughter_startP->at(temp_selpindex),temp_pmom);
                 //_event_histo_1d->h_PiAbs_sig_energeticproton_truevsreco_costheta->Fill(TMath::Cos(t->true_beam_daughter_Theta->at(temp_selpindex)),temp_pcostheta);
                 _event_histo_1d->h_PiAbs_sig_energeticproton_truevsreco_costheta->Fill(t->true_beam_daughter_startPz->at(temp_selpindex)/t->true_beam_daughter_startP->at(temp_selpindex),temp_pcostheta);
                 if(t->true_beam_daughter_startPy->at(temp_selpindex)>0){
                 _event_histo_1d->h_PiAbs_sig_energeticproton_truevsreco_phi->Fill(TMath::ACos(t->true_beam_daughter_startPx->at(temp_selpindex)/TMath::Sqrt(t->true_beam_daughter_startPx->at(temp_selpindex)*t->true_beam_daughter_startPx->at(temp_selpindex)+t->true_beam_daughter_startPy->at(temp_selpindex)*t->true_beam_daughter_startPy->at(temp_selpindex))),temp_pphi);
                 } else {
                 _event_histo_1d->h_PiAbs_sig_energeticproton_truevsreco_phi->Fill(-TMath::ACos(t->true_beam_daughter_startPx->at(temp_selpindex)/TMath::Sqrt(t->true_beam_daughter_startPx->at(temp_selpindex)*t->true_beam_daughter_startPx->at(temp_selpindex)+t->true_beam_daughter_startPy->at(temp_selpindex)*t->true_beam_daughter_startPy->at(temp_selpindex))),temp_pphi);

                 }
        }

 




        else if(t->true_daughter_nPi0>0) {Nchx++;
             _event_histo_1d->h_PiAbs_chxbac_pmult->Fill(t->reco_daughter_allTrack_ID->size());
             int Nreco_pion0=0;
             temp_plen = -999.0; temp_pmom=-999.0;  temp_pcostheta=-999.0;  temp_pphi=-999.0;
             for(unsigned int tmk=0; tmk<t->reco_daughter_PFP_true_byHits_PDG->size(); tmk++){
                 //_event_histo_1d->h_PiAbs_chxbac_pmom->Fill(t->reco_daughter_allTrack_momByRange_proton->at(tmk)); 
                 _event_histo_1d->h_PiAbs_chxbac_pmom->Fill(t->reco_daughter_allTrack_len->at(tmk)); 
                 _event_histo_1d->h_PiAbs_chxbac_pcostheta->Fill(TMath::Cos(t->reco_daughter_allTrack_Theta->at(tmk)));
                 _event_histo_1d->h_PiAbs_chxbac_pphi->Fill(t->reco_daughter_allTrack_Phi->at(tmk));
                 _event_histo_1d->h_PiAbs_chxbac_ptheta->Fill(t->reco_daughter_allTrack_Theta->at(tmk));

                 double temp_KE = 0.0; double temp_Momentum = 0.0; double trkrange = 0.0;
                 if(t->reco_daughter_allTrack_len->at(tmk) > 0 && t->reco_daughter_allTrack_len->at(tmk) < 80){
                     trkrange = t->reco_daughter_allTrack_len->at(tmk); 
                     temp_KE = 29.9317 * std::pow(trkrange, 0.586304);
                 } else if (t->reco_daughter_allTrack_len->at(tmk) > 80 && t->reco_daughter_allTrack_len->at(tmk) < 3.022E3){

                     trkrange = t->reco_daughter_allTrack_len->at(tmk);
                   temp_KE = 149.904 + (3.34146 * trkrange) + (-0.00318856 * trkrange * trkrange) +
                   (4.34587E-6 * trkrange * trkrange * trkrange) +
                   (-3.18146E-9 * trkrange * trkrange * trkrange * trkrange) +
                   (1.17854E-12 * trkrange * trkrange * trkrange * trkrange * trkrange) +
                   (-1.71763E-16 * trkrange * trkrange * trkrange * trkrange * trkrange *trkrange);
                 }
                
                 double M = ProtonMass*1000;
                 temp_Momentum = TMath::Sqrt((temp_KE*temp_KE) + 2*M*temp_KE); 
                 temp_Momentum = temp_Momentum/1000.0;
 



                 if(t->reco_daughter_allTrack_ID->size()){
                      _event_histo_1d->h_chxbac_Ptmissing->Fill(reco_Ptmissing);  
                      _event_histo_1d->h_chxbac_Plongit->Fill(reco_Pz);
                 }
                 if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(tmk) == 22)){
                     Nreco_pion0++;
                     double chi2_perdof= t->reco_daughter_allTrack_Chi2_proton->at(tmk)/t->reco_daughter_allTrack_Chi2_ndof->at(tmk);
                     std::cout<<"BKG_REA "<<t->event<<"  "<<t->run;
                     std::cout<<" "<<chi2_perdof<<"   ";
                     std::cout<<t->reco_daughter_allTrack_Chi2_proton->at(tmk)<<"  "<<t->reco_daughter_allTrack_Chi2_ndof->at(tmk)<<std::endl; 
                     
                 }
                 //get the most energetic proton momentum and angle
                 if(t->reco_daughter_allTrack_len->at(tmk) > temp_plen){
                                temp_plen=t->reco_daughter_allTrack_len->at(tmk);
                                temp_pmom=temp_Momentum;
                                temp_pcostheta=TMath::Cos(t->reco_daughter_allTrack_Theta->at(tmk));
                                temp_pphi=t->reco_daughter_allTrack_Phi->at(tmk);
                 }

             }//end of loop over all the final state particles

             _event_histo_1d->h_PiAbs_chxbac_energeticproton_mom->Fill(temp_pmom);
             _event_histo_1d->h_PiAbs_chxbac_energeticproton_costheta->Fill(temp_pcostheta);
             _event_histo_1d->h_PiAbs_chxbac_energeticproton_phi->Fill(temp_pphi);



             _event_histo_1d->h_reco_photon_rea->Fill(Nreco_pion0);
 

        }
        else if(t->true_daughter_nPiPlus > 0 || t->true_daughter_nPiMinus >0) {Nrea++;
             _event_histo_1d->h_PiAbs_reabac_pmult->Fill(t->reco_daughter_allTrack_ID->size());
             if(t->reco_daughter_allTrack_ID->size()){
                      _event_histo_1d->h_reabac_Ptmissing->Fill(reco_Ptmissing);
                      _event_histo_1d->h_reabac_Plongit->Fill(reco_Pz);  
             }
             int Nreco_pionpm=0;
             temp_plen=-999.0; temp_pmom=-999.0;  temp_pcostheta=-999.0;  temp_pphi=-999.0;
             for(unsigned int tmk=0; tmk<t->reco_daughter_PFP_true_byHits_PDG->size(); tmk++){
                 //_event_histo_1d->h_PiAbs_reabac_pmom->Fill(t->reco_daughter_allTrack_momByRange_proton->at(tmk)); 
                 _event_histo_1d->h_PiAbs_reabac_pmom->Fill(t->reco_daughter_allTrack_len->at(tmk)); 
                 _event_histo_1d->h_PiAbs_reabac_pcostheta->Fill(TMath::Cos(t->reco_daughter_allTrack_Theta->at(tmk)));
                 _event_histo_1d->h_PiAbs_reabac_pphi->Fill(t->reco_daughter_allTrack_Phi->at(tmk));
                 _event_histo_1d->h_PiAbs_reabac_ptheta->Fill(t->reco_daughter_allTrack_Theta->at(tmk));


                 if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(tmk) == 211)){
                     Nreco_pionpm++;
                     double chi2_perdof= t->reco_daughter_allTrack_Chi2_proton->at(tmk)/t->reco_daughter_allTrack_Chi2_ndof->at(tmk);
                     std::cout<<"BKG_REA "<<t->event<<"  "<<t->run;
                     std::cout<<" "<<chi2_perdof<<"   ";
                     std::cout<<t->reco_daughter_allTrack_Chi2_proton->at(tmk)<<"  "<<t->reco_daughter_allTrack_Chi2_ndof->at(tmk)<<std::endl; 
                     
                 }
                 //get the most energetic proton momentum and angle
                 double temp_KE = 0.0; double temp_Momentum = 0.0; double trkrange = 0.0;
                 if(t->reco_daughter_allTrack_len->at(tmk) > 0 && t->reco_daughter_allTrack_len->at(tmk) < 80){
                     trkrange = t->reco_daughter_allTrack_len->at(tmk); 
                     temp_KE = 29.9317 * std::pow(trkrange, 0.586304);
                 } else if (t->reco_daughter_allTrack_len->at(tmk) > 80 && t->reco_daughter_allTrack_len->at(tmk) < 3.022E3){

                     trkrange = t->reco_daughter_allTrack_len->at(tmk);
                   temp_KE = 149.904 + (3.34146 * trkrange) + (-0.00318856 * trkrange * trkrange) +
                   (4.34587E-6 * trkrange * trkrange * trkrange) +
                   (-3.18146E-9 * trkrange * trkrange * trkrange * trkrange) +
                   (1.17854E-12 * trkrange * trkrange * trkrange * trkrange * trkrange) +
                   (-1.71763E-16 * trkrange * trkrange * trkrange * trkrange * trkrange *trkrange);
                 }
                
                 double M = ProtonMass*1000;
                 temp_Momentum = TMath::Sqrt((temp_KE*temp_KE) + 2*M*temp_KE); 
                 temp_Momentum = temp_Momentum/1000.0;
 
                 if(t->reco_daughter_allTrack_len->at(tmk) > temp_plen){
                                temp_plen=t->reco_daughter_allTrack_len->at(tmk);
                                temp_pmom=temp_Momentum;
                                temp_pcostheta=TMath::Cos(t->reco_daughter_allTrack_Theta->at(tmk));
                                temp_pphi=t->reco_daughter_allTrack_Phi->at(tmk);
                 }

             }//end of loop over all the final state particles
             _event_histo_1d->h_PiAbs_reabac_energeticproton_mom->Fill(temp_pmom);
             _event_histo_1d->h_PiAbs_reabac_energeticproton_costheta->Fill(temp_pcostheta);
             _event_histo_1d->h_PiAbs_reabac_energeticproton_phi->Fill(temp_pphi);

             _event_histo_1d->h_reco_pionpm_rea->Fill(Nreco_pionpm);
        }
        //-----------------------------------------------------------------------------------
        //------------------------------------------------------------

        std::cout<<"End of processing a certain event "<<std::endl;

	//-----------------------------------------------------------  

       

  } //end of loop over all the events


  
  file_out->cd();

  //file_out->WriteObject(_event_histo_1d, "PDEventHisto1D");


  /*
  */
  TH1D  *selected_signal_percut = new TH1D("selected_signal_percut", "selected_percut", 4, 0, 4);

  TH1D  *selected_percut = new TH1D("selected_percut", "selected_percut", 4, 0, 4);

  TH1D  *generated_signal_percut = new TH1D("generated_signal_percut", "generated_signal_percut", 4, 0, 4);



  /*_event_histo_1d->*/selected_percut->SetBinContent(1, Noriabs+Norichx+Norirea);
  /*_event_histo_1d->*/selected_percut->SetBinContent(2, Nshwcutabs+Nshwcutchx+Nshwcutrea);
  /*_event_histo_1d->*/selected_percut->SetBinContent(3, Ntmcutabs+Ntmcutchx+Ntmcutrea);
  /*_event_histo_1d->*/selected_percut->SetBinContent(4, Nabs+Nchx+Nrea);



  /*_event_histo_1d->*/selected_signal_percut->SetBinContent(1, Noriabs);
  /*_event_histo_1d->*/selected_signal_percut->SetBinContent(2, Nshwcutabs);
  /*_event_histo_1d->*/selected_signal_percut->SetBinContent(3, Ntmcutabs);
  /*_event_histo_1d->*/selected_signal_percut->SetBinContent(4, Nabs);

  /*_event_histo_1d->*/generated_signal_percut->SetBinContent(1, Noriabs);
  /*_event_histo_1d->*/generated_signal_percut->SetBinContent(2, Noriabs);
  /*_event_histo_1d->*/generated_signal_percut->SetBinContent(3, Noriabs);
  /*_event_histo_1d->*/generated_signal_percut->SetBinContent(4, Noriabs);


  std::vector<std::string> cut_names = {"initial", "TrackScore", "TrunMeandQdx", "Chi2"};
  for(int i=0; i<4; i++){
      std::cout<<cut_names.at(i)<<" & "<</*_event_histo_1d->*/selected_signal_percut->GetBinContent(i+1)<<std::endl;
  }


  TCanvas * canvas_eff_pur_graph_percut = new TCanvas();

  canvas_eff_pur_graph_percut->SetLeftMargin(0.05157593);
  canvas_eff_pur_graph_percut->SetRightMargin(0.1475645);
  canvas_eff_pur_graph_percut->SetTopMargin(0.04210526);
  canvas_eff_pur_graph_percut->SetBottomMargin(0.1578947);

  TH1F *h = new TH1F("h","",4, 0, 4);
  h->SetMaximum(1);
  h->GetXaxis()->SetBinLabel(1,"Initial");
  h->GetXaxis()->SetBinLabel(2,"TrackScore");
  h->GetXaxis()->SetBinLabel(3,"TrunMeandQdx");
  h->GetXaxis()->SetBinLabel(4,"Chi2");

  h->GetXaxis()->SetLabelOffset(0.009);
  h->GetXaxis()->SetLabelSize(0.06);

  h->Draw();  


  TEfficiency* pEff_percut = new TEfficiency(*selected_signal_percut,*generated_signal_percut);
  pEff_percut->SetTitle("EfficiencyPerCut;Cut index;Efficiency");
  pEff_percut->SetLineColor(kGreen+3);
  pEff_percut->SetMarkerColor(kGreen+3);
  pEff_percut->SetMarkerStyle(20);
  pEff_percut->SetMarkerSize(0.6);
  TGraphAsymmErrors * pEff_percut_graph = pEff_percut->CreateGraph();
  for (int i = 0; i < 4; i++) {
    pEff_percut_graph->SetPointEXhigh(i, 0.);
    pEff_percut_graph->SetPointEXlow(i, 0.);
  }
  auto axis = pEff_percut_graph->GetYaxis();
  axis->SetLimits(0.,1.); 

  pEff_percut_graph->GetHistogram()->GetXaxis()->SetBinLabel(1,"Initial");
  pEff_percut_graph->GetHistogram()->GetXaxis()->SetBinLabel(2,"TrackScore");
  pEff_percut_graph->GetHistogram()->GetXaxis()->SetBinLabel(3,"TrunMeandQdx");
  pEff_percut_graph->GetHistogram()->GetXaxis()->SetBinLabel(4,"Chi2");
 
  pEff_percut_graph->Draw("PL"); 



	  
  TEfficiency* pPur_percut = new TEfficiency(/*_event_histo_1d->*/*selected_signal_percut, /*_event_histo_1d->*/*selected_percut);
  pPur_percut->SetTitle("Purity Per Cut; Cut index; Purity");
  pPur_percut->SetLineColor(kRed+3);
  pPur_percut->SetMarkerColor(kRed+3);
  pPur_percut->SetMarkerStyle(20);
  pPur_percut->SetMarkerSize(0.6);
  TGraphAsymmErrors * pPur_percut_graph = pPur_percut->CreateGraph();
  for (int i = 0; i < 4; i++) {
    pPur_percut_graph->SetPointEXhigh(i, 0.);
    pPur_percut_graph->SetPointEXlow(i, 0.);
  }

  pPur_percut_graph->Draw("PL");
  TLegend* l = new TLegend(0.4842407,0.8168421,0.777937,0.9221053,NULL,"brNDC");
  l->AddEntry(pEff_percut_graph,"Efficiency");
  l->AddEntry(pPur_percut_graph,"Purity");
  //leg->AddEntry(gr3,"Neutrino MCFlash","l");
  //  
  l->Draw();

  canvas_eff_pur_graph_percut->SaveAs("eff_pur_test.png");

   

  file_out->Write();
	 
  file_out->Close();

	  
  std::cout<<"Is FV cut on?? "<<FVcuton<<std::endl;

  std::cout<<"Noriabs = "<<Noriabs<<std::endl;
  std::cout<<"Norichx = "<<Norichx<<std::endl;
  std::cout<<"Norirea = "<<Norirea<<std::endl;
  std::cout<<"TestGenSig= "<<TestGenSig<<std::endl;

  std::cout<<"Nshwcutabs = "<<Nshwcutabs<<std::endl;
  std::cout<<"Nshwcutchx = "<<Nshwcutchx<<std::endl;
  std::cout<<"Nshwcutrea = "<<Nshwcutrea<<std::endl;

  std::cout<<"Ntmcutabs = "<<Ntmcutabs<<std::endl;
  std::cout<<"Ntmcutchx = "<<Ntmcutchx<<std::endl;
  std::cout<<"Ntmcutrea = "<<Ntmcutrea<<std::endl;


  std::cout<<"Nabs = "<<Nabs<<std::endl;
  std::cout<<"Nchx = "<<Nchx<<std::endl;
  std::cout<<"Nrea = "<<Nrea<<std::endl;
  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  std::cout << std::endl << std::endl;
  std::cout << "Computing time: " << elapsed_secs << " seconds = " << elapsed_secs/60. << " minutes." << endl << endl;
  return; 

}
#endif
