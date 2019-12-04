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
  chain_pdxsec -> Add(pattern_add.c_str());

  PDEventPro * t=new PDEventPro(chain_pdxsec);

  _event_histo_1d = new PDEventHisto1D();
  _event_histo_1d->InitializeBootstraps();

  int evts = chain_pdxsec -> GetEntries();
  std::cout<<"total number of events is "<<evts<<std::endl;
  //loop over all the events
  int barWidth = 70;
  for(int i= _initial_entry; i<evts; i++){
	if(i !=0) DrawProgressBar((double)i/(double)evts, barWidth);
	chain_pdxsec->GetEntry(i);
	//-----------------------------------------------
        if(abs(t->true_beam_PDG) != 211) continue; 	
        
        //Fill true proton momentum
        int ngen_proton=0;
	for(unsigned int ind=0; ind<t->true_beam_daughter_PDG->size(); ind++){
           if(abs(t->true_beam_daughter_PDG->at(ind))==2212) {
	     ngen_proton++;
             _event_histo_1d->h_mom_gentruep->Fill(t->true_beam_daughter_startP->at(ind));
             // check if this proton is really reconstructed or not
             // loop over the reco particles and compare the track ID


              
             //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
           }
        }



        //std::cout<<t->reco_daughter_PFP_true_byHits_endProcess->size()<<" "<<t->reco_daughter_allTrack_Chi2_proton->size()<<std::endl;
        //----------------------------------------------------------
	int nsim_proton=0;
	//loop over all the daughter particles and get the chi2
        for(unsigned int ipfp=0; ipfp<t->reco_daughter_PFP_true_byHits_endProcess->size(); ipfp++){ 
           if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(ipfp))==2212) {
	     nsim_proton++;
             _event_histo_1d->h_chi2_phypo_proton->Fill(t->reco_daughter_allTrack_Chi2_proton->at(ipfp)/t->reco_daughter_allTrack_Chi2_ndof->at(ipfp));}
           else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(ipfp))==211) {_event_histo_1d->h_chi2_phypo_pionpm->Fill(t->reco_daughter_allTrack_Chi2_proton->at(ipfp)/t->reco_daughter_allTrack_Chi2_ndof->at(ipfp));}
           else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(ipfp))==13) {_event_histo_1d->h_chi2_phypo_muon->Fill(t->reco_daughter_allTrack_Chi2_proton->at(ipfp)/t->reco_daughter_allTrack_Chi2_ndof->at(ipfp));}
           else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(ipfp))==11) {_event_histo_1d->h_chi2_phypo_electron->Fill(t->reco_daughter_allTrack_Chi2_proton->at(ipfp)/t->reco_daughter_allTrack_Chi2_ndof->at(ipfp));}
           else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(ipfp))==321) {_event_histo_1d->h_chi2_phypo_kaon->Fill(t->reco_daughter_allTrack_Chi2_proton->at(ipfp)/t->reco_daughter_allTrack_Chi2_ndof->at(ipfp));}
           else{_event_histo_1d->h_chi2_phypo_other->Fill(t->reco_daughter_allTrack_Chi2_proton->at(ipfp)/t->reco_daughter_allTrack_Chi2_ndof->at(ipfp));}
        }//end of loop over all the daughter PDF particles  
        //--------------------------------------------------------------------------------
        int nsim_showerp=0;
	for(unsigned int ishw=0; ishw<t->reco_daughter_shower_true_byHits_PDG->size(); ishw++){       
           if(abs(t->reco_daughter_shower_true_byHits_PDG->at(ishw))==2212) {
	    nsim_showerp++;
            _event_histo_1d->h_chi2_phypo_proton->Fill(t->reco_daughter_shower_Chi2_proton->at(ishw)/t->reco_daughter_shower_Chi2_ndof->at(ishw));}
           else if(abs(t->reco_daughter_shower_true_byHits_PDG->at(ishw))==211) {_event_histo_1d->h_chi2_phypo_pionpm->Fill(t->reco_daughter_shower_Chi2_proton->at(ishw)/t->reco_daughter_shower_Chi2_ndof->at(ishw));}
           else if(abs(t->reco_daughter_shower_true_byHits_PDG->at(ishw))==13) {_event_histo_1d->h_chi2_phypo_muon->Fill(t->reco_daughter_shower_Chi2_proton->at(ishw)/t->reco_daughter_shower_Chi2_ndof->at(ishw));}
           else if(abs(t->reco_daughter_shower_true_byHits_PDG->at(ishw))==11) {_event_histo_1d->h_chi2_phypo_electron->Fill(t->reco_daughter_shower_Chi2_proton->at(ishw)/t->reco_daughter_shower_Chi2_ndof->at(ishw));}
           else if(abs(t->reco_daughter_shower_true_byHits_PDG->at(ishw))==321) {_event_histo_1d->h_chi2_phypo_kaon->Fill(t->reco_daughter_shower_Chi2_proton->at(ishw)/t->reco_daughter_shower_Chi2_ndof->at(ishw));}
           else{_event_histo_1d->h_chi2_phypo_other->Fill(t->reco_daughter_shower_Chi2_proton->at(ishw)/t->reco_daughter_shower_Chi2_ndof->at(ishw));}
        }  
	//------------------------------------------------
        if(nsim_proton + nsim_showerp > ngen_proton){
				std::cout<<"sim number of proton > gen number of proton "<<ngen_proton<<" "<<nsim_proton+nsim_showerp<<std::endl;
        }
        //-------------------------------------------------




        //std::cout<<"Calculate Truncated Mean dQdx"<<std::endl;
	std::vector<double> trunmeandqdx_test;
	std::vector<double> poop_par;
        //std::cout<<t->reco_daughter_PFP_true_byHits_PDG->size()<<" "<<t->reco_daughter_allTrack_dQdX->size()<<std::endl;
	for(unsigned int jj=0; jj<t->reco_daughter_PFP_true_byHits_PDG->size(); jj++){
            if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==2212){
	            _event_histo_1d->h_trackscore_proton->Fill(t->reco_daughter_PFP_trackScore->at(jj));       
                    if(t->reco_daughter_PFP_true_byHits_parID->at(jj)==t->true_beam_ID ){ 
                        _event_histo_1d->h_mom_recotruep->Fill(t->reco_daughter_PFP_true_byHits_startP->at(jj));
                    }
            } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==211){
        	    _event_histo_1d->h_trackscore_pionpm->Fill(t->reco_daughter_PFP_trackScore->at(jj));  
            } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==11){
	            _event_histo_1d->h_trackscore_electron->Fill(t->reco_daughter_PFP_trackScore->at(jj));       
            } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==22){
		    _event_histo_1d->h_trackscore_photon->Fill(t->reco_daughter_PFP_trackScore->at(jj));       
            } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==13){
        	    _event_histo_1d->h_trackscore_muon->Fill(t->reco_daughter_PFP_trackScore->at(jj));       
            } else if(abs(t->reco_daughter_PFP_true_byHits_PDG->at(jj))==111){
	            _event_histo_1d->h_trackscore_pion0->Fill(t->reco_daughter_PFP_trackScore->at(jj));       
            } else {
        	    _event_histo_1d->h_trackscore_other->Fill(t->reco_daughter_PFP_trackScore->at(jj));       
            }

            if(t->reco_daughter_PFP_trackScore->at(jj)<0.2) continue;
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
	}// end of loop over all the PFP_true_
        //-----------------------------------------------------------------------------------
        //Fill the distance between PFP daughters and end point of beam particles
	TVector3 PFP_daughter_start_v3;
	TVector3 beam_end_v3;

	beam_end_v3.SetXYZ(t->reco_beam_endX, t->reco_beam_endY, t->reco_beam_endZ);
	for(unsigned int kk=0; kk<t->reco_daughter_allTrack_ID->size(); kk++){
		PFP_daughter_start_v3.SetXYZ(t->reco_daughter_allTrack_startX->at(kk), t->reco_daughter_allTrack_startY->at(kk), t->reco_daughter_allTrack_startZ->at(kk));	

                _event_histo_1d->h_trktovtx ->Fill( (PFP_daughter_start_v3 - beam_end_v3).Mag()); 
                //check if the daughter particle is really daughter from beam particle
                if(t->reco_daughter_PFP_true_byHits_parID->at(kk)==t->true_beam_ID ){ 
                _event_histo_1d->h_trktovtx_true ->Fill( (PFP_daughter_start_v3 - beam_end_v3).Mag()); 
                }  
        }
        //------------------------------------------------------------



	//-----------------------------------------------------------  
  } //end of loop over all the events


  file_out->cd();

  file_out->WriteObject(_event_histo_1d, "PDEventHisto1D");

  file_out->Write();
 
  file_out->Close();







  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  std::cout << std::endl << std::endl;
  std::cout << "Computing time: " << elapsed_secs << " seconds = " << elapsed_secs/60. << " minutes." << endl << endl;
  return; 

}
#endif
