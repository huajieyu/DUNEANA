//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Nov 29 18:40:56 2019 by ROOT version 6.18/04
// from TTree beamana/beam analysis tree
// found on file: /pnfs/dune/persistent/users/jiangl/pionana_1128/pionana_mc0to30.root
//////////////////////////////////////////////////////////

#ifndef PDEventPro_h
#define PDEventPro_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
using namespace std;
class PDEventPro {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           subrun;
   Int_t           event;
   Int_t           MC;
   Int_t           reco_beam_type;
   Double_t        reco_beam_startX;
   Double_t        reco_beam_startY;
   Double_t        reco_beam_startZ;
   Double_t        reco_beam_endX;
   Double_t        reco_beam_endY;
   Double_t        reco_beam_endZ;
   Double_t        reco_beam_len;
   Double_t        reco_beam_trackDirX;
   Double_t        reco_beam_trackDirY;
   Double_t        reco_beam_trackDirZ;
   Double_t        reco_beam_trackEndDirX;
   Double_t        reco_beam_trackEndDirY;
   Double_t        reco_beam_trackEndDirZ;
   Double_t        reco_beam_vtxX;
   Double_t        reco_beam_vtxY;
   Double_t        reco_beam_vtxZ;
   Int_t           reco_beam_trackID;
   vector<double>  *reco_beam_dQdX;
   vector<double>  *reco_beam_dEdX;
   vector<double>  *reco_beam_calibrated_dEdX;
   vector<double>  *reco_beam_resRange;
   vector<double>  *reco_beam_TrkPitch;
   vector<double>  *reco_beam_calo_wire;
   vector<double>  *reco_beam_calo_tick;
   Int_t           reco_beam_nTrackDaughters;
   Int_t           reco_beam_nShowerDaughters;
   Bool_t          reco_beam_flipped;
   Bool_t          reco_beam_passes_beam_cuts;
   vector<int>     *reco_daughter_trackID;
   vector<double>  *reco_daughter_true_byE_completeness;
   vector<double>  *reco_daughter_true_byE_purity;
   vector<int>     *reco_daughter_true_byE_PDG;
   vector<int>     *reco_daughter_true_byE_ID;
   vector<int>     *reco_daughter_true_byE_origin;
   vector<int>     *reco_daughter_true_byE_parID;
   vector<string>  *reco_daughter_true_byE_process;
   vector<int>     *reco_daughter_true_byHits_PDG;
   vector<int>     *reco_daughter_true_byHits_ID;
   vector<int>     *reco_daughter_true_byHits_origin;
   vector<int>     *reco_daughter_true_byHits_parID;
   vector<string>  *reco_daughter_true_byHits_process;
   vector<double>  *reco_daughter_true_byHits_purity;
   vector<unsigned long> *reco_daughter_true_byHits_sharedHits;
   vector<unsigned long> *reco_daughter_true_byHits_emHits;
   vector<double>  *reco_daughter_true_byHits_len;
   vector<double>  *reco_daughter_true_byHits_startX;
   vector<double>  *reco_daughter_true_byHits_startY;
   vector<double>  *reco_daughter_true_byHits_startZ;
   vector<double>  *reco_daughter_true_byHits_endX;
   vector<double>  *reco_daughter_true_byHits_endY;
   vector<double>  *reco_daughter_true_byHits_endZ;
   vector<double>  *reco_daughter_true_byHits_startPx;
   vector<double>  *reco_daughter_true_byHits_startPy;
   vector<double>  *reco_daughter_true_byHits_startPz;
   vector<double>  *reco_daughter_true_byHits_startP;
   vector<double>  *reco_daughter_true_byHits_startE;
   vector<int>     *reco_daughter_PFP_true_byHits_PDG;
   vector<int>     *reco_daughter_PFP_true_byHits_ID;
   vector<int>     *reco_daughter_PFP_true_byHits_origin;
   vector<int>     *reco_daughter_PFP_true_byHits_parID;
   vector<string>  *reco_daughter_PFP_true_byHits_process;
   vector<unsigned long> *reco_daughter_PFP_true_byHits_sharedHits;
   vector<unsigned long> *reco_daughter_PFP_true_byHits_emHits;
   vector<double>  *reco_daughter_PFP_true_byHits_len;
   vector<double>  *reco_daughter_PFP_true_byHits_startX;
   vector<double>  *reco_daughter_PFP_true_byHits_startY;
   vector<double>  *reco_daughter_PFP_true_byHits_startZ;
   vector<double>  *reco_daughter_PFP_true_byHits_endX;
   vector<double>  *reco_daughter_PFP_true_byHits_endY;
   vector<double>  *reco_daughter_PFP_true_byHits_endZ;
   vector<double>  *reco_daughter_PFP_true_byHits_startPx;
   vector<double>  *reco_daughter_PFP_true_byHits_startPy;
   vector<double>  *reco_daughter_PFP_true_byHits_startPz;
   vector<double>  *reco_daughter_PFP_true_byHits_startP;
   vector<double>  *reco_daughter_PFP_true_byHits_startE;
   vector<string>  *reco_daughter_PFP_true_byHits_endProcess;
   vector<int>     *reco_daughter_allTrack_ID;
   vector<vector<double> > *reco_daughter_allTrack_dEdX;
   vector<vector<double> > *reco_daughter_allTrack_dQdX;
   vector<vector<double> > *reco_daughter_allTrack_resRange;
   vector<vector<double> > *reco_daughter_allTrack_dQdX_SCE;
   vector<vector<double> > *reco_daughter_allTrack_dEdX_SCE;
   vector<vector<double> > *reco_daughter_allTrack_resRange_SCE;
   vector<vector<double> > *reco_daughter_allTrack_calibrated_dEdX;
   vector<vector<double> > *reco_daughter_allTrack_calibrated_dEdX_SCE;
   vector<double>  *reco_daughter_allTrack_Chi2_proton;
   vector<int>     *reco_daughter_allTrack_Chi2_ndof;
   vector<double>  *reco_daughter_allTrack_startX;
   vector<double>  *reco_daughter_allTrack_startY;
   vector<double>  *reco_daughter_allTrack_startZ;
   vector<double>  *reco_daughter_allTrack_endX;
   vector<double>  *reco_daughter_allTrack_endY;
   vector<double>  *reco_daughter_allTrack_endZ;
   vector<double>  *reco_daughter_allTrack_dR;
   vector<double>  *reco_daughter_allTrack_to_vertex;
   vector<int>     *reco_daughter_shower_true_byE_PDG;
   vector<int>     *reco_daughter_shower_true_byE_ID;
   vector<int>     *reco_daughter_shower_true_byE_origin;
   vector<int>     *reco_daughter_shower_true_byE_parID;
   vector<double>  *reco_daughter_shower_true_byE_startPx;
   vector<double>  *reco_daughter_shower_true_byE_startPy;
   vector<double>  *reco_daughter_shower_true_byE_startPz;
   vector<double>  *reco_daughter_shower_true_byE_startP;
   vector<int>     *reco_daughter_shower_true_byHits_PDG;
   vector<int>     *reco_daughter_shower_true_byHits_ID;
   vector<int>     *reco_daughter_shower_true_byHits_origin;
   vector<int>     *reco_daughter_shower_true_byHits_parID;
   vector<string>  *reco_daughter_shower_true_byHits_process;
   vector<double>  *reco_daughter_shower_true_byHits_purity;
   vector<double>  *reco_daughter_shower_true_byHits_startPx;
   vector<double>  *reco_daughter_shower_true_byHits_startPy;
   vector<double>  *reco_daughter_shower_true_byHits_startPz;
   vector<double>  *reco_daughter_shower_true_byHits_startP;
   vector<int>     *reco_daughter_showerID;
   vector<vector<double> > *reco_daughter_dQdX;
   vector<vector<double> > *reco_daughter_dEdX;
   vector<vector<double> > *reco_daughter_resRange;
   vector<vector<double> > *reco_daughter_shower_dQdX;
   vector<vector<double> > *reco_daughter_shower_dEdX;
   vector<vector<double> > *reco_daughter_shower_resRange;
   vector<double>  *reco_daughter_len;
   vector<double>  *reco_daughter_startX;
   vector<double>  *reco_daughter_startY;
   vector<double>  *reco_daughter_startZ;
   vector<double>  *reco_daughter_endX;
   vector<double>  *reco_daughter_endY;
   vector<double>  *reco_daughter_endZ;
   vector<double>  *reco_daughter_deltaR;
   vector<double>  *reco_daughter_dR;
   vector<double>  *reco_daughter_to_vertex;
   vector<int>     *reco_daughter_slice;
   vector<double>  *reco_daughter_shower_to_vertex;
   vector<double>  *reco_daughter_shower_startX;
   vector<double>  *reco_daughter_shower_startY;
   vector<double>  *reco_daughter_shower_startZ;
   vector<double>  *reco_daughter_shower_len;
   vector<int>     *reco_daughter_PFP_ID;
   vector<double>  *reco_daughter_PFP_trackScore;
   vector<double>  *reco_daughter_PFP_emScore;
   vector<double>  *reco_daughter_PFP_michelScore;
   Int_t           true_beam_PDG;
   Int_t           true_beam_ID;
   string          *true_beam_endProcess;
   Double_t        true_beam_endX;
   Double_t        true_beam_endY;
   Double_t        true_beam_endZ;
   Double_t        true_beam_startX;
   Double_t        true_beam_startY;
   Double_t        true_beam_startZ;
   Double_t        true_beam_startPx;
   Double_t        true_beam_startPy;
   Double_t        true_beam_startPz;
   Double_t        true_beam_startP;
   Double_t        true_beam_endPx;
   Double_t        true_beam_endPy;
   Double_t        true_beam_endPz;
   Double_t        true_beam_endP;
   Double_t        true_beam_startDirX;
   Double_t        true_beam_startDirY;
   Double_t        true_beam_startDirZ;
   Int_t           true_beam_nElasticScatters;
   vector<double>  *true_beam_elastic_costheta;
   vector<double>  *true_beam_elastic_X;
   vector<double>  *true_beam_elastic_Y;
   vector<double>  *true_beam_elastic_Z;
   Double_t        true_beam_IDE_totalDep;
   Bool_t          true_beam_IDE_found_in_recoVtx;
   Int_t           true_daughter_nPi0;
   Int_t           true_daughter_nPiPlus;
   Int_t           true_daughter_nProton;
   Int_t           true_daughter_nNeutron;
   Int_t           true_daughter_nPiMinus;
   Int_t           true_daughter_nNucleus;
   Int_t           reco_beam_vertex_slice;
   vector<vector<double> > *reco_beam_vertex_dRs;
   vector<int>     *reco_beam_vertex_hits_slices;
   vector<int>     *true_beam_daughter_PDG;
   vector<int>     *true_beam_daughter_ID;
   vector<double>  *true_beam_daughter_len;
   vector<double>  *true_beam_daughter_startX;
   vector<double>  *true_beam_daughter_startY;
   vector<double>  *true_beam_daughter_startZ;
   vector<double>  *true_beam_daughter_startPx;
   vector<double>  *true_beam_daughter_startPy;
   vector<double>  *true_beam_daughter_startPz;
   vector<double>  *true_beam_daughter_startP;
   vector<double>  *true_beam_daughter_endX;
   vector<double>  *true_beam_daughter_endY;
   vector<double>  *true_beam_daughter_endZ;
   vector<string>  *true_beam_daughter_Process;
   vector<int>     *true_beam_Pi0_decay_ID;
   vector<int>     *true_beam_Pi0_decay_PDG;
   vector<double>  *true_beam_Pi0_decay_startP;
   vector<int>     *true_beam_grand_daughter_ID;
   vector<int>     *true_beam_grand_daughter_parID;
   vector<int>     *true_beam_grand_daughter_PDG;
   string          *reco_beam_true_byE_endProcess;
   string          *reco_beam_true_byE_process;
   Int_t           reco_beam_true_byE_origin;
   Int_t           reco_beam_true_byE_PDG;
   Int_t           reco_beam_true_byE_ID;
   string          *reco_beam_true_byHits_endProcess;
   string          *reco_beam_true_byHits_process;
   Int_t           reco_beam_true_byHits_origin;
   Int_t           reco_beam_true_byHits_PDG;
   Int_t           reco_beam_true_byHits_ID;
   Bool_t          reco_beam_true_byE_matched;
   Bool_t          reco_beam_true_byHits_matched;
   Double_t        reco_beam_true_byHits_purity;
   vector<string>  *true_beam_processes;
   Bool_t          reco_daughter_true_byE_isPrimary;
   Double_t        data_BI_P;
   Double_t        data_BI_X;
   Double_t        data_BI_Y;
   Double_t        data_BI_Z;
   Int_t           data_BI_nFibersP1;
   Int_t           data_BI_nFibersP2;
   Int_t           data_BI_nFibersP3;
   vector<int>     *data_BI_PDG_candidates;
   Bool_t          quality_reco_view_0_hits_in_TPC5;
   Bool_t          quality_reco_view_1_hits_in_TPC5;
   Bool_t          quality_reco_view_2_hits_in_TPC5;
   Double_t        quality_reco_max_lateral;
   Double_t        quality_reco_max_segment;
   Double_t        quality_reco_view_0_max_segment;
   Double_t        quality_reco_view_1_max_segment;
   Double_t        quality_reco_view_2_max_segment;
   Double_t        quality_reco_view_0_wire_backtrack;
   Double_t        quality_reco_view_1_wire_backtrack;
   Double_t        quality_reco_view_2_wire_backtrack;
   vector<double>  *quality_reco_view_0_wire;
   vector<double>  *quality_reco_view_1_wire;
   vector<double>  *quality_reco_view_2_wire;
   vector<double>  *quality_reco_view_2_z;
   vector<double>  *quality_reco_view_0_tick;
   vector<double>  *quality_reco_view_1_tick;
   vector<double>  *quality_reco_view_2_tick;
   Double_t        reco_beam_Chi2_proton;
   Int_t           reco_beam_Chi2_ndof;
   vector<double>  *reco_daughter_Chi2_proton;
   vector<int>     *reco_daughter_Chi2_ndof;
   vector<double>  *reco_daughter_momByRange_proton;
   vector<double>  *reco_daughter_momByRange_muon;
   vector<double>  *reco_daughter_allTrack_momByRange_proton;
   vector<double>  *reco_daughter_allTrack_momByRange_muon;
   vector<double>  *reco_daughter_shower_Chi2_proton;
   vector<int>     *reco_daughter_shower_Chi2_ndof;
   vector<double>  *reco_daughter_trackScore;
   vector<double>  *reco_daughter_emScore;
   vector<double>  *reco_daughter_michelScore;
   vector<double>  *reco_daughter_shower_trackScore;
   vector<double>  *reco_daughter_shower_emScore;
   vector<double>  *reco_daughter_shower_michelScore;
   Double_t        reco_beam_true_byE_endPx;
   Double_t        reco_beam_true_byE_endPy;
   Double_t        reco_beam_true_byE_endPz;
   Double_t        reco_beam_true_byE_endE;
   Double_t        reco_beam_true_byE_endP;
   Double_t        reco_beam_true_byE_startPx;
   Double_t        reco_beam_true_byE_startPy;
   Double_t        reco_beam_true_byE_startPz;
   Double_t        reco_beam_true_byE_startE;
   Double_t        reco_beam_true_byE_startP;
   Double_t        reco_beam_true_byHits_endPx;
   Double_t        reco_beam_true_byHits_endPy;
   Double_t        reco_beam_true_byHits_endPz;
   Double_t        reco_beam_true_byHits_endE;
   Double_t        reco_beam_true_byHits_endP;
   Double_t        reco_beam_true_byHits_startPx;
   Double_t        reco_beam_true_byHits_startPy;
   Double_t        reco_beam_true_byHits_startPz;
   Double_t        reco_beam_true_byHits_startE;
   Double_t        reco_beam_true_byHits_startP;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_subrun;   //!
   TBranch        *b_event;   //!
   TBranch        *b_MC;   //!
   TBranch        *b_reco_beam_type;   //!
   TBranch        *b_reco_beam_startX;   //!
   TBranch        *b_reco_beam_startY;   //!
   TBranch        *b_reco_beam_startZ;   //!
   TBranch        *b_reco_beam_endX;   //!
   TBranch        *b_reco_beam_endY;   //!
   TBranch        *b_reco_beam_endZ;   //!
   TBranch        *b_reco_beam_len;   //!
   TBranch        *b_reco_beam_trackDirX;   //!
   TBranch        *b_reco_beam_trackDirY;   //!
   TBranch        *b_reco_beam_trackDirZ;   //!
   TBranch        *b_reco_beam_trackEndDirX;   //!
   TBranch        *b_reco_beam_trackEndDirY;   //!
   TBranch        *b_reco_beam_trackEndDirZ;   //!
   TBranch        *b_reco_beam_vtxX;   //!
   TBranch        *b_reco_beam_vtxY;   //!
   TBranch        *b_reco_beam_vtxZ;   //!
   TBranch        *b_reco_beam_trackID;   //!
   TBranch        *b_reco_beam_dQdX;   //!
   TBranch        *b_reco_beam_dEdX;   //!
   TBranch        *b_reco_beam_calibrated_dEdX;   //!
   TBranch        *b_reco_beam_resRange;   //!
   TBranch        *b_reco_beam_TrkPitch;   //!
   TBranch        *b_reco_beam_calo_wire;   //!
   TBranch        *b_reco_beam_calo_tick;   //!
   TBranch        *b_reco_beam_nTrackDaughters;   //!
   TBranch        *b_reco_beam_nShowerDaughters;   //!
   TBranch        *b_reco_beam_flipped;   //!
   TBranch        *b_reco_beam_passes_beam_cuts;   //!
   TBranch        *b_reco_daughter_trackID;   //!
   TBranch        *b_reco_daughter_true_byE_completeness;   //!
   TBranch        *b_reco_daughter_true_byE_purity;   //!
   TBranch        *b_reco_daughter_true_byE_PDG;   //!
   TBranch        *b_reco_daughter_true_byE_ID;   //!
   TBranch        *b_reco_daughter_true_byE_origin;   //!
   TBranch        *b_reco_daughter_true_byE_parID;   //!
   TBranch        *b_reco_daughter_true_byE_process;   //!
   TBranch        *b_reco_daughter_true_byHits_PDG;   //!
   TBranch        *b_reco_daughter_true_byHits_ID;   //!
   TBranch        *b_reco_daughter_true_byHits_origin;   //!
   TBranch        *b_reco_daughter_true_byHits_parID;   //!
   TBranch        *b_reco_daughter_true_byHits_process;   //!
   TBranch        *b_reco_daughter_true_byHits_purity;   //!
   TBranch        *b_reco_daughter_true_byHits_sharedHits;   //!
   TBranch        *b_reco_daughter_true_byHits_emHits;   //!
   TBranch        *b_reco_daughter_true_byHits_len;   //!
   TBranch        *b_reco_daughter_true_byHits_startX;   //!
   TBranch        *b_reco_daughter_true_byHits_startY;   //!
   TBranch        *b_reco_daughter_true_byHits_startZ;   //!
   TBranch        *b_reco_daughter_true_byHits_endX;   //!
   TBranch        *b_reco_daughter_true_byHits_endY;   //!
   TBranch        *b_reco_daughter_true_byHits_endZ;   //!
   TBranch        *b_reco_daughter_true_byHits_startPx;   //!
   TBranch        *b_reco_daughter_true_byHits_startPy;   //!
   TBranch        *b_reco_daughter_true_byHits_startPz;   //!
   TBranch        *b_reco_daughter_true_byHits_startP;   //!
   TBranch        *b_reco_daughter_true_byHits_startE;   //!
   TBranch        *b_reco_daughter_PFP_true_byHits_PDG;   //!
   TBranch        *b_reco_daughter_PFP_true_byHits_ID;   //!
   TBranch        *b_reco_daughter_PFP_true_byHits_origin;   //!
   TBranch        *b_reco_daughter_PFP_true_byHits_parID;   //!
   TBranch        *b_reco_daughter_PFP_true_byHits_process;   //!
   TBranch        *b_reco_daughter_PFP_true_byHits_sharedHits;   //!
   TBranch        *b_reco_daughter_PFP_true_byHits_emHits;   //!
   TBranch        *b_reco_daughter_PFP_true_byHits_len;   //!
   TBranch        *b_reco_daughter_PFP_true_byHits_startX;   //!
   TBranch        *b_reco_daughter_PFP_true_byHits_startY;   //!
   TBranch        *b_reco_daughter_PFP_true_byHits_startZ;   //!
   TBranch        *b_reco_daughter_PFP_true_byHits_endX;   //!
   TBranch        *b_reco_daughter_PFP_true_byHits_endY;   //!
   TBranch        *b_reco_daughter_PFP_true_byHits_endZ;   //!
   TBranch        *b_reco_daughter_PFP_true_byHits_startPx;   //!
   TBranch        *b_reco_daughter_PFP_true_byHits_startPy;   //!
   TBranch        *b_reco_daughter_PFP_true_byHits_startPz;   //!
   TBranch        *b_reco_daughter_PFP_true_byHits_startP;   //!
   TBranch        *b_reco_daughter_PFP_true_byHits_startE;   //!
   TBranch        *b_reco_daughter_PFP_true_byHits_endProcess;   //!
   TBranch        *b_reco_daughter_allTrack_ID;   //!
   TBranch        *b_reco_daughter_allTrack_dEdX;   //!
   TBranch        *b_reco_daughter_allTrack_dQdX;   //!
   TBranch        *b_reco_daughter_allTrack_resRange;   //!
   TBranch        *b_reco_daughter_allTrack_dQdX_SCE;   //!
   TBranch        *b_reco_daughter_allTrack_dEdX_SCE;   //!
   TBranch        *b_reco_daughter_allTrack_resRange_SCE;   //!
   TBranch        *b_reco_daughter_allTrack_calibrated_dEdX;   //!
   TBranch        *b_reco_daughter_allTrack_calibrated_dEdX_SCE;   //!
   TBranch        *b_reco_daughter_allTrack_Chi2_proton;   //!
   TBranch        *b_reco_daughter_allTrack_Chi2_ndof;   //!
   TBranch        *b_reco_daughter_allTrack_startX;   //!
   TBranch        *b_reco_daughter_allTrack_startY;   //!
   TBranch        *b_reco_daughter_allTrack_startZ;   //!
   TBranch        *b_reco_daughter_allTrack_endX;   //!
   TBranch        *b_reco_daughter_allTrack_endY;   //!
   TBranch        *b_reco_daughter_allTrack_endZ;   //!
   TBranch        *b_reco_daughter_allTrack_dR;   //!
   TBranch        *b_reco_daughter_allTrack_to_vertex;   //!
   TBranch        *b_reco_daughter_shower_true_byE_PDG;   //!
   TBranch        *b_reco_daughter_shower_true_byE_ID;   //!
   TBranch        *b_reco_daughter_shower_true_byE_origin;   //!
   TBranch        *b_reco_daughter_shower_true_byE_parID;   //!
   TBranch        *b_reco_daughter_shower_true_byE_startPx;   //!
   TBranch        *b_reco_daughter_shower_true_byE_startPy;   //!
   TBranch        *b_reco_daughter_shower_true_byE_startPz;   //!
   TBranch        *b_reco_daughter_shower_true_byE_startP;   //!
   TBranch        *b_reco_daughter_shower_true_byHits_PDG;   //!
   TBranch        *b_reco_daughter_shower_true_byHits_ID;   //!
   TBranch        *b_reco_daughter_shower_true_byHits_origin;   //!
   TBranch        *b_reco_daughter_shower_true_byHits_parID;   //!
   TBranch        *b_reco_daughter_shower_true_byHits_process;   //!
   TBranch        *b_reco_daughter_shower_true_byHits_purity;   //!
   TBranch        *b_reco_daughter_shower_true_byHits_startPx;   //!
   TBranch        *b_reco_daughter_shower_true_byHits_startPy;   //!
   TBranch        *b_reco_daughter_shower_true_byHits_startPz;   //!
   TBranch        *b_reco_daughter_shower_true_byHits_startP;   //!
   TBranch        *b_reco_daughter_showerID;   //!
   TBranch        *b_reco_daughter_dQdX;   //!
   TBranch        *b_reco_daughter_dEdX;   //!
   TBranch        *b_reco_daughter_resRange;   //!
   TBranch        *b_reco_daughter_shower_dQdX;   //!
   TBranch        *b_reco_daughter_shower_dEdX;   //!
   TBranch        *b_reco_daughter_shower_resRange;   //!
   TBranch        *b_reco_daughter_len;   //!
   TBranch        *b_reco_daughter_startX;   //!
   TBranch        *b_reco_daughter_startY;   //!
   TBranch        *b_reco_daughter_startZ;   //!
   TBranch        *b_reco_daughter_endX;   //!
   TBranch        *b_reco_daughter_endY;   //!
   TBranch        *b_reco_daughter_endZ;   //!
   TBranch        *b_reco_daughter_deltaR;   //!
   TBranch        *b_reco_daughter_dR;   //!
   TBranch        *b_reco_daughter_to_vertex;   //!
   TBranch        *b_reco_daughter_slice;   //!
   TBranch        *b_reco_daughter_shower_to_vertex;   //!
   TBranch        *b_reco_daughter_shower_startX;   //!
   TBranch        *b_reco_daughter_shower_startY;   //!
   TBranch        *b_reco_daughter_shower_startZ;   //!
   TBranch        *b_reco_daughter_shower_len;   //!
   TBranch        *b_reco_daughter_PFP_ID;   //!
   TBranch        *b_reco_daughter_PFP_trackScore;   //!
   TBranch        *b_reco_daughter_PFP_emScore;   //!
   TBranch        *b_reco_daughter_PFP_michelScore;   //!
   TBranch        *b_true_beam_PDG;   //!
   TBranch        *b_true_beam_ID;   //!
   TBranch        *b_true_beam_endProcess;   //!
   TBranch        *b_true_beam_endX;   //!
   TBranch        *b_true_beam_endY;   //!
   TBranch        *b_true_beam_endZ;   //!
   TBranch        *b_true_beam_startX;   //!
   TBranch        *b_true_beam_startY;   //!
   TBranch        *b_true_beam_startZ;   //!
   TBranch        *b_true_beam_startPx;   //!
   TBranch        *b_true_beam_startPy;   //!
   TBranch        *b_true_beam_startPz;   //!
   TBranch        *b_true_beam_startP;   //!
   TBranch        *b_true_beam_endPx;   //!
   TBranch        *b_true_beam_endPy;   //!
   TBranch        *b_true_beam_endPz;   //!
   TBranch        *b_true_beam_endP;   //!
   TBranch        *b_true_beam_startDirX;   //!
   TBranch        *b_true_beam_startDirY;   //!
   TBranch        *b_true_beam_startDirZ;   //!
   TBranch        *b_true_beam_nElasticScatters;   //!
   TBranch        *b_true_beam_elastic_costheta;   //!
   TBranch        *b_true_beam_elastic_X;   //!
   TBranch        *b_true_beam_elastic_Y;   //!
   TBranch        *b_true_beam_elastic_Z;   //!
   TBranch        *b_true_beam_IDE_totalDep;   //!
   TBranch        *b_true_beam_IDE_found_in_recoVtx;   //!
   TBranch        *b_true_daughter_nPi0;   //!
   TBranch        *b_true_daughter_nPiPlus;   //!
   TBranch        *b_true_daughter_nProton;   //!
   TBranch        *b_true_daughter_nNeutron;   //!
   TBranch        *b_true_daughter_nPiMinus;   //!
   TBranch        *b_true_daughter_nNucleus;   //!
   TBranch        *b_reco_beam_vertex_slice;   //!
   TBranch        *b_reco_beam_vertex_dRs;   //!
   TBranch        *b_reco_beam_vertex_hits_slices;   //!
   TBranch        *b_true_beam_daughter_PDG;   //!
   TBranch        *b_true_beam_daughter_ID;   //!
   TBranch        *b_true_beam_daughter_len;   //!
   TBranch        *b_true_beam_daughter_startX;   //!
   TBranch        *b_true_beam_daughter_startY;   //!
   TBranch        *b_true_beam_daughter_startZ;   //!
   TBranch        *b_true_beam_daughter_startPx;   //!
   TBranch        *b_true_beam_daughter_startPy;   //!
   TBranch        *b_true_beam_daughter_startPz;   //!
   TBranch        *b_true_beam_daughter_startP;   //!
   TBranch        *b_true_beam_daughter_endX;   //!
   TBranch        *b_true_beam_daughter_endY;   //!
   TBranch        *b_true_beam_daughter_endZ;   //!
   TBranch        *b_true_beam_daughter_Process;   //!
   TBranch        *b_true_beam_Pi0_decay_ID;   //!
   TBranch        *b_true_beam_Pi0_decay_PDG;   //!
   TBranch        *b_true_beam_Pi0_decay_startP;   //!
   TBranch        *b_true_beam_grand_daughter_ID;   //!
   TBranch        *b_true_beam_grand_daughter_parID;   //!
   TBranch        *b_true_beam_grand_daughter_PDG;   //!
   TBranch        *b_reco_beam_true_byE_endProcess;   //!
   TBranch        *b_reco_beam_true_byE_process;   //!
   TBranch        *b_reco_beam_true_byE_origin;   //!
   TBranch        *b_reco_beam_true_byE_PDG;   //!
   TBranch        *b_reco_beam_true_byE_ID;   //!
   TBranch        *b_reco_beam_true_byHits_endProcess;   //!
   TBranch        *b_reco_beam_true_byHits_process;   //!
   TBranch        *b_reco_beam_true_byHits_origin;   //!
   TBranch        *b_reco_beam_true_byHits_PDG;   //!
   TBranch        *b_reco_beam_true_byHits_ID;   //!
   TBranch        *b_reco_beam_true_byE_matched;   //!
   TBranch        *b_reco_beam_true_byHits_matched;   //!
   TBranch        *b_reco_beam_true_byHits_purity;   //!
   TBranch        *b_true_beam_processes;   //!
   TBranch        *b_reco_daughter_true_byE_isPrimary;   //!
   TBranch        *b_data_BI_P;   //!
   TBranch        *b_data_BI_X;   //!
   TBranch        *b_data_BI_Y;   //!
   TBranch        *b_data_BI_Z;   //!
   TBranch        *b_data_BI_nFibersP1;   //!
   TBranch        *b_data_BI_nFibersP2;   //!
   TBranch        *b_data_BI_nFibersP3;   //!
   TBranch        *b_data_BI_PDG_candidates;   //!
   TBranch        *b_quality_reco_view_0_hits_in_TPC5;   //!
   TBranch        *b_quality_reco_view_1_hits_in_TPC5;   //!
   TBranch        *b_quality_reco_view_2_hits_in_TPC5;   //!
   TBranch        *b_quality_reco_max_lateral;   //!
   TBranch        *b_quality_reco_max_segment;   //!
   TBranch        *b_quality_reco_view_0_max_segment;   //!
   TBranch        *b_quality_reco_view_1_max_segment;   //!
   TBranch        *b_quality_reco_view_2_max_segment;   //!
   TBranch        *b_quality_reco_view_0_wire_backtrack;   //!
   TBranch        *b_quality_reco_view_1_wire_backtrack;   //!
   TBranch        *b_quality_reco_view_2_wire_backtrack;   //!
   TBranch        *b_quality_reco_view_0_wire;   //!
   TBranch        *b_quality_reco_view_1_wire;   //!
   TBranch        *b_quality_reco_view_2_wire;   //!
   TBranch        *b_quality_reco_view_2_z;   //!
   TBranch        *b_quality_reco_view_0_tick;   //!
   TBranch        *b_quality_reco_view_1_tick;   //!
   TBranch        *b_quality_reco_view_2_tick;   //!
   TBranch        *b_reco_beam_Chi2_proton;   //!
   TBranch        *b_reco_beam_Chi2_ndof;   //!
   TBranch        *b_reco_daughter_Chi2_proton;   //!
   TBranch        *b_reco_daughter_Chi2_ndof;   //!
   TBranch        *b_reco_daughter_momByRange_proton;   //!
   TBranch        *b_reco_daughter_momByRange_muon;   //!
   TBranch        *b_reco_daughter_allTrack_momByRange_proton;   //!
   TBranch        *b_reco_daughter_allTrack_momByRange_muon;   //!
   TBranch        *b_reco_daughter_shower_Chi2_proton;   //!
   TBranch        *b_reco_daughter_shower_Chi2_ndof;   //!
   TBranch        *b_reco_daughter_trackScore;   //!
   TBranch        *b_reco_daughter_emScore;   //!
   TBranch        *b_reco_daughter_michelScore;   //!
   TBranch        *b_reco_daughter_shower_trackScore;   //!
   TBranch        *b_reco_daughter_shower_emScore;   //!
   TBranch        *b_reco_daughter_shower_michelScore;   //!
   TBranch        *b_reco_beam_true_byE_endPx;   //!
   TBranch        *b_reco_beam_true_byE_endPy;   //!
   TBranch        *b_reco_beam_true_byE_endPz;   //!
   TBranch        *b_reco_beam_true_byE_endE;   //!
   TBranch        *b_reco_beam_true_byE_endP;   //!
   TBranch        *b_reco_beam_true_byE_startPx;   //!
   TBranch        *b_reco_beam_true_byE_startPy;   //!
   TBranch        *b_reco_beam_true_byE_startPz;   //!
   TBranch        *b_reco_beam_true_byE_startE;   //!
   TBranch        *b_reco_beam_true_byE_startP;   //!
   TBranch        *b_reco_beam_true_byHits_endPx;   //!
   TBranch        *b_reco_beam_true_byHits_endPy;   //!
   TBranch        *b_reco_beam_true_byHits_endPz;   //!
   TBranch        *b_reco_beam_true_byHits_endE;   //!
   TBranch        *b_reco_beam_true_byHits_endP;   //!
   TBranch        *b_reco_beam_true_byHits_startPx;   //!
   TBranch        *b_reco_beam_true_byHits_startPy;   //!
   TBranch        *b_reco_beam_true_byHits_startPz;   //!
   TBranch        *b_reco_beam_true_byHits_startE;   //!
   TBranch        *b_reco_beam_true_byHits_startP;   //!

   PDEventPro(TTree *tree=0);
   virtual ~PDEventPro();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef PDEventPro_cxx
PDEventPro::PDEventPro(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/pnfs/dune/persistent/users/jiangl/pionana_1128/pionana_mc0to30.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/pnfs/dune/persistent/users/jiangl/pionana_1128/pionana_mc0to30.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/pnfs/dune/persistent/users/jiangl/pionana_1128/pionana_mc0to30.root:/pionana");
      dir->GetObject("beamana",tree);

   }
   Init(tree);
}

PDEventPro::~PDEventPro()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PDEventPro::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PDEventPro::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void PDEventPro::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   reco_beam_dQdX = 0;
   reco_beam_dEdX = 0;
   reco_beam_calibrated_dEdX = 0;
   reco_beam_resRange = 0;
   reco_beam_TrkPitch = 0;
   reco_beam_calo_wire = 0;
   reco_beam_calo_tick = 0;
   reco_daughter_trackID = 0;
   reco_daughter_true_byE_completeness = 0;
   reco_daughter_true_byE_purity = 0;
   reco_daughter_true_byE_PDG = 0;
   reco_daughter_true_byE_ID = 0;
   reco_daughter_true_byE_origin = 0;
   reco_daughter_true_byE_parID = 0;
   reco_daughter_true_byE_process = 0;
   reco_daughter_true_byHits_PDG = 0;
   reco_daughter_true_byHits_ID = 0;
   reco_daughter_true_byHits_origin = 0;
   reco_daughter_true_byHits_parID = 0;
   reco_daughter_true_byHits_process = 0;
   reco_daughter_true_byHits_purity = 0;
   reco_daughter_true_byHits_sharedHits = 0;
   reco_daughter_true_byHits_emHits = 0;
   reco_daughter_true_byHits_len = 0;
   reco_daughter_true_byHits_startX = 0;
   reco_daughter_true_byHits_startY = 0;
   reco_daughter_true_byHits_startZ = 0;
   reco_daughter_true_byHits_endX = 0;
   reco_daughter_true_byHits_endY = 0;
   reco_daughter_true_byHits_endZ = 0;
   reco_daughter_true_byHits_startPx = 0;
   reco_daughter_true_byHits_startPy = 0;
   reco_daughter_true_byHits_startPz = 0;
   reco_daughter_true_byHits_startP = 0;
   reco_daughter_true_byHits_startE = 0;
   reco_daughter_PFP_true_byHits_PDG = 0;
   reco_daughter_PFP_true_byHits_ID = 0;
   reco_daughter_PFP_true_byHits_origin = 0;
   reco_daughter_PFP_true_byHits_parID = 0;
   reco_daughter_PFP_true_byHits_process = 0;
   reco_daughter_PFP_true_byHits_sharedHits = 0;
   reco_daughter_PFP_true_byHits_emHits = 0;
   reco_daughter_PFP_true_byHits_len = 0;
   reco_daughter_PFP_true_byHits_startX = 0;
   reco_daughter_PFP_true_byHits_startY = 0;
   reco_daughter_PFP_true_byHits_startZ = 0;
   reco_daughter_PFP_true_byHits_endX = 0;
   reco_daughter_PFP_true_byHits_endY = 0;
   reco_daughter_PFP_true_byHits_endZ = 0;
   reco_daughter_PFP_true_byHits_startPx = 0;
   reco_daughter_PFP_true_byHits_startPy = 0;
   reco_daughter_PFP_true_byHits_startPz = 0;
   reco_daughter_PFP_true_byHits_startP = 0;
   reco_daughter_PFP_true_byHits_startE = 0;
   reco_daughter_PFP_true_byHits_endProcess = 0;
   reco_daughter_allTrack_ID = 0;
   reco_daughter_allTrack_dEdX = 0;
   reco_daughter_allTrack_dQdX = 0;
   reco_daughter_allTrack_resRange = 0;
   reco_daughter_allTrack_dQdX_SCE = 0;
   reco_daughter_allTrack_dEdX_SCE = 0;
   reco_daughter_allTrack_resRange_SCE = 0;
   reco_daughter_allTrack_calibrated_dEdX = 0;
   reco_daughter_allTrack_calibrated_dEdX_SCE = 0;
   reco_daughter_allTrack_Chi2_proton = 0;
   reco_daughter_allTrack_Chi2_ndof = 0;
   reco_daughter_allTrack_startX = 0;
   reco_daughter_allTrack_startY = 0;
   reco_daughter_allTrack_startZ = 0;
   reco_daughter_allTrack_endX = 0;
   reco_daughter_allTrack_endY = 0;
   reco_daughter_allTrack_endZ = 0;
   reco_daughter_allTrack_dR = 0;
   reco_daughter_allTrack_to_vertex = 0;
   reco_daughter_shower_true_byE_PDG = 0;
   reco_daughter_shower_true_byE_ID = 0;
   reco_daughter_shower_true_byE_origin = 0;
   reco_daughter_shower_true_byE_parID = 0;
   reco_daughter_shower_true_byE_startPx = 0;
   reco_daughter_shower_true_byE_startPy = 0;
   reco_daughter_shower_true_byE_startPz = 0;
   reco_daughter_shower_true_byE_startP = 0;
   reco_daughter_shower_true_byHits_PDG = 0;
   reco_daughter_shower_true_byHits_ID = 0;
   reco_daughter_shower_true_byHits_origin = 0;
   reco_daughter_shower_true_byHits_parID = 0;
   reco_daughter_shower_true_byHits_process = 0;
   reco_daughter_shower_true_byHits_purity = 0;
   reco_daughter_shower_true_byHits_startPx = 0;
   reco_daughter_shower_true_byHits_startPy = 0;
   reco_daughter_shower_true_byHits_startPz = 0;
   reco_daughter_shower_true_byHits_startP = 0;
   reco_daughter_showerID = 0;
   reco_daughter_dQdX = 0;
   reco_daughter_dEdX = 0;
   reco_daughter_resRange = 0;
   reco_daughter_shower_dQdX = 0;
   reco_daughter_shower_dEdX = 0;
   reco_daughter_shower_resRange = 0;
   reco_daughter_len = 0;
   reco_daughter_startX = 0;
   reco_daughter_startY = 0;
   reco_daughter_startZ = 0;
   reco_daughter_endX = 0;
   reco_daughter_endY = 0;
   reco_daughter_endZ = 0;
   reco_daughter_deltaR = 0;
   reco_daughter_dR = 0;
   reco_daughter_to_vertex = 0;
   reco_daughter_slice = 0;
   reco_daughter_shower_to_vertex = 0;
   reco_daughter_shower_startX = 0;
   reco_daughter_shower_startY = 0;
   reco_daughter_shower_startZ = 0;
   reco_daughter_shower_len = 0;
   reco_daughter_PFP_ID = 0;
   reco_daughter_PFP_trackScore = 0;
   reco_daughter_PFP_emScore = 0;
   reco_daughter_PFP_michelScore = 0;
   true_beam_endProcess = 0;
   true_beam_elastic_costheta = 0;
   true_beam_elastic_X = 0;
   true_beam_elastic_Y = 0;
   true_beam_elastic_Z = 0;
   reco_beam_vertex_dRs = 0;
   reco_beam_vertex_hits_slices = 0;
   true_beam_daughter_PDG = 0;
   true_beam_daughter_ID = 0;
   true_beam_daughter_len = 0;
   true_beam_daughter_startX = 0;
   true_beam_daughter_startY = 0;
   true_beam_daughter_startZ = 0;
   true_beam_daughter_startPx = 0;
   true_beam_daughter_startPy = 0;
   true_beam_daughter_startPz = 0;
   true_beam_daughter_startP = 0;
   true_beam_daughter_endX = 0;
   true_beam_daughter_endY = 0;
   true_beam_daughter_endZ = 0;
   true_beam_daughter_Process = 0;
   true_beam_Pi0_decay_ID = 0;
   true_beam_Pi0_decay_PDG = 0;
   true_beam_Pi0_decay_startP = 0;
   true_beam_grand_daughter_ID = 0;
   true_beam_grand_daughter_parID = 0;
   true_beam_grand_daughter_PDG = 0;
   reco_beam_true_byE_endProcess = 0;
   reco_beam_true_byE_process = 0;
   reco_beam_true_byHits_endProcess = 0;
   reco_beam_true_byHits_process = 0;
   true_beam_processes = 0;
   data_BI_PDG_candidates = 0;
   quality_reco_view_0_wire = 0;
   quality_reco_view_1_wire = 0;
   quality_reco_view_2_wire = 0;
   quality_reco_view_2_z = 0;
   quality_reco_view_0_tick = 0;
   quality_reco_view_1_tick = 0;
   quality_reco_view_2_tick = 0;
   reco_daughter_Chi2_proton = 0;
   reco_daughter_Chi2_ndof = 0;
   reco_daughter_momByRange_proton = 0;
   reco_daughter_momByRange_muon = 0;
   reco_daughter_allTrack_momByRange_proton = 0;
   reco_daughter_allTrack_momByRange_muon = 0;
   reco_daughter_shower_Chi2_proton = 0;
   reco_daughter_shower_Chi2_ndof = 0;
   reco_daughter_trackScore = 0;
   reco_daughter_emScore = 0;
   reco_daughter_michelScore = 0;
   reco_daughter_shower_trackScore = 0;
   reco_daughter_shower_emScore = 0;
   reco_daughter_shower_michelScore = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("subrun", &subrun, &b_subrun);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("MC", &MC, &b_MC);
   fChain->SetBranchAddress("reco_beam_type", &reco_beam_type, &b_reco_beam_type);
   fChain->SetBranchAddress("reco_beam_startX", &reco_beam_startX, &b_reco_beam_startX);
   fChain->SetBranchAddress("reco_beam_startY", &reco_beam_startY, &b_reco_beam_startY);
   fChain->SetBranchAddress("reco_beam_startZ", &reco_beam_startZ, &b_reco_beam_startZ);
   fChain->SetBranchAddress("reco_beam_endX", &reco_beam_endX, &b_reco_beam_endX);
   fChain->SetBranchAddress("reco_beam_endY", &reco_beam_endY, &b_reco_beam_endY);
   fChain->SetBranchAddress("reco_beam_endZ", &reco_beam_endZ, &b_reco_beam_endZ);
   fChain->SetBranchAddress("reco_beam_len", &reco_beam_len, &b_reco_beam_len);
   fChain->SetBranchAddress("reco_beam_trackDirX", &reco_beam_trackDirX, &b_reco_beam_trackDirX);
   fChain->SetBranchAddress("reco_beam_trackDirY", &reco_beam_trackDirY, &b_reco_beam_trackDirY);
   fChain->SetBranchAddress("reco_beam_trackDirZ", &reco_beam_trackDirZ, &b_reco_beam_trackDirZ);
   fChain->SetBranchAddress("reco_beam_trackEndDirX", &reco_beam_trackEndDirX, &b_reco_beam_trackEndDirX);
   fChain->SetBranchAddress("reco_beam_trackEndDirY", &reco_beam_trackEndDirY, &b_reco_beam_trackEndDirY);
   fChain->SetBranchAddress("reco_beam_trackEndDirZ", &reco_beam_trackEndDirZ, &b_reco_beam_trackEndDirZ);
   fChain->SetBranchAddress("reco_beam_vtxX", &reco_beam_vtxX, &b_reco_beam_vtxX);
   fChain->SetBranchAddress("reco_beam_vtxY", &reco_beam_vtxY, &b_reco_beam_vtxY);
   fChain->SetBranchAddress("reco_beam_vtxZ", &reco_beam_vtxZ, &b_reco_beam_vtxZ);
   fChain->SetBranchAddress("reco_beam_trackID", &reco_beam_trackID, &b_reco_beam_trackID);
   fChain->SetBranchAddress("reco_beam_dQdX", &reco_beam_dQdX, &b_reco_beam_dQdX);
   fChain->SetBranchAddress("reco_beam_dEdX", &reco_beam_dEdX, &b_reco_beam_dEdX);
   fChain->SetBranchAddress("reco_beam_calibrated_dEdX", &reco_beam_calibrated_dEdX, &b_reco_beam_calibrated_dEdX);
   fChain->SetBranchAddress("reco_beam_resRange", &reco_beam_resRange, &b_reco_beam_resRange);
   fChain->SetBranchAddress("reco_beam_TrkPitch", &reco_beam_TrkPitch, &b_reco_beam_TrkPitch);
   fChain->SetBranchAddress("reco_beam_calo_wire", &reco_beam_calo_wire, &b_reco_beam_calo_wire);
   fChain->SetBranchAddress("reco_beam_calo_tick", &reco_beam_calo_tick, &b_reco_beam_calo_tick);
   fChain->SetBranchAddress("reco_beam_nTrackDaughters", &reco_beam_nTrackDaughters, &b_reco_beam_nTrackDaughters);
   fChain->SetBranchAddress("reco_beam_nShowerDaughters", &reco_beam_nShowerDaughters, &b_reco_beam_nShowerDaughters);
   fChain->SetBranchAddress("reco_beam_flipped", &reco_beam_flipped, &b_reco_beam_flipped);
   fChain->SetBranchAddress("reco_beam_passes_beam_cuts", &reco_beam_passes_beam_cuts, &b_reco_beam_passes_beam_cuts);
   fChain->SetBranchAddress("reco_daughter_trackID", &reco_daughter_trackID, &b_reco_daughter_trackID);
   fChain->SetBranchAddress("reco_daughter_true_byE_completeness", &reco_daughter_true_byE_completeness, &b_reco_daughter_true_byE_completeness);
   fChain->SetBranchAddress("reco_daughter_true_byE_purity", &reco_daughter_true_byE_purity, &b_reco_daughter_true_byE_purity);
   fChain->SetBranchAddress("reco_daughter_true_byE_PDG", &reco_daughter_true_byE_PDG, &b_reco_daughter_true_byE_PDG);
   fChain->SetBranchAddress("reco_daughter_true_byE_ID", &reco_daughter_true_byE_ID, &b_reco_daughter_true_byE_ID);
   fChain->SetBranchAddress("reco_daughter_true_byE_origin", &reco_daughter_true_byE_origin, &b_reco_daughter_true_byE_origin);
   fChain->SetBranchAddress("reco_daughter_true_byE_parID", &reco_daughter_true_byE_parID, &b_reco_daughter_true_byE_parID);
   fChain->SetBranchAddress("reco_daughter_true_byE_process", &reco_daughter_true_byE_process, &b_reco_daughter_true_byE_process);
   fChain->SetBranchAddress("reco_daughter_true_byHits_PDG", &reco_daughter_true_byHits_PDG, &b_reco_daughter_true_byHits_PDG);
   fChain->SetBranchAddress("reco_daughter_true_byHits_ID", &reco_daughter_true_byHits_ID, &b_reco_daughter_true_byHits_ID);
   fChain->SetBranchAddress("reco_daughter_true_byHits_origin", &reco_daughter_true_byHits_origin, &b_reco_daughter_true_byHits_origin);
   fChain->SetBranchAddress("reco_daughter_true_byHits_parID", &reco_daughter_true_byHits_parID, &b_reco_daughter_true_byHits_parID);
   fChain->SetBranchAddress("reco_daughter_true_byHits_process", &reco_daughter_true_byHits_process, &b_reco_daughter_true_byHits_process);
   fChain->SetBranchAddress("reco_daughter_true_byHits_purity", &reco_daughter_true_byHits_purity, &b_reco_daughter_true_byHits_purity);
   fChain->SetBranchAddress("reco_daughter_true_byHits_sharedHits", &reco_daughter_true_byHits_sharedHits, &b_reco_daughter_true_byHits_sharedHits);
   fChain->SetBranchAddress("reco_daughter_true_byHits_emHits", &reco_daughter_true_byHits_emHits, &b_reco_daughter_true_byHits_emHits);
   fChain->SetBranchAddress("reco_daughter_true_byHits_len", &reco_daughter_true_byHits_len, &b_reco_daughter_true_byHits_len);
   fChain->SetBranchAddress("reco_daughter_true_byHits_startX", &reco_daughter_true_byHits_startX, &b_reco_daughter_true_byHits_startX);
   fChain->SetBranchAddress("reco_daughter_true_byHits_startY", &reco_daughter_true_byHits_startY, &b_reco_daughter_true_byHits_startY);
   fChain->SetBranchAddress("reco_daughter_true_byHits_startZ", &reco_daughter_true_byHits_startZ, &b_reco_daughter_true_byHits_startZ);
   fChain->SetBranchAddress("reco_daughter_true_byHits_endX", &reco_daughter_true_byHits_endX, &b_reco_daughter_true_byHits_endX);
   fChain->SetBranchAddress("reco_daughter_true_byHits_endY", &reco_daughter_true_byHits_endY, &b_reco_daughter_true_byHits_endY);
   fChain->SetBranchAddress("reco_daughter_true_byHits_endZ", &reco_daughter_true_byHits_endZ, &b_reco_daughter_true_byHits_endZ);
   fChain->SetBranchAddress("reco_daughter_true_byHits_startPx", &reco_daughter_true_byHits_startPx, &b_reco_daughter_true_byHits_startPx);
   fChain->SetBranchAddress("reco_daughter_true_byHits_startPy", &reco_daughter_true_byHits_startPy, &b_reco_daughter_true_byHits_startPy);
   fChain->SetBranchAddress("reco_daughter_true_byHits_startPz", &reco_daughter_true_byHits_startPz, &b_reco_daughter_true_byHits_startPz);
   fChain->SetBranchAddress("reco_daughter_true_byHits_startP", &reco_daughter_true_byHits_startP, &b_reco_daughter_true_byHits_startP);
   fChain->SetBranchAddress("reco_daughter_true_byHits_startE", &reco_daughter_true_byHits_startE, &b_reco_daughter_true_byHits_startE);
   fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_PDG", &reco_daughter_PFP_true_byHits_PDG, &b_reco_daughter_PFP_true_byHits_PDG);
   fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_ID", &reco_daughter_PFP_true_byHits_ID, &b_reco_daughter_PFP_true_byHits_ID);
   fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_origin", &reco_daughter_PFP_true_byHits_origin, &b_reco_daughter_PFP_true_byHits_origin);
   fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_parID", &reco_daughter_PFP_true_byHits_parID, &b_reco_daughter_PFP_true_byHits_parID);
   fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_process", &reco_daughter_PFP_true_byHits_process, &b_reco_daughter_PFP_true_byHits_process);
   fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_sharedHits", &reco_daughter_PFP_true_byHits_sharedHits, &b_reco_daughter_PFP_true_byHits_sharedHits);
   fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_emHits", &reco_daughter_PFP_true_byHits_emHits, &b_reco_daughter_PFP_true_byHits_emHits);
   fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_len", &reco_daughter_PFP_true_byHits_len, &b_reco_daughter_PFP_true_byHits_len);
   fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_startX", &reco_daughter_PFP_true_byHits_startX, &b_reco_daughter_PFP_true_byHits_startX);
   fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_startY", &reco_daughter_PFP_true_byHits_startY, &b_reco_daughter_PFP_true_byHits_startY);
   fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_startZ", &reco_daughter_PFP_true_byHits_startZ, &b_reco_daughter_PFP_true_byHits_startZ);
   fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_endX", &reco_daughter_PFP_true_byHits_endX, &b_reco_daughter_PFP_true_byHits_endX);
   fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_endY", &reco_daughter_PFP_true_byHits_endY, &b_reco_daughter_PFP_true_byHits_endY);
   fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_endZ", &reco_daughter_PFP_true_byHits_endZ, &b_reco_daughter_PFP_true_byHits_endZ);
   fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_startPx", &reco_daughter_PFP_true_byHits_startPx, &b_reco_daughter_PFP_true_byHits_startPx);
   fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_startPy", &reco_daughter_PFP_true_byHits_startPy, &b_reco_daughter_PFP_true_byHits_startPy);
   fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_startPz", &reco_daughter_PFP_true_byHits_startPz, &b_reco_daughter_PFP_true_byHits_startPz);
   fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_startP", &reco_daughter_PFP_true_byHits_startP, &b_reco_daughter_PFP_true_byHits_startP);
   fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_startE", &reco_daughter_PFP_true_byHits_startE, &b_reco_daughter_PFP_true_byHits_startE);
   fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_endProcess", &reco_daughter_PFP_true_byHits_endProcess, &b_reco_daughter_PFP_true_byHits_endProcess);
   fChain->SetBranchAddress("reco_daughter_allTrack_ID", &reco_daughter_allTrack_ID, &b_reco_daughter_allTrack_ID);
   fChain->SetBranchAddress("reco_daughter_allTrack_dEdX", &reco_daughter_allTrack_dEdX, &b_reco_daughter_allTrack_dEdX);
   fChain->SetBranchAddress("reco_daughter_allTrack_dQdX", &reco_daughter_allTrack_dQdX, &b_reco_daughter_allTrack_dQdX);
   fChain->SetBranchAddress("reco_daughter_allTrack_resRange", &reco_daughter_allTrack_resRange, &b_reco_daughter_allTrack_resRange);
   fChain->SetBranchAddress("reco_daughter_allTrack_dQdX_SCE", &reco_daughter_allTrack_dQdX_SCE, &b_reco_daughter_allTrack_dQdX_SCE);
   fChain->SetBranchAddress("reco_daughter_allTrack_dEdX_SCE", &reco_daughter_allTrack_dEdX_SCE, &b_reco_daughter_allTrack_dEdX_SCE);
   fChain->SetBranchAddress("reco_daughter_allTrack_resRange_SCE", &reco_daughter_allTrack_resRange_SCE, &b_reco_daughter_allTrack_resRange_SCE);
   fChain->SetBranchAddress("reco_daughter_allTrack_calibrated_dEdX", &reco_daughter_allTrack_calibrated_dEdX, &b_reco_daughter_allTrack_calibrated_dEdX);
   fChain->SetBranchAddress("reco_daughter_allTrack_calibrated_dEdX_SCE", &reco_daughter_allTrack_calibrated_dEdX_SCE, &b_reco_daughter_allTrack_calibrated_dEdX_SCE);
   fChain->SetBranchAddress("reco_daughter_allTrack_Chi2_proton", &reco_daughter_allTrack_Chi2_proton, &b_reco_daughter_allTrack_Chi2_proton);
   fChain->SetBranchAddress("reco_daughter_allTrack_Chi2_ndof", &reco_daughter_allTrack_Chi2_ndof, &b_reco_daughter_allTrack_Chi2_ndof);
   fChain->SetBranchAddress("reco_daughter_allTrack_startX", &reco_daughter_allTrack_startX, &b_reco_daughter_allTrack_startX);
   fChain->SetBranchAddress("reco_daughter_allTrack_startY", &reco_daughter_allTrack_startY, &b_reco_daughter_allTrack_startY);
   fChain->SetBranchAddress("reco_daughter_allTrack_startZ", &reco_daughter_allTrack_startZ, &b_reco_daughter_allTrack_startZ);
   fChain->SetBranchAddress("reco_daughter_allTrack_endX", &reco_daughter_allTrack_endX, &b_reco_daughter_allTrack_endX);
   fChain->SetBranchAddress("reco_daughter_allTrack_endY", &reco_daughter_allTrack_endY, &b_reco_daughter_allTrack_endY);
   fChain->SetBranchAddress("reco_daughter_allTrack_endZ", &reco_daughter_allTrack_endZ, &b_reco_daughter_allTrack_endZ);
   fChain->SetBranchAddress("reco_daughter_allTrack_dR", &reco_daughter_allTrack_dR, &b_reco_daughter_allTrack_dR);
   fChain->SetBranchAddress("reco_daughter_allTrack_to_vertex", &reco_daughter_allTrack_to_vertex, &b_reco_daughter_allTrack_to_vertex);
   fChain->SetBranchAddress("reco_daughter_shower_true_byE_PDG", &reco_daughter_shower_true_byE_PDG, &b_reco_daughter_shower_true_byE_PDG);
   fChain->SetBranchAddress("reco_daughter_shower_true_byE_ID", &reco_daughter_shower_true_byE_ID, &b_reco_daughter_shower_true_byE_ID);
   fChain->SetBranchAddress("reco_daughter_shower_true_byE_origin", &reco_daughter_shower_true_byE_origin, &b_reco_daughter_shower_true_byE_origin);
   fChain->SetBranchAddress("reco_daughter_shower_true_byE_parID", &reco_daughter_shower_true_byE_parID, &b_reco_daughter_shower_true_byE_parID);
   fChain->SetBranchAddress("reco_daughter_shower_true_byE_startPx", &reco_daughter_shower_true_byE_startPx, &b_reco_daughter_shower_true_byE_startPx);
   fChain->SetBranchAddress("reco_daughter_shower_true_byE_startPy", &reco_daughter_shower_true_byE_startPy, &b_reco_daughter_shower_true_byE_startPy);
   fChain->SetBranchAddress("reco_daughter_shower_true_byE_startPz", &reco_daughter_shower_true_byE_startPz, &b_reco_daughter_shower_true_byE_startPz);
   fChain->SetBranchAddress("reco_daughter_shower_true_byE_startP", &reco_daughter_shower_true_byE_startP, &b_reco_daughter_shower_true_byE_startP);
   fChain->SetBranchAddress("reco_daughter_shower_true_byHits_PDG", &reco_daughter_shower_true_byHits_PDG, &b_reco_daughter_shower_true_byHits_PDG);
   fChain->SetBranchAddress("reco_daughter_shower_true_byHits_ID", &reco_daughter_shower_true_byHits_ID, &b_reco_daughter_shower_true_byHits_ID);
   fChain->SetBranchAddress("reco_daughter_shower_true_byHits_origin", &reco_daughter_shower_true_byHits_origin, &b_reco_daughter_shower_true_byHits_origin);
   fChain->SetBranchAddress("reco_daughter_shower_true_byHits_parID", &reco_daughter_shower_true_byHits_parID, &b_reco_daughter_shower_true_byHits_parID);
   fChain->SetBranchAddress("reco_daughter_shower_true_byHits_process", &reco_daughter_shower_true_byHits_process, &b_reco_daughter_shower_true_byHits_process);
   fChain->SetBranchAddress("reco_daughter_shower_true_byHits_purity", &reco_daughter_shower_true_byHits_purity, &b_reco_daughter_shower_true_byHits_purity);
   fChain->SetBranchAddress("reco_daughter_shower_true_byHits_startPx", &reco_daughter_shower_true_byHits_startPx, &b_reco_daughter_shower_true_byHits_startPx);
   fChain->SetBranchAddress("reco_daughter_shower_true_byHits_startPy", &reco_daughter_shower_true_byHits_startPy, &b_reco_daughter_shower_true_byHits_startPy);
   fChain->SetBranchAddress("reco_daughter_shower_true_byHits_startPz", &reco_daughter_shower_true_byHits_startPz, &b_reco_daughter_shower_true_byHits_startPz);
   fChain->SetBranchAddress("reco_daughter_shower_true_byHits_startP", &reco_daughter_shower_true_byHits_startP, &b_reco_daughter_shower_true_byHits_startP);
   fChain->SetBranchAddress("reco_daughter_showerID", &reco_daughter_showerID, &b_reco_daughter_showerID);
   fChain->SetBranchAddress("reco_daughter_dQdX", &reco_daughter_dQdX, &b_reco_daughter_dQdX);
   fChain->SetBranchAddress("reco_daughter_dEdX", &reco_daughter_dEdX, &b_reco_daughter_dEdX);
   fChain->SetBranchAddress("reco_daughter_resRange", &reco_daughter_resRange, &b_reco_daughter_resRange);
   fChain->SetBranchAddress("reco_daughter_shower_dQdX", &reco_daughter_shower_dQdX, &b_reco_daughter_shower_dQdX);
   fChain->SetBranchAddress("reco_daughter_shower_dEdX", &reco_daughter_shower_dEdX, &b_reco_daughter_shower_dEdX);
   fChain->SetBranchAddress("reco_daughter_shower_resRange", &reco_daughter_shower_resRange, &b_reco_daughter_shower_resRange);
   fChain->SetBranchAddress("reco_daughter_len", &reco_daughter_len, &b_reco_daughter_len);
   fChain->SetBranchAddress("reco_daughter_startX", &reco_daughter_startX, &b_reco_daughter_startX);
   fChain->SetBranchAddress("reco_daughter_startY", &reco_daughter_startY, &b_reco_daughter_startY);
   fChain->SetBranchAddress("reco_daughter_startZ", &reco_daughter_startZ, &b_reco_daughter_startZ);
   fChain->SetBranchAddress("reco_daughter_endX", &reco_daughter_endX, &b_reco_daughter_endX);
   fChain->SetBranchAddress("reco_daughter_endY", &reco_daughter_endY, &b_reco_daughter_endY);
   fChain->SetBranchAddress("reco_daughter_endZ", &reco_daughter_endZ, &b_reco_daughter_endZ);
   fChain->SetBranchAddress("reco_daughter_deltaR", &reco_daughter_deltaR, &b_reco_daughter_deltaR);
   fChain->SetBranchAddress("reco_daughter_dR", &reco_daughter_dR, &b_reco_daughter_dR);
   fChain->SetBranchAddress("reco_daughter_to_vertex", &reco_daughter_to_vertex, &b_reco_daughter_to_vertex);
   fChain->SetBranchAddress("reco_daughter_slice", &reco_daughter_slice, &b_reco_daughter_slice);
   fChain->SetBranchAddress("reco_daughter_shower_to_vertex", &reco_daughter_shower_to_vertex, &b_reco_daughter_shower_to_vertex);
   fChain->SetBranchAddress("reco_daughter_shower_startX", &reco_daughter_shower_startX, &b_reco_daughter_shower_startX);
   fChain->SetBranchAddress("reco_daughter_shower_startY", &reco_daughter_shower_startY, &b_reco_daughter_shower_startY);
   fChain->SetBranchAddress("reco_daughter_shower_startZ", &reco_daughter_shower_startZ, &b_reco_daughter_shower_startZ);
   fChain->SetBranchAddress("reco_daughter_shower_len", &reco_daughter_shower_len, &b_reco_daughter_shower_len);
   fChain->SetBranchAddress("reco_daughter_PFP_ID", &reco_daughter_PFP_ID, &b_reco_daughter_PFP_ID);
   fChain->SetBranchAddress("reco_daughter_PFP_trackScore", &reco_daughter_PFP_trackScore, &b_reco_daughter_PFP_trackScore);
   fChain->SetBranchAddress("reco_daughter_PFP_emScore", &reco_daughter_PFP_emScore, &b_reco_daughter_PFP_emScore);
   fChain->SetBranchAddress("reco_daughter_PFP_michelScore", &reco_daughter_PFP_michelScore, &b_reco_daughter_PFP_michelScore);
   fChain->SetBranchAddress("true_beam_PDG", &true_beam_PDG, &b_true_beam_PDG);
   fChain->SetBranchAddress("true_beam_ID", &true_beam_ID, &b_true_beam_ID);
   fChain->SetBranchAddress("true_beam_endProcess", &true_beam_endProcess, &b_true_beam_endProcess);
   fChain->SetBranchAddress("true_beam_endX", &true_beam_endX, &b_true_beam_endX);
   fChain->SetBranchAddress("true_beam_endY", &true_beam_endY, &b_true_beam_endY);
   fChain->SetBranchAddress("true_beam_endZ", &true_beam_endZ, &b_true_beam_endZ);
   fChain->SetBranchAddress("true_beam_startX", &true_beam_startX, &b_true_beam_startX);
   fChain->SetBranchAddress("true_beam_startY", &true_beam_startY, &b_true_beam_startY);
   fChain->SetBranchAddress("true_beam_startZ", &true_beam_startZ, &b_true_beam_startZ);
   fChain->SetBranchAddress("true_beam_startPx", &true_beam_startPx, &b_true_beam_startPx);
   fChain->SetBranchAddress("true_beam_startPy", &true_beam_startPy, &b_true_beam_startPy);
   fChain->SetBranchAddress("true_beam_startPz", &true_beam_startPz, &b_true_beam_startPz);
   fChain->SetBranchAddress("true_beam_startP", &true_beam_startP, &b_true_beam_startP);
   fChain->SetBranchAddress("true_beam_endPx", &true_beam_endPx, &b_true_beam_endPx);
   fChain->SetBranchAddress("true_beam_endPy", &true_beam_endPy, &b_true_beam_endPy);
   fChain->SetBranchAddress("true_beam_endPz", &true_beam_endPz, &b_true_beam_endPz);
   fChain->SetBranchAddress("true_beam_endP", &true_beam_endP, &b_true_beam_endP);
   fChain->SetBranchAddress("true_beam_startDirX", &true_beam_startDirX, &b_true_beam_startDirX);
   fChain->SetBranchAddress("true_beam_startDirY", &true_beam_startDirY, &b_true_beam_startDirY);
   fChain->SetBranchAddress("true_beam_startDirZ", &true_beam_startDirZ, &b_true_beam_startDirZ);
   fChain->SetBranchAddress("true_beam_nElasticScatters", &true_beam_nElasticScatters, &b_true_beam_nElasticScatters);
   fChain->SetBranchAddress("true_beam_elastic_costheta", &true_beam_elastic_costheta, &b_true_beam_elastic_costheta);
   fChain->SetBranchAddress("true_beam_elastic_X", &true_beam_elastic_X, &b_true_beam_elastic_X);
   fChain->SetBranchAddress("true_beam_elastic_Y", &true_beam_elastic_Y, &b_true_beam_elastic_Y);
   fChain->SetBranchAddress("true_beam_elastic_Z", &true_beam_elastic_Z, &b_true_beam_elastic_Z);
   fChain->SetBranchAddress("true_beam_IDE_totalDep", &true_beam_IDE_totalDep, &b_true_beam_IDE_totalDep);
   fChain->SetBranchAddress("true_beam_IDE_found_in_recoVtx", &true_beam_IDE_found_in_recoVtx, &b_true_beam_IDE_found_in_recoVtx);
   fChain->SetBranchAddress("true_daughter_nPi0", &true_daughter_nPi0, &b_true_daughter_nPi0);
   fChain->SetBranchAddress("true_daughter_nPiPlus", &true_daughter_nPiPlus, &b_true_daughter_nPiPlus);
   fChain->SetBranchAddress("true_daughter_nProton", &true_daughter_nProton, &b_true_daughter_nProton);
   fChain->SetBranchAddress("true_daughter_nNeutron", &true_daughter_nNeutron, &b_true_daughter_nNeutron);
   fChain->SetBranchAddress("true_daughter_nPiMinus", &true_daughter_nPiMinus, &b_true_daughter_nPiMinus);
   fChain->SetBranchAddress("true_daughter_nNucleus", &true_daughter_nNucleus, &b_true_daughter_nNucleus);
   fChain->SetBranchAddress("reco_beam_vertex_slice", &reco_beam_vertex_slice, &b_reco_beam_vertex_slice);
   fChain->SetBranchAddress("reco_beam_vertex_dRs", &reco_beam_vertex_dRs, &b_reco_beam_vertex_dRs);
   fChain->SetBranchAddress("reco_beam_vertex_hits_slices", &reco_beam_vertex_hits_slices, &b_reco_beam_vertex_hits_slices);
   fChain->SetBranchAddress("true_beam_daughter_PDG", &true_beam_daughter_PDG, &b_true_beam_daughter_PDG);
   fChain->SetBranchAddress("true_beam_daughter_ID", &true_beam_daughter_ID, &b_true_beam_daughter_ID);
   fChain->SetBranchAddress("true_beam_daughter_len", &true_beam_daughter_len, &b_true_beam_daughter_len);
   fChain->SetBranchAddress("true_beam_daughter_startX", &true_beam_daughter_startX, &b_true_beam_daughter_startX);
   fChain->SetBranchAddress("true_beam_daughter_startY", &true_beam_daughter_startY, &b_true_beam_daughter_startY);
   fChain->SetBranchAddress("true_beam_daughter_startZ", &true_beam_daughter_startZ, &b_true_beam_daughter_startZ);
   fChain->SetBranchAddress("true_beam_daughter_startPx", &true_beam_daughter_startPx, &b_true_beam_daughter_startPx);
   fChain->SetBranchAddress("true_beam_daughter_startPy", &true_beam_daughter_startPy, &b_true_beam_daughter_startPy);
   fChain->SetBranchAddress("true_beam_daughter_startPz", &true_beam_daughter_startPz, &b_true_beam_daughter_startPz);
   fChain->SetBranchAddress("true_beam_daughter_startP", &true_beam_daughter_startP, &b_true_beam_daughter_startP);
   fChain->SetBranchAddress("true_beam_daughter_endX", &true_beam_daughter_endX, &b_true_beam_daughter_endX);
   fChain->SetBranchAddress("true_beam_daughter_endY", &true_beam_daughter_endY, &b_true_beam_daughter_endY);
   fChain->SetBranchAddress("true_beam_daughter_endZ", &true_beam_daughter_endZ, &b_true_beam_daughter_endZ);
   fChain->SetBranchAddress("true_beam_daughter_Process", &true_beam_daughter_Process, &b_true_beam_daughter_Process);
   fChain->SetBranchAddress("true_beam_Pi0_decay_ID", &true_beam_Pi0_decay_ID, &b_true_beam_Pi0_decay_ID);
   fChain->SetBranchAddress("true_beam_Pi0_decay_PDG", &true_beam_Pi0_decay_PDG, &b_true_beam_Pi0_decay_PDG);
   fChain->SetBranchAddress("true_beam_Pi0_decay_startP", &true_beam_Pi0_decay_startP, &b_true_beam_Pi0_decay_startP);
   fChain->SetBranchAddress("true_beam_grand_daughter_ID", &true_beam_grand_daughter_ID, &b_true_beam_grand_daughter_ID);
   fChain->SetBranchAddress("true_beam_grand_daughter_parID", &true_beam_grand_daughter_parID, &b_true_beam_grand_daughter_parID);
   fChain->SetBranchAddress("true_beam_grand_daughter_PDG", &true_beam_grand_daughter_PDG, &b_true_beam_grand_daughter_PDG);
   fChain->SetBranchAddress("reco_beam_true_byE_endProcess", &reco_beam_true_byE_endProcess, &b_reco_beam_true_byE_endProcess);
   fChain->SetBranchAddress("reco_beam_true_byE_process", &reco_beam_true_byE_process, &b_reco_beam_true_byE_process);
   fChain->SetBranchAddress("reco_beam_true_byE_origin", &reco_beam_true_byE_origin, &b_reco_beam_true_byE_origin);
   fChain->SetBranchAddress("reco_beam_true_byE_PDG", &reco_beam_true_byE_PDG, &b_reco_beam_true_byE_PDG);
   fChain->SetBranchAddress("reco_beam_true_byE_ID", &reco_beam_true_byE_ID, &b_reco_beam_true_byE_ID);
   fChain->SetBranchAddress("reco_beam_true_byHits_endProcess", &reco_beam_true_byHits_endProcess, &b_reco_beam_true_byHits_endProcess);
   fChain->SetBranchAddress("reco_beam_true_byHits_process", &reco_beam_true_byHits_process, &b_reco_beam_true_byHits_process);
   fChain->SetBranchAddress("reco_beam_true_byHits_origin", &reco_beam_true_byHits_origin, &b_reco_beam_true_byHits_origin);
   fChain->SetBranchAddress("reco_beam_true_byHits_PDG", &reco_beam_true_byHits_PDG, &b_reco_beam_true_byHits_PDG);
   fChain->SetBranchAddress("reco_beam_true_byHits_ID", &reco_beam_true_byHits_ID, &b_reco_beam_true_byHits_ID);
   fChain->SetBranchAddress("reco_beam_true_byE_matched", &reco_beam_true_byE_matched, &b_reco_beam_true_byE_matched);
   fChain->SetBranchAddress("reco_beam_true_byHits_matched", &reco_beam_true_byHits_matched, &b_reco_beam_true_byHits_matched);
   fChain->SetBranchAddress("reco_beam_true_byHits_purity", &reco_beam_true_byHits_purity, &b_reco_beam_true_byHits_purity);
   fChain->SetBranchAddress("true_beam_processes", &true_beam_processes, &b_true_beam_processes);
   fChain->SetBranchAddress("reco_daughter_true_byE_isPrimary", &reco_daughter_true_byE_isPrimary, &b_reco_daughter_true_byE_isPrimary);
   fChain->SetBranchAddress("data_BI_P", &data_BI_P, &b_data_BI_P);
   fChain->SetBranchAddress("data_BI_X", &data_BI_X, &b_data_BI_X);
   fChain->SetBranchAddress("data_BI_Y", &data_BI_Y, &b_data_BI_Y);
   fChain->SetBranchAddress("data_BI_Z", &data_BI_Z, &b_data_BI_Z);
   fChain->SetBranchAddress("data_BI_nFibersP1", &data_BI_nFibersP1, &b_data_BI_nFibersP1);
   fChain->SetBranchAddress("data_BI_nFibersP2", &data_BI_nFibersP2, &b_data_BI_nFibersP2);
   fChain->SetBranchAddress("data_BI_nFibersP3", &data_BI_nFibersP3, &b_data_BI_nFibersP3);
   fChain->SetBranchAddress("data_BI_PDG_candidates", &data_BI_PDG_candidates, &b_data_BI_PDG_candidates);
   fChain->SetBranchAddress("quality_reco_view_0_hits_in_TPC5", &quality_reco_view_0_hits_in_TPC5, &b_quality_reco_view_0_hits_in_TPC5);
   fChain->SetBranchAddress("quality_reco_view_1_hits_in_TPC5", &quality_reco_view_1_hits_in_TPC5, &b_quality_reco_view_1_hits_in_TPC5);
   fChain->SetBranchAddress("quality_reco_view_2_hits_in_TPC5", &quality_reco_view_2_hits_in_TPC5, &b_quality_reco_view_2_hits_in_TPC5);
   fChain->SetBranchAddress("quality_reco_max_lateral", &quality_reco_max_lateral, &b_quality_reco_max_lateral);
   fChain->SetBranchAddress("quality_reco_max_segment", &quality_reco_max_segment, &b_quality_reco_max_segment);
   fChain->SetBranchAddress("quality_reco_view_0_max_segment", &quality_reco_view_0_max_segment, &b_quality_reco_view_0_max_segment);
   fChain->SetBranchAddress("quality_reco_view_1_max_segment", &quality_reco_view_1_max_segment, &b_quality_reco_view_1_max_segment);
   fChain->SetBranchAddress("quality_reco_view_2_max_segment", &quality_reco_view_2_max_segment, &b_quality_reco_view_2_max_segment);
   fChain->SetBranchAddress("quality_reco_view_0_wire_backtrack", &quality_reco_view_0_wire_backtrack, &b_quality_reco_view_0_wire_backtrack);
   fChain->SetBranchAddress("quality_reco_view_1_wire_backtrack", &quality_reco_view_1_wire_backtrack, &b_quality_reco_view_1_wire_backtrack);
   fChain->SetBranchAddress("quality_reco_view_2_wire_backtrack", &quality_reco_view_2_wire_backtrack, &b_quality_reco_view_2_wire_backtrack);
   fChain->SetBranchAddress("quality_reco_view_0_wire", &quality_reco_view_0_wire, &b_quality_reco_view_0_wire);
   fChain->SetBranchAddress("quality_reco_view_1_wire", &quality_reco_view_1_wire, &b_quality_reco_view_1_wire);
   fChain->SetBranchAddress("quality_reco_view_2_wire", &quality_reco_view_2_wire, &b_quality_reco_view_2_wire);
   fChain->SetBranchAddress("quality_reco_view_2_z", &quality_reco_view_2_z, &b_quality_reco_view_2_z);
   fChain->SetBranchAddress("quality_reco_view_0_tick", &quality_reco_view_0_tick, &b_quality_reco_view_0_tick);
   fChain->SetBranchAddress("quality_reco_view_1_tick", &quality_reco_view_1_tick, &b_quality_reco_view_1_tick);
   fChain->SetBranchAddress("quality_reco_view_2_tick", &quality_reco_view_2_tick, &b_quality_reco_view_2_tick);
   fChain->SetBranchAddress("reco_beam_Chi2_proton", &reco_beam_Chi2_proton, &b_reco_beam_Chi2_proton);
   fChain->SetBranchAddress("reco_beam_Chi2_ndof", &reco_beam_Chi2_ndof, &b_reco_beam_Chi2_ndof);
   fChain->SetBranchAddress("reco_daughter_Chi2_proton", &reco_daughter_Chi2_proton, &b_reco_daughter_Chi2_proton);
   fChain->SetBranchAddress("reco_daughter_Chi2_ndof", &reco_daughter_Chi2_ndof, &b_reco_daughter_Chi2_ndof);
   fChain->SetBranchAddress("reco_daughter_momByRange_proton", &reco_daughter_momByRange_proton, &b_reco_daughter_momByRange_proton);
   fChain->SetBranchAddress("reco_daughter_momByRange_muon", &reco_daughter_momByRange_muon, &b_reco_daughter_momByRange_muon);
   fChain->SetBranchAddress("reco_daughter_allTrack_momByRange_proton", &reco_daughter_allTrack_momByRange_proton, &b_reco_daughter_allTrack_momByRange_proton);
   fChain->SetBranchAddress("reco_daughter_allTrack_momByRange_muon", &reco_daughter_allTrack_momByRange_muon, &b_reco_daughter_allTrack_momByRange_muon);
   fChain->SetBranchAddress("reco_daughter_shower_Chi2_proton", &reco_daughter_shower_Chi2_proton, &b_reco_daughter_shower_Chi2_proton);
   fChain->SetBranchAddress("reco_daughter_shower_Chi2_ndof", &reco_daughter_shower_Chi2_ndof, &b_reco_daughter_shower_Chi2_ndof);
   fChain->SetBranchAddress("reco_daughter_trackScore", &reco_daughter_trackScore, &b_reco_daughter_trackScore);
   fChain->SetBranchAddress("reco_daughter_emScore", &reco_daughter_emScore, &b_reco_daughter_emScore);
   fChain->SetBranchAddress("reco_daughter_michelScore", &reco_daughter_michelScore, &b_reco_daughter_michelScore);
   fChain->SetBranchAddress("reco_daughter_shower_trackScore", &reco_daughter_shower_trackScore, &b_reco_daughter_shower_trackScore);
   fChain->SetBranchAddress("reco_daughter_shower_emScore", &reco_daughter_shower_emScore, &b_reco_daughter_shower_emScore);
   fChain->SetBranchAddress("reco_daughter_shower_michelScore", &reco_daughter_shower_michelScore, &b_reco_daughter_shower_michelScore);
   fChain->SetBranchAddress("reco_beam_true_byE_endPx", &reco_beam_true_byE_endPx, &b_reco_beam_true_byE_endPx);
   fChain->SetBranchAddress("reco_beam_true_byE_endPy", &reco_beam_true_byE_endPy, &b_reco_beam_true_byE_endPy);
   fChain->SetBranchAddress("reco_beam_true_byE_endPz", &reco_beam_true_byE_endPz, &b_reco_beam_true_byE_endPz);
   fChain->SetBranchAddress("reco_beam_true_byE_endE", &reco_beam_true_byE_endE, &b_reco_beam_true_byE_endE);
   fChain->SetBranchAddress("reco_beam_true_byE_endP", &reco_beam_true_byE_endP, &b_reco_beam_true_byE_endP);
   fChain->SetBranchAddress("reco_beam_true_byE_startPx", &reco_beam_true_byE_startPx, &b_reco_beam_true_byE_startPx);
   fChain->SetBranchAddress("reco_beam_true_byE_startPy", &reco_beam_true_byE_startPy, &b_reco_beam_true_byE_startPy);
   fChain->SetBranchAddress("reco_beam_true_byE_startPz", &reco_beam_true_byE_startPz, &b_reco_beam_true_byE_startPz);
   fChain->SetBranchAddress("reco_beam_true_byE_startE", &reco_beam_true_byE_startE, &b_reco_beam_true_byE_startE);
   fChain->SetBranchAddress("reco_beam_true_byE_startP", &reco_beam_true_byE_startP, &b_reco_beam_true_byE_startP);
   fChain->SetBranchAddress("reco_beam_true_byHits_endPx", &reco_beam_true_byHits_endPx, &b_reco_beam_true_byHits_endPx);
   fChain->SetBranchAddress("reco_beam_true_byHits_endPy", &reco_beam_true_byHits_endPy, &b_reco_beam_true_byHits_endPy);
   fChain->SetBranchAddress("reco_beam_true_byHits_endPz", &reco_beam_true_byHits_endPz, &b_reco_beam_true_byHits_endPz);
   fChain->SetBranchAddress("reco_beam_true_byHits_endE", &reco_beam_true_byHits_endE, &b_reco_beam_true_byHits_endE);
   fChain->SetBranchAddress("reco_beam_true_byHits_endP", &reco_beam_true_byHits_endP, &b_reco_beam_true_byHits_endP);
   fChain->SetBranchAddress("reco_beam_true_byHits_startPx", &reco_beam_true_byHits_startPx, &b_reco_beam_true_byHits_startPx);
   fChain->SetBranchAddress("reco_beam_true_byHits_startPy", &reco_beam_true_byHits_startPy, &b_reco_beam_true_byHits_startPy);
   fChain->SetBranchAddress("reco_beam_true_byHits_startPz", &reco_beam_true_byHits_startPz, &b_reco_beam_true_byHits_startPz);
   fChain->SetBranchAddress("reco_beam_true_byHits_startE", &reco_beam_true_byHits_startE, &b_reco_beam_true_byHits_startE);
   fChain->SetBranchAddress("reco_beam_true_byHits_startP", &reco_beam_true_byHits_startP, &b_reco_beam_true_byHits_startP);
   Notify();
}

Bool_t PDEventPro::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PDEventPro::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PDEventPro::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef PDEventPro_cxx
