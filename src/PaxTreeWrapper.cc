#ifndef paxtreewrapper_cc
#define paxtreewrapper_cc

#include "PaxTreeWrapper.hh"
#include "PaxClasses.hh"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


void PaxTreeWrapper::Loop()
{
//   In a ROOT session, you can do:
//      root> .L tree.C
//      root> tree t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}


PaxTreeWrapper::PaxTreeWrapper(string filename) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
	fInit = false;
	
	//events = new Event;
	
	TFile *f = TFile::Open(filename.c_str(), "read");
	if (!f || !f->IsOpen())
	{
		return;
	}
	
	TTree *tree = (TTree*)f->Get("tree");
	if(!tree)
	{
		delete f;
		return;
	}
	
	Init(tree);
}



PaxTreeWrapper::PaxTreeWrapper(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
	
	fInit = false;
	
	//events = new Event;
	
	if(tree) {
		Init(tree);
	}
}


PaxTreeWrapper::~PaxTreeWrapper()
{
	if (fChain)
	{
		delete fChain->GetCurrentFile();
	}
}


Int_t PaxTreeWrapper::GetEntry(Long64_t entry)
{
// Read contents of entry.
	if (!fChain) return 0;
	return fChain->GetEntry(entry);
}


Long64_t PaxTreeWrapper::LoadTree(Long64_t entry)
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

void PaxTreeWrapper::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
	if (!tree)
	{
		fInit=false;
		return;
	}
	
	fChain = tree;
	fCurrent = -1;
	fChain->SetMakeClass(1);
	
	//fChain->SetBranchAddress("events", &events, &b_events);
	
    fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_events_fUniqueID);
    fChain->SetBranchAddress("fBits", &fBits, &b_events_fBits);
    fChain->SetBranchAddress("all_hits", &all_hits_, &b_events_all_hits_);
    fChain->SetBranchAddress("all_hits.fUniqueID", &all_hits_fUniqueID, &b_all_hits_fUniqueID);
    fChain->SetBranchAddress("all_hits.fBits", &all_hits_fBits, &b_all_hits_fBits);
    fChain->SetBranchAddress("all_hits.area", &all_hits_area, &b_all_hits_area);
    fChain->SetBranchAddress("all_hits.center", &all_hits_center, &b_all_hits_center);
    fChain->SetBranchAddress("all_hits.channel", &all_hits_channel, &b_all_hits_channel);
    fChain->SetBranchAddress("all_hits.found_in_pulse", &all_hits_found_in_pulse, &b_all_hits_found_in_pulse);
    fChain->SetBranchAddress("all_hits.height", &all_hits_height, &b_all_hits_height);
    fChain->SetBranchAddress("all_hits.index_of_maximum", &all_hits_index_of_maximum, &b_all_hits_index_of_maximum);
    fChain->SetBranchAddress("all_hits.is_rejected", &all_hits_is_rejected, &b_all_hits_is_rejected);
    fChain->SetBranchAddress("all_hits.left", &all_hits_left, &b_all_hits_left);
    fChain->SetBranchAddress("all_hits.left_central", &all_hits_left_central, &b_all_hits_left_central);
    fChain->SetBranchAddress("all_hits.n_saturated", &all_hits_n_saturated, &b_all_hits_n_saturated);
    fChain->SetBranchAddress("all_hits.noise_sigma", &all_hits_noise_sigma, &b_all_hits_noise_sigma);
    fChain->SetBranchAddress("all_hits.right", &all_hits_right, &b_all_hits_right);
    fChain->SetBranchAddress("all_hits.right_central", &all_hits_right_central, &b_all_hits_right_central);
    fChain->SetBranchAddress("all_hits.sum_absolute_deviation", &all_hits_sum_absolute_deviation, &b_all_hits_sum_absolute_deviation);
    fChain->SetBranchAddress("block_id", &block_id, &b_events_block_id);
    fChain->SetBranchAddress("dataset_name", &dataset_name, &b_events_dataset_name);
    fChain->SetBranchAddress("event_number", &event_number, &b_events_event_number);
    fChain->SetBranchAddress("interactions", &interactions_, &b_events_interactions_);
    fChain->SetBranchAddress("interactions.fUniqueID", interactions_fUniqueID, &b_interactions_fUniqueID);
    fChain->SetBranchAddress("interactions.fBits", interactions_fBits, &b_interactions_fBits);
    fChain->SetBranchAddress("interactions.drift_time", interactions_drift_time, &b_interactions_drift_time);
    fChain->SetBranchAddress("interactions.r_correction", interactions_r_correction, &b_interactions_r_correction);
    fChain->SetBranchAddress("interactions.s1", interactions_s1, &b_interactions_s1);
    fChain->SetBranchAddress("interactions.s1_area_correction", interactions_s1_area_correction, &b_interactions_s1_area_correction);
    fChain->SetBranchAddress("interactions.s1_area_fraction_top_probability", interactions_s1_area_fraction_top_probability, &b_interactions_s1_area_fraction_top_probability);
    fChain->SetBranchAddress("interactions.s1_pattern_fit", interactions_s1_pattern_fit, &b_interactions_s1_pattern_fit);
    fChain->SetBranchAddress("interactions.s1_pattern_fit_hits", interactions_s1_pattern_fit_hits, &b_interactions_s1_pattern_fit_hits);
    fChain->SetBranchAddress("interactions.s1_saturation_correction", interactions_s1_saturation_correction, &b_interactions_s1_saturation_correction);
    fChain->SetBranchAddress("interactions.s1_spatial_correction", interactions_s1_spatial_correction, &b_interactions_s1_spatial_correction);
    fChain->SetBranchAddress("interactions.s2", interactions_s2, &b_interactions_s2);
    fChain->SetBranchAddress("interactions.s2_area_correction", interactions_s2_area_correction, &b_interactions_s2_area_correction);
    fChain->SetBranchAddress("interactions.s2_lifetime_correction", interactions_s2_lifetime_correction, &b_interactions_s2_lifetime_correction);
    fChain->SetBranchAddress("interactions.x", interactions_x, &b_interactions_x);
    fChain->SetBranchAddress("interactions.xy_posrec_algorithm", interactions_xy_posrec_algorithm, &b_interactions_xy_posrec_algorithm);
    fChain->SetBranchAddress("interactions.xy_posrec_goodness_of_fit", interactions_xy_posrec_goodness_of_fit, &b_interactions_xy_posrec_goodness_of_fit);
    fChain->SetBranchAddress("interactions.xy_posrec_ndf", interactions_xy_posrec_ndf, &b_interactions_xy_posrec_ndf);
    fChain->SetBranchAddress("interactions.y", interactions_y, &b_interactions_y);
    fChain->SetBranchAddress("interactions.z", interactions_z, &b_interactions_z);
    fChain->SetBranchAddress("interactions.z_correction", interactions_z_correction, &b_interactions_z_correction);
    fChain->SetBranchAddress("is_channel_suspicious[260]", is_channel_suspicious, &b_events_is_channel_suspicious);
    fChain->SetBranchAddress("lone_hits_per_channel[260]", lone_hits_per_channel, &b_events_lone_hits_per_channel);
    fChain->SetBranchAddress("lone_hits_per_channel_before[260]", lone_hits_per_channel_before, &b_events_lone_hits_per_channel_before);
    fChain->SetBranchAddress("n_channels", &n_channels, &b_events_n_channels);
    fChain->SetBranchAddress("n_hits_rejected[260]", n_hits_rejected, &b_events_n_hits_rejected);
    fChain->SetBranchAddress("n_pulses", &n_pulses, &b_events_n_pulses);
    fChain->SetBranchAddress("n_pulses_per_channel[260]", n_pulses_per_channel, &b_events_n_pulses_per_channel);
    fChain->SetBranchAddress("noise_pulses_in[260]", noise_pulses_in, &b_events_noise_pulses_in);
    fChain->SetBranchAddress("peaks", &peaks_, &b_events_peaks_);
    fChain->SetBranchAddress("peaks.fUniqueID", peaks_fUniqueID, &b_peaks_fUniqueID);
    fChain->SetBranchAddress("peaks.fBits", peaks_fBits, &b_peaks_fBits);
    fChain->SetBranchAddress("peaks.area", peaks_area, &b_peaks_area);
    fChain->SetBranchAddress("peaks.area_decile_from_midpoint[11]", peaks_area_decile_from_midpoint, &b_peaks_area_decile_from_midpoint);
    fChain->SetBranchAddress("peaks.area_fraction_top", peaks_area_fraction_top, &b_peaks_area_fraction_top);
    fChain->SetBranchAddress("peaks.area_midpoint", peaks_area_midpoint, &b_peaks_area_midpoint);
    fChain->SetBranchAddress("peaks.area_per_channel[260]", peaks_area_per_channel, &b_peaks_area_per_channel);
    fChain->SetBranchAddress("peaks.birthing_split_fraction", peaks_birthing_split_fraction, &b_peaks_birthing_split_fraction);
    fChain->SetBranchAddress("peaks.birthing_split_goodness", peaks_birthing_split_goodness, &b_peaks_birthing_split_goodness);
    fChain->SetBranchAddress("peaks.bottom_hitpattern_spread", peaks_bottom_hitpattern_spread, &b_peaks_bottom_hitpattern_spread);
    fChain->SetBranchAddress("peaks.center_time", peaks_center_time, &b_peaks_center_time);
    fChain->SetBranchAddress("peaks.detector", peaks_detector, &b_peaks_detector);
    fChain->SetBranchAddress("peaks.height", peaks_height, &b_peaks_height);
    fChain->SetBranchAddress("peaks.hit_time_mean", peaks_hit_time_mean, &b_peaks_hit_time_mean);
    fChain->SetBranchAddress("peaks.hit_time_std", peaks_hit_time_std, &b_peaks_hit_time_std);
    fChain->SetBranchAddress("peaks.hits", peaks_hits, &b_peaks_hits);
    fChain->SetBranchAddress("peaks.hits_fraction_top", peaks_hits_fraction_top, &b_peaks_hits_fraction_top);
    fChain->SetBranchAddress("peaks.hits_per_channel[260]", peaks_hits_per_channel, &b_peaks_hits_per_channel);
    fChain->SetBranchAddress("peaks.index_of_maximum", peaks_index_of_maximum, &b_peaks_index_of_maximum);
    fChain->SetBranchAddress("peaks.interior_split_fraction", peaks_interior_split_fraction, &b_peaks_interior_split_fraction);
    fChain->SetBranchAddress("peaks.interior_split_goodness", peaks_interior_split_goodness, &b_peaks_interior_split_goodness);
    fChain->SetBranchAddress("peaks.largest_hit_area", peaks_largest_hit_area, &b_peaks_largest_hit_area);
    fChain->SetBranchAddress("peaks.largest_hit_channel", peaks_largest_hit_channel, &b_peaks_largest_hit_channel);
    fChain->SetBranchAddress("peaks.left", peaks_left, &b_peaks_left);
    fChain->SetBranchAddress("peaks.lone_hit_channel", peaks_lone_hit_channel, &b_peaks_lone_hit_channel);
    fChain->SetBranchAddress("peaks.mean_amplitude_to_noise", peaks_mean_amplitude_to_noise, &b_peaks_mean_amplitude_to_noise);
    fChain->SetBranchAddress("peaks.n_contributing_channels", peaks_n_contributing_channels, &b_peaks_n_contributing_channels);
    fChain->SetBranchAddress("peaks.n_contributing_channels_top", peaks_n_contributing_channels_top, &b_peaks_n_contributing_channels_top);
    fChain->SetBranchAddress("peaks.n_hits", peaks_n_hits, &b_peaks_n_hits);
    fChain->SetBranchAddress("peaks.n_noise_pulses", peaks_n_noise_pulses, &b_peaks_n_noise_pulses);
    fChain->SetBranchAddress("peaks.n_saturated_channels", peaks_n_saturated_channels, &b_peaks_n_saturated_channels);
    fChain->SetBranchAddress("peaks.n_saturated_per_channel[260]", peaks_n_saturated_per_channel, &b_peaks_n_saturated_per_channel);
    fChain->SetBranchAddress("peaks.n_saturated_samples", peaks_n_saturated_samples, &b_peaks_n_saturated_samples);
    fChain->SetBranchAddress("peaks.range_area_decile[11]", peaks_range_area_decile, &b_peaks_range_area_decile);
    fChain->SetBranchAddress("peaks.right", peaks_right, &b_peaks_right);
    fChain->SetBranchAddress("peaks.s2_bottom_spatial_correction", peaks_s2_bottom_spatial_correction, &b_peaks_s2_bottom_spatial_correction);
    fChain->SetBranchAddress("peaks.s2_saturation_correction", peaks_s2_saturation_correction, &b_peaks_s2_saturation_correction);
    fChain->SetBranchAddress("peaks.s2_spatial_correction", peaks_s2_spatial_correction, &b_peaks_s2_spatial_correction);
    fChain->SetBranchAddress("peaks.s2_top_spatial_correction", peaks_s2_top_spatial_correction, &b_peaks_s2_top_spatial_correction);
    fChain->SetBranchAddress("peaks.sum_waveform[251]", peaks_sum_waveform, &b_peaks_sum_waveform);
    fChain->SetBranchAddress("peaks.sum_waveform_top[251]", peaks_sum_waveform_top, &b_peaks_sum_waveform_top);
    fChain->SetBranchAddress("peaks.tight_coincidence", peaks_tight_coincidence, &b_peaks_tight_coincidence);
    fChain->SetBranchAddress("peaks.top_hitpattern_spread", peaks_top_hitpattern_spread, &b_peaks_top_hitpattern_spread);
    fChain->SetBranchAddress("peaks.type", peaks_type, &b_peaks_type);
    fChain->SetBranchAddress("pulses", &pulses_, &b_events_pulses_);
    fChain->SetBranchAddress("pulses.fUniqueID", pulses_fUniqueID, &b_pulses_fUniqueID);
    fChain->SetBranchAddress("pulses.fBits", pulses_fBits, &b_pulses_fBits);
    fChain->SetBranchAddress("pulses.baseline", pulses_baseline, &b_pulses_baseline);
    fChain->SetBranchAddress("pulses.baseline_increase", pulses_baseline_increase, &b_pulses_baseline_increase);
    fChain->SetBranchAddress("pulses.channel", pulses_channel, &b_pulses_channel);
    fChain->SetBranchAddress("pulses.hitfinder_threshold", pulses_hitfinder_threshold, &b_pulses_hitfinder_threshold);
    fChain->SetBranchAddress("pulses.left", pulses_left, &b_pulses_left);
    fChain->SetBranchAddress("pulses.maximum", pulses_maximum, &b_pulses_maximum);
    fChain->SetBranchAddress("pulses.minimum", pulses_minimum, &b_pulses_minimum);
    fChain->SetBranchAddress("pulses.n_hits_found", pulses_n_hits_found, &b_pulses_n_hits_found);
    fChain->SetBranchAddress("pulses.noise_sigma", pulses_noise_sigma, &b_pulses_noise_sigma);
    fChain->SetBranchAddress("pulses.right", pulses_right, &b_pulses_right);
    fChain->SetBranchAddress("sample_duration", &sample_duration, &b_events_sample_duration);
    fChain->SetBranchAddress("start_time", &start_time, &b_events_start_time);
    fChain->SetBranchAddress("stop_time", &stop_time, &b_events_stop_time);
    fChain->SetBranchAddress("trigger_signals", &trigger_signals_, &b_events_trigger_signals_);
    fChain->SetBranchAddress("trigger_signals.fUniqueID", trigger_signals_fUniqueID, &b_trigger_signals_fUniqueID);
    fChain->SetBranchAddress("trigger_signals.fBits", trigger_signals_fBits, &b_trigger_signals_fBits);
    fChain->SetBranchAddress("trigger_signals.area", trigger_signals_area, &b_trigger_signals_area);
    fChain->SetBranchAddress("trigger_signals.left_time", trigger_signals_left_time, &b_trigger_signals_left_time);
    fChain->SetBranchAddress("trigger_signals.n_contributing_channels", trigger_signals_n_contributing_channels, &b_trigger_signals_n_contributing_channels);
    fChain->SetBranchAddress("trigger_signals.n_pulses", trigger_signals_n_pulses, &b_trigger_signals_n_pulses);
    fChain->SetBranchAddress("trigger_signals.right_time", trigger_signals_right_time, &b_trigger_signals_right_time);
    fChain->SetBranchAddress("trigger_signals.time_mean", trigger_signals_time_mean, &b_trigger_signals_time_mean);
    fChain->SetBranchAddress("trigger_signals.time_rms", trigger_signals_time_rms, &b_trigger_signals_time_rms);
    fChain->SetBranchAddress("trigger_signals.trigger", trigger_signals_trigger, &b_trigger_signals_trigger);
    fChain->SetBranchAddress("trigger_signals.type", trigger_signals_type, &b_trigger_signals_type);
    fChain->SetBranchAddress("s1s", &s1s, &b_events_s1s);
    fChain->SetBranchAddress("s2s", &s2s, &b_events_s2s);
	
	
	Notify();
	
	fInit = true;
}


Bool_t PaxTreeWrapper::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


void PaxTreeWrapper::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}


Int_t PaxTreeWrapper::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

#endif
