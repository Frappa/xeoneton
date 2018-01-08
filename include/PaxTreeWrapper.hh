//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Dec  6 10:53:53 2017 by ROOT version 6.06/08
// from TTree tree/Tree with XENON1T events from pax
// found on file: /project2/lgrandi/xenon1t/processed/pax_v6.8.0/170807_1104.root
//////////////////////////////////////////////////////////

#ifndef PAXTREEWRAPPER_HH
#define PAXTREEWRAPPER_HH

#include "PaxClasses.hh"


// Header file for the classes stored in the TTree if any.
#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TObject.h"

#include <vector>



using namespace std;

class PaxTreeWrapper {
	
private:
	bool fInit;//True when the class is initialised
	TTree* fChain;   //!pointer to the analyzed TTree or TChain
	
	Int_t fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
	static const Int_t kMaxall_hits = 10000;
	static const Int_t kMaxinteractions = 20;
	static const Int_t kMaxpeaks = 10000;
	static const Int_t kMaxpulses = 10000;
	static const Int_t kMaxtrigger_signals = 10000;
	

public:
	
	// Declaration of leaf types
	// Declaration of leaf types
	//Event           *events;
    UInt_t          fUniqueID;
    UInt_t          fBits;
    Int_t           all_hits_;
    UInt_t          all_hits_fUniqueID[kMaxall_hits];   //[all_hits_]
    UInt_t          all_hits_fBits[kMaxall_hits];   //[all_hits_]
    Float_t         all_hits_area[kMaxall_hits];   //[all_hits_]
    Float_t         all_hits_center[kMaxall_hits];   //[all_hits_]
    Int_t           all_hits_channel[kMaxall_hits];   //[all_hits_]
    Int_t           all_hits_found_in_pulse[kMaxall_hits];   //[all_hits_]
    Float_t         all_hits_height[kMaxall_hits];   //[all_hits_]
    Int_t           all_hits_index_of_maximum[kMaxall_hits];   //[all_hits_]
    Bool_t          all_hits_is_rejected[kMaxall_hits];   //[all_hits_]
    Int_t           all_hits_left[kMaxall_hits];   //[all_hits_]
    Int_t           all_hits_left_central[kMaxall_hits];   //[all_hits_]
    Int_t           all_hits_n_saturated[kMaxall_hits];   //[all_hits_]
    Float_t         all_hits_noise_sigma[kMaxall_hits];   //[all_hits_]
    Int_t           all_hits_right[kMaxall_hits];   //[all_hits_]
    Int_t           all_hits_right_central[kMaxall_hits];   //[all_hits_]
    Float_t         all_hits_sum_absolute_deviation[kMaxall_hits];   //[all_hits_]
    Int_t           block_id;
    TString         dataset_name; //This branch cannot be retrieved and produces segmentation fault
    Int_t           event_number;
    Int_t           interactions_;
    UInt_t          interactions_fUniqueID[kMaxinteractions];   //[interactions_]
    UInt_t          interactions_fBits[kMaxinteractions];   //[interactions_]
    Float_t         interactions_drift_time[kMaxinteractions];   //[interactions_]
    Float_t         interactions_r_correction[kMaxinteractions];   //[interactions_]
    Int_t           interactions_s1[kMaxinteractions];   //[interactions_]
    Float_t         interactions_s1_area_correction[kMaxinteractions];   //[interactions_]
    Float_t         interactions_s1_area_fraction_top_probability[kMaxinteractions];   //[interactions_]
    Float_t         interactions_s1_pattern_fit[kMaxinteractions];   //[interactions_]
    Float_t         interactions_s1_pattern_fit_hits[kMaxinteractions];   //[interactions_]
    Float_t         interactions_s1_saturation_correction[kMaxinteractions];   //[interactions_]
    Float_t         interactions_s1_spatial_correction[kMaxinteractions];   //[interactions_]
    Int_t           interactions_s2[kMaxinteractions];   //[interactions_]
    Float_t         interactions_s2_area_correction[kMaxinteractions];   //[interactions_]
    Float_t         interactions_s2_lifetime_correction[kMaxinteractions];   //[interactions_]
    Float_t         interactions_x[kMaxinteractions];   //[interactions_]
    TString         interactions_xy_posrec_algorithm[kMaxinteractions];
    Float_t         interactions_xy_posrec_goodness_of_fit[kMaxinteractions];   //[interactions_]
    Float_t         interactions_xy_posrec_ndf[kMaxinteractions];   //[interactions_]
    Float_t         interactions_y[kMaxinteractions];   //[interactions_]
    Float_t         interactions_z[kMaxinteractions];   //[interactions_]
    Float_t         interactions_z_correction[kMaxinteractions];   //[interactions_]
    Bool_t          is_channel_suspicious[260];
    Short_t         lone_hits_per_channel[260];
    Short_t         lone_hits_per_channel_before[260];
    Int_t           n_channels;
    Short_t         n_hits_rejected[260];
    Int_t           n_pulses;
    Short_t         n_pulses_per_channel[260];
    Short_t         noise_pulses_in[260];
    Int_t           peaks_;
    UInt_t          peaks_fUniqueID[kMaxpeaks];   //[peaks_]
    UInt_t          peaks_fBits[kMaxpeaks];   //[peaks_]
    Float_t         peaks_area[kMaxpeaks];   //[peaks_]
    Double_t        peaks_area_decile_from_midpoint[kMaxpeaks][11];   //[peaks_]
    Float_t         peaks_area_fraction_top[kMaxpeaks];   //[peaks_]
    Float_t         peaks_area_midpoint[kMaxpeaks];   //[peaks_]
    Double_t        peaks_area_per_channel[kMaxpeaks][260];   //[peaks_]
    Float_t         peaks_birthing_split_fraction[kMaxpeaks];   //[peaks_]
    Float_t         peaks_birthing_split_goodness[kMaxpeaks];   //[peaks_]
    Float_t         peaks_bottom_hitpattern_spread[kMaxpeaks];   //[peaks_]
    Float_t         peaks_center_time[kMaxpeaks];   //[peaks_]
    TString         peaks_detector[kMaxpeaks];
    Float_t         peaks_height[kMaxpeaks];   //[peaks_]
    Float_t         peaks_hit_time_mean[kMaxpeaks];   //[peaks_]
    Float_t         peaks_hit_time_std[kMaxpeaks];   //[peaks_]
    vector<Hit>     peaks_hits[kMaxpeaks];
    Float_t         peaks_hits_fraction_top[kMaxpeaks];   //[peaks_]
    Short_t         peaks_hits_per_channel[kMaxpeaks][260];   //[peaks_]
    Int_t           peaks_index_of_maximum[kMaxpeaks];   //[peaks_]
    Float_t         peaks_interior_split_fraction[kMaxpeaks];   //[peaks_]
    Float_t         peaks_interior_split_goodness[kMaxpeaks];   //[peaks_]
    Float_t         peaks_largest_hit_area[kMaxpeaks];   //[peaks_]
    Int_t           peaks_largest_hit_channel[kMaxpeaks];   //[peaks_]
    Int_t           peaks_left[kMaxpeaks];   //[peaks_]
    Int_t           peaks_lone_hit_channel[kMaxpeaks];   //[peaks_]
    Float_t         peaks_mean_amplitude_to_noise[kMaxpeaks];   //[peaks_]
    Int_t           peaks_n_contributing_channels[kMaxpeaks];   //[peaks_]
    Int_t           peaks_n_contributing_channels_top[kMaxpeaks];   //[peaks_]
    Int_t           peaks_n_hits[kMaxpeaks];   //[peaks_]
    Int_t           peaks_n_noise_pulses[kMaxpeaks];   //[peaks_]
    Int_t           peaks_n_saturated_channels[kMaxpeaks];   //[peaks_]
    Short_t         peaks_n_saturated_per_channel[kMaxpeaks][260];   //[peaks_]
    Int_t           peaks_n_saturated_samples[kMaxpeaks];   //[peaks_]
    Double_t        peaks_range_area_decile[kMaxpeaks][11];   //[peaks_]
  //vector<ReconstructedPosition> peaks_reconstructed_positions[kMaxpeaks];
    Int_t           peaks_right[kMaxpeaks];   //[peaks_]
    Float_t         peaks_s2_bottom_spatial_correction[kMaxpeaks];   //[peaks_]
    Float_t         peaks_s2_saturation_correction[kMaxpeaks];   //[peaks_]
    Float_t         peaks_s2_spatial_correction[kMaxpeaks];   //[peaks_]
    Float_t         peaks_s2_top_spatial_correction[kMaxpeaks];   //[peaks_]
    Float_t         peaks_sum_waveform[kMaxpeaks][251];   //[peaks_]
    Float_t         peaks_sum_waveform_top[kMaxpeaks][251];   //[peaks_]
    Int_t           peaks_tight_coincidence[kMaxpeaks];   //[peaks_]
    Float_t         peaks_top_hitpattern_spread[kMaxpeaks];   //[peaks_]
    TString         peaks_type[kMaxpeaks];
    Int_t           pulses_;
    UInt_t          pulses_fUniqueID[kMaxpulses];   //[pulses_]
    UInt_t          pulses_fBits[kMaxpulses];   //[pulses_]
    Float_t         pulses_baseline[kMaxpulses];   //[pulses_]
    Float_t         pulses_baseline_increase[kMaxpulses];   //[pulses_]
    Int_t           pulses_channel[kMaxpulses];   //[pulses_]
    Int_t           pulses_hitfinder_threshold[kMaxpulses];   //[pulses_]
    Int_t           pulses_left[kMaxpulses];   //[pulses_]
    Float_t         pulses_maximum[kMaxpulses];   //[pulses_]
    Float_t         pulses_minimum[kMaxpulses];   //[pulses_]
    Int_t           pulses_n_hits_found[kMaxpulses];   //[pulses_]
    Float_t         pulses_noise_sigma[kMaxpulses];   //[pulses_]
    Int_t           pulses_right[kMaxpulses];   //[pulses_]
    Int_t           sample_duration;
    Long64_t        start_time;
    Long64_t        stop_time;
    Int_t           trigger_signals_;
    UInt_t          trigger_signals_fUniqueID[kMaxtrigger_signals];   //[trigger_signals_]
    UInt_t          trigger_signals_fBits[kMaxtrigger_signals];   //[trigger_signals_]
    Float_t         trigger_signals_area[kMaxtrigger_signals];   //[trigger_signals_]
    Int_t           trigger_signals_left_time[kMaxtrigger_signals];   //[trigger_signals_]
    Int_t           trigger_signals_n_contributing_channels[kMaxtrigger_signals];   //[trigger_signals_]
    Int_t           trigger_signals_n_pulses[kMaxtrigger_signals];   //[trigger_signals_]
    Int_t           trigger_signals_right_time[kMaxtrigger_signals];   //[trigger_signals_]
    Float_t         trigger_signals_time_mean[kMaxtrigger_signals];   //[trigger_signals_]
    Float_t         trigger_signals_time_rms[kMaxtrigger_signals];   //[trigger_signals_]
    Bool_t          trigger_signals_trigger[kMaxtrigger_signals];   //[trigger_signals_]
    Int_t           trigger_signals_type[kMaxtrigger_signals];   //[trigger_signals_]
    vector<int>     s1s;
    vector<int>     s2s;

	// List of branches
	TBranch        *b_events;
	TBranch        *b_events_fUniqueID;   //!
	TBranch        *b_events_fBits;   //!
	TBranch        *b_events_all_hits_;   //!
	TBranch        *b_all_hits_fUniqueID;   //!
	TBranch        *b_all_hits_fBits;   //!
	TBranch        *b_all_hits_area;   //!
	TBranch        *b_all_hits_center;   //!
	TBranch        *b_all_hits_channel;   //!
	TBranch        *b_all_hits_found_in_pulse;   //!
	TBranch        *b_all_hits_height;   //!
	TBranch        *b_all_hits_index_of_maximum;   //!
	TBranch        *b_all_hits_is_rejected;   //!
	TBranch        *b_all_hits_left;   //!
	TBranch        *b_all_hits_left_central;   //!
	TBranch        *b_all_hits_n_saturated;   //!
	TBranch        *b_all_hits_noise_sigma;   //!
	TBranch        *b_all_hits_right;   //!
	TBranch        *b_all_hits_right_central;   //!
	TBranch        *b_all_hits_sum_absolute_deviation;   //!
	TBranch        *b_events_block_id;   //!
	TBranch        *b_events_dataset_name;   //!
	TBranch        *b_events_event_number;   //!
	TBranch        *b_events_interactions_;   //!
	TBranch        *b_interactions_fUniqueID;   //!
	TBranch        *b_interactions_fBits;   //!
	TBranch        *b_interactions_drift_time;   //!
	TBranch        *b_interactions_r_correction;   //!
	TBranch        *b_interactions_s1;   //!
	TBranch        *b_interactions_s1_area_correction;   //!
	TBranch        *b_interactions_s1_area_fraction_top_probability;   //!
	TBranch        *b_interactions_s1_pattern_fit;   //!
	TBranch        *b_interactions_s1_pattern_fit_hits;   //!
	TBranch        *b_interactions_s1_saturation_correction;   //!
	TBranch        *b_interactions_s1_spatial_correction;   //!
	TBranch        *b_interactions_s2;   //!
	TBranch        *b_interactions_s2_area_correction;   //!
	TBranch        *b_interactions_s2_lifetime_correction;   //!
	TBranch        *b_interactions_x;   //!
	TBranch        *b_interactions_xy_posrec_algorithm;   //!
	TBranch        *b_interactions_xy_posrec_goodness_of_fit;   //!
	TBranch        *b_interactions_xy_posrec_ndf;   //!
	TBranch        *b_interactions_y;   //!
	TBranch        *b_interactions_z;   //!
	TBranch        *b_interactions_z_correction;   //!
	TBranch        *b_events_is_channel_suspicious;   //!
	TBranch        *b_events_lone_hits_per_channel;   //!
	TBranch        *b_events_lone_hits_per_channel_before;   //!
	TBranch        *b_events_n_channels;   //!
	TBranch        *b_events_n_hits_rejected;   //!
	TBranch        *b_events_n_pulses;   //!
	TBranch        *b_events_n_pulses_per_channel;   //!
	TBranch        *b_events_noise_pulses_in;   //!
	TBranch        *b_events_peaks_;   //!
	TBranch        *b_peaks_fUniqueID;   //!
	TBranch        *b_peaks_fBits;   //!
	TBranch        *b_peaks_area;   //!
	TBranch        *b_peaks_area_decile_from_midpoint;   //!
	TBranch        *b_peaks_area_fraction_top;   //!
	TBranch        *b_peaks_area_midpoint;   //!
	TBranch        *b_peaks_area_per_channel;   //!
	TBranch        *b_peaks_birthing_split_fraction;   //!
	TBranch        *b_peaks_birthing_split_goodness;   //!
	TBranch        *b_peaks_bottom_hitpattern_spread;   //!
	TBranch        *b_peaks_center_time;   //!
	TBranch        *b_peaks_detector;   //!
	TBranch        *b_peaks_height;   //!
	TBranch        *b_peaks_hit_time_mean;   //!
	TBranch        *b_peaks_hit_time_std;   //!
	TBranch        *b_peaks_hits;   //!
	TBranch        *b_peaks_hits_fraction_top;   //!
	TBranch        *b_peaks_hits_per_channel;   //!
	TBranch        *b_peaks_index_of_maximum;   //!
	TBranch        *b_peaks_interior_split_fraction;   //!
	TBranch        *b_peaks_interior_split_goodness;   //!
	TBranch        *b_peaks_largest_hit_area;   //!
	TBranch        *b_peaks_largest_hit_channel;   //!
	TBranch        *b_peaks_left;   //!
	TBranch        *b_peaks_lone_hit_channel;   //!
	TBranch        *b_peaks_mean_amplitude_to_noise;   //!
	TBranch        *b_peaks_n_contributing_channels;   //!
	TBranch        *b_peaks_n_contributing_channels_top;   //!
	TBranch        *b_peaks_n_hits;   //!
	TBranch        *b_peaks_n_noise_pulses;   //!
	TBranch        *b_peaks_n_saturated_channels;   //!
	TBranch        *b_peaks_n_saturated_per_channel;   //!
	TBranch        *b_peaks_n_saturated_samples;   //!
	TBranch        *b_peaks_range_area_decile;   //!
	TBranch        *b_peaks_reconstructed_positions;   //!
	TBranch        *b_peaks_right;   //!
	TBranch        *b_peaks_s2_bottom_spatial_correction;   //!
	TBranch        *b_peaks_s2_saturation_correction;   //!
	TBranch        *b_peaks_s2_spatial_correction;   //!
	TBranch        *b_peaks_s2_top_spatial_correction;   //!
	TBranch        *b_peaks_sum_waveform;   //!
	TBranch        *b_peaks_sum_waveform_top;   //!
	TBranch        *b_peaks_tight_coincidence;   //!
	TBranch        *b_peaks_top_hitpattern_spread;   //!
	TBranch        *b_peaks_type;   //!
	TBranch        *b_events_pulses_;   //!
	TBranch        *b_pulses_fUniqueID;   //!
	TBranch        *b_pulses_fBits;   //!
	TBranch        *b_pulses_baseline;   //!
	TBranch        *b_pulses_baseline_increase;   //!
	TBranch        *b_pulses_channel;   //!
	TBranch        *b_pulses_hitfinder_threshold;   //!
	TBranch        *b_pulses_left;   //!
	TBranch        *b_pulses_maximum;   //!
	TBranch        *b_pulses_minimum;   //!
	TBranch        *b_pulses_n_hits_found;   //!
	TBranch        *b_pulses_noise_sigma;   //!
	TBranch        *b_pulses_right;   //!
	TBranch        *b_events_sample_duration;   //!
	TBranch        *b_events_start_time;   //!
	TBranch        *b_events_stop_time;   //!
	TBranch        *b_events_trigger_signals_;   //!
	TBranch        *b_trigger_signals_fUniqueID;   //!
	TBranch        *b_trigger_signals_fBits;   //!
	TBranch        *b_trigger_signals_area;   //!
	TBranch        *b_trigger_signals_left_time;   //!
	TBranch        *b_trigger_signals_n_contributing_channels;   //!
	TBranch        *b_trigger_signals_n_pulses;   //!
	TBranch        *b_trigger_signals_right_time;   //!
	TBranch        *b_trigger_signals_time_mean;   //!
	TBranch        *b_trigger_signals_time_rms;   //!
	TBranch        *b_trigger_signals_trigger;   //!
	TBranch        *b_trigger_signals_type;   //!
	TBranch        *b_events_s1s;   //!
	TBranch        *b_events_s2s;   //!


	PaxTreeWrapper(TTree *tree=0);
	PaxTreeWrapper(string filename);
	virtual ~PaxTreeWrapper();
	
	virtual Int_t    Cut(Long64_t entry);
	virtual Int_t    GetEntry(Long64_t entry);
	virtual Long64_t LoadTree(Long64_t entry);
	virtual void     Init(TTree *tree);
	virtual void     Loop();
	virtual Bool_t   Notify();
	virtual void     Show(Long64_t entry = -1);
	
	TTree* GetTree(){return fChain;};
	bool IsInit(){return fInit;};
};

#endif

#ifdef paxtreewrapper_cc
#endif // #ifdef paxtreewrapper_cc
