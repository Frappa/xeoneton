//#include "xe1t_analysis.hh"
#include "PaxClasses.hh"
#include "PaxTreeWrapper.hh"

#include "TString.h"
#include "TTree.h"
#include "TChain.h"
//#include "TCanvas.h"
//#include "TApplication.h"
#include "TF1.h"
#include "TColor.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TPaveStats.h"
#include "TText.h"
#include "TLatex.h"
#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TError.h"
//#include "TCint.h"
#include "TInterpreter.h"

//#ifdef AP_MPI
//	#include <mpi.h>
//#endif 

#include <sstream>      
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "math.h"
#include <sys/time.h> // for gettimeofday()
#include <unistd.h>


using namespace std;

void MakeMiniTree(string MiniTreeName);

namespace minitrees
{
	//This class is used to address the branches of the TTree
	class EventData{
	public:
		EventData(){
			S1s_vec= new vector<vector<double> >;
			S1sTot_vec = new vector<double>;
			S1sTop_vec = new vector<double>;
			S1sBot_vec = new vector<double>;
			S1sPeakMaxPos_vec = new vector<int>;
			S1sPeakMax_vec = new vector<double>;
			S1sMidPoint_vec = new vector<double>;
			S1sCenterTime_vec = new vector<double>;
			S1sWidth_vec = new vector<double>;
			S1sLowWidth_vec = new vector<double>;
			S1sNcoins_vec = new vector<int>;
			S1sNcoinsTop_vec = new vector<int>;
			S1sNcoinsBot_vec = new vector<int>;
			S1sCoin_vec = new vector<int>;
			S1sCoinTop_vec = new vector<int>;
			S1sCoinBot_vec = new vector<int>;
			S1sPmtOrder_vec= new vector<vector<int> >;
			
			S2s_vec= new vector<vector<double> >;
			S2sTot_vec = new vector<double>;
			S2sTop_vec = new vector<double>;
			S2sBot_vec = new vector<double>;
			S2sPeakMaxPos_vec = new vector<int>;
			S2sPeakMax_vec = new vector<double>;
			S2sMidPoint_vec = new vector<double>;
			S2sCenterTime_vec = new vector<double>;
			S2sWidth_vec = new vector<double>;
			S2sLowWidth_vec = new vector<double>;
			S2sNcoins_vec = new vector<int>;
			S2sNcoinsTop_vec = new vector<int>;
			S2sNcoinsBot_vec = new vector<int>;
			S2sCoin_vec = new vector<int>;
			S2sCoinTop_vec = new vector<int>;
			S2sCoinBot_vec = new vector<int>;
			S2sPmtOrder_vec= new vector<vector<int> >;
		};
		
		void ResetS1peaks(int nS1peaks)
		{
			S1s_vec->resize(nS1peaks);
			S1sTot_vec->resize(nS1peaks);
			S1sTop_vec->resize(nS1peaks);
			S1sBot_vec->resize(nS1peaks);
			S1sPeakMaxPos_vec->resize(nS1peaks);
			S1sPeakMax_vec->resize(nS1peaks);
			S1sMidPoint_vec->resize(nS1peaks);
			S1sCenterTime_vec->resize(nS1peaks);
			S1sWidth_vec->resize(nS1peaks);
			S1sLowWidth_vec->resize(nS1peaks);
			S1sNcoins_vec->resize(nS1peaks);
			S1sNcoinsTop_vec->resize(nS1peaks);
			S1sNcoinsBot_vec->resize(nS1peaks);
			S1sCoin_vec->resize(nS1peaks);
			S1sCoinTop_vec->resize(nS1peaks);
			S1sCoinBot_vec->resize(nS1peaks);
			S1sPmtOrder_vec->resize(nS1peaks);
		};
		
		void ResetS2peaks(int nS2peaks)
		{
			S2s_vec->resize(nS2peaks);
			S2sTot_vec->resize(nS2peaks);
			S2sTop_vec->resize(nS2peaks);
			S2sBot_vec->resize(nS2peaks);
			S2sPeakMaxPos_vec->resize(nS2peaks);
			S2sPeakMax_vec->resize(nS2peaks);
			S2sMidPoint_vec->resize(nS2peaks);
			S2sCenterTime_vec->resize(nS2peaks);
			S2sWidth_vec->resize(nS2peaks);
			S2sLowWidth_vec->resize(nS2peaks);
			S2sNcoins_vec->resize(nS2peaks);
			S2sNcoinsTop_vec->resize(nS2peaks);
			S2sNcoinsBot_vec->resize(nS2peaks);
			S2sCoin_vec->resize(nS2peaks);
			S2sCoinTop_vec->resize(nS2peaks);
			S2sCoinBot_vec->resize(nS2peaks);
			S2sPmtOrder_vec->resize(nS2peaks);
		};
		
		
		vector<vector<double> > *S1s_vec;
		vector<double> *S1sTot_vec;
		vector<double> *S1sTop_vec;
		vector<double> *S1sBot_vec;
		vector<int> *S1sPeakMaxPos_vec;
		vector<double> *S1sPeakMax_vec;
		vector<double> *S1sMidPoint_vec;
		vector<double> *S1sCenterTime_vec;
		vector<double> *S1sWidth_vec;
		vector<double> *S1sLowWidth_vec;
		vector<int> *S1sNcoins_vec;
		vector<int> *S1sNcoinsTop_vec;
		vector<int> *S1sNcoinsBot_vec;
		vector<int> *S1sCoin_vec;
		vector<int> *S1sCoinTop_vec;
		vector<int> *S1sCoinBot_vec;
		vector<vector<int> > *S1sPmtOrder_vec;
	
		vector<vector<double> > *S2s_vec;
		vector<double> *S2sTot_vec;
		vector<double> *S2sTop_vec;
		vector<double> *S2sBot_vec;
		vector<int> *S2sPeakMaxPos_vec;
		vector<double> *S2sPeakMax_vec;
		vector<double> *S2sMidPoint_vec;
		vector<double> *S2sCenterTime_vec;
		vector<double> *S2sWidth_vec;
		vector<double> *S2sLowWidth_vec;
		vector<int> *S2sNcoins_vec;
		vector<int> *S2sNcoinsTop_vec;
		vector<int> *S2sNcoinsBot_vec;
		vector<int> *S2sCoin_vec;
		vector<int> *S2sCoinTop_vec;
		vector<int> *S2sCoinBot_vec;
		vector<vector<int> > *S2sPmtOrder_vec;
	};
	
	
	class PmtContrib
	{
	public:
		PmtContrib(int _index, double _area){
			index = _index;
			area = _area;
		}
		
		int index;
		double area;
		
		inline bool operator< (const PmtContrib& rhs) const{
			return (area<rhs.area);
		};
		
		inline bool operator> (const PmtContrib& rhs) const{
			return (area>rhs.area);
		};
	};
	
	
	
	
	//Functions
	void MakeBranches(TTree* mt);
	void ProcessEvents(TTree *mt);
	void ExtractS1info(Peak *peak, int peakindex);
	void ExtractS2info(Peak *peak, int peakindex);
	void ReadS1peakData(int iPeak);
	void ReadS2peakData(int iPeak);
	
	
	//Variables
	
	EventData data; //This contains the variables to be put in the new minitree
	
	PaxTreeWrapper *paxdata = NULL; //This is used to set the branch addresses of the pax tree
	
	set<string> datasets;
	
	const string procdir = "/home/fpiastra/xenon/processed/";
	
	const string my_minitreedir = "/home/fpiastra/xenon/minitrees_local/";
	
	
}//End of namespace minitrees


//using namespace minitrees;
void MakeMiniTree(string MiniTreeName)
{
	if(!minitrees::datasets.size()) return;
	
	
	//Load the std classes
	gROOT->ProcessLine("#include <vector>");
	
	//Load the PaxClasses (Events, etc.)
	//gInterpreter->AddIncludePath("/home/fpiastra/xenon/analysis/code/include");
	//gROOT->ProcessLine("#include \"PaxClasses.hh\"");
	
	if(gSystem->Load("libPaxClasses") < 0)
	{
		cerr << "\nERROR: Could not find the library \"libPaxClasses\"" << endl << endl;
		exit(-1);
	}
	
	
	set<string>::iterator It;
	for(It=minitrees::datasets.begin(); It!=minitrees::datasets.end(); It++)
	{
		string infilename = minitrees::procdir + (*It) + string(".root");
	
		//Control that the file exists
		if( access(infilename.c_str(), R_OK|X_OK)!=0 )
		{
			cerr << "\nWARNING: Cannot open or find input file <" << infilename << ">" << endl;
			continue;
		}
	
		TFile *inFile = TFile::Open(infilename.c_str());
	
		if( (!inFile)||(!inFile->IsOpen()) ){
			continue;
		}
	
		TTree *tree = (TTree*)inFile->Get("tree");
	
		if( !tree )
		{
			cerr << "\nWARNING: Cannot open or find the TTree \"tree\" inside the input file <" << infilename << ">" << endl;
			delete inFile;
			continue;
		}
		
		if(!minitrees::paxdata)
		{
			minitrees::paxdata = new PaxTreeWrapper(tree);
		}else{
			minitrees::paxdata->Init(tree);
		}
		
		
		if( !minitrees::paxdata->IsInit() )
		{
			cerr << "\nERROR: Could not initialise the PaxTreeWrapper object with the TTree \"tree\" from file <" << infilename << ">" << endl;
			delete inFile;
			if(minitrees::paxdata)
			{
				delete minitrees::paxdata;
				minitrees::paxdata = NULL;
			}
			continue;
		}
		
		
		string outfilename = minitrees::my_minitreedir + (*It) + string("_") + MiniTreeName + string(".root");
	
		TFile *outFile = TFile::Open( outfilename.c_str(), "recreate" );
	
		if( (!outFile)||(!outFile->IsOpen()) )
		{
			delete inFile;
			tree = NULL;
			continue;
		}
		
		
		TTree *mt = NULL;
		if(!mt) mt = new TTree( MiniTreeName.c_str(), "" );
		
		minitrees::MakeBranches(mt);
		
		minitrees::ProcessEvents(mt);
		
		
	
		outFile->WriteTObject(mt, 0, "overwrite");
	
		delete outFile;
		mt = NULL;
		
	}//End of loop over datasets
}//End of mitrees namespace



//SERVICE FUNCTIONS
void AddDataset(string dsname)
{
	minitrees::datasets.insert(dsname);
};


void ResetDataSets()
{
	minitrees::datasets.clear();
};


void help();




int main(int argc, char **argv)
{
	if(argc<3)
	{
		help();
		return -1;
	}
	
	ResetDataSets();
	
	for(int iArg=2; iArg<argc; iArg++)
	{
		AddDataset(argv[iArg]);
	}
	
	MakeMiniTree(argv[1]);
	
	return 0;
}



void help()
{
	cout << "\n\nUsage:" << endl << endl;
	cout << "MakeMiniTree <MiniTreeName> <ListOfDatasets>" << endl << endl;
	cout << "The list of datasets should be separated by a simple space." << endl << endl;
};





//Definition of the functions
namespace minitrees
{
	void MakeBranches(TTree* mt)
	{
		if(!mt){
			cerr << "\nERROR: Null pointer to the minitree" << endl;
			return;
		}
		
		//Define the branches of the minitree
		mt->Branch("S1s", "vector<vector<double> >", &data.S1s_vec);
		mt->Branch("S1sTot", "vector<double>", &data.S1sTot_vec);
		mt->Branch("S1sTop", "vector<double>", &data.S1sTop_vec);
		mt->Branch("S1sBot", "vector<double>", &data.S1sBot_vec);
		mt->Branch("S1sPeak", "vector<int>", &data.S1sPeakMaxPos_vec);
		mt->Branch("S1sPeakMax", "vector<double>", &data.S1sPeakMax_vec);
		mt->Branch("S1sMidPoint", "vector<double>", &data.S1sMidPoint_vec);
		mt->Branch("S1sCenterTime", "vector<double>", &data.S1sCenterTime_vec);
		mt->Branch("S1sWidth", "vector<double>", &data.S1sWidth_vec);
		mt->Branch("S1sLowWidth", "vector<double>", &data.S1sLowWidth_vec);
		mt->Branch("S1sCoins", "vector<int>", &data.S1sNcoins_vec);
		mt->Branch("S1sCoinsTop", "vector<int>", &data.S1sNcoinsTop_vec);
		mt->Branch("S1sCoinsBot", "vector<int>", &data.S1sNcoinsBot_vec);
		mt->Branch("S1sPmtOrder", "vector<vector<int> >", &data.S1sPmtOrder_vec);
	
		mt->Branch("S2s", "vector<vector<double> >", &data.S2s_vec);
		mt->Branch("S2sTot", "vector<double>", &data.S2sTot_vec);
		mt->Branch("S2sTop", "vector<double>", &data.S2sTop_vec);
		mt->Branch("S2sBot", "vector<double>", &data.S2sBot_vec);
		mt->Branch("S2sPeak", "vector<int>", &data.S2sPeakMaxPos_vec);
		mt->Branch("S2sPeakMax", "vector<double>", &data.S2sPeakMax_vec);
		mt->Branch("S2sMidPoint", "vector<double>", &data.S2sMidPoint_vec);
		mt->Branch("S2sCenterTime", "vector<double>", &data.S2sCenterTime_vec);
		mt->Branch("S2sWidth", "vector<double>", &data.S2sWidth_vec);
		mt->Branch("S2sLowWidth", "vector<double>", &data.S2sLowWidth_vec);
		mt->Branch("S2sCoins", "vector<int>", &data.S2sNcoins_vec);
		mt->Branch("S2sCoinsTop", "vector<int>", &data.S2sNcoinsTop_vec);
		mt->Branch("S2sCoinsBot", "vector<int>", &data.S2sNcoinsBot_vec);
		mt->Branch("S2sPmtOrder", "vector<vector<int> >", &data.S2sPmtOrder_vec);
		
		/*
		//Not well recognised peaks 
		mt->Branch("SXs", "vector<vector<double> >", &data.S1s_vec);
		mt->Branch("SXsTot", "vector<double>", &data.S1sTot_vec);
		mt->Branch("SXsTop", "vector<double>", &data.S1sTop_vec);
		mt->Branch("SXsBot", "vector<double>", &data.S1sBot_vec);
		mt->Branch("SXsPeak", "vector<int>", &data.S1sPeakMaxPos_vec);
		mt->Branch("SXsPeakMax", "vector<double>", &data.S1sPeakMax_vec);
		mt->Branch("SXsMidPoint", "vector<double>", &data.S1sMidPoint_vec);
		mt->Branch("SXsCenterTime", "vector<double>", &data.S1sCenterTime_vec);
		mt->Branch("SXsWidth", "vector<double>", &data.S1sWidth_vec);
		mt->Branch("SXsLowWidth", "vector<double>", &data.S1sLowWidth_vec);
		mt->Branch("SXsCoins", "vector<int>", &data.S1sNcoins_vec);
		mt->Branch("SXsCoinsTop", "vector<int>", &data.S1sNcoinsTop_vec);
		mt->Branch("SXsCoinsBot", "vector<int>", &data.S1sNcoinsBot_vec);
		mt->Branch("SXsPmtOrder", "vector<vector<int> >", &data.S1sPmtOrder_vec);
		*/
	}
	
	
	void ProcessEvents(TTree *mt)
	{
		if(!paxdata)
		{
			cerr << "\nProcessEvents: ERROR: The PaxTreeWrapper pointer is zero!" << endl << endl;
			return;
		}
		
		if( !paxdata->IsInit() )
		{
			cerr << "\nProcessEvents: ERROR: The PaxTreeWrapper is not initialised!" << endl << endl;
			return;
		}
		
		
		TTree *tree = paxdata->GetTree();
		if(!tree)
		{
			cerr << "\nProcessEvents: ERROR: \"tree\" pointer is 0" << endl << endl;
			return;
		}
		
		if(!mt)
		{
			cerr << "\nProcessEvents: ERROR: \"mt\" pointer is 0" << endl << endl;
			return;
		}
		
		tree->SetBranchStatus("*", 0);
		
		tree->SetBranchStatus("all_hits", 1);
		
		tree->SetBranchStatus("peaks", 1);
		tree->SetBranchStatus("peaks.area", 1);
		tree->SetBranchStatus("peaks.area_fraction_top", 1);
		tree->SetBranchStatus("peaks.area_midpoint", 1);
		tree->SetBranchStatus("peaks.area_per_channel[260]", 1);
		tree->SetBranchStatus("peaks.center_time", 1);
		tree->SetBranchStatus("peaks.height", 1);
		tree->SetBranchStatus("peaks.index_of_maximum", 1);
		tree->SetBranchStatus("peaks.n_contributing_channels", 1);
		tree->SetBranchStatus("peaks.n_contributing_channels_top", 1);
		tree->SetBranchStatus("peaks.range_area_decile[11]", 1);
		tree->SetBranchStatus("peaks.s2_bottom_spatial_correction", 1);
		tree->SetBranchStatus("peaks.s2_spatial_correction", 1);
		tree->SetBranchStatus("peaks.s2_top_spatial_correction", 1);
		tree->SetBranchStatus("s1s", 1);
		tree->SetBranchStatus("s2s", 1);
		
		
		
		int nEvents = tree->GetEntries();
		
		
		for(int iEv=0; iEv<nEvents; iEv++)
		//for(int iEv=0; iEv<500; iEv++)
		{
			Long64_t nbytes = tree->GetEntry(iEv);
			if(!(nbytes>0))
			{
				cerr << "\nProcessEvents: WARNING: Read 0 bytes at event #" << iEv << endl << endl;
				continue;
			}
			
			
			if(false){//This is for debug only
				if(paxdata->peaks_)
				{
					std::cout << "\nEntry " << iEv << std::endl;
					std::cout << "peaks_ = " << paxdata->peaks_ << std::endl;
					std::cout << "nS1s = " << (paxdata->s1s).size() << std::endl;
					std::cout << "nS2s = " << (paxdata->s2s).size() << std::endl;
					//std::cout << "peaks_ = " << tw->peaks_ << std::endl;
					//std::cout << "peaks_ = " << tw->peaks_ << std::endl;
					//std::cout << "peaks_ = " << tw->peaks_ << std::endl;
				}
				continue;
			}
			
			int nPeaks = paxdata->peaks_;
			int nS1s=(paxdata->s1s).size();
			int nS2s=(paxdata->s2s).size();
			
			if(nPeaks<(nS1s+nS2s)){//Something is clearly wrong!
				cerr << "\nProcessEvents: WARNING: nPeaks=" << nPeaks << ", while (nS1s+nS2s)=" << (nS1s+nS2s) << " at event #" << iEv << " " << endl << endl;
				//continue;
			}
			
			//Uncomment for debugging
			//cout << "Debug--> ProcessEvents:   Event=" << iEv << "   nPeaks=" << nPeaks << ",   nS1s=" << nS1s << ",   nS2s=" << nS2s << endl;
			
			
			//continue; //Uncomment for debug purpose.
			
			data.ResetS1peaks(0);
			data.ResetS2peaks(0);
			
			
			//Cycle over the S1s
			for(int iS1=0; iS1<nS1s; iS1++)
			{
				int s1index = (paxdata->s1s).at(iS1);
				/*
				if(paxdata.peaks_type[s1index] != TString("s1") )
				{
					cerr << "\nProcessEvents: WARNING: read peak type \"" << paxdata.peaks_type[s1index] << "\" for an expected S1 signal. Peak index: " << s1index << ". Event: " << iEv << endl << endl;
					continue;
				}
				*/
				ReadS1peakData(s1index);
			}
			
			//Cycle over the S2s
			for(int iS2=0; iS2<nS2s; iS2++)
			{
				int s2index = (paxdata->s2s).at(iS2);
				/*
				if(paxdata.peaks_type[s2index] != TString("s2") )
				{
					cerr << "\nProcessEvents: WARNING: read peak type \"" << paxdata.peaks_type[s2index] << "\" for an expected S2 signal. Peak index: " << s2index << ". Event: " << iEv << endl << endl;
					continue;
				}
				*/
				ReadS2peakData(s2index);
			}
			
			
			mt->Fill();
			
		}//End of loop over the tree entries
		
	}
	
	
	void ReadS1peakData(int iPeak)
	{
		
		//Uncomment for debugging
		//cout << "Debug--> ReadS1peakData:   Reading data of peak " << iPeak << endl;
		//cout << "Debug--> ReadS1peakData:   Vector size= " << data.S2sTot_vec->size() << endl;
		
		
		data.S1sTot_vec->push_back( (paxdata->peaks_area)[iPeak] );
		data.S1sTop_vec->push_back( (paxdata->peaks_area)[iPeak] * (paxdata->peaks_area_fraction_top)[iPeak] );
		data.S1sBot_vec->push_back( (paxdata->peaks_area)[iPeak] * (1.-(paxdata->peaks_area_fraction_top)[iPeak]) );
		data.S1sPeakMaxPos_vec->push_back( (paxdata->peaks_index_of_maximum)[iPeak] );
		data.S1sPeakMax_vec->push_back( (paxdata->peaks_height)[iPeak] );
		data.S1sCenterTime_vec->push_back( (paxdata->peaks_center_time)[iPeak] );
		data.S1sWidth_vec->push_back( (paxdata->peaks_range_area_decile)[iPeak][5] ); //Replaces FWHM
		data.S1sLowWidth_vec->push_back( (paxdata->peaks_range_area_decile)[iPeak][9] ); //Replaces FWTM
		data.S1sNcoins_vec->push_back( (paxdata->peaks_n_contributing_channels)[iPeak] );
		data.S1sNcoinsTop_vec->push_back( (paxdata->peaks_n_contributing_channels_top)[iPeak] );
		data.S1sNcoinsBot_vec->push_back( (paxdata->peaks_n_contributing_channels)[iPeak] - (paxdata->peaks_n_contributing_channels_top)[iPeak] );
		
		//Uncomment for debugging
		//cout << "Debug--> ReadS1peakData:   Reading coincidences for peak " << iPeak << endl;
		
		//Contribution from each channel to the total area
		vector<double> s1s( &((paxdata->peaks_area_per_channel)[iPeak][0]), &((paxdata->peaks_area_per_channel)[iPeak][248]) );
		
		//This vector of structures is only used for sorting the peak indexex from the largest to the smallest S1
		vector<PmtContrib> pmtscontrib;
		for(int iCh=0; iCh<s1s.size(); iCh++)
		{
			if(s1s.at(iCh)>0)
			{
				pmtscontrib.push_back( PmtContrib( iCh, s1s.at(iCh)) );
			}
		}
		
		//Uncomment for debugging
		//cout << "Debug--> ReadS1peakData:   Sorting pmtscontrib for peak " << iPeak << endl;
		
		//SORT THE "pmtscontrib" by area and save the sorted pmt indexes
		sort( pmtscontrib.begin(), pmtscontrib.end(), greater<PmtContrib>() );
		vector<int> S1sPmtOrder( pmtscontrib.size() );
		for(int iPmt=0; iPmt<pmtscontrib.size(); iPmt++)
		{
			S1sPmtOrder.at(iPmt) = pmtscontrib.at(iPmt).index;
		}
		
		
		data.S1s_vec->push_back( s1s );
	}
	
	
	void ReadS2peakData(int iPeak)
	{
		//Uncomment for debugging
		//cout << "Debug--> ReadS2peakData:   Reading data of of peak " << iPeak << endl;
		//cout << "Debug--> ReadS2peakData:   Vector size= " << data.S2sTot_vec->size() << endl;
		
		data.S2sTot_vec->push_back( (paxdata->peaks_area)[iPeak] );
		data.S2sTop_vec->push_back( (paxdata->peaks_area)[iPeak] * (paxdata->peaks_area_fraction_top)[iPeak] );
		data.S2sBot_vec->push_back( (paxdata->peaks_area)[iPeak] * (1.-(paxdata->peaks_area_fraction_top)[iPeak]) );
		data.S2sPeakMaxPos_vec->push_back( (paxdata->peaks_index_of_maximum)[iPeak] );
		data.S2sPeakMax_vec->push_back( (paxdata->peaks_height)[iPeak] );
		data.S2sCenterTime_vec->push_back( (paxdata->peaks_center_time)[iPeak] );
		data.S2sWidth_vec->push_back( (paxdata->peaks_range_area_decile)[iPeak][5] ); //Replaces FWHM
		data.S2sLowWidth_vec->push_back( (paxdata->peaks_range_area_decile)[iPeak][9] ); //Replaces FWTM
		data.S2sNcoins_vec->push_back( (paxdata->peaks_n_contributing_channels)[iPeak] );
		data.S2sNcoinsTop_vec->push_back( (paxdata->peaks_n_contributing_channels_top)[iPeak] );
		data.S2sNcoinsBot_vec->push_back( (paxdata->peaks_n_contributing_channels)[iPeak] - (paxdata->peaks_n_contributing_channels_top)[iPeak] );
		
		//Uncomment for debugging
		//cout << "Debug--> ReadS2peakData:   Reading coincidences for peak " << iPeak << endl;
		
		//Contribution from each channel to the total area
		vector<double> s2s( &((paxdata->peaks_area_per_channel)[iPeak][0]), &((paxdata->peaks_area_per_channel)[iPeak][248]) );
		
		//This vector of structures is only used for sorting the peak indexex from the largest to the smallest S1
		vector<PmtContrib> pmtscontrib;
		for(int iCh=0; iCh<s2s.size(); iCh++)
		{
			if(s2s.at(iCh)>0)
			{
				pmtscontrib.push_back( PmtContrib( iCh, s2s.at(iCh)) );
			}
		}
		
		//Uncomment for debugging
		//cout << "Debug--> ReadS2peakData:   Sorting pmtscontrib for peak " << iPeak << endl;
		
		//SORT THE "pmtscontrib" by area and save the sorted pmt indexes
		sort( pmtscontrib.begin(), pmtscontrib.end(), greater<PmtContrib>() );
		vector<int> S2sPmtOrder(pmtscontrib.size());
		for(int iPmt=0; iPmt<pmtscontrib.size(); iPmt++)
		{
			S2sPmtOrder.at(iPmt) = pmtscontrib.at(iPmt).index;
		}
		
		
		data.S2s_vec->push_back( s2s );
		
	}
}