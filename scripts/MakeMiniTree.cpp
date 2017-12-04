#include "xe1t_analysis.hh"
#include "PaxTree.hh"

#include "TString.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TF1.h"
#include "TColor.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TPaveStats.h"
#include "TText.h"
#include "TLatex.h"
#include "TFile.h"
#include "TError.h"

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
	set<string> datasets;
	
	const string procdir = "/home/fpiastra/xenon/processed/";
	
	const string my_minitreedir = "/home/fpiastra/xenon/minitrees_local/";
	
	
	void MakeBranches(TTree* mt);
	
	void ProcessEvent(Event *ev);
	void ExtractS1info(Peak *peak, int peakindex);
	void ExtractS2info(Peak *peak, int peakindex);
	
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
			S1sCenterTime_vec = new vector<double>;
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
			S2sCenterTime_vec = new vector<double>;
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
			S1sCenterTime_vec->resize(nS1peaks);
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
			S2sCenterTime_vec->resize(nS2peaks);
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
		vector<double> *S1sCenterTime_vec;
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
		vector<double> *S2sCenterTime_vec;
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
	
	
	
	EventData data;
	
	
	void MakeBrances(TTree* mt)
	{
		if(!mt) return;
		
		//Define the branches of the minitree
		mt->Branch("S1s", "vector<vector<double> >", &data.S1s_vec);
		mt->Branch("S1sTot", "vector<double>", &data.S1sTot_vec);
		mt->Branch("S1sTop", "vector<double>", &data.S1sTop_vec);
		mt->Branch("S1sBot", "vector<double>", &data.S1sBot_vec);
		mt->Branch("S1sPeak", "vector<int>", &data.S1sPeakMaxPos_vec);
		mt->Branch("S1sPeakMax", "vector<double>", &data.S1sPeakMax_vec);
		mt->Branch("S1sCenterTime", "vector<double>", &data.S1sCenterTime_vec);
		mt->Branch("S1sCoins", "vector<int>", &data.S1sNcoins_vec);
		mt->Branch("S1sCoinsTop", "vector<int>", &data.S1sNcoinsTop_vec);
		mt->Branch("S1sCoinsBot", "vector<int>", &data.S1sNcoinsBot_vec);
		//mt->Branch("S1sCoinTop", "vector<int>", &data.S1sCoinTop_vec);
		//mt->Branch("S1sCoinBot", "vector<int>", &data.S1sCoinBot_vec);
		mt->Branch("S1sPmtOrder", "vector<vector<int> >", &data.S1sPmtOrder_vec);
	
		mt->Branch("S2s", "vector<vector<double> >", &data.S2s_vec);
		mt->Branch("S2sTot", "vector<double>", &data.S2sTot_vec);
		mt->Branch("S2sTop", "vector<double>", &data.S2sTop_vec);
		mt->Branch("S2sBot", "vector<double>", &data.S2sBot_vec);
		mt->Branch("S2sPeak", "vector<int>", &data.S2sPeakMaxPos_vec);
		mt->Branch("S2sPeakMax", "vector<double>", &data.S1sPeakMax_vec);
		mt->Branch("S2sCenterTime", "vector<double>", &data.S2sCenterTime_vec);
		mt->Branch("S2sCoins", "vector<int>", &data.S2sNcoins_vec);
		mt->Branch("S2sCoinsTop", "vector<int>", &data.S2sNcoinsTop_vec);
		mt->Branch("S2sCoinsBot", "vector<int>", &data.S2sNcoinsBot_vec);
		//mt->Branch("S2sCoinTop", "vector<int>", &data.S2sCoin_vec);
		//mt->Branch("S2sCoinBot", "vector<int>", &data.S2sCoin_vec);
		mt->Branch("S2sPmtOrder", "vector<vector<int> >", &data.S2sPmtOrder_vec);
	};
	
	
	void ProcessEvent(Event *ev)
	{
		if(!ev) return;
		
		int nPeaks = (ev->peaks).size();
		int nS1s = (ev->s1s).size();
		int nS2s = (ev->s2s).size();
		
		data.ResetS1peaks(nS1s);
		data.ResetS2peaks(nS2s);
		
		if(nPeaks<nS1s+nS2s){//Something is clearly wrong!
			return;
		}
		
		
		//Cycle over the S1s
		int counter = 0;
		vector<Int_t>::iterator ItS1;
		for(ItS1=(ev->s1s).begin(); ItS1!=(ev->s1s).end(); ItS1++)
		{
			Peak peak = (ev->peaks).at(*ItS1);
			
			if(peak.type != TString("s1") )
			{
				counter++;
				continue;
			}
			
			ExtractS1info(&peak, counter);
			counter++;
		}
		
		
		//Cycle over the S1s
		counter = 0;
		vector<Int_t>::iterator ItS2;
		for(ItS2=(ev->s2s).begin(); ItS2!=(ev->s2s).end(); ItS2++)
		{
			Peak peak = (ev->peaks).at(*ItS2);
			if(peak.type != TString("s2") )
			{
				counter++;
				continue;
			}
			
			ExtractS2info(&peak, counter);
			counter++;
		}
		
	};
	
	
	
	void ExtractS1info(Peak* peak, int peakindex)
	{
		//Get the data from the summed waveform
		vector<double> s1s(peak->area_per_channel, peak->area_per_channel+259);//Contribution from each channel to the total  area
		double s1tot = peak->area;
		double s1top = (peak->area)*(peak->area_fraction_top);
		double s1bot = s1tot-s1top;
		
		int peakmaxpos = peak->index_of_maximum;
		double peakmax = peak->height;
		
		double centertime = peak->center_time;
		
		int ncoins = peak->n_contributing_channels;
		int ncoinstop = peak->n_contributing_channels_top;
		int ncoinsbot = ncoins - ncoinstop;
		
		vector<PmtContrib> pmtscontrib;//This must be sorted from the largest to the smallest
		for(int iPmt=0; iPmt<s1s.size(); iPmt++)
		{
			if(s1s.at(iPmt)>0)
			{
				pmtscontrib.push_back( PmtContrib( iPmt, s1s.at(iPmt)) );
			}
		}
		
		//SORT THE "pmtscontrib" by area and save the sorted pmt indexes
		sort(pmtscontrib.begin(), pmtscontrib.end(), greater<PmtContrib>());
		vector<int> S1sPmtOrder(pmtscontrib.size());
		for(int iPmt=0; iPmt<pmtscontrib.size(); iPmt++)
		{
			S1sPmtOrder.at(iPmt) = pmtscontrib.at(iPmt).index;
		}
		
		
		data.S1s_vec->at(peakindex) = s1s;
		data.S1sTot_vec->at(peakindex) = s1tot;
		data.S1sTop_vec->at(peakindex) = s1top;
		data.S1sBot_vec->at(peakindex) = s1bot;
		data.S1sPeakMaxPos_vec->at(peakindex) = peakmaxpos;
		data.S1sPeakMax_vec->at(peakindex) = peakmax;
		data.S1sCenterTime_vec->at(peakindex) = centertime;
		data.S1sNcoins_vec->at(peakindex) = ncoins;
		data.S1sNcoinsTop_vec->at(peakindex) = ncoinstop;
		data.S1sNcoinsBot_vec->at(peakindex) = ncoinsbot;
		data.S1sPmtOrder_vec->at(peakindex) = S1sPmtOrder;
	};
	
	
	void ExtractS2info(Peak* peak, int peakindex)
	{
		//Get the data from the summed waveform
		vector<double> s2s(peak->area_per_channel, peak->area_per_channel+259);//Contribution from each channel to the total  area
		double s2tot = peak->area;
		double s2top = (peak->area)*(peak->area_fraction_top);
		double s2bot = s2tot-s2top;
		
		int peakmaxpos = peak->index_of_maximum;
		double peakmax = peak->height;
		
		double centertime = peak->center_time;
		
		int ncoins = peak->n_contributing_channels;
		int ncoinstop = peak->n_contributing_channels_top;
		int ncoinsbot = ncoins - ncoinstop;
		
		vector<PmtContrib> pmtscontrib;//This must be sorted from the largest to the smallest
		for(int iPmt=0; iPmt<s2s.size(); iPmt++)
		{
			if(s2s.at(iPmt)>0)
			{
				pmtscontrib.push_back( PmtContrib( iPmt, s2s.at(iPmt) ) );
			}
		}
		
		//SORT THE "pmtscontrib" by area and save the sorted pmt indexes
		sort(pmtscontrib.begin(), pmtscontrib.end(), greater<PmtContrib>());
		vector<int> S2sPmtOrder( pmtscontrib.size() );
		for(int iPmt=0; iPmt<pmtscontrib.size(); iPmt++)
		{
			S2sPmtOrder.at(iPmt) = pmtscontrib.at(iPmt).index;
		}
		
		
		data.S2s_vec->at(peakindex) = s2s;
		data.S2sTot_vec->at(peakindex) = s2tot;
		data.S2sTop_vec->at(peakindex) = s2top;
		data.S2sBot_vec->at(peakindex) = s2bot;
		data.S2sPeakMaxPos_vec->at(peakindex) = peakmaxpos;
		data.S2sPeakMax_vec->at(peakindex) = peakmax;
		data.S2sCenterTime_vec->at(peakindex) = centertime;
		data.S2sNcoins_vec->at(peakindex) = ncoins;
		data.S2sNcoinsTop_vec->at(peakindex) = ncoinstop;
		data.S2sNcoinsBot_vec->at(peakindex) = ncoinsbot;
		data.S2sPmtOrder_vec->at(peakindex) = S2sPmtOrder;
	};
	
	
	
	
}//End of namespace minitrees


//using namespace minitrees;
void MakeMiniTree(string MiniTreeName)
{
	Event* paxEvent = NULL;


	if(!minitrees::datasets.size()) return;


	set<string>::iterator It;
	for(It=minitrees::datasets.begin(); It!=minitrees::datasets.end(); It++)
	{
		string infilename = minitrees::procdir + (*It) + string(".root");
	
		//Control that the file exists
		if( access(infilename.c_str(), R_OK|X_OK)!=0 ) continue;
	
		TFile *inFile = TFile::Open(infilename.c_str());
	
		if( (!inFile)||(!inFile->IsOpen()) ) continue;
	
		TTree *tree = (TTree*)inFile->Get("tree");
	
		if( !tree )
		{
			delete inFile;
			continue;
		}
	
		tree->SetBranchStatus("*", 0);
		tree->SetBranchStatus("events", 1);
	
		tree->SetBranchAddress("events", &paxEvent);
	
		if(!paxEvent)
		{
			delete inFile;
			tree = NULL;
			continue;
		}
	
		int nEvents = tree->GetEntries();
	
	
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
		
		
		//Loop over the events
		for(int iEv=0; iEv<nEvents; iEv++)
		{
			tree->GetEntry(iEv);
		
			minitrees::ProcessEvent(paxEvent);
		
			mt->Fill();
		
		}//End of the loop over the events
	
		outFile->WriteTObject(mt, 0, "overwrite");
	
		delete outFile;
		mt = NULL;
	
		delete inFile;
		tree = NULL;
	
	}//End of loop over datasets
}


//SERVICE FUNCTIONS
void AddDataset(string dsname)
{
	minitrees::datasets.insert(dsname);
};

void ResetDataSets()
{
	minitrees::datasets.clear();
};
