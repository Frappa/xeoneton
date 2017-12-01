#ifndef XE1T_ANALYSIS_HH
#define XE1T_ANALYSIS_HH

#include "TTree.h"
#include "TChain.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraph.h"
#include "TProfile2D.h"
#include "TString.h"

#include <iostream>
#include <string>


using namespace std;

class xe1t_analysis{

 public:
  xe1t_analysis(){};
  virtual ~xe1t_analysis(){};

  virtual void Loop();
  virtual void SetOutputFileName(TString f){OutputFileName = f;}
  virtual void initialize();
  virtual void execute();
  virtual void finalize();


 private:
  TFile *m_file;
  TString OutputFileName;

  // Event No
  Int_t nentries;
  Int_t jentry;
  Int_t ientry;
  Int_t nbytes, nb;
  

  TTree *out_tree;
  
  Int_t EventNumber;
  Long64_t run_start_time;
  Long64_t time_since_run_start;
  
};


#endif
