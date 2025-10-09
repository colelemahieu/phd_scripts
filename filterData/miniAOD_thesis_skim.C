//miniAOD_thesis_skim.C

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "TString.h"
#include "TTree.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TBranch.h"
#include "TFile.h"
#include <TTree.h>
#include <TChain.h>
#include "TCanvas.h"
#include "THStack.h"
#include "TLine.h"
#include "TMath.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
using namespace std;

void miniAOD_thesis_skim(string inputfile, string outfile)
{
  // access trees in files
  TChain* calochain = new TChain("ggHiNtuplizer/EventTree");
  //TChain* jetchain = new TChain("akCs2PFJetAnalyzer/t");
  TChain* jetchain = new TChain("akCs4PFJetAnalyzer/t");
  //TChain* jetchain = new TChain("akCs6PFJetAnalyzer/t"); 
  
  TChain* hiEvtchain = new TChain("hiEvtAnalyzer/HiTree");
  TChain* skimchain = new TChain("skimanalysis/HltTree");
  TChain* pfchain = new TChain("particleFlowAnalyser/pftree");
  TChain* zdcchain = new TChain("zdcanalyzer/zdcrechit");

  jetchain->Add(inputfile.c_str());
  hiEvtchain->Add(inputfile.c_str());
  skimchain->Add(inputfile.c_str());
  zdcchain->Add(inputfile.c_str());
  pfchain->Add(inputfile.c_str());

  int jetentries = jetchain->GetEntries();
  int hiEvtEntries = hiEvtchain->GetEntries();
  int skimEntries = skimchain->GetEntries();
  cout << "Entries in jet chain = " << jetentries << endl;
  cout << "Entries in hiEvt chain = " << hiEvtEntries << endl;
  cout << "Entries in skim chain = " << skimEntries << endl;

  // Tree Variables
  int nRef=0, nZDC=0, nTrk=0, nTower=0, n=0, nPFpart=0, vtxFilter=0, clusterFilter=0;
  unsigned int run=0, lumis=0;
  ULong64_t event=0;
  float vtx_z=0, ZDC_minus, ZDC_plus, jetEta1=0, jetEta2=0;
  float trkPt[10000], trkEta[10000];
  float jetPt[200]={0}, jetEta[200]={0}, jetPhi[200]={0}, jetM[200]={0};
  float e[18]={0};
  int zside[18]={0};

  vector<float> *jetPt_vec=0, *jetEta_vec=0, *jetPhi_vec=0;
  vector<float> *caloEta=0, *caloPhi=0, *caloEn=0, *calo_hadE=0, *calo_emE=0;
  vector<float> *trkPt_vec=0, *trkEta_vec=0;
  vector<float> *pfEnergy=0, *pfEta=0, *pfPhi=0, *pfPt=0, *pfEt=0, *pfE=0, *pfTrkEta=0, *pfTrkPt=0;
  vector<int> *pfId=0;

  // miniAOD read Variables
  float hiHE_pfle=0, hiHEMinus_pfle=0, hiHEPlus_pfle=0;
  float hiEE_pfle=0, hiEEMinus_pfle=0, hiEEPlus_pfle=0; 

  hiEvtchain->SetBranchAddress("vz",&vtx_z);
  hiEvtchain->SetBranchAddress("hiHE_pfle", &hiHE_pfle);
  hiEvtchain->SetBranchAddress("hiHEMinus_pfle", &hiHEMinus_pfle);
  hiEvtchain->SetBranchAddress("hiHEPlus_pfle", &hiHEPlus_pfle);
  hiEvtchain->SetBranchAddress("hiEE_pfle", &hiEE_pfle);
  hiEvtchain->SetBranchAddress("hiEEMinus_pfle", &hiEEMinus_pfle);
  hiEvtchain->SetBranchAddress("hiEEPlus_pfle", &hiEEPlus_pfle);
  skimchain->SetBranchAddress("pprimaryVertexFilter", &vtxFilter);
  skimchain->SetBranchAddress("pclusterCompatibilityFilter", &clusterFilter);

  jetchain->SetBranchAddress("nref", &nRef);
  jetchain->SetBranchStatus("*",0);
  jetchain->SetBranchStatus("jtpt",1);
  jetchain->SetBranchStatus("jteta",1);
  jetchain->SetBranchStatus("jtphi",1);
  jetchain->SetBranchStatus("jtm",1);
  jetchain->SetBranchAddress("jtpt", &jetPt);
  jetchain->SetBranchAddress("jteta", &jetEta);
  jetchain->SetBranchAddress("jtphi", &jetPhi);
  jetchain->SetBranchAddress("jtm", &jetM);
  
  pfchain->SetBranchStatus("*",0);
  pfchain->SetBranchStatus("nPF",1);
  pfchain->SetBranchStatus("pfPt",1);
  pfchain->SetBranchStatus("pfEt",1);
  pfchain->SetBranchStatus("pfE",1);
  pfchain->SetBranchStatus("pfEta",1);
  pfchain->SetBranchStatus("pfPhi",1);
  pfchain->SetBranchStatus("pfId",1);
  pfchain->SetBranchStatus("pfTrkEta",1);
  pfchain->SetBranchStatus("pfTrkPt",1);
  pfchain->SetBranchAddress("nPF", &nPFpart);
  pfchain->SetBranchAddress("pfPt", &pfPt);
  pfchain->SetBranchAddress("pfEt", &pfEt);
  pfchain->SetBranchAddress("pfE", &pfE);
  pfchain->SetBranchAddress("pfEta", &pfEta);
  pfchain->SetBranchAddress("pfPhi", &pfPhi);
  pfchain->SetBranchAddress("pfId", &pfId);
  pfchain->SetBranchAddress("pfTrkEta", &pfTrkEta);
  pfchain->SetBranchAddress("pfTrkPt", &pfTrkPt);

  zdcchain->SetBranchStatus("*",0);
  zdcchain->SetBranchStatus("n",1);
  zdcchain->SetBranchStatus("e",1);
  zdcchain->SetBranchStatus("zside",1);
  zdcchain->SetBranchAddress("n", &nZDC);
  zdcchain->SetBranchAddress("e", &e);
  zdcchain->SetBranchAddress("zside", &zside);

  // Clone Trees
  TFile newfile(outfile.c_str(),"recreate");
  TTree *pfClone = pfchain->CloneTree(0);
  TTree *jetClone = jetchain->CloneTree(0);
  TTree *zdcClone = zdcchain->CloneTree(0);

  // New ROOT File
  TTree *hiEvent = new TTree("hiEvent","hiEvent");
  hiEvent->Branch("Vertex_Z", &vtx_z);
  hiEvent->Branch("hiHE_pfle", &hiHE_pfle); 
  hiEvent->Branch("hiHEMinus_pfle", &hiHEMinus_pfle);
  hiEvent->Branch("hiHEPlus_pfle", &hiHEPlus_pfle);
  hiEvent->Branch("hiEE_pfle", &hiEE_pfle);
  hiEvent->Branch("hiEEMinus_pfle", &hiEEMinus_pfle);
  hiEvent->Branch("hiEEPlus_pfle", &hiEEPlus_pfle);
  hiEvent->Branch("primaryVertexFilter", &vtxFilter);
  hiEvent->Branch("clusterCompatibilityFilter", &clusterFilter);

  int goodEvt=0;
  int pfGoodTrk, pfHadronYes;

  
  // Event Loop
  for (int i=0; i<skimEntries; i++)
    {
      jetchain->GetEntry(i);
      hiEvtchain->GetEntry(i);
      skimchain->GetEntry(i);
      zdcchain->GetEntry(i);
      pfchain->GetEntry(i);
      
      ZDC_minus=0;
      ZDC_plus=0; 
      pfGoodTrk=0, pfHadronYes=0;

      // Cuts
      if (vtxFilter==0) continue;
      if (clusterFilter==0) continue;

      for (int pf=0; pf<nPFpart; pf++)
      {
        if (pfId->at(pf)==1) pfHadronYes = pfHadronYes + 1;
      }

      if (pfHadronYes<2) continue;
      if (fabs(vtx_z) > 20) continue;

     
      // ZDC Info
      for (int j=0; j<nZDC; j++)
        {
          if (zside[j]>0)
            ZDC_plus += e[j];
          if (zside[j]<0)
            ZDC_minus += e[j];
        }
      
      //if (ZDC_minus>1200 || ZDC_plus<1200) continue;  // 0nXn
      //if (ZDC_minus<1200 || ZDC_plus>1200) continue; // Xn0n
      if ((ZDC_minus>1200 && ZDC_plus>1200) || (ZDC_minus<1200 && ZDC_plus<1200)) continue;
      //if (ZDC_minus>1200 || ZDC_plus>1200) continue;   // 0n0n

      if (nRef==0) continue;
      if (nRef==1) continue;
      if ((fabs(jetEta[0]) > 3.0) || (fabs(jetEta[1]) > 3.0)) continue;
      
      // Fill Branches
      pfClone->Fill();     
      jetClone->Fill();
      zdcClone->Fill();
      hiEvent->Fill();
   
      goodEvt=  goodEvt+ 1;
    }
 
  cout << "goodEvt=" << goodEvt << endl;
  newfile.Write();
  newfile.Close();
}
