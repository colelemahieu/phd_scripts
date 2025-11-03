// getThesis_histograms.C

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
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
#include "TProfile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "../commonUtility.h"
using namespace std;

int return_z(TH2 *h2_dummy, float qt, float phi)
{
 int x_bin = h2_dummy->GetXaxis()->FindBin(qt);
 int y_bin = h2_dummy->GetYaxis()->FindBin(phi);
 int z = y_bin + 12*(x_bin-1);

 return z;
}


// calculate rapidity gap in the negative eta region
float backRapGap(int nTrk, vector<float> *trkPt, vector<float> *trkEta, float ptThreshold)
{
  // FOR USE INSIDE A LOOP OVER EVENTS

  float etaEdg = 2.5;
  float etaBinNumber = 25;
  float etaBin = 2*etaEdg/etaBinNumber;
  float rapGap = 0;
  bool hit = false;

  // loop over eta bins
  for(int j=1; j<etaBinNumber; j++)
    {
      // loop over tracks
      for (int k=0; k<nTrk; k++)
        {
          // For hit=true, track must be above threshold, and within the current bin, starting from the left
          if (trkPt->at(k) > ptThreshold && (trkEta->at(k) > -1.0*etaEdg && trkEta->at(k) < (-1.0*etaEdg+etaBin*j))) hit=true;
        }

      if( hit == false)
        {
          rapGap = rapGap+etaBin;
        }
    }

  return rapGap;
}



// calculate rapidity gap in the positive eta region
float forwardRapGap(int nTrk, vector<float> *trkPt, vector<float> *trkEta, float ptThreshold)
{
  // FOR USE INSIDE A LOOP OVER EVENTS

  float etaEdg = 2.5;
  float etaBinNumber = 25;
  float etaBin = 2*etaEdg/etaBinNumber;
  float rapGap = 0;
  bool hit = false;

  // loop over eta bins
  for(int j=1; j<etaBinNumber; j++)
     {
      //cout << " " << endl;
      // loop over tracks
      for (int k=0; k<nTrk; k++)
        {
          // For hit=true, track must be above threshold, and within the current bin, starting from the right
	  //cout << trkPt[k] << endl;
          if (trkPt->at(k) > ptThreshold && (trkEta->at(k) < etaEdg && trkEta->at(k) > (etaEdg-etaBin*j))) hit=true;
        }

      if( hit == false)
        {
          rapGap = rapGap+etaBin;
        }
    }

  return rapGap;
}


float backRapGap_ext(int nTrk, vector<float> *trkPt, vector<float> *trkEta, float ptThreshold)
{
  // FOR USE INSIDE A LOOP OVER EVENTS

  float etaEdg = 3.0;
  float etaBinNumber = 30;
  float etaBin = 2*etaEdg/etaBinNumber;
  float rapGap = 0;
  bool hit = false;

  // loop over eta bins
  for(int j=1; j<etaBinNumber; j++)
    {
      // loop over tracks
      for (int k=0; k<nTrk; k++)
        {
          // For hit=true, track must be above threshold, and within the current bin, starting from the left
          if (trkPt->at(k) > ptThreshold && (trkEta->at(k) > -1.0*etaEdg && trkEta->at(k) < (-1.0*etaEdg+etaBin*j))) hit=true;
        }

      if( hit == false)
        {
          rapGap = rapGap+etaBin;
        }
    }

  return rapGap;
}



// calculate rapidity gap in the positive eta region
float forwardRapGap_ext(int nTrk, vector<float> *trkPt, vector<float> *trkEta, float ptThreshold)
{
  // FOR USE INSIDE A LOOP OVER EVENTS

  float etaEdg = 3;
  float etaBinNumber = 30;
  float etaBin = 2*etaEdg/etaBinNumber;
  float rapGap = 0;
  bool hit = false;

  // loop over eta bins
  for(int j=1; j<etaBinNumber; j++)
    {
      //cout << " " << endl;
      // loop over tracks
      for (int k=0; k<nTrk; k++)
        {
          // For hit=true, track must be above threshold, and within the current bin, starting from the right
	  //cout << trkPt[k] << endl;
          if (trkPt->at(k) > ptThreshold && (trkEta->at(k) < etaEdg && trkEta->at(k) > (etaEdg-etaBin*j))) hit=true;
        }

      if( hit == false)
        {
          rapGap = rapGap+etaBin;
        }
    }

  return rapGap;
}




void getThesis_histograms()
{
  bool thirdJetCut=true;
  bool dPhiCut=false;
  
  // new data
  //TFile *file = new TFile("../miniAOD_r2_newCuts.root");
  TFile *file = new TFile("../miniAOD_r4_newCuts.root");
  //TFile *file = new TFile("../miniAOD_r6_newCuts.root");

  // read data
  TTree *jetTree = (TTree*)file->Get("t");
  TTree *zdcTree = (TTree*)file->Get("zdcrechit");
  TTree *evtTree = (TTree*)file->Get("hiEvent");
  TTree *pfTree = (TTree*)file->Get("pftree");


  // noise thresholds
  float heM_noiseTH=1.1, heP_noiseTH=1.06, eeM_noiseTH=19, eeP_noiseTH=16.2, hfM_noiseTH=11.75, hfP_noiseTH=11.75, trackGap=1.8;
  // vary HE by 5%
  //float heM_noiseTH=0.78, heP_noiseTH=0.78, eeM_noiseTH=19, eeP_noiseTH=16.2, hfM_noiseTH=11.75, hfP_noiseTH=11.75, trackGap=1.8;
  // vary EE by 5%
  //float heM_noiseTH=1.1, heP_noiseTH=1.06, eeM_noiseTH=0.2, eeP_noiseTH=0.2, hfM_noiseTH=11.75, hfP_noiseTH=11.75, trackGap=1.8;
  // vary HF by 5%
  //float heM_noiseTH=1.1, heP_noiseTH=1.06, eeM_noiseTH=19, eeP_noiseTH=16.2, hfM_noiseTH=6.75, hfP_noiseTH=6.75, trackGap=1.8;
  // vary rapGap up
  //float heM_noiseTH=1.1, heP_noiseTH=1.06, eeM_noiseTH=19, eeP_noiseTH=16.2, hfM_noiseTH=11.75, hfP_noiseTH=11.75, trackGap=2.0;
  // vary rapGap down
  //float heM_noiseTH=1.1, heP_noiseTH=1.06, eeM_noiseTH=19, eeP_noiseTH=16.2, hfM_noiseTH=11.75, hfP_noiseTH=11.75, trackGap=1.6;
  

  // read variables
  int nTrk=0, nZDC=0, nRef=0;
  float vtx_z=0;
  unsigned int run=0, lumis=0;
  ULong64_t event=0;
  float jetEnt = jetTree->GetEntries();
  vector<float> *trkPt=0, *trkEta=0;
  float jetPt[200]={0}, jetEta[200]={0}, jetPhi[200]={0}, jetM[200]={0};
  float genPt[200]={0}, genEta[200]={0}, genPhi[200]={0}, genM[200]={0};
  vector<float> *caloEta=0, *caloPhi=0, *caloEn=0, *calo_hadE=0, *calo_emE=0;
  float e[18]={0};
  int zside[18]={0};
  float hiHE_pfle=0, hiEE_pfle=0;

  int nPF=0;
  vector<int> *pfId=0;
  vector<float> *pfEta=0, *pfE=0, *pfEt=0, *pfPhi=0, *pfTrkEta=0, *pfTrkPt=0;

  // data
  jetTree->SetBranchAddress("nref", &nRef);
  jetTree->SetBranchAddress("jteta", &jetEta);
  jetTree->SetBranchAddress("jtpt", &jetPt);
  jetTree->SetBranchAddress("jtphi", &jetPhi);
  jetTree->SetBranchAddress("jtm", &jetM);

  zdcTree->SetBranchAddress("n", &nZDC);
  zdcTree->SetBranchAddress("e", &e);
  zdcTree->SetBranchAddress("zside", &zside);
  evtTree->SetBranchAddress("Vertex_Z", &vtx_z);
  evtTree->SetBranchAddress("hiEE_pfle", &hiEE_pfle);
  evtTree->SetBranchAddress("hiHE_pfle", &hiHE_pfle);
  pfTree->SetBranchAddress("nPF", &nPF);
  pfTree->SetBranchAddress("pfEt", &pfEt);
  pfTree->SetBranchAddress("pfE", &pfE);
  pfTree->SetBranchAddress("pfEta", &pfEta);
  pfTree->SetBranchAddress("pfPhi", &pfPhi);
  pfTree->SetBranchAddress("pfTrkEta", &pfTrkEta);
  pfTree->SetBranchAddress("pfTrkPt", &pfTrkPt);
  pfTree->SetBranchAddress("pfId", &pfId);

  // Analysis
  // Vectors for Calculating v2
  float rapPrime=0;
  float pi = TMath::Pi();
  float px_1=0, px_2=0, py_1=0, py_2=0, avgRap=0;
  TVector2 Q_T, P_T;
  float QT_norm, PT_norm, QT_phi, PT_phi, angleDiff;
  TVector2 QT_unit, PT_unit;
  float cos12, sin12, angle12, cos_2phi=0, cos_4phi=0, cos_6phi=0, cos_8phi=0, cos_10phi=0;
  float avg_v1=0, avg_v2=0, avg_v4=0, avg_v6=0, avg_v8=0, avg_v10=0;
  float v1_total=0, v2_total=0, v4_total=0, v6_total=0, v8_total=0, v10_total=0;
  int goodEvt=0;
  float heBad, eeBad, ZDC_minus, ZDC_plus;
  float eeMax_0nXn, eeMax_Xn0n;

  // New Root File
  TFile newfile("allPlots_test.root","recreate");
  TH1F *h_zvtx = new TH1F("h_zvtx","h_zvtx",100,-16,16);
  TH1F *h_zdcMinus = new TH1F("h_zdcMinus", "h_zdcMinus", 500, 0, 155000);
  TH1F *h_zdcPlus = new TH1F("h_zdcPlus", "h_zdcPlus", 500, 0, 155000);
  TH1F *h_ee_0nXn = new TH1F("h_ee_0nXn","h_ee_0nXn",100,0,20);
  TH1F *h_ee_Xn0n = new TH1F("h_ee_Xn0n","h_ee_Xn0n",100,0,20);
  TH1F *h_nRef = new TH1F("h_nRef","nRef",10,0,10);
  TH1F *h_jtPt_L = new TH1F("h_jtPt_L","h_jtPt_L",40,0,80);
  TH1F *h_jtPt_S = new TH1F("h_jtPt_S","h_jtPt_S",40,0,80);
  TH1F *h_dijetMass = new TH1F("h_dijetMass","h_dijetMass",50,0,150);
  TH2F *h2_rap1_rap2 = new TH2F("h2_rap1_rap2","h2_rap1_rap2",20,-2.4,2.4,20,-2.4,2.4);
  TH1F *h_dijetRap = new TH1F("h_dijetRap","h_dijetRap",30,-2.5,2.5);
  TH1F *h_dijetRapPrime = new TH1F("h_dijetRapPrime","h_dijetRapPrime",10,-2.5,2.5);
  TH1F *h_PT = new TH1F("h_PT","h_PT",40,0,80);
  TH1F *h_QT = new TH1F("h_QT","h_QT",40,0,80);
  TH1F *h_PT_10 = new TH1F("h_PT_10","h_PT_10",10,0,80);
  TH1F *h_QT_10 = new TH1F("h_QT_10","h_QT_10",10,0,80);
  TH1D *h_QT_meas = new TH1D("h_QT_meas","h_QT_meas",40,0,40);
  TH1D *h_QT10_meas = new TH1D("h_QT10_meas","h_QT10_meas",10,0,40);
  TH1D *h_phiQTPT_12 = new TH1D("h_phiQTPT_12","h_phiQTPT_12",12,-TMath::Pi(),TMath::Pi());
  TH1D *h_phiQTPT_12_abs = new TH1D("h_phiQTPT_12_abs","h_phiQTPT_12_abs",12,0,TMath::Pi());
  float minHist=-TMath::Pi(), maxHist=TMath::Pi();
  TH1D *h_phiQTPT_bin1 = new TH1D("h_phiQTPT_bin1","h_phiQTPT_bin1",12,-TMath::Pi(),TMath::Pi());
  TH1D *h_phiQTPT_bin2 = new TH1D("h_phiQTPT_bin2","h_phiQTPT_bin2",12,-TMath::Pi(),TMath::Pi());
  TH1D *h_phiQTPT_bin3 = new TH1D("h_phiQTPT_bin3","h_phiQTPT_bin3",12,-TMath::Pi(),TMath::Pi());
  TH1D *h_phiQTPT_bin4 = new TH1D("h_phiQTPT_bin4","h_phiQTPT_bin4",12,-TMath::Pi(),TMath::Pi());
  TH1D *h_phiQTPT_bin5 = new TH1D("h_phiQTPT_bin5","h_phiQTPT_bin5",12,-TMath::Pi(),TMath::Pi());
  TH1D *h_phiQTPT_6_bin1 = new TH1D("h_phiQTPT_6_bin1","h_phiQTPT_6_bin1",6,-TMath::Pi(),TMath::Pi());
  TH1D *h_phiQTPT_6_bin2 = new TH1D("h_phiQTPT_6_bin2","h_phiQTPT_6_bin2",6,-TMath::Pi(),TMath::Pi());
  TH1D *h_phiQTPT_6_bin3 = new TH1D("h_phiQTPT_6_bin3","h_phiQTPT_6_bin3",6,-TMath::Pi(),TMath::Pi());
  TH1D *h_phiQTPT_6_bin4 = new TH1D("h_phiQTPT_6_bin4","h_phiQTPT_6_bin4",6,-TMath::Pi(),TMath::Pi());
  TH1D *h_phiQTPT_6_bin5 = new TH1D("h_phiQTPT_6_bin5","h_phiQTPT_6_bin5",6,-TMath::Pi(),TMath::Pi());
  TH1D *h_cos2phi_meas = new TH1D("h_cos2phi_meas","h_cos2phi_meas",40,-1,1);
  TH1D *h_cos2phi_meas_12 = new TH1D("h_cos2phi_meas_12","h_cos2phi_meas_12",12,-1,1);
  TH2D *h2_cos12_QT5_meas = new TH2D("h2_cos12_QT5_meas","h2_cos12_QT5_meas",5,0,40,12,-1,1);
  TH2F *h2_QTPT = new TH2F("h2_QTPT", "h2_QTPT",100,0,40,100,20,100);
  TH1F *h_deltaPhi_pre = new TH1F("h_deltaPhi_pre","h_deltaPhi_pre",100,-2*(TMath::Pi()),2*(TMath::Pi()));
  TH2D *h2_QT_cos = new TH2D("h2_QT_cos","h2_QT_cos",40,0.0,40.0,40,-1.0,1.0);
  TH2D *h2_QT8_cos = new TH2D("h2_QT8_cos","h2_QT8_cos",8,0.0,40.0,40,-1.0,1.0);
  TH2D *h2_QT8_cos20 = new TH2D("h2_QT8_cos20","h2_QT8_cos20",8,0.0,40.0,20,-1.0,1.0);
  TH2D *h2_QT8_cos15 = new TH2D("h2_QT8_cos15","h2_QT8_cos15",8,0.0,40.0,15,-1.0,1.0);
  TH1F *h_deltaPhi = new TH1F("h_deltaPhi","h_deltaPhi",100,0,2*(TMath::Pi()));
  TH1F *h_angle12_pre = new TH1F("h_angle12_pre","h_angle12_pre",100,0,2*(TMath::Pi()));
  TH1F *h_angle12 = new TH1F("h_angle12","h_angle12",80,-TMath::Pi(),TMath::Pi());
  TH2F *h2_rapgaps_ext = new TH2F("h2_rapgaps_ext", "h2_rapgaps_ext",25,0,5,25,0,5);
  TH2F *h2_rapgaps_0nXn = new TH2F("h2_rapgaps_0nXn", "h2_rapgaps_0nXn",25,0,5,25,0,5);
  TH2F *h2_rapgaps_Xn0n = new TH2F("h2_rapgaps_Xn0n", "h2_rapgaps_Xn0n",25,0,5,25,0,5);
  TH1F *h_rapgap_forward = new TH1F("h_rapgaps_forward", "h_rapgaps_forward",25,0,5);
  TH1F *h_rapgap_backward = new TH1F("h_rapgaps_backward", "h_rapgaps_backward",25,0,5);
  TH1F *h_dijetRap_0nXn = new TH1F("h_dijetRap_0nXn", "h_dijetRap_0nXn",25,-3,3);
  TH1F *h_dijetRap_Xn0n = new TH1F("h_dijetRap_Xn0n", "h_dijetRap_Xn0n",25,-3,3);
  TProfile *h_profPT_5  = new TProfile("h_profPT_5","Profile PT 5",5,0,40,0,150);
  TProfile *h_profV1  = new TProfile("h_profV1","Profile v1",12,0,40,-1,1);
  TProfile *h_profV2  = new TProfile("h_profV2","Profile v2",12,0,40,-1,1);
  TProfile *h_profV2_5  = new TProfile("h_profV2_5","Profile v2 5",5,0,40,-1,1);
  TProfile *h_profV2_8  = new TProfile("h_profV2_8","Profile v2 8",8,0,40,-1,1);
  TProfile *h_profV2_10  = new TProfile("h_profV2_10","Profile v2 10",10,0,40,-1,1);
  TProfile *h_profV4  = new TProfile("h_profV4","Profile v4",12,0,40,-1,1);
  TProfile *h_profV6  = new TProfile("h_profV6","Profile v6",12,0,40,-1,1);
  TProfile *h_profV8  = new TProfile("h_profV8","Profile v8",12,0,40,-1,1);
  TProfile *h_profV10  = new TProfile("h_profV10","Profile v10",12,0,40,-1,1);
  TProfile *h_profV2_rap1  = new TProfile("h_profV2_rap1","Profile v2 rap1",12,0,40,-1,1);
  TProfile *h_profV2_rap2  = new TProfile("h_profV2_rap2","Profile v2 rap2",12,0,40,-1,1);
  TProfile *h_profV2_rap3  = new TProfile("h_profV2_rap3","Profile v2 rap3",12,0,40,-1,1);
  TProfile *h_profV2_rap4  = new TProfile("h_profV2_rap4","Profile v2 rap4",12,0,40,-1,1);
  TProfile *h_profV2_rap5  = new TProfile("h_profV2_rap5","Profile v2 rap5",12,0,40,-1,1);
  TProfile *h_profV2_rap6  = new TProfile("h_profV2_rap6","Profile v2 rap6",12,0,40,-1,1);
  TProfile *h_profV2_rap7  = new TProfile("h_profV2_rap7","Profile v2 rap7",12,0,40,-1,1);
  TProfile *h_profV2_rap8  = new TProfile("h_profV2_rap8","Profile v2 rap8",12,0,40,-1,1);
  TProfile *h_profV2_rap1_4  = new TProfile("h_profV2_rap1_4","Profile v2 rap1_4",4,0,40,-1,1);
  TProfile *h_profV2_rap2_4  = new TProfile("h_profV2_rap2_4","Profile v2 rap2_4",4,0,40,-1,1);
  TProfile *h_profV2_rap3_4  = new TProfile("h_profV2_rap3_4","Profile v2 rap3_4",4,0,40,-1,1);
  TProfile *h_profV2_rap4_4  = new TProfile("h_profV2_rap4_4","Profile v2 rap4_4",4,0,40,-1,1);
  TProfile *h_profV2_rap5_4  = new TProfile("h_profV2_rap5_4","Profile v2 rap5_4",4,0,40,-1,1);
  TProfile *h_profV2_rap6_4  = new TProfile("h_profV2_rap6_4","Profile v2 rap6_4",4,0,40,-1,1);
  TProfile *h_profV2_rap7_4  = new TProfile("h_profV2_rap7_4","Profile v2 rap7_4",4,0,40,-1,1);
  TProfile *h_profV2_rap8_4  = new TProfile("h_profV2_rap8_4","Profile v2 rap8_4",4,0,40,-1,1);
  TProfile *h_profV2_vtx1  = new TProfile("h_profV2_vtx1","Profile v2 vtx1",8,0,40,-1,1);
  TProfile *h_profV2_vtx2  = new TProfile("h_profV2_vtx2","Profile v2 vtx2",8,0,40,-1,1);
  TProfile *h_profV2_vtx3  = new TProfile("h_profV2_vtx3","Profile v2 vtx3",8,0,40,-1,1);
  TProfile *h_profV2_vtx4  = new TProfile("h_profV2_vtx4","Profile v2 vtx4",8,0,40,-1,1);
  TProfile *h_profV2_vtx5  = new TProfile("h_profV2_vtx5","Profile v2 vtx5",8,0,40,-1,1);
  TProfile *h_profV2_vtx6  = new TProfile("h_profV2_vtx6","Profile v2 vtx6",8,0,40,-1,1);
  // QT Bins of 3
  TProfile *h_profV2_rapPrime1_3  = new TProfile("h_profV2_rapPrime1_3","h_profV2_rapPrime1_3",3,0,40,-1,1);
  TProfile *h_profV2_rapPrime2_3  = new TProfile("h_profV2_rapPrime2_3","h_profV2_rapPrime2_3",3,0,40,-1,1);
  TProfile *h_profV2_rapPrime3_3  = new TProfile("h_profV2_rapPrime3_3","h_profV2_rapPrime3_3",3,0,40,-1,1);
  TProfile *h_profV2_rapPrime4_3  = new TProfile("h_profV2_rapPrime4_3","h_profV2_rapPrime4_3",3,0,40,-1,1);
  TProfile *h_profV2_rapPrime5_3  = new TProfile("h_profV2_rapPrime5_3","h_profV2_rapPrime5_3",3,0,40,-1,1);
  // QT Bins of 2
  double prime_bins[] = {0, 12, 40};
  TProfile *h_profV2_rapPrime1_2  = new TProfile("h_profV2_rapPrime1_2","h_profV2_rapPrime1_2",2,prime_bins,-1,1);
  TProfile *h_profV2_rapPrime2_2  = new TProfile("h_profV2_rapPrime2_2","h_profV2_rapPrime2_2",2,prime_bins,-1,1);
  TProfile *h_profV2_rapPrime3_2  = new TProfile("h_profV2_rapPrime3_2","h_profV2_rapPrime3_2",2,prime_bins,-1,1);
  TProfile *h_profV2_rapPrime4_2  = new TProfile("h_profV2_rapPrime4_2","h_profV2_rapPrime4_2",2,prime_bins,-1,1);
  TProfile *h_profV2_rapPrime5_2  = new TProfile("h_profV2_rapPrime5_2","h_profV2_rapPrime5_2",2,prime_bins,-1,1);

  // 2D Phi vs Rapidity
  const int nBins=5;
  double rapidityBins[nBins+1]= {-2.5,0.0,0.5,1.0,1.5,2.5};
  TH2D *h2_phi12_rap_qt0 = new TH2D("h2_phi12_rap_qt0","h2_phi12_rap_qt0",nBins,rapidityBins, 12,-TMath::Pi(),TMath::Pi());
  TH2D *h2_phi12_rap_qt1 = new TH2D("h2_phi12_rap_qt1","h2_phi12_rap_qt1",nBins,rapidityBins, 12,-TMath::Pi(),TMath::Pi());

  // Phi vs Rapidity Distributions for Unfolding
  TH1F *h_phi_rapidity[2][5];
  for (int qtBin=0; qtBin<2; qtBin++)
    {
      for (int rapBin=0; rapBin<5; rapBin++)
	{
	  h_phi_rapidity[qtBin][rapBin] = new TH1F(Form("h_phi_rapidity_%i_%i",qtBin,rapBin),Form("h_phi_rapidity_%i_%i",qtBin,rapBin),12,-TMath::Pi(),TMath::Pi());
	} // end rap loop
    } // end qtBin loop

  // HE and EE histograms
  TH1F *h_heP_2p8 = new TH1F("h_heP_2p8","h_heP_2p8",100,0,4);
  TH1F *h_heM_2p8 = new TH1F("h_heM_2p8","h_heM_2p8",100,0,4);
  TH1F *h_eeP_2p8 = new TH1F("h_eeP_2p8","h_eeP_2p8",100,0,40);
  TH1F *h_eeM_2p8 = new TH1F("h_eeM_2p8","h_eeM_2p8",100,0,40);
  TH1F *h_chH_eta_0nXn = new TH1F("h_chH_eta_0nXn","h_chH_eta_0nXn",50,-5,5);
  TH1F *h_chH_eta_Xn0n = new TH1F("h_chH_eta_Xn0n","h_chH_eta_Xn0n",50,-5,5);

  // PF histograms
  TH1F *h_radius = new TH1F("h_radius","h_radius",100,0,5);
  TH1F *h_pfEt = new TH1F("h_pfEt","h_pfEt",100,0,30);
  TH1F *h_pfEta = new TH1F("h_pfEta","h_pfEta",100,-3.2,3.2);
  TH1F *h_pfEtaW = new TH1F("h_pfEtaW","h_pfEtaW",100,-3.2,3.2);

  TH1F *h_pfEt_2p8_15 = new TH1F("h_pfEt_2p8_15","h_pfEt_2p8_15",100,0,5);
  TH1F *h_pfEtsum_2p8_15 = new TH1F("h_pfEtsum_2p8_15","h_pfEtsum_2p8_15",100,0,5);
  TH1F *h_pfEsum_2p8_15 = new TH1F("h_pfEsum_2p8_15","h_pfEsum_2p8_15",100,0,30);
  TH1F *h_pfEt_2p8_4 = new TH1F("h_pfEt_2p8_4","h_pfEt_2p8_4",100,0,5);
  TH1F *h_pfE_2p8_4 = new TH1F("h_pfE_2p8_4","h_pfE_2p8_4",100,0,5);
  TH1F *h_pfId_2p8 = new TH1F("h_pfId_2p8", "h_pfId_2p8",8,0,8);

  // New Unfolding
  TH2D *h2_phi6_QT5 = new TH2D("h2_phi6_QT5","h2_phiQTPT6_QT5",5,0,40,6,0,2*(TMath::Pi()));
  TH2D *h2_phiQTPT_QT = new TH2D("h2_phiQTPT_QT","h2_phiQTPT_QT ",5,0,40,12,0,2*(TMath::Pi()));
  TH2D *h2_phiQTPT_QT_mpi = new TH2D("h2_phiQTPT_QT_mpi","h2_phiQTPT_QT_mpi ",5,0,40,12,-TMath::Pi(),TMath::Pi());
  TH2D *h2_cosQTPT_QT_mpi = new TH2D("h2_cosQTPT_QT_mpi","h2_cosQTPT_QT_mpi ",30,0,40,30,-1,1);
  TH2D *h2_cosQTPT_rap0 = new TH2D("h2_cosQTPT_rap0","h2_cosQTPT_rap0",30,-2.5,2.5,30,-1,1);
  TH2D *h2_cosQTPT_rap1 = new TH2D("h2_cosQTPT_rap1","h2_cosQTPT_rap1",30,-2.5,2.5,30,-1,1);
  TH2D *h2_phi9_QT5 = new TH2D("h2_phi9_QT5","h2_phiQTPT9_QT5",5,0,40,9,0,2*(TMath::Pi()));

  // angles
  TH1F *h_trigPhi = new TH1F("h_trigPhi", "h_trigPhi", 100, -3.4, 3.4);
  // for twiki
  TH1F *h_nPF = new TH1F("h_nPF","h_nPF",50,0,200);
  TH1F *h_nPF_H = new TH1F("h_nPF_H","h_nPF_H",60,0,60);
  TH1F *h_nPF_ch = new TH1F("h_nPF_ch","h_nPF_ch",40,0,40);
  TH2F *h_ntrk_vtx = new TH2F("h_ntrk_vtx","h_ntrk_vtx",20,50,200,20,-20,20);
  TH2F *h_nPF_vtx = new TH2F("h_nPF_vtx","h_nPF_vtx",20,50,200,20,-20,20);
  TH2F *h_chH_vtx = new TH2F("h_chH_vtx","h_chH_vtx",30,0,30,20,-20,20);
  TH1F *h_zvtx_3chH = new TH1F("h_zvtx_3chH","h_zvtx_3chH",100,-16,16);
  TH1F *h_zvtx_4chH_up = new TH1F("h_zvtx_4chH_up","h_zvtx_4chH_up",100,-16,16);
  TH1F *h_PF_pt = new TH1F("h_PF_pt", "h_PF_pt", 100, 0, 20);
  TH1F *h_zvtx_4chH_up_300 = new TH1F("h_zvtx_4chH_up_300","h_zvtx_4chH_up_300",100,-16,16);

  // forward folding
  TH2F *h2_dummy = new TH2F("h2_dummy","h2_dummy",5,0,40,12,-TMath::Pi(), TMath::Pi());
  TH1D *h_measZ = new TH1D("h_measZ","h_measZ",60,0,60);
  TH2F *h2_dummy_rap0 = new TH2F("h2_dummy_rap0","h2_dummy_rap0",nBins,rapidityBins,12,-TMath::Pi(), TMath::Pi());
  TH1D *h_measZ_rap0 = new TH1D("h_measZ_rap0","h_measZ_rap0",60,0,60);
  TH2F *h2_dummy_rap1 = new TH2F("h2_dummy_rap1","h2_dummy_rap1",nBins,rapidityBins,12,-TMath::Pi(), TMath::Pi());
  TH1D *h_measZ_rap1 = new TH1D("h_measZ_rap1","h_measZ_rap1",60,0,60);

  // pf leading energies
  TH1F *h_pf_hfP_max_0nXn = new TH1F("h_pf_hfP_max_0nXn","h_pf_hfP_max_0nXn",100,0,50);
  TH1F *h_pf_hfP_max_Xn0n = new TH1F("h_pf_hfP_max_Xn0n","h_pf_hfP_max_Xn0n",100,0,50);
  TH1F *h_pf_hfM_max_0nXn = new TH1F("h_pf_hfM_max_0nXn","h_pf_hfM_max_0nXn",100,0,50);
  TH1F *h_pf_hfM_max_Xn0n = new TH1F("h_pf_hfM_max_Xn0n","h_pf_hfM_max_Xn0n",100,0,50);
  


  // pf leading energies
  float pfle_had, pfle_em;
  float pftotal_had, pftotal_hadE, pftotal_hadE_plus;
  float pf_hfM, pf_hfP, pf_hfM_max, pf_hfP_max;
  float pf_heM, pf_heP, pf_heM_max, pf_heP_max;
  float pf_eeM, pf_eeP, pf_eeM_max, pf_eeP_max;
  float pfId_H_sum=0, pfId_chH_sum=0, pfId_chH_sum_less300=0;
  int chH_0nXn, chH_Xn0n;
  // pf eta max
  float forward_rapgap=0, backward_rapgap=0, pf_etaMax, pf_etaMin;
  int badJets=0;
  

  // event loop
  for (int i=0; i<jetEnt; i++)
    {
      jetTree->GetEntry(i);
      zdcTree->GetEntry(i);
      evtTree->GetEntry(i);
      pfTree->GetEntry(i);
      heBad=0;
      eeBad=0;
      ZDC_minus=0;
      ZDC_plus=0;
      eeMax_0nXn=0;
      eeMax_Xn0n=0;
      pfId_H_sum=0, pfId_chH_sum=0, pfId_chH_sum_less300=0;
      pfle_had=0, pfle_em=0, pftotal_had=0, pftotal_hadE=0, pftotal_hadE_plus=0;
      pf_hfM=0, pf_hfP=0, pf_hfM_max=0, pf_hfP_max=0;
      pf_heM=0, pf_heP=0, pf_heM_max=0, pf_heP_max=0;
      pf_eeM=0, pf_eeP=0, pf_eeM_max=0, pf_eeP_max=0;
      chH_0nXn=0, chH_Xn0n=0;
      pf_etaMax=-1000, pf_etaMin=1000;
      

      // get zdc information
      for (int j=0; j<nZDC; j++)
        {
          if (zside[j]>0)
            ZDC_plus += e[j];
          if (zside[j]<0)
            ZDC_minus += e[j];
        }


      // get max pfE for hf, HE
      // PF Loop
      for (int iPF=0; iPF<nPF; iPF++)
	{
	  // HF+
	  if (pfEta->at(iPF)>3 && pfEta->at(iPF)<5.2)
	    {
	      pf_hfP = pfE->at(iPF);
	    }
	  if (pf_hfP >= pf_hfP_max) pf_hfP_max = pf_hfP;

	  // HF-
	  if (pfEta->at(iPF)>-5.2 && pfEta->at(iPF)<-3)
	    {
	      pf_hfM = pfE->at(iPF);
	    }
	  if (pf_hfM >= pf_hfM_max) pf_hfM_max = pf_hfM;

	  // HE+
	  if (pfEta->at(iPF)>2.8 && pfEta->at(iPF)<3.0)
	    {
	      if (pfId->at(iPF)==1 || pfId->at(iPF)==5)
		{
		  if (pfE->at(iPF) >= pf_heP_max) pf_heP_max = pfE->at(iPF);
		}
	    }

	  // HE-
	  if (pfEta->at(iPF)>(-3.0) && pfEta->at(iPF)<(-2.8))
	    {
	      if (pfId->at(iPF)==1 || pfId->at(iPF)==5)
		{
		  if (pfE->at(iPF) >= pf_heM_max) pf_heM_max = pfE->at(iPF);
		}
	    }
	  // EE+
	  if (pfEta->at(iPF)>2.8 && pfEta->at(iPF)<3.0)
	    {
	      if (pfId->at(iPF)==4)
		{
		  if (pfE->at(iPF) >= pf_eeP_max) pf_eeP_max = pfE->at(iPF);
		}
	    }
	  // EE-
	  if (pfEta->at(iPF)>(-3.0) && pfEta->at(iPF)<(-2.8))
	    {
	      if (pfId->at(iPF)==4)
		{
		  if (pfE->at(iPF) >= pf_eeM_max) pf_eeM_max = pfE->at(iPF);
		}
	    }
	 
	  
	} // end PF Loop

     
      


      // RAPGAP CUT
      for (int iPF=0; iPF<nPF; iPF++)
	{
	  if (pfTrkEta->at(iPF)==-99) continue;
	  if (pfTrkEta->at(iPF)>pf_etaMax) pf_etaMax=pfTrkEta->at(iPF);
	  if (pfTrkEta->at(iPF)<pf_etaMin) pf_etaMin=pfTrkEta->at(iPF);
	}
      forward_rapgap=3-pf_etaMax;
      backward_rapgap=-3-pf_etaMin;
      if (ZDC_minus>1200) h2_rapgaps_Xn0n->Fill(fabs(backward_rapgap), fabs(forward_rapgap));
      if (ZDC_plus>1200) h2_rapgaps_0nXn->Fill(fabs(backward_rapgap), fabs(forward_rapgap));
      if (ZDC_minus>1200 && (fabs(forward_rapgap)<trackGap)) continue;
      if (ZDC_plus>1200 && (fabs(backward_rapgap)<trackGap)) continue;
      if (ZDC_minus>1200 && forward_rapgap<fabs(backward_rapgap)) continue;
      if (ZDC_plus>1200 && fabs(backward_rapgap)<forward_rapgap) continue;
      
   

      // Calculate px1, py1, px2, py2
      px_1=(jetPt[0])*(cos(jetPhi[0]));
      py_1=(jetPt[0])*(sin(jetPhi[0]));
      px_2=(jetPt[1])*(cos(jetPhi[1]));
      py_2=(jetPt[1])*(sin(jetPhi[1]));

      // Define QT and PT 2-vectors
      float random=0;
      random=rand()%(2)+1;
      Q_T.Set(px_1+px_2, py_1+py_2);
      if (random==1) P_T.Set(0.5*(px_1-px_2), 0.5*(py_1-py_2));
      if (random==2) P_T.Set(0.5*(px_2-px_1), 0.5*(py_2-py_1));

      // Compute the norm of QT, PT
      QT_norm = sqrt(Q_T.X()*Q_T.X() + Q_T.Y()*Q_T.Y());
      PT_norm = sqrt(P_T.X()*P_T.X() + P_T.Y()*P_T.Y());

      // Further selections
      h2_QTPT->Fill(QT_norm,PT_norm);
      if (PT_norm<QT_norm) continue;

      h_nRef->Fill(nRef);

      // get rid of 3rd jets
      if (thirdJetCut)
	{
	  if (jetPt[2] > 20) continue;
	}

      // count number ch hadrons
      for (int iPF=0; iPF<nPF; iPF++)
	{
	  h_PF_pt->Fill(pfE->at(iPF));
	  
	  if (pfId->at(iPF)==1) pfId_chH_sum = pfId_chH_sum+1;
	  if (pfId->at(iPF)==5) pfId_H_sum = pfId_H_sum+1;
	  if ((pfId->at(iPF)==1) && (pfE->at(iPF)<0.3)) pfId_chH_sum_less300 = pfId_chH_sum_less300+1;
	}
	

      // fill rapgap hist 0nXn
      if (ZDC_minus<1200)
	{
	  h_pf_hfM_max_0nXn->Fill(pf_hfM_max);
	  h_pf_hfP_max_0nXn->Fill(pf_hfP_max);
	  
	  if (pf_hfM_max>hfM_noiseTH) continue;
	  h_heM_2p8->Fill(pf_heM_max);
	  h_eeM_2p8->Fill(pf_eeM_max);
	  if (pf_heM_max>heM_noiseTH) continue;
	  if (pf_eeM_max>eeM_noiseTH) continue;

	  // PF Loop
	  for (int iPF=0; iPF<nPF; iPF++)
	    {
	      if (pfId->at(iPF)==1)
		{
		  h_chH_eta_0nXn->Fill(pfEta->at(iPF));
		  if (pfEta->at(iPF)>(-3.0) && pfEta->at(iPF)<(-2.4))
		    {
		      chH_0nXn=chH_0nXn+1;
		    }
		}
	    } // end PF loop
	
	}
	  
      // fill rapgap hist Xn0n
      if (ZDC_minus>1200)
	{
	  h_pf_hfM_max_Xn0n->Fill(pf_hfM_max);
	  h_pf_hfP_max_Xn0n->Fill(pf_hfP_max);
	  if (pf_hfP_max>hfP_noiseTH) continue;
	  h_heP_2p8->Fill(pf_heP_max);
	  h_eeP_2p8->Fill(pf_eeP_max);
	  if (pf_heP_max>heP_noiseTH) continue;
	  if (pf_eeP_max>eeP_noiseTH) continue;

	  // PF Loop
	  for (int iPF=0; iPF<nPF; iPF++)
	    {
	      if (pfId->at(iPF)==1)
		{
		  h_chH_eta_Xn0n->Fill(pfEta->at(iPF));
		  if (pfEta->at(iPF)>(2.4) && pfEta->at(iPF)<(3.0))
		    {
		      chH_Xn0n=chH_Xn0n+1;
		    }
		}
	    } // end PF loop
	  }
      


      
      // make TLorentz Dijet objects for the good evts
      TLorentzVector jet1, jet2;
      jet1.SetPtEtaPhiM(jetPt[0],jetEta[0],jetPhi[0],jetM[0]);
      jet2.SetPtEtaPhiM(jetPt[1],jetEta[1],jetPhi[1],jetM[1]);
      float dijetRap = (jet1+jet2).Rapidity();
      float dPhi=fabs(jet1.DeltaPhi(jet2));
      if (ZDC_minus<1200) rapPrime = dijetRap*(1);
      if (ZDC_minus>=1200) rapPrime = dijetRap*(-1);
      // 5pi/6 cut
      if (dPhiCut)
	{
	  if (dPhi<5*TMath::Pi()/6) continue;
	}

     
      goodEvt = goodEvt + 1;


      // Calculate v2 stuff
      // Make unit vectors of QT, PT
      QT_unit.Set(Q_T.X()/QT_norm, Q_T.Y()/QT_norm);
      PT_unit.Set(P_T.X()/PT_norm, P_T.Y()/PT_norm);

      // cos(phi) = Dot Product of the unit vectors
      cos12 = QT_unit.X()*PT_unit.X() + QT_unit.Y()*PT_unit.Y();
      // sin(phi) = Cross Product of the unit vectors
      sin12 = QT_unit.X()*PT_unit.Y() - QT_unit.Y()*PT_unit.X();

      // Qt and Pt angles
      QT_phi = atan2(QT_unit.Y(), QT_unit.X());
      PT_phi = atan2(PT_unit.Y(), PT_unit.X());
      angleDiff = PT_phi - QT_phi;
      if (angleDiff < 0) angleDiff = angleDiff + 2*pi;

      
      // Compute the angle by using arctan2 function
      angle12 = atan2(sin12, cos12);
      float angle12_orig = angle12;
      h_trigPhi->Fill(angle12);
   

      // cos(n*phi) values
      cos_2phi = cos(2*angle12);
      float cos_2phi_orig = cos(2*angle12_orig);
      cos_4phi = cos(4*angle12);
      cos_6phi = cos(6*angle12);
      cos_8phi = cos(8*angle12);
      cos_10phi = cos(10*angle12);

      v1_total = v1_total + cos12;
      v2_total = v2_total + cos_2phi;
      v4_total = v4_total + cos_4phi;
      v6_total = v6_total + cos_6phi;
      v8_total = v8_total + cos_8phi;
      v10_total = v10_total + cos_10phi;


      h2_rap1_rap2->Fill(jet1.Rapidity(), jet2.Rapidity());
      h_dijetRapPrime->Fill(rapPrime);
      if (rapPrime>-2.4 && rapPrime<0)
	{
	  h_profV2_rapPrime1_3->Fill(QT_norm,cos_2phi);
	  h_profV2_rapPrime1_2->Fill(QT_norm,cos_2phi);
	  if (QT_norm>0 && QT_norm<12) h_phi_rapidity[0][0]->Fill(angle12);
	  if (QT_norm>12 && QT_norm<40) h_phi_rapidity[1][0]->Fill(angle12);
	}
      if (rapPrime>0 && rapPrime<0.5)
	{
	  h_profV2_rapPrime2_3->Fill(QT_norm,cos_2phi);
	  h_profV2_rapPrime2_2->Fill(QT_norm,cos_2phi);
	  if (QT_norm>0 && QT_norm<12) h_phi_rapidity[0][1]->Fill(angle12);
	  if (QT_norm>12 && QT_norm<40) h_phi_rapidity[1][1]->Fill(angle12);
	}
      if (rapPrime>0.5 && rapPrime<1.0)
	{
	  h_profV2_rapPrime3_3->Fill(QT_norm,cos_2phi);
	  h_profV2_rapPrime3_2->Fill(QT_norm,cos_2phi);
	  if (QT_norm>0 && QT_norm<12) h_phi_rapidity[0][2]->Fill(angle12);
	  if (QT_norm>12 && QT_norm<40) h_phi_rapidity[1][2]->Fill(angle12);
	}
      if (rapPrime>1.0 && rapPrime<1.5)
	{
	  h_profV2_rapPrime4_3->Fill(QT_norm,cos_2phi);
	  h_profV2_rapPrime4_2->Fill(QT_norm,cos_2phi);
	  if (QT_norm>0 && QT_norm<12) h_phi_rapidity[0][3]->Fill(angle12);
	  if (QT_norm>12 && QT_norm<40) h_phi_rapidity[1][3]->Fill(angle12);
	}
      if (rapPrime>1.5 && rapPrime<2.4)
	{
	  h_profV2_rapPrime5_3->Fill(QT_norm,cos_2phi);
	  h_profV2_rapPrime5_2->Fill(QT_norm,cos_2phi);
	  if (QT_norm>0 && QT_norm<12) h_phi_rapidity[0][4]->Fill(angle12);
	  if (QT_norm>12 && QT_norm<40) h_phi_rapidity[1][4]->Fill(angle12);
	}

      if (QT_norm>0 && QT_norm<12) h2_phi12_rap_qt0->Fill(rapPrime, angle12);
      if (QT_norm>12 && QT_norm<40) h2_phi12_rap_qt1->Fill(rapPrime, angle12);


      
	
     	

      // Fill Histograms
      h2_rapgaps_ext->Fill(backRapGap_ext(nPF, pfTrkPt, pfTrkEta, 0.2),forwardRapGap_ext(nPF, pfTrkPt, pfTrkEta, 0.2));

      /// TEEEEEEEEEEEEEEEEEEEEEEEEST /////////////
      if (ZDC_plus<1200)
	{
	  //h_pf_hfM_max_Xn0n->Fill(pf_hfM_max);
	  //h_pf_hfP_max_Xn0n->Fill(pf_hfP_max);
	  h2_rapgaps_Xn0n->Fill(backRapGap(nPF, pfTrkPt, pfTrkEta, 0.2),forwardRapGap(nPF, pfTrkPt, pfTrkEta, 0.2));
	  h_rapgap_backward->Fill(backRapGap(nPF, pfTrkPt, pfTrkEta, 0.2));
	  h_rapgap_forward->Fill(forwardRapGap(nPF, pfTrkPt, pfTrkEta, 0.2));
	  }

    
      
      h_measZ->Fill(return_z(h2_dummy, QT_norm, angle12));
      if (QT_norm>0 && QT_norm<12) h_measZ_rap0->Fill(return_z(h2_dummy_rap0, rapPrime, angle12));
      if (QT_norm>12 && QT_norm<40) h_measZ_rap1->Fill(return_z(h2_dummy_rap1, rapPrime, angle12));
      h_zvtx->Fill(vtx_z);
      h_zdcMinus->Fill(ZDC_minus);
      h_zdcPlus->Fill(ZDC_plus);
      h_jtPt_L->Fill(jetPt[0]);
      h_jtPt_S->Fill(jetPt[1]);
      h_dijetMass->Fill((jet1 + jet2).M());
      h_dijetRap->Fill(dijetRap);
      if (ZDC_plus>1200) h_dijetRap_0nXn->Fill(dijetRap);
      if (ZDC_minus>1200) h_dijetRap_Xn0n->Fill(dijetRap);
      h_PT->Fill(PT_norm);
      h_QT->Fill(QT_norm);
      h_PT_10->Fill(PT_norm);
      h_QT_10->Fill(QT_norm);
      h_QT_meas->Fill(QT_norm);
      h_QT10_meas->Fill(QT_norm);
      h_phiQTPT_12->Fill(angle12_orig);
      h_phiQTPT_12_abs->Fill(fabs(angle12_orig));
      if (0<QT_norm && QT_norm<8) h_phiQTPT_bin1->Fill(angle12_orig);
      if (8<QT_norm && QT_norm<16) h_phiQTPT_bin2->Fill(angle12_orig);
      if (16<QT_norm && QT_norm<24) h_phiQTPT_bin3->Fill(angle12_orig);
      if (24<QT_norm && QT_norm<32) h_phiQTPT_bin4->Fill(angle12_orig);
      if (32<QT_norm && QT_norm<40) h_phiQTPT_bin5->Fill(angle12_orig);
      if (0<QT_norm && QT_norm<8) h_phiQTPT_6_bin1->Fill(angle12_orig);
      if (8<QT_norm && QT_norm<16) h_phiQTPT_6_bin2->Fill(angle12_orig);
      if (16<QT_norm && QT_norm<24) h_phiQTPT_6_bin3->Fill(angle12_orig);
      if (24<QT_norm && QT_norm<32) h_phiQTPT_6_bin4->Fill(angle12_orig);
      if (32<QT_norm && QT_norm<40) h_phiQTPT_6_bin5->Fill(angle12_orig);
      h_cos2phi_meas->Fill(cos_2phi_orig);
      h_cos2phi_meas_12->Fill(cos_2phi);
      h2_cos12_QT5_meas->Fill(QT_norm, cos_2phi_orig);
      h2_QT_cos->Fill(QT_norm, cos_2phi);
      h2_QT8_cos->Fill(QT_norm, cos_2phi);
      h2_QT8_cos20->Fill(QT_norm, cos_2phi);
      h2_QT8_cos15->Fill(QT_norm, cos_2phi);
      h_profPT_5->Fill(QT_norm, PT_norm);
      h_profV1->Fill(QT_norm, cos12);
      h_profV2->Fill(QT_norm,cos_2phi);
      h_profV2_5->Fill(QT_norm,cos_2phi);
      h_profV2_8->Fill(QT_norm,cos_2phi);
      h_profV2_10->Fill(QT_norm,cos_2phi);
      h_profV4->Fill(QT_norm,cos_4phi);
      h_profV6->Fill(QT_norm,cos_6phi);
      h_profV8->Fill(QT_norm,cos_8phi);
      h_profV10->Fill(QT_norm,cos_10phi);
      h_deltaPhi->Fill(angleDiff);
      h_deltaPhi_pre->Fill(PT_phi - QT_phi);
      h_angle12->Fill(angle12);
      h_angle12_pre->Fill(atan2(sin12, cos12));
      h2_phi6_QT5->Fill(QT_norm, angle12);
      h2_phi9_QT5->Fill(QT_norm, angle12);
      h2_phiQTPT_QT->Fill(QT_norm, angle12);
      h2_phiQTPT_QT_mpi->Fill(QT_norm, angle12_orig);
      h2_cosQTPT_QT_mpi->Fill(QT_norm, cos(2*angle12_orig));
      if (QT_norm>0 && QT_norm<12) h2_cosQTPT_rap0->Fill(rapPrime, cos(2*angle12_orig));
      if (QT_norm>12 && QT_norm<40) h2_cosQTPT_rap1->Fill(rapPrime, cos(2*angle12_orig));

      if (dijetRap>-2.5 && dijetRap<-1.5)
	{
	  h_profV2_rap1->Fill(QT_norm,cos_2phi);
	  h_profV2_rap1_4->Fill(QT_norm,cos_2phi);
	}
      if (dijetRap>-1.5 && dijetRap<-1.0)
	{
	  h_profV2_rap2->Fill(QT_norm,cos_2phi);
	  h_profV2_rap2_4->Fill(QT_norm,cos_2phi);
	}
      if (dijetRap>-1.0 && dijetRap<-0.5)
	{
	  h_profV2_rap3->Fill(QT_norm,cos_2phi);
	  h_profV2_rap3_4->Fill(QT_norm,cos_2phi);
	}
      if (dijetRap>-0.5 && dijetRap<0)
	{
	  h_profV2_rap4->Fill(QT_norm,cos_2phi);
	  h_profV2_rap4_4->Fill(QT_norm,cos_2phi);
	}
      if (dijetRap>0 && dijetRap<0.5)
	{
	  h_profV2_rap5->Fill(QT_norm,cos_2phi);
	  h_profV2_rap5_4->Fill(QT_norm,cos_2phi);
	}
      if (dijetRap>0.5 && dijetRap<1.0)
	{
	  h_profV2_rap6->Fill(QT_norm,cos_2phi);
	  h_profV2_rap6_4->Fill(QT_norm,cos_2phi);
	}
      if (dijetRap>1.0 && dijetRap<1.5)
	{
	  h_profV2_rap7->Fill(QT_norm,cos_2phi);
	  h_profV2_rap7_4->Fill(QT_norm,cos_2phi);
	}
      if (dijetRap>1.5 && dijetRap<2.5)
	{
	  h_profV2_rap8->Fill(QT_norm,cos_2phi);
	  h_profV2_rap8_4->Fill(QT_norm,cos_2phi);
	}
      
      if (vtx_z>-20 && vtx_z<-5) h_profV2_vtx1->Fill(QT_norm,cos_2phi);
      if (vtx_z>-5 && vtx_z<-2.5) h_profV2_vtx2->Fill(QT_norm,cos_2phi);
      if (vtx_z>-2.5 && vtx_z<0) h_profV2_vtx3->Fill(QT_norm,cos_2phi);
      if (vtx_z>0 && vtx_z<2.5) h_profV2_vtx4->Fill(QT_norm,cos_2phi);
      if (vtx_z>2.5 && vtx_z<5) h_profV2_vtx5->Fill(QT_norm,cos_2phi);
      if (vtx_z>5 && vtx_z<20) h_profV2_vtx6->Fill(QT_norm,cos_2phi);

      h_nPF->Fill(nPF);
      h_nPF_H->Fill(pfId_chH_sum+pfId_H_sum);
      h_nPF_ch->Fill(pfId_chH_sum);
      h_ntrk_vtx->Fill(pfTrkPt->size(), vtx_z);
      h_nPF_vtx->Fill(nPF, vtx_z);
      
      h_chH_vtx->Fill(pfId_chH_sum, vtx_z);
      if (pfId_chH_sum<4) h_zvtx_3chH->Fill(vtx_z);
      if (pfId_chH_sum>=4) h_zvtx_4chH_up->Fill(vtx_z);
      if (pfId_chH_sum>=4 && pfId_chH_sum_less300==0) h_zvtx_4chH_up_300->Fill(vtx_z); 
         
  
      
    } // end event loop 

  avg_v2 = v2_total / goodEvt;
  avg_v4 = v4_total / goodEvt;
  avg_v6 = v6_total / goodEvt;
  avg_v8 = v8_total / goodEvt;
  avg_v10 = v10_total / goodEvt;
  cout << "overall v2=" << avg_v2 << endl;
  cout << "overall v4=" << avg_v4 << endl;
  cout << "overall v6=" << avg_v6 << endl;
  cout << "overall v8=" << avg_v8 << endl;
  cout << "overall v10=" << avg_v10 << endl;
  cout << goodEvt << endl;

  // write histograms
  h_measZ->Write();
  h_measZ_rap0->Write();
  h_measZ_rap1->Write();
  h_zvtx->Write();
  h_zdcMinus->Write();
  h_zdcPlus->Write();
  h_trigPhi->Write();
  h_ee_0nXn->Write();
  h_ee_Xn0n->Write();
  h_nRef->Write();
  h_jtPt_L->Write();
  h_jtPt_S->Write();
  h_dijetMass->Write();
  h2_rap1_rap2->Write();
  h_dijetRap->Write();
  h_dijetRap_0nXn->Write();
  h_dijetRap_Xn0n->Write();
  h_PT->Write();
  h_QT->Write();
  h_PT_10->Write();
  h_QT_10->Write();
  h_QT_meas->Write();
  h_QT10_meas->Write();
  h_phiQTPT_12->Write();
  h_phiQTPT_12_abs->Write();
  h_cos2phi_meas_12->Write();
  h_cos2phi_meas->Write();
  h2_cos12_QT5_meas->Write();
  h2_QTPT->Write();
  h_deltaPhi->Write();
  h_deltaPhi_pre->Write();
  h_angle12->Write();
  h_angle12_pre->Write();
  h2_rapgaps_ext->Write();
  h2_rapgaps_0nXn->Write();
  h2_rapgaps_Xn0n->Write();
  h_rapgap_backward->Write();
  h_rapgap_forward->Write();
  
  h2_QT_cos->Write();
  h2_QT8_cos->Write();
  h2_QT8_cos20->Write();
  h2_QT8_cos15->Write();
  h_profPT_5->Write();
  h_profV2->Write();
  h_profV2_5->Write();
  h_profV2_8->Write();
  h_profV2_10->Write();
  h_profV4->Write();
  h_profV6->Write();
  h_profV8->Write();
  h_profV10->Write();
  h_profV2_rap1->Write();
  h_profV2_rap2->Write();
  h_profV2_rap3->Write();
  h_profV2_rap4->Write();
  h_profV2_rap5->Write();
  h_profV2_rap6->Write();
  h_profV2_rap7->Write();
  h_profV2_rap8->Write();
  h_profV2_vtx1->Write();
  h_profV2_vtx2->Write();
  h_profV2_vtx3->Write();
  h_profV2_vtx4->Write();
  h_profV2_vtx5->Write();
  h_profV2_vtx6->Write();
  h_profV2_rapPrime1_3->Write();
  h_profV2_rapPrime2_3->Write();
  h_profV2_rapPrime3_3->Write();
  h_profV2_rapPrime4_3->Write();
  h_profV2_rapPrime5_3->Write();
  h_profV2_rapPrime1_2->Write();
  h_profV2_rapPrime2_2->Write();
  h_profV2_rapPrime3_2->Write();
  h_profV2_rapPrime4_2->Write();
  h_profV2_rapPrime5_2->Write();
  h2_phi12_rap_qt0->Write();
  h2_phi12_rap_qt1->Write();
  h2_phi6_QT5->Write();
  h2_phi9_QT5->Write();
  h2_phiQTPT_QT->Write();
  h2_phiQTPT_QT_mpi->Write();
  h2_cosQTPT_QT_mpi->Write();
  h2_cosQTPT_rap0->Write();
  h2_cosQTPT_rap1->Write();
  h_phiQTPT_bin1->Write();
  h_phiQTPT_bin2->Write();
  h_phiQTPT_bin3->Write();
  h_phiQTPT_bin4->Write();
  h_phiQTPT_bin5->Write();
  h_phiQTPT_6_bin1->Write();
  h_phiQTPT_6_bin2->Write();
  h_phiQTPT_6_bin3->Write();
  h_phiQTPT_6_bin4->Write();
  h_phiQTPT_6_bin5->Write();
  h_pf_hfM_max_0nXn->Write();
  h_pf_hfM_max_Xn0n->Write();
  h_pf_hfP_max_0nXn->Write();
  h_pf_hfP_max_Xn0n->Write();
  for (int qtBin=0; qtBin<2; qtBin++)
    {
      for (int rapBin=0; rapBin<5; rapBin++)
	{
	  h_phi_rapidity[qtBin][rapBin]->Write();
	} // end rap loop
    } // end qtBin loop
  h_heP_2p8->Write();
  h_heM_2p8->Write();
  h_eeP_2p8->Write();
  h_eeM_2p8->Write();
  h_chH_eta_0nXn->Write();
  h_chH_eta_Xn0n->Write();
  newfile.Write();
  h_nPF->Write();
  h_nPF_H->Write();
  h_nPF_ch->Write();
  h_ntrk_vtx->Write();
  h_nPF_vtx->Write();
  h_chH_vtx->Write();
  h_zvtx_3chH->Write();
  h_zvtx_4chH_up->Write();
  h_PF_pt->Write();
  h_zvtx_4chH_up_300->Write();

} // end code
