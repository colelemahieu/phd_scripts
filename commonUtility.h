// commonUtility.h

#include <iostream>
#include <fstream>
#include <ctime>
#include "Riostream.h"
#include <TGraphAsymmErrors.h>
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
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TVirtualFitter.h"
#include <TBox.h>

void drawText(float xp, float yp, const char *text, int textColor=kBlack, double textSize=22){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(43);
  tex->SetTextSize(textSize);
  tex->SetTextColor(textColor);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw();
}

void setMargins(TCanvas *canvas)
{
  canvas->SetRightMargin(0.125);
  canvas->SetLeftMargin(0.125);
  canvas->SetBottomMargin(0.125);
}

void drawCMS_p(TCanvas *canvas)
{
  canvas->SetRightMargin(0.125);
  canvas->SetLeftMargin(0.125);
  canvas->SetBottomMargin(0.125);
  drawText(0.1325, 0.915, "#bf{CMS} #it{work in progress}", kBlack, 26);
  drawText(0.53, 0.915, "PbPb #sqrt{s}=5.02 TeV; 2018", kBlack, 26);
}

void drawCMS_lumi(TCanvas *canvas)
{
  canvas->SetRightMargin(0.125);
  canvas->SetLeftMargin(0.125);
  canvas->SetBottomMargin(0.125);
  drawText(0.1325, 0.915, "#bf{CMS} #it{Preliminary}", kBlack, 26);
  drawText(0.49, 0.915, "PbPb #sqrt{s}=5.02 TeV; 1.67 nb^{-1}", kBlack, 26);
}

void lumiLabel()
{
  drawText(0.13, 0.93, "#bf{CMS Preliminary}");
  drawText(0.71, 0.93, "5.02 TeV (2018)");
}

void drawErrorBox(float x1,float y1, float x2, float y2, int theColor=kSpring)
{
	TBox* tt = new TBox(x1,y1,x2,y2);
	tt->SetFillColor(theColor);
	tt->SetFillStyle(3001);
	tt->Draw();
}

float find_hist_peak(TH1 *hist)
{
  int maxBin = hist->GetMaximumBin();
  float peak = hist->GetBinCenter(maxBin);
  
  return (peak);
}

float find_hist_chi2(TH1 *h1, TH1 *h2)
{
  int binNumber = h1->GetNbinsX();
  float h1_value=0, h2_value=0, binDiff=0;
  float h1_error=0, h2_error=0, totError=0;
  float chi2=0, chi_num=0, chi_denom=0;

  for (int iBin=1; iBin<binNumber+1; iBin++)
    {
      h1_value=h1->GetBinContent(iBin);
      h2_value=h2->GetBinContent(iBin);
      binDiff=fabs(h1_value-h2_value);

      h1_error=h1->GetBinError(iBin);
      h2_error=h2->GetBinError(iBin);
      totError=sqrt(h1_error*h1_error + h2_error*h2_error);

      chi_num = binDiff*binDiff;
      chi_denom = totError*totError;
      chi2 = chi2 + chi_num/chi_denom;
    } // end bins loop

  return chi2;
}


TH1F* makeResidual_histogram(TH1 *h1, TH1 *h2)
{
  const int binNumber = h1->GetNbinsX();
  const int xMin = h1->GetXaxis()->GetXmin();
  const int xMax = h2->GetXaxis()->GetXmax();
  float h1_content=0, h2_content=0;
  float h1_error=0, h2_error=0, residual_error=0;
  
  TH1F *h_residual = new TH1F("h_residual", "", binNumber, xMin, xMax);

  // loop thru Bins
  for (int iBin=1; iBin<binNumber+1; iBin++)
    {
      h1_content=h1->GetBinContent(iBin);
      h2_content=h2->GetBinContent(iBin);

      h1_error=h1->GetBinError(iBin);
      h2_error=h2->GetBinError(iBin);
      residual_error = sqrt(h1_error*h1_error + h2_error*h2_error);
      
      
      h_residual->SetBinContent(iBin, (h2_content - h1_content));
      h_residual->SetBinError(iBin, residual_error);
    } // end bin loop

  return h_residual;
  
}

int jetsProperOrder(float genPhi_0, float genEta_0, float recoPhi_0, float recoEta_0, float recoPhi_1, float recoEta_1)
{
  // find distance between gen0 and reco0
  float phiDiff_00 = genPhi_0 - recoPhi_0;
  float etaDiff_00 = genEta_0 - recoEta_0;
  float jetDist_00 = sqrt(phiDiff_00*phiDiff_00 + etaDiff_00*etaDiff_00);

  // find distance between gen0 and reco1
  float phiDiff_01 = genPhi_0 - recoPhi_1;
  float etaDiff_01 = genEta_0 - recoEta_1;
  float jetDist_01 = sqrt(phiDiff_01*phiDiff_01 + etaDiff_01*etaDiff_01);

  int indexMin;
  if (jetDist_00 < jetDist_01) indexMin=1;
  if (jetDist_00 > jetDist_01) indexMin=0;
      
  return indexMin;
}



float integralThreshold(TH1 *h1, float threshold)
{
  h1->Scale(1/h1->Integral());
  int binNumber = h1->GetNbinsX();
  float currentIntegral;
  float xValue=0;

  for (int bin=1; bin<binNumber+1; bin++)
    {
      currentIntegral = h1->Integral(1,bin);
      xValue = h1->GetXaxis()->GetBinCenter(bin);

      if (currentIntegral>threshold) break;
    }
    
  return xValue;
}




bool check_HFp_Ex(int nTower, vector<float> *energy, vector<float> *caloEta)
{
  // noise threshold
  float hf_P = 7.2;
  int flag=1;

  for (int tower=0; tower<nTower; tower++)
    {

      // HFp region
      if (caloEta->at(tower)<3.0) continue;
      
      if (energy->at(tower) > hf_P)
	{
	  flag=0;
	  break;
	}
      if (energy->at(tower) <= hf_P)
	{
	  flag=1;
	}
    }

  return flag;
}



bool check_HFm_Ex(int nTower, vector<float> *energy, vector<float> *caloEta)
{
  // noise threshold
  float hf_M = 7.5;  
  int flag=1;

  for (int tower=0; tower<nTower; tower++)
    {

      // HFm region
      if (caloEta->at(tower)>-3.0) continue;
      
      if (energy->at(tower) > hf_M)
	{
	  flag=0;
	  break;
	}
      if (energy->at(tower) <= hf_M)
	{
	  flag=1;
	}
    }

  return flag;
}

bool check_HE_Ex(int nTower, vector<float> *energy, vector<float> *caloEta)
{
  // noise threshold
  float he = 2.4;
  int flag=1;

  for (int tower=0; tower<nTower; tower++)
    {

      // HE region
      if (caloEta->at(tower)<1.3 || caloEta->at(tower)>3) continue;
      
      if (energy->at(tower) > he)
	{
	  flag=0;
	  break;
	}
      if (energy->at(tower) <= he)
	{
	  flag=1;
	}
    }

  return flag;
}

bool check_EE_Ex(int nTower, vector<float> *energy, vector<float> *caloEta)
{
  // noise threshold
  float ee = 7.5;
  int flag=1;

  for (int tower=0; tower<nTower; tower++)
    {

      // EE region
      if (caloEta->at(tower)<1.48 || caloEta->at(tower)>3.0) continue;
      
      if (energy->at(tower) > ee)
	{
	  flag=0;
	  break;
	}
      if (energy->at(tower) <= ee)
	{
	  flag=1;
	}
    }

  return flag;
}


bool checkEx(float energy, float noiseTh)
{
  int flag=1;

  if (energy > noiseTh)
    {
      flag=0;
    }

  return flag;
}




float getYmax(TH1 *h1, TH1 *h2, TH1 *h3, TH1 *h4){
  float yMax1 = h1->GetMaximum();
  float yMax2 = h2->GetMaximum();
  float yMax3 = h3->GetMaximum();
  float yMax4 = h4->GetMaximum();
  float yMax = yMax1;
  if (yMax2>yMax1) yMax = yMax2;
  if (yMax3>yMax2) yMax = yMax3;
  if (yMax4>yMax3) yMax = yMax4;

  return (1.1*yMax);  
}

void makeDotHist(TH1 *h1, int markerColor=kBlack){
  h1->SetMarkerStyle(kFullDotLarge);
  h1->SetMarkerColor(markerColor);
  h1->SetStats(0);
}

void modDotHist(TH1 *h1, int markerColor=kBlack){
  h1->SetMarkerStyle(kFullDotLarge);
  h1->SetMarkerColor(markerColor);
  h1->Scale(1/(h1->Integral()));
  h1->SetStats(0);
}

void modBarHist(TH1 *h1, int markerColor=kBlack){
  h1->SetLineColor(markerColor);
  h1->Scale(1/(h1->Integral()));
  h1->SetLineWidth(1);
  h1->SetStats(0);
}

void makeTriHist(TH1 *h1, int markerColor=kBlack){
  h1->SetMarkerStyle(kFullTriangleUp);
  h1->SetMarkerColor(markerColor);
  h1->SetStats(0);
}

void drawLine(float x1, float y1, float x2, float y2){
  TLine *l = new TLine(x1,y1,x2,y2);
  //l->SetLineWidth(3);
  l->SetLineColor(15);
  l->SetLineWidth(4);
  l->Draw();
}

void drawDashLine(float x1, float y1, float x2, float y2){
  TLine *l = new TLine(x1,y1,x2,y2);
  //l->SetLineColor(kGray+1);
  l->SetLineColor(kBlack);
  l->SetLineWidth(2);
  l->SetLineStyle(9);
  l->Draw();
}

void drawRedDashLine(float x1, float y1, float x2, float y2){
  TLine *l = new TLine(x1,y1,x2,y2);
  //l->SetLineColor(kGray+1);
  l->SetLineColor(kRed);
  l->SetLineWidth(4);
  l->SetLineStyle(9);
  l->Draw();
}

void drawLegend(TH1 *h1, TH1 *h2, float x1, float x2, float y1, float y2, const char *text1, const char *text2, int nCol=1){
  TLegend* legend = new TLegend(x1,y1,x2,y2);
  legend->SetNColumns(nCol);
  legend->AddEntry(h1, text1, "lpe");
  legend->AddEntry(h2, text2, "lpe");
  legend->SetTextSize(0.025);
  legend->SetBorderSize(0);
  legend->Draw("SAME");
  }


void drawLegendAsymm(TH1 *h1, TGraphAsymmErrors *a1, float x1, float x2, float y1, float y2, const char *text1, const char *text2, int nCol=1){
  TLegend* legend = new TLegend(x1,y1,x2,y2);
  legend->SetNColumns(nCol);
  legend->AddEntry(h1, text1, "pe");
  legend->AddEntry(a1, text2, "pe");
  legend->SetTextSize(0.045);
  legend->SetBorderSize(0);
  legend->Draw("SAME");
  }



// lumi legend
void lumiLegend(TH1 *h1, TH1 *h2, float x1, float x2, float y1, float y2, string text1, string text2, int nCol=1){
  TLegend* legend = new TLegend(x2,y2,x1,y1);
  legend->SetNColumns(nCol);
  legend->AddEntry(h1, text1.c_str(), "f");   // l f p e
  legend->AddEntry(h2, text2.c_str(), "lep");
  legend->SetTextSize(0.05);
  legend->SetBorderSize(0);
  legend->Draw("SAME");
}

void drawLegendStr(TH1 *h1, TH1 *h2, float x1, float x2, float y1, float y2, string text1, string text2, int nCol=1){
  TLegend* legend = new TLegend(x2,y2,x1,y1);
  legend->SetNColumns(nCol);
  legend->AddEntry(h1, text1.c_str(), "lep");   // l f p e
  legend->AddEntry(h2, text2.c_str(), "lep");
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->Draw("SAME");
  }


void drawLegend_TH1_TProfile(TH1 *h1, TProfile *p1, float x1, float x2, float y1, float y2, string text1, string text2, int nCol=1){
  TLegend* legend = new TLegend(x2,y2,x1,y1);
  legend->SetNColumns(nCol);
  legend->AddEntry(h1, text1.c_str(), "pe");   // l f p 
  legend->AddEntry(p1, text2.c_str(), "pe");
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->Draw("SAME");
}



void drawLegend_TProfile_TH1_TProfile(TProfile *p1, TH1 *h2, TProfile *p2, float x1, float x2, float y1, float y2, string text1, string text2, string text3, int nCol=1){
  TLegend* legend = new TLegend(x2,y2,x1,y1);
  legend->SetNColumns(nCol);
  legend->AddEntry(p1, text1.c_str(), "p");   // l f p
  legend->AddEntry(h2, text2.c_str(), "p");
  legend->AddEntry(p2, text3.c_str(), "p");
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->Draw("SAME");
}



void drawLegendTF1(TF1 *f1, TF1 *f2, float x1, float x2, float y1, float y2, string text1, string text2, int nCol=1){
  TLegend* legend = new TLegend(x2,y2,x1,y1);
  legend->SetNColumns(nCol);
  legend->AddEntry(f1, text1.c_str(), "l");   // l f p 
  legend->AddEntry(f2, text2.c_str(), "l");
  legend->SetTextSize(0.05);
  legend->SetBorderSize(0);
  legend->Draw("SAME");
}

void drawLegend_TH1_3TF1(TH1 *h1, TF1 *fit, TF1 *f1, TF1 *f2, float x1, float x2, float y1, float y2, string text1, string text2, string text3, string text4, int nCol=1){
  TLegend* legend = new TLegend(x2,y2,x1,y1);
  legend->SetNColumns(nCol);
  legend->AddEntry(h1, text1.c_str(), "pe");   // l f p 
  legend->AddEntry(fit, text2.c_str(), "l");
  legend->AddEntry(f1, text3.c_str(), "l");   
  legend->AddEntry(f2, text4.c_str(), "l");
  legend->SetTextSize(0.05);
  legend->SetBorderSize(0);
  legend->Draw("SAME");
}

void drawLegendTGraph(TGraph *g1, TGraph *g2, float x1, float x2, float y1, float y2, string text1, string text2, int nCol=1){
  TLegend* legend = new TLegend(x2,y2,x1,y1);
  legend->SetNColumns(nCol);
  legend->AddEntry(g1, text1.c_str(), "P");   // l f p 
  legend->AddEntry(g2, text2.c_str(), "P");
  legend->SetTextSize(0.025);
  legend->SetBorderSize(0);
  legend->Draw("SAME");
}

void drawLegendProfile(TProfile *p1, TProfile *p2, float x1, float x2, float y1, float y2, string text1, string text2, int nCol=1){
  TLegend* legend = new TLegend(x2,y2,x1,y1);
  legend->SetNColumns(nCol);
  legend->AddEntry(p1, text1.c_str(), "F");   // l f p 
  legend->AddEntry(p2, text2.c_str(), "lep");
  legend->SetTextSize(0.055);
  legend->SetBorderSize(0);
  legend->Draw("SAME");
}

void drawLegendProfile_3(TProfile *p1, TProfile *p2, TProfile *p3, float x1, float x2, float y1, float y2, string text1, string text2, string text3, int nCol=1){
  TLegend* legend = new TLegend(x2,y2,x1,y1);
  legend->SetNColumns(nCol);
  legend->AddEntry(p1, text1.c_str(), "f");   // l f p 
  legend->AddEntry(p2, text2.c_str(), "f");
  legend->AddEntry(p3, text3.c_str(), "P");
  legend->SetTextSize(0.035);
  legend->SetBorderSize(0);
  legend->Draw("SAME");
}

void drawLegendProfile_4(TProfile *p1, TProfile *p2, TProfile *p3, TProfile *p4, float x1, float x2, float y1, float y2, string text1, string text2, string text3, string text4, int nCol=1){
  TLegend* legend = new TLegend(x2,y2,x1,y1);
  legend->SetNColumns(nCol);
  legend->AddEntry(p1, text1.c_str(), "leP");   // l f p 
  legend->AddEntry(p2, text2.c_str(), "leP");
  legend->AddEntry(p3, text3.c_str(), "leP");
  legend->AddEntry(p4, text4.c_str(), "leP");
  legend->SetTextSize(0.035);
  legend->SetBorderSize(0);
  legend->Draw("SAME");
}

void drawLegendFit(TH1 *h1, TF1 *fit, float x1, float x2, float y1, float y2, string text1, string text2, int nCol=1){
  TLegend* legend = new TLegend(x2,y2,x1,y1);
  legend->SetNColumns(nCol);
  legend->AddEntry(h1, text1.c_str(), "f");   // l f p 
  legend->AddEntry(fit, text2.c_str(), "l");
  legend->SetTextSize(0.05);
  legend->SetBorderSize(0);
  legend->Draw("SAME");
}

void drawLegendStr3(TH1 *h1, TH1 *h2, TH1 *h3, float x1, float x2, float y1, float y2, string text1, string text2, string text3, int nCol=1){
  TLegend* legend = new TLegend(x2,y2,x1,y1);
  legend->SetNColumns(nCol);
  legend->AddEntry(h1, text1.c_str(), "lep"); // l f p
  legend->AddEntry(h2, text2.c_str(), "le");
  legend->AddEntry(h3, text3.c_str(), "le");
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->Draw("SAME");
}


void drawLegendStr4(TH1 *h1, string leg1, TH1 *h2, string leg2, TH1 *h3, string leg3, TH1 *h4, string leg4, float x1, float x2, float y1, float y2, string text1, string text2, string text3, string text4, int nCol=1){
  TLegend* legend = new TLegend(x2,y2,x1,y1);
  legend->SetNColumns(nCol);
  legend->AddEntry(h1, text1.c_str(), leg1.c_str()); // l p f
  legend->AddEntry(h2, text2.c_str(), leg2.c_str());
  legend->AddEntry(h3, text3.c_str(), leg3.c_str());
  legend->AddEntry(h4, text4.c_str(), leg4.c_str());
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->Draw("SAME");
}




void drawLegendStr5(TH1 *h1, string leg1, TH1 *h2, string leg2, TH1 *h3, string leg3, TH1 *h4, string leg4, TH1 *h5, string leg5, float x1, float x2, float y1, float y2, string text1, string text2, string text3, string text4, string text5, int nCol=1){
  TLegend* legend = new TLegend(x2,y2,x1,y1);
  legend->SetNColumns(nCol);
  legend->AddEntry(h1, text1.c_str(), leg1.c_str()); // l p f
  legend->AddEntry(h2, text2.c_str(), leg2.c_str());
  legend->AddEntry(h3, text3.c_str(), leg3.c_str());
  legend->AddEntry(h4, text4.c_str(), leg4.c_str());
  legend->AddEntry(h5, text5.c_str(), leg5.c_str());
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->Draw("SAME");
}




void histComp(TH1 *h1, TH1 *h2, const char *title, TCanvas *canvas, const char *label1, const char *label2, float x1=0.3, float x2=0.7, float y1=0.84, float y2=0.85){
  canvas->cd();
  float ymax = h1->GetMaximum();
  h1->SetMaximum(ymax*1.15);
  h1->SetTitle(title);
  h1->Draw("hist");
  h2->Draw("histsame");

  drawLegend(h1, h2, x1, x2, y1, y2, label1, label2, 2);
}


// input scaled histograms to make nice plot
void makeRatioPlot(TH1F *h1, TH1F *h2, const char *titleX, const char *titleY, TCanvas *canvas, string label1, string label2, float x1=0.60, float x2=0.80, float y1=0.65, float y2=0.80){
  
  TH1F *h_ratio = new TH1F(*h2);
  h_ratio->Divide(h1);
  h_ratio->SetMarkerColor(kBlack);
  h_ratio->SetMarkerStyle(kFullDotLarge);

  TPad *main = new TPad("main", "new pad", 0.01, 0.315, 0.99, 0.99);
  TPad *ratio = new TPad("ratio", "new pad 2", 0.01, 0.1, 0.99, 0.38);

  TPad *titlePad = new TPad("all","all",0,0,1,1);
  titlePad->SetFillStyle(4000); 
  titlePad->Draw();

  canvas->cd();
  main->Draw();
  ratio->Draw();

  main->cd();
  main->SetLeftMargin(0.125);
  main->SetBottomMargin(0.11);
  main->SetLogy();
  main->SetTickx(1);
  main->SetTicky(1);
  float ymax = h2->GetMaximum();
  h1->SetStats(0);
  h1->SetMaximum(ymax*1.3);
  // log scale y
  h1->SetMaximum(ymax*5);
  //h1->GetYaxis()->SetRangeUser(2,ymax*5);
  //h1->GetXaxis()->SetRangeUser(0,20); 
  // dirap
  //h1->SetMaximum(ymax*21);
  h1->GetXaxis()->SetLabelSize(0.0000000000001);
  h1->GetYaxis()->SetTitle(titleY);
  h1->GetYaxis()->SetTitleSize(0.07);
  h1->GetYaxis()->SetLabelSize(0.04);
  // rapidity
  h1->GetYaxis()->SetTitleOffset(0.9);
  h1->SetFillColor(kBlue-9);
  h2->SetLineColor(kBlack);
  h2->SetLineWidth(2);
  // plot overflow bin data
  /*int lastBin_1 = h1->GetNbinsX();
  double overflow_1 = h1->GetBinContent(lastBin_1 + 1);
  h1->SetBinContent(lastBin_1, h1->GetBinContent(lastBin_1) + overflow_1);
  // plot overflow bin mc
  int lastBin_2 = h2->GetNbinsX();
  double overflow_2 = h2->GetBinContent(lastBin_2 + 1);
  h2->SetBinContent(lastBin_2, h2->GetBinContent(lastBin_2) + overflow_2);*/
  h1->Draw("hist");
  h2->Draw("same pE");
  drawText(.13, .915, "#bf{CMS} #it{Preliminary}", kBlack, 30);
  drawText(0.5, 0.915, "PbPb 1.43 nb^{#minus1} (5.02 TeV)", kBlack, 30);
  
  // drawLegend 
  //lumiLegend(h1, h2, x1+0.1, x2+0.075, y1+0.05, y2+0.05, label1, label2);
  // tick marks
  //lumiLegend(h1, h2, x1+0.04, x2+0.04, y1+0.03, y2+0.03, label1, label2);
  // tick marks rap1, rap2
  lumiLegend(h1, h2, x1+0.06, x2+0.06, y1+0.03, y2+0.03, label1, label2);

  ratio->cd();
  ratio->SetTopMargin(0.035);
  ratio->SetLeftMargin(0.125);
  ratio->SetBottomMargin(0.15);
  ratio->SetTickx(1);
  ratio->SetTicky(1);
  //ratio->SetLogx();
  h_ratio->SetStats(0);
  h_ratio->SetLineColor(kBlack);
  h_ratio->SetTitle("");
  h_ratio->GetXaxis()->SetLabelSize(0.09);
  //h_ratio->GetXaxis()->SetTitleSize(0.08);
  //h_ratio->GetXaxis()->SetTitleOffset(1.0);
  //h_ratio->GetXaxis()->SetTitle(titleX);
  h_ratio->GetYaxis()->SetLabelSize(0.09);
  h_ratio->GetYaxis()->SetNdivisions(505);
  string ratio_label = label1;
  //ratio_label = ratio_label + " / " + label2;
  ratio_label = "Data/MC";
  float ymax_r = h_ratio->GetMaximum();
  float ymin_r = h_ratio->GetMinimum();
  // range for rap
  //h_ratio->GetYaxis()->SetRangeUser(0.7,1.8);
  // range for dirap
  //h_ratio->GetYaxis()->SetRangeUser(0,7);
  // pt1, pt2
  //h_ratio->GetYaxis()->SetRangeUser(0.2,1.8);
  // dimass range
  //h_ratio->GetYaxis()->SetRangeUser(0.7,1.3);
  //h_ratio->GetXaxis()->SetRangeUser(0,20);
  // range for diPt
  //h_ratio->GetYaxis()->SetRangeUser(0,12);
  // range for aco
  //h_ratio->GetYaxis()->SetRangeUser(0.4,50);
  // ext dipt
  h_ratio->GetYaxis()->SetRangeUser(0,11);
  h_ratio->GetYaxis()->SetTitle(ratio_label.c_str());
  h_ratio->GetYaxis()->SetTitleOffset(0.38);
  h_ratio->GetYaxis()->SetTitleSize(0.155);
  h_ratio->SetMarkerSize(1.2);
  ratio->SetBottomMargin(0.15);
  h_ratio->Draw("same hist p");

  float xmax = h_ratio->GetXaxis()->GetXmax();
  float xmin = h_ratio->GetXaxis()->GetXmin();
  TLine *line = new TLine(xmin,1,xmax,1);
  line->SetLineColor(kBlack);
  line->SetLineWidth(3);
  line->Draw("LSAME");
  h_ratio->Draw("same hist pE");

  titlePad->cd();
  TLatex *xlabel = new TLatex();
  xlabel->SetTextSize(0.045);
  // position for aco, pt1, pt2, diPt, dimass
  xlabel->DrawLatex(.7,.05, titleX);
  // position for rap1, rap2, dirap
  //xlabel->DrawLatex(.825,.05, titleX);
  // position for aco
  //xlabel->DrawLatex(.675,.055, titleX);
  

  canvas->SaveAs("test.png");
}


void makePtPlot(TH1F *h1, TH1F *h2, const char *titleX, const char *titleY, TCanvas *canvas, string label1, string label2, float x1=0.60, float x2=0.80, float y1=0.65, float y2=0.80){
  
  TH1F *h_ratio = new TH1F(*h2);
  h_ratio->Divide(h1);
  h_ratio->SetMarkerColor(kBlack);
  h_ratio->SetMarkerStyle(kFullDotLarge);

  TH1F *h_blank = new TH1F("h_blank","h_blank",10,3,20);
  h_blank->SetTitle("");
  TH1F *h_blank2 = new TH1F("h_blank2","h_blank2",10,3,20);
  h_blank2->SetTitle("");

  TPad *main = new TPad("main", "new pad", 0.01, 0.315, 0.99, 0.99);
  TPad *ratio = new TPad("ratio", "new pad 2", 0.01, 0.1, 0.99, 0.38);

  TPad *titlePad = new TPad("all","all",0,0,1,1);
  titlePad->SetFillStyle(4000); 
  titlePad->Draw();

  canvas->cd();
  main->Draw();
  ratio->Draw();

  main->cd();
  main->SetLeftMargin(0.125);
  main->SetBottomMargin(0.11);
  main->SetLogy();
  main->SetTickx(1);
  main->SetTicky(1);
  float ymax = h2->GetMaximum();
  h_blank->SetStats(0);
  // log scale y
  h_blank->GetYaxis()->SetRangeUser(10, ymax*5);
  h_blank->GetXaxis()->SetLabelSize(0.0000000000001);
  h_blank->GetYaxis()->SetTitle(titleY);
  h_blank->GetYaxis()->SetTitleSize(0.07);
  h_blank->GetYaxis()->SetLabelSize(0.04);
  // rapidity
  h1->GetYaxis()->SetTitleOffset(0.9);
  h1->SetFillColor(kBlue-9);
  h2->SetLineColor(kBlack);
  h2->SetLineWidth(2);
  // plot overflow bin data
  int lastBin_1 = h1->GetNbinsX();
  double overflow_1 = h1->GetBinContent(lastBin_1 + 1);
  h1->SetBinContent(lastBin_1, h1->GetBinContent(lastBin_1) + overflow_1);
  // plot overflow bin mc
  int lastBin_2 = h2->GetNbinsX();
  double overflow_2 = h2->GetBinContent(lastBin_2 + 1);
  h2->SetBinContent(lastBin_2, h2->GetBinContent(lastBin_2) + overflow_2);
  h_blank->Draw("hist");
  h1->Draw("same hist");
  h2->Draw("same pE");
  drawText(.13, .915, "#bf{CMS} #it{Preliminary}", kBlack, 30);
  drawText(0.5, 0.915, "PbPb 1.43 nb^{#minus1} (5.02 TeV)", kBlack, 30);
  
  // drawLegend 
  lumiLegend(h1, h2, x1+0.06, x2+0.06, y1+0.03, y2+0.03, label1, label2);

  ratio->cd();
  ratio->SetTopMargin(0.035);
  ratio->SetLeftMargin(0.125);
  ratio->SetBottomMargin(0.15);
  ratio->SetTickx(1);
  ratio->SetTicky(1);
  h_blank2->SetStats(0);
  h_ratio->SetLineColor(kBlack);
  h_blank2->SetTitle("");
  h_blank2->GetXaxis()->SetLabelSize(0.09);
    h_blank2->GetYaxis()->SetLabelSize(0.09);
  h_blank2->GetYaxis()->SetNdivisions(505);
  string ratio_label = label1;
  //ratio_label = ratio_label + " / " + label2;
  ratio_label = "Data/MC";
  float ymax_r = h_ratio->GetMaximum();
  float ymin_r = h_ratio->GetMinimum();
  // pt1, pt2
  h_blank2->GetYaxis()->SetRangeUser(0.6,1.4);
  h_blank2->GetYaxis()->SetTitle(ratio_label.c_str());
  h_blank2->GetYaxis()->SetTitleOffset(0.38);
  h_blank2->GetYaxis()->SetTitleSize(0.155);
  h_ratio->SetMarkerSize(1.2);
  ratio->SetBottomMargin(0.15);
  h_blank2->Draw();
  h_ratio->Draw("same hist p");

  float xmax = h_blank2->GetXaxis()->GetXmax();
  float xmin = h_blank2->GetXaxis()->GetXmin();
  TLine *line = new TLine(xmin,1,xmax,1);
  line->SetLineColor(kBlack);
  line->SetLineWidth(3);
  line->Draw("LSAME");
  h_ratio->Draw("same hist pE");

  titlePad->cd();
  TLatex *xlabel = new TLatex();
  xlabel->SetTextSize(0.045);
  // position for pt1, pt2
  xlabel->DrawLatex(.7,.05, titleX);

  canvas->SaveAs("test.png");
}


void makeRapPlot(TH1F *h1, TH1F *h2, const char *titleX, const char *titleY, TCanvas *canvas, string label1, string label2, float x1=0.60, float x2=0.80, float y1=0.65, float y2=0.80){
  
  TH1F *h_ratio = new TH1F(*h2);
  h_ratio->Divide(h1);
  h_ratio->SetMarkerColor(kBlack);
  h_ratio->SetMarkerStyle(kFullDotLarge);

  TPad *main = new TPad("main", "new pad", 0.01, 0.315, 0.99, 0.99);
  TPad *ratio = new TPad("ratio", "new pad 2", 0.01, 0.1, 0.99, 0.38);

  TPad *titlePad = new TPad("all","all",0,0,1,1);
  titlePad->SetFillStyle(4000); 
  titlePad->Draw();

  TH1F *h_blankrap = new TH1F("h_blankrap","h_blankrap",10,-2.6,2.6);
  h_blankrap->SetTitle("");
  TH1F *h_blank2rap = new TH1F("h_blank2rap","h_blank2rap",10,-2.6,2.6);
  h_blank2rap->SetTitle("");

  canvas->cd();
  main->Draw();
  ratio->Draw();

  main->cd();
  main->SetLeftMargin(0.125);
  main->SetBottomMargin(0.11);
  main->SetTickx(1);
  main->SetTicky(1);
  float ymax = h2->GetMaximum();
  h_blankrap->SetStats(0);
  h_blankrap->SetMaximum(ymax*1.3);
  h_blankrap->GetXaxis()->SetLabelSize(0.0000000000001);
  h_blankrap->GetYaxis()->SetTitle(titleY);
  h_blankrap->GetYaxis()->SetTitleSize(0.07);
  h_blankrap->GetYaxis()->SetLabelSize(0.04);
  // rapidity
  h_blankrap->GetYaxis()->SetTitleOffset(0.9);
  h1->SetFillColor(kBlue-9);
  h2->SetLineColor(kBlack);
  h2->SetLineWidth(2);
  h_blankrap->Draw();
  h1->Draw("same hist");
  h2->Draw("same pE");
  drawText(.13, .915, "#bf{CMS} #it{Preliminary}", kBlack, 30);
  drawText(0.5, 0.915, "PbPb 1.43 nb^{#minus1} (5.02 TeV)", kBlack, 30);
  
  lumiLegend(h1, h2, x1+0.06, x2+0.06, y1+0.03, y2+0.03, label1, label2);

  ratio->cd();
  ratio->SetTopMargin(0.035);
  ratio->SetLeftMargin(0.125);
  ratio->SetBottomMargin(0.15);
  ratio->SetTickx(1);
  ratio->SetTicky(1);
  h_blank2rap->SetStats(0);
  h_ratio->SetLineColor(kBlack);
  h_blank2rap->SetTitle("");
  h_blank2rap->GetXaxis()->SetLabelSize(0.09);
  h_blank2rap->GetYaxis()->SetLabelSize(0.09);
  h_blank2rap->GetYaxis()->SetNdivisions(505);
  string ratio_label = label1;
  //ratio_label = ratio_label + " / " + label2;
  ratio_label = "Data/MC";
  float ymax_r = h_ratio->GetMaximum();
  float ymin_r = h_ratio->GetMinimum();
  // range for rap
  h_blank2rap->GetYaxis()->SetRangeUser(0,1.8);
  h_blank2rap->GetYaxis()->SetTitle(ratio_label.c_str());
  h_blank2rap->GetYaxis()->SetTitleOffset(0.38);
  h_blank2rap->GetYaxis()->SetTitleSize(0.155);
  h_ratio->SetMarkerSize(1.2);
  ratio->SetBottomMargin(0.15);
  h_blank2rap->Draw();
  h_ratio->Draw("same hist p");

  float xmax = h_blankrap->GetXaxis()->GetXmax();
  float xmin = h_blankrap->GetXaxis()->GetXmin();
  TLine *line = new TLine(xmin,1,xmax,1);
  line->SetLineColor(kBlack);
  line->SetLineWidth(3);
  line->Draw("LSAME");
  h_ratio->Draw("same hist pE");

  titlePad->cd();
  TLatex *xlabel = new TLatex();
  xlabel->SetTextSize(0.045);
  // position for rap1, rap2, dirap
  xlabel->DrawLatex(.825,.05, titleX);

  canvas->SaveAs("test.png");
}


void makeDimassPlot(TH1F *h1, TH1F *h2, const char *titleX, const char *titleY, TCanvas *canvas, string label1, string label2, float x1=0.60, float x2=0.80, float y1=0.65, float y2=0.80){
  
  TH1F *h_ratio = new TH1F(*h2);
  h_ratio->Divide(h1);
  h_ratio->SetMarkerColor(kBlack);
  h_ratio->SetMarkerStyle(kFullDotLarge);

  TPad *main = new TPad("main", "new pad", 0.01, 0.315, 0.99, 0.99);
  TPad *ratio = new TPad("ratio", "new pad 2", 0.01, 0.1, 0.99, 0.38);

  TH1F *h_blankDimass = new TH1F("h_blankDimass","h_blankDimass",10,7,25);
  h_blankDimass->SetTitle("");
  TH1F *h_blank2Dimass = new TH1F("h_blank2Dimass","h_blank2Dimass",10,7,25);
  h_blank2Dimass->SetTitle("");

  TPad *titlePad = new TPad("all","all",0,0,1,1);
  titlePad->SetFillStyle(4000); 
  titlePad->Draw();

  canvas->cd();
  main->Draw();
  ratio->Draw();

  main->cd();
  main->SetLeftMargin(0.125);
  main->SetBottomMargin(0.11);
  main->SetTickx(1);
  main->SetTicky(1);
  float ymax = h2->GetMaximum();
  h_blankDimass->SetStats(0);
  h_blankDimass->SetMaximum(ymax*1.3);
  h_blankDimass->GetXaxis()->SetLabelSize(0.0000000000001);
  h_blankDimass->GetYaxis()->SetTitle(titleY);
  h_blankDimass->GetYaxis()->SetTitleSize(0.07);
  h_blankDimass->GetYaxis()->SetLabelSize(0.04);
  // rapidity
  h_blankDimass->GetYaxis()->SetTitleOffset(0.9);
  h1->SetFillColor(kBlue-9);
  h2->SetLineColor(kBlack);
  h2->SetLineWidth(2);
  // plot overflow bin data
  /*int lastBin_1 = h1->GetNbinsX();
  double overflow_1 = h1->GetBinContent(lastBin_1 + 1);
  h1->SetBinContent(lastBin_1, h1->GetBinContent(lastBin_1) + overflow_1);
  // plot overflow bin mc
  int lastBin_2 = h2->GetNbinsX();
  double overflow_2 = h2->GetBinContent(lastBin_2 + 1);
  h2->SetBinContent(lastBin_2, h2->GetBinContent(lastBin_2) + overflow_2);*/
  h_blankDimass->Draw();
  h1->Draw("same hist");
  h2->Draw("same pE");
  drawText(.13, .915, "#bf{CMS} #it{Preliminary}", kBlack, 30);
  drawText(0.5, 0.915, "PbPb 1.43 nb^{#minus1} (5.02 TeV)", kBlack, 30);
  
  // drawLegend 
  lumiLegend(h1, h2, x1+0.06, x2+0.06, y1+0.03, y2+0.03, label1, label2);

  ratio->cd();
  ratio->SetTopMargin(0.035);
  ratio->SetLeftMargin(0.125);
  ratio->SetBottomMargin(0.15);
  ratio->SetTickx(1);
  ratio->SetTicky(1);
  h_blank2Dimass->SetStats(0);
  h_ratio->SetLineColor(kBlack);
  h_blank2Dimass->SetTitle("");
  h_blank2Dimass->GetXaxis()->SetLabelSize(0.09);
  h_blank2Dimass->GetYaxis()->SetLabelSize(0.09);
  h_blank2Dimass->GetYaxis()->SetNdivisions(505);
  string ratio_label = label1;
  //ratio_label = ratio_label + " / " + label2;
  ratio_label = "Data/MC";
  float ymax_r = h_ratio->GetMaximum();
  float ymin_r = h_ratio->GetMinimum();
  // dimass range
  h_blank2Dimass->GetYaxis()->SetRangeUser(0.7,1.3);
  h_blank2Dimass->GetYaxis()->SetTitle(ratio_label.c_str());
  h_blank2Dimass->GetYaxis()->SetTitleOffset(0.38);
  h_blank2Dimass->GetYaxis()->SetTitleSize(0.155);
  h_ratio->SetMarkerSize(1.2);
  ratio->SetBottomMargin(0.15);
  h_blank2Dimass->Draw();
  h_ratio->Draw("same hist p");

  float xmax = h_blankDimass->GetXaxis()->GetXmax();
  float xmin = h_blankDimass->GetXaxis()->GetXmin();
  TLine *line = new TLine(xmin,1,xmax,1);
  line->SetLineColor(kBlack);
  line->SetLineWidth(3);
  line->Draw("LSAME");
  h_ratio->Draw("same hist pE");

  titlePad->cd();
  TLatex *xlabel = new TLatex();
  xlabel->SetTextSize(0.045);
  // position for aco, pt1, pt2, diPt, dimass
  xlabel->DrawLatex(.7,.05, titleX);

  canvas->SaveAs("test.png");
}






void pad4(TH1F *h1, TH1F *h2, TH1F *h3, TH1F *h4, TH1F *h1b, TH1F *h2b, TH1F *h3b, TH1F *h4b, string leg1, string leg2)
{
  h1b->SetLineColor(kBlack);
  h2b->SetLineColor(kBlack);
  h3b->SetLineColor(kBlack);
  h4b->SetLineColor(kBlack);

  // Plotting
  gStyle->SetOptStat(0);
  auto *c1 = new TCanvas("c1","Title",900,900);
  TPad* p1 = new TPad("p1","p1",0.01,0.5,0.5,0.99,0); p1->Draw();
  TPad* p2 = new TPad("p2","p2",0.5,0.5,0.09,0.99,0); p2->Draw();
  TPad* p3 = new TPad("p3","p3",0.01,0.01,0.5,0.5,0); p3->Draw();
  TPad* p4 = new TPad("p4","p4",0.5,0.01,0.99,0.5,0); p4->Draw();
  p1->SetRightMargin(0);
  p2->SetLeftMargin(0);
  p3->SetRightMargin(0);
  p4->SetLeftMargin(0);

  TPad *titlePad = new TPad("all","all",0,0,1,1);
  titlePad->SetFillStyle(4000); 
  titlePad->Draw();

  p1->cd();
  float xmax = 180;
  float xmax1 = h4->GetXaxis()->GetXmax();
  float ymax1 = h1b->GetMaximum();
  h1->GetXaxis()->SetRangeUser(0,xmax);
  h1->GetYaxis()->SetRangeUser(0,1.1*ymax1);
  h1->SetTitle(";;");
  h1->Draw("hist");
  h1b->Draw("SAME hist");
  drawLegendStr(h1, h1b, 0.4, 0.9, 0.6, 0.85, leg1, leg2);
  drawText(.15, .85, "#bf{#eta: 3.15-3.5}", kBlack, 22);

  p2->cd();
  h2->GetXaxis()->SetRangeUser(0,xmax);
  h2->GetYaxis()->SetRangeUser(0,1.1*ymax1);
  h2->SetTitle(";;");
  h2->Draw("hist");
  h2b->Draw("SAME hist");
  drawText(.05, .85, "#bf{#eta: 3.5-4}", kBlack, 22);

  p3->cd();
  h3->GetXaxis()->SetRangeUser(0,xmax);
  h3->GetYaxis()->SetRangeUser(0,1.1*ymax1);
  h3->SetTitle(";;");
  h3->Draw("hist");
  h3b->Draw("SAME hist");
  drawText(.15, .85, "#bf{#eta: 4-4.5}", kBlack, 22);

  p4->cd();
  h4->GetXaxis()->SetRangeUser(0,xmax);
  h4->GetYaxis()->SetRangeUser(0,1.1*ymax1);
  h4->SetTitle(";;");
  h4->Draw("hist");
  h4b->Draw("SAME hist");
  drawText(.05, .85, "#bf{#eta: 4.5-5.2}", kBlack, 22);

  titlePad->cd();
  TLatex *title = new TLatex();
  title->SetTextSize(0.0375);
  title->DrawLatex(.28,.97,"HF+ Avg Energy Tower");
  TLatex *xlabel = new TLatex();
  xlabel->SetTextSize(0.03);
  xlabel->DrawLatex(.4,.005, "Energy [GeV]");
  TLatex *ylabel = new TLatex();
  ylabel->SetTextSize(0.03);
  ylabel->SetTextAngle(90);
  ylabel->DrawLatex(.0125,.475, "Normalized Counts");

  c1->SaveAs("hfp_partialSums.png");
}



void pad5(TH1F *h1, TH1F *h2, TH1F *h3, TH1F *h4, TH1F *h5)
{
  // Plotting
  gStyle->SetOptStat(0);
  auto *c1 = new TCanvas("c1","Title",1600,800);
  TPad* p1 = new TPad("p1","p1",0.01,0.5,0.33,0.99,0); p1->Draw();
  TPad* p2 = new TPad("p2","p2",0.33,0.5,0.66,0.99,0); p2->Draw();
  TPad* p3 = new TPad("p3","p3",0.66,0.5,0.99,0.99,0); p3->Draw();
  TPad* p4 = new TPad("p4","p4",0.17,0.01,0.5,0.5,0); p4->Draw();
  TPad* p5 = new TPad("p5","p5",0.5,0.01,0.83,0.5,0); p5->Draw();
  p1->SetRightMargin(0);
  p2->SetLeftMargin(0);
  p2->SetRightMargin(0);
  p3->SetLeftMargin(0);
  p4->SetTopMargin(0);
  p4->SetRightMargin(0);
  p5->SetTopMargin(0);
  p5->SetLeftMargin(0);

  TPad *titlePad = new TPad("all","all",0,0,1,1);
  titlePad->SetFillStyle(4000); 
  titlePad->Draw();

  p1->cd();
  float xmin = -3.2;
  float xmax = 3.2;
  float xmax1 = h1->GetXaxis()->GetXmax();
  float ymax1 = h2->GetMaximum();
  h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->GetYaxis()->SetRangeUser(0,1.1*ymax1);
  h1->SetTitle(";;");
  h1->Draw("hist");
  drawText(.15, .8, "#bf{0 GeV < Q_{T} < 5 GeV}", kBlack, 22);

  p2->cd();
  h2->GetXaxis()->SetRangeUser(xmin,xmax);
  h2->GetYaxis()->SetRangeUser(0,1.1*ymax1);
  h2->SetTitle(";;");
  h2->Draw("hist");
  drawText(.05, .8, "#bf{5 GeV < Q_{T} < 10 GeV}", kBlack, 22);

  p3->cd();
  h3->GetXaxis()->SetRangeUser(xmin,xmax);
  h3->GetYaxis()->SetRangeUser(0,1.1*ymax1);
  h3->SetTitle(";;");
  h3->Draw("hist");
  drawText(.05, .8, "#bf{10 GeV < Q_{T} < 15 GeV}", kBlack, 22);

  p4->cd();
  h4->GetXaxis()->SetRangeUser(xmin,xmax);
  h4->GetYaxis()->SetRangeUser(0,1.1*ymax1);
  h4->SetTitle(";;");
  h4->Draw("hist");
  drawText(.15, .9, "#bf{15 GeV < Q_{T} < 20 GeV}", kBlack, 22);

  p5->cd();
  h5->GetXaxis()->SetRangeUser(xmin,xmax);
  h5->GetYaxis()->SetRangeUser(0,1.1*ymax1);
  h5->SetTitle(";;");
  h5->Draw("hist");
  drawText(.05, .9, "#bf{20 GeV < Q_{T} < 25 GeV}", kBlack, 22);

  titlePad->cd();
  TLatex *title = new TLatex();
  title->SetTextSize(0.0375);
  title->DrawLatex(.4,.97,"Q_{T} #phi (HF Cuts)");
  TLatex *xlabel = new TLatex();
  xlabel->SetTextSize(0.03);
  xlabel->DrawLatex(.45,.005, "#phi [radians]");
  TLatex *ylabel = new TLatex();
  ylabel->SetTextSize(0.03);
  ylabel->SetTextAngle(90);
  ylabel->DrawLatex(.0125,.475, "Counts");

  c1->SaveAs("qt_angles_HF.png");
}





void pad5_5(TH1F *h1, TH1F *h2, TH1F *h3, TH1F *h4, TH1F *h5, TH1F *h1b, TH1F *h2b, TH1F *h3b, TH1F *h4b, TH1F *h5b, string leg1, string leg2)
{
  h1b->SetLineColor(kRed);
  h2b->SetLineColor(kRed);
  h3b->SetLineColor(kRed);
  h4b->SetLineColor(kRed);
  h5b->SetLineColor(kRed);

  // Plotting
  gStyle->SetOptStat(0);
  auto *c1 = new TCanvas("c1","Title",1600,800);
  TPad* p1 = new TPad("p1","p1",0.01,0.5,0.33,0.99,0); p1->Draw();
  TPad* p2 = new TPad("p2","p2",0.33,0.5,0.66,0.99,0); p2->Draw();
  TPad* p3 = new TPad("p3","p3",0.66,0.5,0.99,0.99,0); p3->Draw();
  TPad* p4 = new TPad("p4","p4",0.17,0.01,0.5,0.5,0); p4->Draw();
  TPad* p5 = new TPad("p5","p5",0.5,0.01,0.83,0.5,0); p5->Draw();
  p1->SetRightMargin(0);
  p2->SetLeftMargin(0);
  p2->SetRightMargin(0);
  p3->SetLeftMargin(0);
  p4->SetTopMargin(0);
  p4->SetRightMargin(0);
  p5->SetTopMargin(0);
  p5->SetLeftMargin(0);

  TPad *titlePad = new TPad("all","all",0,0,1,1);
  titlePad->SetFillStyle(4000); 
  titlePad->Draw();

  p1->cd();
  float xmax1 = h1->GetXaxis()->GetXmax();
  float ymax1 = h5->GetMaximum();
  h1->GetXaxis()->SetRangeUser(0,200);
  h1->GetYaxis()->SetRangeUser(0,1.1*ymax1);
  h1->SetTitle(";;");
  h1->Draw("hist");
  h1b->Draw("SAME hist");
  drawLegendStr(h1, h1b, 0.7, 0.9, 0.5, 0.7, leg1, leg2);
  drawText(.15, .8, "#bf{EM1}", kBlack, 22);

  p2->cd();
  h2->GetXaxis()->SetRangeUser(0,200);
  h2->GetYaxis()->SetRangeUser(0,1.1*ymax1);
  h2->SetTitle(";;");
  h2->Draw("hist");
  h2b->Draw("SAME hist");
  drawText(.05, .8, "#bf{EM2}", kBlack, 22);

  p3->cd();
  h3->GetXaxis()->SetRangeUser(0,200);
  h3->GetYaxis()->SetRangeUser(0,1.1*ymax1);
  h3->SetTitle(";;");
  h3->Draw("hist");
  h3b->Draw("SAME hist");
  drawText(.05, .8, "#bf{EM3}", kBlack, 22);

  p4->cd();
  h4->GetXaxis()->SetRangeUser(0,200);
  h4->GetYaxis()->SetRangeUser(0,1.1*ymax1);
  h4->SetTitle(";;");
  h4->Draw("hist");
  h4b->Draw("SAME hist");
  drawText(.15, .9, "#bf{EM4}", kBlack, 22);

  p5->cd();
  h5->GetXaxis()->SetRangeUser(0,200);
  h5->GetYaxis()->SetRangeUser(0,1.1*ymax1);
  h5->SetTitle(";;");
  h5->Draw("hist");
  h5b->Draw("SAME hist");
  drawText(.05, .9, "#bf{EM5}", kBlack, 22);

  titlePad->cd();
  TLatex *title = new TLatex();
  title->SetTextSize(0.0375);
  title->DrawLatex(.4,.97,"Max Charge EM Channels 1-5");
  TLatex *xlabel = new TLatex();
  xlabel->SetTextSize(0.03);
  xlabel->DrawLatex(.45,.005, "Max Charge (fC)");
  TLatex *ylabel = new TLatex();
  ylabel->SetTextSize(0.03);
  ylabel->SetTextAngle(90);
  ylabel->DrawLatex(.0125,.475, "Counts");

  c1->SaveAs("allEM_overlay1M.png");
}


void pad5_4(TH1F *h1, TH1F *h2, TH1F *h3, TH1F *h4, TH1F *h5, TH1F *h1b, TH1F *h2b, TH1F *h3b, TH1F *h4b, TH1F *h5b, TH1F *h1c, TH1F *h2c, TH1F *h3c, TH1F *h4c, TH1F *h5c, TH1F *h1d, TH1F *h2d, TH1F *h3d, TH1F *h4d, TH1F *h5d, string leg1, string leg2, string leg3, string leg4)
{
  string str_0 = leg1;
  string str_1000 = leg2;
  string str_null = leg3;
  string str_null2 = leg4;

  h1b->SetLineColor(kRed);
  h2b->SetLineColor(kRed);
  h3b->SetLineColor(kRed);
  h4b->SetLineColor(kRed);
  h5b->SetLineColor(kRed);
  
  h1->SetLineColor(kGreen+2);
  h2->SetLineColor(kGreen+2);
  h3->SetLineColor(kGreen+2);
  h4->SetLineColor(kGreen+2);
  h5->SetLineColor(kGreen+2);

  h1d->SetLineColor(kViolet);
  h2d->SetLineColor(kViolet);
  h3d->SetLineColor(kViolet);
  h4d->SetLineColor(kViolet);
  h5d->SetLineColor(kViolet);

 

  // Plotting
  gStyle->SetOptStat(0);
  auto *c1 = new TCanvas("c1","Title",1600,800);
  TPad* p1 = new TPad("p1","p1",0.01,0.5,0.33,0.99,0); p1->Draw();
  TPad* p2 = new TPad("p2","p2",0.33,0.5,0.66,0.99,0); p2->Draw();
  TPad* p3 = new TPad("p3","p3",0.66,0.5,0.99,0.99,0); p3->Draw();
  TPad* p4 = new TPad("p4","p4",0.17,0.01,0.5,0.5,0); p4->Draw();
  TPad* p5 = new TPad("p5","p5",0.5,0.01,0.83,0.5,0); p5->Draw();
  p1->SetRightMargin(0);
  p2->SetLeftMargin(0);
  p2->SetRightMargin(0);
  p3->SetLeftMargin(0);
  p4->SetTopMargin(0);
  p4->SetRightMargin(0);
  p5->SetTopMargin(0);
  p5->SetLeftMargin(0);

  TPad *titlePad = new TPad("all","all",0,0,1,1);
  titlePad->SetFillStyle(4000); 
  titlePad->Draw();

  p1->cd();
  p1->SetLogy();
  float xmax1 = h1->GetXaxis()->GetXmax();
  float xmax2 = 500;
  float ymax1 = h1d->GetMaximum();
  h1->GetXaxis()->SetRangeUser(0,xmax2);
  //h1->GetYaxis()->SetRangeUser(0,1.1*ymax1);
  h1->GetYaxis()->SetRangeUser(1,10.1*ymax1);
  h1->SetTitle(";;");
  h1->Draw("hist");
  h1b->Draw("SAME hist");
  h1c->Draw("SAME hist");
  h1d->Draw("SAME hist");
  //drawLegendStr4(h1c, h1d, h1, h1b, 0.55, 0.85, 0.65, 0.85, str_null, str_null2, str_0, str_1000);
  drawText(.25, .8, "#bf{EM1}", kBlack, 22);

  p2->cd();
  p2->SetLogy();
  h2->GetXaxis()->SetRangeUser(0,xmax2);
  //h2->GetYaxis()->SetRangeUser(0,1.1*ymax1);
  h2->GetYaxis()->SetRangeUser(1,10.1*ymax1);
  h2->SetTitle(";;");
  h2->Draw("hist");
  h2b->Draw("SAME hist");
  h2c->Draw("SAME hist");
  h2d->Draw("SAME hist");
  drawText(.15, .8, "#bf{EM2}", kBlack, 22);

  p3->cd();
  p3->SetLogy();
  h3->GetXaxis()->SetRangeUser(0,xmax2);
  //h3->GetYaxis()->SetRangeUser(0,1.1*ymax1);
  h3->GetYaxis()->SetRangeUser(1,10.1*ymax1);
  h3->SetTitle(";;");
  h3->Draw("hist");
  h3b->Draw("SAME hist");
  h3c->Draw("SAME hist");
  h3d->Draw("SAME hist");
  drawText(.15, .8, "#bf{EM3}", kBlack, 22);

  p4->cd();
  p4->SetLogy();
  h4->GetXaxis()->SetRangeUser(0,xmax2);
  //h4->GetYaxis()->SetRangeUser(0,1.1*ymax1);
  h4->GetYaxis()->SetRangeUser(1,10.1*ymax1);
  h4->SetTitle(";;");
  h4->Draw("hist");
  h4b->Draw("SAME hist");
  h4c->Draw("SAME hist");
  h4d->Draw("SAME hist");
  drawText(.25, .9, "#bf{EM4}", kBlack, 22);

  p5->cd();
  p5->SetLogy();
  h5->GetXaxis()->SetRangeUser(0,xmax2);
  //h5->GetYaxis()->SetRangeUser(0,1.1*ymax1);
  h5->GetYaxis()->SetRangeUser(1,10.1*ymax1);
  h5->SetTitle(";;");
  h5->Draw("hist");
  h5b->Draw("SAME hist");
  h5c->Draw("SAME hist");
  h5d->Draw("SAME hist");
  drawText(.15, .9, "#bf{EM5}", kBlack, 22);

  titlePad->cd();
  TLatex *title = new TLatex();
  title->SetTextSize(0.0375);
  title->DrawLatex(.4,.97,"Max Charge EM Channels 1-5");
  TLatex *xlabel = new TLatex();
  xlabel->SetTextSize(0.03);
  xlabel->DrawLatex(.42,.005, "Max Charge (fC)");
  TLatex *ylabel = new TLatex();
  ylabel->SetTextSize(0.03);
  ylabel->SetTextAngle(90);
  ylabel->DrawLatex(.0125,.475, "Counts");

  c1->SaveAs("masterCompEM_logYext.png");
  
}


void pad4_4(TH1F *h1, TH1F *h2, TH1F *h3, TH1F *h4, TH1F *h1b, TH1F *h2b, TH1F *h3b, TH1F *h4b, TH1F *h1c, TH1F *h2c, TH1F *h3c, TH1F *h4c, TH1F *h1d, TH1F *h2d, TH1F *h3d, TH1F *h4d, string leg1, string leg2, string leg3, string leg4)
{
  string str_0 = leg1;
  string str_1000 = leg2;
  string str_null = leg3;
  string str_null2 = leg4;

  h1b->SetLineColor(kRed);
  h2b->SetLineColor(kRed);
  h3b->SetLineColor(kRed);
  h4b->SetLineColor(kRed);
  
  h1->SetLineColor(kGreen+2);
  h2->SetLineColor(kGreen+2);
  h3->SetLineColor(kGreen+2);
  h4->SetLineColor(kGreen+2);

  h1d->SetLineColor(kViolet);
  h2d->SetLineColor(kViolet);
  h3d->SetLineColor(kViolet);
  h4d->SetLineColor(kViolet);

  // Plotting
  gStyle->SetOptStat(0);
  auto *c1 = new TCanvas("c1","Title",900,900);
  TPad* p1 = new TPad("p1","p1",0.01,0.5,0.5,0.99,0); p1->Draw();
  TPad* p2 = new TPad("p2","p2",0.5,0.5,0.99,0.99,0); p2->Draw();
  TPad* p3 = new TPad("p3","p3",0.01,0.01,0.5,0.5,0); p3->Draw();
  TPad* p4 = new TPad("p4","p4",0.5,0.01,0.99,0.5,0); p4->Draw();
  p1->SetRightMargin(0);
  p2->SetLeftMargin(0);
  p3->SetRightMargin(0);
  p4->SetLeftMargin(0);
  
  TPad *titlePad = new TPad("all","all",0,0,1,1);
  titlePad->SetFillStyle(4000); 
  titlePad->Draw();

  p1->cd();
  p1->SetLogy();
  float xmax1 = h1->GetXaxis()->GetXmax();
  float xmax2 = 450;
  float ymax1 = h2c->GetMaximum();
  h1->GetXaxis()->SetRangeUser(0,xmax2);
  //h1->GetYaxis()->SetRangeUser(0,1.1*ymax1);
  h1->GetYaxis()->SetRangeUser(1,10.1*ymax1);
  h1->SetTitle(";;");
  h1->Draw("hist");
  h1b->Draw("SAME hist");
  h1c->Draw("SAME hist");
  h1d->Draw("SAME hist");
  //drawLegendStr4(h1c, h1d, h1, h1b, 0.55, 0.9, 0.6, 0.85, str_null, str_null2, str_0, str_1000);
  drawText(.15, .8, "#bf{HAD1}", kBlack, 22);

  p2->cd();
  p2->SetLogy();
  h2->GetXaxis()->SetRangeUser(0,xmax2);
  //h2->GetYaxis()->SetRangeUser(0,1.1*ymax1);
  h2->GetYaxis()->SetRangeUser(1,10.1*ymax1);
  h2->SetTitle(";;");
  h2->Draw("hist");
  h2b->Draw("SAME hist");
  h2c->Draw("SAME hist");
  h2d->Draw("SAME hist");
  drawText(.05, .8, "#bf{HAD2}", kBlack, 22);

  p3->cd();
  p3->SetLogy();
  h3->GetXaxis()->SetRangeUser(0,xmax2);
  //h3->GetYaxis()->SetRangeUser(0,1.1*ymax1);
  h3->GetYaxis()->SetRangeUser(1,10.1*ymax1);
  h3->SetTitle(";;");
  h3->Draw("hist");
  h3b->Draw("SAME hist");
  h3c->Draw("SAME hist");
  h3d->Draw("SAME hist");
  drawText(.15, .8, "#bf{HAD3}", kBlack, 22);

  p4->cd();
  p4->SetLogy();
  h4->GetXaxis()->SetRangeUser(0,xmax2);
  //h4->GetYaxis()->SetRangeUser(0,1.1*ymax1);
  h4->GetYaxis()->SetRangeUser(1,10.1*ymax1);
  h4->SetTitle(";;");
  h4->Draw("hist");
  h4b->Draw("SAME hist");
  h4c->Draw("SAME hist");
  h4d->Draw("SAME hist");
  drawText(.05, .8, "#bf{HAD4}", kBlack, 22);

  titlePad->cd();
  TLatex *title = new TLatex();
  title->SetTextSize(0.0375);
  title->DrawLatex(.25,.97,"Max Charge HAD Channels 1-5");
  TLatex *xlabel = new TLatex();
  xlabel->SetTextSize(0.03);
  xlabel->DrawLatex(.42,.005, "Max Charge (fC)");
  TLatex *ylabel = new TLatex();
  ylabel->SetTextSize(0.03);
  ylabel->SetTextAngle(90);
  ylabel->DrawLatex(.0125,.475, "Counts");

  c1->SaveAs("masterCompHad_logYext.png");
  
}




void pad6(TH1F *h1, TH1F *h2, TH1F *h3, TH1F *h4, TH1F *h5, TH1F *h6)
{

  // Plotting
  gStyle->SetOptStat(0);
  auto *c1 = new TCanvas("c1","Title",1600,800);
  TPad* p1 = new TPad("p1","p1",0.01,0.5,0.33,0.99,0); p1->Draw();
  TPad* p2 = new TPad("p2","p2",0.33,0.5,0.66,0.99,0); p2->Draw();
  TPad* p3 = new TPad("p3","p3",0.66,0.5,0.99,0.99,0); p3->Draw();
  TPad* p4 = new TPad("p4","p4",0.01,0.01,0.33,0.5,0); p4->Draw();
  TPad* p5 = new TPad("p5","p5",0.33,0.01,0.66,0.5,0); p5->Draw();
  TPad* p6 = new TPad("p6","p6",0.66,0.01,0.99,0.5,0); p6->Draw();
  p1->SetBottomMargin(0);
  p1->SetRightMargin(0);
  p2->SetBottomMargin(0);
  p2->SetLeftMargin(0);
  p2->SetRightMargin(0);
  p3->SetBottomMargin(0);
  p3->SetLeftMargin(0);
  p4->SetTopMargin(0);
  p4->SetRightMargin(0);
  p5->SetTopMargin(0);
  p5->SetLeftMargin(0);
  p5->SetRightMargin(0);
  p6->SetTopMargin(0);
  p6->SetLeftMargin(0);

  TPad *titlePad = new TPad("all","all",0,0,1,1);
  titlePad->SetFillStyle(4000); 
  titlePad->Draw();

  p1->cd();
  //p1->SetLogx();
  //p1->SetLogy();
  float ymin = 1;
  float xmax1 = h1->GetXaxis()->GetXmax();
  float ymax1 = 0.6*(h4->GetMaximum());
  //float ymax1 = 0.025;
  float xmax = 10500;
  h1->GetXaxis()->SetRangeUser(0,xmax);
  h1->GetYaxis()->SetRangeUser(ymin,ymax1);
  h1->SetTitle(";;");
  h1->Draw("hist");
  //drawText(.85, .8, "#bf{TS 1}", kBlack, 22);
  drawText(.5, .85, "#bf{evt 1, evt 2: p_{T} cuts}", kBlack, 22);
  drawText(.5, .775, "#bf{evt 1, evt 2: P_{T} > Q_{T}}", kBlack, 22);

  p2->cd();
  //p2->SetLogx();
  //p2->SetLogy();
  h2->GetXaxis()->SetRangeUser(0,xmax);
  h2->GetYaxis()->SetRangeUser(ymin,ymax1);
  h2->SetTitle(";;");
  h2->Draw("hist");
  //drawText(.85, .8, "#bf{TS 2}", kBlack, 22);
  drawText(.5, .85, "#bf{evt 1, evt 2: p_{T} cuts}", kBlack, 22);
  drawText(.5, .775, "#bf{evt 1, evt 2: P_{T} > Q_{T}}", kBlack, 22);
  drawText(.5, .7, "#bf{Mixed evt: P_{T} > Q_{T}}", kBlack, 22);

  p3->cd();
  //p3->SetLogx();
  //p3->SetLogy();
  h3->GetXaxis()->SetRangeUser(0,xmax);
  h3->GetYaxis()->SetRangeUser(ymin,ymax1);
  h3->SetTitle(";;");
  h3->Draw("hist");
  //drawText(.75, .8, "#bf{TS 3}", kBlack, 22);
  drawText(.5, .85, "#bf{evt 1, evt 2: p_{T} cuts}", kBlack, 22);
  drawText(.5, .775, "#bf{Mixed evt: P_{T} > Q_{T}}", kBlack, 22);

  p4->cd();
  //p4->SetLogx();
  //p4->SetLogy();
  h4->GetXaxis()->SetRangeUser(0,xmax);
  h4->GetYaxis()->SetRangeUser(ymin,ymax1);
  h4->SetTitle(";;");
  h4->Draw("hist");
  //drawText(.85, .9, "#bf{TS 4}", kBlack, 22);
  drawText(.5, .875, "#bf{Mixed evt: p_{T} cuts}", kBlack, 22);
  drawText(.5, .8, "#bf{evt 1, evt 2: P_{T} > Q_{T}}", kBlack, 22);

  p5->cd();
  //p5->SetLogx();
  //p5->SetLogy();
  h5->GetXaxis()->SetRangeUser(0,xmax);
  h5->GetYaxis()->SetRangeUser(ymin,ymax1);
  h5->SetTitle(";;");
  h5->Draw("hist");
  //drawText(.85, .9, "#bf{TS 5}", kBlack, 22);
  drawText(.5, .875, "#bf{Mixed evt: p_{T} cuts}", kBlack, 22);
  drawText(.5, .8, "#bf{evt 1, evt 2: P_{T} > Q_{T}}", kBlack, 22);
  drawText(.5, .725, "#bf{Mixed evt: P_{T} > Q_{T}}", kBlack, 22);

  p6->cd();
  //p6->SetLogx();
  //p6->SetLogy();
  h6->GetXaxis()->SetRangeUser(0,xmax);
  h6->GetYaxis()->SetRangeUser(ymin,ymax1);
  h6->SetTitle(";;");
  h6->Draw("hist");
  //drawText(.75, .9, "#bf{TS 6}", kBlack, 22);
  drawText(.5, .875, "#bf{Mixed evt: p_{T} cuts}", kBlack, 22);
  drawText(.5, .8, "#bf{Mixed evt: P_{T} > Q_{T}}", kBlack, 22);

  titlePad->cd();
  TLatex *title = new TLatex();
  title->SetTextSize(0.0375);
  title->DrawLatex(.5,.97,"#Delta#phi");
  TLatex *xlabel = new TLatex();
  xlabel->SetTextSize(0.03);
  xlabel->DrawLatex(.45,.005, "#Delta#phi [radians]");
  TLatex *ylabel = new TLatex();
  ylabel->SetTextSize(0.03);
  ylabel->SetTextAngle(90);
  ylabel->DrawLatex(.0125,.475, "Normalized Counts");

  c1->SaveAs("oct26ped_allTS_fib17_ch2.png");
}




void pad6_6(TH1F *h1, TH1F *h2, TH1F *h3, TH1F *h4, TH1F *h5, TH1F *h6, TH1F *h1_b, TH1F *h2_b, TH1F *h3_b, TH1F *h4_b, TH1F *h5_b, TH1F *h6_b)
{

  // Plotting
  gStyle->SetOptStat(0);
  auto *c1 = new TCanvas("c1","Title",1600,800);
  TPad* p1 = new TPad("p1","p1",0.01,0.5,0.33,0.99,0); p1->Draw();
  TPad* p2 = new TPad("p2","p2",0.33,0.5,0.66,0.99,0); p2->Draw();
  TPad* p3 = new TPad("p3","p3",0.66,0.5,0.99,0.99,0); p3->Draw();
  TPad* p4 = new TPad("p4","p4",0.01,0.01,0.33,0.5,0); p4->Draw();
  TPad* p5 = new TPad("p5","p5",0.33,0.01,0.66,0.5,0); p5->Draw();
  TPad* p6 = new TPad("p6","p6",0.66,0.01,0.99,0.5,0); p6->Draw();
  p1->SetBottomMargin(0);
  p1->SetRightMargin(0);
  p2->SetBottomMargin(0);
  p2->SetLeftMargin(0);
  p2->SetRightMargin(0);
  p3->SetBottomMargin(0);
  p3->SetLeftMargin(0);
  p4->SetTopMargin(0);
  p4->SetRightMargin(0);
  p5->SetTopMargin(0);
  p5->SetLeftMargin(0);
  p5->SetRightMargin(0);
  p6->SetTopMargin(0);
  p6->SetLeftMargin(0);

  TPad *titlePad = new TPad("all","all",0,0,1,1);
  titlePad->SetFillStyle(4000); 
  titlePad->Draw();

  p1->cd();
  //p1->SetLogx();
  //p1->SetLogy();
  float ymin = 1;
  float xmax1 = h1->GetXaxis()->GetXmax();
  float ymax1 = h4->GetMaximum();
  float xmax = 10500;
  h1->GetXaxis()->SetRangeUser(0,xmax);
  h1->GetYaxis()->SetRangeUser(ymin,ymax1);
  h1->SetTitle(";;");
  h1->Draw("hist");
  h1_b->Draw("hist same");
  //drawText(.85, .8, "#bf{TS 1}", kBlack, 22);
  drawLegendStr(h1, h1_b, 0.5, 0.7, 0.6, 0.7, "leading jet", "subleading jet");
  drawText(.5, .85, "#bf{evt 1, evt 2: p_{T} cuts}", kBlack, 22);
  drawText(.5, .775, "#bf{evt 1, evt 2: P_{T} > Q_{T}}", kBlack, 22);

  p2->cd();
  //p2->SetLogx();
  //p2->SetLogy();
  h2->GetXaxis()->SetRangeUser(0,xmax);
  h2->GetYaxis()->SetRangeUser(ymin,ymax1);
  h2->SetTitle(";;");
  h2->Draw("hist");
  h2_b->Draw("hist same");
  //drawText(.85, .8, "#bf{TS 2}", kBlack, 22);
  drawText(.5, .85, "#bf{evt 1, evt 2: p_{T} cuts}", kBlack, 22);
  drawText(.5, .775, "#bf{evt 1, evt 2: P_{T} > Q_{T}}", kBlack, 22);
  drawText(.5, .7, "#bf{Mixed evt: P_{T} > Q_{T}}", kBlack, 22);

  p3->cd();
  //p3->SetLogx();
  //p3->SetLogy();
  h3->GetXaxis()->SetRangeUser(0,xmax);
  h3->GetYaxis()->SetRangeUser(ymin,ymax1);
  h3->SetTitle(";;");
  h3->Draw("hist");
  h3_b->Draw("hist same");
  //drawText(.75, .8, "#bf{TS 3}", kBlack, 22);
  drawText(.5, .85, "#bf{evt 1, evt 2: p_{T} cuts}", kBlack, 22);
  drawText(.5, .775, "#bf{Mixed evt: P_{T} > Q_{T}}", kBlack, 22);

  p4->cd();
  //p4->SetLogx();
  //p4->SetLogy();
  h4->GetXaxis()->SetRangeUser(0,xmax);
  h4->GetYaxis()->SetRangeUser(ymin,ymax1);
  h4->SetTitle(";;");
  h4->Draw("hist");
  h4_b->Draw("hist same");
  //drawText(.85, .9, "#bf{TS 4}", kBlack, 22);
  drawText(.5, .875, "#bf{Mixed evt: p_{T} cuts}", kBlack, 22);
  drawText(.5, .8, "#bf{evt 1, evt 2: P_{T} > Q_{T}}", kBlack, 22);

  p5->cd();
  //p5->SetLogx();
  //p5->SetLogy();
  h5->GetXaxis()->SetRangeUser(0,xmax);
  h5->GetYaxis()->SetRangeUser(ymin,ymax1);
  h5->SetTitle(";;");
  h5->Draw("hist");
  h5_b->Draw("hist same");
  //drawText(.85, .9, "#bf{TS 5}", kBlack, 22);
  drawText(.5, .875, "#bf{Mixed evt: p_{T} cuts}", kBlack, 22);
  drawText(.5, .8, "#bf{evt 1, evt 2: P_{T} > Q_{T}}", kBlack, 22);
  drawText(.5, .725, "#bf{Mixed evt: P_{T} > Q_{T}}", kBlack, 22);

  p6->cd();
  //p6->SetLogx();
  //p6->SetLogy();
  h6->GetXaxis()->SetRangeUser(0,xmax);
  h6->GetYaxis()->SetRangeUser(ymin,ymax1);
  h6->SetTitle(";;");
  h6->Draw("hist");
  h6_b->Draw("hist same");
  //drawText(.75, .9, "#bf{TS 6}", kBlack, 22);
  drawText(.5, .875, "#bf{Mixed evt: p_{T} cuts}", kBlack, 22);
  drawText(.5, .8, "#bf{Mixed evt: P_{T} > Q_{T}}", kBlack, 22);

  titlePad->cd();
  TLatex *title = new TLatex();
  title->SetTextSize(0.0375);
  title->DrawLatex(.35,.97,"");
  TLatex *xlabel = new TLatex();
  xlabel->SetTextSize(0.03);
  xlabel->DrawLatex(.45,.005, "p_{T} [GeV]");
  TLatex *ylabel = new TLatex();
  ylabel->SetTextSize(0.03);
  ylabel->SetTextAngle(90);
  ylabel->DrawLatex(.0125,.475, "Normalized Counts");

  c1->SaveAs("oct26ped_allTS_fib17_ch2.png");
}




void pad6_2d(TH2F *h1, TH2F *h2, TH2F *h3, TH2F *h4, TH2F *h5, TH2F *h6)
{

  // Plotting
  gStyle->SetOptStat(0);
  auto *c1 = new TCanvas("c1","Title",1600,800);
  TPad* p1 = new TPad("p1","p1",0.01,0.5,0.33,0.99,0); p1->Draw();
  TPad* p2 = new TPad("p2","p2",0.33,0.5,0.66,0.99,0); p2->Draw();
  TPad* p3 = new TPad("p3","p3",0.66,0.5,0.99,0.99,0); p3->Draw();
  TPad* p4 = new TPad("p4","p4",0.01,0.01,0.33,0.5,0); p4->Draw();
  TPad* p5 = new TPad("p5","p5",0.33,0.01,0.66,0.5,0); p5->Draw();
  TPad* p6 = new TPad("p6","p6",0.66,0.01,0.99,0.5,0); p6->Draw();
  p1->SetBottomMargin(0);
  p1->SetRightMargin(0);
  p2->SetBottomMargin(0);
  p2->SetLeftMargin(0);
  p2->SetRightMargin(0);
  p3->SetBottomMargin(0);
  p3->SetLeftMargin(0);
  p4->SetTopMargin(0);
  p4->SetRightMargin(0);
  p5->SetTopMargin(0);
  p5->SetLeftMargin(0);
  p5->SetRightMargin(0);
  p6->SetTopMargin(0);
  p6->SetLeftMargin(0);

  TPad *titlePad = new TPad("all","all",0,0,1,1);
  titlePad->SetFillStyle(4000); 
  titlePad->Draw();

  p1->cd();
  //p1->SetLogz();
  float ymin = 1;
  float xmax1 = h1->GetXaxis()->GetXmax();
  float ymax1 = h4->GetMaximum();
  float xmax = 10500;
  //h1->GetXaxis()->SetRangeUser(0,xmax);
  h1->GetYaxis()->SetRangeUser(ymin,1.1*ymax1);
  h1->SetTitle(";;");
  h1->Draw("colz");
  drawText(.5, .85, "#bf{evt 1, evt 2: p_{T} cuts}", kBlack, 22);
  drawText(.5, .775, "#bf{evt 1, evt 2: P_{T} > Q_{T}}", kBlack, 22);
  

  p2->cd();
  //h2->GetXaxis()->SetRangeUser(0,xmax);
  h2->GetYaxis()->SetRangeUser(ymin,1.1*ymax1);
  h2->SetTitle(";;");
  h2->Draw("colz");
  //drawText(.85, .8, "#bf{TS 2}", kBlack, 22);
  drawText(.5, .85, "#bf{evt 1, evt 2: p_{T} cuts}", kBlack, 22);
  drawText(.5, .775, "#bf{evt 1, evt 2: P_{T} > Q_{T}}", kBlack, 22);
  drawText(.5, .7, "#bf{Mixed evt: P_{T} > Q_{T}}", kBlack, 22);

  p3->cd();
  //h3->GetXaxis()->SetRangeUser(0,xmax);
  h3->GetYaxis()->SetRangeUser(ymin,1.1*ymax1);
  h3->SetTitle(";;");
  h3->Draw("colz");
  //drawText(.75, .8, "#bf{TS 3}", kBlack, 22);
  drawText(.5, .85, "#bf{evt 1, evt 2: p_{T} cuts}", kBlack, 22);
  drawText(.5, .775, "#bf{Mixed evt: P_{T} > Q_{T}}", kBlack, 22);

  p4->cd();
  //h4->GetXaxis()->SetRangeUser(0,xmax);
  h4->GetYaxis()->SetRangeUser(ymin,1.1*ymax1);
  h4->SetTitle(";;");
  h4->Draw("colz");
  //drawText(.85, .9, "#bf{TS 4}", kBlack, 22);
  drawText(.5, .875, "#bf{Mixed evt: p_{T} cuts}", kBlack, 22);
  drawText(.5, .8, "#bf{evt 1, evt 2: P_{T} > Q_{T}}", kBlack, 22);

  p5->cd();
  //h5->GetXaxis()->SetRangeUser(0,xmax);
  h5->GetYaxis()->SetRangeUser(ymin,1.1*ymax1);
  h5->SetTitle(";;");
  h5->Draw("colz");
  //drawText(.85, .9, "#bf{TS 5}", kBlack, 22);
  drawText(.5, .85, "#bf{evt 1, evt 2: p_{T} cuts}", kBlack, 22);
  drawText(.5, .775, "#bf{evt 1, evt 2: P_{T} > Q_{T}}", kBlack, 22);
  drawText(.5, .7, "#bf{Mixed evt: P_{T} > Q_{T}}", kBlack, 22);


  p6->cd();
  //h6->GetXaxis()->SetRangeUser(0,xmax);
  h6->GetYaxis()->SetRangeUser(ymin,1.1*ymax1);
  h6->SetTitle(";;");
  h6->Draw("colz");
  //drawText(.75, .9, "#bf{TS 6}", kBlack, 22);
  drawText(.5, .875, "#bf{Mixed evt: p_{T} cuts}", kBlack, 22);
  drawText(.5, .8, "#bf{Mixed evt: P_{T} > Q_{T}}", kBlack, 22);

  titlePad->cd();
  TLatex *title = new TLatex();
  title->SetTextSize(0.0375);
  title->DrawLatex(.35,.97,"");
  TLatex *xlabel = new TLatex();
  xlabel->SetTextSize(0.03);
  xlabel->DrawLatex(.45,.005, "Q_{T} [GeV]");
  TLatex *ylabel = new TLatex();
  ylabel->SetTextSize(0.03);
  ylabel->SetTextAngle(90);
  ylabel->DrawLatex(.0125,.475, "P_{T} [GeV]");

  c1->SaveAs("oct26ped_allTS_fib17_ch2.png");
}




void pad8_profile(TProfile *h1, TProfile *h2, TProfile *h3, TProfile *h4, TProfile *h5, TProfile *h6, TProfile *h7, TProfile *h8)
{

  // Plotting
  gStyle->SetOptStat(0);
  auto *c1 = new TCanvas("c1","Title",1600,800);
  TPad* p1 = new TPad("p1","p1",0.01,0.5,0.25,0.99,0); p1->Draw();
  TPad* p2 = new TPad("p2","p2",0.25,0.5,0.50,0.99,0); p2->Draw();
  TPad* p3 = new TPad("p3","p3",0.50,0.5,0.75,0.99,0); p3->Draw();
  TPad* p4 = new TPad("p4","p4",0.75,0.5,0.99,0.99,0); p4->Draw();
  TPad* p5 = new TPad("p5","p5",0.01,0.01,0.25,0.5,0); p5->Draw();
  TPad* p6 = new TPad("p6","p6",0.25,0.01,0.50,0.5,0); p6->Draw();
  TPad* p7 = new TPad("p7","p7",0.50,0.01,0.75,0.5,0); p7->Draw();
  TPad* p8 = new TPad("p8","p8",0.75,0.01,0.99,0.5,0); p8->Draw();
  p1->SetBottomMargin(0);
  p1->SetRightMargin(0);
  p2->SetBottomMargin(0);
  p2->SetLeftMargin(0);
  p2->SetRightMargin(0);
  p3->SetBottomMargin(0);
  p3->SetLeftMargin(0);
  p3->SetRightMargin(0);
  p4->SetBottomMargin(0);
  p4->SetLeftMargin(0);
  p5->SetTopMargin(0);
  p5->SetRightMargin(0);
  p6->SetTopMargin(0);
  p6->SetLeftMargin(0);
  p6->SetRightMargin(0);
  p7->SetTopMargin(0);
  p7->SetLeftMargin(0);
  p7->SetRightMargin(0);
  p8->SetLeftMargin(0);
  p8->SetTopMargin(0);

  TPad *titlePad = new TPad("all","all",0,0,1,1);
  titlePad->SetFillStyle(4000); 
  titlePad->Draw();

  /*TLegend* legend = new TLegend(0.8,0.7,0.6,0.55);
  legend->SetNColumns(1);
  legend->AddEntry(h1, leg1.c_str(), "l");   // l f p 
  legend->AddEntry(h1b, leg2.c_str(), "l");
  legend->SetTextSize(0.05);
  legend->SetBorderSize(0);*/

  p1->cd();
  //p1->SetLogx();
  //p1->SetLogy();
  //float ymin = 0;
  //float ymin = -1.1;
  float ymin = -0.25;
  //float xmax1 = h1->GetXaxis()->GetXmax();
  //float ymax1 = 1.1;
  float xmax1 = h1->GetXaxis()->GetXmax();
  float ymax1 = 1.1;
  //float ymax1=0.65;
  float xmax = 40;
  h1->GetXaxis()->SetRangeUser(0,xmax);
  h1->GetYaxis()->SetRangeUser(ymin,ymax1);
  h1->SetTitle(";;");
  h1->GetYaxis()->SetLabelSize(0.045);
  h1->SetMarkerStyle(kFullDotLarge);
  h1->SetMarkerColor(kBlack);
  h1->SetLineColor(kBlack);
  h1->SetMarkerSize(1);
  h1->Draw("hist p E");
  //legend->Draw("same");
  drawText(.45, .8, "#bf{-2.5 < #eta_{dijet} < -1.5}", kBlack, 22);
  //drawText(.45, .8, "#bf{-2.5 < #eta'_{dijet} < -1.5}", kBlack, 22);
  //drawLine(0, 0, 40, 0);

  p2->cd();
  //p2->SetLogx();
  //p2->SetLogy();
  h2->GetXaxis()->SetRangeUser(0,xmax);
  h2->GetYaxis()->SetRangeUser(ymin,ymax1);
  h2->SetTitle(";;");
  h2->SetMarkerStyle(kFullDotLarge);
  h2->SetMarkerColor(kBlack);
  h2->SetLineColor(kBlack);
  h2->SetMarkerSize(1);
  h2->Draw("hist p E");
  drawText(.45, .8, "#bf{-1.5 < #eta_{dijet} < -1.0}", kBlack, 22);
  //drawText(.45, .8, "#bf{-1.5 < #eta'_{dijet} < -1.0}", kBlack, 22);
  //drawLine(0, 0, 40, 0);

  p3->cd();
  //p3->SetLogx();
  //p3->SetLogy();
  h3->GetXaxis()->SetRangeUser(0,xmax);
  h3->GetYaxis()->SetRangeUser(ymin,ymax1);
  h3->SetTitle(";;");
  h3->SetMarkerStyle(kFullDotLarge);
  h3->SetMarkerColor(kBlack);
  h3->SetLineColor(kBlack);
  h3->SetMarkerSize(1);
  h3->Draw("hist p E");
  drawText(.45, .8, "#bf{-1.0 < #eta_{dijet} < -0.5}", kBlack, 22);
  //drawText(.45, .8, "#bf{-1.0 < #eta'_{dijet} < -0.5}", kBlack, 22);
  //drawLine(0, 0, 40, 0);

  p4->cd();
  //p4->SetLogx();
  //p4->SetLogy();
  h4->GetXaxis()->SetRangeUser(0,xmax);
  h4->GetYaxis()->SetRangeUser(ymin,ymax1);
  h4->SetTitle(";;");
  h4->SetMarkerStyle(kFullDotLarge);
  h4->SetMarkerColor(kBlack);
  h4->SetLineColor(kBlack);
  h4->SetMarkerSize(1);
  h4->Draw("hist p E");
  drawText(.25, .8, "#bf{-0.5 < #eta_{dijet} < 0}", kBlack, 22);
  //drawText(.25, .8, "#bf{-0.5 < #eta'_{dijet} < 0}", kBlack, 22);
  //drawLine(0, 0, 40, 0);

  p5->cd();
  //p5->SetLogx();
  //p5->SetLogy();
  h5->GetXaxis()->SetRangeUser(0,xmax);
  h5->GetYaxis()->SetRangeUser(ymin,ymax1);
  h5->SetTitle(";;");
  h5->GetXaxis()->SetLabelSize(0.0575);
  h5->GetYaxis()->SetLabelSize(0.045);
  h5->SetMarkerStyle(kFullDotLarge);
  h5->SetMarkerColor(kBlack);
  h5->SetLineColor(kBlack);
  h5->SetMarkerSize(1);
  h5->Draw("hist p E");
  drawText(.45, .9, "#bf{0 < #eta_{dijet} < 0.5}", kBlack, 22);
  //drawText(.45, .9, "#bf{0 < #eta'_{dijet} < 0.5}", kBlack, 22);
  //drawLine(0, 0, 40, 0);

  p6->cd();
  //p6->SetLogx();
  //p6->SetLogy();
  h6->GetXaxis()->SetRangeUser(0,xmax);
  h6->GetYaxis()->SetRangeUser(ymin,ymax1);
  h6->SetTitle(";;");
  h6->GetXaxis()->SetLabelSize(0.0575);
  h6->SetMarkerStyle(kFullDotLarge);
  h6->SetMarkerColor(kBlack);
  h6->SetLineColor(kBlack);
  h6->SetMarkerSize(1);
  h6->Draw("hist p E");
  drawText(.35, .9, "#bf{0.5 < #eta_{dijet} < 1.0}", kBlack, 22);
  //drawText(.35, .9, "#bf{0.5 < #eta'_{dijet} < 1.0}", kBlack, 22);
  //drawLine(0, 0, 40, 0);

  p7->cd();
  //p7->SetLogx();
  //p7->SetLogy();
  h7->GetXaxis()->SetRangeUser(0,xmax);
  h7->GetYaxis()->SetRangeUser(ymin,ymax1);
  h7->SetTitle(";;");
  h7->GetXaxis()->SetLabelSize(0.0575);
  h7->SetMarkerStyle(kFullDotLarge);
  h7->SetMarkerColor(kBlack);
  h7->SetLineColor(kBlack);
  h7->SetMarkerSize(1);
  h7->Draw("hist p E");
  drawText(.35, .9, "#bf{1.0 < #eta_{dijet} < 1.5}", kBlack, 22);
  //drawText(.35, .9, "#bf{1.0 < #eta'_{dijet} < 1.5}", kBlack, 22);
  //drawLine(0, 0, 40, 0);

  p8->cd();
  //p8->SetLogx();
  //p8->SetLogy();
  h8->GetXaxis()->SetRangeUser(0,xmax);
  h8->GetYaxis()->SetRangeUser(ymin,ymax1);
  h8->SetTitle(";;");
  h8->GetXaxis()->SetLabelSize(0.0575);
  h8->SetMarkerStyle(kFullDotLarge);
  h8->SetMarkerColor(kBlack);
  h8->SetLineColor(kBlack);
  h8->SetMarkerSize(1);
  h8->Draw("hist p E");
  drawText(.35, .9, "#bf{1.5 < #eta_{dijet} < 2.5}", kBlack, 22);
  //drawText(.35, .9, "#bf{1.5 < #eta'_{dijet} < 2.5}", kBlack, 22);
  //drawLine(0, 0, 40, 0);

  titlePad->cd();
  TLatex *title = new TLatex();
  title->SetTextSize(0.0375);
  title->DrawLatex(.35,.97,"Differential v_{2}: 0nXn");
  //title->DrawLatex(.35,.97,"v_{2}: 0nXn - Xn0n");
  //title->DrawLatex(.35,.97,"Avg v_{2}: 0.5(0nXn + Xn0n)");
  TLatex *xlabel = new TLatex();
  xlabel->SetTextSize(0.03);
  xlabel->DrawLatex(.45,.005, "Q_{T} [GeV]");
  TLatex *ylabel = new TLatex();
  ylabel->SetTextSize(0.03);
  ylabel->SetTextAngle(90);
  ylabel->DrawLatex(.0125,.475, "v_{2} = <cos(2#phi)>");

  c1->SaveAs("oct26ped_allTS_fib17_ch2.png");
}







void pad6_profile(TProfile *h1, TProfile *h2, TProfile *h3, TProfile *h4, TProfile *h5, TProfile *h6)
{

  // Plotting
  gStyle->SetOptStat(0);
  auto *c1 = new TCanvas("c1","Title",1600,800);
  TPad* p1 = new TPad("p1","p1",0.01,0.5,0.33,0.99,0); p1->Draw();
  TPad* p2 = new TPad("p2","p2",0.33,0.5,0.66,0.99,0); p2->Draw();
  TPad* p3 = new TPad("p3","p3",0.66,0.5,0.99,0.99,0); p3->Draw();
  TPad* p4 = new TPad("p4","p4",0.01,0.01,0.33,0.5,0); p4->Draw();
  TPad* p5 = new TPad("p5","p5",0.33,0.01,0.66,0.5,0); p5->Draw();
  TPad* p6 = new TPad("p6","p6",0.66,0.01,0.99,0.5,0); p6->Draw();
  p1->SetBottomMargin(0);
  p1->SetRightMargin(0);
  p2->SetBottomMargin(0);
  p2->SetLeftMargin(0);
  p2->SetRightMargin(0);
  p3->SetBottomMargin(0);
  p3->SetLeftMargin(0);
  p4->SetTopMargin(0);
  p4->SetRightMargin(0);
  p5->SetTopMargin(0);
  p5->SetLeftMargin(0);
  p5->SetRightMargin(0);
  p6->SetTopMargin(0);
  p6->SetLeftMargin(0);

  TPad *titlePad = new TPad("all","all",0,0,1,1);
  titlePad->SetFillStyle(4000); 
  titlePad->Draw();

  /*TLegend* legend = new TLegend(0.8,0.7,0.6,0.55);
  legend->SetNColumns(1);
  legend->AddEntry(h1, leg1.c_str(), "l");   // l f p 
  legend->AddEntry(h1b, leg2.c_str(), "l");
  legend->SetTextSize(0.05);
  legend->SetBorderSize(0);*/

  p1->cd();
  //p1->SetLogx();
  //p1->SetLogy();
  float ymin = 0;
  float xmax1 = h1->GetXaxis()->GetXmax();
  float ymax1 = 1.1;
  float xmax = 40;
  h1->GetXaxis()->SetRangeUser(0,xmax);
  //h1->GetYaxis()->SetRangeUser(ymin,ymax1);
  h1->SetTitle(";;");
  h1->GetYaxis()->SetLabelSize(0.045);
  h1->SetMarkerStyle(kFullDotLarge);
  h1->SetMarkerColor(kBlack);
  h1->SetLineColor(kBlack);
  h1->SetMarkerSize(1);
  h1->Draw("hist p E");
  //legend->Draw("same");
  drawText(.5, .85, "#bf{evt 1, evt 2: p_{T} cuts}", kBlack, 22);
  drawText(.5, .775, "#bf{evt 1, evt 2: P_{T} > Q_{T}}", kBlack, 22);

  p2->cd();
  //p2->SetLogx();
  //p2->SetLogy();
  h2->GetXaxis()->SetRangeUser(0,xmax);
  //h2->GetYaxis()->SetRangeUser(ymin,ymax1);
  h2->SetTitle(";;");
  h2->SetMarkerStyle(kFullDotLarge);
  h2->SetMarkerColor(kBlack);
  h2->SetLineColor(kBlack);
  h2->SetMarkerSize(1);
  h2->Draw("hist p E");
  drawText(.5, .85, "#bf{evt 1, evt 2: p_{T} cuts}", kBlack, 22);
  drawText(.5, .775, "#bf{evt 1, evt 2: P_{T} > Q_{T}}", kBlack, 22);
  drawText(.5, .7, "#bf{Mixed evt: P_{T} > Q_{T}}", kBlack, 22);

  p3->cd();
  //p3->SetLogx();
  //p3->SetLogy();
  h3->GetXaxis()->SetRangeUser(0,xmax);
  //h3->GetYaxis()->SetRangeUser(ymin,ymax1);
  h3->SetTitle(";;");
  h3->SetMarkerStyle(kFullDotLarge);
  h3->SetMarkerColor(kBlack);
  h3->SetLineColor(kBlack);
  h3->SetMarkerSize(1);
  h3->Draw("hist p E");
  drawText(.5, .85, "#bf{evt 1, evt 2: p_{T} cuts}", kBlack, 22);
  drawText(.5, .775, "#bf{Mixed evt: P_{T} > Q_{T}}", kBlack, 22);

  p4->cd();
  //p4->SetLogx();
  //p4->SetLogy();
  h4->GetXaxis()->SetRangeUser(0,xmax);
  //h4->GetYaxis()->SetRangeUser(ymin,ymax1);
  h4->SetTitle(";;");
  h4->SetMarkerStyle(kFullDotLarge);
  h4->SetMarkerColor(kBlack);
  h4->SetLineColor(kBlack);
  h4->GetXaxis()->SetLabelSize(0.0575);
  h4->GetYaxis()->SetLabelSize(0.045);
  h4->SetMarkerSize(1);
  h4->Draw("hist p E");
  drawText(.5, .875, "#bf{Mixed evt: p_{T} cuts}", kBlack, 22);
  drawText(.5, .8, "#bf{evt 1, evt 2: P_{T} > Q_{T}}", kBlack, 22);

  p5->cd();
  //p5->SetLogx();
  //p5->SetLogy();
  h5->GetXaxis()->SetRangeUser(0,xmax);
  //h5->GetYaxis()->SetRangeUser(ymin,ymax1);
  h5->SetTitle(";;");
  h5->GetXaxis()->SetLabelSize(0.0575);
  h5->GetYaxis()->SetLabelSize(0.045);
  h5->SetMarkerStyle(kFullDotLarge);
  h5->SetMarkerColor(kBlack);
  h5->SetLineColor(kBlack);
  h5->SetMarkerSize(1);
  h5->Draw("hist p E");
  drawText(.5, .875, "#bf{Mixed evt: p_{T} cuts}", kBlack, 22);
  drawText(.5, .8, "#bf{evt 1, evt 2: P_{T} > Q_{T}}", kBlack, 22);
  drawText(.5, .725, "#bf{Mixed evt: P_{T} > Q_{T}}", kBlack, 22);

  p6->cd();
  //p6->SetLogx();
  //p6->SetLogy();
  h6->GetXaxis()->SetRangeUser(0,xmax);
  //h6->GetYaxis()->SetRangeUser(ymin,ymax1);
  h6->SetTitle(";;");
  h6->GetXaxis()->SetLabelSize(0.0575);
  h6->SetMarkerStyle(kFullDotLarge);
  h6->SetMarkerColor(kBlack);
  h6->SetLineColor(kBlack);
  h6->SetMarkerSize(1);
  h6->Draw("hist p E");
  drawText(.5, .875, "#bf{Mixed evt: p_{T} cuts}", kBlack, 22);
  drawText(.5, .8, "#bf{Mixed evt: P_{T} > Q_{T}}", kBlack, 22);

  titlePad->cd();
  TLatex *title = new TLatex();
  title->SetTextSize(0.0375);
  title->DrawLatex(.5,.97,"v_{2}");
  TLatex *xlabel = new TLatex();
  xlabel->SetTextSize(0.03);
  xlabel->DrawLatex(.45,.005, "Q_{T} [GeV]");
  TLatex *ylabel = new TLatex();
  ylabel->SetTextSize(0.03);
  ylabel->SetTextAngle(90);
  ylabel->DrawLatex(.0125,.475, "v_{2} = <cos(2#phi)>");

  c1->SaveAs("oct26ped_allTS_fib17_ch2.png");
}







void pad8_8(TH1F *h1, TH1F *h2, TH1F *h3, TH1F *h4, TH1F *h5, TH1F *h6, TH1F *h7, TH1F *h8, TH1F *h1b, TH1F *h2b, TH1F *h3b, TH1F *h4b, TH1F *h5b, TH1F *h6b, TH1F *h7b, TH1F *h8b, string leg1, string leg2)
{

  // Plotting
  gStyle->SetOptStat(0);
  auto *c1 = new TCanvas("c1","Title",1600,800);
  TPad* p1 = new TPad("p1","p1",0.01,0.5,0.25,0.99,0); p1->Draw();
  TPad* p2 = new TPad("p2","p2",0.25,0.5,0.50,0.99,0); p2->Draw();
  TPad* p3 = new TPad("p3","p3",0.50,0.5,0.75,0.99,0); p3->Draw();
  TPad* p4 = new TPad("p4","p4",0.75,0.5,0.99,0.99,0); p4->Draw();
  TPad* p5 = new TPad("p5","p5",0.01,0.01,0.25,0.5,0); p5->Draw();
  TPad* p6 = new TPad("p6","p6",0.25,0.01,0.50,0.5,0); p6->Draw();
  TPad* p7 = new TPad("p7","p7",0.50,0.01,0.75,0.5,0); p7->Draw();
  TPad* p8 = new TPad("p8","p8",0.75,0.01,0.99,0.5,0); p8->Draw();
  p1->SetBottomMargin(0);
  p1->SetRightMargin(0);
  p2->SetBottomMargin(0);
  p2->SetLeftMargin(0);
  p2->SetRightMargin(0);
  p3->SetBottomMargin(0);
  p3->SetLeftMargin(0);
  p3->SetRightMargin(0);
  p4->SetBottomMargin(0);
  p4->SetLeftMargin(0);
  p5->SetTopMargin(0);
  p5->SetRightMargin(0);
  p6->SetTopMargin(0);
  p6->SetLeftMargin(0);
  p6->SetRightMargin(0);
  p7->SetTopMargin(0);
  p7->SetLeftMargin(0);
  p7->SetRightMargin(0);
  p8->SetLeftMargin(0);
  p8->SetTopMargin(0);

  TPad *titlePad = new TPad("all","all",0,0,1,1);
  titlePad->SetFillStyle(4000); 
  titlePad->Draw();

  TLegend* legend = new TLegend(0.8,0.7,0.6,0.55);
  legend->SetNColumns(1);
  legend->AddEntry(h1, leg1.c_str(), "l");   // l f p 
  legend->AddEntry(h1b, leg2.c_str(), "l");
  legend->SetTextSize(0.05);
  legend->SetBorderSize(0);

  p1->cd();
  //p1->SetLogx();
  //p1->SetLogy();
  float ymin = 0;
  float xmax1 = h1->GetXaxis()->GetXmax();
  float ymax1 = h5->GetMaximum();
  float xmax = 30;
  h1->GetXaxis()->SetRangeUser(0,xmax);
  h1->GetYaxis()->SetRangeUser(ymin,1.1*ymax1);
  h1->SetTitle(";;");
  h1->GetYaxis()->SetLabelSize(0.045);
  h1->Draw("hist");
  h1b->Draw("hist same");
  legend->Draw("same");
  drawText(.65, .8, "#bf{-5.2 < #eta < -4.5}", kBlack, 22);

  p2->cd();
  //p2->SetLogx();
  //p2->SetLogy();
  h2->GetXaxis()->SetRangeUser(0,xmax);
  h2->GetYaxis()->SetRangeUser(ymin,1.1*ymax1);
  h2->SetTitle(";;");
  h2->Draw("hist");
  h2b->Draw("hist same");
  drawText(.65, .8, "#bf{-4.5 < #eta < -4.0}", kBlack, 22);

  p3->cd();
  //p3->SetLogx();
  //p3->SetLogy();
  h3->GetXaxis()->SetRangeUser(0,xmax);
  h3->GetYaxis()->SetRangeUser(ymin,1.1*ymax1);
  h3->SetTitle(";;");
  h3->Draw("hist");
  h3b->Draw("hist same");
  drawText(.65, .8, "#bf{-4.0 < #eta < -3.5}", kBlack, 22);

  p4->cd();
  //p4->SetLogx();
  //p4->SetLogy();
  h4->GetXaxis()->SetRangeUser(0,xmax);
  h4->GetYaxis()->SetRangeUser(ymin,1.1*ymax1);
  h4->SetTitle(";;");
  h4->Draw("hist");
  h4b->Draw("hist same");
  drawText(.45, .8, "#bf{-3.5 < #eta < -3.15}", kBlack, 22);

  p5->cd();
  //p5->SetLogx();
  //p5->SetLogy();
  h5->GetXaxis()->SetRangeUser(0,xmax);
  h5->GetYaxis()->SetRangeUser(ymin,1.1*ymax1);
  h5->SetTitle(";;");
  h5->GetXaxis()->SetLabelSize(0.0575);
  h5->GetYaxis()->SetLabelSize(0.045);
  h5->Draw("hist");
  h5b->Draw("hist same");
  drawText(.65, .9, "#bf{3.15 < #eta < 3.5}", kBlack, 22);

  p6->cd();
  //p6->SetLogx();
  //p6->SetLogy();
  h6->GetXaxis()->SetRangeUser(0,xmax);
  h6->GetYaxis()->SetRangeUser(ymin,1.1*ymax1);
  h6->SetTitle(";;");
  h6->GetXaxis()->SetLabelSize(0.0575);
  h6->Draw("hist");
  h6b->Draw("hist same");
  drawText(.55, .9, "#bf{3.5 < #eta < 4.0}", kBlack, 22);

  p7->cd();
  //p7->SetLogx();
  //p7->SetLogy();
  h7->GetXaxis()->SetRangeUser(0,xmax);
  h7->GetYaxis()->SetRangeUser(ymin,1.1*ymax1);
  h7->SetTitle(";;");
  h7->GetXaxis()->SetLabelSize(0.0575);
  h7->Draw("hist");
  h7b->Draw("hist same");
  drawText(.55, .9, "#bf{4.0 < #eta < 4.5}", kBlack, 22);

  p8->cd();
  //p8->SetLogx();
  //p8->SetLogy();
  h8->GetXaxis()->SetRangeUser(0,xmax);
  h8->GetYaxis()->SetRangeUser(ymin,1.1*ymax1);
  h8->SetTitle(";;");
  h8->GetXaxis()->SetLabelSize(0.0575);
  h8->Draw("hist");
  h8b->Draw("hist same");
  drawText(.55, .9, "#bf{4.5 < #eta < 5.2}", kBlack, 22);

  titlePad->cd();
  TLatex *title = new TLatex();
  title->SetTextSize(0.0375);
  title->DrawLatex(.35,.97,"HF+ and HF- Partial Energy Max in #eta (Positive Dijets)");
  TLatex *xlabel = new TLatex();
  xlabel->SetTextSize(0.03);
  xlabel->DrawLatex(.45,.005, "Energy [GeV]");
  TLatex *ylabel = new TLatex();
  ylabel->SetTextSize(0.03);
  ylabel->SetTextAngle(90);
  ylabel->DrawLatex(.0125,.475, "Counts Normalized to Unity");

  c1->SaveAs("oct26ped_allTS_fib17_ch2.png");
}

void pad10_10(TH1F *h1, TH1F *h2, TH1F *h3, TH1F *h4, TH1F *h5, TH1F *h6, TH1F *h7, TH1F *h8, TH1F *h9, TH1F *h10,TH1F *h1_b, TH1F *h2_b, TH1F *h3_b, TH1F *h4_b, TH1F *h5_b, TH1F *h6_b, TH1F *h7_b, TH1F *h8_b, TH1F *h9_b, TH1F *h10_b)
{

  // Plotting
  gStyle->SetOptStat(0);
  auto *c1 = new TCanvas("c1","Title",1600,800);
  c1->SetRightMargin(0.125);
  c1->SetLeftMargin(0.125);
  c1->SetBottomMargin(0.125);
  
  TPad* p1 = new TPad("p1","p1",0.01,0.5,0.2,0.99,0); p1->Draw();
  TPad* p2 = new TPad("p2","p2",0.2,0.5,0.4,0.99,0); p2->Draw();
  TPad* p3 = new TPad("p3","p3",0.4,0.5,0.6,0.99,0); p3->Draw();
  TPad* p4 = new TPad("p4","p4",0.6,0.5,0.8,0.99,0); p4->Draw();
  TPad* p5 = new TPad("p5","p5",0.8,0.5,0.99,0.99,0); p5->Draw();

  TPad* p6 = new TPad("p6","p6",0.01,0.01,0.2,0.5,0); p6->Draw();
  TPad* p7 = new TPad("p7","p7",0.2,0.01,0.4,0.5,0); p7->Draw();
  TPad* p8 = new TPad("p8","p8",0.4,0.01,0.6,0.5,0); p8->Draw();
  TPad* p9 = new TPad("p9","p9",0.6,0.01,0.8,0.5,0); p9->Draw();
  TPad* p10 = new TPad("p10","p10",0.8,0.01,0.99,0.5,0); p10->Draw();
  p1->SetBottomMargin(0);
  p1->SetRightMargin(0);
  p2->SetBottomMargin(0);
  p2->SetLeftMargin(0);
  p2->SetRightMargin(0);
  p3->SetBottomMargin(0);
  p3->SetLeftMargin(0);
  p3->SetRightMargin(0);
  p4->SetBottomMargin(0);
  p4->SetLeftMargin(0);
  p4->SetRightMargin(0);
  p5->SetBottomMargin(0);
  p5->SetLeftMargin(0);
  
  p6->SetTopMargin(0);
  p6->SetRightMargin(0);
  p7->SetTopMargin(0);
  p7->SetLeftMargin(0);
  p7->SetRightMargin(0);
  p8->SetTopMargin(0);
  p8->SetLeftMargin(0);
  p8->SetRightMargin(0);
  p9->SetTopMargin(0);
  p9->SetLeftMargin(0);
  p9->SetRightMargin(0);
  p10->SetTopMargin(0);
  p10->SetLeftMargin(0);

  TPad *titlePad = new TPad("all","all",0,0,1,1);
  titlePad->SetFillStyle(4000); 
  titlePad->Draw();

  string leg1 = "HIEmptyBX";
  string leg2 = "HIForward Data";
  TLegend* legend = new TLegend(0.8,0.7,0.6,0.55);
  legend->SetNColumns(1);
  legend->AddEntry(h1, leg1.c_str(), "l");   // l f p 
  legend->AddEntry(h1_b, leg2.c_str(), "l");
  legend->SetTextSize(0.05);
  legend->SetBorderSize(0);

  p1->cd();
  //p1->SetLogx();
  //p1->SetLogy();
  float ymin = 0;
  float xmax1 = h1->GetXaxis()->GetXmax();
  float ymax1 = (h7->GetMaximum())*1.1;
  float tex_y1 = (h7->GetMaximum())*5.7;
  float tex_x1 = 0.3;
  float tex_y = (h7->GetMaximum())*6.5;
  //float ymax1 = 0.35;
  float xmax = 15;
  h1->GetXaxis()->SetRangeUser(0,xmax);
  h1->GetYaxis()->SetRangeUser(ymin,ymax1);
  h1->GetYaxis()->SetLabelSize(0.045);
  h1->SetTitle(";;");
  h1->Draw("hist");
  h1_b->Draw("hist same");
  drawText(tex_x1, tex_y1, "#bf{-2.5 < #eta < -2.0}", kBlack, 22);
  drawDashLine(5.57895, 0, 5.57895, ymax1);
  legend->Draw("same");

  p2->cd();
  //p2->SetLogx();
  //p2->SetLogy();
  h2->GetXaxis()->SetRangeUser(0,xmax);
  h2->GetYaxis()->SetRangeUser(ymin,ymax1);
  h2->SetTitle(";;");
  h2->Draw("hist");
  h2_b->Draw("hist same");
  drawText(tex_x1, tex_y1, "#bf{-2.0 < #eta < -1.5}", kBlack, 22);
  drawDashLine(1.78947, 0, 1.78947, ymax1);

  p3->cd();
  //p3->SetLogx();
  //p3->SetLogy();
  h3->GetXaxis()->SetRangeUser(0,xmax);
  h3->GetYaxis()->SetRangeUser(ymin,ymax1);
  h3->SetTitle(";;");
  h3->Draw("hist");
  h3_b->Draw("hist same");
  drawText(tex_x1, tex_y1, "#bf{-1.5 < #eta < -1.0}", kBlack, 22);
  drawDashLine(9.57895, 0, 9.57895, ymax1);

  p4->cd();
  //p4->SetLogx();
  //p4->SetLogy();
  h4->GetXaxis()->SetRangeUser(0,xmax);
  h4->GetYaxis()->SetRangeUser(ymin,ymax1);
  h4->SetTitle(";;");
  h4->Draw("hist");
  h4_b->Draw("hist same");
  drawText(tex_x1, tex_y1, "#bf{-1.0 < #eta < -0.5}", kBlack, 22);
  drawDashLine(6.21053, 0, 6.21053, ymax1);

  p5->cd();
  //p5->SetLogx();
  //p5->SetLogy();
  h5->GetXaxis()->SetRangeUser(0,xmax);
  h5->GetYaxis()->SetRangeUser(ymin,ymax1);
  h5->SetTitle(";;");
  h5->Draw("hist");
  h5_b->Draw("hist same");
  drawText(tex_x1, tex_y1, "#bf{-0.5 < #eta < 0}", kBlack, 22);
  drawDashLine(6.42105, 0, 6.42105, ymax1);

  p6->cd();
  //p6->SetLogx();
  //p6->SetLogy();
  h6->GetXaxis()->SetRangeUser(0,xmax);
  h6->GetYaxis()->SetRangeUser(ymin,ymax1);
  h6->GetXaxis()->SetLabelSize(0.06);
  h6->GetYaxis()->SetLabelSize(0.045);
  h6->SetTitle(";;");
  h6->Draw("hist");
  h6_b->Draw("hist same");
  drawText(tex_x1, tex_y, "#bf{0 < #eta < 0.5}", kBlack, 22);
  drawDashLine(5.36842, 0, 5.36842, ymax1);

  p7->cd();
  //p6->SetLogx();
  //p6->SetLogy();
  h7->GetXaxis()->SetRangeUser(0,xmax);
  h7->GetYaxis()->SetRangeUser(ymin,ymax1);
  h7->GetXaxis()->SetLabelSize(0.06);
  h7->SetTitle(";;");
  h7->Draw("hist");
  h7_b->Draw("hist same");
  drawText(tex_x1, tex_y, "#bf{0.5 < #eta < 1.0}", kBlack, 22);
  drawDashLine(5.36842, 0, 5.36842, ymax1);

  p8->cd();
  //p6->SetLogx();
  //p6->SetLogy();
  h8->GetXaxis()->SetRangeUser(0,xmax);
  h8->GetYaxis()->SetRangeUser(ymin,ymax1);
  h8->GetXaxis()->SetLabelSize(0.06);
  h8->SetTitle(";;");
  h8->Draw("hist");
  h8_b->Draw("hist same");
  drawText(tex_x1, tex_y, "#bf{1.0 < #eta < 1.5}", kBlack, 22);
  drawDashLine(14.4211, 0, 14.4211, ymax1);

  p9->cd();
  //p6->SetLogx();
  //p6->SetLogy();
  h9->GetXaxis()->SetRangeUser(0,xmax);
  h9->GetYaxis()->SetRangeUser(ymin,ymax1);
  h9->GetXaxis()->SetLabelSize(0.06);
  h9->SetTitle(";;");
  h9->Draw("hist");
  h9_b->Draw("hist same");
  drawText(tex_x1, tex_y, "#bf{1.5 < #eta < 2.0}", kBlack, 22);
  drawDashLine(1.36842, 0, 1.36842, ymax1);

  p10->cd();
  //p6->SetLogx();
  //p6->SetLogy();
  h10->GetXaxis()->SetRangeUser(0,xmax);
  h10->GetYaxis()->SetRangeUser(ymin,ymax1);
  h10->GetXaxis()->SetLabelSize(0.06);
  h10->SetTitle(";;");
  h10->Draw("hist");
  h10_b->Draw("hist same");
  drawText(tex_x1, tex_y, "#bf{2.0 < #eta < 2.5}", kBlack, 22);
  drawDashLine(6.21053, 0, 6.21053, ymax1);

  titlePad->cd();
  TLatex *title = new TLatex();
  title->SetTextSize(0.0375);
  title->DrawLatex(.42,.97,"Noise Particle Flow Sums");
  TLatex *xlabel = new TLatex();
  xlabel->SetTextSize(0.03);
  xlabel->DrawLatex(.45,.005, "Energy [GeV]");
  TLatex *ylabel = new TLatex();
  ylabel->SetTextSize(0.03);
  ylabel->SetTextAngle(90);
  ylabel->DrawLatex(.0125,.475, "Normalized Counts");

  c1->SaveAs("oct26ped_allTS_fib17_ch2.png");
}


void setErrorHistDivision(TH1 *h1, TH1 *h2, TH1 *h3)
{
  float nBins = h1->GetNbinsX();

  for (int iBin=1; iBin<(nBins+1); iBin++)
    {
      float binContent1 = h1->GetBinContent(iBin);
      float binContent2 = h2->GetBinContent(iBin);
      float binError1 = h1->GetBinError(iBin);
      float binError2 = h2->GetBinError(iBin);

      h3->SetBinError(iBin, sqrt(binError1*binError1/(binContent2*binContent2)+(binError2*binError2*binContent1*binContent1/(binContent2*binContent2*binContent2*binContent2))));
      
    } // end bin loop
}




