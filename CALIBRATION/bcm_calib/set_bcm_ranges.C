#include <TString.h>
#include <TLegend.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TCutG.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TProfile.h>
#include <TObjArray.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void set_bcm_ranges(TString basename="none",Int_t nrun=20277) {
  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11);
  gStyle->SetTitleOffset(1.,"Y");
  gStyle->SetTitleOffset(.7,"X");
  gStyle->SetLabelSize(0.04,"XY");
  gStyle->SetTitleSize(0.06,"XY");
  gStyle->SetPadLeftMargin(0.14);
  TFile *fsimc;
  TString inputroot;
  inputroot="ROOTfiles/scalers/deut_shms_replay_scalers_20277_-1.root";
  cout << " infile root = " << inputroot << endl;
  fsimc =  new TFile(inputroot);
  TTree *tsimc = (TTree*) fsimc->Get("TSP");
  Double_t  Unser;
  tsimc->SetBranchAddress("P.Unser.scalerRate",&Unser);
  Double_t  Time;
  tsimc->SetBranchAddress("P.1MHz.scalerTime",&Time);
  //
  Double_t hmax=1000000;
  TH2F *hUnser_Time = new TH2F("hUnser_Time",Form("Run %d ;   ; )",nrun),3000,0,8000, 1000,0,hmax);
   //
Long64_t nentries = tsimc->GetEntries();

	for (int i = 0; i < nentries; i++) {
     		tsimc->GetEntry(i);
	  hUnser_Time->Fill(Time,Unser);
	}
	//
	vector<Double_t> vxlo;
	vector<Double_t> vxhi;
  TCanvas *cUnser;
     cUnser = new TCanvas("cunser","Unser", 700,700);
     cUnser->Divide(1,1);
     Int_t check=0;
     while (check==0) {
     hUnser_Time->Draw("L");
    for (Int_t j=0;j<vxlo.size();j++) {
      TLine *line1 = new TLine(vxlo[j],0,vxlo[j],hmax);
      TLine *line2 = new TLine(vxhi[j],0,vxhi[j],hmax);
      line1->SetLineColor(2);
      line2->SetLineColor(2);
      line1->Draw();
      line2->Draw();
       cout << vxlo[j] << " " << vxhi[j] << endl;
     }
     
      gPad->Update();
    TCutG *tempg = (TCutG*) gPad->WaitPrimitive("CUTG","CutG");
    //     gPad->Update();
    if (!tempg) cout << " no cut" << endl;
    Double_t zxlo=0.;
    Double_t zxhi=100.;
    Double_t xlo=0.;
    Double_t xhi=100.;
    if (tempg)	{
      Double_t ydummy;
      tempg->GetPoint(0,xlo,ydummy);
      tempg->GetPoint(1,xhi,ydummy);
      //      cout << xlo << " " << xhi  << endl;
      TH2F *hclone =(TH2F*)hUnser_Time->Clone();
      hclone->GetXaxis()->SetRangeUser(xlo,xhi);
      hclone->Draw("L");
      //      TCutG *tempg = (TCutG*) gPad->WaitPrimitive("CUTG","CutG");
       gPad->Update();
       TCutG *cutg = (TCutG*) gPad->WaitPrimitive("CUTG","CutG");
        cutg->GetPoint(0,xlo,ydummy);
        cutg->GetPoint(1,xhi,ydummy);
        cutg->GetPoint(2,zxlo,ydummy);
        cutg->GetPoint(3,zxhi,ydummy);
	//cout << xlo << " " << xhi  << endl;
      vxlo.push_back(zxlo);
      vxhi.push_back(zxhi);
       vxlo.push_back(xlo);
      vxhi.push_back(xhi);
     
       }
    //
    //check=1;
    cout << " to quit set cont = 1 , cont ?" << endl;
    cin >> check ;
       }
     //
     for (Int_t j=0;j<vxlo.size();j++) {
       cout << vxlo[j] << " " << vxhi[j] << endl;
     }
 //
}
