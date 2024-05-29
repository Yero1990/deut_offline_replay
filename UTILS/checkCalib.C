//Script to check if the detectors are calibrated properly
#include <sys/stat.h>
#include "checkCalib.h"

void checkCalib(TString filename="", int run=0, TString hms_pid="", TString shms_pid="")
{

  // filename ---> full .root file path to check calibration
  // run ---> run number
  // hms_pid --> "p" (for protons) or "e" (for electrons)
  // shms_pid ---> "p" or "e"
  
  //Prevent plot 
  gROOT->SetBatch(kTRUE);
  

  //Create a directory where to store the plots and output root file
  mkdir(Form("./hms_Calib_%d", run), S_IRWXU);
  mkdir(Form("./shms_Calib_%d", run), S_IRWXU);
	
  //=============================
  //INITIALIZE HISTOGRAM BINNING
  //=============================

  //HMS			                     //SHMS			     
  hxfp_nbins = 100,                          pxfp_nbins = 100;   
  hxfp_xmin = -50,                           pxfp_xmin = -50;
  hxfp_xmax = 50,                            pxfp_xmax = 50; 

  hyfp_nbins = 100,                          pyfp_nbins = 100;
  hyfp_xmin = -50,                           pyfp_xmin = -50;                                                                                               
  hyfp_xmax = 50,                            pyfp_xmax = 50;  

  hdelta_nbins = 100,                        pdelta_nbins = 100;              
  hdelta_xmin = -10,                         pdelta_xmin = -10;                                                                                                 
  hdelta_xmax = 10,                          pdelta_xmax = 22;        

  hcalEtrkNorm_nbins = 100,         	     pcalEtrkNorm_nbins = 100  ;
  hcalEtrkNorm_xmin = 0.001, 		     pcalEtrkNorm_xmin = 0.001 ;
  hcalEtrkNorm_xmax = 2.0,   		     pcalEtrkNorm_xmax = 2.0;   
			     		   			     
  hcalEtot_nbins = 100,      		     pcalEtot_nbins = 100;      
  hcalEtot_xmin = 0.001,     		     pcalEtot_xmin = 0.001;     
  hcalEtot_xmax = 3.0,	     		     pcalEtot_xmax = 3.0;	     
			     		   			     
  hcalXtrk_nbins = 100,      		     pcalXtrk_nbins = 100;      
  hcalXtrk_xmin = -50, 	     		     pcalXtrk_xmin = -50; 	     
  hcalXtrk_xmax = 50,  	     		     pcalXtrk_xmax = 50;  	     
  			     		     			     
  hcalYtrk_nbins = 100,      		     pcalYtrk_nbins = 100;      
  hcalYtrk_xmin = -50, 	     		     pcalYtrk_xmin = -50; 	     
  hcalYtrk_xmax = 50, 	     		     pcalYtrk_xmax = 50; 	     
  			     		     			     
  hdcTime_nbins = 300,	     		     pdcTime_nbins = 300;	     
  hdcTime_xmin = -50.,	     		     pdcTime_xmin = -50.;	     
  hdcTime_xmax = 250.,	     		     pdcTime_xmax = 250.;	     
  			     		     			     
  hdcRes_nbins = 100.,	     		     pdcRes_nbins = 100.;	     
  hdcRes_xmin = -0.2,	     		     pdcRes_xmin = -0.2;	     
  hdcRes_xmax = 0.2,	     		     pdcRes_xmax = 0.2;	     
  			     		     			     
  hdcDist_nbins = 50.,	     		     pdcDist_nbins = 50.;	     
  hdcDist_xmin = -0.05,	     		     pdcDist_xmin = -0.05;	     
  hdcDist_xmax = 0.55,	     		     pdcDist_xmax = 0.55;	     
			     		   			     
  hhodBeta_nbins = 100,	     		     phodBeta_nbins = 100;	     
  hhodBeta_xmin = 0.3,	     		     phodBeta_xmin = 0.3;	     
  hhodBeta_xmax = 1.6,	     		     phodBeta_xmax = 1.6;	     
  			     		     			     
  hhodXtrk_nbins = 80,	     		     phodXtrk_nbins = 80;	     
  hhodXtrk_xmin = -60,	     		     phodXtrk_xmin = -60;	     
  hhodXtrk_xmax = 60,	     		     phodXtrk_xmax = 60;	     
			     		   			     
  hhodYtrk_nbins = 80,	     		     phodYtrk_nbins = 80;	     
  hhodYtrk_xmin = -60,	     		     phodYtrk_xmin = -60;	     
  hhodYtrk_xmax = 60,	     		     phodYtrk_xmax = 60;	     
  			     		     			     
  hcer_nbins = 100,	     		     phgcer_nbins = 100,  pngcer_nbins = 100;	     
  hcer_xmin = 0.01,	     		     phgcer_xmin = 0.01,  pngcer_xmin = 0.01;     
  hcer_xmax = 10,           		     phgcer_xmax = 20,    pngcer_xmax = 20;        


  //=========================
  //====OPEN ROOT FILE=======
  //=========================
  //TString filename = "../../../ROOTfiles/coin_replay_coin_all_3288_20000.root";
  //TString filename = "../../../ROOTfiles/coin_replay_hod_calib_3288_50000.root"; 
  //TString filename = "../../../ROOTfiles/coin_replay_hdc_calib_3288_-1_hdcCalib.root";
  //TString filename = Form("../../../ROOTfiles/coin_replay_deep_check_%d_-1.root", run);
  //TString filename = Form("../../../ROOTfiles/D2_heep_updated/Emiss_alignment/coin_replay_heep_check_%d_-1.root", run);
  
  //TString filename = Form("../../../ROOTfiles/good_Heep_hmsElec/coin_replay_heep_check_%d_-1.root", run);
  //TString filename = Form("../../../ROOTfiles/good_Heep_hmsElec/g%d_coin.root", run);
  //TString filename = Form("../../../ROOTfiles/good_Heep_hmsProt/hprot_kg%d.root", run);
  //TString filename = Form("../../../ROOTfiles/coin_replay_trkStudy_%d_-1.root", run);

 


  TFile *data_file = new TFile(filename.Data(), "READ");
  TTree *T = (TTree*)data_file->Get("T");
  
  TFile *houtROOT;
  TFile *poutROOT;
  TFile *outROOT;

  //Create output root file where histograms will be stored
  houtROOT = new TFile(Form("./hms_Calib_%d/hms_Calib_histos%d.root", run, run), "recreate");
  poutROOT = new TFile(Form("./shms_Calib_%d/shms_Calib_histos%d.root", run, run), "recreate");
 

  //===========================
  //===Set Branch Address======
  //===========================

  // focal plane in [cm]
  nhxfp = "H.dc.x_fp";                          npxfp = "P.dc.x_fp"; 
  nhyfp = "H.dc.y_fp";                          npyfp = "P.dc.y_fp";  
  nhdelta = "H.gtr.dp";                        npdelta = "P.gtr.dp";    

  nhdc_res = "H.dc.residualExclPlane";	       npdc_res = "P.dc.residualExclPlane";	     
  nhhod_beta = "H.hod.beta";		       nphod_beta = "P.hod.beta";		     
  nhhod_beta_notrk = "H.hod.betanotrack";      nphod_beta_notrk = "P.hod.betanotrack";    
  					       					     
  nhcer_npesum = "H.cer.npeSum";	       nphgcer_npesum = "P.hgcer.npeSum";    npngcer_npesum = "P.ngcer.npeSum";	     
  nhcer_npe = "H.cer.npe";		       nphgcer_npe = "P.hgcer.npe";	     npngcer_npe = "P.ngcer.npe";	     
					     					     
  nhcal_etot = "H.cal.etot";	     	       npcal_etot = "P.cal.etot";	     	     
  nhcal_etrknorm = "H.cal.etottracknorm"; 	       npcal_etrknorm = "P.cal.etottracknorm"; 	     
  nhcal_xtrack = "H.cal.xtrack";	       npcal_xtrack = "P.cal.xtrack";	     
  nhcal_ytrack = "H.cal.ytrack";               npcal_ytrack = "P.cal.ytrack";             
    


  houtROOT->cd();
  //Initialize Some Histograms
  //===HMS====
  H_hbeta_peak = new TH1F("H_hbeta_peak", "HMS Hodoscope Beta (no tracking)", hhodBeta_nbins, hhodBeta_xmin, hhodBeta_xmax);
  
  H_hhodBeta = new TH1F("hHod_Beta", "HMS Hodoscope Beta w/ Tracking", hhodBeta_nbins, hhodBeta_xmin, hhodBeta_xmax);
  H_hhodBetaNoTrk = new TH1F("hHod_Beta_noTracks", "HMS Hodoscope Beta w/out Tracking", hhodBeta_nbins, hhodBeta_xmin, hhodBeta_xmax);
  H_hcalEtrkNorm = new TH1F("hCal_eTrkNorm", "HMS Calorimeter Normalized Track Energy", hcalEtrkNorm_nbins, hcalEtrkNorm_xmin, hcalEtrkNorm_xmax);
  H_hcalEtot = new TH1F("hCal_eTot", "HMS Calorimeter Total Energy", hcalEtot_nbins, hcalEtot_xmin, hcalEtot_xmax);
  H_hcalEtrkNorm_vs_xtrk = new TH2F("hCal_eTrkNorm_v_xtrack", "HMS Norm. Trk E v. xTrack", hcalXtrk_nbins, hcalXtrk_xmin, hcalXtrk_xmax, hcalEtrkNorm_nbins, hcalEtrkNorm_xmin, hcalEtrkNorm_xmax);
  H_hcalEtrkNorm_vs_ytrk = new TH2F("hCal_eTrkNorm_v_ytrack", "HMS Norm. Trk E v. yTrack", hcalYtrk_nbins, hcalYtrk_xmin, hcalYtrk_xmax, hcalEtrkNorm_nbins, hcalEtrkNorm_xmin, hcalEtrkNorm_xmax);
  H_hcerNpeSum = new TH1F("hCer_npeSum", "HMS Cherenkov Total P.E. Sum", hcer_nbins, hcer_xmin, hcer_xmax);

  H_hcalEtrkNorm_vs_xfp = new TH2F("hCal_eTrkNorm_v_xfp", "HMS Edep/Ptrk vs X_{fp}", hxfp_nbins, hxfp_xmin, hxfp_xmax, hcalEtrkNorm_nbins, hcalEtrkNorm_xmin, hcalEtrkNorm_xmax);
  H_hcalEtrkNorm_vs_yfp = new TH2F("hCal_eTrkNorm_v_yfp", "HMS Edep/Ptrk vs Y_{fp}", hyfp_nbins, hyfp_xmin, hyfp_xmax, hcalEtrkNorm_nbins, hcalEtrkNorm_xmin, hcalEtrkNorm_xmax);  
  H_hcalEtrkNorm_vs_delta = new TH2F("hCal_eTrkNorm_v_delta", "HMS Edep/Ptrk vs #delta", hdelta_nbins, hdelta_xmin, hdelta_xmax, hcalEtrkNorm_nbins, hcalEtrkNorm_xmin, hcalEtrkNorm_xmax);  

  H_hhodBeta_vs_xfp = new TH2F("hhodBeta_vs_xfp", "HMS Beta (w/Trk) vs X_{fp}", hxfp_nbins, hxfp_xmin, hxfp_xmax, hhodBeta_nbins, hhodBeta_xmin, hhodBeta_xmax);                   
  H_hhodBeta_vs_yfp = new TH2F("hhodBeta_vs_yfp", "HMS Beta (w/Trk) vs Y_{fp}", hyfp_nbins, hyfp_xmin, hyfp_xmax, hhodBeta_nbins, hhodBeta_xmin, hhodBeta_xmax);                           
  H_hhodBeta_vs_delta = new TH2F("hhodBeta_vs_delta", "HMS Beta (w/Trk) vs #delta", hdelta_nbins, hdelta_xmin, hdelta_xmax, hhodBeta_nbins, hhodBeta_xmin, hhodBeta_xmax);                           
  

  poutROOT->cd();
  //===SHMS===
  H_pbeta_peak = new TH1F("H_pbeta_peak", "SHMS Hodoscope Beta (no tracking)", phodBeta_nbins, phodBeta_xmin, phodBeta_xmax);

  H_phodBeta = new TH1F("pHod_Beta", "SHMS Hodoscope Beta w/ Tracking", phodBeta_nbins, phodBeta_xmin, phodBeta_xmax);
  H_phodBetaNoTrk = new TH1F("pHod_Beta_noTracks", "SHMS Hodoscope Beta w/out Tracking", phodBeta_nbins, phodBeta_xmin, phodBeta_xmax);
  H_pcalEtrkNorm = new TH1F("pCal_eTrkNorm", "SHMS Calorimeter Normalized Track Energy", pcalEtrkNorm_nbins, pcalEtrkNorm_xmin, pcalEtrkNorm_xmax);
  H_pcalEtot = new TH1F("pCal_eTot", "SHMS Calorimeter Total Energy", pcalEtot_nbins, pcalEtot_xmin, pcalEtot_xmax);
  H_pcalEtrkNorm_vs_xtrk = new TH2F("pCal_eTrkNorm_v_xtrack", "SHMS Norm. Trk E v. xTrack", pcalXtrk_nbins, pcalXtrk_xmin, pcalXtrk_xmax, pcalEtrkNorm_nbins, pcalEtrkNorm_xmin, pcalEtrkNorm_xmax);
  H_pcalEtrkNorm_vs_ytrk = new TH2F("pCal_eTrkNorm_v_ytrack", "SHMS Norm. Trk E v. yTrack", pcalYtrk_nbins, pcalYtrk_xmin, pcalYtrk_xmax, pcalEtrkNorm_nbins, pcalEtrkNorm_xmin, pcalEtrkNorm_xmax);  
  H_phgcerNpeSum = new TH1F("phgCer_npeSum", "SHMS HGCER Total P.E. Sum", phgcer_nbins, phgcer_xmin, phgcer_xmax);
  H_pngcerNpeSum = new TH1F("pngCer_npeSum", "SHMS NGCER Total P.E. Sum", pngcer_nbins, pngcer_xmin, pngcer_xmax);

  H_pcalEtrkNorm_vs_xfp = new TH2F("pCal_eTrkNorm_v_xfp", "SHMS Edep/Ptrk vs X_{fp}", pxfp_nbins, pxfp_xmin, pxfp_xmax, pcalEtrkNorm_nbins, pcalEtrkNorm_xmin, pcalEtrkNorm_xmax);                
  H_pcalEtrkNorm_vs_yfp = new TH2F("pCal_eTrkNorm_v_yfp", "SHMS Edep/Ptrk vs Y_{fp}", pyfp_nbins, pyfp_xmin, pyfp_xmax, pcalEtrkNorm_nbins, pcalEtrkNorm_xmin, pcalEtrkNorm_xmax);                                            
  H_pcalEtrkNorm_vs_delta = new TH2F("pCal_eTrkNorm_v_delta", "SHMS Edep/Ptrk vs #delta", pdelta_nbins, pdelta_xmin, pdelta_xmax, pcalEtrkNorm_nbins, pcalEtrkNorm_xmin, pcalEtrkNorm_xmax);                             
                                                                                                                                                                                                                           
  H_phodBeta_vs_xfp = new TH2F("phodBeta_vs_xfp", "SHMS Beta (w/Trk) vs X_{fp}", pxfp_nbins, pxfp_xmin, pxfp_xmax, phodBeta_nbins, phodBeta_xmin, phodBeta_xmax);                                                               
  H_phodBeta_vs_yfp = new TH2F("phodBeta_vs_yfp", "SHMS Beta (w/Trk) vs Y_{fp}", pyfp_nbins, pyfp_xmin, pyfp_xmax, phodBeta_nbins, phodBeta_xmin, phodBeta_xmax);                                                                
  H_phodBeta_vs_delta = new TH2F("phodBeta_vs_delta", "SHMS Beta (w/Trk) vs #delta", pdelta_nbins, pdelta_xmin, pdelta_xmax, phodBeta_nbins, phodBeta_xmin, phodBeta_xmax); 


  //Set Branch Address for Hodo, Calo, and Cherenkov
  
  //HMS
  T->SetBranchAddress(nhxfp, &hxfp);  
  T->SetBranchAddress(nhyfp, &hyfp);
  T->SetBranchAddress(nhdelta, &hdelta);

  T->SetBranchAddress(nhhod_beta, &hhod_beta);
  T->SetBranchAddress(nhhod_beta_notrk, &hhod_beta_notrk);
  T->SetBranchAddress(nhcer_npesum, &hcer_npesum);
  T->SetBranchAddress(nhcal_etot, &hcal_etot);
  T->SetBranchAddress(nhcal_etrknorm, &hcal_etrknorm);
  T->SetBranchAddress(nhcal_xtrack, &hcal_xtrack);
  T->SetBranchAddress(nhcal_ytrack, &hcal_ytrack);
  
 
  //SHMS
  T->SetBranchAddress(npxfp, &pxfp);                                                                                                                                                                          
  T->SetBranchAddress(npyfp, &pyfp);                                                                                                                                                                                         
  T->SetBranchAddress(npdelta, &pdelta);

  T->SetBranchAddress(nphod_beta, &phod_beta);
  T->SetBranchAddress(nphod_beta_notrk, &phod_beta_notrk);
  T->SetBranchAddress(nphgcer_npesum, &phgcer_npesum);
  T->SetBranchAddress(npngcer_npesum, &pngcer_npesum);
  T->SetBranchAddress(npcal_etot, &pcal_etot);
  T->SetBranchAddress(npcal_etrknorm, &pcal_etrknorm);
  T->SetBranchAddress(npcal_xtrack, &pcal_xtrack);
  T->SetBranchAddress(npcal_ytrack, &pcal_ytrack);



  //Loop over DC Planes
  for (Int_t npl = 0; npl < dc_PLANES; npl++ )
    {
      x[npl] = npl+1;  //set x-axis for use with TGraph
      x_err[npl] = 0; // set x-axis error (none)


      houtROOT->cd();    
      //Initialize DC Histograms
      
      //HMS DC
      H_hdcDist[npl] = new TH1F(Form("hDC_%s_DriftDist", hdc_pl_names[npl].c_str()), Form("HMS DC Drift Distance, Plane %s", hdc_pl_names[npl].c_str()), hdcDist_nbins, hdcDist_xmin, hdcDist_xmax);
      H_hdcDist[npl]->GetXaxis()->SetTitle("Drift Distance (cm) ");
      H_hdcDist[npl]->GetXaxis()->CenterTitle();
      H_hdcDist[npl]->GetYaxis()->SetTitle("Counts");
      H_hdcDist[npl]->GetYaxis()->CenterTitle();
      
      H_hdcTime[npl] = new TH1F(Form("hDC_%s_DriftTime", hdc_pl_names[npl].c_str()), Form("HMS DC Drift Time, Plane %s", hdc_pl_names[npl].c_str()), hdcTime_nbins, hdcTime_xmin, hdcTime_xmax);
      H_hdcTime[npl]->GetXaxis()->SetTitle("Drift Time (ns) ");
      H_hdcTime[npl]->GetXaxis()->CenterTitle();
      H_hdcTime[npl]->GetYaxis()->SetTitle("Counts");
      H_hdcTime[npl]->GetYaxis()->CenterTitle();
      
      H_hdcRes[npl] = new TH1F(Form("hDC_%s_DriftResiduals", hdc_pl_names[npl].c_str()), Form("HMS DC Residuals, Plane %s", hdc_pl_names[npl].c_str()), hdcRes_nbins, hdcRes_xmin, hdcRes_xmax);
      H_hdcRes[npl]->GetXaxis()->SetTitle("Drift Residuals (cm) ");
      H_hdcRes[npl]->GetXaxis()->CenterTitle();
      H_hdcRes[npl]->GetYaxis()->SetTitle("Counts");
      H_hdcRes[npl]->GetYaxis()->CenterTitle();
      
      poutROOT->cd();    
      //SHMS DC
      H_pdcDist[npl] = new TH1F(Form("pDC_%s_DriftDist", pdc_pl_names[npl].c_str()), Form("SHMS DC Drift Distance, Plane %s", pdc_pl_names[npl].c_str()), pdcDist_nbins, pdcDist_xmin, pdcDist_xmax);
      H_pdcDist[npl]->GetXaxis()->SetTitle("Drift Distance (cm) ");
      H_pdcDist[npl]->GetXaxis()->CenterTitle();
      H_pdcDist[npl]->GetYaxis()->SetTitle("Counts");
      H_pdcDist[npl]->GetYaxis()->CenterTitle();
      
      H_pdcTime[npl] = new TH1F(Form("pDC_%s_DriftTime", pdc_pl_names[npl].c_str()), Form("SHMS DC Drift Time, Plane %s", pdc_pl_names[npl].c_str()), pdcTime_nbins, pdcTime_xmin, pdcTime_xmax);
      H_pdcTime[npl]->GetXaxis()->SetTitle("Drift Time (ns) ");
      H_pdcTime[npl]->GetXaxis()->CenterTitle();
      H_pdcTime[npl]->GetYaxis()->SetTitle("Counts");
      H_pdcTime[npl]->GetYaxis()->CenterTitle();
      
      H_pdcRes[npl] = new TH1F(Form("pDC_%s_DriftResiduals", pdc_pl_names[npl].c_str()), Form("SHMS DC Residuals, Plane %s", pdc_pl_names[npl].c_str()), pdcRes_nbins, pdcRes_xmin, pdcRes_xmax);
      H_pdcRes[npl]->GetXaxis()->SetTitle("Drift Residuals (cm) ");
      H_pdcRes[npl]->GetXaxis()->CenterTitle();
      H_pdcRes[npl]->GetYaxis()->SetTitle("Counts");
      H_pdcRes[npl]->GetYaxis()->CenterTitle();

      houtROOT->cd();    
      //2D Histos
      //HMS
      H_hres_vs_wire[npl] = new TH2F(Form("hRes_vs_Wire, %s", hdc_pl_names[npl].c_str()), Form("HMS DC Residuals vs. Wire, Plane %s", hdc_pl_names[npl].c_str()), hnwires[npl], 0., hnwires[npl], hdcRes_nbins, hdcRes_xmin, hdcRes_xmax);
      H_hres_vs_wire[npl]->GetXaxis()->SetTitle("Wire Number ");
      H_hres_vs_wire[npl]->GetXaxis()->CenterTitle();
      H_hres_vs_wire[npl]->GetYaxis()->SetTitle("Drift Residuals (cm)");
      H_hres_vs_wire[npl]->GetYaxis()->CenterTitle();
      
      H_htime_vs_wire[npl] = new TH2F(Form("hTime_vs_Wire, %s", hdc_pl_names[npl].c_str()), Form("HMS DC Time vs. Wire, Plane %s", hdc_pl_names[npl].c_str()), hnwires[npl], 0., hnwires[npl], hdcTime_nbins, hdcTime_xmin, hdcTime_xmax);
      H_htime_vs_wire[npl]->GetXaxis()->SetTitle("Wire Number ");
      H_htime_vs_wire[npl]->GetXaxis()->CenterTitle();
      H_htime_vs_wire[npl]->GetYaxis()->SetTitle("Drift Time (ns)");
      H_htime_vs_wire[npl]->GetYaxis()->CenterTitle();
      
      H_hdist_vs_wire[npl] = new TH2F(Form("hDist_vs_Wire, %s", hdc_pl_names[npl].c_str()), Form("HMS DC Distance vs. Wire, Plane %s", hdc_pl_names[npl].c_str()), hnwires[npl], 0., hnwires[npl], hdcDist_nbins, hdcDist_xmin, hdcDist_xmax);
      H_hdist_vs_wire[npl]->GetXaxis()->SetTitle("Wire Number ");
      H_hdist_vs_wire[npl]->GetXaxis()->CenterTitle();
      H_hdist_vs_wire[npl]->GetYaxis()->SetTitle("Drift Distance (cm)");
      H_hdist_vs_wire[npl]->GetYaxis()->CenterTitle();
      
      poutROOT->cd();    
      //SHMS
      H_pres_vs_wire[npl] = new TH2F(Form("pRes_vs_Wire, %s", pdc_pl_names[npl].c_str()), Form("SHMS DC Residuals vs. Wire, Plane %s", pdc_pl_names[npl].c_str()), pnwires[npl], 0., pnwires[npl], pdcRes_nbins, pdcRes_xmin, pdcRes_xmax);
      H_pres_vs_wire[npl]->GetXaxis()->SetTitle("Wire Number ");
      H_pres_vs_wire[npl]->GetXaxis()->CenterTitle();
      H_pres_vs_wire[npl]->GetYaxis()->SetTitle("Drift Residuals (cm)");
      H_pres_vs_wire[npl]->GetYaxis()->CenterTitle();
      
      H_ptime_vs_wire[npl] = new TH2F(Form("pTime_vs_Wire, %s", pdc_pl_names[npl].c_str()), Form("SHMS DC Time vs. Wire, Plane %s", pdc_pl_names[npl].c_str()), pnwires[npl], 0., pnwires[npl], pdcTime_nbins, pdcTime_xmin, pdcTime_xmax);
      H_ptime_vs_wire[npl]->GetXaxis()->SetTitle("Wire Number ");
      H_ptime_vs_wire[npl]->GetXaxis()->CenterTitle();
      H_ptime_vs_wire[npl]->GetYaxis()->SetTitle("Drift Time (ns)");
      H_ptime_vs_wire[npl]->GetYaxis()->CenterTitle();
      
      H_pdist_vs_wire[npl] = new TH2F(Form("pDist_vs_Wire, %s", pdc_pl_names[npl].c_str()), Form("SHMS DC Distance vs. Wire, Plane %s", pdc_pl_names[npl].c_str()), pnwires[npl], 0., pnwires[npl], pdcDist_nbins, pdcDist_xmin, pdcDist_xmax);
      H_pdist_vs_wire[npl]->GetXaxis()->SetTitle("Wire Number ");
      H_pdist_vs_wire[npl]->GetXaxis()->CenterTitle();
      H_pdist_vs_wire[npl]->GetYaxis()->SetTitle("Drift Distance (cm)");
      H_pdist_vs_wire[npl]->GetYaxis()->CenterTitle();

      //===HMS===
      //----Define TTree Leaf Names-----
      base = "H.dc." + hdc_pl_names[npl];
      nhdc_time = base + "." + "time";
      nhdc_wire = base + "." + "wirenum";
      nhdc_dist = base + "." + "dist";
      nhdc_nhit = base + "." + "nhit";
      nhdc_ndata = "Ndata." + base + "." + "time";
      
      //------Set Branch Address-------
      T->SetBranchAddress(nhdc_ndata, &hdc_ndata[npl]);
      T->SetBranchAddress(nhdc_time, hdc_time[npl]);
      T->SetBranchAddress(nhdc_wire, hdc_wire[npl]);
      T->SetBranchAddress(nhdc_dist, hdc_dist[npl]);
      T->SetBranchAddress(nhdc_nhit, &hdc_nhit[npl]);
      T->SetBranchAddress("H.dc.residualExclPlane", &hdc_res[0]);
         
      //===SHMS===
      //----Define TTree Leaf Names-----
      base = "P.dc." + pdc_pl_names[npl];
      npdc_time = base + "." + "time";
      npdc_wire = base + "." + "wirenum";
      npdc_dist = base + "." + "dist";
      npdc_nhit = base + "." + "nhit";
      npdc_ndata = "Ndata." + base + "." + "time";
      
      //------Set Branch Address-------
      T->SetBranchAddress(npdc_ndata, &pdc_ndata[npl]);
      T->SetBranchAddress(npdc_time, pdc_time[npl]);
      T->SetBranchAddress(npdc_wire, pdc_wire[npl]);
      T->SetBranchAddress(npdc_dist, pdc_dist[npl]);
      T->SetBranchAddress(npdc_nhit, &pdc_nhit[npl]);
      T->SetBranchAddress("P.dc.residualExclPlane", &pdc_res[0]);

      

    } //end dc plane loop
  
  //Loop over HODO Planes
  for (Int_t npl = 0; npl < hod_PLANES; npl++ )
    {
      houtROOT->cd();    
      //Initialize HMS HODO Histograms
      H_hhodBeta_v_Xtrk[npl] = new TH2F(Form("hHod_%s_Beta_v_Xtrk", hod_pl_names[npl].c_str()), Form("HMS Hodo Beta vs. X-track, Plane %s", hod_pl_names[npl].c_str()), hhodXtrk_nbins, hhodXtrk_xmin, hhodXtrk_xmax,  hhodBeta_nbins, hhodBeta_xmin, hhodBeta_xmax);
      H_hhodBeta_v_Xtrk[npl]->GetXaxis()->SetTitle("Hodoscope X-Track (cm)");
      H_hhodBeta_v_Xtrk[npl]->GetXaxis()->CenterTitle();
      H_hhodBeta_v_Xtrk[npl]->GetYaxis()->SetTitle("Hodoscope Beta");
      H_hhodBeta_v_Xtrk[npl]->GetYaxis()->CenterTitle();
      
      H_hhodBeta_v_Ytrk[npl] = new TH2F(Form("hHod_%s_Beta_v_Ytrk", hod_pl_names[npl].c_str()), Form("HMS Hodo Beta vs. Y-track, Plane %s", hod_pl_names[npl].c_str()), hhodYtrk_nbins, hhodYtrk_xmin, hhodYtrk_xmax,  hhodBeta_nbins, hhodBeta_xmin, hhodBeta_xmax);
      H_hhodBeta_v_Ytrk[npl]->GetXaxis()->SetTitle("Hodoscope Y-Track (cm)");
      H_hhodBeta_v_Ytrk[npl]->GetXaxis()->CenterTitle();
      H_hhodBeta_v_Ytrk[npl]->GetYaxis()->SetTitle("Hodoscope Beta");
      H_hhodBeta_v_Ytrk[npl]->GetYaxis()->CenterTitle();
      
      poutROOT->cd();    
      //Initialize SHMS HODO Histograms
      H_phodBeta_v_Xtrk[npl] = new TH2F(Form("pHod_%s_Beta_v_Xtrk", hod_pl_names[npl].c_str()), Form("SHMS Hodo Beta vs. X-track, Plane %s", hod_pl_names[npl].c_str()), phodXtrk_nbins, phodXtrk_xmin, phodXtrk_xmax,  phodBeta_nbins, phodBeta_xmin, phodBeta_xmax);
      H_phodBeta_v_Xtrk[npl]->GetXaxis()->SetTitle("Hodoscope X-Track (cm)");
      H_phodBeta_v_Xtrk[npl]->GetXaxis()->CenterTitle();
      H_phodBeta_v_Xtrk[npl]->GetYaxis()->SetTitle("Hodoscope Beta");
      H_phodBeta_v_Xtrk[npl]->GetYaxis()->CenterTitle();
      
      H_phodBeta_v_Ytrk[npl] = new TH2F(Form("pHod_%s_Beta_v_Ytrk", hod_pl_names[npl].c_str()), Form("SHMS Hodo Beta vs. Y-track, Plane %s", hod_pl_names[npl].c_str()), phodYtrk_nbins, phodYtrk_xmin, phodYtrk_xmax,  phodBeta_nbins, phodBeta_xmin, phodBeta_xmax);
      H_phodBeta_v_Ytrk[npl]->GetXaxis()->SetTitle("Hodoscope Y-Track (cm)");
      H_phodBeta_v_Ytrk[npl]->GetXaxis()->CenterTitle();
      H_phodBeta_v_Ytrk[npl]->GetYaxis()->SetTitle("Hodoscope Beta");
      H_phodBeta_v_Ytrk[npl]->GetYaxis()->CenterTitle();

      //----Define HMS TTree Leaf Names-----
      base = "H.hod." + hod_pl_names[npl];
      nhhod_xtrack = "H.hod." + hod_pl_names[npl] + ".TrackXPos";
      nhhod_ytrack = "H.hod." + hod_pl_names[npl] + ".TrackYPos";
      
      //------Set HMS Branch Address-------
      T->SetBranchAddress(nhhod_xtrack, &hhod_xtrack[npl]);
      T->SetBranchAddress(nhhod_ytrack, &hhod_ytrack[npl]);
      
      
      //----Define SHMS TTree Leaf Names-----
      base = "P.hod." + hod_pl_names[npl];
      nphod_xtrack = "P.hod." + hod_pl_names[npl] + ".TrackXPos";
      nphod_ytrack = "P.hod." + hod_pl_names[npl] + ".TrackYPos";
      
      //------Set SHMS Branch Address-------
      T->SetBranchAddress(nphod_xtrack, &phod_xtrack[npl]);
      T->SetBranchAddress(nphod_ytrack, &phod_ytrack[npl]);
    

      if(npl==0)
	{
	  //Loop over Cherenkov PMTs
	  for (int ipmt = 0; ipmt < 4; ipmt++)
	    {
	      poutROOT->cd();    
	      //Initialize HGC/NGC Histos
	      H_phgcerNpe[ipmt] = new TH1F(Form("pHGCER_pmt%d", ipmt+1), Form("SHMS HGC PMT %d", ipmt+1), phgcer_nbins, phgcer_xmin, phgcer_xmax);
	      H_pngcerNpe[ipmt] = new TH1F(Form("pNGCER_pmt%d", ipmt+1), Form("SHMS NGC PMT %d", ipmt+1), pngcer_nbins, pngcer_xmin, pngcer_xmax);
     

	      //------Set Branch Address-------
	      T->SetBranchAddress(nphgcer_npe, phgcer_npe);
	      T->SetBranchAddress(npngcer_npe, pngcer_npe);

	      //HMS Gas Cherenkov
	      if(ipmt<2)
		{
		  houtROOT->cd();    
		  H_hcerNpe[ipmt] = new TH1F(Form("hCER_pmt%d", ipmt+1), Form("HMS Cherenkov PMT %d", ipmt+1), hcer_nbins, hcer_xmin, hcer_xmax);

		  //------Set Branch Address-------
		  T->SetBranchAddress(nhcer_npe, hcer_npe);

		} //end HMS Cer PMT Loop

	      
	    } //end SHMS Cer PMT Loop
	  
	

	} //end Plane==0 requirement, for cer pmt loop

    } //end hodo plane loop
  
  
  Bool_t hnhit;
  Bool_t pnhit;


  // pid cuts used in calibration (must also be applied accordingly here)
  Bool_t beta_cut = 1;
  Bool_t cal_elec = 1;
  Bool_t cer_elec = 1;
  
  Bool_t hms_pid_cut = 1;
  Bool_t shms_pid_cut = 1;

  //==============================
  //=====LOOP OVER ALL EVENTS=====
  //==============================
  
  Long64_t nentries = T->GetEntries();


  //Loop over 50k sample entries (to get beta peak)
  for(Long64_t i=0; i<50000; i++)
    {
      T->GetEntry(i);  
      
      // Get Hod Beta Peak to make pid cuts
      H_hbeta_peak->Fill(hhod_beta_notrk);
      H_pbeta_peak->Fill(phod_beta_notrk);
	

    }

  // bin number corresponding to maximum bin content (this is assumed to be the main beta peak the user is interested in, may be e- or protons)
  int hbinmax = H_hbeta_peak->GetMaximumBin();
  int pbinmax = H_pbeta_peak->GetMaximumBin();
  
  // x-value corresponding to bin number with max content (i.e., peak)
  double hbeta_central = H_hbeta_peak->GetXaxis()->GetBinCenter(hbinmax);
  double pbeta_central = H_pbeta_peak->GetXaxis()->GetBinCenter(pbinmax);
  
  cout << "HMS Beta Peak: " <<  hbeta_central << endl;
  cout << "SHMS Beta Peak: " <<  pbeta_central << endl;
  
  
  //Loop over all entries
  for(Long64_t i=0; i<nentries; i++)
    {
      
      T->GetEntry(i);  
      
      if(hms_pid=="p"){

	//cout << "hms pid: protons " << endl;    
	beta_cut= abs(hbeta_central-hhod_beta_notrk)<0.1;
	hms_pid_cut = beta_cut;
	//cout << "beta_cut = " << beta_cut << endl;
	//cout << Form("(beta_central-beta)=(%.3f, %.3f)", hbeta_central, hhod_beta_notrk) << endl;
      }
      
      else if(hms_pid=="e"){

	//cout << "hms pid: electrons " << endl;    
	cal_elec = hcal_etot > 0.1;  
	cer_elec = hcer_npesum>0.5;
	beta_cut= abs(hbeta_central-hhod_beta_notrk)<0.1;
	hms_pid_cut = cal_elec && cer_elec && beta_cut;
      }

    
      if(shms_pid=="p"){
           
	//cout << "shms pid: protons " << endl; 
       	beta_cut= abs(pbeta_central-phod_beta_notrk)<0.1;
	shms_pid_cut = beta_cut;

      }
      else if (shms_pid=="e"){
	
	//cout << "shms pid: electrons " << endl;
	cal_elec = pcal_etot > 0.1;  
	cer_elec = pngcer_npesum>1.0;
	beta_cut= abs(pbeta_central-phod_beta_notrk)<0.1;
	shms_pid_cut = cal_elec && cer_elec && beta_cut;
      
      }

      

      //======HMS DRIFT CHAMBERS=====
      //Loop over all DC planes
      for (Int_t npl = 0; npl < dc_PLANES; npl++ )
	{
	  

	  //Require single hit per plane
	  hnhit = hdc_nhit[0]==1&&hdc_nhit[1]==1&&hdc_nhit[2]==1&&hdc_nhit[3]==1&&hdc_nhit[4]==1&&hdc_nhit[5]==1&& 
	    hdc_nhit[6]==1&&hdc_nhit[7]==1&&hdc_nhit[8]==1&&hdc_nhit[9]==1&&hdc_nhit[10]==1&&hdc_nhit[11]==1;
	 

	  //Loop over hits
	  for (Int_t j=0; j < hdc_ndata[npl]; j++)
	    {
	      
	      if(hms_pid_cut&&hnhit){
		//Fill Histograms
		H_hdcTime[npl]->Fill(hdc_time[npl][j]);
		H_hdcDist[npl]->Fill(hdc_dist[npl][j]);
		
		//Fill 2D Histos
		H_hres_vs_wire[npl]->Fill(hdc_wire[npl][j], hdc_res[npl]);
		H_htime_vs_wire[npl]->Fill(hdc_wire[npl][j], hdc_time[npl][j]);
		H_hdist_vs_wire[npl]->Fill(hdc_wire[npl][j], hdc_dist[npl][j]);
		
		//Fill Residual
		H_hdcRes[npl]->Fill(hdc_res[npl]);
	      } // end single hit requirement
	      
	    } //end loop over hits
	  
	} //end DC Plane loop
      
      
      //====HMS CALORIMETER=====
      H_hcalEtrkNorm->Fill(hcal_etrknorm);
      H_hcalEtot->Fill(hcal_etot); 
      H_hcalEtrkNorm_vs_xtrk->Fill(hcal_xtrack, hcal_etrknorm);
      H_hcalEtrkNorm_vs_ytrk->Fill(hcal_ytrack, hcal_etrknorm);
      
      H_hcalEtrkNorm_vs_xfp->Fill(hxfp,     hcal_etrknorm);
      H_hcalEtrkNorm_vs_yfp->Fill(hyfp,     hcal_etrknorm);
      H_hcalEtrkNorm_vs_delta->Fill(hdelta, hcal_etrknorm);
      
      

      //====HMS HODOSCOPES=====
      H_hhodBeta->Fill(hhod_beta);
      H_hhodBetaNoTrk->Fill(hhod_beta_notrk);
      
      H_hhodBeta_vs_xfp->Fill(hxfp,     hhod_beta);
      H_hhodBeta_vs_yfp->Fill(hyfp,     hhod_beta);
      H_hhodBeta_vs_delta->Fill(hdelta, hhod_beta);

      //Loop over all HODO planes
      for (Int_t npl = 0; npl < hod_PLANES; npl++ )
	{
	  H_hhodBeta_v_Xtrk[npl]->Fill(hhod_xtrack[npl], hhod_beta);
	  H_hhodBeta_v_Ytrk[npl]->Fill(hhod_ytrack[npl], hhod_beta);


	} //end hodo plane loop      

      
      //---------------------------------------------------------------------
      
      //======sHMS DRIFT CHAMBERS=====
      //Loop over all DC planes
      for (Int_t npl = 0; npl < dc_PLANES; npl++ )
	{
	  //Require single hit per plane
	  pnhit = pdc_nhit[0]==1&&pdc_nhit[1]==1&&pdc_nhit[2]==1&&pdc_nhit[3]==1&&pdc_nhit[4]==1&&pdc_nhit[5]==1&& 
	    pdc_nhit[6]==1&&pdc_nhit[7]==1&&pdc_nhit[8]==1&&pdc_nhit[9]==1&&pdc_nhit[10]==1&&pdc_nhit[11]==1;
	  
	  //Loop over hits
	  for (Int_t j=0; j < pdc_ndata[npl]; j++)
	    {
	      
	      if(shms_pid_cut&&pnhit){
		//Fill Histograms
		H_pdcTime[npl]->Fill(pdc_time[npl][j]);
		H_pdcDist[npl]->Fill(pdc_dist[npl][j]);
		//Fill 2D Histos
		H_pres_vs_wire[npl]->Fill(pdc_wire[npl][j], pdc_res[npl]);
		H_ptime_vs_wire[npl]->Fill(pdc_wire[npl][j], pdc_time[npl][j]);
		H_pdist_vs_wire[npl]->Fill(pdc_wire[npl][j], pdc_dist[npl][j]);
		
		//Fill Residual
		H_pdcRes[npl]->Fill(pdc_res[npl]);
	      } // end single hit requirement
	      
	    } //end loop over hits
	  
	} //end DC Plane loop
      
      
      //====SHMS CALORIMETER=====
      H_pcalEtrkNorm->Fill(pcal_etrknorm);
      H_pcalEtot->Fill(pcal_etot); 
      H_pcalEtrkNorm_vs_xtrk->Fill(pcal_xtrack, pcal_etrknorm);
      H_pcalEtrkNorm_vs_ytrk->Fill(pcal_ytrack, pcal_etrknorm);
      
      H_pcalEtrkNorm_vs_xfp->Fill(pxfp,     pcal_etrknorm);                                                                                                                                                                                       
      H_pcalEtrkNorm_vs_yfp->Fill(pyfp,     pcal_etrknorm);                                                                                                                                                                                               
      H_pcalEtrkNorm_vs_delta->Fill(pdelta, pcal_etrknorm);


      //====SHMS HODOSCOPES=====
      H_phodBeta->Fill(phod_beta);
      H_phodBetaNoTrk->Fill(phod_beta_notrk);
      
      H_phodBeta_vs_xfp->Fill(pxfp,     phod_beta);                                                                                                                                                                                           
      H_phodBeta_vs_yfp->Fill(pyfp,     phod_beta);                                                                                                                                                                                                      
      H_phodBeta_vs_delta->Fill(pdelta, phod_beta);   

      //Loop over all HODO planes
      for (Int_t npl = 0; npl < hod_PLANES; npl++ )
	{
	  H_phodBeta_v_Xtrk[npl]->Fill(phod_xtrack[npl], phod_beta);
	  H_phodBeta_v_Ytrk[npl]->Fill(phod_ytrack[npl], phod_beta);
	  
	  
	} //end hodo plane loop      
      
      //---------------------------------------------------------------------



      //===CHERENKOVS===
      for (int ipmt = 0; ipmt<4; ipmt++)
	{
	  //HMS
	  if(ipmt<2)
	    {
	      H_hcerNpe[ipmt]->Fill(hcer_npe[ipmt]);
	      
	    }//end hms cherenkov pmt loop

	  H_phgcerNpe[ipmt]->Fill(phgcer_npe[ipmt]);
	  H_pngcerNpe[ipmt]->Fill(pngcer_npe[ipmt]);
	}//end shms cherenkov pmt loop

      H_hcerNpeSum->Fill(hcer_npesum);
      H_phgcerNpeSum->Fill(phgcer_npesum);
      H_pngcerNpeSum->Fill(pngcer_npesum);

    }// end loop over entries
  

 
  
  //========================
  //===== DRAW CANVAS ======
  //========================
  
  //===HMS Drift Chambers===
  
  hdcTimeCanv = new TCanvas("hDC Times", "HMS DC TIMES",  2500, 1000);
  hdcTimeCanv->Divide(6,2);
  
  hdcDistCanv = new TCanvas("hDC Dist", "HMS DC Distance",  2500, 1000);
  hdcDistCanv->Divide(6,2);
  
  hdcResCanv = new TCanvas("hDC Residuals", "HMS DC Residuals",  2500, 1000);
  hdcResCanv->SetLeftMargin(0.25);
  hdcResCanv->Divide(6,2);
  
  hdcResGraphCanv = new TCanvas("hDC Residuals Graph", "HMS DC Residuals Graph",  2500, 1000);
  hdcResGraphCanv->Divide(2,1);
  
  //2d histos 
  hdcResCanv2D = new TCanvas("hDC Residuals vs. Wire", "HMS DC Residuals vs. Wire",  2500, 1000);
  hdcResCanv2D->Divide(6,2);
  
  hdcTimeCanv2D = new TCanvas("hDC Time vs. Wire", "HMS DC Time vs. Wire",  2500, 1000);
  hdcTimeCanv2D->Divide(6,2);
  
  hdcDistCanv2D = new TCanvas("hDC Dist vs. Wire", "HMS DC Dist vs. Wire",  2500, 1000);
  hdcDistCanv2D->Divide(6,2);
  
  //Profile Histograms
  hdcResCanvProf = new TCanvas("hDC Residuals vs. Wire: Profile", "HMS DC Residuals vs. Wire, Profile",  2500, 1000);
  hdcResCanvProf->Divide(6,2);
      
  //Loop over DC planes
  for (Int_t npl = 0; npl < dc_PLANES; npl++ )
    {
      
      hdcTimeCanv->cd(npl+1);
      gPad->SetLogy();
      H_hdcTime[npl]->Draw();
      
      hdcDistCanv->cd(npl+1);
      binmax = H_hdcDist[npl]->GetMaximumBin();
      upperlim =  H_hdcDist[npl]->GetBinContent(binmax) + 200.;
      H_hdcDist[npl]->GetYaxis()->SetRangeUser(0., upperlim);
      H_hdcDist[npl]->Draw();
      
      
      // determine gaussian fit limits for residuals
      binmax = H_hdcRes[npl]->GetMaximumBin(); 	  
      stdev  = H_hdcRes[npl]->GetStdDev(); 
      amplitude = H_hdcRes[npl]->GetBinContent(binmax);
      center = H_hdcRes[npl]->GetBinCenter(binmax);
      xmin_fit = center - (1.5*stdev);
      xmax_fit = center + (1.5*stdev); 
      
      cout << std::setprecision(3) << "HMS DC xmin_fit, xmax_fit = " << xmin_fit << ", " << xmax_fit << endl;
      
      
      //change to pad of canvas
      hdcResCanv->cd(npl+1);
      
      // define fit fucntion
      gaus_fit = new TF1("gaus_fit","gaus", xmin_fit, xmax_fit);
      gaus_fit->SetParameters(amplitude, center, stdev );
      
      // Fit gaussian to the residuals
      H_hdcRes[npl]->Fit("gaus_fit", "R");	  
      H_hdcRes[npl]->Draw();
      
      //Get Mean/Sigma for residuals and conver to microns                                                                                                                                       
      mean[npl]      = gaus_fit->GetParameter(1)*1e4;                                                                                                                       
      mean_err[npl]  = gaus_fit->GetParError(1)*1e4; 
      sigma[npl]     = gaus_fit->GetParameter(2)*1e4; 
      sigma_err[npl] = gaus_fit->GetParError(2)*1e4;
      
      //2D and Profile Histograms
      hdcResCanv2D->cd(npl+1);
      H_hres_vs_wire[npl]->Draw("COLZ");
      
      hdcTimeCanv2D->cd(npl+1);
      H_htime_vs_wire[npl]->Draw("COLZ");
      
      hdcDistCanv2D->cd(npl+1);
      H_hdist_vs_wire[npl]->Draw("COLZ");
      
      
      houtROOT->cd(); 
      hdcResCanvProf->cd(npl+1);
      hdcResProf[npl] = H_hres_vs_wire[npl]->ProfileX(Form("HMS Profile of Residuals, Plane %s", hdc_pl_names[npl].c_str()), 0., hnwires[npl]);
      hdcResProf[npl]->Draw();
      
    } //END LOOP OVER DC PLANES
  
  
  TGraph *hgr_mean = new TGraphErrors(12, x, mean, x_err, mean_err);
  
  //Change to SupPad 2 to plot mean
  hdcResGraphCanv->cd(1);
  //dcResGraphCanv->SetGrid();
  hgr_mean->SetTitle("HMS DC Residuals Mean");
  hgr_mean->SetMarkerStyle(22);
  hgr_mean->SetMarkerColor(kBlue);
  hgr_mean->SetMarkerSize(1);
  hgr_mean->GetYaxis()->SetRangeUser(-250, 250);
  
  //Set Axis Titles
  hgr_mean->GetXaxis()->SetTitle("HMS DC Planes Residuals");
  hgr_mean->GetXaxis()->CenterTitle();
  hgr_mean->GetYaxis()->SetTitle("HMS DC Residual Mean (#mum)");
  hgr_mean->GetYaxis()->CenterTitle();
  hgr_mean->SetTitle("HMS DC Plane Residuals Mean");
  
  hgr_mean->Draw("AP");
  
  
  //Change to SubPad 1 to plot sigma
  hdcResGraphCanv->cd(2);
  TGraph *hgr_residual = new TGraphErrors(12, x, sigma, x_err, sigma_err);
  //dcResGraphCanv->SetGrid();
  hgr_residual->SetTitle("HMS DC Residuals Sigma");
  hgr_residual->SetMarkerStyle(22);
  hgr_residual->SetMarkerColor(kRed);
  hgr_residual->SetMarkerSize(1);
  hgr_residual->GetYaxis()->SetRangeUser(0, 1000.);
  
  //Set Axis Titles
  hgr_residual->GetXaxis()->SetTitle("HMS DC Planes Residuals");
  hgr_residual->GetXaxis()->CenterTitle();
  hgr_residual->GetYaxis()->SetTitle("HMS DC Residual #sigma (#mum)");
  hgr_residual->GetYaxis()->CenterTitle();
  hgr_residual->SetTitle("HMS DC Plane Residuals Sigma");
  hgr_residual->Draw("AP");
  hdcResGraphCanv->Update();
  
  hdcTimeCanv->SaveAs(Form("./hms_Calib_%d/hDC_Times.pdf", run));
  hdcDistCanv->SaveAs(Form("./hms_Calib_%d/hDC_Dist.pdf", run));
  hdcResCanv->SaveAs(Form("./hms_Calib_%d/hDC_Res.pdf", run));
  
  hdcResCanv2D->SaveAs(Form("./hms_Calib_%d/hDC_Res2D.pdf", run));
  hdcTimeCanv2D->SaveAs(Form("./hms_Calib_%d/hDC_Time2D.pdf", run));
  hdcDistCanv2D->SaveAs(Form("./hms_Calib_%d/hDC_Dist2D.pdf", run));
  
  
  hdcResCanvProf->SaveAs(Form("./hms_Calib_%d/hDC_ResProfile.pdf", run));
  hdcResGraphCanv->SaveAs(Form("./hms_Calib_%d/hDC_ResPlot.pdf", run));
  
  //====CALORIMETERS====
  hcalCanv = new TCanvas("HMS Calorimeter Canv" , "HMS Calorimeter Plots", 1500, 500);
  hcalCanv->Divide(3,2);
  
  hcalCanv->cd(1);
  H_hcalEtrkNorm->Draw();
  hcalCanv->Update();
  
  hcalCanv->cd(2);
  H_hcalEtot->Draw();
  
  hcalCanv->cd(3);
  H_hcalEtrkNorm_vs_xtrk->Draw("COLZ");
  
  hcalCanv->cd(4);
  H_hcalEtrkNorm_vs_ytrk->Draw("COLZ");
  
  
  hcalCanv->SaveAs(Form("./hms_Calib_%d/hCal_CalibPlots.pdf", run));

  
  //======HODOSCOPOES=====
  hhodCanv = new TCanvas("HMS Hodoscope Beta Canv" , "HMS Hodoscope Beta Plots", 1500, 500);
  hhodCanv->Divide(2,1);
  hhodCanv->cd(1);
  H_hhodBetaNoTrk->Draw();
  hhodCanv->cd(2);
  H_hhodBeta->Draw();    
  
  
  hhodCanv2D = new TCanvas("HMS Hodoscope 2D Beta Canvas" , "2D HMS Hodoscope Beta Plots", 1500, 500);
  hhodCanv2D->Divide(4,2);
  
  hhodProfCanv = new TCanvas("Hodoscope Beta Profile Canvas", "Hodoscope Beta ProfileX", 1500, 500);
  hhodProfCanv->Divide(4,2);
  
  
  //Loop over Hodo planes
  for (Int_t npl = 0; npl < hod_PLANES; npl++ )
    {
      
      hhodCanv2D->cd(npl + 1);
      H_hhodBeta_v_Xtrk[npl]->Draw("COLZ");
      
      hhodCanv2D->cd(npl + 5);
      H_hhodBeta_v_Ytrk[npl]->Draw("COLZ");
      
      houtROOT->cd();
      //X-Profile of 2D Histos beta vs xtrk (or ytrk)
      hhod_xProfX[npl] = new TProfile();
      hhod_yProfX[npl] = new TProfile();
      
      hhodProfCanv->cd(npl + 1);
      hhod_xProfX[npl] = H_hhodBeta_v_Xtrk[npl]->ProfileX(Form("h%s_betaXtrk_ProfileX", hod_pl_names[npl].c_str()), 0, 200);
      hhod_xProfX[npl]->Draw();
      hhodProfCanv->cd(npl + 5);
      hhod_yProfX[npl] = H_hhodBeta_v_Ytrk[npl]->ProfileX(Form("h%s_betaYtrk_ProfileX", hod_pl_names[npl].c_str()), 0, 200);
      hhod_yProfX[npl]->Draw();
      
    } // END LOOP OVER HODO PLANES
  
  hhodCanv->SaveAs(Form("./hms_Calib_%d/hHodBetaPlots.pdf", run));
  hhodCanv2D->SaveAs(Form("./hms_Calib_%d/hHodBeta2DPlots.pdf", run));
  hhodProfCanv->SaveAs(Form("./hms_Calib_%d/hHodBetaProfilePlots.pdf", run));

  
  //====CHERENKOV====
  hcerCanv = new TCanvas("HMS Cherenkov", "HMS Cherenkov Calib. Plots", 1500, 500);
  hcerCanv->Divide(3,1);
  hcerCanv->cd(1);
  H_hcerNpe[0]->Draw();
  hcerCanv->cd(2);
  H_hcerNpe[1]->Draw();
  hcerCanv->cd(3);
  H_hcerNpeSum->Draw();
  
  hcerCanv->SaveAs(Form("./hms_Calib_%d/hCer.pdf", run));
      
  //Write Histograms to ROOT file 
  houtROOT->cd();
  houtROOT->Write(); 
  houtROOT->WriteTObject(hgr_mean);
  houtROOT->WriteTObject(hgr_residual);
  houtROOT->Close();

  
  //===SHMS Drift Chambers===
  
  pdcTimeCanv = new TCanvas("pDC Times", "SHMS DC TIMES",  1500, 500);
  pdcTimeCanv->Divide(6,2);
  
  pdcDistCanv = new TCanvas("pDC Dist", "SHMS DC Distance",  1500, 500);
  pdcDistCanv->Divide(6,2);
  
  pdcResCanv = new TCanvas("pDC Residuals", "SHMS DC Residuals",  1500, 500);
  pdcResCanv->Divide(6,2);
  
  pdcResGraphCanv = new TCanvas("pDC Residuals Graph", "SHMS DC Residuals Graph",  1500, 500);
  pdcResGraphCanv->Divide(2,1);
  
  //2d histos 
  pdcResCanv2D = new TCanvas("pDC Residuals vs. Wire", "SHMS DC Residuals vs. Wire",  1500, 500);
  pdcResCanv2D->Divide(6,2);
  
  pdcTimeCanv2D = new TCanvas("pDC Time vs. Wire", "SHMS DC Time vs. Wire",  1500, 500);
  pdcTimeCanv2D->Divide(6,2);
  
  pdcDistCanv2D = new TCanvas("pDC Dist vs. Wire", "SHMS DC Dist vs. Wire",  1500, 500);
  pdcDistCanv2D->Divide(6,2);
  
  //Profile Histograms
  pdcResCanvProf = new TCanvas("pDC Residuals vs. Wire: Profile", "SHMS DC Residuals vs. Wire, Profile",  1500, 500);
  pdcResCanvProf->Divide(6,2);
  
  //Loop over DC planes
  for (Int_t npl = 0; npl < dc_PLANES; npl++ )
    {
      
      pdcTimeCanv->cd(npl+1);
      gPad->SetLogy();
      H_pdcTime[npl]->Draw();
      
      pdcDistCanv->cd(npl+1);
      binmax = H_pdcDist[npl]->GetMaximumBin();
      upperlim =  H_pdcDist[npl]->GetBinContent(binmax) + 200.;
      H_pdcDist[npl]->GetYaxis()->SetRangeUser(0., upperlim);
      H_pdcDist[npl]->Draw();
      
      
      // determine gaussian fit limits for residuals
      binmax = H_pdcRes[npl]->GetMaximumBin(); 	  
      stdev  = H_pdcRes[npl]->GetStdDev(); 
      amplitude = H_pdcRes[npl]->GetBinContent(binmax);
      center = H_pdcRes[npl]->GetBinCenter(binmax);
      xmin_fit = center - (1.5*stdev);
      xmax_fit = center + (1.5*stdev); 
      
      cout << std::setprecision(3) << "SHMS DC xmin_fit, xmax_fit = " << xmin_fit << ", " << xmax_fit << endl;
      
      //change to pad of canvas
      pdcResCanv->cd(npl+1);
      
      // define fit fucntion
      gaus_fit = new TF1("gaus_fit","gaus", xmin_fit, xmax_fit);
      gaus_fit->SetParameters(amplitude, center, stdev );
      
      // Fit gaussian to the residuals
      H_pdcRes[npl]->Fit("gaus_fit", "R");	  
      H_pdcRes[npl]->Draw();
      
      
      //Get Mean/Sigma for residuals and conver to microns                                                                                                                                       
      mean[npl]      = gaus_fit->GetParameter(1)*1e4;            
      mean_err[npl]  = gaus_fit->GetParError(1)*1e4; 
      sigma[npl]     = gaus_fit->GetParameter(2)*1e4; 
      sigma_err[npl] = gaus_fit->GetParError(2)*1e4; 
      
      //cout << "fit[um]: mean, err = " <<  mean[npl] << ", "<<  mean_err[npl] << endl; 
      //cout << "fit[um]: sigma, err = " <<  sigma[npl] << ", "<<  sigma_err[npl] << endl; 
      
      //2D and Profile Histograms
      pdcResCanv2D->cd(npl+1);
      H_pres_vs_wire[npl]->Draw("COLZ");
      
      pdcTimeCanv2D->cd(npl+1);
      H_ptime_vs_wire[npl]->Draw("COLZ");
      
      pdcDistCanv2D->cd(npl+1);
      H_pdist_vs_wire[npl]->Draw("COLZ");
      
      poutROOT->cd(); 
      pdcResCanvProf->cd(npl+1);
      pdcResProf[npl] = H_pres_vs_wire[npl]->ProfileX(Form("SHMS Profile of Residuals, Plane %s", pdc_pl_names[npl].c_str()), 0., pnwires[npl]);
      pdcResProf[npl]->Draw();
      
    } //END LOOP OVER DC PLANES
  
  
  
  TGraph *pgr_mean = new TGraphErrors(12, x, mean, x_err, mean_err);
  
  //Change to SupPad 2 to plot mean
  pdcResGraphCanv->cd(1);
  //dcResGraphCanv->SetGrid();
  pgr_mean->SetTitle("SHMS DC Residuals Mean");   
  pgr_mean->SetMarkerStyle(22);
  pgr_mean->SetMarkerColor(kBlue);
  pgr_mean->SetMarkerSize(1);
  pgr_mean->GetYaxis()->SetRangeUser(-250, 250);
  
  //Set Axis Titles
  pgr_mean->GetXaxis()->SetTitle("SHMS DC Planes Residuals");
  pgr_mean->GetXaxis()->CenterTitle();
  pgr_mean->GetYaxis()->SetTitle("SHMS DC Residual Mean (#mum)");
  pgr_mean->GetYaxis()->CenterTitle();
  pgr_mean->SetTitle("SHMS DC Plane Residuals Mean");
  
  pgr_mean->Draw("AP");
  
  
  //Change to SubPad 1 to plot sigma
  pdcResGraphCanv->cd(2);
  TGraph *pgr_residual = new TGraphErrors(12, x, sigma, x_err, sigma_err);
  //dcResGraphCanv->SetGrid();
  pgr_residual->SetTitle("SHMS DC Residuals Sigma");
  pgr_residual->SetMarkerStyle(22);
  pgr_residual->SetMarkerColor(kRed);
  pgr_residual->SetMarkerSize(1);
  pgr_residual->GetYaxis()->SetRangeUser(0, 1000.);
  
  //Set Axis Titles
  pgr_residual->GetXaxis()->SetTitle("SHMS DC Planes Residuals");
  pgr_residual->GetXaxis()->CenterTitle();
  pgr_residual->GetYaxis()->SetTitle("SHMS DC Residual #sigma (#mum)");
  pgr_residual->GetYaxis()->CenterTitle();
  pgr_residual->SetTitle("SHMS DC Plane Residuals Sigma");
  pgr_residual->Draw("AP");
  pdcResGraphCanv->Update();
  
  pdcTimeCanv->SaveAs(Form("./shms_Calib_%d/pDC_Times.pdf",  run));
  pdcDistCanv->SaveAs(Form("./shms_Calib_%d/pDC_Dist.pdf",  run));
  pdcResCanv->SaveAs(Form("./shms_Calib_%d/pDC_Res.pdf",  run));
  
  pdcResCanv2D->SaveAs(Form("./shms_Calib_%d/pDC_Res2D.pdf",  run));
  pdcTimeCanv2D->SaveAs(Form("./shms_Calib_%d/pDC_Time2D.pdf",  run));
  pdcDistCanv2D->SaveAs(Form("./shms_Calib_%d/pDC_Dist2D.pdf",  run));
  
  
  pdcResCanvProf->SaveAs(Form("./shms_Calib_%d/pDC_ResProfile.pdf",  run));
  pdcResGraphCanv->SaveAs(Form("./shms_Calib_%d/pDC_ResPlot.pdf",  run));

  
  //====CALORIMETERS====
  pcalCanv = new TCanvas("SHMS Calorimeter Canv" , "SHMS Calorimeter Plots", 1500, 500);
  pcalCanv->Divide(3,2);
      
  pcalCanv->cd(1);
  H_pcalEtrkNorm->Draw();
  pcalCanv->Update();
  
  pcalCanv->cd(2);
  H_pcalEtot->Draw();
  
  pcalCanv->cd(3);
  H_pcalEtrkNorm_vs_xtrk->Draw("COLZ");
  
  pcalCanv->cd(4);
  H_pcalEtrkNorm_vs_ytrk->Draw("COLZ");
  

  pcalCanv->SaveAs(Form("./shms_Calib_%d/pCal_CalibPlots.pdf",  run));

  
  //======HODOSCOPOES=====
  phodCanv = new TCanvas("SHMS Hodoscope Beta Canv" , "SHMS Hodoscope Beta Plots", 1500, 500);
  phodCanv->Divide(2,1);
  phodCanv->cd(1);
  H_phodBetaNoTrk->Draw();
  phodCanv->cd(2);
  H_phodBeta->Draw();    
  
  
  phodCanv2D = new TCanvas("SHMS Hodoscope 2D Beta Canvas" , "2D SHMS Hodoscope Beta Plots", 1500, 500);
  phodCanv2D->Divide(4,2);
  
  phodProfCanv = new TCanvas("SHMS Hodoscope Beta Profile Canvas", "SHMS Hodoscope Beta ProfileX", 1500, 500);
  phodProfCanv->Divide(4,2);
  
  
  //Loop over Hodo planes
  for (Int_t npl = 0; npl < hod_PLANES; npl++ )
    {
      
      phodCanv2D->cd(npl + 1);
      H_phodBeta_v_Xtrk[npl]->Draw("COLZ");
	  
      phodCanv2D->cd(npl + 5);
      H_phodBeta_v_Ytrk[npl]->Draw("COLZ");
      
      poutROOT->cd();
      //X-Profile of 2D Histos beta vs xtrk (or ytrk)
      phod_xProfX[npl] = new TProfile();
      phod_yProfX[npl] = new TProfile();
      
      phodProfCanv->cd(npl + 1);
      phod_xProfX[npl] = H_phodBeta_v_Xtrk[npl]->ProfileX(Form("p%s_betaXtrk_ProfileX", hod_pl_names[npl].c_str()), 0, 200);
      phod_xProfX[npl]->Draw();
      phodProfCanv->cd(npl + 5);
      phod_yProfX[npl] = H_phodBeta_v_Ytrk[npl]->ProfileX(Form("p%s_betaYtrk_ProfileX", hod_pl_names[npl].c_str()), 0, 200);
      phod_yProfX[npl]->Draw();
      
    } // END LOOP OVER HODO PLANES
    
  phodCanv->SaveAs(Form("./shms_Calib_%d/pHodBetaPlots.pdf",  run));
  phodCanv2D->SaveAs(Form("./shms_Calib_%d/pHodBeta2DPlots.pdf",  run));
  phodProfCanv->SaveAs(Form("./shms_Calib_%d/pHodBetaProfilePlots.pdf",  run));

  
  //====CHERENKOVS====
  phgcerCanv = new TCanvas("SHMS HGCER", "SHMS HGCER Calib. Plots", 1500, 1500);
  pngcerCanv = new TCanvas("SHMS NGCER", "SHMS NGCER Calib. Plots", 1500, 1500);
  
  phgcerSumCanv = new TCanvas("SHMS HGCER SUM", "SHMS HGCER SUM Calib. Plots", 1500, 500);
  pngcerSumCanv = new TCanvas("SHMS NGCER SUM", "SHMS NGCER SUM Calib. Plots", 1500, 500);

  phgcerCanv->Divide(2,2);
  pngcerCanv->Divide(2,2);

  //HGCER
  phgcerCanv->cd(1);
  H_phgcerNpe[0]->Draw();
  phgcerCanv->cd(2);
  H_phgcerNpe[1]->Draw();
  phgcerCanv->cd(3);
  H_phgcerNpe[2]->Draw();
  phgcerCanv->cd(4);
  H_phgcerNpe[3]->Draw();

  phgcerSumCanv->cd();
  H_phgcerNpeSum->Draw();

  //NGCER
  pngcerCanv->cd(1);
  H_pngcerNpe[0]->Draw();
  pngcerCanv->cd(2);
  H_pngcerNpe[1]->Draw();
  pngcerCanv->cd(3);
  H_pngcerNpe[2]->Draw();
  pngcerCanv->cd(4);
  H_pngcerNpe[3]->Draw();
  
  pngcerSumCanv->cd();
  H_pngcerNpeSum->Draw();
  
  phgcerCanv->SaveAs(Form("./shms_Calib_%d/pHGCER_PMTs.pdf",  run));
  pngcerCanv->SaveAs(Form("./shms_Calib_%d/pNGCER_PMTs.pdf",  run));

  phgcerSumCanv->SaveAs(Form("./shms_Calib_%d/pHGCER_SUM.pdf",  run));
  pngcerSumCanv->SaveAs(Form("./shms_Calib_%d/pNGCER_SUM.pdf",  run));


  //Write Histograms to ROOT file
  poutROOT->cd();
  poutROOT->Write(); 
  poutROOT->WriteTObject(pgr_mean);
  poutROOT->WriteTObject(pgr_residual);    
  poutROOT->Close();
	


}




