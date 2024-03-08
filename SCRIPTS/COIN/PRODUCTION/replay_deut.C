#include <chrono>
#include <iostream>
#include "../../../UTILS/parse_utils.h" //useful C++ string parsing utilities

using namespace std;
using namespace std::chrono;


void replay_deut(Int_t RunNumber = 0, Int_t MaxEvent = 0, TString ftype="") {

  
  // Get RunNumber and MaxEvent if not provided.
  if(RunNumber == 0) {
    cout << "Enter a Run Number (-1 to exit): ";
    cin >> RunNumber;
    if( RunNumber<=0 ) return;
  }
  if(MaxEvent == 0) {
    cout << "\nNumber of Events to analyze: ";
    cin >> MaxEvent;
    if(MaxEvent == 0) {
      cerr << "...Invalid entry\n";
      //exit;
    }
  }
  
  if(ftype==""){
    cout  << "\nEnter analysis type to use (e.g., shms50k, hms50k, sample, prod): \n " << endl;
    cin >> ftype;
    if(ftype==""){
      cerr << "...Invalid file type\n";
      //exit;
    }
  }
  
  // Get starting timepoint
  auto start = high_resolution_clock::now();
    
  // Create file name patterns.
  //const char* RunFileNamePattern = "coin_all_%05d.dat";
  const char* RunFileNamePattern = "shms_all_%05d.dat";

  vector<TString> pathList;
  pathList.push_back(".");
  pathList.push_back("./raw");
  pathList.push_back("./raw.copiedtotape");
  pathList.push_back("./CACHE_LINKS/cache_pionlt");
  pathList.push_back("./CACHE_LINKS/cache_cafe"); 
  pathList.push_back("./CACHE_LINKS/cache_deut");

  //const char* RunFileNamePattern = "raw/coin_all_%05d.dat";

  // Create dir. to store monitoring histos
  TString cmd = ""; // Form("mkdir -p ROOTfiles/%s", ftype.Data());

  if((ftype=="prod") && (MaxEvent==-1)){
    ftype="prod";
    cmd = Form("mkdir -p ROOTfiles/%s", ftype.Data());
  }
  else if( (ftype=="prod") && (MaxEvent!=-1) ){
    ftype="sample";
    cmd = Form("mkdir -p ROOTfiles/%s", ftype.Data());
  }
  else{
    cmd = Form("mkdir -p ROOTfiles/%s", ftype.Data());
  }


  gSystem->Exec(cmd); // create study type dir. if it doesn't exist

  if((ftype=="shms50k") || (ftype=="hms50k")){
    
    cmd = Form("mkdir -p HISTOGRAMS/%s", ftype.Data());
    gSystem->Exec(cmd); // create study type dir. if it doesn't exist
    
    cmd = Form("mkdir -p HISTOGRAMS/%s/PDF", ftype.Data());
    gSystem->Exec(cmd); // create study type dir. if it doesn't exist
    
    cmd = Form("mkdir -p HISTOGRAMS/%s/ROOT", ftype.Data());
    gSystem->Exec(cmd); // create study type dir. if it doesn't exist
  }

  
  const char* ROOTFileNamePattern = "ROOTfiles/%s/deut_replay_%s_%d_%d.root";

  
  // Load global parameters
  gHcParms->Define("gen_run_number", "Run Number", RunNumber);
  gHcParms->AddString("g_ctp_database_filename", "DBASE/COIN/standard.database");
  gHcParms->Load(gHcParms->GetString("g_ctp_database_filename"), RunNumber);  // load the standard.database
  gHcParms->Load(gHcParms->GetString("g_ctp_parm_filename"));                 // load the general param
  gHcParms->Load(gHcParms->GetString("g_ctp_calib_filename"));                // load the detector calib param
  gHcParms->Load(gHcParms->GetString("g_ctp_cuts_filename"));                 // load the detector cuts param 
  gHcParms->Load(gHcParms->GetString("g_ctp_kinematics_filename"), RunNumber); // load the standard.kinematics file

  //Load params for coin. trigger configuration (up-to-date)
  //gHcParms->Load("PARAM/TRIG/tcoin.param"); 

  // Load the Hall C detector map
  gHcDetectorMap = new THcDetectorMap();

  gHcDetectorMap->Load("MAPS/COIN/DETEC/coin.map");     
  
  //use spring18 config det. map  / params
  //All runs before coin 4361 did NOT have hDCREF_5 added in ROC3
  if(RunNumber<=4361){
    
    //Load params for coin. trigger configuration (spring 2018, for testing using 2018 data) 
    //gHcParms->Load("PARAM/TRIG/archive/spring18/tcoin_spring18.param");
    
    // Load 2018 map
    gHcDetectorMap->Load("MAPS/COIN/DETEC/coin_comm18.map");
  }
  
  else if(RunNumber>=14933 && RunNumber<=15079){
    
    //Load map which has SHMS Cal ROC 04 SLot 9 offset by +1 channels (channel 0 empty)
    gHcDetectorMap->Load("MAPS/COIN/DETEC/coin_FADC_ROC4_Slot9_Offset.map");
  
  }
  
  //=:=:=:=
  // SHMS 
  //=:=:=:=
     // Dec data
  gHaApps->Add(new Podd::DecData("D","Decoder raw data"));
 
  // Set up the equipment to be analyzed.
  THcHallCSpectrometer* SHMS = new THcHallCSpectrometer("P", "SHMS");
  SHMS->SetEvtType(1);
  SHMS->AddEvtType(4);
  SHMS->AddEvtType(5);
  SHMS->AddEvtType(6);
  SHMS->AddEvtType(7);
  gHaApps->Add(SHMS);
  // Add Noble Gas Cherenkov to SHMS apparatus
  THcCherenkov* pngcer = new THcCherenkov("ngcer", "Noble Gas Cherenkov");
  SHMS->AddDetector(pngcer);
  // Add drift chambers to SHMS apparatus
  THcDC* pdc = new THcDC("dc", "Drift Chambers");
  SHMS->AddDetector(pdc);
  // Add hodoscope to SHMS apparatus
  THcHodoscope* phod = new THcHodoscope("hod", "Hodoscope");
  SHMS->AddDetector(phod);
  // Add Heavy Gas Cherenkov to SHMS apparatus
  THcCherenkov* phgcer = new THcCherenkov("hgcer", "Heavy Gas Cherenkov");
  SHMS->AddDetector(phgcer);
  // Add Aerogel Cherenkov to SHMS apparatus
  THcAerogel* paero = new THcAerogel("aero", "Aerogel");
  SHMS->AddDetector(paero);
  // Add calorimeter to SHMS apparatus
  THcShower* pcal = new THcShower("cal", "Calorimeter");
  SHMS->AddDetector(pcal);

  // Add rastered beam apparatus
  THaApparatus* pbeam = new THcRasteredBeam("P.rb", "Rastered Beamline");
  gHaApps->Add(pbeam);
  // Add physics modules
  // Calculate reaction point
  THcReactionPoint* prp = new THcReactionPoint("P.react", "SHMS reaction point", "P", "P.rb");
  gHaPhysics->Add(prp);
  // Calculate extended target corrections
  THcExtTarCor* pext = new THcExtTarCor("P.extcor", "HMS extended target corrections", "P", "P.react");
  gHaPhysics->Add(pext);
  // Calculate golden track quantites
  THaGoldenTrack* pgtr = new THaGoldenTrack("P.gtr", "SHMS Golden Track", "P");
  gHaPhysics->Add(pgtr);
  // Calculate the hodoscope efficiencies
  THcHodoEff* peff = new THcHodoEff("phodeff", "SHMS hodo efficiency", "P.hod");
  gHaPhysics->Add(peff);  

  // Add event handler for scaler events
  THcScalerEvtHandler* pscaler = new THcScalerEvtHandler("P", "Hall C scaler event type 1");
  pscaler->AddEvtType(1);
  pscaler->AddEvtType(4);
  pscaler->AddEvtType(5);
  pscaler->AddEvtType(6);
  pscaler->AddEvtType(7);
  pscaler->AddEvtType(129);
  pscaler->SetDelayedType(129);
  pscaler->SetUseFirstEvent(kTRUE);
  gHaEvtHandlers->Add(pscaler);

  /*
  //Add SHMS event handler for helicity scalers
  THcHelicityScaler *phelscaler = new THcHelicityScaler("P", "Hall C helicity scaler");
  //phelscaler->SetDebugFile("PHelScaler.txt");
  phelscaler->SetROC(8);   
  phelscaler->SetUseFirstEvent(kTRUE); 
  gHaEvtHandlers->Add(phelscaler); 
  */
  
  //=:=:=
  // HMS 
  //=:=:=

  // Set up the equipment to be analyzed.
  THcHallCSpectrometer* HMS = new THcHallCSpectrometer("H", "HMS");
  HMS->SetEvtType(2);
  HMS->AddEvtType(4);
  HMS->AddEvtType(5);
  HMS->AddEvtType(6);
  HMS->AddEvtType(7);
  gHaApps->Add(HMS);
  // Add drift chambers to HMS apparatus
  THcDC* hdc = new THcDC("dc", "Drift Chambers");
  HMS->AddDetector(hdc);
  // Add hodoscope to HMS apparatus
  THcHodoscope* hhod = new THcHodoscope("hod", "Hodoscope");
  HMS->AddDetector(hhod);
  // Add Cherenkov to HMS apparatus
  THcCherenkov* hcer = new THcCherenkov("cer", "Heavy Gas Cherenkov");
  HMS->AddDetector(hcer);
  // Add Aerogel Cherenkov to HMS apparatus
  // THcAerogel* haero = new THcAerogel("aero", "Aerogel");
  // HMS->AddDetector(haero);
  // Add calorimeter to HMS apparatus
  THcShower* hcal = new THcShower("cal", "Calorimeter");
  HMS->AddDetector(hcal);

  // Add rastered beam apparatus
  THaApparatus* hbeam = new THcRasteredBeam("H.rb", "Rastered Beamline");
  gHaApps->Add(hbeam);  
  // Add physics modules
  // Calculate reaction point
  THaReactionPoint* hrp = new THaReactionPoint("H.react", "HMS reaction point", "H", "H.rb");
  gHaPhysics->Add(hrp);
  // Calculate extended target corrections
  THcExtTarCor* hext = new THcExtTarCor("H.extcor", "HMS extended target corrections", "H", "H.react");
  gHaPhysics->Add(hext);
  // Calculate golden track quantities
  THaGoldenTrack* hgtr = new THaGoldenTrack("H.gtr", "HMS Golden Track", "H");
  gHaPhysics->Add(hgtr);
  // Calculate the hodoscope efficiencies
  THcHodoEff* heff = new THcHodoEff("hhodeff", "HMS hodo efficiency", "H.hod");
  gHaPhysics->Add(heff);

  // Add event handler for scaler events
  THcScalerEvtHandler *hscaler = new THcScalerEvtHandler("H", "Hall C scaler event type 4");  
  hscaler->AddEvtType(2);
  hscaler->AddEvtType(4);
  hscaler->AddEvtType(5);
  hscaler->AddEvtType(6);
  hscaler->AddEvtType(7);
  hscaler->AddEvtType(129);
  hscaler->SetDelayedType(129);
  hscaler->SetUseFirstEvent(kTRUE);
  gHaEvtHandlers->Add(hscaler);

  /*
  // Add HMS event handler for helicity scalers                                                                                             
  THcHelicityScaler *hhelscaler = new THcHelicityScaler("H", "Hall C helicity scaler"); 
  //hhelscaler->SetDebugFile("HHelScaler.txt");                                                                
  hhelscaler->SetROC(5);
  hhelscaler->SetUseFirstEvent(kTRUE); 
  gHaEvtHandlers->Add(hhelscaler);
  */
  
  //=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=
  // Kinematics Modules
  //=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=

  // ---------------------------------
  // electrons in SHMS, protons in HMS
  // ---------------------------------
  // Add physics module to calculate primary (scattered electrons) beam kinematics
  THcPrimaryKine* pkin_primary = new THcPrimaryKine("P.kin.primary", "SHMS Single Arm Kinematics", "P", "P.rb");
  gHaPhysics->Add(pkin_primary);
  // Add physics module to calculate secondary (scattered hadrons) beam kinematics
  THcSecondaryKine* hkin_secondary = new THcSecondaryKine("H.kin.secondary", "HMS Single Arm Kinematics", "H", "P.kin.primary");
  gHaPhysics->Add(hkin_secondary);
 
  //=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=
  // Global Objects & Event Handlers
  //=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=

  // Add trigger apparatus
  THaApparatus* TRG = new THcTrigApp("T", "TRG");
  gHaApps->Add(TRG);
  // Add trigger detector to trigger apparatus
  THcTrigDet* coin = new THcTrigDet("coin", "Coincidence Trigger Information");

  //Add coin physics module
  THcCoinTime* coinTime = new THcCoinTime("CTime", "Coincidende Time Determination", "H", "P", "T.coin");
  gHaPhysics->Add(coinTime);

  // Suppress missing reference time warnings for these event types
  coin->SetEvtType(1);
  coin->AddEvtType(2);
  TRG->AddDetector(coin); 

  //Add RF physics module THcRFTime::THcRFTime (const char *name, const char* description, const char* hadArmName, 
  // const char* elecArmName, const char* RFname) :
  THcRFTime* RFTime = new THcRFTime("RFTime", "RF Time Determination", "P", "H", "T.coin");
  gHaPhysics->Add(RFTime);

  
  // Add event handler for prestart event 125.
  THcConfigEvtHandler* ev125 = new THcConfigEvtHandler("HC", "Config Event type 125");
  gHaEvtHandlers->Add(ev125);
  // Add event handler for EPICS events
  THaEpicsEvtHandler* hcepics = new THaEpicsEvtHandler("epics", "HC EPICS event type 181");
  gHaEvtHandlers->Add(hcepics);
 
  // Set up the analyzer - we use the standard one,
  // but this could be an experiment-specific one as well.
  // The Analyzer controls the reading of the data, executes
  // tests/cuts, loops over Acpparatus's and PhysicsModules,
  // and executes the output routines.
  THcAnalyzer* analyzer = new THcAnalyzer;

  // A simple event class to be output to the resulting tree.
  // Creating your own descendant of THaEvent is one way of
  // defining and controlling the output.
  THaEvent* event = new THaEvent;

  // Define the run(s) that we want to analyze.
  // We just set up one, but this could be many.
  THcRun* run = new THcRun( pathList, Form(RunFileNamePattern, RunNumber) );
  run->Init();  // initialize run (necessary to get start_of_run) NOTE: run is a TDatime ROOT cern class reference, so TDatime methods can be invoked
  
  // Set to read in Hall C run database parameters
  run->SetRunParamClass("THcRunParameters");
  
  // Eventually need to learn to skip over, or properly analyze the pedestal events
  run->SetEventRange(1, MaxEvent); // Physics Event number, does not include scaler or control events.
  run->SetNscan(1);
  run->SetDataRequired(0x7);
  run->Print();

  // Define the analysis parameters
  TString ROOTFileName = Form(ROOTFileNamePattern, ftype.Data(), ftype.Data(), RunNumber, MaxEvent);


  analyzer->SetCountMode(2);  // 0 = counter is # of physics triggers
                              // 1 = counter is # of all decode reads
                              // 2 = counter is event number

  analyzer->SetEvent(event);

  // original (epics180), does not exits, so beam positions force to (0,0)
  //analyzer->SetEpicsEvtType(180);

  // Set EPICS event type (Feb 24, 2023 added these, suggested by M. Jones)
  // TO get non-zero beam position
  analyzer->SetEpicsEvtType(181);
  analyzer->AddEpicsEvtType(182);  

  // Define crate map
  analyzer->SetCrateMapFileName("MAPS/db_cratemap.dat");

  // Define output ROOT file
  analyzer->SetOutFile(ROOTFileName.Data());

  // Define DEF-file+
  TString DefTreeFile=Form("DEF-files/deut_%s.def",ftype.Data());
  analyzer->SetOdefFile(DefTreeFile);

  // Define cuts file
  TString DefCutTreeFile;

  if(ftype.Data()=="shms50k"){
    DefCutTreeFile="DEF-files/CUTS/deut_cuts_shms.def";
  }
  else if(ftype.Data()=="hms50k"){
    DefCutTreeFile="DEF-files/CUTS/deut_cuts_hms.def";
  }
  else{
    DefCutTreeFile="DEF-files/CUTS/deut_cuts.def";
  }
  
  analyzer->SetCutFile(DefCutTreeFile);  // optional

  
  // File to record accounting information for cuts
  cmd = Form("mkdir -p REPORT_OUTPUT/%s", ftype.Data());  
  gSystem->Exec(cmd); // create study type dir. if it doesn't exist
  analyzer->SetSummaryFile(Form("REPORT_OUTPUT/%s/summary_deut_%s_%d_%d.report", ftype.Data(), ftype.Data(), RunNumber, MaxEvent));  // optional

  // Start the actual analysis.
  analyzer->Process(run);

  // Create report file from template
  // C.Y. (for now we just have 1 template file, but this can be expanded to: deut_ftype.template (e.g. deut_heep.template, etc)
  TString REPORT_FileName=Form("REPORT_OUTPUT/%s/deut_%s_%d_%d.report", ftype.Data(), ftype.Data(), RunNumber, MaxEvent);
  TString TEMPLATE_FileName="TEMPLATES/deut_prod.template";

  //testing (to use deut template file to mimic what the shift crew would see during Deut data-taking, commnent the 'RunNumber<3400' below)
  //if(RunNumber<=3400){
  //  REPORT_FileName=Form("REPORT_OUTPUT/%s/deut_%s_spring18.report", ftype.Data(), ftype.Data() );
  //  TEMPLATE_FileName="TEMPLATES/deut_spring18.template";
  //}
  
  analyzer->PrintReport( TEMPLATE_FileName, REPORT_FileName );  // optional


  //----------------------------------------
  
  // C.Y. Jun 25, 2022  Added functionality to get start_of_run time, run_length, and hence calculate end_of_run time)
  // This calculation reads the REPORT_FILE output and hence assumes certain syntax, which MUST NOT be changed for this to work

  //  Write timestamp to REPORT_FILE to keep track of run start_time 
  ofstream file_object;
  file_object.open(REPORT_FileName.Data(), std::ios_base::app); 
  file_object << "" << endl;
  file_object << "start_of_run = " << (run->GetDate()).AsSQLString() << endl; // written in format: yyyy-mm-dd HH:MM:SS )
 
  
  // get the year, month day, hour, min sec of the start_of_run
  int year = (run->GetDate()).GetYear();
  int month = (run->GetDate()).GetMonth();
  int day = (run->GetDate()).GetDay();
  int hour = (run->GetDate()).GetHour();
  int minute = (run->GetDate()).GetMinute();
  int second = (run->GetDate()).GetSecond();

  // set TTimeStamp for start_of_run
  TTimeStamp *start_of_run_timestamp = new TTimeStamp(year, month, day, hour, minute, second);

  // get the run length from the REPORT FILE (NOTE this is the run length corresponding to the number of evts replayed.
  // if the Full Replay was done (all events),  end_of_run time will correspond to the HCLOG end_of_run time, otherwise
  // it will correspond to the run lenght of the number of events replayed)
  double run_len = stod(split(FindString("SHMS_Run_Length_sec", REPORT_FileName.Data())[0], ':')[1]); // in sec

  // add the run length in sec.
  start_of_run_timestamp->Add(run_len);

  // convert output to string in format: yyyy-mm-dd HH:MM:SS, with an implied UTC
  string end_of_run_time = start_of_run_timestamp->AsString("s");

  file_object << "end_of_run = " << end_of_run_time.c_str() << endl;
  file_object.close(); 
  

    
  // Get ending timepoint
  auto stop = high_resolution_clock::now();
  
  // Get duration. Substart timepoints to
  // get duration. To cast it to proper unit
  // use duration cast method
  auto duration = duration_cast<seconds>(stop - start);
  
  cout << "Time taken by replay_deut: "
       << duration.count() << " sec." << endl;


 stop:
  cout << "Exiting Now . . ." << endl;

}
