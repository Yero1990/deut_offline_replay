//Main Calibration Code
#include "DC_calib.h"
#include "DC_calib.C"
#include <iostream>
#include <ctime>
using namespace std;

int main_calib()
{

  //prevent root from displaying graphs while executing
  gROOT->SetBatch(1);


  //measure execution time
  clock_t cl;
  cl = clock();
  
  //template arguments
  //DC_calib obj("spec", "path/to/rootfile.root", runNUM, eventNUm, "pid_flag", "calib_mode"); pid_flag: "pid_elec" or "pid_kFALSE", calib_mode: "wire" or "card"
               
  DC_calib obj("HMS", "ROOTfiles/dc_calib/cafe_replay_dc_UnCalib_16962_-1.root ", 16962, -1, "pid_kFALSE", "card");
  
  obj.setup_Directory();
  obj.SetPlaneNames();
  obj.GetDCLeafs();
  obj.AllocateDynamicArrays();
  obj.SetTdcOffset();
  obj.CreateHistoNames();
  obj.EventLoop("FillUncorrectedTimes");
  obj.Calculate_tZero();
  obj.EventLoop("ApplyT0Correction");
  obj.WriteTZeroParam();
  obj.WriteLookUpTable();
  obj.WriteToFile(1);  //set argument to (1) for debugging
 

  //stop clock
 cl = clock() - cl;
 cout << "execution time: " << cl/(double)CLOCKS_PER_SEC << " sec" << endl;

  return 0;
}
