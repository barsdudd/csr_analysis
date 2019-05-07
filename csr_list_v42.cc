#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>

#include <map>

#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include "Event.h"
#include "Variables.h"
#include "Selection.h"

#include "Fitter2DPol1.h"
#include "Fitter2DPol2.h"
#include "FitterCommonPol2.h"
#include "Fitter2DCos.h"


#define FILE_EXIST 1
#define DIR_EXIST 2
#define NO_FILE_DIR_EXIST 0

vector<string> fileList;

using namespace std;

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
/////////////////////
//    VARIABLES    //
/////////////////////
///// INPUT /////
TFile*   file;
TTree*   tree;
Event*  event;

TFile*   saveFile;

double   spillID;
double targetPos;
double    RF[17];
double     G2SEM;
double    QIEsum;
double inh_thres;
double     runID;

vector<int>rs;
vector<int>rs_cont;
int        rs_current;
int        rs_current_cont;
string outDir;
string list;
///// DATA /////
TFile* saveCSRFile;
TTree* tree_for_dbub;
int    b_tgt;
int    b_ix2;
double b_mass;
double b_x1;
double b_x2;
double b_weight;
vector<TH1D*>         vec_th1_save;
vector<TGraphErrors*> vec_tge_save;
vector<TGraphAsymmErrors*> vec_csr_save;

/////// x2 ////////
TH1D* h1_rf_x2[4][8]; // [target type (dummy, LH2, Empty, LD2)][x2 bins]
TH1D* h1_avg_x2;
TH1D* h1_raw_x2;
TH1D* h1_avg_inte_x2[8]; // [x2]
TGraphErrors* gr_rf_ratio_x2[8]; // [x2 bins]
/////// x2 ////////

/////// x2 syst ////////
TH1D* h1_rf_x2_ped[2][4][8]; // systematics for pedestal: [low/high][target type (dummy, LH2, Empty, LD2)][x2 bins]
TH1D* h1_avg_inte_x2_ped[2][8]; // [low/high][x2]
TGraphErrors* gr_rf_ratio_x2_ped[2][8]; // [low/high][x2 bins]
/////// x2 syst ////////

/////// xF ////////
TH1D* h1_rf_xF[4][20]; // [target type (dummy, LH2, Empty, LD2)][x2 bins]
TH1D* h1_avg_xF;
TH1D* h1_raw_xF;
TH1D* h1_avg_inte_xF[20];
TGraphErrors* gr_rf_ratio_xF[20];
/////// xF ////////

/////// mass ////////
TH1D* h1_rf_mass[4][20]; // [target type (dummy, LH2, Empty, LD2)][x2 bins]
TH1D* h1_avg_mass;
TH1D* h1_raw_mass;
TH1D* h1_avg_inte_mass[20];
TGraphErrors* gr_rf_ratio_mass[20];
/////// mass ////////

/////// x1 ////////
TH1D* h1_rf_x1[4][20]; // [target type (dummy, LH2, Empty, LD2)][x2 bins]
TH1D* h1_avg_x1;
TH1D* h1_raw_x1;
TH1D* h1_avg_inte_x1[20];
TGraphErrors* gr_rf_ratio_x1[20];
/////// x1 ////////

/////// pT ////////
TH1D* h1_rf_pT[4][20]; // [target type (dummy, LH2, Empty, LD2)][x2 bins]
TH1D* h1_avg_pT;
TH1D* h1_raw_pT;
TH1D* h1_avg_inte_pT[20];
TGraphErrors* gr_rf_ratio_pT[20];
/////// pT ////////


TH1D* h1_raw    [8]; // [x2]


vector<double> nDimuons;
double nDimuonsRS;
bool potMode;

///// ADJESTABLE /////
const char* saveType = "png";

int binRange;
int    nbins;

int    nBinsrf;
double   maxrf;
double  fitMax;

double    cutMass;
bool    looseCuts;
bool    tightCuts;
bool      recCuts;
bool beamDumpCuts;

int bugFinder = 0;

bool isValidNumbers = true;

bool fittingMode = false;

/////////////////////
double THH;
double THD;
double TDD;

double  AH;
double  AD;
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////

//////////////////////
//    INITIALIZE    //
//////////////////////
int fileE(const char* filename){
   struct stat st;
   if(stat(filename, &st) != 0){
      return NO_FILE_DIR_EXIST;
   }else{
      mode_t m = st.st_mode;
      if(S_ISDIR(m)){
         return NO_FILE_DIR_EXIST;
      }else{
         return FILE_EXIST;
      }
   }
}

void initInputRS(const char* list){
   ostringstream oss;
   int rs_tmp;
   ifstream ifs;
   ifs.open(list);
   while( ifs >> rs_tmp ){
      oss.str("");
      oss << fileDir << rs_tmp << "/" << rs_tmp << "_with_cuts.root";
      fileList.push_back(oss.str().c_str());

      rs.push_back(rs_tmp);

      rs_cont.push_back(rs_tmp);
      if( rs_tmp == 67 ) rs_cont.push_back(68);
   }
   return;
}

bool initInputFile(int index){
   ostringstream oss;
   
   cout << " preparing " << fileList[index] << "......" << endl;
   if( fileE(fileList[index].c_str()) == NO_FILE_DIR_EXIST ){
      cout << fileList[index] << " DOES NOT EXIST. GO TO NEXT RUN..." << endl;
      return false;
   }
   file = new TFile(fileList[index].c_str(), "READ");
   tree = (TTree*)file->Get("save");

   tree->SetBranchAddress("event", &event);
   cout << "  Analyze " << fileList[index].c_str() << "..." << endl;
   nDimuonsRS = 0.0;
   return true;
}

void initTH1(TH1D*& h1, const char* name, const char* title, int nBin, double min, double max, int color, bool save=false){
   h1 = new TH1D(name, title, nBin, min, max);
   h1->Sumw2();
   h1->SetStats(0);
   h1->SetLineColor(color);
   h1->SetLineWidth(2);
   if( save )
      vec_th1_save.push_back(h1);
}

void initTH1(){
   ostringstream oss;
   ostringstream title;

   for( int ix = 0 ; ix < nBinsX2 ; ix++ ){
      for( int it = 1 ; it < 4 ; it++ ){
         oss.str("");
         oss << "h1_rf_x2" << "_" << it << "_" << ix;
         title.str("");
         title << "Counts (" << x2_bin[ix] << " < x_{T} < " << x2_bin[ix+1] 
               << ");Intensity;Counts";
         initTH1(h1_rf_x2[it][ix], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, colorList[it], true);
      }

      oss.str("");
      oss << "h1_avg_inte_x2_" << ix;
      title.str("");
      title << "Average (" << x2_bin[ix] << " < x_{2} < " << x2_bin[ix+1] 
            << ");Intensity;Counts";
      initTH1(h1_avg_inte_x2[ix], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, 2);
   }

   h1_avg_x2 = new TH1D("h1_avg_x2", "x_{2} average;x_{2};x_{2} average", nBinsX2, x2_bin);
   h1_avg_x2->Sumw2();
   h1_avg_x2->SetStats(0);
   h1_avg_x2->SetLineColor(2);
   h1_avg_x2->SetLineWidth(2);
   vec_th1_save.push_back(h1_avg_x2);
   h1_raw_x2     = new TH1D("h1_raw_x2"    , "x_{2} raw;x_{2};x_{2} raw"        , nBinsX2, x2_bin);
   h1_raw_x2->Sumw2();
   h1_raw_x2->SetStats(0);
   h1_raw_x2->SetLineColor(2);
   h1_raw_x2->SetLineWidth(2);

//////// x2 syst /////////
   for( int ip = 0 ; ip < 2 ; ip++ ){
      for( int ix = 0 ; ix < nBinsX2 ; ix++ ){
         for( int it = 1 ; it < 4 ; it++ ){
            oss.str("");
            oss << "h1_rf_x2_ped_" << ip << "_" << it << "_" << ix;
            title.str("");
            title << "Counts (" << x2_bin[ix] << " < x_{T} < " << x2_bin[ix+1] 
                  << ");Intensity;Counts";
            initTH1(h1_rf_x2_ped[ip][it][ix], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, colorList[it], true);
         }

         oss.str("");
         oss << "h1_avg_inte_x2_ped_" << ip << "_" << ix;
         title.str("");
         title << "Average (" << x2_bin[ix] << " < x_{2} < " << x2_bin[ix+1] 
               << ");Intensity;Counts";
         initTH1(h1_avg_inte_x2_ped[ip][ix], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, 2);
      }
   }
//////// x2 syst /////////

/////////////////////// mass ////////////////////////
   h1_avg_mass = new TH1D("h1_avg_mass", "mass average;mass;mass average", nBinsMass, mass_bin);
   h1_avg_mass->Sumw2();
   h1_avg_mass->SetStats(0);
   h1_avg_mass->SetLineColor(2);
   h1_avg_mass->SetLineWidth(2);
   vec_th1_save.push_back(h1_avg_mass);
   h1_raw_mass = new TH1D("h1_raw_mass"    , "mass raw;mass;mass raw"    , nBinsMass, mass_bin);
   h1_raw_mass->Sumw2();
   h1_raw_mass->SetStats(0);
   h1_raw_mass->SetLineColor(2);
   h1_raw_mass->SetLineWidth(2);
/////////////////////// mass ////////////////////////

/////////////////////// x1 ////////////////////////
   h1_avg_x1 = new TH1D("h1_avg_x1", "x1 average;x1;x1 average", nBinsX1, x1_bin);
   h1_avg_x1->Sumw2();
   h1_avg_x1->SetStats(0);
   h1_avg_x1->SetLineColor(2);
   h1_avg_x1->SetLineWidth(2);
   vec_th1_save.push_back(h1_avg_x1);
   h1_raw_x1 = new TH1D("h1_raw_x1"    , "x1 raw;x1;x1 raw"    , nBinsX1, x1_bin);
   h1_raw_x1->Sumw2();
   h1_raw_x1->SetStats(0);
   h1_raw_x1->SetLineColor(2);
   h1_raw_x1->SetLineWidth(2);
/////////////////////// x1 ////////////////////////

/////////////////////// xF ////////////////////////
   h1_avg_xF = new TH1D("h1_avg_xF", "xF average;xF;xF average", nBinsXF, xF_bin);
   h1_avg_xF->Sumw2();
   h1_avg_xF->SetStats(0);
   h1_avg_xF->SetLineColor(2);
   h1_avg_xF->SetLineWidth(2);
   vec_th1_save.push_back(h1_avg_xF);
   h1_raw_xF = new TH1D("h1_raw_xF"    , "xF raw;xF;xF raw"    , nBinsXF, xF_bin);
   h1_raw_xF->Sumw2();
   h1_raw_xF->SetStats(0);
   h1_raw_xF->SetLineColor(2);
   h1_raw_xF->SetLineWidth(2);
/////////////////////// xF ////////////////////////

/////////////////////// pT ////////////////////////
   h1_avg_pT = new TH1D("h1_avg_pT", "pT average;pT;pT average", nBinsPT, pT_bin);
   h1_avg_pT->Sumw2();
   h1_avg_pT->SetStats(0);
   h1_avg_pT->SetLineColor(2);
   h1_avg_pT->SetLineWidth(2);
   vec_th1_save.push_back(h1_avg_pT);
   h1_raw_pT = new TH1D("h1_raw_pT"    , "pT raw;pT;pT raw"    , nBinsPT, pT_bin);
   h1_raw_pT->Sumw2();
   h1_raw_pT->SetStats(0);
   h1_raw_pT->SetLineColor(2);
   h1_raw_pT->SetLineWidth(2);
/////////////////////// pT ////////////////////////

   for( int iMass = 0 ; iMass < nBinsMass ; iMass++ ){
      for( int it = 1 ; it < 4 ; it++ ){
         oss.str("");
         oss << "h1_rf_mass" << "_" << it << "_" << iMass;
         title.str("");
         title << "Counts (" << mass_bin[iMass] << " < Mass < " << mass_bin[iMass+1] 
               << ");Intensity;Counts";
         initTH1(h1_rf_mass[it][iMass], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, colorList[it], true);
      }
      oss.str("");
      oss << "h1_avg_mass_" << iMass;
      title.str("");
      title << "Counts (" << mass_bin[iMass] << " < Mass < " << mass_bin[iMass+1] 
            << ");Intensity;Avg. Intensity";
      initTH1(h1_avg_inte_mass[iMass], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, colorList[1], true);
   }
   for( int iXF = 0 ; iXF < nBinsXF ; iXF++ ){
      for( int it = 1 ; it < 4 ; it++ ){
         oss.str("");
         oss << "h1_rf_xF" << "_" << it << "_" << iXF;
         title.str("");
         title << "Counts (" << xF_bin[iXF] << " < XF < " << xF_bin[iXF+1] 
               << ");Intensity;Counts";
         initTH1(h1_rf_xF[it][iXF], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, colorList[it], true);
      }
      oss.str("");
      oss << "h1_avg_xF_" << iXF;
      title.str("");
      title << "Counts (" << xF_bin[iXF] << " < XF < " << xF_bin[iXF+1] 
            << ");Intensity;Avg. Intensity";
      initTH1(h1_avg_inte_xF[iXF], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, colorList[1], true);
   }
   for( int iX1 = 0 ; iX1 < nBinsX1 ; iX1++ ){
      for( int it = 1 ; it < 4 ; it++ ){
         oss.str("");
         oss << "h1_rf_x1" << "_" << it << "_" << iX1;
         title.str("");
         title << "Counts (" << x1_bin[iX1] << " < X1 < " << x1_bin[iX1+1] 
               << ");Intensity;Counts";
         initTH1(h1_rf_x1[it][iX1], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, colorList[it], true);
      }
      oss.str("");
      oss << "h1_avg_x1_" << iX1;
      title.str("");
      title << "Counts (" << x1_bin[iX1] << " < X1 < " << x1_bin[iX1+1] 
            << ");Intensity;Avg. Intensity";
      initTH1(h1_avg_inte_x1[iX1], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, colorList[1], true);
   }
   for( int iPT = 0 ; iPT < nBinsPT ; iPT++ ){
      for( int it = 1 ; it < 4 ; it++ ){
         oss.str("");
         oss << "h1_rf_pT" << "_" << it << "_" << iPT;
         title.str("");
         title << "Counts (" << pT_bin[iPT] << " < PT < " << pT_bin[iPT+1] 
               << ");Intensity;Counts";
         initTH1(h1_rf_pT[it][iPT], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, colorList[it], true);
      }
      oss.str("");
      oss << "h1_avg_pT_" << iPT;
      title.str("");
      title << "Counts (" << pT_bin[iPT] << " < PT < " << pT_bin[iPT+1] 
            << ");Intensity;Avg. Intensity";
      initTH1(h1_avg_inte_pT[iPT], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, colorList[1], true);
      
   }
}

void initOutput(bool mode){
   if( !mode ){
      saveFile = new TFile("results.root", "recreate");
      if( !saveFile->IsOpen() ){
         cout << "saveFile was not created" << endl;
         exit(0);
      }
   }
   saveCSRFile = new TFile("results_csr.root", "recreate");
   tree_for_dbub = new TTree("tr_for_dbub", "Event data for calc-dy-ken.");
   tree_for_dbub->Branch("tgt"   , &b_tgt   ,    "tgt/I");
   tree_for_dbub->Branch("ix2"   , &b_ix2   ,    "ix2/I");
   tree_for_dbub->Branch("mass"  , &b_mass  ,   "mass/D");
   tree_for_dbub->Branch("x1"    , &b_x1    ,     "x1/D");
   tree_for_dbub->Branch("x2"    , &b_x2    ,     "x2/D");
   tree_for_dbub->Branch("weight", &b_weight, "weight/D");
}

void initMain(int argc, char* argv[]){
   ostringstream oss;

   if( argc < 5 ){
      cout << "!! NEED at least 5 ARGUMENTS !!: ROADSET NBINS FITRANGE POTMODE[POT/DIMUON] CUTMODE[LC/RC/TC/BC] ONLYFITTING[FIT/ANY]" << endl;
      exit(0);
   }

   initInputRS   (argv[1]);
   nbins   = atoi(argv[2]);
   fitMax  = atof(argv[3]);
   potMode = strcmp(argv[4], "POT" ) == 0 ? true : false;

   binRange = fitMax;

   looseCuts    = strcmp(argv[5], "LC") == 0 ? true : false;
   recCuts      = strcmp(argv[5], "RC") == 0 ? true : false;
   tightCuts    = strcmp(argv[5], "TC") == 0 ? true : false;
   beamDumpCuts = strcmp(argv[5], "BC") == 0 ? true : false;

   int cutTypeIsValid = 0;
   if(    looseCuts ) cutTypeIsValid++;
   if(      recCuts ) cutTypeIsValid++;
   if(    tightCuts ) cutTypeIsValid++;
   if( beamDumpCuts ) cutTypeIsValid++;

   if ( cutTypeIsValid != 1 ){
      cout << "!! YOU NEED TO GIVE A \"CORRECT\" TYPE OF CUTS !! [LC, RC, TC, BC]" << endl;
      exit(1);
   }

   if( argc > 6 )
      if( strcmp(argv[6], "FIT") == 0 ) fittingMode = true;

   event = new Event();

   const char* potApp = potMode ? "ON" : "OFF";

   oss.str("");
   oss << "results_list_v42_" << argv[1];
   oss << "/nbins_" << nbins << "_range_" << fitMax << "_POTMODE_" << potApp << "_cuts_" << argv[5];// << tightApp;
   outDir = oss.str();

   gSystem->mkdir(outDir.c_str(), true);
   gSystem->cd   (outDir.c_str()      );

   initOutput( fittingMode );

   if( ! fittingMode ) initTH1();
}

///////////////////
//    EXTRACT    //
///////////////////
int findTrackIndex(int trackID){
   int index = -1;
   for( int it = 0 ; it < (int)event->tracks.size() ; it++ ){
      if( event->tracks[it].trackID == trackID ){
         index = it;
         break;
      }
   }
   if( index == -1 ){
      cout << "!! NO SUCH A TRACK !!" << endl;
      exit(0);
   }
   return index;
}

void extData(){

   for( int id = 0 ; id < (int)event->dimuons.size() ; id++ ){
      Dimuon dimuon = event->dimuons[id];
      Track  trackP = event->tracks[findTrackIndex(dimuon.trackID_pos)];
      Track  trackN = event->tracks[findTrackIndex(dimuon.trackID_neg)];

//      if(   dimuon.xF <= 0.1 ) continue;

      if( !   trackIsValid_2111_v42( trackP, rs[rs_current] ) ) continue;
      if( !   trackIsValid_2111_v42( trackN, rs[rs_current] ) ) continue;
      if( !  dimuonIsValid_2111_v42( dimuon, rs[rs_current], looseCuts ) ) continue;
      if( ! tracksAreValid_2111_v42( trackP, trackN, dimuon ) ) continue;

      if( tightCuts )
         if( ! tightMode_2111_v42(trackP, trackN, dimuon) ) continue;

      if( beamDumpCuts ) 
         if( ! beamDump_2111_v42(trackP, trackN, dimuon, rs[rs_current]) ) continue;

      if( event->targetPos == 3 ) nDimuonsRS++;
      int ix = getIX2(dimuon.x2);
      if( ix == -1 ) continue;
      if( event->targetPos > 3 ) continue;

//      event->inte_t[16] *= 1. - 0.0016;

      h1_rf_x2[event->targetPos][ix]->Fill(event->inte_t[16]);
      double inte_low  = (event->RF[16] - 34+4)*event->G2SEM/(event->QIEsum - turns*(34-4)*buckets);
      double inte_high = (event->RF[16] - 34-4)*event->G2SEM/(event->QIEsum - turns*(34+4)*buckets);
      h1_rf_x2_ped[0][event->targetPos][ix]->Fill(inte_low );
      h1_rf_x2_ped[1][event->targetPos][ix]->Fill(inte_high);

//////////////////////////////////////
      int iMass = getIMass(dimuon.mass);
      if( iMass != -1 ){
         h1_rf_mass      [event->targetPos][iMass]->Fill(event->inte_t[16]                   );
         h1_avg_inte_mass                  [iMass]->Fill(event->inte_t[16], event->inte_t[16]);
      }
      int iX1 = getIX1(dimuon.x1);
      if( iX1   != -1 ){
         h1_rf_x1      [event->targetPos][iX1]->Fill(event->inte_t[16]                   );
         h1_avg_inte_x1                  [iX1]->Fill(event->inte_t[16], event->inte_t[16]);
      }
      int iXF = getIXF(dimuon.xF);
      if( iXF   != -1 ){
         h1_rf_xF      [event->targetPos][iXF]->Fill(event->inte_t[16]                   );
         h1_avg_inte_xF                  [iXF]->Fill(event->inte_t[16], event->inte_t[16]);
      }
      double pT = sqrt( pow(dimuon.vtx_mom.X(), 2) + pow(dimuon.vtx_mom.Y(), 2) );
      int iPT = getIPT(pT);
      if( iPT   != -1 ){
         h1_rf_pT      [event->targetPos][iPT]->Fill(event->inte_t[16]                   );
         h1_avg_inte_pT                  [iPT]->Fill(event->inte_t[16], event->inte_t[16]);
      }
/////////////////////////////////////
      h1_avg_inte_x2[ix]->Fill(event->inte_t[16], event->inte_t[16]);

      h1_avg_inte_x2_ped[0][ix]->Fill(inte_low , inte_low );
      h1_avg_inte_x2_ped[1][ix]->Fill(inte_high, inte_high);

      if( fitMax > event->inte_t[16] ){
         if( event->targetPos != 2 ){
            h1_avg_x2->Fill(dimuon.x2, dimuon.x2);
            h1_raw_x2->Fill(dimuon.x2);

            if( iMass != -1 ){
               h1_avg_mass->Fill(dimuon.mass, dimuon.mass);
               h1_raw_mass->Fill(dimuon.mass);
            }
            if( iX1 != -1 ){
               h1_avg_x1->Fill(dimuon.x1, dimuon.x1);
               h1_raw_x1->Fill(dimuon.x1);
            }
            if( iXF != -1 ){
               h1_avg_xF->Fill(dimuon.xF, dimuon.xF);
               h1_raw_xF->Fill(dimuon.xF);
            }
            if( iPT != -1 ){
               h1_avg_pT->Fill(pT, pT);
               h1_raw_pT->Fill(pT);
            }
         }
         b_tgt    = event->targetPos;
         b_ix2    = ix;
         b_mass   = dimuon.mass;
         b_x1     = dimuon.x1;
         b_x2     = dimuon.x2;
         b_weight = 1;
         tree_for_dbub->Fill();
      }
   }

}

void extMain(){
   int prev = -1;
   int curr =  0;
   bool cont_change = false;
   for( int ie = 0 ; ie < tree->GetEntries() ; ++ie){    
      tree->GetEntry(ie);
      curr = (ie+1) * 100 / tree->GetEntries();
      if( curr != prev ){
         cout << "\r" << (ie+1) << " / " << tree->GetEntries() << " = " << curr << " %" << flush;
         prev = curr;
      }

      if( !occupancyIsValid_2111_v42(event) ) continue;
      if( event->runID >= 14653 && rs[rs_current] == 67 && !cont_change ){
         rs_current_cont++;
         nDimuons.push_back(nDimuonsRS);
         nDimuonsRS = 0;
         cout << "\n" << "RUNID: " << event->runID << " | CHANGE TO \"RS68\"" << endl;
         cont_change = true;
      }
      extData();

   }
   cout << endl;
   nDimuons.push_back(nDimuonsRS);
}

///////////////////
//    ANALYZE    //
///////////////////
void normHist(TH1D*& h1, int tgt){
   h1->Scale(1/getPoT(rs, tgt));
}

void getRatioAndError(double& ratio, double& error, double nume, double deno, double subt, double eNume, double eDeno, double eSubt){

   double eNumeCal = (THH * aMass_d / 2. / TDD / aMass_h) * 1. / AD / ( deno / AH - subt ) ;
   double eDenoCal = (THH * aMass_d / 2. / TDD / aMass_h) * ( nume / AD - subt) / pow( deno / AH - subt, 2 ) / AH ;
   double eSubtCal = (THH * aMass_d / 2. / TDD / aMass_h) * ( nume / AD - deno / AH ) / pow( deno / AH - subt, 2);

   ratio = (THH * aMass_d / 2. / TDD / aMass_h) * ( (nume / AD - subt) / ( deno / AH - subt ) - THD / THH ) ;   
   error = 
      sqrt(
         pow( eNume * eNumeCal, 2) +
         pow( eDeno * eDenoCal, 2) +
         pow( eSubt * eSubtCal, 2)
         )
      ;
   
}

TGraphErrors* getRatio(TH1D* h1_nume, TH1D* h1_deno, TH1D* h1_subt, TH1D* h1_ave, const char* gn){
   TGraphErrors* tge = new TGraphErrors();
   tge->SetName(gn);
   for( int ib = 1 ; ib <= h1_nume->GetNbinsX() ; ib++ ){
      double nume = h1_nume->GetBinContent(ib);
      double deno = h1_deno->GetBinContent(ib);
      double subt = h1_subt->GetBinContent(ib);
      
      if( deno - subt < 0 ) { cout << h1_nume->GetBinCenter(ib) << ": !! deno ( " << deno << " ) - subt ( " << subt << " ) < 0 !! " << endl; isValidNumbers = false;; }
      if( nume - subt < 0 ) { cout << h1_nume->GetBinCenter(ib) << ": !! nume ( " << nume << " ) - subt ( " << subt << " ) < 0 !! " << endl; isValidNumbers = false;; }
      double eNume = h1_nume->GetBinError(ib);
      double eDeno = h1_deno->GetBinError(ib);
      double eSubt = h1_subt->GetBinError(ib);

      double ratio, error;

      getRatioAndError(ratio, error, nume, deno, subt, eNume, eDeno, eSubt);

      tge->SetPoint     ( tge->GetN()  , h1_ave->GetBinContent(ib), ratio );
      tge->SetPointError( tge->GetN()-1,                         0, error );
      tge->SetMarkerColor(2);
      tge->SetMarkerSize(1.2);
      tge->SetMarkerStyle(20);

   }
   vec_tge_save.push_back(tge);
   return tge;
}

void anaMain(){
   ostringstream oss;

   TDD = targetLength * rho_d / getRDAverage(rs_cont, nDimuons, potMode) * (getFracDAverage(rs_cont, nDimuons, potMode)+getFracHDAverage(rs_cont, nDimuons, potMode)/2.);
   THD = targetLength * rho_d / getRDAverage(rs_cont, nDimuons, potMode) *                                              getFracHDAverage(rs_cont, nDimuons, potMode)/2. ;
   THH = targetLength * rho_h;

   lambda_d = getLambdaLD2TAverage(rs_cont, nDimuons, potMode);
   AH = lambda_h * ( 1 - exp( - targetLength / lambda_h ) )/targetLength;
   AD = lambda_d * ( 1 - exp( - targetLength / lambda_d ) )/targetLength;

   TH1D* h1_raw_temp;
/////// x2 /////////
   for( int ix = 0 ; ix < nBinsX2 ; ix++ ){
      for( int it = 1 ; it <= 3 ; it++ ){
         if( it == 1 ) h1_raw_temp = (TH1D*)h1_rf_x2[it][ix]->Clone();
         else          h1_raw_temp->Add(h1_rf_x2[it][ix]);
//         h1_rf_x2[it][ix]->Scale(1/getRawPoT(rs, it));
         h1_rf_x2[it][ix]->Scale(1/getPoT(rs, it));
      }

      h1_avg_inte_x2[ix]->Divide(h1_raw_temp);
      oss.str("");
      oss << "gr_rf_ratio_x2_" << ix;
      gr_rf_ratio_x2[ix] = getRatio(h1_rf_x2[3][ix], h1_rf_x2[1][ix], h1_rf_x2[2][ix], h1_avg_inte_x2[ix], oss.str().c_str());
   }
/////// x2 /////////

   h1_avg_x2    ->Divide(h1_raw_x2  );
   h1_avg_mass  ->Divide(h1_raw_mass);
   h1_avg_x1    ->Divide(h1_raw_x1  );
   h1_avg_xF    ->Divide(h1_raw_xF  );
   h1_avg_pT    ->Divide(h1_raw_pT  );

/////// x2 syst /////////
   for( int ip = 0 ; ip < 2 ; ip++ ){
      for( int ix = 0 ; ix < nBinsX2 ; ix++ ){
         for( int it = 1 ; it <= 3 ; it++ ){
            if( it == 1 ) h1_raw_temp = (TH1D*)h1_rf_x2_ped[ip][it][ix]->Clone();
            else          h1_raw_temp ->   Add(h1_rf_x2_ped[ip][it][ix]);
            h1_rf_x2_ped[ip][it][ix]->Scale(1/getRawPoT(rs, it));
         }
         h1_avg_inte_x2_ped[ip][ix]->Divide(h1_raw_temp);
         oss.str("");
         oss << "gr_rf_ratio_x2_ped_" << ip << "_" << ix;
         gr_rf_ratio_x2_ped[ip][ix] = getRatio(h1_rf_x2_ped[ip][3][ix], h1_rf_x2_ped[ip][1][ix], h1_rf_x2_ped[ip][2][ix], h1_avg_inte_x2_ped[ip][ix], oss.str().c_str());
      }
   }
/////// x2 syst /////////


//////////////////////////////////////////////////
   for( int iMass = 0 ; iMass < nBinsMass ; iMass++ ){
      for( int it = 1 ; it <= 3 ; it++ ){
         if( it == 1 ) h1_raw_temp = (TH1D*)h1_rf_mass[it][iMass]->Clone();
         else          h1_raw_temp->Add(h1_rf_mass[it][iMass]);
         h1_rf_mass[it][iMass]->Scale(1/getRawPoT(rs, it));
      }
      h1_avg_inte_mass[iMass]->Divide(h1_raw_temp);
      oss.str("");
      oss << "gr_rf_ratio_mass_" << iMass;
      gr_rf_ratio_mass[iMass] = getRatio(h1_rf_mass[3][iMass], h1_rf_mass[1][iMass], h1_rf_mass[2][iMass], h1_avg_inte_mass[iMass], oss.str().c_str());
   }

   for( int iX1 = 0 ; iX1 < nBinsX1 ; iX1++ ){
      for( int it = 1 ; it <= 3 ; it++ ){
         if( it == 1 ) h1_raw_temp = (TH1D*)h1_rf_x1[it][iX1]->Clone();
         else          h1_raw_temp->Add(h1_rf_x1[it][iX1]);
         h1_rf_x1[it][iX1]->Scale(1/getRawPoT(rs, it));
      }
      h1_avg_inte_x1[iX1]->Divide(h1_raw_temp);
      oss.str("");
      oss << "gr_rf_ratio_x1_" << iX1;
      gr_rf_ratio_x1[iX1] = getRatio(h1_rf_x1[3][iX1], h1_rf_x1[1][iX1], h1_rf_x1[2][iX1], h1_avg_inte_x1[iX1], oss.str().c_str());
   }

   for( int iXF = 0 ; iXF < nBinsXF ; iXF++ ){
      for( int it = 1 ; it <= 3 ; it++ ){
         if( it == 1 ) h1_raw_temp = (TH1D*)h1_rf_xF[it][iXF]->Clone();
         else          h1_raw_temp->Add(h1_rf_xF[it][iXF]);
         h1_rf_xF[it][iXF]->Scale(1/getRawPoT(rs, it));
      }
      h1_avg_inte_xF[iXF]->Divide(h1_raw_temp);
      oss.str("");
      oss << "gr_rf_ratio_xF_" << iXF;
      gr_rf_ratio_xF[iXF] = getRatio(h1_rf_xF[3][iXF], h1_rf_xF[1][iXF], h1_rf_xF[2][iXF], h1_avg_inte_xF[iXF], oss.str().c_str());
   }

   for( int iPT = 0 ; iPT < nBinsPT ; iPT++ ){
      for( int it = 1 ; it <= 3 ; it++ ){
         if( it == 1 ) h1_raw_temp = (TH1D*)h1_rf_pT[it][iPT]->Clone();
         else          h1_raw_temp->Add(h1_rf_pT[it][iPT]);
         h1_rf_pT[it][iPT]->Scale(1/getRawPoT(rs, it));
      }
      h1_avg_inte_pT[iPT]->Divide(h1_raw_temp);
      oss.str("");
      oss << "gr_rf_ratio_pT_" << iPT;
      gr_rf_ratio_pT[iPT] = getRatio(h1_rf_pT[3][iPT], h1_rf_pT[1][iPT], h1_rf_pT[2][iPT], h1_avg_inte_pT[iPT], oss.str().c_str());
   }
//////////////////////////////////////////////////

}

////////////////
//    DRAW    //
////////////////
void writeAverage(){
   ofstream ofs;
   ofs.open("x2_average.txt");
   for( int ib = 1; ib <= h1_avg_x2->GetNbinsX() ; ib++ ){
      ofs << h1_avg_x2->GetBinContent(ib) << "\t";
   }
   ofs << endl;
   ofs.close();
}

void drawTwoTGraphErrors(TGraphErrors* tge_1, TGraphErrors* tge_2, double xmin, double xmax, double ymin, double ymax, const char* label_1, const char* label_2, const char* title, const char* name){
   ostringstream oss;
   TCanvas* c1 = new TCanvas("c1", "", 800, 600);
   c1->SetGrid();

   TH1* frame = c1->DrawFrame(xmin, ymin, xmax, ymax);
   frame->SetTitle(title);

   TGraphErrors* tmp_1 = (TGraphErrors*) tge_1->Clone();
   TGraphErrors* tmp_2 = (TGraphErrors*) tge_2->Clone();

   tmp_1->SetMarkerColor(2);
   tmp_1->SetLineColor  (2);
   tmp_1->SetMarkerStyle(20);
   tmp_1->SetMarkerSize(1.2);
   tmp_2->SetMarkerColor(4);
   tmp_2->SetLineColor  (4);
   tmp_2->SetMarkerStyle(21);
   tmp_2->SetMarkerSize(1.2);

   frame->Draw();
   tmp_1->Draw("PE SAME");
   tmp_2->Draw("PE SAME");

   TLegend* leg = new TLegend(0.1, 0.7, 0.3, 0.9);
   leg->AddEntry(tmp_1, label_1, "pl");
   leg->AddEntry(tmp_2, label_2, "pl");
   leg->Draw();

   oss.str("");
   oss << name << "." << saveType;
   c1->SaveAs(oss.str().c_str());
   delete c1;
}


void drawTGraphErrors(TGraphErrors* tge, const char* title, const char* name){
   ostringstream oss;

   TCanvas* c1 = new TCanvas("c1", "", 800, 600);
   c1->SetGrid();

   TH1* frame = c1->DrawFrame(0, 0, binRange, 2);
   frame->SetTitle(title);
   frame->Draw();

   tge->Draw("SAME P");

   oss.str("");
   oss << name << "." << saveType;  
   c1->SaveAs(oss.str().c_str());

   delete c1;
}

void drawTGraphErrors(TGraphErrors* tge, const char* title){
   ostringstream oss;

   TCanvas* c1 = new TCanvas("c1", "", 800, 600);
   c1->SetGrid();

   TH1* frame = c1->DrawFrame(0, 0, binRange, 2);
   frame->SetTitle(title);
   frame->Draw();

   tge->Draw("SAME P");

   oss.str("");
   oss << tge->GetName() << "." << saveType;  
   c1->SaveAs(oss.str().c_str());

   delete c1;
}

void fitTGraphErrors(TGraphErrors* tge, const char* title, 
                     const char* fitFunc, double* paras, double* errors,
                     double& chisq, int& ndf,
                     const char* name,
                     double fitMax = binRange){
   gStyle->SetOptFit(1111);
   ostringstream oss;

   TCanvas* c1 = new TCanvas("c1", "", 800, 600);
   c1->SetGrid();

   TH1* frame = c1->DrawFrame(0, 0, binRange, 5);
   frame->SetTitle(title);
   frame->Draw();

   for( int ip = 0 ; ip < tge->GetN() ; ip++){
      double x, y;
      tge->GetPoint(ip, x, y);
      cout << x << "\t" << y << endl;
   }

   tge->Draw("SAME P");

   TF1* f1 = new TF1("f1", fitFunc, 0, fitMax);
   f1->SetParameters(paras);
   for( int ip = 0 ; ip < f1->GetNpar() ; ip++ ){
      f1->SetParLimits(ip, 0, 10);
   }
   tge->Fit("f1", "", "", 0, fitMax);

   for( int ip = 0 ; ip < f1->GetNumberFreeParameters() ; ip++ ){
      paras [ip] = f1->GetParameter(ip);
      errors[ip] = f1->GetParError (ip);
   }
   chisq = f1->GetChisquare();
   ndf   = f1->GetNDF();

   oss.str("");
   oss << "fit_";
   if( binRange != fitMax ) oss << fitMax << "_";
   oss << name << "." << saveType;  
   c1->SaveAs(oss.str().c_str());

   gStyle->SetOptFit(0);
   delete c1;
}

void drawTwoTH1s(TH1D* h1_1, TH1D* h1_2, const char* label_1, const char* label_2, const char* title, const char* name){
   ostringstream oss;
   TCanvas* c1 = new TCanvas("c1", "", 800, 600);
   c1->SetGrid();

   h1_1->SetTitle(title);
   h1_1->Draw("PE");
   h1_2->Draw("PE SAME");

   TLegend* leg = new TLegend(0.1, 0.7, 0.3, 0.9);
   leg->AddEntry(h1_1, label_1, "pl");
   leg->AddEntry(h1_2, label_2, "pl");
   leg->Draw();

   oss.str("");
   oss << name << "." << saveType;
   c1->SaveAs(oss.str().c_str());
   delete c1;
}

void printRatio(TGraphErrors* tge, const char* name){
   ofstream ofs;
   ofs.open(name);
   for( int ip = 0 ; ip < tge->GetN() ; ip++ ){
      double x, y,ye;
      tge->GetPoint(ip, x, y);
      ye = tge->GetErrorY(ip);
      ofs << y << "\t" << ye << "\t" << x << endl;
   }
   ofs.close();
}

void printRatioTH1(TH1D* h1, const char* name){
   ofstream ofs;
   ofs.open(name);
   for( int ip = 1 ; ip <= h1->GetNbinsX() ; ip++ ){
      ofs << h1->GetBinCenter(ip) << "\t" << h1->GetBinContent(ip) << "\t" 
          << h1->GetBinError (ip) << endl;
   }
   ofs.close();
}

void saveHists(){
   saveFile->cd();
   for( int ih = 0 ; ih < (int)vec_th1_save.size() ; ih++ )
      vec_th1_save[ih]->Write();
   for( int ig = 0 ; ig < (int)vec_tge_save.size() ; ig++ )
      vec_tge_save[ig]->Write();
   saveFile->Close();
}

void saveCSR(){
   //TFile* saveCSRFile = new TFile("results_csr.root", "recreate");
   saveCSRFile->cd();
   for( int ic = 0 ; ic < (int)vec_csr_save.size() ; ic++ )
      vec_csr_save[ic]->Write();
   tree_for_dbub->Write();
   saveCSRFile->Close();
}

void drawComp(){
   ostringstream oss;
   TCanvas* c1 = new TCanvas("c1", "", 800, 600);
   c1->SetGrid();
   TLegend* leg;
   for( int ix = 0 ; ix < nBinsX2 ; ix++ ){
      leg = new TLegend(0.1, 0.7, 0.3, 0.9);
      double max = 0 ;
      for( int it = 1 ; it <= 3 ; it++ ) 
         max = max > h1_rf_x2[it][ix]->GetMaximum() ? max : h1_rf_x2[it][ix]->GetMaximum();
      for( int it = 1 ; it <= 3 ; it++ ) {
         if( it == 1 ){
            h1_rf_x2[it][ix]->GetYaxis()->SetRangeUser(0, max*1.05);
            h1_rf_x2[it][ix]->Draw("hist");
         }
         else
            h1_rf_x2[it][ix]->Draw("hist same");
         leg->AddEntry(h1_rf_x2[it][ix], targetType[it], "l");
      }
      leg->Draw();
      oss.str("");
      oss << "h1_rf_x2_comp_" << ix << "." << saveType;
      c1->SaveAs(oss.str().c_str());
   }
   delete c1;
   delete leg;
}

void drawCSR(){
   TCanvas* c1 = new TCanvas("c1", "", 800, 600);
   c1->SetGrid();
   TH1* frame;
   double margin = 0.05;
   for( int ic = 0 ; ic < (int)vec_csr_save.size() ; ic++ ){
      TGraphAsymmErrors* gn = vec_csr_save[ic];
      double xErr;
      double x, y;
      gn->GetPoint(0, x, y);
      xErr = gn->GetErrorXlow(0);
      double xLow = x - xErr - margin;
      gn->GetPoint(gn->GetN()-1, x, y);
      xErr = gn->GetErrorXhigh(gn->GetN()-1);
      double xHigh = x + xErr + margin;

      frame = c1->DrawFrame(xLow, 0.7, xHigh, 1.4);
      frame->SetTitle(gn->GetTitle());
      gn->Draw("PSAME");
      c1->Print( ((string)(gn->GetName())+"."+saveType).c_str() );
   }
   delete c1;
}

double calcSystSquare(TGraphAsymmErrors* center, TGraphAsymmErrors* syst, int ip){
   double xC, yC, xS, yS;
   center->GetPoint(ip, xC, yC);
   syst  ->GetPoint(ip, xS, yS);
   return pow(yC-yS, 2);
}

void drawMain(){
   ostringstream oss;
   ostringstream title;

   if( ! fittingMode ){
      writeAverage();
      drawComp();
      
      
      for( int ix = 0 ; ix < nBinsX2 ; ix++ ){
         oss.str("");
         oss << "ratio_intensity_sum2_" << ix << ".txt";
         printRatio(gr_rf_ratio_x2[ix], oss.str().c_str());
      }
            
      for( int ig = 0 ; ig < (int)vec_tge_save.size() ; ig++ ){
         drawTGraphErrors(vec_tge_save[ig], ";Intensity;CSR");
      }
      
      saveHists();
   }


///////////////////
///// FITTING /////
///////////////////

   cout << "\n///////////////////"   << endl;
   cout <<   "///// FITTING /////"   << endl;
   cout <<   "///////////////////\n" << endl;

   FitterMultiHist* fit;
   ostringstream oss_avg;
   vector<string> valSet;
   valSet.push_back("x2");
   valSet.push_back("mass");
   valSet.push_back("x1");
   valSet.push_back("xF");
   valSet.push_back("pT");
   TGraphAsymmErrors* cent_csr;
   TGraphAsymmErrors* syst_fit[2];

   for( int iv = 0 ; iv < (int)valSet.size() ; iv++ ){
      ////////////// POL2 //////////////
      fit = new Fitter2DPol2(); 
      oss.str(""); oss << "gr_rf_ratio_" << valSet[iv];
      oss_avg.str(""); oss_avg << "h1_avg_" << valSet[iv];
      fit->Init("results.root", oss.str().c_str(), oss_avg.str().c_str());
      fit->DoFit();
      oss.str(""); oss << "result_pol2_fit_" << valSet[iv];
      fit->PrintFit(oss.str().c_str());
      fit->PrintFit(cout);
      fit->DrawFit();
      oss.str(""); oss << "csr_pol2_as_func_" << valSet[iv];
      if( iv == 0 ){
         cent_csr = fit->getCSRPlot(oss.str().c_str(), valSet[iv], "#sigma_{pd}/2#sigma_{pd}");
         vec_csr_save.push_back(cent_csr);
      }
      else vec_csr_save.push_back(fit->getCSRPlot(oss.str().c_str(), valSet[iv], "#sigma_{pd}/2#sigma_{pd}"));
      ////////////// POL2 //////////////

      ////////////// POL1 //////////////
      fit = new Fitter2DPol1(); 
      oss.str(""); oss << "gr_rf_ratio_" << valSet[iv];
      oss_avg.str(""); oss_avg << "h1_avg_" << valSet[iv];
      fit->Init("results.root", oss.str().c_str(), oss_avg.str().c_str());
      fit->DoFit();
      oss.str(""); oss << "result_pol1_fit_" << valSet[iv];
      fit->PrintFit(oss.str().c_str());
      fit->PrintFit(cout);
      fit->DrawFit();
      oss.str(""); oss << "csr_pol1_as_func_" << valSet[iv];
      vec_csr_save.push_back(fit->getCSRPlot(oss.str().c_str(), valSet[iv], "#sigma_{pd}/2#sigma_{pd}"));
      ////////////// POL1 //////////////

      ////////////// COS //////////////
      fit = new Fitter2DCos(); 
      oss.str(""); oss << "gr_rf_ratio_" << valSet[iv];
      oss_avg.str(""); oss_avg << "h1_avg_" << valSet[iv];
      fit->Init("results.root", oss.str().c_str(), oss_avg.str().c_str());
      fit->DoFit();
      oss.str(""); oss << "result_cos_fit_" << valSet[iv];
      fit->PrintFit(oss.str().c_str());
      fit->PrintFit(cout);
      fit->DrawFit();
      oss.str(""); oss << "csr_cos_as_func_" << valSet[iv];
      if( iv == 0 ){
         syst_fit[0] = fit->getCSRPlot(oss.str().c_str(), valSet[iv], "#sigma_{pd}/2#sigma_{pd}");
         vec_csr_save.push_back(syst_fit[0]);
      }
      else vec_csr_save.push_back(fit->getCSRPlot(oss.str().c_str(), valSet[iv], "#sigma_{pd}/2#sigma_{pd}"));
      ////////////// COS //////////////

      ////////////// COMMON POL2 //////////////
      fit = new FitterCommonPol2(); 
      oss.str(""); oss << "gr_rf_ratio_" << valSet[iv];
      oss_avg.str(""); oss_avg << "h1_avg_" << valSet[iv];
      fit->Init("results.root", oss.str().c_str(), oss_avg.str().c_str());
      fit->DoFit();
      oss.str(""); oss << "result_common_pol2_fit_" << valSet[iv];
      fit->PrintFit(oss.str().c_str());
      fit->PrintFit(cout);
      fit->DrawFit();
      oss.str(""); oss << "csr_common_pol2_as_func_" << valSet[iv];
      if( iv == 0 ){
         syst_fit[1] = fit->getCSRPlot(oss.str().c_str(), valSet[iv], "#sigma_{pd}/2#sigma_{pd}");
         vec_csr_save.push_back(syst_fit[1]);
      }
      else vec_csr_save.push_back(fit->getCSRPlot(oss.str().c_str(), valSet[iv], "#sigma_{pd}/2#sigma_{pd}"));
      ////////////// COMMON POL2 //////////////
   }

/////// pedestal ////////
   TGraphAsymmErrors* syst_ped[2];
   for( int ip = 0 ; ip < 2 ; ip++){
      fit = new Fitter2DPol2(); 
      oss.str(""); oss << "gr_rf_ratio_x2_ped_" << ip;
      oss_avg.str(""); oss_avg << "h1_avg_x2";
      fit->Init("results.root", oss.str().c_str(), oss_avg.str().c_str());
      fit->DoFit();
      oss.str(""); oss << "result_pol2_fit_x2_ped_" << ip;
      fit->PrintFit(oss.str().c_str());
      fit->PrintFit(cout);
      fit->DrawFit();
      oss.str(""); oss << "csr_pol2_as_func_x2_ped_" << ip;
      syst_ped[ip] = fit->getCSRPlot(oss.str().c_str(), "x2", "#sigma_{pd}/2#sigma_{pd}");
   }
/////// pedestal ////////

   double syst  [20] = {0.};
   double x2_avg[20] = {0.};
   cout << setprecision(4) << fixed;
   for( int ix = 0 ; ix < nBinsX2 ; ix++ ){
      double temp_syst;
      temp_syst  = 0;
      temp_syst += calcSystSquare(cent_csr, syst_fit[0], ix);
      temp_syst += calcSystSquare(cent_csr, syst_fit[1], ix);
      syst[ix]  += temp_syst;
      cout << sqrt(temp_syst) << "\t";

      temp_syst = 0;
      temp_syst += calcSystSquare(cent_csr, syst_ped[0], ix);
      temp_syst += calcSystSquare(cent_csr, syst_ped[1], ix);
      temp_syst  = temp_syst / 2;
      cout << sqrt(temp_syst) << "\t";

      syst[ix]  += temp_syst;      

      syst[ix]  = sqrt(syst[ix]);
      double temp;
      cout << syst[ix] << endl;
      cent_csr->GetPoint(ix, x2_avg[ix], temp);
   }

   TGraphAsymmErrors* tge_syst = new TGraphAsymmErrors(nBinsX2, x2_avg, syst);
   tge_syst->SetName("csr_pol2_as_func_x2_sys");
   tge_syst->SetTitle(";x_{2};Systematics");
   tge_syst->SetMarkerSize(1.2);
   tge_syst->SetMarkerColor(4);
   tge_syst->SetLineColor(1);   
   tge_syst->SetLineWidth(2);
   tge_syst->SetMarkerStyle(20);

   vec_csr_save.push_back(tge_syst);

   drawCSR();
   saveCSR();
}

////////////////
//    MAIN    //
////////////////
void deleteMemory(){
   delete   file;
}

int main(int argc, char* argv[]){
   gROOT ->SetStyle("Plain");
   gStyle->SetOptStat(1);
   gStyle->SetPadGridX(kTRUE);
   gStyle->SetPadGridY(kTRUE);
   
   initMain(argc, argv);

   if( !fittingMode ){
      rs_current_cont = 0;
      for( int ir = 0 ; ir < (int)fileList.size() ; ir++ ){
         if( !initInputFile (ir) ) continue;
         rs_current = ir;
         rs_current_cont++;
         extMain ();
         deleteMemory();
      }
      anaMain ();
   }

   drawMain();      

   if( ! fittingMode ) 
      cout << "TDD = " << TDD << " | THD = " << THD << " | THH = " << THH << " / AD = " << AD << " | AH = " << AH 
           << " / RD = " << getRDAverage(rs_cont, nDimuons, potMode) << " | lambda_d = " << lambda_d 
           << " | C = "  << getFracHDAverage(rs_cont, nDimuons, potMode) << endl;
   cout << "OUTPUT DIR: " << gSystem->pwd() << endl;
}
