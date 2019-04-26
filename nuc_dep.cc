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
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
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

string variableType;
int iInteBin;
int iValBin;

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
vector<TH1D*>         vec_th1_save;
vector<TGraphErrors*> vec_tge_save;
vector<TGraphAsymmErrors*> vec_csr_save;

/////// x2 ////////
TH1D* h1_rf_x2[8][8]; // [target type (dummy, LH2, Empty, LD2, None, Fe, C, W)][x2 bins]
TH1D* h1_avg_x2[8];
TH1D* h1_raw_x2[8];
TH1D* h1_avg_inte_x2[8][8]; // [x2]
TGraphErrors* gr_rf_ratio_x2[8][8]; // [x2 bins]
/////// x2 ////////

/////// xF ////////
TH1D* h1_rf_xF[8][20]; // [target type (dummy, LH2, Empty, LD2, None, Fe, C, W)][x2 bins]
TH1D* h1_avg_xF[8];
TH1D* h1_raw_xF[8];
TH1D* h1_avg_inte_xF[8][20];
TGraphErrors* gr_rf_ratio_xF[8][20];
/////// xF ////////

/////// mass ////////
TH1D* h1_rf_mass[8][20]; // [target type (dummy, LH2, Empty, LD2, None, Fe, C, W)][mass bins]
TH1D* h1_avg_mass[8];
TH1D* h1_raw_mass[8];
TH1D* h1_avg_inte_mass[8][20];
TGraphErrors* gr_rf_ratio_mass[8][20];
/////// mass ////////

/////// x1 ////////
TH1D* h1_rf_x1[8][20]; // [target type (dummy, LH2, Empty, LD2, None, Fe, C, W)][x1 bins]
TH1D* h1_avg_x1[8];
TH1D* h1_raw_x1[8];
TH1D* h1_avg_inte_x1[8][20];
TGraphErrors* gr_rf_ratio_x1[8][20];
/////// x1 ////////

/////// pT ////////
TH1D* h1_rf_pT[8][20]; // [target type (dummy, LH2, Empty, LD2, None, Fe, C, W)][pT bins]
TH1D* h1_avg_pT[8];
TH1D* h1_raw_pT[8];
TH1D* h1_avg_inte_pT[8][20];
TGraphErrors* gr_rf_ratio_pT[8][20];

TH1D* h1_pT_dist[8];
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
double T[8];
double A[8];


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


   for( int it = 1 ; it < 8 ; it++ ){
      for( int ix = 0 ; ix < 8 ; ix++ ){
         oss.str("");
         oss << "h1_rf_x2" << "_" << it << "_" << ix;
         title.str("");
         title << "Counts (" << x2_bin[ix] << " < x_{2} < " << x2_bin[ix+1] 
               << ");Intensity;Counts";
         initTH1(h1_rf_x2[it][ix], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, colorList[it], true);
         
	 oss.str("");
	 oss << "h1_avg_inte_x2_" << it << "_" << ix;
	 title.str("");
	 title << "Average (" << x2_bin[ix] << " < x_{2} < " << x2_bin[ix+1] 
	       << ");Intensity;Counts";
	 initTH1(h1_avg_inte_x2[it][ix], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, 2);
	 
      }
     
      oss.str("");
      oss << "h1_avg_x2_" << it;
      h1_avg_x2[it] = new TH1D(oss.str().c_str(), "x_{2} average;x_{2};x_{2} average", nBinsX2, x2_bin);
      h1_avg_x2[it]->Sumw2();
      h1_avg_x2[it]->SetStats(0);
      h1_avg_x2[it]->SetLineColor(2);
      h1_avg_x2[it]->SetLineWidth(2);
      vec_th1_save.push_back(h1_avg_x2[it]);
      oss.str("");
      oss << "h1_raw_x2_" << it;
      h1_raw_x2[it] = new TH1D(oss.str().c_str(), "x_{2} raw;x_{2};x_{2} raw"        , nBinsX2, x2_bin);
      h1_raw_x2[it]->Sumw2();
      h1_raw_x2[it]->SetStats(0);
      h1_raw_x2[it]->SetLineColor(2);
      h1_raw_x2[it]->SetLineWidth(2);
   
/////////////////////// mass ////////////////////////
      for( int iMass = 0 ; iMass < nBinsMass ; iMass++ ){
         oss.str("");
         oss << "h1_rf_mass" << "_" << it << "_" << iMass;
         title.str("");
         title << "Counts (" << mass_bin[iMass] << " < Mass < " << mass_bin[iMass+1] 
               << ");Intensity;Counts";
         initTH1(h1_rf_mass[it][iMass], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, colorList[it], true);
       
         oss.str("");
         oss << "h1_avg_mass_" << it << "_" << iMass;
         title.str("");
         title << "Counts (" << mass_bin[iMass] << " < Mass < " << mass_bin[iMass+1] 
               << ");Intensity;Avg. Intensity";
         initTH1(h1_avg_inte_mass[it][iMass], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, colorList[1], true);
      }
     
      oss.str("");
      oss << "h1_avg_mass_" << it;
      h1_avg_mass[it] = new TH1D(oss.str().c_str(), "mass average;mass;mass average", nBinsMass, mass_bin);
      h1_avg_mass[it]->Sumw2();
      h1_avg_mass[it]->SetStats(0);
      h1_avg_mass[it]->SetLineColor(2);
      h1_avg_mass[it]->SetLineWidth(2);
      vec_th1_save.push_back(h1_avg_mass[it]);
      oss.str("");
      oss << "h1_raw_mass_" << it;
      h1_raw_mass[it] = new TH1D(oss.str().c_str(), "mass raw;mass;mass raw"    , nBinsMass, mass_bin);
      h1_raw_mass[it]->Sumw2();
      h1_raw_mass[it]->SetStats(0);
      h1_raw_mass[it]->SetLineColor(2);
      h1_raw_mass[it]->SetLineWidth(2);
      /////////////////////// mass ////////////////////////
     
      /////////////////////// x1 ////////////////////////
      for( int iX1 = 0 ; iX1 < nBinsX1 ; iX1++ ){
         oss.str("");
         oss << "h1_rf_x1" << "_" << it << "_" << iX1;
         title.str("");
         title << "Counts (" << x1_bin[iX1] << " < X1 < " << x1_bin[iX1+1] 
               << ");Intensity;Counts";
         initTH1(h1_rf_x1[it][iX1], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, colorList[it], true);
       
         oss.str("");
         oss << "h1_avg_x1_" << it << "_" << iX1;
         title.str("");
         title << "Counts (" << x1_bin[iX1] << " < X1 < " << x1_bin[iX1+1] 
               << ");Intensity;Avg. Intensity";
         initTH1(h1_avg_inte_x1[it][iX1], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, colorList[1], true);
      }

      oss.str("");
      oss << "h1_avg_x1" << "_" << it;
      h1_avg_x1[it] = new TH1D(oss.str().c_str(), "x1 average;x1;x1 average", nBinsX1, x1_bin);
      h1_avg_x1[it]->Sumw2();
      h1_avg_x1[it]->SetStats(0);
      h1_avg_x1[it]->SetLineColor(2);
      h1_avg_x1[it]->SetLineWidth(2);
      vec_th1_save.push_back(h1_avg_x1[it]);
     
      oss.str("");
      oss << "h1_raw_x1" << "_" << it;
      h1_raw_x1[it] = new TH1D(oss.str().c_str(), "x1 raw;x1;x1 raw"    , nBinsX1, x1_bin);
      h1_raw_x1[it]->Sumw2();
      h1_raw_x1[it]->SetStats(0);
      h1_raw_x1[it]->SetLineColor(2);
      h1_raw_x1[it]->SetLineWidth(2);
      /////////////////////// x1 ////////////////////////

      /////////////////////// xF ////////////////////////
      for( int iXF = 0 ; iXF < nBinsXF ; iXF++ ){
         oss.str("");
         oss << "h1_rf_xF" << "_" << it << "_" << iXF;
         title.str("");
         title << "Counts (" << xF_bin[iXF] << " < XF < " << xF_bin[iXF+1] 
               << ");Intensity;Counts";
         initTH1(h1_rf_xF[it][iXF], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, colorList[it], true);
       
         oss.str("");
         oss << "h1_avg_xF_" << it << "_" << iXF;
         title.str("");
         title << "Counts (" << xF_bin[iXF] << " < XF < " << xF_bin[iXF+1] 
               << ");Intensity;Avg. Intensity";
         initTH1(h1_avg_inte_xF[it][iXF], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, colorList[1], true);
      }
     
      oss.str("");
      oss << "h1_avg_xF" << "_" << it;
      h1_avg_xF[it] = new TH1D(oss.str().c_str(), "xF average;xF;xF average", nBinsXF, xF_bin);
      h1_avg_xF[it]->Sumw2();
      h1_avg_xF[it]->SetStats(0);
      h1_avg_xF[it]->SetLineColor(2);
      h1_avg_xF[it]->SetLineWidth(2);
      vec_th1_save.push_back(h1_avg_xF[it]);

      oss.str("");
      oss << "h1_raw_xF" << "_" << it;
      h1_raw_xF[it] = new TH1D(oss.str().c_str(), "xF raw;xF;xF raw"    , nBinsXF, xF_bin);
      h1_raw_xF[it]->Sumw2();
      h1_raw_xF[it]->SetStats(0);
      h1_raw_xF[it]->SetLineColor(2);
      h1_raw_xF[it]->SetLineWidth(2);
/////////////////////// xF ////////////////////////
     
/////////////////////// pT ////////////////////////
      for( int iPT = 0 ; iPT < nBinsPT ; iPT++ ){
         oss.str("");
         oss << "h1_rf_pT" << "_" << it << "_" << iPT;
         title.str("");
         title << "Counts (" << pT_bin[iPT] << " < PT < " << pT_bin[iPT+1] 
               << ");Intensity;Counts";
         initTH1(h1_rf_pT[it][iPT], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, colorList[it], true);
     
         oss.str("");
         oss << "h1_avg_pT_" << it << "_" << iPT;
         title.str("");
         title << "Counts (" << pT_bin[iPT] << " < PT < " << pT_bin[iPT+1] 
               << ");Intensity;Avg. Intensity";
         initTH1(h1_avg_inte_pT[it][iPT], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, colorList[1], true);
      }    
     
      oss.str("");
      oss << "h1_avg_pT" << "_" << it;
      h1_avg_pT[it] = new TH1D(oss.str().c_str(), "pT average;pT;pT average", nBinsPT, pT_bin);
      h1_avg_pT[it]->Sumw2();
      h1_avg_pT[it]->SetStats(0);
      h1_avg_pT[it]->SetLineColor(2);
      h1_avg_pT[it]->SetLineWidth(2);
      vec_th1_save.push_back(h1_avg_pT[it]);
      oss.str("");
      oss << "h1_raw_pT" << "_" << it;
      h1_raw_pT[it] = new TH1D(oss.str().c_str(), "pT raw;pT;pT raw"    , nBinsPT, pT_bin);
      h1_raw_pT[it]->Sumw2();
      h1_raw_pT[it]->SetStats(0);
      h1_raw_pT[it]->SetLineColor(2);
      h1_raw_pT[it]->SetLineWidth(2);
      /////////////////////// pT ////////////////////////

      oss << "h1_pT_dist" << "_" << it;
      h1_pT_dist[it] = new TH1D(oss.str().c_str(), "pT raw;pT;pT raw", 100, 0, 2.5);
      h1_pT_dist[it]->Sumw2();
      h1_pT_dist[it]->SetStats(0);
      h1_pT_dist[it]->SetLineColor(2);
      h1_pT_dist[it]->SetLineWidth(2);
     
   }
}

void initOutput(){
   saveFile = new TFile("results.root", "recreate");
   if( !saveFile->IsOpen() ){
      cout << "saveFile was not created" << endl;
      exit(0);
   }
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
   oss << "results_nuc_list_v42_" << argv[1];
   oss << "/nbins_" << nbins << "_range_" << fitMax << "_POTMODE_" << potApp << "_cuts_" << argv[5];// << tightApp;
   outDir = oss.str();

   gSystem->mkdir(outDir.c_str(), true);
   gSystem->cd   (outDir.c_str()      );

   if( ! fittingMode ) initOutput();
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

//      if(   dimuon.xF <= 0. ) continue;

      if( !   trackIsValid_2111_v42( trackP, rs[rs_current] ) ) continue;
      if( !   trackIsValid_2111_v42( trackN, rs[rs_current] ) ) continue;
      if( !  dimuonIsValid_2111_v42( dimuon, rs[rs_current], looseCuts ) ) continue;
      if( ! tracksAreValid_2111_v42( trackP, trackN, dimuon ) ) continue;

      if( tightCuts )
         if( ! tightMode_2111_v42(trackP, trackN, dimuon) ) continue;

      if( beamDumpCuts ) 
         if( ! beamDump_2111_v42(trackP, trackN, dimuon, rs[rs_current]) ) continue;

      if( event->targetPos == 3 ) nDimuonsRS++;//counting the # of dimuons
      int ix = getIX2(dimuon.x2);
      if( ix == -1 ) continue;
      if( event->targetPos > 7 ) continue;

      h1_rf_x2      [event->targetPos][ix]->Fill(event->inte_t[16]);
      h1_avg_inte_x2[event->targetPos][ix]->Fill(event->inte_t[16], event->inte_t[16]);
         
//////////////////////////////////////
      int iMass = getIMass(dimuon.mass);
      if( iMass != -1 ){
         h1_rf_mass      [event->targetPos][iMass]->Fill(event->inte_t[16]                   );
         h1_avg_inte_mass[event->targetPos][iMass]->Fill(event->inte_t[16], event->inte_t[16]);
      }
      int iX1 = getIX1(dimuon.x1);
      if( iX1   != -1 ){
         h1_rf_x1      [event->targetPos][iX1]->Fill(event->inte_t[16]                   );
         h1_avg_inte_x1[event->targetPos][iX1]->Fill(event->inte_t[16], event->inte_t[16]);
      }
      int iXF = getIXF(dimuon.xF);
      if( iXF   != -1 ){
         h1_rf_xF      [event->targetPos][iXF]->Fill(event->inte_t[16]                   );
         h1_avg_inte_xF[event->targetPos][iXF]->Fill(event->inte_t[16], event->inte_t[16]);
      }
      double pT = sqrt( pow(dimuon.vtx_mom.X(), 2) + pow(dimuon.vtx_mom.Y(), 2) );
      int iPT = getIPT(pT);
      if( iPT   != -1 ){
         h1_rf_pT      [event->targetPos][iPT]->Fill(event->inte_t[16]                   );
         h1_avg_inte_pT[event->targetPos][iPT]->Fill(event->inte_t[16], event->inte_t[16]);
      }
      h1_pT_dist[event->targetPos]->Fill(pT);
/////////////////////////////////////


      if( fitMax > event->inte_t[16] ){
         if( event->targetPos != 2 ){
            h1_avg_x2[event->targetPos]->Fill(dimuon.x2, dimuon.x2);
            h1_raw_x2[event->targetPos]->Fill(dimuon.x2);

            if( iMass != -1 ){
               h1_avg_mass[event->targetPos]->Fill(dimuon.mass, dimuon.mass);
               h1_raw_mass[event->targetPos]->Fill(dimuon.mass);
            }
            if( iX1 != -1 ){
               h1_avg_x1[event->targetPos]->Fill(dimuon.x1, dimuon.x1);
               h1_raw_x1[event->targetPos]->Fill(dimuon.x1);
            }
            if( iXF != -1 ){
               h1_avg_xF[event->targetPos]->Fill(dimuon.xF, dimuon.xF);
               h1_raw_xF[event->targetPos]->Fill(dimuon.xF);
            }
            if( iPT != -1 ){
               h1_avg_pT[event->targetPos]->Fill(pT, pT);
               h1_raw_pT[event->targetPos]->Fill(pT);
            }
         }
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
//void normHist(TH1D*& h1, int tgt){
//   h1->Scale(1/getPoT(rs, tgt));
//}
void getQuantiles(TH1D* h1, int into, const char* name){
   ofstream ofs;
   ofs.open(name);
   double xq[100];
   for( int iq = 0 ; iq <= into ; iq++ ) xq[iq] = 1./into * iq;
   double yq[100];
   h1->GetQuantiles(into+1, yq, xq);
   for( int iq = 0 ; iq <= into ; iq++ )
      ofs << yq[iq] << "\t";
   ofs << endl;
   ofs.close();
}

void getRatioAndError(double& ratio, double& error, 
                      double  nume, double  deno1, double  deno2, double  subtE, double  subtN, 
                      double eNume, double eDeno1, double eDeno2, double eSubtE, double eSubtN, 
                      int it){

   if( it != 5 && it != 6 && it != 7 ){
      cout << "!!! TARGET TYPE SHOULD BE 5, 6, OR 7 !!!\n" << "YOUR TARGET TYPE IS: " << it << endl;
      exit(0);
   }
   double e_A = 1;
   double e_d = 1;
   double e_h = 1;
   
// double CommonFactor = aMassNum[3] / aMassNum[it] * aMass_a[it] / avogadro / e_A / T[it] / (aMass_a[3] / TDD / avogadro / e_d);
// double CommonFactor = aMassNum[3] / aMassNum[it] * aMass_a[it] /            e_A / T[it] / (aMass_a[3] / TDD /            e_d);
   double CommonFactor =  1          /                                         e_A / T[it] / (    1      / TDD /            e_d);

// double sigma_h = ( deno2 / AH    - subtE ) * aMass_a [ 1] / avogadro / e_h / THH  ;
   double sigma_h = ( deno2 / AH    - subtE ) * aMassNum[ 1] / avogadro / e_h / THH  ;
   double sigma_A = ( nume  / A[it] - subtN );// * aMass_a[it] / avogadro / e_A / T[it];
// double sigma_d = ( deno1 / AD    - subtE - THD * avogadro * sigma_h * e_d / aMass_a [1]);// * aMass_a[3] / TDD / avogadro / e_d;
   double sigma_d = ( deno1 / AD    - subtE - THD * avogadro * sigma_h * e_d / aMassNum[1]);// * aMass_a[3] / TDD / avogadro / e_d;

   if( sigma_A < 0 || sigma_d < 0 || sigma_A / sigma_d * CommonFactor > 2.0 )
      cout << "!!! OUTSIDE OF THE RANGE !!! | " << variableType << " | TARGET TYPE: " << it << " | INTENSITY_BIN: " << iInteBin << " | Val. Bin: " << iValBin << " : sigma_A = " << sigma_A << " / sigma_d = " << sigma_d << " / Ratio = " << sigma_A / sigma_d * CommonFactor << endl; 

   double  eNumeCal = 1 / sigma_d / A[it];
   double eDeno1Cal = sigma_A / pow(sigma_d, 2) / AD;
// double eDeno2Cal = sigma_A / pow(sigma_d, 2) * (THD * avogadro * e_d / aMass_a [1]) * ( 1 / AH * aMass_a [ 1] / avogadro / e_h / THH);
   double eDeno2Cal = sigma_A / pow(sigma_d, 2) * (THD * avogadro * e_d / aMassNum[1]) * ( 1 / AH * aMassNum[ 1] / avogadro / e_h / THH);
   double eSubtECal = sigma_A / pow(sigma_d, 2) * (THD * e_d / e_h / THH - 1 );
   double eSubtNCal = 1 / sigma_d;

   ratio = CommonFactor * sigma_A / sigma_d;
   error = CommonFactor * 
      sqrt(
         pow(  eNume *  eNumeCal, 2) +
         pow( eDeno1 * eDeno1Cal, 2) +
         pow( eDeno2 * eDeno2Cal, 2) +
         pow( eSubtE * eSubtECal, 2) +
         pow( eSubtN * eSubtNCal, 2)
         )
      ;   
}

TGraphErrors* getRatio(TH1D* h1_nume, TH1D* h1_deno1, TH1D* h1_deno2, TH1D* h1_subtE, TH1D* h1_subtN, 
                       TH1D* h1_ave, const char* gn, int it){
   TGraphErrors* tge = new TGraphErrors();
   tge->SetName(gn);
   for( int ib = 1 ; ib <= h1_nume->GetNbinsX() ; ib++ ){
      iInteBin = ib;
      double  nume = h1_nume ->GetBinContent(ib);
      double deno1 = h1_deno1->GetBinContent(ib);
      double deno2 = h1_deno2->GetBinContent(ib);
      double subtE = h1_subtE->GetBinContent(ib);
      double subtN = h1_subtN->GetBinContent(ib);
      
      // if( deno - subt < 0 ) { cout << h1_nume->GetBinCenter(ib) << ": !! deno ( " << deno << " ) - subt ( " << subt << " ) < 0 !! " << endl; isValidNumbers = false;; }
      // if( nume - subt < 0 ) { cout << h1_nume->GetBinCenter(ib) << ": !! nume ( " << nume << " ) - subt ( " << subt << " ) < 0 !! " << endl; isValidNumbers = false;; }

      double  eNume = h1_nume ->GetBinError(ib);
      double eDeno1 = h1_deno1->GetBinError(ib);
      double eDeno2 = h1_deno2->GetBinError(ib);
      double eSubtE = h1_subtE->GetBinError(ib);
      double eSubtN = h1_subtN->GetBinError(ib);

      double ratio, error;

      getRatioAndError(ratio, error, nume, deno1, deno2, subtE, subtN, eNume, eDeno1, eDeno2, eSubtE, eSubtN, it);

      tge->SetPoint     ( tge->GetN()  , h1_ave->GetBinContent(ib), ratio );
      tge->SetPointError( tge->GetN()-1,                         0, error );
      tge->SetMarkerColor(2);
      tge->SetMarkerSize(1.2);
      tge->SetMarkerStyle(20);

   }
   vec_tge_save.push_back(tge);
   return tge;
}

void printYields(TH1D* h1_rf[], TH1D* h1_avg_inte[], TH1D* h1_avg_val, double val[], const char* valName, int nBins, const char* name, bool texMode=false){
   string fileName = (string)("yields_") + name;
   if( texMode ) fileName += ".tex";
   else          fileName += ".txt";
   ofstream ofs;
   ofs.open(fileName.c_str());
   string sep = texMode ? " & " : "\t";
   if( texMode ){
      ofs << "\\begin{landscape}\n\\begin{table}\n\\center\n\\begin{tabular}{c@{-}c";
      for( int iBins = 0 ; iBins < nBins ; iBins++ ) ofs << "||c|r";
      ofs << "}\n";
      ofs << "\\multicolumn{2}{c||}{}";
      for( int iBins = 0 ; iBins < nBins ; iBins++ ){
         ofs << sep << "\\multicolumn{2}{c";
         if( iBins < nBins - 1 ) ofs << "||";
         ofs << "}{$" << val[iBins] << " < " << valName << " < " << val[iBins+1] << "$}";
      }
      ofs << "\\\\\\cline{3-" << 2 * (nBins+1) << "}" << endl;
      ofs << "\\multicolumn{2}{c||}{Int. Range}";
      for( int iBins = 0 ; iBins < nBins ; iBins++ ){
         ofs << sep << "\\multicolumn{2}{c";
         if( iBins < nBins - 1 ) ofs << "||";
         ofs <<"}{$\\langle " << valName << "\\rangle = " << h1_avg_val->GetBinContent(iBins+1) << "$}";
      }
      ofs << "\\\\\\cline{3-" << 2 * (nBins+1) << "}" << endl;
      ofs << "\\multicolumn{2}{c||}{}";
      for( int iBins = 0 ; iBins < nBins ; iBins++ ){
         ofs << sep << "Int. Avg" << sep << "\\multicolumn{1}{c";
         if( iBins < nBins - 1 ) ofs << "||";
         ofs <<"}{Yield}";
      }
      ofs << "\\\\\\hline\\hline" << endl;
   } 
   double step = fitMax / nbins;
   ofs << 0 << sep << fitMax << sep;
   for( int iBins = 0 ; iBins < nBins ; iBins++ ){
      ofs << h1_rf[iBins]->GetMean() << sep << h1_rf[iBins]->Integral(1, h1_rf[iBins]->GetNbinsX());
      if( iBins < nBins - 1 ) ofs <<  sep;
      else 
         if( texMode ) ofs << "\\\\\\hline" << endl;
         else          ofs                  << endl;
   }
   
   for( int iI = 1 ; iI <= h1_rf[0]->GetNbinsX() ; iI++ ){
      ofs << step * (iI-1) << sep << step * iI << sep;
      for( int iBins = 0 ; iBins < nBins ; iBins++ ){
         ofs << h1_avg_inte[iBins]->GetBinContent(iI) << sep << h1_rf[iBins]->GetBinContent(iI);
         if( iBins < nBins - 1 ) ofs <<  sep;
         else 
            if( texMode ) ofs << "\\\\" << endl;
            else          ofs           << endl;
      }
   }
   if( texMode )
      ofs << "\\end{tabular}\n\\end{table}\n\\end{landscape}" << endl;      

   ofs.close();

   if( !texMode ){
      fileName = (string)"average_" + name + ".txt";
      ofs.open( fileName.c_str() );
      for( int iBins = 0 ; iBins < nBins ; iBins++ ){
         ofs <<  h1_avg_val->GetBinContent(iBins+1);
         if( iBins < nBins - 1 ) ofs <<  sep;
         else                    ofs << endl;
      }
      ofs.close();
   }

}

void anaMain(){
   ostringstream oss;

   TDD = targetLength * rho_d / getRDAverage(rs_cont, nDimuons, potMode) * (getFracDAverage(rs_cont, nDimuons, potMode)+getFracHDAverage(rs_cont, nDimuons, potMode)/2.);
   THD = targetLength * rho_d / getRDAverage(rs_cont, nDimuons, potMode) *                                              getFracHDAverage(rs_cont, nDimuons, potMode)/2. ;
   THH = targetLength * rho_h;

   for( int it = 5 ; it < 8 ; it++){
      T[it] = thickness[it]*rho[it];
   }
   //T5 = Fe, T6 = C, T7 = W
   
   lambda_d = getLambdaLD2TAverage(rs_cont, nDimuons, potMode);
   AH = lambda_h * ( 1 - exp( - targetLength / lambda_h ) )/targetLength;
   AD = lambda_d * ( 1 - exp( - targetLength / lambda_d ) )/targetLength;
   
   for( int it = 5 ; it < 8 ; it++){
      A[it] = lambda[it]*( 1 - exp( - thickness[it] / lambda[it] ) )/thickness[it];
   }

   //A5 = Fe, A6 = C, A7 = W 
   TH1D* h1_raw_temp[8];
   TH1D* h1_raw_yield[8][20];
   /////// x2 /////////
   variableType = "x2";
   for( int ix = 0 ; ix < nBinsX2 ; ix++ ){
      for( int it = 1 ; it < 8 ; it++ ){
         if( it != 4 )
            h1_raw_temp [it] = (TH1D*)h1_rf_x2[ 3][ix]->Clone();
         h1_raw_yield[it][ix] = (TH1D*)h1_rf_x2[it][ix]->Clone();
      }
      for( int it = 1 ; it <=7 ; it++ ){
         if( it != 4 )
            h1_raw_temp[it]->Add(h1_rf_x2[it][ix]);
         else
            h1_raw_temp [it] = (TH1D*)h1_rf_x2[it][ix]->Clone();
         h1_rf_x2[it][ix]->Scale(1/getRawPoT(rs, it));

         if( it == 4 ){
            h1_avg_inte_x2[it][ix]->Divide(h1_raw_temp[it]);
         }
         if (it < 5) continue;
         h1_avg_inte_x2[it][ix]->Add(h1_avg_inte_x2[ 3][ix]);
         h1_avg_inte_x2[it][ix]->Divide(h1_raw_temp[it]);
         oss.str("");
         oss << "gr_rf_ratio_x2_" << it << "_" << ix;
         gr_rf_ratio_x2[it][ix] = getRatio(h1_rf_x2[it][ix], h1_rf_x2[3][ix], h1_rf_x2[1][ix], h1_rf_x2[2][ix], h1_rf_x2[4][ix], h1_avg_inte_x2[it][ix], oss.str().c_str(), it);
         
         if( ix != 0  ) continue;
         h1_avg_x2[it]->Add   (h1_avg_x2[ 3]);
         h1_raw_x2[it]->Add   (h1_raw_x2[ 3]);
         h1_avg_x2[it]->Divide(h1_raw_x2[it]);
      }
   }

   string yields_file = "yields_";
   for( int it = 4 ; it <= 7 ; it++ ){
      printYields(h1_raw_yield[it], h1_avg_inte_x2[it], h1_avg_x2[it], x2_bin, "x2", nBinsX2, 
                  ((string)((string)targetType[it]+"_x2")).c_str());
      printYields(h1_raw_yield[it], h1_avg_inte_x2[it], h1_avg_x2[it], x2_bin, "x_2", nBinsX2, 
                  ((string)((string)targetType[it]+"_x2")).c_str(), true);
   }

   /////// x2 /////////

   /////// mass /////////
   variableType = "mass";
   for( int iMass = 0 ; iMass < nBinsMass ; iMass++ ){
      for( int it = 1 ; it < 8 ; it++ ){
         h1_raw_temp[it] = (TH1D*)h1_rf_mass[3][iMass]->Clone();
         h1_raw_yield[it][iMass] = (TH1D*)h1_rf_mass[it][iMass]->Clone();
      }
      for( int it = 1 ; it <=7 ; it++ ){
         h1_raw_temp[it]->Add(h1_rf_mass[it][iMass]);
         h1_rf_mass[it][iMass]->Scale(1/getRawPoT(rs, it));
         if (it < 5) continue;
         h1_avg_inte_mass[it][iMass]->Add(h1_avg_inte_mass[3][iMass]);
         h1_avg_inte_mass[it][iMass]->Divide(h1_raw_temp[it]);
         oss.str("");
         oss << "gr_rf_ratio_mass_" << it << "_" << iMass;
         gr_rf_ratio_mass[it][iMass] = getRatio(h1_rf_mass[it][iMass], h1_rf_mass[3][iMass], h1_rf_mass[1][iMass], h1_rf_mass[2][iMass], h1_rf_mass[4][iMass], h1_avg_inte_mass[it][iMass], oss.str().c_str(), it);
         if( iMass != 0  ) continue;
         h1_avg_mass[it]->Add   (h1_avg_mass[ 3]);
         h1_raw_mass[it]->Add   (h1_raw_mass[ 3]);
         h1_avg_mass[it]->Divide(h1_raw_mass[it]);
      }
   }
   for( int it = 5 ; it <= 7 ; it++ ){
      printYields(h1_raw_yield[it], h1_avg_inte_mass[it], h1_avg_mass[it], mass_bin, "mass", nBinsMass, 
                  ((string)((string)targetType[it]+"_mass")).c_str());
      printYields(h1_raw_yield[it], h1_avg_inte_mass[it], h1_avg_mass[it], mass_bin, "x_2", nBinsMass, 
                  ((string)((string)targetType[it]+"_mass")).c_str(), true);
   }
   /////// mass /////////

   /////// x1 /////////
   variableType = "x1";
   for( int iX1 = 0 ; iX1 < nBinsX1 ; iX1++ ){
      for( int it = 1 ; it < 8 ; it++ ){
         h1_raw_temp[it] = (TH1D*)h1_rf_x1[3][iX1]->Clone();
         h1_raw_yield[it][iX1] = (TH1D*)h1_rf_x1[it][iX1]->Clone();
      }
      for( int it = 1 ; it <=7 ; it++ ){
         h1_raw_temp[it]->Add(h1_rf_x1[it][iX1]);
         h1_rf_x1[it][iX1]->Scale(1/getRawPoT(rs, it));
         if (it < 5) continue;
         h1_avg_inte_x1[it][iX1]->Add(h1_avg_inte_x1[3][iX1]);
         h1_avg_inte_x1[it][iX1]->Divide(h1_raw_temp[it]);
         oss.str("");
         oss << "gr_rf_ratio_x1_" << it << "_" << iX1;
         gr_rf_ratio_x1[it][iX1] = getRatio(h1_rf_x1[it][iX1], h1_rf_x1[3][iX1], h1_rf_x1[1][iX1], h1_rf_x1[2][iX1], h1_rf_x1[4][iX1], h1_avg_inte_x1[it][iX1], oss.str().c_str(), it);
         if( iX1 != 0  ) continue;
         h1_avg_x1[it]->Add   (h1_avg_x1[ 3]);
         h1_raw_x1[it]->Add   (h1_raw_x1[ 3]);
         h1_avg_x1[it]->Divide(h1_raw_x1[it]);
      }
   }
   for( int it = 5 ; it <= 7 ; it++ ){
      printYields(h1_raw_yield[it], h1_avg_inte_x1[it], h1_avg_x1[it], x1_bin, "x1", nBinsX1, 
                  ((string)((string)targetType[it]+"_x1")).c_str());
      printYields(h1_raw_yield[it], h1_avg_inte_x1[it], h1_avg_x1[it], x1_bin, "x_2", nBinsX1, 
                  ((string)((string)targetType[it]+"_x1")).c_str(), true);
   }
   /////// x1 /////////

   /////// xF /////////
   variableType = "xF";
   for( int iXF = 0 ; iXF < nBinsXF ; iXF++ ){
      for( int it = 1 ; it < 8 ; it++ ){
         h1_raw_temp[it] = (TH1D*)h1_rf_xF[3][iXF]->Clone();
         h1_raw_yield[it][iXF] = (TH1D*)h1_rf_xF[it][iXF]->Clone();
      }
      for( int it = 1 ; it <=7 ; it++ ){
         h1_raw_temp[it]->Add(h1_rf_xF[it][iXF]);
         h1_rf_xF[it][iXF]->Scale(1/getRawPoT(rs, it));
         if (it < 5) continue;
         h1_avg_inte_xF[it][iXF]->Add(h1_avg_inte_xF[3][iXF]);
         h1_avg_inte_xF[it][iXF]->Divide(h1_raw_temp[it]);
         oss.str("");
         oss << "gr_rf_ratio_xF_" << it << "_" << iXF;
         gr_rf_ratio_xF[it][iXF] = getRatio(h1_rf_xF[it][iXF], h1_rf_xF[3][iXF], h1_rf_xF[1][iXF], h1_rf_xF[2][iXF], h1_rf_xF[4][iXF], h1_avg_inte_xF[it][iXF], oss.str().c_str(), it);
         if( iXF != 0  ) continue;
         h1_avg_xF[it]->Add   (h1_avg_xF[ 3]);
         h1_raw_xF[it]->Add   (h1_raw_xF[ 3]);
         h1_avg_xF[it]->Divide(h1_raw_xF[it]);
      }
   }
   for( int it = 5 ; it <= 7 ; it++ ){
      printYields(h1_raw_yield[it], h1_avg_inte_xF[it], h1_avg_xF[it], xF_bin, "xF", nBinsXF, 
                  ((string)((string)targetType[it]+"_xF")).c_str());
      printYields(h1_raw_yield[it], h1_avg_inte_xF[it], h1_avg_xF[it], xF_bin, "x_2", nBinsXF, 
                  ((string)((string)targetType[it]+"_xF")).c_str(), true);
   }
   /////// xF /////////

   /////// pT /////////
   variableType = "pT";
   for( int iPT = 0 ; iPT < nBinsPT ; iPT++ ){
      iValBin = iPT;
      for( int it = 1 ; it < 8 ; it++ ){
         h1_raw_temp[it] = (TH1D*)h1_rf_pT[3][iPT]->Clone();
         h1_raw_yield[it][iPT] = (TH1D*)h1_rf_pT[it][iPT]->Clone();
      }
      for( int it = 1 ; it <=7 ; it++ ){
         h1_raw_temp[it]->Add(h1_rf_pT[it][iPT]);
         h1_rf_pT[it][iPT]->Scale(1/getRawPoT(rs, it));
         if (it < 5) continue;
         h1_avg_inte_pT[it][iPT]->Add(h1_avg_inte_pT[3][iPT]);
         h1_avg_inte_pT[it][iPT]->Divide(h1_raw_temp[it]);
         oss.str("");
         oss << "gr_rf_ratio_pT_" << it << "_" << iPT;
         gr_rf_ratio_pT[it][iPT] = getRatio(h1_rf_pT[it][iPT], h1_rf_pT[3][iPT], h1_rf_pT[1][iPT], h1_rf_pT[2][iPT], h1_rf_pT[4][iPT], h1_avg_inte_pT[it][iPT], oss.str().c_str(), it);
         if( iPT != 0  ) continue;
         h1_avg_pT[it]->Add   (h1_avg_pT[ 3]);
         h1_raw_pT[it]->Add   (h1_raw_pT[ 3]);
         h1_avg_pT[it]->Divide(h1_raw_pT[it]);
      }
   }
   for( int it = 5 ; it <= 7 ; it++ ){
      printYields(h1_raw_yield[it], h1_avg_inte_pT[it], h1_avg_pT[it], pT_bin, "pT", nBinsPT, 
                  ((string)((string)targetType[it]+"_pT")).c_str());
      printYields(h1_raw_yield[it], h1_avg_inte_pT[it], h1_avg_pT[it], pT_bin, "x_2", nBinsPT, 
                  ((string)((string)targetType[it]+"_pT")).c_str(), true);
   }
   /////// pT /////////
}

////////////////
//    DRAW    //
////////////////
void writeAverage(){
   // ofstream ofs;
   // ofs.open("x2_average.txt");
   // for( int ib = 1; ib <= h1_avg_x2->GetNbinsX() ; ib++ ){
   //    ofs << h1_avg_x2->GetBinContent(ib) << "\t";
   // }
   // ofs << endl;
   // ofs.close();
}

void drawTH1(TH1D* h1){
   ostringstream oss;
   TCanvas* c1 = new TCanvas("c1", "", 800, 600);
   c1->SetGrid();
   h1->Draw("HIST");
   oss.str("");
   oss << h1->GetName() << "." << saveType;
   c1->SaveAs(oss.str().c_str());
   delete c1;
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
   TFile* saveCSRFile = new TFile("results_csr.root", "recreate");
   saveCSRFile->cd();
   for( int ic = 0 ; ic < (int)vec_csr_save.size() ; ic++ )
      vec_csr_save[ic]->Write();
   saveCSRFile->Close();
}

void drawComp(){
   ostringstream oss;
   TCanvas* c1 = new TCanvas("c1", "", 800, 600);
   c1->SetGrid();
   TLegend* leg;
   for( int ix = 0 ; ix < 8 ; ix++ ){
      leg = new TLegend(0.1, 0.7, 0.3, 0.9);
      double max = 0 ;
      for( int it = 1 ; it <= 8 ; it++ ) 
         max = max > h1_rf_x2[it][ix]->GetMaximum() ? max : h1_rf_x2[it][ix]->GetMaximum();
      for( int it = 1 ; it <= 8 ; it++ ) {
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

void drawMain(){
   ostringstream oss;
   ostringstream title;

   if( ! fittingMode ){
      writeAverage();
//      drawComp();
      for( int it = 1 ; it < 8 ; it++){
         oss.str("");
         oss << "quantiles_pt_" << it << ".txt";
         getQuantiles(h1_pT_dist[it], nBinsPT, oss.str().c_str());
      }
      
      
      for( int ix = 0 ; ix < 8 ; ix++ ){
         oss.str("");
         oss << "ratio_intensity_sum2_" << ix << ".txt";
//         printRatio(gr_rf_ratio_x2[ix], oss.str().c_str());
      }
            
      for( int ig = 0 ; ig < (int)vec_tge_save.size() ; ig++ ){
         drawTGraphErrors(vec_tge_save[ig], ";Intensity;CSR");
      }



      // for( int ih = 0 ; ih < (int)vec_th1_save.size() ; ih++ ){
      //    drawTH1(vec_th1_save[ih]);
      // }
      
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

   for( int it = 5 ; it <= 7 ; it++ ){
      for( int iv = 0 ; iv < (int)valSet.size() ; iv++ ){
         ////////////// POL2 //////////////
         fit = new Fitter2DPol2(); 
         oss.str(""); oss << "gr_rf_ratio_" << valSet[iv] << "_" << it;
         oss_avg.str(""); oss_avg << "h1_avg_" << valSet[iv] << "_" << it;
         fit->Init("results.root", oss.str().c_str(), oss_avg.str().c_str());
         fit->DoFit();
         oss.str(""); oss << "result_pol2_fit_" << it << "_" << valSet[iv] << ".txt";
         fit->PrintFit(oss.str().c_str());
         fit->PrintFit(cout);
         fit->DrawFit();
         oss.str(""); oss << "csr_pol2_" << it << "_as_func_" << valSet[iv];
         vec_csr_save.push_back(fit->getCSRPlot(oss.str().c_str(), valSet[iv], (string)(targetType[it])+"/LD_{2}"));
         ////////////// POL2 //////////////
      
         ////////////// POL1 //////////////
         fit = new Fitter2DPol1(); 
         oss.str(""); oss << "gr_rf_ratio_" << valSet[iv] << "_" << it;
         oss_avg.str(""); oss_avg << "h1_avg_" << valSet[iv] << "_" << it;
         fit->Init("results.root", oss.str().c_str(), oss_avg.str().c_str());
         fit->DoFit();
         oss.str(""); oss << "result_pol1_fit_" << it << "_" << valSet[iv] << ".txt";
         fit->PrintFit(oss.str().c_str());
         fit->PrintFit(cout);
         fit->DrawFit();
         oss.str(""); oss << "csr_pol1_" << it << "_as_func_" << valSet[iv];
         vec_csr_save.push_back(fit->getCSRPlot(oss.str().c_str(), valSet[iv], (string)(targetType[it])+"/LD_{2}"));
         ////////////// POL1 //////////////

         ////////////// COS //////////////
         fit = new Fitter2DCos(); 
         oss.str(""); oss << "gr_rf_ratio_" << valSet[iv] << "_" << it;
         oss_avg.str(""); oss_avg << "h1_avg_" << valSet[iv] << "_" << it;
         fit->Init("results.root", oss.str().c_str(), oss_avg.str().c_str());
         fit->DoFit();
         oss.str(""); oss << "result_cos_fit_" << it << "_" << valSet[iv] << ".txt";
         fit->PrintFit(oss.str().c_str());
         fit->PrintFit(cout);
         fit->DrawFit();
         oss.str(""); oss << "csr_cos_" << it << "_as_func_" << valSet[iv];
         vec_csr_save.push_back(fit->getCSRPlot(oss.str().c_str(), valSet[iv], (string)(targetType[it])+"/LD_{2}"));
         ////////////// COS //////////////

         ////////////// COMMON POL2 //////////////
         fit = new FitterCommonPol2(); 
         oss.str(""); oss << "gr_rf_ratio_" << valSet[iv] << "_" << it;
         oss_avg.str(""); oss_avg << "h1_avg_" << valSet[iv] << "_" << it;
         fit->Init("results.root", oss.str().c_str(), oss_avg.str().c_str());
         fit->DoFit();
         oss.str(""); oss << "result_common_pol2_fit_" << it << "_" << valSet[iv] << ".txt";
         fit->PrintFit(oss.str().c_str());
         fit->PrintFit(cout);
         fit->DrawFit();
         oss.str(""); oss << "csr_common_pol2_" << it << "_as_func_" << valSet[iv];
         vec_csr_save.push_back(fit->getCSRPlot(oss.str().c_str(), valSet[iv], (string)(targetType[it])+"/LD_{2}"));
        ////////////// COMMON POL2 //////////////
      }
   }

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
