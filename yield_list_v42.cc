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
#include <TH2.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include <TPaveStats.h>
#include "Event.h"
#include "Variables.h"
#include "Selection.h"

#define FILE_EXIST 1
#define DIR_EXIST 2
#define NO_FILE_DIR_EXIST 0

vector<string> fileList;
vector<string> fileListNim3;

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

double     runID;
double   spillID;
double targetPos;
double    RF[17];
double     G2SEM;
double    QIEsum;
double inh_thres;

vector<int> rs;
int         rs_current;
string outDir;
string list;
///// DATA /////
// [target type (dummy, LH2, Empty, LD2)][x2 bins][rs]
TH1D* h1_rf_x2     [4][8][5];
TH1D* h1_rf_x2_sum1[4][8]; 
TH1D* h1_rf_x2_sum2[4][8]; 
TGraphErrors*tge_rf[4][8];
TH1D* h1_rf_average[4][8];

TH1D* h1_rf_x2_kEff     [4][8][5];
TH1D* h1_rf_x2_kEff_sum1[4][8];
TH1D* h1_rf_x2_kEff_sum2[4][8];
TGraphErrors*tge_rf_kEff[4][8];
TH1D* h1_rf_average_kEff[4][8];

TH1D* h1_rf_x2_fine   [4][8]; // [target type (dummy, LH2, Empty, LD2)][x2 bins]
TGraphErrors* tge_fine   [8]; // [x2 bins]

TH1D* h1_raw_tgt[4][8]; // [x2]

TH1D* h1_average[8]; // [x2]
TH1D* h1_raw    [8]; // [x2]

TH1D* h1_average_fine[8]; // [x2]
TH1D* h1_raw_fine    [8]; // [x2]

TH1D* h1_average_x2;
TH1D* h1_raw_x2    ;
TH1D* h1_average_x2_tgt[8]; // targetType
TH1D* h1_raw_x2_tgt    [8]; // targetType
vector<double> nDimuons;
double nDimuonsRS;

TH1D* h1_nim3    [4][5];
TH1D* h1_nim3_sum[4];

TH2D* h2_qiesum_g2sem;

TH1D* h1_rf_all_kEff        [4][5];
TH1D* h1_rf_all_kEff_sum1   [4]   ;
TGraphErrors*tge_rf_all_kEff[4]   ;
TH1D* h1_rf_all_average_kEff[4]   ;
///// FIXED VARIABLES /////
int color[4] = {0, kBlack, kRed, kBlue};
double factor = 2.3105;
double     purity;
double effAtoms_h;
double effAtoms_d;

///// ADJESTABLE /////
const char* saveType = "png";

int binRange = 60000;
int   nbins =    15;

int    nBinsrf =     300;
double   maxrf =   90000;
double  fitMax = binRange;

double cutMass;
bool   tightCuts;

int bugFinder = 0;

double fineBin[] = {0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 16000, 20000, 24000, 28000, 32000, 36000, 40000, 44000, 48000, 52000, 56000, 60000};
int fineBinNum = sizeof fineBin / sizeof fineBin[0] - 1;

double adjustment;
int prevSpill = 0;
// int prevRun   = 0;
// int runCount  = 0;
// int histID    = 0;
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

      // oss.str("");
      // oss << fileDirNim3 << rs_tmp << "/NIM3_profiles.root";
      // fileListNim3.push_back(oss.str().c_str());

      rs.push_back(rs_tmp);
   }

//   fileListNim3.push_back("/seaquest/users/arunts/NIM3/NIM3_profiles/special_runs/NIM3_profiles.root");
   fileListNim3.push_back("/seaquest/users/knagai/test_arun/extract/beam/special/NIM3_profiles.root");
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

bool initInputFileNim3(int index){
   ostringstream oss;
   
   cout << " preparing " << fileListNim3[index] << "......" << endl;
   if( fileE(fileListNim3[index].c_str()) == NO_FILE_DIR_EXIST ){
      cout << fileListNim3[index] << " DOES NOT EXIST. GO TO NEXT RUN..." << endl;
      return false;
   }
   file = new TFile(fileListNim3[index].c_str(), "READ");
   tree = (TTree*)file->Get("NIM3");

   tree->SetBranchAddress(    "runID", &    runID);
   tree->SetBranchAddress(  "spillID", &  spillID);
   tree->SetBranchAddress("targetPos", &targetPos);
   tree->SetBranchAddress(    "G2SEM", &    G2SEM);
   tree->SetBranchAddress(   "QIEsum", &   QIEsum);
   tree->SetBranchAddress("inh_thres", &inh_thres);

   for( int irf = 0 ; irf <= 8 ; irf++ ){
      oss.str("");
      oss << "RF0" << irf;
      tree->SetBranchAddress(oss.str().c_str(), & RF[irf+8]);
      
      if( irf == 0 )continue;
      oss.str("");
      oss << "RFm0" << irf;
      tree->SetBranchAddress(oss.str().c_str(), & RF[8-irf]);
   }

   cout << "  Analyze " << fileListNim3[index].c_str() << "..." << endl;
   return true;
}

// void initFiles(){
//    ostringstream oss;
//    oss.str("");
//    oss << fileDir << rs << "/" << rs << "_with_cuts.root";
//    file  = new TFile(oss.str().c_str(), "READ");
//    tree = (TTree*)file->Get("save");
//    tree->SetBranchAddress("event", &event);
// }

void initTH1(TH1D*& h1, const char* name, const char* title, int nBin, double min, double max, int color){
   h1 = new TH1D(name, title, nBin, min, max);
   h1->Sumw2();
   h1->SetStats(0);
   h1->SetLineColor(color);
   h1->SetLineWidth(2);
}

void initTH1(){
   ostringstream oss;
   ostringstream title;

   for( int ix = 0 ; ix < 8 ; ix++ ){
      for( int it = 1 ; it < 4 ; it++ ){
         for( int irs = 0 ; irs < (int)fileList.size() ; irs++ ){
            oss.str("");
            oss << "h1_rf_x2_" << rs[irs] << "_" << it << "_" << ix;
            title.str("");
            title << "FPGA1 profile (" << x2_bin[ix] << " < x_{T} < " << x2_bin[ix+1] 
                  << ";Intensity;Counts";
            initTH1(h1_rf_x2[it][ix][irs], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, color[it]);

            oss.str("");
            oss << "h1_rf_x2_kEff" << rs[irs] << "_" << it << "_" << ix;
            title.str("");
            title << "FPGA1 Profile / kEff (" << x2_bin[ix] << " < x_{T} < " << x2_bin[ix+1] 
                  << ";Intensity;Counts";
            initTH1(h1_rf_x2_kEff[it][ix][irs], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, color[it]);

            if( ix != 0 ) continue;

            oss.str("");
            oss << "h1_rf_all_kEff" << rs[irs] << "_" << it;
            title.str("");
            title << "FPGA1 Profile / kEff (ALL x_{2});Intensity;Counts";
            initTH1(h1_rf_all_kEff[it][irs], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, color[it]);            
         }

         oss.str("");
         oss << "h1_rf_average_"<< it << "_" << ix;
         title.str("");
         title << "Intensity Average (" << x2_bin[ix] << " < x_{T} < " << x2_bin[ix+1] 
               << ";Intensity;Counts";
         initTH1(h1_rf_average[it][ix], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, color[it]);

         oss.str("");
         oss << "h1_rf_average_kEff_" << it << "_" << ix;
         title.str("");
         title << "Intensity Average (" << x2_bin[ix] << " < x_{T} < " << x2_bin[ix+1] 
               << ";Intensity;Counts";
         initTH1(h1_rf_average_kEff[it][ix], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, color[it]);

         if( ix != 0 ) continue;
         oss.str("");
         oss << "h1_rf_all_average_kEff_" << it;
         title.str("");
         title << "Intensity Average (ALL x_{2});Intensity;Counts";
         initTH1(h1_rf_all_average_kEff[it], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, color[it]);

      }

      oss.str("");
      oss << "h1_average_" << ix;
      title.str("");
      title << "Average (" << x2_bin[ix] << " < x_{2} < " << x2_bin[ix+1] 
            << ");Intensity;Counts";
      initTH1(h1_average[ix], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, 2);

      oss.str("");
      oss << "h1_raw_" << ix;
      title.str("");
      title << "raw (" << x2_bin[ix] << " < x_{2} < " << x2_bin[ix+1] 
            << ");Intensity;Counts";
      initTH1(h1_raw[ix], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, 2);


      oss.str("");
      oss << "h1_average_fine_" << ix;
      title.str("");
      title << "Average (" << x2_bin[ix] << " < x_{2} < " << x2_bin[ix+1] 
            << ");Intensity;Counts";
      h1_average_fine[ix] = new TH1D(oss.str().c_str(), title.str().c_str(), fineBinNum, fineBin);
      h1_average_fine[ix]->Sumw2();
      h1_average_fine[ix]->SetStats(0);
      h1_average_fine[ix]->SetLineColor(2);
      h1_average_fine[ix]->SetLineWidth(2);

      oss.str("");
      oss << "h1_raw_fine_" << ix;
      title.str("");
      title << "Raw (" << x2_bin[ix] << " < x_{2} < " << x2_bin[ix+1] 
            << ");Intensity;Counts";
      h1_raw_fine[ix] = new TH1D(oss.str().c_str(), title.str().c_str(), fineBinNum, fineBin);
      h1_raw_fine[ix]->Sumw2();
      h1_raw_fine[ix]->SetStats(0);
      h1_raw_fine[ix]->SetLineColor(2);
      h1_raw_fine[ix]->SetLineWidth(2);

   }

   for( int it = 1 ; it < 4 ; it++ ){
      oss.str("");
      oss << "h1_average_x2_tgt_" << it;
      h1_average_x2_tgt[it] = new TH1D(oss.str().c_str(), "", nBinsX2, x2_bin);
      h1_average_x2_tgt[it]->Sumw2();
      h1_average_x2_tgt[it]->SetStats(0);
      h1_average_x2_tgt[it]->SetMarkerColor( colorList[it]);
      h1_average_x2_tgt[it]->SetMarkerStyle(markerList[it]);
      h1_average_x2_tgt[it]->SetMarkerSize(1.2);

      oss.str("");
      oss << "h1_raw_x2_tgt_" << it;
      h1_raw_x2_tgt[it] = new TH1D(oss.str().c_str(), "", nBinsX2, x2_bin);
      h1_raw_x2_tgt[it]->Sumw2();
      h1_raw_x2_tgt[it]->SetStats(0);
      h1_raw_x2_tgt[it]->SetMarkerColor( colorList[it]);
      h1_raw_x2_tgt[it]->SetMarkerStyle(markerList[it]);
      h1_raw_x2_tgt[it]->SetMarkerSize(1.2);

      for( int irs = 0 ; irs < (int)fileList.size() ; irs++ ){
         oss.str("");
         oss << "h1_nim3_" << it << "_" << irs;
         title.str("");
         title << "NIM3;Intensity;Yield";
         initTH1(h1_nim3[it][irs], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, color[it]);
      }
   }


   h1_average_x2 = new TH1D("h1_average_x2", "x_{2} average;x_{2};x_{2} average", nBinsX2, x2_bin);
   h1_average_x2->Sumw2();
   h1_average_x2->SetStats(0);
   h1_average_x2->SetLineColor(2);
   h1_average_x2->SetLineWidth(2);
   h1_raw_x2     = new TH1D("h1_raw_x2"    , "x_{2} raw;x_{2};x_{2} raw"        , nBinsX2, x2_bin);
   h1_raw_x2->Sumw2();
   h1_raw_x2->SetStats(0);
   h1_raw_x2->SetLineColor(2);
   h1_raw_x2->SetLineWidth(2);

   h2_qiesum_g2sem = new TH2D("h2_qiesum_g2sem", "", 100, 2e12, 8e12, 100, 0, 3.5e11);
   h2_qiesum_g2sem->SetStats(0);
   h2_qiesum_g2sem->Sumw2();

}

void initMain(int argc, char* argv[]){
   ostringstream oss;

   if( argc < 5 ){
      cout << "!! NEED 5 ARGUMENTS !!: ROADSET NBINS FITRANGE CUTOFF-MASS [CA cut apply (CA or others), optional]" << endl;
      exit(0);
   }

   initInputRS   (argv[1]);
   nbins   = atoi(argv[2]);
   fitMax  = atof(argv[3]);
   cutMass = atof(argv[4]);

   if( argc >= 5 ) tightCuts  = strcmp(argv[5], "TC") == 0 ? true : false;

   event = new Event();

   const char* tightApp = tightCuts ? "app" : "not";

   oss.str("");
   oss << "results_yield_list_v42_" << argv[1];
   oss << "/nbins_" << nbins << "_range_" << fitMax << "_cutoffMass_" << cutMass << "_tight_" << tightApp;
   outDir = oss.str();

//   initFiles();
   initTH1();

//   shiftX2(0.02);
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

   if( event->targetPos > 3 ) return;

   if( event->runID >= 12598 && event->runID <= 12734 )
      return;



   for( int id = 0 ; id < (int)event->dimuons.size() ; id++ ){
      Dimuon dimuon = event->dimuons[id];
      Track  trackP = event->tracks[findTrackIndex(dimuon.trackID_pos)];
      Track  trackN = event->tracks[findTrackIndex(dimuon.trackID_neg)];
      if( dimuon.mass <= cutMass ) continue;
      if( tightCuts ){
         if( !   trackIsValid_2111_v42( trackP ) ) continue;
         if( !   trackIsValid_2111_v42( trackN ) ) continue;
         if( !  dimuonIsValid_2111_v42( dimuon ) ) continue;
         if( ! tracksAreValid_2111_v42( trackP, trackN, dimuon ) ) continue;
      }

      if( event->targetPos == 3 ) nDimuonsRS++;

      double pedestal = getPedestal(spillID) - adjustment;
      double intensity = ( event->RF[16] - pedestal ) * event->G2SEM / ( event->QIEsum - turns * buckets * pedestal );

      double kEff = getKEff(dimuon.x2, event->occChams[0], event->targetPos);
      h1_rf_all_kEff        [event->targetPos][rs_current]->Fill(intensity,        1./kEff);
      h1_rf_all_average_kEff[event->targetPos]            ->Fill(intensity, intensity/kEff);

      int ix = getIX2(dimuon.x2);
      if( ix == -1 ) continue;

      h1_rf_x2     [event->targetPos][ix][rs_current]->Fill(intensity);
      h1_rf_x2_kEff[event->targetPos][ix][rs_current]->Fill(intensity, 1./kEff);

      h1_rf_average     [event->targetPos][ix]->Fill(intensity, intensity);
      h1_rf_average_kEff[event->targetPos][ix]->Fill(intensity, intensity/kEff);

      // if( event->targetPos == 3 ){
      //    cout << event->inte_t[16] << endl;
      //    exit(0);
      // }

//      h1_rf_x2_fine[event->targetPos][ix]->Fill(event->inte_t[16]);

      // h1_average[ix]->Fill(event->inte_t[16], event->inte_t[16]);
      // h1_raw    [ix]->Fill(event->inte_t[16]);

      // h1_average_fine[ix]->Fill(event->inte_t[16], event->inte_t[16]);
      // h1_raw_fine    [ix]->Fill(event->inte_t[16]);

      // if( fitMax > event->inte_t[16] ){
      //    h1_average_x2_tgt[event->targetPos]->Fill(dimuon.x2, dimuon.x2);
      //    h1_raw_x2_tgt    [event->targetPos]->Fill(dimuon.x2);
      //    if( event->targetPos != 2 ){
      //       h1_average_x2->Fill(dimuon.x2, dimuon.x2);
      //       h1_raw_x2    ->Fill(dimuon.x2);
      //    }
      // }
   }
}

void extNIM3(){

   if( (int)targetPos > 3 ) return;

   if( (int)runID >= 12598 && (int)runID <= 12734 )
      return;

//////// ONLY SPECIAL RUNS ////////
   bool execute = false;
   for( int ir = 0 ; ir < 4 ; ir++){
      if( (int)runID == NIM3SpecialRuns[ir] ){
         execute = true;
         break;
      }
   }
   if( !execute ) return;
//////// ONLY SPECIAL RUNS ////////

   for( int irf = -8 ; irf <= 8 ; irf++ )
      if( RF[irf+8] >= inh_thres ) return;

   if( ! (  G2SEM > 2e12 &&  G2SEM < 1e13 ) ) return;
   if( ! ( QIEsum > 4e10 && QIEsum < 1e12 ) ) return;
      
   double  pedestal = getPedestal(spillID) - adjustment;
   double intensity = ( RF[8] - pedestal ) * G2SEM / ( QIEsum - turns * buckets * pedestal );

   h1_nim3[(int)targetPos][rs_current]->Fill(intensity);

   if( prevSpill == (int)spillID ) return;
   h2_qiesum_g2sem->Fill(G2SEM, QIEsum);
   prevSpill = (int)spillID;
}

void extMain(bool NIM3mode = false){
   int prev = -1;
   int curr =  0;

   for( int ie = 0 ; ie < tree->GetEntries() ; ++ie){    
      tree->GetEntry(ie);
      curr = (ie+1) * 100 / tree->GetEntries();
      if( curr != prev ){
         cout << "\r" << (ie+1) << " / " 
              << tree->GetEntries() << " = " 
              << curr << " %" << flush;
         prev = curr;
      }


      if( ! NIM3mode )
         extData();
      else
         extNIM3();
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

TGraphErrors* getRatio(TH1D* h1_nume, TH1D* h1_deno, TH1D* h1_subt, TH1D* h1_ave){
   TGraphErrors* tge = new TGraphErrors();
   for( int ib = 1 ; ib <= h1_nume->GetNbinsX() ; ib++ ){
      double nume = h1_nume->GetBinContent(ib);
      double deno = h1_deno->GetBinContent(ib);
      double subt = h1_subt->GetBinContent(ib);
      
      if( deno - subt == 0 ) continue;
      double eNume = h1_nume->GetBinError(ib);
      double eDeno = h1_deno->GetBinError(ib);
      double eSubt = h1_subt->GetBinError(ib);
//      double ePuri = getPurityError(rs);

      double ratio = ( ( nume - subt ) * effAtoms_h / ( deno - subt ) / effAtoms_d + purity - 1 ) / purity / 2.;
      double eNumeCal = effAtoms_h / effAtoms_d / purity / ( deno - subt );
      double eDenoCal = effAtoms_h / effAtoms_d / purity * ( nume - subt ) / pow( deno - subt, 2 );
      double eSubtCal = effAtoms_h / effAtoms_d / purity * ( deno - nume ) / pow( deno - subt, 2 );
//      double ePuriCal = ( 1 - ( ( nume - subt ) * effAtoms_h / ( deno - subt ) / effAtoms_d + purity - 1 ) / purity ) / purity;
//      double ePuriCal = ( ( nume - subt ) * effAtoms_h / ( deno - subt ) / effAtoms_d + 1 ) / pow(purity, 2);
      double error = 
         sqrt(
            pow( eNume * eNumeCal, 2) +
            pow( eDeno * eDenoCal, 2) +
            pow( eSubt * eSubtCal, 2)
//            pow( ePuri * ePuriCal, 2)
            ) / 2.;

//      cout << h1_ave->GetBinContent(ib) << endl;

      tge->SetPoint     ( tge->GetN()  , h1_ave->GetBinContent(ib), ratio );
      tge->SetPointError( tge->GetN()-1,                         0, error );
      tge->SetMarkerColor(2);
      tge->SetMarkerSize(1.2);
      tge->SetMarkerStyle(20);
      // double x, y;
      // tge->GetPoint(tge->GetN()-1, x, y);
      // cout << x << "\t" << y << endl;
   }
   return tge;
}

TGraphErrors* convTH12TgeAve(TH1D* h1, TH1D* h1_ave, const char* name){
   TGraphErrors* tge = new TGraphErrors();
   tge->SetName(name);
   for( int ibin = 1 ; ibin <= h1->GetNbinsX() ; ibin++ ){
      tge->SetPoint     (ibin-1, h1_ave->GetBinContent(ibin), h1->GetBinContent(ibin));
      tge->SetPointError(ibin-1,                           0, h1->GetBinError  (ibin));
   }
   return tge;
}

void anaMain(){
   ostringstream oss;

   purity = getPurityAverage(rs, nDimuons);
   effAtoms_h = effAtoms(aMass_h, rho_h, lambda_h);
   effAtoms_d = effAtoms(effLinear(aMass_d, aMass_h, purity), effLinear(rho_d, rho_h, purity), effInverse(lambda_d, lambda_h, purity));
   cout << purity << " | " << effAtoms_h << " | " << effAtoms_d << endl;

////////////////// NIM3 /////////////////
   for( int it = 1 ; it <= 3 ; it++ ){
      for( int irs = 0 ; irs < (int)fileList.size() ; irs++ ){
         if( irs == 0 ){
            h1_nim3_sum[it] = (TH1D*) h1_nim3[it][irs]->Clone();
         }
         else{
            h1_nim3_sum[it]->Add(h1_nim3[it][irs]);
         }

         h1_nim3[it][irs]->Scale(1/h1_nim3[it][irs]->Integral());
      }
   }
////////////////// NIM3 /////////////////

//    for( int ix = 0 ; ix < 8 ; ix++ ){
//       for( int irs = 0 ; irs < (int)fileList.size() ; irs++ ){
//          TH1D* h1_lh2 = (TH1D*) h1_rf_x2[1][ix][irs]->Clone();
//          TH1D* h1_ld2 = (TH1D*) h1_rf_x2[3][ix][irs]->Clone();
//          h1_rf_x2[1][ix][irs]->Divide(h1_nim3[1][irs]);
//          h1_rf_x2[3][ix][irs]->Divide(h1_nim3[3][irs]);

//          if( irs == 0 ){
// //-------------------------------------------------------------------------------------//
//             h1_rf_x2_sum1[1][ix] = (TH1D*)h1_rf_x2[1][ix][0]->Clone();
//             h1_rf_x2_sum1[3][ix] = (TH1D*)h1_rf_x2[3][ix][0]->Clone();
//             oss.str("");
//             oss << "h1_rf_x2_sum1_1_" << ix;
//             h1_rf_x2_sum1[1][ix]->SetName(oss.str().c_str());
//             oss.str("");
//             oss << "h1_rf_x2_sum1_3_" << ix;
//             h1_rf_x2_sum1[3][ix]->SetName(oss.str().c_str());
// //-------------------------------------------------------------------------------------//
//             h1_rf_x2_sum2[1][ix] = (TH1D*)h1_lh2->Clone();
//             h1_rf_x2_sum2[3][ix] = (TH1D*)h1_ld2->Clone();
//             oss.str("");
//             oss << "h1_rf_x2_sum2_1_" << ix;
//             h1_rf_x2_sum2[1][ix]->SetName(oss.str().c_str());
//             oss.str("");
//             oss << "h1_rf_x2_sum2_3_" << ix;
//             h1_rf_x2_sum2[3][ix]->SetName(oss.str().c_str());
// //-------------------------------------------------------------------------------------//
//          }
//          else{
//             h1_rf_x2_sum1[1][ix]->Add(h1_rf_x2[1][ix][irs]);
//             h1_rf_x2_sum1[3][ix]->Add(h1_rf_x2[3][ix][irs]);
//             h1_rf_x2_sum2[1][ix]->Add(h1_lh2              );
//             h1_rf_x2_sum2[3][ix]->Add(h1_ld2              );
//          }
//          delete h1_lh2;
//          delete h1_ld2;
//       }


   TH1D* h1_tmp;
   TH1D* h1_tmp_ave[4][8];
/////////// WITHOUT CORRECTION //////////
   for( int iit = 1 ; iit <= 3 ; iit++ ){
      int it = iit % 3 + 1;
      cout << it << endl;
      for( int ix = 0 ; ix < 8 ; ix++ ){
         for( int irs = 0 ; irs < (int)fileList.size() ; irs++ ){
            h1_tmp = (TH1D*) h1_rf_x2[it][ix][irs]->Clone();
            h1_rf_x2[it][ix][irs]->Divide(h1_nim3[it][irs]);

            if( irs == 0 ){
//-------------------------------------------------------------------------------------//
               h1_rf_x2_sum1[it][ix] = (TH1D*)h1_rf_x2[it][ix][0]->Clone();
               oss.str("");
               oss << "h1_rf_x2_sum1_" << it << "_" << ix;
               h1_rf_x2_sum1[it][ix]->SetName(oss.str().c_str());
//-------------------------------------------------------------------------------------//
               h1_rf_x2_sum2[it][ix] = (TH1D*)h1_tmp->Clone();
               oss.str("");
               oss << "h1_rf_x2_sum2_" << it << "_" << ix;
               h1_rf_x2_sum2[it][ix]->SetName(oss.str().c_str());
//-------------------------------------------------------------------------------------//
            }
            else{
               h1_rf_x2_sum1[it][ix]->Add(h1_rf_x2[it][ix][irs]);
               h1_rf_x2_sum2[it][ix]->Add(h1_tmp               );
            }
            delete h1_tmp;
         }

         if( it != 2 ){
            h1_tmp_ave[it][ix] = (TH1D*)h1_rf_x2_sum2[it][ix]->Clone();
            h1_tmp_ave[it][ix]->Add(h1_rf_x2_sum2[2][ix]);
         }

         h1_rf_x2_sum1[it][ix]->Scale(1/getPoT(rs, it));
         h1_rf_x2_sum2[it][ix]->Divide(h1_nim3_sum[it]);
         h1_rf_x2_sum2[it][ix]->Scale(h1_rf_x2_sum1[it][ix]->Integral() / h1_rf_x2_sum2[it][ix]->Integral());
      }
   }
/////////// WITHOUT CORRECTION //////////
   TH1D* h1_tmp_ave_kEff[4][8];
/////////// WITH    CORRECTION //////////
   for( int ix = 0 ; ix < 8 ; ix++ ){
      for( int iit = 1 ; iit <= 3 ; iit++ ){
         int it = iit % 3 + 1;
         for( int irs = 0 ; irs < (int)fileList.size() ; irs++ ){
            h1_tmp = (TH1D*) h1_rf_x2_kEff[it][ix][irs]->Clone();
            h1_rf_x2_kEff[it][ix][irs]->Divide(h1_nim3[it][irs]);

            if( irs == 0 ){
//-------------------------------------------------------------------------------------//
               h1_rf_x2_kEff_sum1[it][ix] = (TH1D*)h1_rf_x2_kEff[it][ix][0]->Clone();
               oss.str("");
               oss << "h1_rf_x2_kEff_sum1_" << it << "_" << ix;
               h1_rf_x2_kEff_sum1[it][ix]->SetName(oss.str().c_str());
//-------------------------------------------------------------------------------------//
               h1_rf_x2_kEff_sum2[it][ix] = (TH1D*)h1_tmp->Clone();
               oss.str("");
               oss << "h1_rf_x2_kEff_sum2_" << it << "_" << ix;
               h1_rf_x2_kEff_sum2[it][ix]->SetName(oss.str().c_str());
//-------------------------------------------------------------------------------------//
            }
            else{
               h1_rf_x2_kEff_sum1[it][ix]->Add(h1_rf_x2_kEff[it][ix][irs]);
               h1_rf_x2_kEff_sum2[it][ix]->Add(h1_tmp                    );
            }
            delete h1_tmp;
         }

         if( it != 2 ){
            h1_tmp_ave_kEff[it][ix] = (TH1D*)h1_rf_x2_kEff_sum2[it][ix]->Clone();
            h1_tmp_ave_kEff[it][ix]->Add(h1_rf_x2_kEff_sum2[2][ix]);
         }

         h1_rf_x2_kEff_sum1[it][ix]->Scale(1/getPoT(rs, it));
         h1_rf_x2_kEff_sum2[it][ix]->Divide(h1_nim3_sum[it]);
//         h1_rf_x2_kEff_sum2[it][ix]->Scale(h1_rf_x2_kEff_sum1[it][ix]->Integral() / h1_rf_x2_kEff_sum2[it][ix]->Integral());
      }
      // h1_rf_x2_kEff_sum1[1][ix]->Add(h1_rf_x2_kEff_sum1[2][ix], -1);
      // h1_rf_x2_kEff_sum1[3][ix]->Add(h1_rf_x2_kEff_sum1[2][ix], -1);
      // h1_rf_x2_kEff_sum2[1][ix]->Add(h1_rf_x2_kEff_sum2[2][ix], -1);
      // h1_rf_x2_kEff_sum2[3][ix]->Add(h1_rf_x2_kEff_sum2[2][ix], -1);

      h1_rf_x2_kEff_sum2[1][ix]->Scale(h1_rf_x2_kEff_sum1[1][ix]->Integral() / h1_rf_x2_kEff_sum2[1][ix]->Integral());
      h1_rf_x2_kEff_sum2[3][ix]->Scale(h1_rf_x2_kEff_sum1[3][ix]->Integral() / h1_rf_x2_kEff_sum2[3][ix]->Integral());
   }
/////////// WITH    CORRECTION //////////
   TH1D* h1_tmp_ave_all_kEff[4];
/////////// ALL X2 //////////
   for( int iit = 1 ; iit <= 3 ; iit++ ){
      int it = iit % 3 + 1;
      for( int irs = 0 ; irs < (int)fileList.size() ; irs++ ){
         h1_tmp = (TH1D*) h1_rf_all_kEff[it][irs]->Clone();
         h1_rf_all_kEff[it][irs]->Divide(h1_nim3[it][irs]);

         if( irs == 0 ){
//-------------------------------------------------------------------------------------//
            h1_rf_all_kEff_sum1[it] = (TH1D*)h1_rf_all_kEff[it][0]->Clone();
            oss.str("");
            oss << "h1_rf_all_kEff_sum1_" << it;
            h1_rf_all_kEff_sum1[it]->SetName(oss.str().c_str());
// //-------------------------------------------------------------------------------------//
            h1_tmp_ave_all_kEff[it] = (TH1D*)h1_tmp->Clone();
//             h1_rf_all_kEff_sum2[it] = (TH1D*)h1_tmp->Clone();
//             oss.str("");
//             oss << "h1_rf_all_kEff_sum2_" << it << "_" << ix;
//             h1_rf_all_kEff_sum2[it]->SetName(oss.str().c_str());
//-------------------------------------------------------------------------------------//
         }
         else{
            h1_rf_all_kEff_sum1[it]->Add(h1_rf_all_kEff[it][irs]);
            h1_tmp_ave_all_kEff[it]->Add(h1_tmp                    );
         }
//         delete h1_tmp;
      }

      if( it != 2 ){
//         h1_tmp_ave_all_kEff[it] = (TH1D*)h1_rf_all_kEff_sum2[it]->Clone();
         h1_tmp_ave_all_kEff[it]->Add(h1_tmp_ave_all_kEff[2]);
      }

      h1_rf_all_kEff_sum1[it]->Scale(1/getPoT(rs, it));
//      h1_rf_all_kEff_sum2[it]->Divide(h1_nim3_sum[it]);
//         h1_rf_all_kEff_sum2[it]->Scale(h1_rf_all_kEff_sum1[it]->Integral() / h1_rf_all_kEff_sum2[it]->Integral());
   }
   h1_rf_all_kEff_sum1[1]->Add(h1_rf_all_kEff_sum1[2], -1);
   h1_rf_all_kEff_sum1[3]->Add(h1_rf_all_kEff_sum1[2], -1);
   // h1_rf_all_kEff_sum2[1]->Add(h1_rf_all_kEff_sum1[2], -1);
   // h1_rf_all_kEff_sum2[3]->Add(h1_rf_all_kEff_sum1[2], -1);

   // h1_rf_all_kEff_sum2[1]->Scale(h1_rf_all_kEff_sum1[1]->Integral() / h1_rf_all_kEff_sum2[1]->Integral());
   // h1_rf_all_kEff_sum2[3]->Scale(h1_rf_all_kEff_sum1[3]->Integral() / h1_rf_all_kEff_sum2[3]->Integral());
/////////// ALL X2 //////////

/////////// INTENSITY AVERAGE //////////
   for( int it = 1 ; it <= 3 ; it++ ){
      for( int ix = 0 ; ix <  8 ; ix++ ){
         if( it == 2 ) continue;
         h1_rf_average     [it][ix]->Add(h1_rf_average     [2][ix]);
         h1_rf_average_kEff[it][ix]->Add(h1_rf_average_kEff[2][ix]);

         h1_rf_average     [it][ix]->Divide(h1_tmp_ave     [it][ix]);
         h1_rf_average_kEff[it][ix]->Divide(h1_tmp_ave_kEff[it][ix]);
                  
      }
      if( it == 2 ) continue;
      h1_rf_all_average_kEff[it]->Add(h1_rf_all_average_kEff[ 2]);
      h1_rf_all_average_kEff[it]->Divide(h1_tmp_ave_all_kEff[it]);
   }
/////////// INTENSITY AVERAGE //////////
   for( int it = 1 ; it <= 3 ; it++ ){
      for( int ix = 0 ; ix <  8 ; ix++ ){
         if( it == 2 ) continue;
         oss.str("");
         oss << "tge_rf_" << it << "_" << ix;
         tge_rf     [it][ix] = convTH12TgeAve(h1_rf_x2_sum1     [it][ix], h1_rf_average     [it][ix], oss.str().c_str());
         oss.str("");
         oss << "tge_rf_kEff_" << it << "_" << ix;
         tge_rf_kEff[it][ix] = convTH12TgeAve(h1_rf_x2_kEff_sum1[it][ix], h1_rf_average_kEff[it][ix], oss.str().c_str());
      }
      oss.str("");
      oss << "tge_rf_all_kEff_" << it;
      tge_rf_all_kEff[it] = convTH12TgeAve(h1_rf_all_kEff_sum1[it], h1_rf_all_average_kEff[it], oss.str().c_str());
   }

}

////////////////
//    DRAW    //
////////////////
void writeAverage(){
   ofstream ofs;
   ofs.open("x2_average.txt");
   for( int ib = 1; ib <= h1_average_x2->GetNbinsX() ; ib++ ){
      ofs << h1_average_x2->GetBinContent(ib) << "\t";
   }
   ofs << endl;
   ofs.close();
}

void printYields(){
   ostringstream oss;
   ofstream ofs;
   for( int ix = 0 ; ix < nBinsX2 ; ix++ ){
      oss.str("");
      oss << "table_yields_" << ix << ".txt";
      ofs.open(oss.str().c_str());
      for( int ib = 1 ; ib <= h1_raw_tgt[1][ix]->GetNbinsX() ; ib++ ){
         ofs << h1_raw_tgt[1][ix]->GetBinLowEdge(ib)   << "\t" 
             << h1_raw_tgt[1][ix]->GetBinLowEdge(ib+1) << "\t";
         for( int it = 1 ; it <= 3 ; it++ )
            ofs << h1_raw_tgt[it][ix]->GetBinContent(ib) << "\t";
         ofs << endl;
      }
      ofs.close();
   }
}

void drawTGraphErrors(TGraphErrors* tge, const char* title, const char* name, double maxRange=binRange){
   ostringstream oss;

   TCanvas* c1 = new TCanvas("c1", "", 800, 600);
   c1->SetGrid();

   double x, y;
   double max = -100;
   for( int ip = 0 ; ip < tge->GetN() ; ip++ ){
      tge->GetPoint(ip, x, y);
      if( max < y ) max = y;
   }

   TH1* frame = c1->DrawFrame(0, 0, maxRange, max * 1.05);
   frame->SetTitle(title);
   frame->Draw();

   tge->SetMarkerColor(2);
   tge->SetMarkerStyle(20);
   tge->SetMarkerSize(1.2);
   tge->Draw("SAME P");

   if( strcmp(name, "") == 0 )
      name = tge->GetName();

   oss.str("");
   oss << name << "." << saveType;  
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

void drawTH1(TH1D* h1, const char* mode="HIST"){
   ostringstream oss;
   TCanvas* c1 = new TCanvas("c1", "", 800, 600);
   c1->SetGrid();

   h1->GetYaxis()->SetRangeUser(0, h1->GetMaximum() * 1.05);

   h1->Draw(mode);

   oss.str("");
   oss << h1->GetName() << "." << saveType;

   c1->SaveAs(oss.str().c_str());
   delete c1;
}

void drawTwoTH1s(TH1D* h1_1, TH1D* h1_2, const char* label_1, const char* label_2, const char* title, const char* name){
   ostringstream oss;
   TCanvas* c1 = new TCanvas("c1", "", 800, 600);
   c1->SetGrid();

   TH1D* tmp_1 = (TH1D*) h1_1->Clone();
   TH1D* tmp_2 = (TH1D*) h1_2->Clone();

   tmp_1->SetMarkerColor(2);
   tmp_1->SetLineColor  (2);
   tmp_1->SetMarkerStyle(20);
   tmp_1->SetMarkerSize(1.2);
   tmp_2->SetMarkerColor(4);
   tmp_2->SetLineColor  (4);
   tmp_2->SetMarkerStyle(21);
   tmp_2->SetMarkerSize(1.2);

   tmp_1->SetTitle(title);
   tmp_1->Draw("PE");
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

void drawMain(){
   ostringstream oss;
   ostringstream title;
   gSystem->mkdir(outDir.c_str(), true);
   gSystem->cd   (outDir.c_str()      );

   for( int ix = 0 ; ix < 8 ; ix++ ){
      for( int irs = 0 ; irs < (int)fileList.size() ; irs++ ){
         drawTH1(h1_rf_x2[1][ix][irs]);
         drawTH1(h1_rf_x2[3][ix][irs]);

         drawTH1(h1_rf_x2_kEff[1][ix][irs]);
         drawTH1(h1_rf_x2_kEff[3][ix][irs]);
      }
      drawTH1(h1_rf_x2_sum1[1][ix], "PE");
      drawTH1(h1_rf_x2_sum1[3][ix], "PE");
      drawTH1(h1_rf_x2_sum2[1][ix], "PE");
      drawTH1(h1_rf_x2_sum2[3][ix], "PE");

      drawTH1(h1_rf_x2_kEff_sum1[1][ix], "PE");
      drawTH1(h1_rf_x2_kEff_sum1[3][ix], "PE");
      drawTH1(h1_rf_x2_kEff_sum2[1][ix], "PE");
      drawTH1(h1_rf_x2_kEff_sum2[3][ix], "PE");

      drawTGraphErrors(tge_rf     [1][ix], ";Intensity;Yield;", "");
      drawTGraphErrors(tge_rf     [3][ix], ";Intensity;Yield;", "");
      drawTGraphErrors(tge_rf_kEff[1][ix], ";Intensity;Yield;", "");
      drawTGraphErrors(tge_rf_kEff[3][ix], ";Intensity;Yield;", "");

      oss.str("");
      oss << "h1_rf_x2_sum_comp_1_" << ix;
      drawTwoTH1s(h1_rf_x2_sum1[1][ix], h1_rf_x2_sum2[1][ix], "Sum v1", "Sum v2", 
                  ";Intensity;Yield", oss.str().c_str());

      oss.str("");
      oss << "h1_rf_x2_sum_comp_3_" << ix;
      drawTwoTH1s(h1_rf_x2_sum1[3][ix], h1_rf_x2_sum2[3][ix], "Sum v1", "Sum v2", 
                  ";Intensity;Yield", oss.str().c_str());


      oss.str("");
      oss << "h1_rf_x2_kEff_sum_comp_1_" << ix;
      drawTwoTH1s(h1_rf_x2_kEff_sum1[1][ix], h1_rf_x2_kEff_sum2[1][ix], "Sum v1", "Sum v2", 
                  ";Intensity;Yield", oss.str().c_str());

      oss.str("");
      oss << "h1_rf_x2_kEff_sum_comp_3_" << ix;
      drawTwoTH1s(h1_rf_x2_kEff_sum1[3][ix], h1_rf_x2_kEff_sum2[3][ix], "Sum v1", "Sum v2", 
                  ";Intensity;Yield", oss.str().c_str());
   }

   drawTGraphErrors(tge_rf_all_kEff[1], ";Intensity;Yield;", "");
   drawTGraphErrors(tge_rf_all_kEff[3], ";Intensity;Yield;", "");

//////////////////////////////////////////////////////////
   TCanvas* c1 = new TCanvas("c1", "", 800, 600);
   c1->SetGrid();
   TH1* h1_frame;
   TF1* f1[4][8];
   gStyle->SetOptFit(1111);
   TGraphErrors* tge_ratio = new TGraphErrors();
   tge_ratio->SetName("tge_ratio");
   for( int ix = 0 ; ix < nBinsX2 ; ix++ ){
      for( int it = 1 ; it < 4;  it++){
         if( it == 2 ) continue;
         double x, y;
         double max = -100;
         for( int ip = 0 ; ip < tge_rf_kEff[it][ix]->GetN() ; ip++ ){
            tge_rf_kEff[it][ix]->GetPoint(ip, x, y);
            if( max < y ) max = y;
         }
         h1_frame = c1->DrawFrame(0, 0, binRange, max * 1.05);
         oss.str("");
         oss << x2_bin[ix] << " < x_{2} < " << x2_bin[ix+1] << ";Intensity;Normalized Yield";
         h1_frame->SetTitle(oss.str().c_str());
         h1_frame->Draw();
         tge_rf_kEff[it][ix]->Draw("SAME P");
      
         oss.str("");
         oss << "f1_" << it << "_" << ix;
         f1[it][ix] = new TF1(oss.str().c_str(), "[0] + [1]*x + [2]*x*x", 0, binRange);
         f1[it][ix]->SetParameter(0, 0.1);
         f1[it][ix]->SetParameter(1, 0.1);
         f1[it][ix]->SetParameter(2, 0.1);
//         f1[it][ix]->SetParameter(3, 0.1);
//         f1[it][ix]->FixParameter(0, 0.0);

         tge_rf_kEff[it][ix]->Fit(f1[it][ix], "EX0", "P", 0, fitMax);

         c1->Modified();
         c1->Update();
         TPaveStats* stats = (TPaveStats*) c1->GetPrimitive("stats");
         stats->SetName("h1_stats");
         stats->SetY1NDC(.7);
         stats->SetY2NDC(.9);
         double width = fabs( stats->GetX1NDC() - stats->GetX2NDC() );
         stats->SetX1NDC(.1);
         stats->SetX2NDC(.1+width);
         stats->SetTextColor(1);
         c1->Update();

         oss.str("");
         oss << "fit_" << tge_rf_kEff[it][ix]->GetName() << "." << saveType;  
         c1->SaveAs(oss.str().c_str());
         delete stats;
      }
      double ratio = ( f1[3][ix]->GetParameter(1) * effAtoms_h / f1[1][ix]->GetParameter(1) / effAtoms_d + purity - 1 ) / purity / 2.;
      tge_ratio->SetPoint(ix, ( x2_bin[ix] + x2_bin[ix+1] )/2., ratio);
   }
   gStyle->SetOptFit(0);
//////////////////////////////////////////////////////////

   h2_qiesum_g2sem->Draw("COLZ");
   c1->SaveAs("h2_qiesum_g2sem.png");

   delete c1;

   drawTGraphErrors(tge_ratio, "Ratio;x_{2};Ratio", "tge_ratio", 0.6 );
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

   adjustment = 0.0;
   
   initMain(argc, argv);

   for( int ir = 0 ; ir < (int)fileList.size() ; ir++ ){
      if( !initInputFile (ir) ) continue;
      rs_current = ir;
      extMain ();
      deleteMemory();
   }

   for( int ir = 0 ; ir < (int)fileListNim3.size() ; ir++ ){
      if( !initInputFileNim3 (ir) ) continue;
      rs_current = ir;
      extMain (true);
      deleteMemory();
   }

   anaMain ();
   drawMain();      

//   cout << purity << " | " << effAtoms_h << " | " << effAtoms_d << endl;
}
