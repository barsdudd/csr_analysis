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
#include <TLegend.h>
#include <TLatex.h>
#include "Event.h"
#include "Variables.h"
#include "Selection.h"

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
int        rs;
string outDir;

///// DATA /////
TH1D* h1_rf_x2   [4][8]; // [target type (dummy, LH2, Empty, LD2)][x2 bins]
TGraphErrors* tge   [8]; // [x2 bins]

TH1D* h1_average[8]; // [x2]
TH1D* h1_raw    [8]; // [x2]
TH1D* h1_average_x2;
TH1D* h1_raw_x2    ;
TH1D* h1_average_x2_tgt[8]; // targetType
TH1D* h1_raw_x2_tgt    [8]; // targetType

vector<TH1D*> vec_th1s;

///// FIXED VARIABLES /////
int color[3] = {kBlack, kRed, kBlue};
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
bool    caCuts;

int bugFinder = 0;

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////

//////////////////////
//    INITIALIZE    //
//////////////////////
void initFiles(){
   ostringstream oss;
   oss.str("");
   oss << fileDir << rs << "/" << rs << "_with_cuts.root";
   cout << "FILE: " << oss.str() << endl;
   file  = new TFile(oss.str().c_str(), "READ");
   tree = (TTree*)file->Get("save");
   tree->SetBranchAddress("event", &event);
}

void initTH1(TH1D*& h1, const char* name, const char* title, int nBin, double min, double max, int color){
   h1 = new TH1D(name, title, nBin, min, max);
   h1->Sumw2();
   h1->SetStats(0);
   h1->SetLineColor(color);
   h1->SetLineWidth(2);
   vec_th1s.push_back(h1);
}

void initTH1(){
   ostringstream oss;
   ostringstream title;

   for( int ix = 0 ; ix < 8 ; ix++ ){
      for( int it = 1 ; it < 6 ; it++ ){
         oss.str("");
         oss << "h1_rf_" << it << "_" << ix;
         title.str("");
         title << "Trigger Intensity Ratio Plot (" << x2_bin[ix] << " < x_{T} < " << x2_bin[ix+1] 
               << ";Trigger Intensity;Counts";
         initTH1(h1_rf_x2[it][ix], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, color[it]);
      }

      oss.str("");
      oss << "h1_average_" << ix;
      title.str("");
      title << "Average (" << x2_bin[ix] << " < x_{2} < " << x2_bin[ix+1] 
            << ");Trigger Intensity;Counts";
      initTH1(h1_average[ix], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, 2);

      oss.str("");
      oss << "h1_raw_" << ix;
      title.str("");
      title << "raw (" << x2_bin[ix] << " < x_{2} < " << x2_bin[ix+1] 
            << ");Trigger Intensity;Counts";
      initTH1(h1_raw[ix], oss.str().c_str(), title.str().c_str(), nbins, 0, binRange, 2);

   }

   for( int it = 1 ; it < 6 ; it++ ){
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

}

void initMain(int argc, char* argv[]){
   ostringstream oss;

   if( argc < 5 ){
      cout << "!! NEED 5 ARGUMENTS !!: ROADSET NBINS FITRANGE CUTOFF-MASS [CA cut apply (CA or others), optional]" << endl;
      exit(0);
   }

   rs      = atoi(argv[1]);
   nbins   = atoi(argv[2]);
   fitMax  = atof(argv[3]);
   cutMass = atof(argv[4]);

   if( argc >= 5 ) caCuts  = strcmp(argv[5], "CA") == 0 ? true : false;

   event = new Event();

   const char* caApplied = caCuts ? "app" : "not";

   oss.str("");
   oss << "results_" << rs;
   oss << "/nbins_" << nbins << "_range_" << fitMax << "_cutoffMass_" << cutMass << "_CA_" << caApplied;
   outDir = oss.str();

   initFiles();
   initTH1();

   purity = getPurity(rs);
   effAtoms_h = effAtoms(aMass_h, rho_h, lambda_h);
   effAtoms_d = effAtoms(effLinear(aMass_d, aMass_h, rs), effLinear(rho_d, rho_h, rs), effInverse(lambda_d, lambda_h, rs));

   cout << effAtoms_h << " | " << effAtoms_d << endl;

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

void extMain(){
   int prev = -1;
   int curr =  0;
   for( int ie = 0 ; ie < tree->GetEntries() ; ++ie){    
      tree->GetEntry(ie);
      curr = (ie+1) * 100 / tree->GetEntries();
      if( curr != prev ){
         cout << "\r" << (ie+1) << " / " << tree->GetEntries() << " = " << curr << " %" << flush;
         prev = curr;
      }
      for( int id = 0 ; id < (int)event->dimuons.size() ; id++ ){
         Dimuon dimuon = event->dimuons[id];
         Track  trackP = event-> tracks[findTrackIndex(dimuon.trackID_pos)];
         Track  trackN = event-> tracks[findTrackIndex(dimuon.trackID_neg)];
         if( dimuon.mass <= cutMass ) continue;
         if( caCuts && !CAisSatisfied_2111_v32(trackP, trackN, rs) ) continue;
         int ix = getIX2(dimuon.x2);
         if( ix == -1 ) continue;
         if( event->targetPos > 5 ) continue;
         h1_rf_x2[event->targetPos][ix]->Fill(event->inte_t[16]);

         cout << event->inte_t[16] << endl;

         if( bugFinder > 100 ) exit(1);
         bugFinder++;

         h1_average[ix]->Fill(event->inte_t[16], event->inte_t[16]);
         h1_raw    [ix]->Fill(event->inte_t[16]);

         if( event->targetPos == 5 ) cout << event->targetPos << " " << ix << endl;

         if( ++bugFinder > 100 ) break;

         if( fitMax > event->inte_t[16] ){
            h1_average_x2_tgt[event->targetPos]->Fill(dimuon.x2, dimuon.x2);
            h1_raw_x2_tgt    [event->targetPos]->Fill(dimuon.x2);
            if( event->targetPos != 2 ){
               h1_average_x2->Fill(dimuon.x2, dimuon.x2);
               h1_raw_x2    ->Fill(dimuon.x2);
            }
         }
      }
   }
   cout << endl;
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

      tge->SetPoint     ( tge->GetN()  , h1_ave->GetBinContent(ib), ratio );
      tge->SetPointError( tge->GetN()-1,                         0, error );
      tge->SetMarkerColor(2);
      tge->SetMarkerSize(1.2);
      tge->SetMarkerStyle(20);
   }
   return tge;
}

void anaMain(){
   for( int ix = 0 ; ix < 8 ; ix++ ){
      normHist(h1_rf_x2[1][ix], 1);
      normHist(h1_rf_x2[2][ix], 2);
      normHist(h1_rf_x2[3][ix], 3);
      normHist(h1_rf_x2[4][ix], 4);
      normHist(h1_rf_x2[5][ix], 5);
      h1_average[ix]->Divide(h1_raw[ix]);
      tge[ix] = getRatio(h1_rf_x2[3][ix], h1_rf_x2[1][ix], h1_rf_x2[2][ix], h1_average[ix]);
   }
   h1_average_x2_tgt[1]->Divide(h1_raw_x2_tgt[1]);
   h1_average_x2_tgt[3]->Divide(h1_raw_x2_tgt[3]);

   h1_average_x2->Divide(h1_raw_x2);
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

void drawTH1s(){
   ostringstream oss;
   TCanvas* c1 = new TCanvas("c1", "", 800, 600);
   c1->SetGrid();
   for( int ith = 0 ; ith < (int)vec_th1s.size() ; ith++ ){
      cout << vec_th1s[ith]->GetBinContent(1) << endl;
      TH1D* h1 = (TH1D*)vec_th1s[ith]->Clone();
      h1->SetLineColor(2);
      h1->SetLineWidth(2);
      h1->Draw("HIST");
      oss.str("");
      oss << vec_th1s[ith]->GetName() << "." << saveType;
      c1->SaveAs(oss.str().c_str());
   }
   delete c1;
}

void drawTGraphErrors(TGraphErrors* tge, const char* title, const char* name){
   ostringstream oss;

   TCanvas* c1 = new TCanvas("c1", "", 800, 600);
   c1->SetGrid();

   TH1* frame = c1->DrawFrame(0, 0, binRange, 5);
   frame->SetTitle(title);
   frame->Draw();

   tge->Draw("SAME P");

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
   ostringstream oss;

   TCanvas* c1 = new TCanvas("c1", "", 800, 600);
   c1->SetGrid();

   TH1* frame = c1->DrawFrame(0, 0, binRange, 5);
   frame->SetTitle(title);
   frame->Draw();

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

void drawMain(){
   ostringstream oss;
   ostringstream title;
   gSystem->mkdir(outDir.c_str(), true);
   gSystem->cd   (outDir.c_str()      );

   drawTH1s();

   writeAverage();

   drawTwoTH1s(h1_average_x2_tgt[1], h1_average_x2_tgt[3], "LH_{2}", "LD_{2}", 
               "Average x_{2};x_{2};Average x_{2}", "h1_average_x2_comp");

   ofstream ofs_lin;
   ofs_lin.open("fit_results_linear.txt");
   ofstream ofs_qua;
   ofs_qua.open("fit_results_quadra.txt");

   ofstream ofs_full_lin;
   ofs_full_lin.open("fit_results_full_linear.txt");
   ofstream ofs_full_qua;
   ofs_full_qua.open("fit_results_full_quadra.txt");

   double parameters[2];
   double errors    [2];
   double chisq        ;
   int    ndf          ;

   double parameters_q[3];
   double errors_q    [3];
   for( int ix = 0 ;ix < 8 ; ix++ ){
      parameters[0] = 2;
      parameters[1] = 0.001;
      oss.str("");
      oss << "tge_ratio_" << ix;
      title.str("");
      title << "Trigger Intensity Ratio (" << x2_bin[ix] << " < x < " << x2_bin[ix+1] << ")"
            << ";Trigger Intensity;(LD_{2} - Empty)/(LH_{2} - Empty)";
      drawTGraphErrors(tge[ix], title.str().c_str(), oss.str().c_str());
      oss.str("");
      oss << "lin_tge_ratio_" << ix;
      fitTGraphErrors (tge[ix], title.str().c_str(), "[0]+[1]*x", parameters, errors, chisq, ndf, oss.str().c_str(), fitMax);
      ofs_lin << parameters[0] << "\t" << errors[0] << "\t"
//              << parameters[1] << "\t" << errors[1] << "\t"
              << chisq         << "\t" << ndf       << "\t"
              << chisq              /     ndf       << endl;
      ofs_full_lin << parameters[0] << "\t" << errors[0] << "\t"
                   << parameters[1] << "\t" << errors[1] << "\t"
                   << chisq         << "\t" << ndf       << "\t"
                   << chisq              /     ndf       << endl;

      parameters_q[0] = 2;
      parameters_q[1] = 0.001;
      parameters_q[2] = 0.001;
      oss.str("");
      oss << "qua_tge_ratio_" << ix;
      fitTGraphErrors(tge[ix], title.str().c_str(), "[0]+[1]*x+[2]*x*x", parameters_q, errors_q, chisq, ndf, oss.str().c_str(), fitMax);
      ofs_qua << parameters_q[0] << "\t" << errors_q[0] << "\t"
              // << parameters_q[1] << "\t" << errors_q[1] << "\t"
              // << parameters_q[2] << "\t" << errors_q[2] << "\t"
              << chisq           << "\t" << ndf         << "\t"
              << chisq                /     ndf         << endl;
      ofs_full_qua << parameters_q[0] << "\t" << errors_q[0] << "\t"
                   << parameters_q[1] << "\t" << errors_q[1] << "\t"
//                   << parameters_q[2] << "\t" << errors_q[2] << "\t"
                   << chisq           << "\t" << ndf         << "\t"
                   << chisq                /     ndf         << endl;
   }
   ofs_lin.close();
   ofs_qua.close();
   ofs_full_lin.close();
   ofs_full_qua.close();
}

////////////////
//    MAIN    //
////////////////
int main(int argc, char* argv[]){
   gROOT ->SetStyle("Plain");
   gStyle->SetOptStat(1);
   gStyle->SetPadGridX(kTRUE);
   gStyle->SetPadGridY(kTRUE);
   
   initMain(argc, argv);
   extMain ();
   anaMain ();
   drawMain();      
}
