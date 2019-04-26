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

using namespace std;

int        rs;
string outDir;

vector<TGraphErrors*> tge_csr;

vector<double> intercept [8];
vector<double> interError[8];
vector<double> slope     [8];
vector<double> slopeError[8];
vector<string> redChisq  [8];

vector<double> paraAerr[8][10]; // [ix][condition]

double nBins[] = {15, 20, 25, 30};
double range[] = {20e3, 30e3, 40e3, 50e3, 60e3};
string fitMode[] = {"linear", "quadra"};
string CA[] = {"app", "not"};
double mass[] = {4.2, 4.7, 5.0, 5.3};

string changeVal[] = {"nBins", "range", "mass", "CA", "fitMode"};

vector<double> x2Average[100];

int nNBins = sizeof nBins / sizeof nBins[0];
int nRange = sizeof range / sizeof range[0];
int nMode  = 2;
int nCA    = 2;
int nMass  = sizeof mass / sizeof mass[0];

void initOutput(){
   ostringstream oss;
   oss << "results_" << rs << "/compare";
   outDir = oss.str();
   gSystem->mkdir(outDir.c_str(), true);
}

void getX2Average(vector<double>& x2_average, double nBins, double range, double cutMass, const char* caApplied){
   if( (int)x2_average.size() > 0 ) 
      x2_average.clear();
   ostringstream oss;
   ifstream ifs;
   oss.str("");
   oss << "results_" << rs;
   oss << "/nbins_" << nBins << "_range_" << range << "_cutoffMass_" << cutMass << "_CA_" << caApplied << "/x2_average.txt";
   ifs.open(oss.str().c_str());
   double x2_ave;
   while( ifs >> x2_ave )
      x2_average.push_back(x2_ave);
}

void writeTable(const char* name){
   ostringstream oss;
   ofstream ofs;
   for( int ix = 0 ;ix < 8 ; ix++ ){
      oss.str("");
      oss << outDir << "/" << name << "_" << ix << ".txt";
      ofs.open(oss.str().c_str());
      for( int ic = 0 ; ic < (int)intercept [ix].size() ; ic++ ) ofs << intercept [ix][ic] << "\t";
      ofs << endl;
      for( int ic = 0 ; ic < (int)interError[ix].size() ; ic++ ) ofs << interError[ix][ic] << "\t";
      ofs << endl;
      for( int ic = 0 ; ic < (int)redChisq  [ix].size() ; ic++ ) ofs << redChisq  [ix][ic] << "\t";
      ofs << endl;
      ofs.close();
   }
}

void writeLinQuaFull(const char* name){
   ostringstream oss;
   ofstream ofs;
   oss.str("");
   oss << outDir << "/" << name << ".txt";
   ofs.open(oss.str().c_str());
   for( int ix = 0 ; ix < nBinsX2 ; ix++ ){
      ofs << x2_bin[ix] << "\t" << x2_bin[ix+1] << "\t";
      for( int in = 0 ; in < 10 ; in++ ){
         for( int id = 0 ; id < (int)paraAerr[ix][in].size() ; id++ )
            ofs << paraAerr[ix][in][id] << "\t";
      }
      ofs << endl;
   }
   ofs.close();
}

void drawCSR(int loop, int nConditions, const char* title, const char* fitModeF){

   if( nConditions != (int)intercept[0].size() ){
      cout << "!!! # of conditions and size of intercept are not the same !!! : " 
           << nConditions << " | " <<  (int)intercept[0].size()
           << endl;
      exit(1);
   }
   ostringstream oss;
   TCanvas* c1 = new TCanvas("c1", "", 800, 600);
   c1->SetGrid();
   TH1* frame = c1->DrawFrame(0, 0, 0.6, 2.5);
   oss.str("");
   oss << "CSR | " << title << "CSR";
   frame->SetTitle(oss.str().c_str());
   frame->Draw();

   oss.str("");
   oss << "different_" << changeVal[loop];
   if( loop != 4 ) oss << "_" << fitModeF;
   string name = oss.str();

   TLegend* leg = new TLegend(0.1, 0.6, 0.25, 0.9);
   TGraphErrors* tge[100];
   double width = 0.02;
   double step  = width/nConditions;
   double start = width/2.;
   for( int ic = 0 ; ic < nConditions ; ic++ ){
      tge[ic] = new TGraphErrors();
      tge[ic]->SetMarkerSize(1.2);
      tge[ic]->SetMarkerStyle(markerList[ic]);
      tge[ic]->SetMarkerColor(colorList[ic]);
      tge[ic+nConditions] = new TGraphErrors();
      tge[ic+nConditions]->SetMarkerSize(1.2);
      tge[ic+nConditions]->SetMarkerStyle(markerList[ic]);
      tge[ic+nConditions]->SetMarkerColor(colorList[ic]);
      double offset = start + step * ic;
      for( int ix = 0 ; ix < 8 ; ix++ ){

         tge[ic]->SetPoint     (ix, x2Average[ic][ix]+offset, intercept [ix][ic]);
         tge[ic]->SetPointError(ix,                        0, interError[ix][ic]);

         tge[ic+nConditions]->SetPoint     (ix, x2Average[ic][ix]+offset, slope     [ix][ic]);
         tge[ic+nConditions]->SetPointError(ix,                        0, slopeError[ix][ic]);
      }

      tge[ic]->Draw("SAME P");
      oss.str("");
      if( loop == 0 ) oss << nBins  [ic];
      if( loop == 1 ) oss << range  [ic];
      if( loop == 2 ) oss << mass   [ic];
      if( loop == 3 ) oss << CA     [ic];
      if( loop == 4 ) oss << fitMode[ic];
//      cout << "Loop: " << loop << " / ic: " << ic << " / value: " << oss.str().c_str() << endl;
      leg->AddEntry(tge[ic], oss.str().c_str(), "pl");
   }
   leg->Draw();
   oss.str("");
   oss << outDir << "/" << name << ".pdf";
   c1->SaveAs(oss.str().c_str());

   oss.str("");
   oss << "Slope | " << title << "Slope";
   frame->GetYaxis()->SetRangeUser(-1e-4, 1e-4);
   frame->SetTitle(oss.str().c_str());
   frame->Draw();
   for( int ic = 0 ; ic < nConditions ; ic++ )
      tge[ic+nConditions]->Draw("SAME P");
   leg->Draw();
   oss.str("");
   oss << outDir << "/slope_" << name << ".pdf";
   c1->SaveAs(oss.str().c_str());

   delete c1;
}


void pushBackParas(const char* fileName){
   ifstream ifs;
   ifs.open(fileName);
   double params[100];
   double errors[100];
   string redCh;
   int ix = 0;
   if( ifs.fail() ){
      cout << " !!! FAILED READING DATA !!! : " << fileName << endl;
      exit(1);
   }
   while( ifs >> params[0] >> errors[0] ){
      ifs >> params[1] >> errors[1];
      ifs >> redCh >> redCh >> redCh;
      intercept [ix].push_back(params[0]);
      interError[ix].push_back(errors[0]);
      slope     [ix].push_back(params[1]);
      slopeError[ix].push_back(errors[1]);
      redChisq  [ix].push_back(redCh);
      ix++;
   }   
   ifs.close();
}

const char* getFileName(double nBinF, double rangeF, double cutMassF, const char* caAppliedF, const char* fitModeF, 
                        int loop, int ic){
   ostringstream oss;
   oss.str("");
   oss << "results_" << rs << "/";

   if( loop == 0 ) oss << "nbins_" << nBins[ic] << "_";
   else            oss << "nbins_" << nBinF      << "_";

   if( loop == 1 ) oss << "range_" << range [ic] << "_";
   else            oss << "range_" << rangeF     << "_";

   if( loop == 2 ) oss << "cutoffMass_" << mass[ic] << "_";
   else            oss << "cutoffMass_" << cutMassF << "_";

   if( loop == 3 ) oss << "CA_" << CA[ic]    ;
   else            oss << "CA_" << caAppliedF;

   oss << "/fit_results_full_";

   if( loop == 4 ) oss << fitMode[ic];
   else            oss << fitModeF;

   oss << ".txt";

   return oss.str().c_str();
}


void collectData(double nBinF, double rangeF, double cutMassF, const char* caAppliedF, const char* fitModeF, 
                 int loop, int nConditions){
   ostringstream oss;
   ostringstream title;

   for( int ic = 0 ; ic < nConditions ; ic++ ){
      pushBackParas(    getFileName(nBinF, rangeF, cutMassF, caAppliedF, fitModeF, loop, ic) );
      getX2Average ( x2Average[ic], nBinF, rangeF, cutMassF, caAppliedF                      );
   }

   title.str("");
   title << "nBins = ";// << nBinF << " | FitMax = " << rangeF << " | Lin. & Qua. ;x_{2};";
   if( loop != 0 ) title << nBinF;
   else            title << "XXX";
   title << " | Fit Max = ";
   if( loop != 1 ) title << rangeF;
   else            title << "XXX";
   title << " | Mass Cut = ";
   if( loop != 2 ) title << cutMassF;
   else            title << "XXX";
   title << " | CA = ";
   if( loop != 3 ) title << caAppliedF;
   else            title << "XXX";
   title << " | Fit = ";
   if( loop != 4 ) title << fitModeF;
   else            title << "Lin./Qua";
   title << ";x_{2};";

   drawCSR   (loop, nConditions, title.str().c_str(), fitModeF);

   for( int ix = 0 ; ix < 8 ; ix++ ){
      intercept [ix].clear();
      interError[ix].clear();
      slope     [ix].clear();
      slopeError[ix].clear();
      redChisq  [ix].clear();
   }
}


int main(int argc, char* argv[]){
   gROOT->SetStyle("Plain");
   if( argc < 2 ){
      cout << "!! NEED 1 ARGUMENT !!: ROADSET" << endl;
      exit(0);
   }
   rs = atoi(argv[1]);
   initOutput();

   collectData(15, 60000, 4.2, "not", "linear", 0, nNBins);
   collectData(15, 60000, 4.2, "not", "quadra", 0, nNBins);

   collectData(15, 60000, 4.2, "not", "linear", 1, nRange);
   collectData(15, 60000, 4.2, "not", "quadra", 1, nRange);

   collectData(15, 60000, 4.2, "not", "linear", 2, nMass);
   collectData(15, 60000, 4.2, "not", "quadra", 2, nMass);

   collectData(15, 60000, 4.2, "not", "linear", 3, nCA);
   collectData(15, 60000, 4.2, "not", "quadra", 3, nCA);

   collectData(15, 60000, 4.2, "not", "linear", 4, nMode);

   cout << "output: " << gSystem->pwd() << "/" << outDir << endl;
}
