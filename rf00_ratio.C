#include <map>

double x2_bin[9] = {0.1, 0.13, 0.16, 0.20, 0.24, 0.29, 0.35, 0.45, 0.58};

int color[3] = {kBlack, kRed, kBlue};
const char* targetType[4] = {"dummy_variable","LH_{2}", "Empty", "LD_{2}"};
//double PoT_data[3] = {5.99079e17, 5.55621e17, 3.66565e16};
double PoT_data[4] = {0, 1.55804e17, 3.66565e16, 7.2706e16};

const char* rf_type     [4] = {"D1", "D2", "D3", "D1+D2+D3"};
const char* rf_type_name[4] = {"D1", "D2", "D3", "sum"};
double turns = 369000;
//double pedestal = 36.2;
double buckets = 588;
int binsize = 60000;
int    nbins = 15;
double  factor = 2.3105;

int nBinsrf = 300;
double maxrf = 90000;

const char* filePlace = "/seaquest/users/arunts/data/";

TFile* file;
TTree* tree;
int D1, D2, D3, targetPos;
double mass, xT, xF, costh, RF00, G2SEM, QIEsum;
double x1_mu_p, x1_mu_m, pz1_kTrack_mu_p, pz1_kTrack_mu_m, dz, Intensity_p;

//TH1D* h1_rf[4][3]; //]
TH1D* h1_rf_x2   [4][8]; // [target type (dummy, LH2, Empty, LD2)][x2 bins]
TGraphErrors* tge    [8]; // [x2 bins]

using namespace std;

map<int, double> purity;
map<int, double> pedestal;

int getIX2(double x2){
   if( x2 <= x2_bin[0] || x2 > x2_bin[8] )
      return -1;
   int ix = -1;
   while( x2 > x2_bin[ix+1] ) ix++;
   return ix;
}

void initFiles(int rs){
   ostringstream oss;
   oss.str("");
   oss << filePlace << rs << "/" << rs << "_with_cuts.root";
   file  = new TFile(oss.str().c_str(), "READ");
   oss.str("");
   oss << "roadset_" << rs;
   tree = (TTree*)file->Get(oss.str().c_str());
   tree->SetBranchAddress("mass",&mass);
   tree->SetBranchAddress("xT",&xT);
   tree->SetBranchAddress("D1",&D1);
   tree->SetBranchAddress("D2",&D2);
   tree->SetBranchAddress("D3",&D3);
   tree->SetBranchAddress("xF",&xF);
   tree->SetBranchAddress("dz",&dz);
   tree->SetBranchAddress("targetPos",&targetPos);
   tree->SetBranchAddress("RF00",&RF00);
   tree->SetBranchAddress("G2SEM",&G2SEM);
   tree->SetBranchAddress("Intensity_p",&Intensity_p);
   tree->SetBranchAddress("QIEsum",&QIEsum);
   tree->SetBranchAddress("costh",&costh);
   tree->SetBranchAddress("x1_mu_p",&x1_mu_p);
   tree->SetBranchAddress("x1_mu_m",&x1_mu_m);
   tree->SetBranchAddress("pz1_kTrack_mu_p",&pz1_kTrack_mu_p);
   tree->SetBranchAddress("pz1_kTrack_mu_m",&pz1_kTrack_mu_m);
}

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
   for( int it = 1 ; it < 4 ; it++ ){
//      for( int io = 0 ; io < 4 ; io++ ){
//         oss.str("");
//         oss << "h1_rf_" << rf_type_name[io] << "_" << it;
//         title.str("");
//         title << rf_type[io] << " RF Plot;" << rf_type[io] << ";Counts";
//         initTH1(h1_rf[io][it], oss.str().c_str(), title.str().c_str(), nBinsrf, 0, maxrf, color[it]);
//      }
      
      for( int ix = 0 ; ix < 8 ; ix++ ){ // "ix" = "ith x2 bin"
         oss.str("");
         oss << "h1_rf_" << it << "_" << ix;
         title.str("");
         title << "Trigger Intensity Ratio Plot (" << x2_bin[ix] << " < x_{T} < " << x2_bin[ix+1] << ";Trigger Intensity;Counts";
         initTH1(h1_rf_x2[it][ix], oss.str().c_str(), title.str().c_str(), nbins, 0, binsize, color[it]);
      }
   }
}

void initMain(int rs){
   initFiles(rs);
   initTH1();
   
   purity[57] = (100-9.6) / 100.;
   purity[59] = (100-9.6) / 100.;
   purity[62] = (100-8.7) / 100.;
   purity[67] = (100-4.7) / 100.;
   purity[70] = 1.;
   
   pedestal[57] = 36.2;
   pedestal[59] = 36.2;
   pedestal[62] = 36.2;
   pedestal[67] = 32.6;
   pedestal[70] = 32.6;
   
}

void anaMain(int rs){
   for( int ie = 0 ; ie < tree->GetEntries() ; ++ie){    
      tree->GetEntry(ie);
      if( mass <= 4.2) continue;         
      int ix = getIX2(xT);
      if( ix == -1 ) continue;
      if(targetPos>3) continue;
      h1_rf_x2[targetPos][ix]->Fill((RF00 - pedestal[rs])*G2SEM/(QIEsum - turns*pedestal[rs]*buckets));
   }   
}

void rf00_ratio(int rs){
   gROOT->SetStyle("Plain");
   gStyle->SetOptStat(1);
   gStyle->SetPadGridX(kTRUE);
   gStyle->SetPadGridY(kTRUE);
   
   initMain(rs);
   anaMain(rs);
   
   ostringstream oss;
   
   oss.str("");
   oss << "results_" << rs;
   gSystem->mkdir(oss.str().c_str(), true);
   gSystem->cd(oss.str().c_str());
   
   TCanvas *c_x2[2];
   c_x2[0] = new TCanvas("c_x2_0","multipads", 1000, 800);
   c_x2[1] = new TCanvas("c_x2_1","multipads", 1000, 800);
   c_x2[0]->Divide(2,2);
   c_x2[1]->Divide(2,2);
   
   ofstream table_intercepts;
   table_intercepts.open("table_intercepts_errors_rf.txt");
   table_intercepts<<"p0, p0_err, p1, p1_err"<<"\n";
   
   for( int ix = 0 ; ix < 8 ; ix++ ){
      int iCanvas = ix / 4;
      int iPad    = ix % 4 + 1;
      TPad* pad = c_x2[iCanvas]->cd(iPad);
      
      tge[ix] = new TGraphErrors();
      tge[ix]->SetMarkerColor(2);
      tge[ix]->SetMarkerSize(1.2);
      tge[ix]->SetMarkerStyle(20);
      
      //h1_occ_x2 is the histogram from FPGA1
      h1_rf_x2[1][ix]->Scale(1/PoT_data[1]);
      h1_rf_x2[2][ix]->Scale(1/PoT_data[2]);
      h1_rf_x2[3][ix]->Scale(1/PoT_data[3]);
      
      for (int iBin = 1 ; iBin <= h1_rf_x2[1][ix]->GetNbinsX() ; iBin++ ){
         double content_lh2 = h1_rf_x2[1][ix]->GetBinContent(iBin);
         double content_emp = h1_rf_x2[2][ix]->GetBinContent(iBin);
         double content_ld2 = h1_rf_x2[3][ix]->GetBinContent(iBin);
         
         if(content_lh2 - content_emp == 0)
            continue;
         double error_lh2 = h1_rf_x2[1][ix]->GetBinError(iBin);
         double error_emp = h1_rf_x2[2][ix]->GetBinError(iBin);
         double error_ld2 = h1_rf_x2[3][ix]->GetBinError(iBin);
         
         double error_ratio = sqrt((1/(content_lh2 - content_emp))*(1/(content_lh2 - content_emp))*error_ld2*error_ld2 + (content_ld2 - content_emp)*(content_ld2 - content_emp)*(1/(content_lh2 - content_emp))*(1/(content_lh2 - content_emp))*(1/(content_lh2 - content_emp))*(1/(content_lh2 - content_emp))*error_lh2*error_lh2 + (content_ld2 - content_lh2)*(content_ld2 - content_lh2)*(1/(content_lh2 - content_emp))*(1/(content_lh2 - content_emp))*(1/(content_lh2 - content_emp))*(1/(content_lh2 - content_emp))*error_emp*error_emp)/factor;
         double ratio = ((content_ld2 - content_emp)/((content_lh2 - content_emp)*(factor)) + (purity[rs] - 1)/2. ) / purity[rs];
         tge[ix]->SetPoint     (tge[ix]->GetN()  , h1_rf_x2[1][ix]->GetBinCenter(iBin), ratio);
         tge[ix]->SetPointError(tge[ix]->GetN()-1, 0, error_ratio);
      }
            //TH1* frame = pad->DrawFrame(0, 0, 300, 5);
      TH1* frame = pad->DrawFrame(0, 0, binsize, 5);

      oss.str("");
      oss << "Trigger Intensity Ratio (" << x2_bin[ix] << " < x < " << x2_bin[ix+1] << ")"
          << ";Trigger Intensity (RF00);(LD_{2} - Empty)/(LH_{2} - Empty)";
      frame->SetTitle(oss.str().c_str());
      frame->Draw();
      tge[ix]->Draw("SAME P");
      leg = new TLegend (0.5, 0.77, 0.99, 0.9);
      leg->SetTextFont(32);
      leg->SetTextSize(0.045);
      leg->AddEntry(tge[ix], "LD_{2} - empty/LH_{2} - empty", "pl");
      leg-> Draw();
      TF1 * f1 = new TF1("f1", "[0] + [1]*(x) ", 0, binsize);
      f1->SetParameter(0,2);
      f1->SetParameter(1,0.001);
      tge[ix]->Fit("f1");
      
      double p0      = f1->GetParameter(0);
      double p0_err  = f1->GetParError(0);
      double p1      = f1->GetParameter(1);
      double p1_err  = f1->GetParError(1);
      double chisq   = f1->GetChisquare();
      double ndf     = f1->GetNDF();
      double rdchisq = chisq/ndf;
      
      double xpos = 0.12;     
      TLatex * text = new TLatex();
      oss.str("");
      oss << "chisq: " << chisq;
      text->DrawLatexNDC(xpos, 0.85,oss.str().c_str());
      oss.str("");
      oss << "NDF: " << ndf;
      text->DrawLatexNDC(xpos, 0.80,oss.str().c_str());
      oss.str("");
      oss << "red_chisq: " << rdchisq;
      text->DrawLatexNDC(xpos, 0.75,oss.str().c_str());
      oss.str("");
      oss << "p0: " << p0;
      text->DrawLatexNDC(xpos, 0.70,oss.str().c_str());
      oss.str("");
      oss << "p0_err: " << p0_err;
      text->DrawLatexNDC(xpos, 0.65,oss.str().c_str());
      oss.str("");
      oss << "p1: " << p1;
      text->DrawLatexNDC(xpos, 0.60,oss.str().c_str());
      oss.str("");
      oss << "p1_err: " << p1_err;
      text->DrawLatexNDC(xpos, 0.55,oss.str().c_str());
      
      
      table_intercepts<< p0 <<"," << p0_err << "," << p1 << "," << p1_err << "\n";
      
      delete f1;     
      
      if( ix % 4 == 3 ){
         oss.str("");
         oss << "plots_rf00_ratio_" << iCanvas << ".png";
         c_x2[iCanvas]->SaveAs(oss.str().c_str());
      }
   }
   table_intercepts.close();
}
