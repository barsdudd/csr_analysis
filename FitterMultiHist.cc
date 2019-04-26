#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TSystem.h>
#include <TH1.h>
#include "FitterMultiHist.h"
using namespace std;

// const double FitterMultiHist::X2_LOW[FitterMultiHist::N_HIST+1] = {
//   0.10, 0.13, 0.16, 0.20, 0.24, 0.29, 0.35, 0.45 // , 0.58
// };

// /// Taken from doc4450 v1 p3 (not exactly correct)
// const double FitterMultiHist::X2_AVG[FitterMultiHist::N_HIST] = {
//   0.119, 0.146, 0.179, 0.218, 0.262, 0.315, 0.385 // , 0.483
// };

FitterMultiHist::FitterMultiHist()
{
   m_minuit = 0;
   m_label  = "dummy";
   m_n_par  = 0; 
}

void FitterMultiHist::Init(const char* fn_hist, const char* gn, const char* avg_hist)
{
   cout << "FitterMultiHist::Init():\n" << gn << "\n";
   hist_name = gn;
   if (m_n_par == 0) {
      cerr << "ERROR: The number of parameters (m_n_par) hasn't been set, which must be set in the constructor of your class.  Abort." << endl;
      exit(1);

   }
   m_minuit = new TFitterMinuit(m_n_par);
   m_minuit->SetMinuitFCN(this);
   m_minuit->CreateMinimizer();
//  m_minuit->SetPrintLevel(0); // Comment out this line to debug the fitting result

   TFile* file = new TFile(fn_hist);
   if (! file) {
      cout << "  ERROR:  Failed at opening '" << fn_hist << "'.  Abort." << endl;
      exit(1);
   }
   ostringstream oss;
   h1_avg = (TH1D*)file->Get(avg_hist);
   N_HIST = h1_avg->GetNbinsX();
   for (int ih = 0; ih < N_HIST; ih++) {
      oss.str("");
      oss << gn << "_" << ih;
      m_list_gr[ih] = (TGraphErrors*)file->Get(oss.str().c_str());
      if (! m_list_gr[ih]) {
         cout << "  ERROR:  Failed at getting '" << oss.str() << "'.  Abort." << endl;
         exit(1);
      }
      cout << "  ";
      m_list_gr[ih]->Print();
   }
   cout << endl;

   InitPar();
}

void FitterMultiHist::ChangeDir()
{
   string dir_out = "result/" + Label();
   cout << "Change directory to '" << dir_out << "'." << endl;
   gSystem->mkdir(dir_out.c_str(), true);
   gSystem->cd   (dir_out.c_str());
}

void FitterMultiHist::DoFit()
{
   if (! m_minuit) {
      cerr << "ERROR: Called DoFit() before Init().  Abort." << endl;
      exit(1);
   }   

   const int nfcn = 1000000; // N of max function calls per minimization
   int ret = m_minuit->Minimize(nfcn);
   if (ret != 0) {
      cerr << "!!ERROR!!  FitterMultiHist::DoFit():  Bad fit result (ret = " << ret << ").\n"
           << "Please fix the problem.  Abort.\n";
      exit(1);
   }
   int n_par = m_minuit->GetNumberTotalParameters();
   m_pars.clear();
   m_errs.clear();
   for (int ii = 0; ii < n_par; ii++) {
      m_pars.push_back(m_minuit->GetParameter(ii));
      m_errs.push_back(m_minuit->GetParError (ii));
   }
   CalcChi2(m_pars, m_chi2, m_ndf);
   m_cov_mat = m_minuit->GetCovarianceMatrix();
}

void FitterMultiHist::PrintFit(std::string fname)
{
   ofstream ofs(fname.c_str());
   PrintFit(ofs);
   ofs.close();
}

void FitterMultiHist::PrintFit(std::ostream& os)
{
   int n_pars = m_pars.size();

   os << "\nFit parameters:\n";
   for (int i_par = 0; i_par < n_pars; i_par++) {
      os << "  par " << i_par << " :  " << m_minuit->GetParName(i_par) << "\n";
   }
   os << "\n";

   os << "Fit result:\n"
      << "  chi2 / ndf = " << m_chi2 << " / " << m_ndf << " = " <<  m_chi2 / m_ndf << "\n";
   for (int i_par = 0; i_par < n_pars; i_par++) {
      os << "  par " << i_par << " :  " << m_pars[i_par] << " +- " << m_errs[i_par] << "\n";
   }
}

void FitterMultiHist::DrawFit()
{
   cerr << "ERROR:  The virtual DrawFit() function was called.  Abort." << endl;
   exit(1);
}

void FitterMultiHist::InitPar()
{
   cerr << "ERROR:  The virtual InitPar() function was called.  Abort." << endl;
   exit(1);
}

void FitterMultiHist::CalcChi2(const std::vector<double>& pars, double& chi2, int& ndf) const
{
   cerr << "ERROR:  The virtual CalcChi2() function was called.  Abort." << endl;
   exit(1);
}

double FitterMultiHist::operator()(const std::vector<double>& pars) const
{
   double chi2;
   int ndf;
   CalcChi2(pars, chi2, ndf);
   return chi2;
}

TGraphAsymmErrors* FitterMultiHist::getCSRPlot(const char* name, string xTitle, string yTitle){
   TGraphAsymmErrors* tge = new TGraphAsymmErrors();
   string title = (string)(";") + xTitle + ";" + yTitle;
   tge->SetName(name);
   tge->SetTitle(title.c_str());
   tge->SetMarkerSize(1.2);
   tge->SetMarkerColor(4);
   tge->SetLineColor(1);   
   tge->SetLineWidth(2);
   tge->SetMarkerStyle(20);
   for( int ip = 0 ; ip < h1_avg->GetNbinsX() ; ip++ ){
      double x     = h1_avg->GetBinContent(ip+1);
      double y     = m_pars[ip+m_n_common];
      double err   = m_errs[ip+m_n_common];
      double xErrL =   x - h1_avg->GetBinLowEdge(ip+1);
      double xErrH = - x + h1_avg->GetBinLowEdge(ip+1) + h1_avg->GetBinWidth(ip+1);
      tge->SetPoint     (ip, x,   y);
      tge->SetPointError(ip, xErrL, xErrH, err, err);
   }
   return tge;
}
