#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <TFitterMinuit.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include "FitterIndividual.h"
using namespace std;

FitterIndividual::FitterIndividual()
{
  m_label  = "indi";
  m_n_par  = 9999; // not used
}

void FitterIndividual::DrawFit()
{
  ofstream ofs("result_fit.txt");
  double chi2_tot = 0;
  int     ndf_tot = 0;

  ostringstream oss;
  TCanvas* c1 = new TCanvas("c1", "");
  c1->SetGrid();
  for (int ih = 0; ih < N_HIST; ih++) {
    TGraphErrors* gr = m_list_gr[ih];
    gr->SetLineColor  (kRed);
    gr->SetMarkerColor(kRed);
    gr->SetMarkerStyle(8);
    gr->GetYaxis()->SetRangeUser(0, 3);
    gr->Fit("pol2", "EM", "E1");

    TF1* f1     = gr->GetFunction("pol2");
    double chi2 = f1->GetChisquare();
    int    ndf  = f1->GetNDF();
    double p0   = f1->GetParameter(0);
    double p1   = f1->GetParameter(1);
    double p2   = f1->GetParameter(1);
    double e0   = f1->GetParError (0);
    double e1   = f1->GetParError (1);
    double e2   = f1->GetParError (1);
    f1->SetLineColor(kBlue);

    ofs << chi2 << " " << ndf << " "
        << p0 << " \\pm " << e0 << " " << p1 << " \\pm " << e1 << " " << p2 << " \\pm " << e2 << "\n";
    chi2_tot += chi2;
    ndf_tot  += ndf;

    TLatex tex;
    tex.SetNDC(true);
    oss.str("");
    oss << "p0 = " << p0 << " #pm " << e0;
    tex.DrawLatex(0.15, 0.85, oss.str().c_str());
    oss.str("");
    oss << "p1 = " << p1 << " #pm " << e1;
    tex.DrawLatex(0.15, 0.80, oss.str().c_str());
    oss.str("");
    oss << "p2 = " << p2 << " #pm " << e2;
    tex.DrawLatex(0.15, 0.75, oss.str().c_str());

    oss.str("");
    oss << "hist_fitted_" << ih << ".png";
    c1->SaveAs(oss.str().c_str());
  }
  delete c1;
  ofs << chi2_tot << "/" << ndf_tot << endl;
  ofs.close();
}

void FitterIndividual::CalcChi2(const std::vector<double>& pars, double& chi2, int& ndf) const
{
  cerr << "ERROR:  This CalcChi2() function must not be called.  Abort." << endl;
  exit(1);
}
