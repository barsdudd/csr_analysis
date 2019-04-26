#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <TFitterMinuit.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include "FitterCommonPol2.h"
using namespace std;

FitterCommonPol2::FitterCommonPol2()
{
   m_label    = "common_pol2";
   m_n_common = 2;
   m_n_par    = N_HIST + m_n_common;
}

void FitterCommonPol2::DrawFit()
{
   ostringstream oss;
   TCanvas* c1 = new TCanvas("c1", "", 800, 600);
   c1->SetGrid();

   double avg;
   for (int ih = 0; ih < N_HIST; ih++) {
      TGraphErrors* gr = m_list_gr[ih];
      gr->SetLineColor  (kRed);
      gr->SetMarkerColor(kRed);
      gr->SetMarkerStyle(8);
      gr->GetYaxis()->SetRangeUser(0, 3);
      gr->Draw("APE1");
      avg = h1_avg->GetBinContent(ih+1);
      TF1 func("func", "[0] + [1] * x  + [2] * x*x", 0, 1e8);
      double p2 = m_pars[0];
      double p1 = m_pars[1];
      double p0 = m_pars[ih + m_n_common];
      double e2 = m_errs[0];
      double e1 = m_errs[1];
      double e0  = m_errs[ih + m_n_common];
      //func.SetParameters(p0, p10, p11, p20, p21, ih);
//    func.SetParameters(p0, p10, p11, p20, p21, X2_AVG[ih]);
      func.SetParameters(p0, p1, p2, avg );
      func.SetLineColor(kBlue);
      func.Draw("Csame");

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
      oss << hist_name << "_common_pol2_fitted_" << ih << ".png";
      c1->SaveAs(oss.str().c_str());
   }
   delete c1;
}

void FitterCommonPol2::InitPar()
{
   ostringstream oss;
   m_minuit->SetParameter(0, "p2 (common)" , 0.0, 1e-6, -1.0, +1.0);
   m_minuit->SetParameter(1, "p1 (common)" , 0.0, 1e-6, -1.0, +1.0);
   for (int ih = 0; ih < N_HIST; ih++) {
      oss.str("");
      oss << "p0 for hist #" << ih;
      m_minuit->SetParameter(ih + m_n_common, oss.str().c_str(), 1.0, 0.01, 0.0, 2.0);
   }
}

void FitterCommonPol2::CalcChi2(const std::vector<double>& pars, double& chi2, int& ndf) const
{
   chi2 = 0;
   ndf  = 0;
//   cout << "N_HIST = " << N_HIST << "\n";
   for (int ih = 0; ih < N_HIST; ih++) {
      double p2 = pars[0]; // common slope
      double p1 = pars[1]; // common slope
      double p0 = pars[ih+m_n_common];
      TGraphErrors* gr = m_list_gr[ih];
      int n_pt = gr->GetN();
      for (int i_pt = 0; i_pt < n_pt; i_pt++) {
         double err = gr->GetErrorY(i_pt);
         if (err == 0) continue;

         double cent, cont;
         gr->GetPoint(i_pt, cent, cont);

         //double pred = p0 + (p10 + p11 * ih) * cent  + (p20 + p21 * ih) * pow(cent, 2);
         double pred = p0 + p1 * cent  + p2 *  pow(cent, 2);
         chi2 += pow((cont - pred) / err, 2);
         ndf++;
      }
   }
   ndf -= m_minuit->GetNumberFreeParameters();
}
