#ifndef __FIT_MULTI_HIST_H__
#define __FIT_MULTI_HIST_H__
#include <string>
#include <vector>
#include <TFitterMinuit.h>
#include <Minuit2/FCNBase.h>
using namespace std;
class TGraphAsymmErrors;
class TGraphErrors;
class TH1D;
/// Class for fitting multiple histograms simultaneously
class FitterMultiHist : public ROOT::Minuit2::FCNBase {
protected:
   int N_HIST;
   std::string hist_name;
   /* static const double X2_LOW[N_HIST+1]; */
   /* static const double X2_AVG[N_HIST]; */
   /* double VAL_AVG[20]; */
   TGraphErrors* m_list_gr[20];
   TFitterMinuit* m_minuit;
   TH1D* h1_avg;
   
   std::string m_label;
   int    m_n_par;
   int    m_ndf;
   int    m_n_common;
   double m_chi2;
   
   std::vector<double> m_pars;
   std::vector<double> m_errs;
   double*             m_cov_mat;
   
   FitterMultiHist(); ///< Not "public" since this class is abstract.
   virtual ~FitterMultiHist() {;}
   
public:
   void Init(const char* fn_hist, const char* gn, const char* avg_hist);
   void ChangeDir();
   void DoFit();
   void PrintFit(std::string fname);
   virtual void PrintFit(std::ostream& os);
   virtual void DrawFit();
   
   std::string Label() { return m_label; }
   int    NPar() { return m_n_par; }
   int    NDF () { return m_ndf ; }
   double Chi2() { return m_chi2; }
   double Par(const int i) { return m_pars[i]; }
   double Err(const int i) { return m_errs[i]; }
   double* CovMat() { return m_cov_mat; }

   TGraphAsymmErrors* getCSRPlot(const char* name, string xTitle="", string yTitle="Cross Section Ratio");
  
protected:
   virtual void InitPar();
   virtual void CalcChi2(const std::vector<double>& pars, double& chi2, int& ndf) const;

   double operator()(const std::vector<double>& r) const;
   double Up() const { return 1.0; }
};

#endif // __FIT_MULTI_HIST_H__
