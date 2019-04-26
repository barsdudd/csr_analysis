#ifndef __FITTER_2D_COS_H__
#define __FITTER_2D_COS_H__
#include "FitterMultiHist.h"
class TH1;

/// Class for fitting multiple histograms simultaneously
class Fitter2DCos : public FitterMultiHist {
 public:
  Fitter2DCos();
  virtual ~Fitter2DCos() {;}
  
  virtual void DrawFit();

 protected:
  virtual void InitPar();
  virtual void CalcChi2(const std::vector<double>& pars, double& chi2, int& ndf) const;
};

#endif // __FITTER_2D_COS_H__
