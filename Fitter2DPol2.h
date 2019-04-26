#ifndef __FITTER_2D_POL2_H__
#define __FITTER_2D_POL2_H__
#include "FitterMultiHist.h"
class TH1;

/// Class for fitting multiple histograms simultaneously
class Fitter2DPol2 : public FitterMultiHist {
 public:
  Fitter2DPol2();
  virtual ~Fitter2DPol2() {;}
  
  virtual void DrawFit();

 protected:
  virtual void InitPar();
  virtual void CalcChi2(const std::vector<double>& pars, double& chi2, int& ndf) const;
};

#endif // __FITTER_2D_POL2_H__
