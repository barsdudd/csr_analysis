#ifndef __FITTER_COMMON_POL2_H__
#define __FITTER_COMMON_POL2_H__
#include "FitterMultiHist.h"
class TH1;

/// Class for fitting multiple histograms simultaneously
class FitterCommonPol2 : public FitterMultiHist {
 public:
  FitterCommonPol2();
  virtual ~FitterCommonPol2() {;}
  
  virtual void DrawFit();

 protected:
  virtual void InitPar();
  virtual void CalcChi2(const std::vector<double>& pars, double& chi2, int& ndf) const;
};

#endif // __FITTER_COMMON_POL2_H__
