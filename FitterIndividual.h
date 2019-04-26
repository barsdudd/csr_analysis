#ifndef __FITTER_INDIVIDUAL_H__
#define __FITTER_INDIVIDUAL_H__
#include "FitterMultiHist.h"
class TH1;

/// Class for fitting multiple histograms individually
/** Actually this class doesn't use most part of FitterMultiHist...
 */
class FitterIndividual : public FitterMultiHist {
 public:
  FitterIndividual();
  virtual ~FitterIndividual() {;}
  
  virtual void DrawFit();

 protected:
  virtual void InitPar() {;}
  virtual void CalcChi2(const std::vector<double>& pars, double& chi2, int& ndf) const;
};

#endif // __FITTER_INDIVIDUAL_H__
