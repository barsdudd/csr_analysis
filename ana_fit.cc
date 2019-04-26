#include "Fitter2DPol2.h"
#include "Fitter2DCos.h"
using namespace std;

/** Main function to run the fitting of multiple "ratio vs intensity" histograms simultaneously.
 */
int main(int argc, char *argv[])
{
  FitterMultiHist* fit = new Fitter2DPol2();
  //FitterMultiHist* fit = new Fitter2DCos();

  /// Choose one of the data files (cf. data/00note.txt)
  //fit->Init("data/results.root");
//  fit->Init("data/plotRF00_2111v42_arun_v2.root");

  fit->ChangeDir();
  fit->DoFit();
  fit->PrintFit("result_fit.txt");
  fit->PrintFit(cout);
  fit->DrawFit();
}
