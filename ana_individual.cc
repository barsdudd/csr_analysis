#include "FitterIndividual.h"
using namespace std;

int main(int argc, char *argv[])
{
  FitterMultiHist* fit = new FitterIndividual();

  /// Choose one of the data files (cf. data/00note.txt)
  //fit->Init("data/results.root");
//  fit->Init("data/plotRF00_2111v42_arun_v2.root");
  fit->ChangeDir();
  fit->DrawFit();
}
