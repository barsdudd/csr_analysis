#include <stdlib.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <math.h>

#include "Variables.h"

///// ROOT /////
int markerList[] = {20, 21, 22, 23, 29, 33, 34, 39, 41, 43, 45, 47, 48, 49};
int  colorList[] = {1, 2, 3, 6, 8, 9};
///// ROOT /////

//// FILE ////
const char* fileDir = "/seaquest/users/knagai/R008_analysis/csr_analysis/data/";
const char* fileDirNim3 = "/seaquest/users/arunts/NIM3/NIM3_profiles/";
//// FILE ////

///// TARGET /////
const char* targetType[] = {"All", "LH2", "Empty", "LD2", "None", "Fe", "C", "W"};
double targetLength = 50.8;
double thickness[]  = {0.0, 50.8, 50.8, 50.8, 0.0, 1.905, 3.322, 0.953};
double getRD(int rs){
   if( rs < 68 )
      return 1.007708;
   return 1;
}

double getFracD(int rs){
// This is the molecular fraction of D in Deuterium
   if( rs < 68 )
      return 0.918;
   return 1;
}

double getFracHD(int rs){
// This is the molecular fraction of HD in sample
   if( rs < 68 )
      return 0.082;
   return 0;
}

double getRDAverage(vector<int> rs, vector<double> nDimuons, bool potMode){
   //RD = relative volume of D compared to normal Deuterium
   if( rs.size() != nDimuons.size() ){ 
      cout << rs.size() << "\t" << nDimuons.size() << endl;
      exit(0);
   }
   double average = 0;
   double total   = 0;
   for( int irs = 0 ; irs < (int)rs.size() ; irs++ ){
//      double weight = potMode ? getPoT(rs[irs], 3,  true, true) : nDimuons[irs];
      double weight = potMode ? getPoT(rs[irs], 3, false, true) : nDimuons[irs];
      average += weight * getRD(rs[irs]);
      total   += weight;
   }
   return average / total;
}
double getFracDAverage(vector<int> rs, vector<double> nDimuons, bool potMode){
   double average = 0;
   double total   = 0;
   for( int irs = 0 ; irs < (int)rs.size() ; irs++ ){
//      double weight = potMode ? getPoT(rs[irs], 3,  true, true) : nDimuons[irs];
      double weight = potMode ? getPoT(rs[irs], 3, false, true) : nDimuons[irs];
      average += weight * getFracD(rs[irs]);
      total   += weight;
   }
   return average / total;
}


double getFracHDAverage(vector<int> rs, vector<double> nDimuons, bool potMode){
   double average = 0;
   double total   = 0;
   for( int irs = 0 ; irs < (int)rs.size() ; irs++ ){
//      double weight = potMode ? getPoT(rs[irs], 3,  true, true) : nDimuons[irs];
      double weight = potMode ? getPoT(rs[irs], 3, false, true) : nDimuons[irs];
      average += weight * getFracHD(rs[irs]);
      total   += weight;
   }
   return average / total;
}

double getLambdaLD2T(int rs){
   if( rs < 68 )
      return 448.2; // mixture LD2 & HD
   return 424.9; // pure LD2

}

double getLambdaLD2TAverage(vector<int> rs, vector<double> nDimuons, bool potMode){
   double c = getFracHDAverage(rs, nDimuons, potMode);

   return 1. / (
      rho_d * sigma_h * avogadro *      c /   2.   /  getRDAverage(rs, nDimuons, potMode) / aMass_d  + 
      rho_d * sigma_d * avogadro * (1 - c /   2. ) /  getRDAverage(rs, nDimuons, potMode) / aMass_d 
      ) ;

}

double getPurity(int rs){
   if( rs < 68 )
      return 0.958;
   return 1.;
}
double getPurityError(int rs){
   if( rs < 68 )
      return 0.002;
   return 0.;
}
double getPurityAverage(vector<int> rs, vector<double> nDimuons, bool potMode){
   double average = 0;
   double total   = 0;
   for( int irs = 0 ; irs < (int)rs.size() ; irs++ ){
//      double weight = potMode ? getPoT(rs[irs], 3,  true, true) : nDimuons[irs];
      double weight = potMode ? getPoT(rs[irs], 3, false, true) : nDimuons[irs];
      average += weight * getPurity(rs[irs]);
      total   += weight;
   }
   return average / total;
}
///// TARGET /////

//// BEAM ////
int turns   = 369000;
int buckets =    588;
double getPoT(int rs, int targetPos, bool rawMode, bool divMode){
   ostringstream oss;
   ifstream ifs;

   string dummy;
   double rawPoT;
   double livePoT;
   int    it = 0;

//    if( rs == 67 ){
//       if( targetPos == 1 )
//          return 1.221932e+16;
//       if( targetPos == 2 )
//          return 5.45744e+15;
//       if( targetPos == 3 )
//          return 1.3623e+16; 
//    }

   oss.str("");
   oss << "/data2/production/list/R008/good_pot_";
   if( rs == 67 && divMode ){
      oss << "67_1";
   }else if( rs == 68 && divMode ){
      oss << "67_2";
   }else{
      oss << rs;
   }
   oss << ".txt";
   ifs.open(oss.str().c_str());
   while( ifs >> dummy >> dummy >> rawPoT >> livePoT ){
      if( it == targetPos ){
         if( rawMode ) return  rawPoT;
         else          return livePoT;
      }
      it++;
   }
   cout << "NO SUCH DATA" << endl;
   exit(0);
   return 0;
}

double getPoT(vector<int> rs, int targetPos){
   double PoT = 0;
   for( int irs = 0 ; irs < (int)rs.size() ; irs++ ){
      PoT += getPoT(rs[irs], targetPos);
   }
   return PoT;
}

double getRawPoT(vector<int> rs, int targetPos){
   double PoT = 0;
   for( int irs = 0 ; irs < (int)rs.size() ; irs++ ){
      PoT += getPoT(rs[irs], targetPos, true);
   }
   return PoT;
}

double getPedestal(int spillID){
   if( spillID < 450000 )
      return 36.2;
   return 32.6;
}

double getBeamOffset(int rs){
   return ( rs == 57 || rs == 59 ) ? 0.4 : 1.6;
}

double shiftIntensity(double intensity, double by){
   return intensity * by;
}
//// BEAM ////

///// X2 /////
double x2_bin[] = {0.1, 0.13, 0.16, 0.195, 0.24, 0.29, 0.35, 0.45};//, 0.58};
int nBinsX2 = sizeof x2_bin / sizeof x2_bin[0] - 1;
void shiftX2(double by){
   for( int ix = 0 ; ix < nBinsX2+1 ; ix++ )
      x2_bin[ix] += by;
}
int getIX2(double x2){
   // if( x2 < x2_bin[0] || x2 > x2_bin[nBinsX2] )
   //    return -1;
   int ix = -1;
   while( x2 >= x2_bin[ix+1] ) ix++;
   if( ix >= nBinsX2 ) ix = -1;
   return ix;
}
///// X2 /////

///// xF /////
double xF_bin[] = {-0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.95};
int nBinsXF = sizeof xF_bin / sizeof xF_bin[0] - 1;
int getIXF(double xF){
   if( xF <= xF_bin[0] || xF > xF_bin[nBinsXF] )
      return -1;
   int ix = -1;
   while( xF > xF_bin[ix+1] ) ix++;
   return ix;
}
///// xF /////

///// mass /////
double mass_bin[] = {4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.5, 6.5, 8.8};
int nBinsMass = sizeof mass_bin / sizeof mass_bin[0] - 1;
int getIMass(double mass){
   if( mass <= mass_bin[0] || mass > mass_bin[nBinsMass] )
      return -1;
   int ix = -1;
   while( mass > mass_bin[ix+1] ) ix++;
   return ix;
}
///// mass /////

///// x1 /////
double x1_bin[] = {0.3, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.95};
int nBinsX1 = sizeof x1_bin / sizeof x1_bin[0] - 1;
int getIX1(double x1){
   if( x1 <= x1_bin[0] || x1 > x1_bin[nBinsX1] )
      return -1;
   int ix = -1;
   while( x1 > x1_bin[ix+1] ) ix++;
   return ix;
}
///// x1 /////

///// pT /////
//double pT_bin[] = {0, 0.3, 0.5, 0.7, 0.9, 1.2, 1.5, 1.8, 2.3};
//double pT_bin[] = {0, 0.35, 0.54, 0.72, 0.99, 1.42};
double pT_bin[] = {0, 0.41, 0.62, 0.82, 1.1, 2.475};
int nBinsPT = sizeof pT_bin / sizeof pT_bin[0] - 1;
int getIPT(double pT){
   if( pT <= pT_bin[0] || pT > pT_bin[nBinsPT] )
      return -1;
   int ix = -1;
   while( pT > pT_bin[ix+1] ) ix++;
   return ix;
}
///// pT /////

////// CSR //////
double aMass_h   = 1.008;
//double aMass_h  = 1.00794;
double aMass_d   = 2.014;
double sigma_h   = 32.2e-27;
double sigma_d   = 46.6e-27;
//double lambda_h = 52.0;
double lambda_d;
double lambda_h  = 734.5;
double rho_h     = 0.0708;
double rho_d     = 0.1634;
double rho   [8] = {0.0, 0.0708, 0.0 , 0.1634, 0.0, 7.87, 1.80, 19.30};
double lambda[8] = {0.0, 734.46, 0.0, 448.2, 0.0, 16.78, 47.61, 9.94};
double avogadro  = 6.0221409e23;

double aMassNum[8] = {0, 1, 0, 2, 0, 56, 12, 184};
double aMass_a [8] = {0, 1.008, 0, 2.0014, 0, 55.845, 12, 183.84};

//double avogadro = 6.022140857e23;
double effAtoms(double aMass, double rho, double lambda){
   double alpha = 1 - exp(-targetLength * rho / lambda);
   return avogadro * alpha * lambda / aMass;
}
double effLinear(double d2, double h2, int rs){
   double purity = getPurity(rs);
   return d2 * purity + h2 * ( 1 - purity );
}
double effInverse(double d2, double h2, int rs){
   double purity = getPurity(rs);
   return 1 / ( purity / d2 + ( 1 - purity ) / h2 );
}
double effLinear(double d2, double h2, double purity){
   return d2 * purity + h2 * ( 1 - purity );
}
double effInverse(double d2, double h2, double purity){
   return 1 / ( purity / d2 + ( 1 - purity ) / h2 );
}
////// CSR //////

////// kEff //////
double p0     [2] = {0.9886, 0.7858};
double slopeP1[2] = {-0.001632, -0.0009815};
double slopeP0[2] = {-0.001852, -0.0018180};
double getKEff(double x2, int occD1, int targetPos){
   int index = targetPos == 2 ? 1 : 0;
   return p0[index] + ( slopeP0[index] + slopeP1[index] * x2 ) * occD1;
}
////// kEff //////

////// SPECIAL RUNS //////
int NIM3SpecialRuns[] = {
   14124,
   14303,
//   14888,
   15509
};
////// SPECIAL RUNS //////
