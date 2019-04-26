#ifndef __VARIABLES_H__
#define __VARIABLES_H__

#include <stdlib.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>

#include <map>

using namespace std;

///// ROOT /////
extern int markerList[];
extern int  colorList[];
///// ROOT /////

//// FILE ////
extern const char* fileDir;
extern const char* fileDirNim3;
//// FILE ////

///// TARGET /////
extern const char* targetType[];
extern double targetLength;
extern double thickness[];
double getRD(int rs);
double getFracD(int rs);
double getFracHD(int rs);
double getRDAverage(vector<int> rs, vector<double> nDimuons, bool potMode=false);
double getFracDAverage(vector<int> rs, vector<double> nDimuons, bool potMode=false);
double getFracHDAverage(vector<int> rs, vector<double> nDimuons, bool potMode=false);
double getLambdaLD2T(int rs);
double getLambdaLD2TAverage(vector<int> rs, vector<double> nDimuons, bool potMode=false);
double getPurity(int rs);
double getPurityError(int rs);
double getPurityAverage(vector<int> rs, vector<double> nDimuons, bool potMode=false);
///// TARGET /////

//// BEAM /////
extern int turns;
extern int buckets;
double getPoT(int rs, int targetPos, bool rawMode=false, bool divMode=false);
double getPoT(vector<int> rs, int targetPos);
double getRawPoT(vector<int> rs, int targetPos);
double getPedestal(int spillID);
double getBeamOffset(int rs);
double shiftIntensity(double intensity, double by);
//// BEAM /////

/////  X2  /////
extern double x2_bin[];
extern int nBinsX2;
void shiftX2(double by);
int getIX2(double x2);
/////  X2  /////

/////  xF  /////
extern double xF_bin[];
extern int nBinsXF;
//void shiftX2(double by);
int getIXF(double xF);
/////  xF  /////

/////  mass  /////
extern double mass_bin[];
extern int nBinsMass;
int getIMass(double mass);
/////  mass  /////

/////  x1  /////
extern double x1_bin[];
extern int nBinsX1;
int getIX1(double x1);
/////  x1  /////

/////  pT  /////
extern double pT_bin[];
extern int nBinsPT;
int getIPT(double pT);
/////  pT  /////

////// CSR //////
extern double aMass_h;
extern double aMass_d;
extern double sigma_d;
extern double sigma_h;
extern double lambda_h;
extern double lambda_d;
extern double rho_h;
extern double rho_d;
extern double rho[];
extern double lambda[];
extern double avogadro;
extern double aMassNum[8];
extern double aMass_a[8];
double effAtoms  (double aMass, double rho, double lambda);
double effLinear (double d2, double h2, int rs);
double effInverse(double d2, double h2, int rs);
double effLinear (double d2, double h2, double purity);
double effInverse(double d2, double h2, double purity);
////// CSR //////

////// kEff //////
extern double p0     [2];
extern double slopeP1[2];
extern double slopeP0[2];
double getKEff(double x2, int occD1, int targetPos);
////// kEff //////

////// SPECIAL RUNS //////
extern int NIM3SpecialRuns[];
////// SPECIAL RUNS //////

#endif // __VARIABLES_H__
