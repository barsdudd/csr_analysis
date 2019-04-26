#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>

#include <map>

#include <TROOT.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include "Event.h"
#include "Variables.h"

using namespace std;

TFile* file;
TTree* tree;
//defining the strings for making querries in sql
//declaring all the variables in kDimuon field in merged roadset 57
int dimuonID, runID, spillID, eventID, posTrackID, negTrackID, targetPos;

double dx, dy, dz, dpx, dpy, dpz, mass, xF, xB, xT, costh, phi, trackSeparation, chisq_dimuon, px1_kDimuon, py1_kDimuon, pz1_kDimuon, px2_kDimuon, py2_kDimuon, pz2_kDimuon, isValid, isTarget, isDump; 
    
//declaring all the variables in kTrack field for mu_p
int trackID_mu_p, runID_mu_p, eventID_mu_p, spillID_mu_p, roadID_mu_p, charge_mu_p, numHits_mu_p, numHitsSt1_mu_p, numHitsSt2_mu_p, numHitsSt3_mu_p, numHitsSt4H_mu_p, numHitsSt4V_mu_p; 
double chisq_mu_p, x0_mu_p, y0_mu_p, z0_mu_p, xD_mu_p, yD_mu_p, xT_mu_p, yT_mu_p, pxD_mu_p, pyD_mu_p, pzD_mu_p, pxT_mu_p, pyT_mu_p, pzT_mu_p, z0x_mu_p, z0y_mu_p, px0_mu_p, py0_mu_p, pz0_mu_p, x1_mu_p, y1_mu_p, z1_mu_p, px1_kTrack_mu_p, py1_kTrack_mu_p, pz1_kTrack_mu_p, x3_mu_p, y3_mu_p, z3_mu_p, px3_kTrack_mu_p, py3_kTrack_mu_p, pz3_kTrack_mu_p, thbend_mu_p, tx_PT_mu_p, ty_PT_mu_p, chisq_target_mu_p, chisq_dump_mu_p, chisq_upstream_mu_p;

//declaring all the variables in kTrack field for mu_m
int trackID_mu_m,runID_mu_m, eventID_mu_m, spillID_mu_m, roadID_mu_m, charge_mu_m, numHits_mu_m, numHitsSt1_mu_m, numHitsSt2_mu_m, numHitsSt3_mu_m, numHitsSt4H_mu_m, numHitsSt4V_mu_m; 

double chisq_mu_m, x0_mu_m, y0_mu_m, z0_mu_m, xD_mu_m, yD_mu_m, xT_mu_m, yT_mu_m, pxD_mu_m, pyD_mu_m, pzD_mu_m, pxT_mu_m, pyT_mu_m, pzT_mu_m, z0x_mu_m, z0y_mu_m, px0_mu_m, py0_mu_m, pz0_mu_m, x1_mu_m, y1_mu_m, z1_mu_m, px1_kTrack_mu_m, py1_kTrack_mu_m, pz1_kTrack_mu_m, x3_mu_m, y3_mu_m, z3_mu_m, px3_kTrack_mu_m, py3_kTrack_mu_m, pz3_kTrack_mu_m, thbend_mu_m, tx_PT_mu_m, ty_PT_mu_m, chisq_target_mu_m, chisq_dump_mu_m, chisq_upstream_mu_m;

//declaring all the occupancy variables
int D1, D2, D3, H1, H2, H3, H4, P1, P2, D1L, D1R, D2L, D2R, D3L, D3R;
  
//declaring all the variables in QIE 
double  Intensity_p;

//declaring all the variables in BeamDAQ
double  QIEsum, liveProton, RF00, G2SEM, inh_thres;


void prepareOutput(int rs, int index){
   ostringstream oss;
   oss.str("");
   oss << "/seaquest/users/knagai/R008_analysis/data/raw_data/" << rs;
   gSystem->mkdir(oss.str().c_str(), true);
   oss << "/R008_without_cuts_" << rs << "_" << index << ".root";
   file = new TFile(oss.str().c_str(), "recreate");
   oss.str("");
   oss << "roadset_" << rs;
   tree = new TTree(oss.str().c_str(), "save");
   
   tree->Branch("Intensity_p",&Intensity_p);
   tree->Branch("trackID_mu_p",&trackID_mu_p);
   tree->Branch("charge_mu_p",&charge_mu_p);
   tree->Branch("runID_mu_p",&runID_mu_p);
   tree->Branch("spillID_mu_p",&spillID_mu_p);
   tree->Branch("eventID_mu_p", &eventID_mu_p);
   tree->Branch("roadID_mu_p", &roadID_mu_p);
   tree->Branch("numHits_mu_p", &numHits_mu_p);
   tree->Branch("numHitsSt1_mu_p",&numHitsSt1_mu_p);
   tree->Branch("numHitsSt2_mu_p",&numHitsSt2_mu_p);
   tree->Branch("numHitsSt3_mu_p", &numHitsSt3_mu_p);
   tree->Branch("numHitsSt4H_mu_p", &numHitsSt4H_mu_p);
   tree->Branch("numHitsSt4V_mu_p", &numHitsSt4V_mu_p);
   tree->Branch("chisq_mu_p", &chisq_mu_p);
   tree->Branch("x0_mu_p", &x0_mu_p);
   tree->Branch("y0_mu_p", &y0_mu_p);
   tree->Branch("z0_mu_p", &z0_mu_p);
   tree->Branch("xD_mu_p", &xD_mu_p);
   tree->Branch("yD_mu_p", &yD_mu_p);
   tree->Branch("xT_mu_p", &xT_mu_p);
   tree->Branch("yT_mu_p", &yT_mu_p);
   tree->Branch("pxD_mu_p", &pxD_mu_p);
   tree->Branch("pyD_mu_p", &pyD_mu_p);
   tree->Branch("pzD_mu_p", &pzD_mu_p);
   tree->Branch("pxT_mu_p", &pxT_mu_p);
   tree->Branch("pyT_mu_p", &pyT_mu_p);
   tree->Branch("pzT_mu_p", &pzT_mu_p);
   tree->Branch("z0x_mu_p", &z0x_mu_p);
   tree->Branch("z0y_mu_p", &z0y_mu_p);
   tree->Branch("px0_mu_p", &px0_mu_p);
   tree->Branch("py0_mu_p", &py0_mu_p);
   tree->Branch("pz0_mu_p", &pz0_mu_p);
   tree->Branch("x1_mu_p", &x1_mu_p);
   tree->Branch("y1_mu_p", &y1_mu_p);
   tree->Branch("z1_mu_p", &z1_mu_p);
   tree->Branch("px1_kTrack_mu_p", &px1_kTrack_mu_p);
   tree->Branch("py1_kTrack_mu_p", &py1_kTrack_mu_p);
   tree->Branch("pz1_kTrack_mu_p", &pz1_kTrack_mu_p);
   tree->Branch("x3_mu_p", &x3_mu_p);
   tree->Branch("y3_mu_p", &y3_mu_p);
   tree->Branch("z3_mu_p", &z3_mu_p);
   tree->Branch("px3_kTrack_mu_p", &px3_kTrack_mu_p);
   tree->Branch("py3_kTrack_mu_p", &py3_kTrack_mu_p);
   tree->Branch("pz3_kTrack_mu_p", &pz3_kTrack_mu_p);
   tree->Branch("thbend_mu_p", &thbend_mu_p);
   tree->Branch("tx_PT_mu_p", &tx_PT_mu_p);
   tree->Branch("ty_PT_mu_p", &ty_PT_mu_p);
   tree->Branch("chisq_target_mu_p", &chisq_target_mu_p);
   tree->Branch("chisq_dump_mu_p", &chisq_dump_mu_p);
   tree->Branch("chisq_upstream_mu_p", &chisq_upstream_mu_p);
   
   
   tree->Branch("trackID_mu_m",&trackID_mu_m);
   tree->Branch("runID_mu_m",&runID_mu_m);
   tree->Branch("spillID_mu_m",&spillID_mu_m);
   tree->Branch("eventID_mu_m", &eventID_mu_m);
   tree->Branch("charge_mu_m",&charge_mu_m);
   tree->Branch("roadID_mu_m",&roadID_mu_m);
   tree->Branch("numHits_mu_m", &numHits_mu_m);
   tree->Branch("numHitsSt1_mu_m", &numHitsSt1_mu_m);
   tree->Branch("numHitsSt2_mu_m", &numHitsSt2_mu_m);
   tree->Branch("numHitsSt3_mu_m", &numHitsSt3_mu_m);
   tree->Branch("numHitsSt4H_mu_m", &numHitsSt4H_mu_m);
   tree->Branch("numHitsSt4V_mu_m", &numHitsSt4V_mu_m);
   tree->Branch("chisq_mu_m", &chisq_mu_m);
   tree->Branch("x0_mu_m", &x0_mu_m);
   tree->Branch("y0_mu_m", &y0_mu_m);
   tree->Branch("z0_mu_m", &z0_mu_m);
   tree->Branch("xD_mu_m", &xD_mu_m);
   tree->Branch("yD_mu_m", &yD_mu_m);
   tree->Branch("xT_mu_m", &xT_mu_m);
   tree->Branch("yT_mu_m", &yT_mu_m);
   tree->Branch("pxD_mu_m", &pxD_mu_m);
   tree->Branch("pyD_mu_m", &pyD_mu_m);
   tree->Branch("pzD_mu_m", &pzD_mu_m);
   tree->Branch("pxT_mu_m", &pxT_mu_m);
   tree->Branch("pyT_mu_m", &pyT_mu_m);
   tree->Branch("pzT_mu_m", &pzT_mu_m);
   tree->Branch("z0x_mu_m", &z0x_mu_m);
   tree->Branch("z0y_mu_m", &z0y_mu_m);
   tree->Branch("px0_mu_m", &px0_mu_m);
   tree->Branch("py0_mu_m", &py0_mu_m);
   tree->Branch("pz0_mu_m", &pz0_mu_m);
   tree->Branch("x1_mu_m", &x1_mu_m);
   tree->Branch("y1_mu_m", &y1_mu_m);
   tree->Branch("z1_mu_m", &z1_mu_m);
   tree->Branch("px1_kTrack_mu_m", &px1_kTrack_mu_m);
   tree->Branch("py1_kTrack_mu_m", &py1_kTrack_mu_m);
   tree->Branch("pz1_kTrack_mu_m", &pz1_kTrack_mu_m);
   tree->Branch("x3_mu_m", &x3_mu_m);
   tree->Branch("y3_mu_m", &y3_mu_m);
   tree->Branch("z3_mu_m", &z3_mu_m);
   tree->Branch("px3_kTrack_mu_m", &px3_kTrack_mu_m);
   tree->Branch("py3_kTrack_mu_m", &py3_kTrack_mu_m);
   tree->Branch("pz3_kTrack_mu_m", &pz3_kTrack_mu_m);
   tree->Branch("thbend_mu_m", &thbend_mu_m);
   tree->Branch("tx_PT_mu_m", &tx_PT_mu_m);
   tree->Branch("ty_PT_mu_m", &ty_PT_mu_m);
   tree->Branch("chisq_target_mu_m", &chisq_target_mu_m);
   tree->Branch("chisq_dump_mu_m", &chisq_dump_mu_m);
   tree->Branch("chisq_upstream_mu_m", &chisq_upstream_mu_m);
   
   
   tree->Branch("dimuonID", &dimuonID);
   tree->Branch("targetPos", &targetPos);
   tree->Branch("posTrackID", &posTrackID);
   tree->Branch("negTrackID", &negTrackID);
   tree->Branch("dx", &dx);
   tree->Branch("dy", &dy);
   tree->Branch("dz", &dz);
   tree->Branch("dpx", &dpx);
   tree->Branch("dpy", &dpy);
   tree->Branch("dpz", &dpz);
   tree->Branch("mass", &mass);
   tree->Branch("xF", &xF);
   tree->Branch("xB", &xB);
   tree->Branch("xT", &xT);
   tree->Branch("costh", &costh);
   tree->Branch("phi", &phi);
   tree->Branch("trackSeparation", &trackSeparation);
   tree->Branch("chisq_dimuon", &chisq_dimuon);
   tree->Branch("px1_kDimuon", &px1_kDimuon);
   tree->Branch("py1_kDimuon", &py1_kDimuon);
   tree->Branch("pz1_kDimuon", &pz1_kDimuon);
   tree->Branch("px2_kDimuon", &px2_kDimuon);
   tree->Branch("py2_kDimuon", &py2_kDimuon);
   tree->Branch("pz2_kDimuon", &pz2_kDimuon);
   tree->Branch("isValid", &isValid);
   tree->Branch("isTarget", &isTarget);
   tree->Branch("isDump", &isDump);
   
   tree->Branch("D1", &D1);
   tree->Branch("D2", &D2);
   tree->Branch("D3", &D3);
   tree->Branch("H1", &H1);
   tree->Branch("H2", &H2);
   tree->Branch("H3", &H3);
   tree->Branch("H4", &H4);
   tree->Branch("P1", &P1);
   tree->Branch("P2", &P2);
   tree->Branch("D1L", &D1L);
   tree->Branch("D1R", &D1R);
   tree->Branch("D2L", &D2L);
   tree->Branch("D2R", &D2R);
   tree->Branch("D3L", &D3L);
   tree->Branch("D3R", &D3R);
   tree->Branch("QIEsum", &QIEsum);
   tree->Branch("G2SEM", &G2SEM);
   //  tree->Branch("inhibit_block_sum",&inhibit_block_sum);
   //  tree->Branch("trigger_sum_no_inhibit",&trigger_sum_no_inhibit);
   tree->Branch("liveProton",&liveProton);  
   tree->Branch("RF00",&RF00);
   tree->Branch("inh_thres",&inh_thres);
   
}

TString getQuery(TString run, TString run_R007){
   TString q1;
   
   q1  = "SELECT ";
   q1 += run_R007 + ".QIE.Intensity_p, ";
   q1 += run      + ".MU_PLUS.*, ";
   q1 += run      + ".MU_MINUS.*, ";
   q1 += run      + ".kDimuon.*, ";
   q1 += run_R007 + ".Occupancy.*, ";
   q1 += run_R007 + ".BeamDAQ.QIEsum, ";
   q1 += run_R007 + ".BeamDAQ.inh_thres, ";
   q1 += run_R007 + ".Spill.liveProton, ";
   q1 += run_R007 + ".QIE.`RF+00`, "; 
   q1 += run_R007 + ".Beam.value as G2SEM, "; 
   q1 += run_R007 + ".Spill.targetPos "; 
   q1 += "FROM " + run_R007 + ".QIE ";
   q1 += "JOIN " + run      + ".kDimuon ON kDimuon.runID = QIE.runID AND kDimuon.eventID = QIE.eventID ";
   q1 += "JOIN " + run      + ".kTrack AS MU_PLUS ON MU_PLUS.trackID = kDimuon.posTrackID AND kDimuon.runID = MU_PLUS.runID ";
   q1 += "JOIN " + run      + ".kTrack AS MU_MINUS ON MU_MINUS.trackID = kDimuon.negTrackID AND kDimuon.runID = MU_MINUS.runID ";
   q1 += "JOIN " + run_R007 + ".Spill ON kDimuon.targetPos = Spill.targetPos AND kDimuon.spillID = Spill.spillID ";
   q1 += "JOIN " + run_R007 + ".Event ON kDimuon.runID = Event.runID AND kDimuon.eventID = Event.eventID ";
   q1 += "JOIN " + run      + ".kEvent ON kDimuon.runID = kEvent.runID AND kDimuon.eventID = kEvent.eventID ";
   q1 += "JOIN " + run_R007 + ".BeamDAQ ON BeamDAQ.spillID = QIE.spillID ";
   q1 += "JOIN " + run_R007 + ".Occupancy ON Occupancy.runID = kDimuon.runID AND Occupancy.spillID = kDimuon.spillID AND kDimuon.eventID = Occupancy.eventID ";
   q1 += "JOIN " + run_R007 + ".Beam ON Beam.name = 'S:G2SEM' AND Beam.spillID = Spill.spillID ";
   q1 += "WHERE Spill.dataQuality = 0 ";
   //  This takes care of all the spill level cuts as given in the redmine site https://cdcvs.fnal.gov/redmine/projects/seaquest-production/wiki/Data_Quality_in_our_Productions
   //  EVENT LEVEL CUTS
   //q1 += "AND MATRIX1 = 1 ";
   //  We are asking only FPGA 1 events
   q1 += "AND kEvent.status = 0. ";
   //  TARGET SELECTION
   //  1 = LIQUID HYDROGEN, 2 = EMPTY, 3 = LIQUID DEUTERIUM, 4 = NONE, 5 = Fe, 6 = C, 7 = W
   q1 += "AND kDimuon.targetPos BETWEEN 1. AND 7. ";
   // or kDimuon.targetPos = 2 or kDimuon.targetPos = 1 or kDimuon.targetPos = 4 ";
   //  TRACK SELECTION CUTS
   //  EVENT LEVEL CUTS                                                                                    
   q1 += "AND Event.MATRIX1 = 1 ";

   q1 += "AND kDimuon.dz BETWEEN -350 AND 50 ";
   // q1 += "AND kDimuon.dpz BETWEEN 37. AND 116. ";
   // q1 += "AND kDimuon.chisq_dimuon < 18. ";
   // q1 += "AND ABS(kDimuon.costh) < 0.5 ";
   // q1 += "AND kDimuon.xF BETWEEN -0.15 AND 0.95 ";

   q1 += "AND ABS(kDimuon.dx) < 1. ";
   q1 += "AND ABS(kDimuon.dy - 1.6) < 1. ";

   // q1 += "AND ABS(kDimuon.dpx) < 2. ";
   // q1 += "AND ABS(kDimuon.dpy) < 2. ";

   // q1 += "AND kDimuon.dx*kDimuon.dx + (kDimuon.dy - 1.6)*(kDimuon.dy - 1.6) < 0.09 ";//RS 62, 67, 70 setting
   // q1 += "AND kDimuon.dpx*kDimuon.dpx + kDimuon.dpy*kDimuon.dpy < 5.5 ";
   // q1 += "AND MU_PLUS.numHits > 13 ";
   // q1 += "AND MU_MINUS.numHits > 13 ";
   // q1 += "AND MU_PLUS.numHits + MU_MINUS.numHits > 29 ";
   // q1 += "AND MU_PLUS.numHitsSt1 + MU_MINUS.numHitsSt1 > 8. ";
   // q1 += "AND MU_PLUS.roadID < 55000. ";// ONLY RS 67 cut
   // q1 += "AND MU_PLUS.chisq/(MU_PLUS.numHits - 5.) < 13. ";
   // q1 += "AND MU_MINUS.chisq/(MU_MINUS.numHits - 5.) < 13. ";
   // q1 += "AND MU_PLUS.pz1 BETWEEN 9. AND 75. ";
   // q1 += "AND MU_MINUS.pz1 BETWEEN 9. AND 75. ";
   // q1 += "AND MU_PLUS.z0 BETWEEN -350. AND 0. ";
   // q1 += "AND MU_MINUS.z0 BETWEEN -350. AND 0. ";
   // q1 += "AND MU_PLUS.xD*MU_PLUS.xD + (MU_PLUS.yD - 1.6)*(MU_PLUS.yD - 1.6)> 12. ";
   // q1 += "AND MU_MINUS.xD*MU_MINUS.xD +(MU_MINUS.yD - 1.6)*(MU_MINUS.yD - 1.6) > 12. ";
   // q1 += "AND MU_PLUS.xD*MU_PLUS.xD + (MU_PLUS.yD - 1.6)*(MU_PLUS.yD - 1.6)< 1200. ";
   // q1 += "AND MU_MINUS.xD*MU_MINUS.xD +(MU_MINUS.yD - 1.6)*(MU_MINUS.yD - 1.6) < 1200. ";
   // q1 += "AND MU_PLUS.xT*MU_PLUS.xT + (MU_PLUS.yT - 1.6)*(MU_PLUS.yT - 1.6) <400. ";
   // q1 += "AND MU_MINUS.xT*MU_MINUS.xT +(MU_MINUS.yT - 1.6)*(MU_MINUS.yT - 1.6) <400. ";


   // q1 += "AND MU_PLUS.chisq_target<1.5*MU_PLUS.chisq_upstream ";
   // q1 += "AND MU_MINUS.chisq_target<1.5*MU_MINUS.chisq_upstream ";
   // q1 += "AND MU_PLUS.chisq_target<1.5*MU_PLUS.chisq_dump ";
   // q1 += "AND MU_MINUS.chisq_target<1.5*MU_MINUS.chisq_dump ";
   // q1 += "AND ABS(MU_PLUS.chisq_target + MU_MINUS.chisq_target - kDimuon.chisq_dimuon) < 2. ";
   // q1 += "AND MU_PLUS.chisq_target<15. ";
   // q1 += "AND MU_MINUS.chisq_target<15. ";
   // q1 += "AND ABS(kDimuon.trackSeparation) < 280. ";


   //  For RS 67 AND 70
   // q1 += "AND ABS(MU_PLUS.px1 - MU_PLUS.px3 + 0.416) < 0.008 ";
   // q1 += "AND ABS(MU_MINUS.px1 - MU_MINUS.px3 - 0.416) < 0.008 ";

   // q1 += "AND ABS(MU_PLUS.py1 - MU_PLUS.py3) < 0.008 ";
   // q1 += "AND ABS(MU_MINUS.py1 - MU_MINUS.py3) < 0.008 ";
   // q1 += "AND ABS(MU_PLUS.pz1 - MU_PLUS.pz3) < 0.08 ";
   // q1 += "AND ABS(MU_MINUS.pz1 - MU_MINUS.pz3) < 0.08 ";
   // //q1 += "AND LEAST(MU_PLUS.pz1,MU_MINUS.pz1) - GREATEST(-MU_PLUS.x1, MU_MINUS.x1)/4.5 > 13 ";//RS 57, 59 AND 62
   // //q1 += "AND LEAST(MU_PLUS.pz1,MU_MINUS.pz1) - GREATEST(MU_PLUS.x1, -MU_MINUS.x1)/4.5 > 13 "; //RS 67 and 70 
   // q1 += "AND MU_MINUS.y1*MU_MINUS.y3 > 0. ";
   // q1 += "AND MU_PLUS.y1*MU_PLUS.y3 > 0. ";
   // q1 += "AND MU_PLUS.y3*MU_MINUS.y3 < 0. ";
   // q1 += "AND D1<400. ";
   // q1 += "AND D2<400. ";
   // q1 += "AND D3<400.";
   //   q1 += "limit 200" ;                              
   //   q1 += "AND triggerCount > 0.";
   //  sometimes this is a negative value meaning there were some issues with the readout
   //  q1 += "AND turnOnset BETWEEN 0 AND 363000 ";
   //  There seem to be some events with turnOnset < 0 and this isnt really correct
   //   q1+= "AND Intensity_p < 60000.";
   
   //Something to keep in mind for dbar/ubar analysis is that turnOnset and rfOnset Kenichi's presentation
   //  use one apostrophe whenever sql requires "" in the query 
   return q1;
}

void extractData(int rs, int index){

   ostringstream oss;
   ifstream ifs;
   oss.str("");
   oss << "/seaquest/users/knagai/R008_analysis/list/R008_roadset_" << rs << "_" << index << ".list";
   ifs.open(oss.str().c_str());
   ifstream ifs_R007;
   oss.str("");
   oss << "/seaquest/users/knagai/R008_analysis/list/R007_roadset_" << rs << "_" << index << ".list";
   ifs_R007.open(oss.str().c_str());
   string host_tmp; TString run; //int rs;
   string host_tmp_R007; TString run_R007; int rs_R007;

   ofstream ofs;
   ofs.open("/seaquest/users/knagai/R008_analysis/list/runs_not_in_R008_1.txt");
   ofs << "RUN_LIST" << endl;

   while( ifs >> host_tmp >> run >> rs && ifs_R007 >> host_tmp_R007 >> run_R007 >> rs_R007 ){
      TString q1 = getQuery(run, run_R007);
      oss.str("");
      oss << "mysql://" << host_tmp << ":";
      if( strcmp(host_tmp.c_str(), "seaquestdb01.fnal.gov") == 0 ) oss << 3310;
      else                                                         oss << 3306;
//      oss << "/" << run;
      const char* host = oss.str().c_str();
      const char* user = "seaguest";
      const char* pass = "qqbar2mu+mu-";

      // cout << host << " / " << user << " / " << pass << " / " << run << endl;

      // cout << "mysql -u " << user << " -p" << pass << " -h " << host << endl;
 
      TSQLServer* db = TSQLServer::Connect(host, user, pass);
      TSQLResult* result;
      TSQLRow* row;

      if ( db->SelectDataBase(run.Data()) != 0 )
      {
         ofs << run << endl;
         continue;
      }  

      result = db->Query(q1.Data());

//filling up the Tree from one row at a time

      int n = result->GetRowCount();
      for (int i = 0; i < n; i++)
      {
         row = result->Next();
         Intensity_p = atof(row->GetField(0));
         trackID_mu_p = atoi(row->GetField(1));
         runID_mu_p = atoi(row->GetField(2));
         spillID_mu_p = atoi(row->GetField(3));
         eventID_mu_p = atoi(row->GetField(4));
         charge_mu_p = atoi(row->GetField(5));
         roadID_mu_p = atoi(row->GetField(6));
         numHits_mu_p = atoi(row->GetField(7));
         numHitsSt1_mu_p = atoi(row->GetField(8));
         numHitsSt2_mu_p = atoi(row->GetField(9));
         numHitsSt3_mu_p = atoi(row->GetField(10));
         numHitsSt4H_mu_p = atoi(row->GetField(11));
         numHitsSt4V_mu_p = atoi(row->GetField(12));
         chisq_mu_p = atof(row->GetField(13));
         x0_mu_p = atof(row->GetField(14));
         y0_mu_p = atof(row->GetField(15));
         z0_mu_p = atof(row->GetField(16));
         xD_mu_p = atof(row->GetField(17));
         yD_mu_p = atof(row->GetField(18));
         xT_mu_p = atof(row->GetField(19));
         yT_mu_p = atof(row->GetField(20));
         pxD_mu_p = atof(row->GetField(21));
         pyD_mu_p = atof(row->GetField(22));
         pzD_mu_p = atof(row->GetField(23));
         pxT_mu_p = atof(row->GetField(24));
         pyT_mu_p = atof(row->GetField(25));
         pzT_mu_p = atof(row->GetField(26));
         z0x_mu_p = atof(row->GetField(27));
         z0y_mu_p = atof(row->GetField(28));
         px0_mu_p = atof(row->GetField(29));
         py0_mu_p = atof(row->GetField(30));
         pz0_mu_p = atof(row->GetField(31));
         x1_mu_p = atof(row->GetField(32));
         y1_mu_p = atof(row->GetField(33));
         z1_mu_p = atof(row->GetField(34));
         px1_kTrack_mu_p = atof(row->GetField(35));
         py1_kTrack_mu_p = atof(row->GetField(36));
         pz1_kTrack_mu_p = atof(row->GetField(37));
         x3_mu_p = atof(row->GetField(38));
         y3_mu_p = atof(row->GetField(39));
         z3_mu_p = atof(row->GetField(40));
         px3_kTrack_mu_p = atof(row->GetField(41));
         py3_kTrack_mu_p = atof(row->GetField(42));
         pz3_kTrack_mu_p = atof(row->GetField(43));
         thbend_mu_p = atof(row->GetField(44));
         tx_PT_mu_p = atof(row->GetField(45));
         ty_PT_mu_p = atof(row->GetField(46));
         chisq_target_mu_p = atof(row->GetField(47));
         chisq_dump_mu_p = atof(row->GetField(48));
         chisq_upstream_mu_p = atof(row->GetField(49));




         trackID_mu_m = atoi(row->GetField(51));
         runID_mu_m = atoi(row->GetField(52));
         spillID_mu_m = atoi(row->GetField(53));
         eventID_mu_m = atof(row->GetField(54));
         charge_mu_m = atoi(row->GetField(55));
         roadID_mu_m = atoi(row->GetField(56));
         numHits_mu_m = atoi(row->GetField(57));
         numHitsSt1_mu_m = atoi(row->GetField(58));
         numHitsSt2_mu_m = atoi(row->GetField(59));
         numHitsSt3_mu_m = atoi(row->GetField(60));
         numHitsSt4H_mu_m = atoi(row->GetField(61));
         numHitsSt4V_mu_m = atoi(row->GetField(62));
         chisq_mu_m = atof(row->GetField(63));
         x0_mu_m = atof(row->GetField(64));
         y0_mu_m = atof(row->GetField(65));
         z0_mu_m = atof(row->GetField(66));
         xD_mu_m = atof(row->GetField(67));
         yD_mu_m = atof(row->GetField(68));
         xT_mu_m = atof(row->GetField(69));
         yT_mu_m = atof(row->GetField(70));
         pxD_mu_m = atof(row->GetField(71));
         pyD_mu_m = atof(row->GetField(72));
         pzD_mu_m = atof(row->GetField(73));
         pxT_mu_m = atof(row->GetField(74));
         pyT_mu_m = atof(row->GetField(75));
         pzT_mu_m = atof(row->GetField(76));
         z0x_mu_m = atof(row->GetField(77));
         z0y_mu_m = atof(row->GetField(78));
         px0_mu_m = atof(row->GetField(79));
         py0_mu_m = atof(row->GetField(80));
         pz0_mu_m = atof(row->GetField(81));
         x1_mu_m = atof(row->GetField(82));
         y1_mu_m = atof(row->GetField(83));
         z1_mu_m = atof(row->GetField(84));
         px1_kTrack_mu_m = atof(row->GetField(85));
         py1_kTrack_mu_m = atof(row->GetField(86));
         pz1_kTrack_mu_m = atof(row->GetField(87));
         x3_mu_m = atof(row->GetField(88));
         y3_mu_m = atof(row->GetField(89));
         z3_mu_m = atof(row->GetField(90));
         px3_kTrack_mu_m = atof(row->GetField(91));
         py3_kTrack_mu_m = atof(row->GetField(92));
         pz3_kTrack_mu_m = atof(row->GetField(93));
         thbend_mu_m = atof(row->GetField(94));
         tx_PT_mu_m = atof(row->GetField(95));
         ty_PT_mu_m = atof(row->GetField(96));
         chisq_target_mu_m = atof(row->GetField(97));
         chisq_dump_mu_m = atof(row->GetField(98));
         chisq_upstream_mu_m = atof(row->GetField(99));
         dimuonID = atoi(row->GetField(101));
         targetPos = atoi(row->GetField(105));
         posTrackID = atoi(row->GetField(106));
         negTrackID = atoi(row->GetField(107));
         dx = atof(row->GetField(108));
         dy = atof(row->GetField(109));
         dz = atof(row->GetField(110));
         dpx = atof(row->GetField(111));
         dpy = atof(row->GetField(112));
         dpz = atof(row->GetField(113));
         mass = atof(row->GetField(114));
         xF = atof(row->GetField(115));
         xB = atof(row->GetField(116));
         xT = atof(row->GetField(117));
         costh = atof(row->GetField(118));
         phi = atof(row->GetField(119));
         trackSeparation = atof(row->GetField(120));
         chisq_dimuon = atof(row->GetField(121));
         px1_kDimuon = atof(row->GetField(122));
         py1_kDimuon = atof(row->GetField(123));
         pz1_kDimuon = atof(row->GetField(124));
         px2_kDimuon = atof(row->GetField(125));
         py2_kDimuon = atof(row->GetField(126));
         pz2_kDimuon = atof(row->GetField(127));
         isValid = atof(row->GetField(128));
         isTarget = atof(row->GetField(129));
         isDump = atof(row->GetField(130));
	 D1 = atoi(row->GetField(134));
	 D2 = atoi(row->GetField(135));
	 D3 = atoi(row->GetField(136));
	 H1 = atoi(row->GetField(137));
	 H2 = atoi(row->GetField(138));
	 H3 = atoi(row->GetField(139));
	 H4 = atoi(row->GetField(140));
	 P1 = atoi(row->GetField(141));
	 P2 = atoi(row->GetField(142));
	 D1L = atoi(row->GetField(143));
	 D1R = atoi(row->GetField(144));
	 D2L = atoi(row->GetField(145));
	 D2R = atoi(row->GetField(146));
	 D3L = atoi(row->GetField(147));
	 D3R = atoi(row->GetField(148));
	 QIEsum = atof(row->GetField(149));
	 inh_thres = atof(row->GetField(150));
	 //	 trigger_sum_no_inhibit = atof(row->GetField(150));
	 //	 inhibit_block_sum = atof(row->GetField(151));
	 liveProton = atof(row->GetField(151));
	 RF00 = atof(row->GetField(152));
	 G2SEM = atof(row->GetField(153));
         tree->Fill();
	 
      }
      
   }
}

int main(int argc, char* argv[]){

   int roadset = atoi(argv[1]);
   int   index = atoi(argv[2]);
   
   prepareOutput(roadset, index);
   
   extractData(roadset, index);
   
   file->cd();
   tree->Write();
   file->Close();
}
