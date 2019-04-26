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
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include "Event.h"
#include "Selection.h"

#define FILE_EXIST 1
#define DIR_EXIST 2
#define NO_FILE_DIR_EXIST 0

using namespace std;

TFile* dataFile;
TTree* dataTree;

TFile* saveFile;
TTree* saveTree;

Event* event;
char* outDir;

vector<int> spillList;
int gl_is = 0;
int   rs;
int targ;
const char* targName[] = {"All", "LH2", "Empty", "LD2", "None", "Fe", "C", "W"};

const char* filePlace = "/seaquest/users/arunts/data/";
const char* file       = "/seaquest/users/knagai/R008_analysis/data/";

double nbuckets =    588;
double   nturns = 369000;
map<int, double> pedestal;

int runID_prev, eventID_prev;

int dimuonID, runID, spillID, eventID, posTrackID, negTrackID, targetPos;

double dx, dy, dz, dpx, dpy, dpz, mass, xF, xB, xT, costh, phi, trackSeparation, chisq_dimuon, px1_kDimuon, py1_kDimuon, pz1_kDimuon, px2_kDimuon, py2_kDimuon, pz2_kDimuon, isValid, isTarget, isDump; 
    
//declaring all the variables in kTrack field for mu_p
int trackID_mu[2], runID_mu[2], eventID_mu[2], spillID_mu[2], roadID_mu[2], charge_mu[2], numHits_mu[2], numHitsSt1_mu[2], numHitsSt2_mu[2], numHitsSt3_mu[2], numHitsSt4H_mu[2], numHitsSt4V_mu[2]; 

double chisq_mu[2], x0_mu[2], y0_mu[2], z0_mu[2], xD_mu[2], yD_mu[2], xT_mu[2], yT_mu[2], pxD_mu[2], pyD_mu[2], pzD_mu[2], pxT_mu[2], pyT_mu[2], pzT_mu[2], z0x_mu[2], z0y_mu[2], px0_mu[2], py0_mu[2], pz0_mu[2], x1_mu[2], y1_mu[2], z1_mu[2], px1_kTrack_mu[2], py1_kTrack_mu[2], pz1_kTrack_mu[2], x3_mu[2], y3_mu[2], z3_mu[2], px3_kTrack_mu[2], py3_kTrack_mu[2], pz3_kTrack_mu[2], thbend_mu[2], tx_PT_mu[2], ty_PT_mu[2], chisq_target_mu[2], chisq_dump_mu[2], chisq_upstream_mu[2];

// //declaring all the variables in kTrack field for mu_m
// int trackID_mu_m,runID_mu_m, eventID_mu_m, spillID_mu_m, roadID_mu_m, charge_mu_m, numHits_mu_m, numHitsSt1_mu_m, numHitsSt2_mu_m, numHitsSt3_mu_m, numHitsSt4H_mu_m, numHitsSt4V_mu_m; 

// double chisq_mu_m, x0_mu_m, y0_mu_m, z0_mu_m, xD_mu_m, yD_mu_m, xT_mu_m, yT_mu_m, pxD_mu_m, pyD_mu_m, pzD_mu_m, pxT_mu_m, pyT_mu_m, pzT_mu_m, z0x_mu_m, z0y_mu_m, px0_mu_m, py0_mu_m, pz0_mu_m, x1_mu_m, y1_mu_m, z1_mu_m, px1_kTrack_mu_m, py1_kTrack_mu_m, pz1_kTrack_mu_m, x3_mu_m, y3_mu_m, z3_mu_m, px3_kTrack_mu_m, py3_kTrack_mu_m, pz3_kTrack_mu_m, thbend_mu_m, tx_PT_mu_m, ty_PT_mu_m, chisq_target_mu_m, chisq_dump_mu_m, chisq_upstream_mu_m;

//declaring all the occupancy variables
int D1, D2, D3, H1, H2, H3, H4, P1, P2, D1L, D1R, D2L, D2R, D3L, D3R;
  
//declaring all the variables in QIE 
double  Intensity_p, PotPerQie;

//declaring all the variables in BeamDAQ
double  QIEsum, liveProton, RF00, G2SEM, inh_thres;


//////////////
//   INIT   //
//////////////
void initInputFile(int rs){
   ostringstream oss;
   oss.str("");
   oss << filePlace << rs << "/" << rs << "_with_cuts.root";
   dataFile = new TFile(oss.str().c_str(), "READ");
   oss.str("");
   oss << "roadset_" << rs;
   dataTree = (TTree*)dataFile->Get(oss.str().c_str());

   dataTree->SetBranchAddress("Intensity_p",&Intensity_p);

   dataTree->SetBranchAddress("trackID_mu_p",&trackID_mu[0]);
   dataTree->SetBranchAddress("charge_mu_p",&charge_mu[0]);
   dataTree->SetBranchAddress("runID_mu_p",&runID_mu[0]);
   dataTree->SetBranchAddress("spillID_mu_p",&spillID_mu[0]);
   dataTree->SetBranchAddress("eventID_mu_p", &eventID_mu[0]);
   dataTree->SetBranchAddress("roadID_mu_p", &roadID_mu[0]);
   dataTree->SetBranchAddress("numHits_mu_p", &numHits_mu[0]);
   dataTree->SetBranchAddress("numHitsSt1_mu_p",&numHitsSt1_mu[0]);
   dataTree->SetBranchAddress("numHitsSt2_mu_p",&numHitsSt2_mu[0]);
   dataTree->SetBranchAddress("numHitsSt3_mu_p", &numHitsSt3_mu[0]);
   dataTree->SetBranchAddress("numHitsSt4H_mu_p", &numHitsSt4H_mu[0]);
   dataTree->SetBranchAddress("numHitsSt4V_mu_p", &numHitsSt4V_mu[0]);
   dataTree->SetBranchAddress("chisq_mu_p", &chisq_mu[0]);
   dataTree->SetBranchAddress("x0_mu_p", &x0_mu[0]);
   dataTree->SetBranchAddress("y0_mu_p", &y0_mu[0]);
   dataTree->SetBranchAddress("z0_mu_p", &z0_mu[0]);
   dataTree->SetBranchAddress("xD_mu_p", &xD_mu[0]);
   dataTree->SetBranchAddress("yD_mu_p", &yD_mu[0]);
   dataTree->SetBranchAddress("xT_mu_p", &xT_mu[0]);
   dataTree->SetBranchAddress("yT_mu_p", &yT_mu[0]);
   dataTree->SetBranchAddress("pxD_mu_p", &pxD_mu[0]);
   dataTree->SetBranchAddress("pyD_mu_p", &pyD_mu[0]);
   dataTree->SetBranchAddress("pzD_mu_p", &pzD_mu[0]);
   dataTree->SetBranchAddress("pxT_mu_p", &pxT_mu[0]);
   dataTree->SetBranchAddress("pyT_mu_p", &pyT_mu[0]);
   dataTree->SetBranchAddress("pzT_mu_p", &pzT_mu[0]);
   dataTree->SetBranchAddress("z0x_mu_p", &z0x_mu[0]);
   dataTree->SetBranchAddress("z0y_mu_p", &z0y_mu[0]);
   dataTree->SetBranchAddress("px0_mu_p", &px0_mu[0]);
   dataTree->SetBranchAddress("py0_mu_p", &py0_mu[0]);
   dataTree->SetBranchAddress("pz0_mu_p", &pz0_mu[0]);
   dataTree->SetBranchAddress("x1_mu_p", &x1_mu[0]);
   dataTree->SetBranchAddress("y1_mu_p", &y1_mu[0]);
   dataTree->SetBranchAddress("z1_mu_p", &z1_mu[0]);
   dataTree->SetBranchAddress("px1_kTrack_mu_p", &px1_kTrack_mu[0]);
   dataTree->SetBranchAddress("py1_kTrack_mu_p", &py1_kTrack_mu[0]);
   dataTree->SetBranchAddress("pz1_kTrack_mu_p", &pz1_kTrack_mu[0]);
   dataTree->SetBranchAddress("x3_mu_p", &x3_mu[0]);
   dataTree->SetBranchAddress("y3_mu_p", &y3_mu[0]);
   dataTree->SetBranchAddress("z3_mu_p", &z3_mu[0]);
   dataTree->SetBranchAddress("px3_kTrack_mu_p", &px3_kTrack_mu[0]);
   dataTree->SetBranchAddress("py3_kTrack_mu_p", &py3_kTrack_mu[0]);
   dataTree->SetBranchAddress("pz3_kTrack_mu_p", &pz3_kTrack_mu[0]);
   dataTree->SetBranchAddress("thbend_mu_p", &thbend_mu[0]);
   dataTree->SetBranchAddress("tx_PT_mu_p", &tx_PT_mu[0]);
   dataTree->SetBranchAddress("ty_PT_mu_p", &ty_PT_mu[0]);
   dataTree->SetBranchAddress("chisq_target_mu_p", &chisq_target_mu[0]);
   dataTree->SetBranchAddress("chisq_dump_mu_p", &chisq_dump_mu[0]);
   dataTree->SetBranchAddress("chisq_upstream_mu_p", &chisq_upstream_mu[0]);
   
   
   dataTree->SetBranchAddress("trackID_mu_m",&trackID_mu[1]);
   dataTree->SetBranchAddress("runID_mu_m",&runID_mu[1]);
   dataTree->SetBranchAddress("spillID_mu_m",&spillID_mu[1]);
   dataTree->SetBranchAddress("eventID_mu_m", &eventID_mu[1]);
   dataTree->SetBranchAddress("charge_mu_m",&charge_mu[1]);
   dataTree->SetBranchAddress("roadID_mu_m",&roadID_mu[1]);
   dataTree->SetBranchAddress("numHits_mu_m", &numHits_mu[1]);
   dataTree->SetBranchAddress("numHitsSt1_mu_m", &numHitsSt1_mu[1]);
   dataTree->SetBranchAddress("numHitsSt2_mu_m", &numHitsSt2_mu[1]);
   dataTree->SetBranchAddress("numHitsSt3_mu_m", &numHitsSt3_mu[1]);
   dataTree->SetBranchAddress("numHitsSt4H_mu_m", &numHitsSt4H_mu[1]);
   dataTree->SetBranchAddress("numHitsSt4V_mu_m", &numHitsSt4V_mu[1]);
   dataTree->SetBranchAddress("chisq_mu_m", &chisq_mu[1]);
   dataTree->SetBranchAddress("x0_mu_m", &x0_mu[1]);
   dataTree->SetBranchAddress("y0_mu_m", &y0_mu[1]);
   dataTree->SetBranchAddress("z0_mu_m", &z0_mu[1]);
   dataTree->SetBranchAddress("xD_mu_m", &xD_mu[1]);
   dataTree->SetBranchAddress("yD_mu_m", &yD_mu[1]);
   dataTree->SetBranchAddress("xT_mu_m", &xT_mu[1]);
   dataTree->SetBranchAddress("yT_mu_m", &yT_mu[1]);
   dataTree->SetBranchAddress("pxD_mu_m", &pxD_mu[1]);
   dataTree->SetBranchAddress("pyD_mu_m", &pyD_mu[1]);
   dataTree->SetBranchAddress("pzD_mu_m", &pzD_mu[1]);
   dataTree->SetBranchAddress("pxT_mu_m", &pxT_mu[1]);
   dataTree->SetBranchAddress("pyT_mu_m", &pyT_mu[1]);
   dataTree->SetBranchAddress("pzT_mu_m", &pzT_mu[1]);
   dataTree->SetBranchAddress("z0x_mu_m", &z0x_mu[1]);
   dataTree->SetBranchAddress("z0y_mu_m", &z0y_mu[1]);
   dataTree->SetBranchAddress("px0_mu_m", &px0_mu[1]);
   dataTree->SetBranchAddress("py0_mu_m", &py0_mu[1]);
   dataTree->SetBranchAddress("pz0_mu_m", &pz0_mu[1]);
   dataTree->SetBranchAddress("x1_mu_m", &x1_mu[1]);
   dataTree->SetBranchAddress("y1_mu_m", &y1_mu[1]);
   dataTree->SetBranchAddress("z1_mu_m", &z1_mu[1]);
   dataTree->SetBranchAddress("px1_kTrack_mu_m", &px1_kTrack_mu[1]);
   dataTree->SetBranchAddress("py1_kTrack_mu_m", &py1_kTrack_mu[1]);
   dataTree->SetBranchAddress("pz1_kTrack_mu_m", &pz1_kTrack_mu[1]);
   dataTree->SetBranchAddress("x3_mu_m", &x3_mu[1]);
   dataTree->SetBranchAddress("y3_mu_m", &y3_mu[1]);
   dataTree->SetBranchAddress("z3_mu_m", &z3_mu[1]);
   dataTree->SetBranchAddress("px3_kTrack_mu_m", &px3_kTrack_mu[1]);
   dataTree->SetBranchAddress("py3_kTrack_mu_m", &py3_kTrack_mu[1]);
   dataTree->SetBranchAddress("pz3_kTrack_mu_m", &pz3_kTrack_mu[1]);
   dataTree->SetBranchAddress("thbend_mu_m", &thbend_mu[1]);
   dataTree->SetBranchAddress("tx_PT_mu_m", &tx_PT_mu[1]);
   dataTree->SetBranchAddress("ty_PT_mu_m", &ty_PT_mu[1]);
   dataTree->SetBranchAddress("chisq_target_mu_m", &chisq_target_mu[1]);
   dataTree->SetBranchAddress("chisq_dump_mu_m", &chisq_dump_mu[1]);
   dataTree->SetBranchAddress("chisq_upstream_mu_m", &chisq_upstream_mu[1]);
   
   
   dataTree->SetBranchAddress("dimuonID", &dimuonID);
   dataTree->SetBranchAddress("targetPos", &targetPos);
   dataTree->SetBranchAddress("posTrackID", &posTrackID);
   dataTree->SetBranchAddress("negTrackID", &negTrackID);
   dataTree->SetBranchAddress("dx", &dx);
   dataTree->SetBranchAddress("dy", &dy);
   dataTree->SetBranchAddress("dz", &dz);
   dataTree->SetBranchAddress("dpx", &dpx);
   dataTree->SetBranchAddress("dpy", &dpy);
   dataTree->SetBranchAddress("dpz", &dpz);
   dataTree->SetBranchAddress("mass", &mass);
   dataTree->SetBranchAddress("xF", &xF);
   dataTree->SetBranchAddress("xB", &xB);
   dataTree->SetBranchAddress("xT", &xT);
   dataTree->SetBranchAddress("costh", &costh);
   dataTree->SetBranchAddress("phi", &phi);
   dataTree->SetBranchAddress("trackSeparation", &trackSeparation);
   dataTree->SetBranchAddress("chisq_dimuon", &chisq_dimuon);
   dataTree->SetBranchAddress("px1_kDimuon", &px1_kDimuon);
   dataTree->SetBranchAddress("py1_kDimuon", &py1_kDimuon);
   dataTree->SetBranchAddress("pz1_kDimuon", &pz1_kDimuon);
   dataTree->SetBranchAddress("px2_kDimuon", &px2_kDimuon);
   dataTree->SetBranchAddress("py2_kDimuon", &py2_kDimuon);
   dataTree->SetBranchAddress("pz2_kDimuon", &pz2_kDimuon);
   dataTree->SetBranchAddress("isValid", &isValid);
   dataTree->SetBranchAddress("isTarget", &isTarget);
   dataTree->SetBranchAddress("isDump", &isDump);
   
   dataTree->SetBranchAddress("D1", &D1);
   dataTree->SetBranchAddress("D2", &D2);
   dataTree->SetBranchAddress("D3", &D3);
   dataTree->SetBranchAddress("H1", &H1);
   dataTree->SetBranchAddress("H2", &H2);
   dataTree->SetBranchAddress("H3", &H3);
   dataTree->SetBranchAddress("H4", &H4);
   dataTree->SetBranchAddress("P1", &P1);
   dataTree->SetBranchAddress("P2", &P2);
   dataTree->SetBranchAddress("D1L", &D1L);
   dataTree->SetBranchAddress("D1R", &D1R);
   dataTree->SetBranchAddress("D2L", &D2L);
   dataTree->SetBranchAddress("D2R", &D2R);
   dataTree->SetBranchAddress("D3L", &D3L);
   dataTree->SetBranchAddress("D3R", &D3R);
   dataTree->SetBranchAddress("QIEsum", &QIEsum);
   dataTree->SetBranchAddress("G2SEM", &G2SEM);
   dataTree->SetBranchAddress("liveProton",&liveProton);  
   dataTree->SetBranchAddress("RF00",&RF00);
   dataTree->SetBranchAddress("inh_thres",&inh_thres);

   dataTree->SetBranchAddress("PotPerQie", &PotPerQie);
   
}

void initOutputFile(int rs){
   ostringstream oss;
   oss.str("");
   oss << file << rs << "/"; 
   gSystem->mkdir(oss.str().c_str(), true);
   oss << rs << "_with_cuts.root";
   saveFile = new TFile(oss.str().c_str(), "RECREATE");
   oss.str("");
   oss << "roadset_" << rs;
   saveTree = new TTree("save", "save");   
   saveTree->Branch("event", &event);
}

void initMain(int argc, char* argv[]){
   if( argc < 2 ){
      cout << "!! NEED 1 ARGUMENT !! : ROADSET" << endl;
      exit(0);
   }
   rs    = atoi(argv[1]);
   event = new Event();

   initInputFile (rs);
   initOutputFile(rs);

   pedestal[57] = 36.2;
   pedestal[59] = 36.2;
   pedestal[62] = 36.2;
   pedestal[67] = 32.6;
   pedestal[70] = 32.6;
}


/////////////////
//   ANALYZE   //
/////////////////
TVector3 getSt2Pos(TVector3 posSt3, TVector3 momSt3){
   double z = 1350.;
   double x = posSt3.X() + (z - posSt3.Z()) * momSt3.X() / momSt3.Z();
   double y = posSt3.Y() + (z - posSt3.Z()) * momSt3.Y() / momSt3.Z();
   return TVector3(x, y, z);
}

void fillEvent(){
   runID_prev   =   runID_mu[0];
   eventID_prev = eventID_mu[0];

   event->roadset = rs;
   event->spillID = spillID_mu[0];
   event->  runID =   runID_mu[0];
   event->eventID = eventID_mu[0];
   
   event->targetPos = targetPos;
   
   event->occChams[0] = D1;
   event->occChams[1] = D2;
   event->occChams[2] = D3;
   
   event->occChamsLR[0][0] = D1L;
   event->occChamsLR[1][0] = D2L;
   event->occChamsLR[2][0] = D3L;
   event->occChamsLR[0][1] = D1R;
   event->occChamsLR[1][1] = D2R;
   event->occChamsLR[2][1] = D3R;
   
   event->occHodos[0] = H1;
   event->occHodos[1] = H2;
   event->occHodos[2] = H3;
   event->occHodos[3] = H4;
   
   event->occProps[0] = 0;
   event->occProps[1] = 0;
   
   for( int irf = -16 ; irf <= 16 ; irf++ ){
      event->RF    [irf+16] = 0;
      event->inte_t[irf+16] = 0;
   }
   event->RF[16] = RF00;
   event->inte_t[16] = (RF00 - 34)*G2SEM/(QIEsum - nturns*34*nbuckets);
//   event->inte_t[16] = (RF00 - pedestal[rs])*PotPerQie;
   event->inh_thres = inh_thres;
   event->inte_p = Intensity_p;
   event->PotPerQie = PotPerQie;
   event->G2SEM = G2SEM;
   event->QIEsum = QIEsum;
}

int findTrackIndex(int trackID){
   for( int it = 0 ; it < (int)event->tracks.size() ; it++ ){
      Track track = event->tracks[it];
      if( trackID == track.trackID ) return it;
   }
   return -1;
}

Dimuon fillDimuon(){
   Dimuon dimuon;
   dimuon.dimuonID = dimuonID;
   dimuon.mass     =     mass;
   dimuon.x1       =       xB;
   dimuon.x2       =       xT;
   dimuon.xF       =       xF;
   dimuon.costh    =    costh;
   dimuon.phi      =      phi;

   dimuon.chisq_dimuon = chisq_dimuon;

   dimuon.vtx_pos = TVector3( dx,  dy,  dz);
   dimuon.vtx_mom = TVector3(dpx, dpy, dpz);

   dimuon.st1_pos = TVector3(0,  0,  0);
   dimuon.st1_mom = TVector3(0,  0,  0);
   dimuon.st3_pos = TVector3(0,  0,  0);
   dimuon.st3_mom = TVector3(0,  0,  0);

   dimuon.trackID_pos = posTrackID;
   dimuon.trackID_neg = negTrackID;
   dimuon.trackIndex_pos = findTrackIndex(posTrackID);
   dimuon.trackIndex_neg = findTrackIndex(negTrackID);

   if( dimuon.trackIndex_pos != -1 && dimuon.trackIndex_neg != -1 ){
      Track trackP = event->tracks[dimuon.trackIndex_pos];
      Track trackN = event->tracks[dimuon.trackIndex_neg];
      dimuon.CAisValid_2111_v32 = CAisSatisfied_2111_v32(trackP, trackN, rs);
      dimuon.  isValid_2111_v32 = dimuonIsValid_2111_v32(dimuon, event) && trackP.isValid_2111_v32 && trackN.isValid_2111_v32;
   }

   return dimuon;
}

Track fillTrack(int pm){
   Track track;
   
   track.trackID = trackID_mu[pm];
   track. roadID =  roadID_mu[pm];
   track. charge =  charge_mu[pm];

   track.numHits = numHits_mu[pm];

   track.numHitsInSt[0] = numHitsSt1_mu[pm];
   track.numHitsInSt[1] = numHitsSt2_mu[pm];
   track.numHitsInSt[2] = numHitsSt3_mu[pm];
   track.numHitsInSt[3] = numHitsSt4H_mu[pm];
   track.numHitsInSt[4] = numHitsSt4V_mu[pm];

   track.chisq          = chisq_mu         [pm];
   track.chisq_target   = chisq_target_mu  [pm];
   track.chisq_dump     = chisq_dump_mu    [pm];
   track.chisq_upstream = chisq_upstream_mu[pm];

   track.posVtx = TVector3( x0_mu[pm],  y0_mu[pm],  z0_mu[pm]);
   track.momVtx = TVector3(px0_mu[pm], py0_mu[pm], pz0_mu[pm]);

   track.posSt[0] = TVector3( x1_mu       [pm],  y1_mu       [pm],  z1_mu       [pm]);
   track.momSt[0] = TVector3(px1_kTrack_mu[pm], py1_kTrack_mu[pm], pz1_kTrack_mu[pm]);
   track.posSt[2] = TVector3( x3_mu       [pm],  y3_mu       [pm],  z3_mu       [pm]);
   track.momSt[2] = TVector3(px3_kTrack_mu[pm], py3_kTrack_mu[pm], pz3_kTrack_mu[pm]);

   track.posSt[1] = getSt2Pos(track.posSt[2], track.momSt[2]);
   track.momSt[1] = TVector3(px3_kTrack_mu[pm], py3_kTrack_mu[pm], pz3_kTrack_mu[pm]);

   track.posTarg = TVector3( xT_mu[pm],  yT_mu[pm],      -129);
   track.momTarg = TVector3(pxT_mu[pm], pyT_mu[pm], pzT_mu[pm]);
   track.posDump = TVector3( xD_mu[pm],  yD_mu[pm],        42);
   track.momDump = TVector3(pxD_mu[pm], pyD_mu[pm], pzD_mu[pm]);

   track.z0x = z0x_mu[pm];
   track.z0y = z0y_mu[pm];

   track.isValid_2111_v32 = trackIsValid_2111_v32(track, rs);

   return track;
}

void anaMain(){
   int nEvtMax = dataTree->GetEntries();
   int current =  0;
   int prev    = -1;
   for(int i = 0; i < nEvtMax ; ++i)
   {
      dataTree->GetEntry(i);
      current = (int)(double (i+1)/nEvtMax * 100);
      if( current != prev ){
         cout << "\r" << (i+1) << " / " << nEvtMax << " = " << (int)(double (i+1)/nEvtMax * 100) << " %" << flush;
         prev = (int)(double (i+1)/nEvtMax * 100);
      }

      if( i != 0 && ( runID_mu[0] != runID_prev || eventID_mu[0] != eventID_prev ) ){
         saveTree->Fill();
         if( saveTree->GetEntries() % 1000 == 1 ) saveTree->AutoSave("SaveSelf");
         delete event;
         event = new Event();
         fillEvent();
      }

      if( i == 0 ) fillEvent();

      Track  trackP = fillTrack(0);
      Track  trackN = fillTrack(1);
      Dimuon dimuon = fillDimuon();
      event->dimuons.push_back(dimuon);
      event-> tracks.push_back(trackP);
      event-> tracks.push_back(trackN);
   }
   saveTree->Fill();
   cout << endl;
}

////////////
//  MAIN  //
////////////
int main(int argc, char* argv[]){
   initMain(argc, argv);
   anaMain();
   saveFile->cd();
   saveTree->Write();
   saveFile->Close();
   return EXIT_SUCCESS;
}

