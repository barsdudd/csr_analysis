#ifndef __EVENT_H__
#define __EVENT_H__
#include <TVector3.h>

using namespace std;

class Dimuon {
public:
   int    dimuonID;
   double     mass;
   double       x1;
   double       x2;
   double       xF;
   double    costh;
   double      phi;
   double chisq_dimuon;
   double trackSeparation;

   TVector3 vtx_pos;
   TVector3 vtx_mom;

   TVector3 st1_pos;
   TVector3 st1_mom;
   TVector3 st3_pos;
   TVector3 st3_mom;

   int trackID_pos;
   int trackID_neg;
   int trackIndex_pos;
   int trackIndex_neg;

   bool isValid_2111_v32;
   bool CAisValid_2111_v32;

   bool operator<( const Dimuon& dimuon) const;
   bool operator==(const Dimuon& dimuon) const;

Dimuon() : mass(-1), x1(-1), x2(-1), isValid_2111_v32(false), CAisValid_2111_v32(false) {;}
   virtual ~Dimuon(){;}

   ClassDef(Dimuon, 1);
};

class Track {
public:
   int     trackID;
   int      roadID;
   int      charge;

   int     numHits;
   int     numHitsInSt[5];

   double chisq;
   double chisq_target;
   double chisq_dump;
   double chisq_upstream;

   TVector3  posVtx;
   TVector3  momVtx;
   TVector3  posSt[3];
   TVector3  momSt[3];

   TVector3 posTarg;
   TVector3 momTarg;
   TVector3 posDump;
   TVector3 momDump;

   double z0x;
   double z0y;

   double tx_PT;
   double ty_PT;

   bool isValid_2111_v32;

   bool operator<( const Track& track) const;
   bool operator==(const Track& track) const;

   Track(){ roadID = -1; posSt[0] = TVector3(-1, -1, -1); posSt[1] = TVector3(-1, -1, -1); posSt[3] = TVector3(-1, -1, -1); isValid_2111_v32 = false;}
   virtual ~Track(){;}
   ClassDef(Track, 1);
};

class Event {
public:
   int   roadset;
   int   spillID;
   int     runID;
   int   eventID;

   int targetPos;

   int occChams[3];
   int occChamsLR[3][2];
   int occHodos[4];
   int occProps[2];
   int RF[33];
   int inh_thres;
   double inte_t[33];
   double inte_p;
   double QIEsum;
   double G2SEM;
   double PotPerQie;
   vector<Dimuon> dimuons;
   vector<Track>   tracks;

   int occTotal(){
      return occChams[0] + occChams[1] + occChams[2];
   };

   void clear(){
      roadset = -1; eventID = -1; spillID = -1 ; targetPos = -1; runID = -1;
      dimuons.clear(); tracks.clear();
   };
   
   Event() { runID = -1; eventID = -1; spillID = -1 ; targetPos = -1; }
   virtual ~Event(){;}
   ClassDef(Event, 1);
};

#endif // __EVENT_H__
