#include <stdlib.h>
#include <cstdlib>
#include <iostream>
#include "Variables.h"
#include "Selection.h"

using namespace std;

bool roughCuts(Event* event){
   if( event->dimuons.size() < 1 ) return false;
   bool validity = true;
   for( int id = 0 ; id < (int)event->dimuons.size() ; id++ ){
      Dimuon dimuon = event->dimuons[id];
      bool trackIsValid[2] = {false, false};
      for( int it = 0 ; it < (int)event->tracks.size() ; it++ ){
         Track track = event->tracks[it];
         if( track.trackID == dimuon.trackID_pos ) trackIsValid[0] = true;
         if( track.trackID == dimuon.trackID_neg ) trackIsValid[1] = true;
      }
      validity = validity && trackIsValid[0] && trackIsValid[1];
   }
   return validity;
}

bool CAisSatisfied_2111_v32(Track trackP, Track trackN, int rs){
  double least = trackP.momSt[0].Z() <  trackN.momSt[0].Z() ? trackP.momSt[0].Z() :  trackN.momSt[0].Z();
  double great = trackP.posSt[0].X() > -trackN.posSt[0].X() ? trackP.posSt[0].X() : -trackN.posSt[0].X();
  if( rs < 67) great = -trackP.posSt[0].X() > trackN.posSt[0].X() ? -trackP.posSt[0].X() : trackN.posSt[0].X();
//  cout << least - great / 4.5 << " | " << least << " | " << great << endl;
  if( ! ( least - great / 4.5 > 13 ) ) return false;
  return true;
}

bool tracksAreValid_2111_v32( Track trackP, Track trackN , Dimuon dimuon){

  if( ! ( trackP.posSt[2].Y() * trackN.posSt[2].Y() < 0 ) ) return false;
  if( ! ( trackP.numHits + trackN.numHits > 29 ) ) return false;
  if( ! ( trackP.numHitsInSt[0] + trackN.numHitsInSt[0] > 8 ) ) return false;
  if( ! ( fabs( trackP.chisq_target + trackN.chisq_target - dimuon.chisq_dimuon ) < 2 ) ) return false;
  return true;
}

bool trackIsValid_2111_v32(Track track, int rs){
  //  return track.isValid;  // STD
  //  cout << fabs( track.momSt[0].X() - track.momSt[2].X() + track.chisq * .416 ) << endl;
   if( ! ( track.chisq_target < 15  ) ) return false;
      
   if( ! ( track.momSt[0].Z() >  9 ) ) return false;
   if( ! ( track.momSt[0].Z() < 75 ) ) return false;
         
   if( ! ( track.numHits > 13 ) ) return false;
   if( ! ( track.posTarg.X() * track.posTarg.X() + (track.posTarg.Y() - 1.6) * (track.posTarg.Y() - 1.6) <  400 ) ) return false;
   if( ! ( track.posDump.X() * track.posDump.X() + (track.posDump.Y() - 1.6) * (track.posDump.Y() - 1.6) >   12 ) ) return false;
   if( ! ( track.posDump.X() * track.posDump.X() + (track.posDump.Y() - 1.6) * (track.posDump.Y() - 1.6) < 1200 ) ) return false;
   if( rs == 67 && track.charge == 1 && ! ( track.roadID < 55000 ) ) return false;
   if( ! ( track.chisq_target < 1.5 * track.chisq_upstream ) ) return false;
   if( ! ( track.chisq_target < 1.5 * track.chisq_dump     ) ) return false;
   if( ! ( track.posVtx.Z() <    0 ) ) return false;
   if( ! ( track.posVtx.Z() > -350 ) ) return false;
   if( ! ( track.chisq / ( track.numHits - 5. ) < 13 ) ) return false;
   if( ! ( fabs( track.momSt[0].X() - track.momSt[2].X() + track.charge * 0.416 ) < 0.008 ) ) return false;
   if( ! ( fabs( track.momSt[0].Y() - track.momSt[2].Y() ) < 0.008 ) ) return false;
   if( ! ( fabs( track.momSt[0].Z() - track.momSt[2].Z() ) < 0.08  ) ) return false;
   if( ! ( track.posSt[0].Y() * track.posSt[2].Y() > 0 ) ) return false;
   return true;
   //  return trackIsValidSimple(track);  // SIMPLE
  //  return true;   // SIMPLEST (NO SELECTION)
}

bool dimuonIsValid_2111_v32(Dimuon dimuon, Event* event){
   Track trackP = event->tracks[dimuon.trackIndex_pos];
   Track trackN = event->tracks[dimuon.trackIndex_neg];
   double trackSeparation = trackP.posVtx.Z() - trackN.posVtx.Z();
     if( ! ( fabs(dimuon.vtx_pos.X()    ) <    0.30 ) ) return false;
     if( ! ( fabs(dimuon.vtx_pos.Y()-1.6) <    0.22 ) ) return false;
     if( ! ( dimuon.vtx_pos.Z()           > -280.00 ) ) return false;
     if( ! ( dimuon.vtx_pos.Z()           < -  5.00 ) ) return false;
     if( ! ( fabs(dimuon.vtx_mom.X()    ) <    2.00 ) ) return false;
     if( ! ( fabs(dimuon.vtx_mom.Y()    ) <    2.00 ) ) return false;
     if( ! ( dimuon.vtx_mom.X() * dimuon.vtx_mom.X() + dimuon.vtx_mom.Y() * dimuon.vtx_mom.Y() < 5.5 ) ) return false;
     if( ! ( dimuon.vtx_mom.Z()           >   37.00 ) ) return false;
     if( ! ( dimuon.vtx_mom.Z()           <  116.00 ) ) return false;
     if( ! ( dimuon.mass                  <    8.80 ) ) return false;
     if( ! ( dimuon.mass                  >    4.20 ) ) return false;
     if( ! ( dimuon.vtx_pos.X() * dimuon.vtx_pos.X() + (dimuon.vtx_pos.Y()-1.6) * (dimuon.vtx_pos.Y()-1.6) < 0.09 ) ) return false;
     if( ! ( dimuon.xF                    <    0.95 ) ) return false;
     if( ! ( dimuon.xF                    > -  0.15 ) ) return false;
     // if( ! ( dimuon.x2                    >    0.05 ) ) return false;
     // if( ! ( dimuon.x2                    <    0.50 ) ) return false;
     if( ! ( fabs(dimuon.costh)           <    0.50 ) ) return false;
     if( ! ( fabs(trackSeparation       ) <  280.00 ) ) return false;
     if( ! ( dimuon.chisq_dimuon          <   18    ) ) return false;
     return true;
  //  return true;  // SIMPLE
}

bool occupancyIsValid_2111_v42(Event* event){
   for( int id = 0 ; id < 3 ; id++ )
      if( event->occChams[id] >= 400 ) return false;
   if( event->occTotal() >= 1000 ) return false;
   return true;
}

bool dimuonIsValid_2111_v42(Dimuon dimuon, int rs, bool looseMode){
    double beamOffset = getBeamOffset(rs);
   if( ! ( fabs(dimuon.vtx_pos.X()    ) < .25 ) ) return false;
   if( ! ( fabs(dimuon.vtx_pos.Y()-beamOffset) < .22 ) ) return false;
   if( ! (      dimuon.vtx_pos.Z()      < -  5) ) return false;
   if( ! (      dimuon.vtx_pos.Z()      > -280) ) return false;
   if( ! ( fabs(dimuon.vtx_mom.X()    ) < 1.8 ) ) return false;
   if( ! ( fabs(dimuon.vtx_mom.Y()    ) < 2.0 ) ) return false;
   if( ! ( fabs(dimuon.costh          ) < 0.5 ) ) return false;
   if( ! (      dimuon.vtx_mom.Z()      < 116 ) ) return false;
   if( ! (      dimuon.vtx_mom.Z()      >  38 ) ) return false;

   if( ! ( pow(dimuon.vtx_mom.X(), 2) + pow(dimuon.vtx_mom.Y()    , 2) < 5  ) ) return false;
      if( ! ( pow(dimuon.vtx_pos.X(), 2) + pow(dimuon.vtx_pos.Y()-beamOffset, 2) < .06) ) return false;

   if( ! ( dimuon.mass < 8.8 ) ) return false;
   if( !looseMode ){
      if( ! ( dimuon.mass > 4.2 ) ) return false;
   }else{
      if( ! ( dimuon.mass > 4.1 ) ) return false;
   }


   if( ! ( dimuon.xF <  .95 ) ) return false;
   if( !looseMode ){
      if( ! ( dimuon.xF > -.10 ) ) return false;
   }else{
      if( ! ( dimuon.xF > -.10 ) ) return false;
   }
   if( ! ( dimuon.x2 >  .05 ) ) return false;
   if( ! ( dimuon.x2 <  .55 ) ) return false;

   if( ! ( fabs(dimuon.costh          ) <  .5 ) ) return false;
   if( ! ( fabs(dimuon.trackSeparation) < 270 ) ) return false;

   if( ! ( dimuon.chisq_dimuon < 18 ) ) return false;

   return true;
}

bool trackIsValid_2111_v42(Track track, int rs){
  double beamOffset = getBeamOffset(rs);
   if( ! ( track.chisq_target < 15 ) ) return false;
  
   if( ! ( track.momSt[0].Z() >  9 ) ) return false;
   if( ! ( track.momSt[0].Z() < 75 ) ) return false;

   if( ! ( track.numHits > 13 ) ) return false;

   if( ! ( pow(track.posTarg.X(), 2) + pow(track.posTarg.Y()-beamOffset, 2) <  320 ) ) return false;
   if( ! ( pow(track.posDump.X(), 2) + pow(track.posDump.Y()-beamOffset, 2) < 1100 ) ) return false;
   if( ! ( pow(track.posDump.X(), 2) + pow(track.posDump.Y()-beamOffset, 2) >   16 ) ) return false;

   if( ! ( track.chisq_target < 1.5 * track.chisq_upstream ) ) return false;
   if( ! ( track.chisq_target < 1.5 * track.chisq_dump     ) ) return false;

   if( ! ( track.posVtx.Z() < -  5 ) ) return false;
   if( ! ( track.posVtx.Z() > -320 ) ) return false;

   if( ! ( track.chisq / ( track.numHits - 5 ) < 12 ) ) return false;

   if( ! ( track.posSt[0].Y() / track.posSt[2].Y() < 1. ) ) return false;

   if( ! fabs(( fabs( track.momSt[0].X() - track.momSt[2].X() ) - .416) < .008 ) ) return false;
   if( ! ( fabs( track.momSt[0].Y() - track.momSt[2].Y() )        < .008 ) ) return false;
   if( ! ( fabs( track.momSt[0].Z() - track.momSt[2].Z() )        < .08  ) ) return false;

   if( ! ( track.posSt[0].Y() * track.posSt[2].Y() > 0. ) ) return false;

   if( ! ( fabs( track.momSt[0].Y() ) > .02 ) ) return false;

   return true;
}

bool tracksAreValid_2111_v42(Track trackP, Track trackN, Dimuon dimuon){
   if( ! ( fabs( trackP.chisq_target + trackN.chisq_target - dimuon.chisq_dimuon ) < 2 ) ) return false;

   if( ! ( trackP.posSt[2].Y() * trackN.posSt[2].Y() < 0 ) ) return false;

   if( ! ( trackP.numHits        + trackN.numHits        > 29 ) ) return false;
   if( ! ( trackP.numHitsInSt[0] + trackN.numHitsInSt[0] >  8 ) ) return false;
   if( ! ( fabs( trackP.posSt[0].X() + trackN.posSt[0].X() ) < 42 ) ) return false; 
   // this should not be in the loose mode

   return true;
}

bool tightMode_2111_v42(Track trackP, Track trackN, Dimuon dimuon){
   if( dimuon.mass        <= 4.3 ) return false;
   if( dimuon.vtx_pos.Z() >= -60 ) return false;
   if( trackP.chisq_target >= trackP.chisq_dump ) return false;   
   if( trackN.chisq_target >= trackN.chisq_dump ) return false;   
   if( fabs(trackP.chisq_target + trackN.chisq_target - dimuon.chisq_dimuon) >= 2 ) return false;
   if( ! ( trackP.posSt[0].X() + trackN.posSt[0].X() < 32 ) ) return false;
   return true;
}

bool beamDump_2111_v42(Track trackP, Track trackN , Dimuon dimuon, int rs){
   if( dimuon.vtx_pos.Z() >= -60 ) return false;
   if( trackP.chisq_target >= 0.6*trackP.chisq_dump ) return false;   
   if( trackN.chisq_target >= 0.6*trackN.chisq_dump ) return false;   
   if( ! ( trackP.posSt[0].X() + trackN.posSt[0].X() < 42 ) ) return false;
   double beamOffset = getBeamOffset(rs);

   if( ! ( pow(trackP.posTarg.X(), 2) + pow(trackP.posTarg.Y()-beamOffset, 2) <  280 ) ) return false;
   if( ! ( pow(trackP.posDump.X(), 2) + pow(trackP.posDump.Y()-beamOffset, 2) < 900 ) ) return false;
   if( ! ( pow(trackP.posDump.X(), 2) + pow(trackP.posDump.Y()-beamOffset, 2) >   36 ) ) return false;

   if( ! ( pow(trackN.posTarg.X(), 2) + pow(trackN.posTarg.Y()-beamOffset, 2) <  280 ) ) return false;
   if( ! ( pow(trackN.posDump.X(), 2) + pow(trackN.posDump.Y()-beamOffset, 2) < 900 ) ) return false;
   if( ! ( pow(trackN.posDump.X(), 2) + pow(trackN.posDump.Y()-beamOffset, 2) >   36 ) ) return false;

   return true;
}
