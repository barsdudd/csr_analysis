#ifndef __SELECTION_H__
#define __SELECTION_H__

#include "Event.h"

using namespace std;

bool roughCuts(Event* event);

bool CAisSatisfied_2111_v32(Track trackP, Track trackN, int rs);
bool tracksAreValid_2111_v32( Track trackP, Track trackN , Dimuon dimuon);
bool trackIsValid_2111_v32(Track track, int rs);
bool dimuonIsValid_2111_v32(Dimuon dimuon, Event* event);

bool occupancyIsValid_2111_v42(Event* event);
bool tracksAreValid_2111_v42(Track trackP, Track trackN , Dimuon dimuon);
bool trackIsValid_2111_v42(Track track, int rs=67);
bool dimuonIsValid_2111_v42(Dimuon dimuon, int rs=67, bool looseMode=false);
// just setting int rs = 67 for the sake of running the code for individual targets which requires an integer value. If we use individual targets, we have to change it.
bool tightMode_2111_v42(Track trackP, Track trackN , Dimuon dimuon);
bool beamDump_2111_v42(Track trackP, Track trackN , Dimuon dimuon, int rs);

#endif // __SELECTION_H__
