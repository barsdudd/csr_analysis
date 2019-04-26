#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <fstream>
#include "Event.h"

bool Dimuon::operator<(const Dimuon& dimuon) const
{
   if( dimuonID <  dimuon.dimuonID) 
      return true; 
   return false;
}

bool Dimuon::operator==(const Dimuon& dimuon) const
{
   if( dimuonID == dimuon.dimuonID) 
      return true; 
   return false;
}

bool Track::operator<(const Track& track) const
{
   if( trackID <  track.trackID) 
      return true; 
   return false;
}

bool Track::operator==(const Track& track) const
{
   if( trackID == track.trackID) 
      return true;
   return false;
}
