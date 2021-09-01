#include "pwmig/pwmigcore/PWMIGmigrated_seismogram.h"
namespace pwmig::pwmigcore {
PWMIGmigrated_seismogram::PWMIGmigrated_seismogram(const int i,
  const int j,const int n) : migrated_data(3,n),elog()
{
  ix1=i;
  ix2=j;
  domega.reserve(n);
  dweight.reserve(n);
  live=false;
}
PWMIGmigrated_seismogram::PWMIGmigrated_seismogram(const PWMIGmigrated_seismogram& parent)
  : migrated_data(parent.migrated_data),elog(parent.elog),domega(parent.domega),dweight(parent.dweight)
{
  ix1=parent.ix1;
  ix2=parent.ix2;
  live=parent.live;
}

PWMIGmigrated_seismogram& PWMIGmigrated_seismogram::operator=(const PWMIGmigrated_seismogram& parent)
{
  if(this != &parent)
  {
    ix1=parent.ix1;
    ix2=parent.ix2;
    live=parent.live;
    migrated_data=parent.migrated_data;
    elog=parent.elog;
    domega=parent.domega;
    dweight=parent.dweight;
  }
  return *this;
}
} // End namespace encapsulation
