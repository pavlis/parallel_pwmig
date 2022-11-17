#include <omp.h>
// used only for debug - remove when finished
//#include <iostream>
//using namespace std;
#include "mspass/utility/Metadata.h"
#include "mspass/seismic/Seismogram.h"
#include "mspass/seismic/Ensemble.h"
#include "pwmig/gclgrid/gclgrid.h"
#include "pwmig/seispp/VelocityModel_1d.h"
#include "pwmig/pwmigcore/PWMIGmigrated_seismogram.h"
#include "pwmig/pwmigcore/pwmig.h"
#include "pwmig/pwmigcore/PWMIGmigrated_seismogram.h"
#include "pwmig/pwmigcore/SlownessVectorMatrix.h"

namespace pwmig::pwmigcore
{
using namespace std;
using namespace pwmig::gclgrid;
using namespace pwmig::seispp;
using namespace pwmig::pwmigcore;
using namespace mspass::seismic;
using namespace mspass::utility;

PWMIGfielddata migrate_component(ThreeComponentEnsemble& d,
  GCLgrid& parent,
    GCLscalarfield3d& TPgrid,
      SlownessVectorMatrix& VPsvm,
         GCLscalarfield3d& Us3d,
           VelocityModel_1d& Vp1d,
             VelocityModel_1d& Vs1d,
               Metadata& control)
{

  int nmembers;
  nmembers = d.member.size();
  /* Perhaps should throw an exception here, but an empty ensemble should
  be allowed and handled seamlessly*/
  if(nmembers==0)
      return PWMIGfielddata();
  double tmax,zmax,dt;
  const double VPVSmax(2.0);
  zmax=control.get_double("maximum_depth");
  tmax=control.get_double("maximum_time_lag");
  dt=control.get_double("data_sample_interval");
  GCLscalarfield3d *raygrid;
  raygrid = Build_GCLraygrid(parent,VPsvm,Vs1d,zmax,VPVSmax*tmax,dt*VPVSmax);
  PWMIGfielddata pwdgrid(*raygrid);
  //omp_set_num_threads(8);
  //#pragma omp parallel for
  for(int m=0;m<d.member.size();++m)
  {
    PWMIGmigrated_seismogram dout;
    dout = migrate_one_seismogram(d.member[m], parent, *raygrid, TPgrid,Us3d,
                      Vp1d, Vs1d, control);
    pwdgrid.accumulate(dout);
  }
  delete raygrid;
  return pwdgrid;
}
}  //end namespace
