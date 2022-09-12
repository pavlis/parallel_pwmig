#include <thread>
#include <mutex>
#include <queue>
/* needed for ref */
#include <functional>
#include "mspass/utility/Metadata.h"
#include "mspass/utility/MsPASSError.h"
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
/* Mutex for managing input list of ensemble members */
std::mutex queue_lock;
std::mutex grid_lock;
/* We put this list at file scope so all the functions in this file
can use it without having to pass it around.   Hope that doesn't cause an issue*/
std::queue<int> member_to_process;
int fill_member_list(const ThreeComponentEnsemble& d)
{
  /* Better to be sure we don't have debris in this before starting */
  if(!member_to_process.empty())
      throw MsPASSError(string("fill_member_list:   ")
          + "member_to_process queue used to synchronize theads is not empty!\n"
          + "This should not happen and could cause a runaway - aborting",
        ErrorSeverity::Fatal);
  for(int i=0;i<d.member.size();++i)
  {
      if(d.member[i].live()) member_to_process.push(i);
  }
  return member_to_process.size();
}
/*! Multithreaded function to process one seismogram in an input ensemble d.
*/
void migrate_members_threaded(ThreeComponentEnsemble& d,
  GCLgrid& parent,
    GCLscalarfield3d& raygrid,
      GCLscalarfield3d& TPgrid,
         GCLscalarfield3d& Us3d,
           VelocityModel_1d& Vp1d,
             VelocityModel_1d& Vs1d,
               Metadata& control,
                 PWMIGfielddata& pwdgrid)
{
  int m;
  while(!member_to_process.empty())
  {
    queue_lock.lock();
    m = member_to_process.front();
    member_to_process.pop();
    queue_lock.unlock();
    PWMIGmigrated_seismogram dout;
    dout = migrate_one_seismogram(d.member[m], parent, raygrid, TPgrid,Us3d,
                      Vp1d, Vs1d, control);
    grid_lock.lock();
    pwdgrid.accumulate(dout);
    grid_lock.unlock();
  }
}

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
  nmembers = fill_member_list(d);
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
  int number_threads;
  number_threads = control.get_int("number_of_threads_per_worker");
  std::vector<std::thread> thread_pool;
  for(unsigned i=0;i<number_threads;++i)
  {
    /* The use of ref on all these arguments is an obscure issue related
    to template handling with std::thread.  The reference I used is
    here:  https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwiyiqPnzIz6AhUiATQIHWqLADsQFnoECBgQAQ&url=https%3A%2F%2Fcs.millersville.edu%2F~wkillian%2Farchive%2F2020%2Fspring%2Ffiles%2Fcsci476%2F05-C%2B%2B-Threads.pptx&usg=AOvVaw0oi1zLohTSsrci5GeowRyb

    This also has a lot thrown together putting the thread constructor
    inside the call to push_back, but using a temporary that would be
    clearer to read does no compiles - some obscure copy issue I think,
    */
    thread_pool.push_back(std::thread(migrate_members_threaded,
         ref(d), ref(parent), ref(*raygrid), ref(TPgrid), ref(Us3d),
           ref(Vp1d), ref(Vs1d), ref(control), ref(pwdgrid)));
  }
  for(unsigned i=0;i<number_threads;++i)
  {
    thread_pool[i].join();
  }
  delete raygrid;
  return pwdgrid;
}
}  //end namespace
