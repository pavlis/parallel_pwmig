#include <sstream>
#include <string>
#include "mspass/utility/MsPASSError.h"
#include "mspass/utility/AntelopePf.h"
#include "pwmig/pwmigcore/DepthDependentAperture.h"
namespace pwmig::pwmigcore
{
using namespace std;
using namespace pwmig::pwmigcore;
using mspass::utility::AntelopePf;

//
// Pf-based constrructor for DepthDependentAperture object.
//
//  Arguments:
//      md - PfStyleMetadata object created from a pf file.
//	tag - name identifying the Tbl defining the designed input
//
//  This routine throws an exception with a message string if the required
// information from the pf cannot be retrieved.
// This should probably normally be fatal.  A more complex return
// to the error was deemed unnesseary since this is a specialized
// object for pwstack
//
//  Author:  Gary L. Pavlis
//

DepthDependentAperture::DepthDependentAperture(const AntelopePf& md,
        const string tag)
{
    try {
        list<string> tlist=md.get_tbl(tag);
        size_t npts=tlist.size();
        t.reserve(npts);
        aperture.reserve(npts);
        cutoff.reserve(npts);
        list<string>::iterator tptr;
        int i;
        for(tptr=tlist.begin(),i=0;tptr!=tlist.end();++tptr,++i)
        {
            stringstream instr(tptr->c_str());
            double val;
            instr >>val; t.push_back(val);
            instr >>val; aperture.push_back(val);
            instr >>val; cutoff.push_back(val);
        }
    }catch(...){throw;};
}
/* Fresnel zone constructor.  See include file for parameter definitions. */
DepthDependentAperture::DepthDependentAperture(const double vs,
        const double vp,
            const double period,
                const double dtau,
                    const int ntau,
                        const double cutoff_multiplier,
                            const bool echoresult)
{
    size_t npts;
    npts = ntau;
    t.reserve(npts);
    aperture.reserve(npts);
    cutoff.reserve(npts);
    for(auto i=0;i<npts;++i) t.push_back(dtau*static_cast<double>(i));
    double lagtots=vp/(vp-vs);  // multiplier from S-P time to S time
    double ts;
    if(echoresult) cout << "Fresnel zone aperture calculation used"<<endl
                        << "lag   aperture(km)    cutoff(km)"<<endl;
    for(auto i=0;i<npts;++i)
    {
        ts=lagtots*t[i];
        double term1=ts+period/2.0;
        //aperture[i]=vs*sqrt(term1*term1 - ts*ts);
        //cutoff[i]=cutoff_multiplier*aperture[i];
        double anow, cutnow;
        anow=vs*sqrt(term1*term1 - ts*ts);
        cutnow=cutoff_multiplier*anow;
        aperture.push_back(anow);
        cutoff.push_back(cutnow);
        if(echoresult)cout<<t[i]<<"  "<<aperture[i]<<"  "<<cutoff[i]<<endl;
    }
}

DepthDependentAperture::DepthDependentAperture(const DepthDependentAperture& parent)
         : t(parent.t),aperture(parent.aperture),cutoff(parent.cutoff)
{
}


double linear_interpolate_irregular(const int npts,
  const vector<double>& x, const vector<double>& y,const double xsearch)
{
    int i1,i2;
    double yout;
    //silently return endpoints if outside the range
    if (xsearch<x[0]) return(y[0]);
    if (xsearch>x[npts-1]) return(y[npts-1]);
    for(i2=1;i2<npts;++i2)
        if(x[i2]>xsearch) break;
    i1 = i2-1;
    yout = y[i1] + (xsearch-x[i1])*(y[i2]-y[i1])/(x[i2]-x[i1]);
    return(yout);
}


double DepthDependentAperture::get_aperture(const double t0) const
{
  /* this avoids a const problem */
  size_t npts=this->npoints();
  return(linear_interpolate_irregular(npts,t,aperture,t0));
}


double DepthDependentAperture::get_cutoff(const double t0) const
{
  /* this avoids a const problem */
  size_t npts=this->npoints();
  return(linear_interpolate_irregular(npts,t,cutoff,t0));
}
DepthDependentAperture& DepthDependentAperture::operator=(const DepthDependentAperture& parent)
{
    if(this!=&parent)
    {
        t=parent.t;
        aperture=parent.aperture;
        cutoff=parent.cutoff;
    }
    return(*this);
}
double DepthDependentAperture::maximum_cutoff() const
{
  double maxval(0.0);
  size_t npts=cutoff.size();
  for(size_t i=0;i<npts;++i)
  {
    if(cutoff[i]>maxval) maxval=cutoff[i];
  }
  return maxval;
}
double DepthDependentAperture::maximum_aperture() const
{
  double maxval(0.0);
  size_t npts=aperture.size();
  for(size_t i=0;i<npts;++i)
  {
    if(aperture[i]>maxval) maxval=aperture[i];
  }
  return maxval;
}
} //end namespace
