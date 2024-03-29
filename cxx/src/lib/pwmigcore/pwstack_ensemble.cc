#include <math.h>
#include <cstdio>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>

#include "misc/blas.h"
#include "pwmig/pwmigcore/pwstack.h"
#include "mspass/seismic/keywords.h"
#include "mspass/utility/ProcessingHistory.h"
#include "mspass/utility/SphericalCoordinate.h"
#include "mspass/utility/dmatrix.h"


#define RADIUS_EARTH 6378.164

namespace pwmig{
namespace pwmigcore
{
using namespace std;
using namespace mspass::seismic;
using namespace pwmig::pwmigcore;
using namespace mspass::utility;
using mspass::algorithms::TopMute;
/*
using mspass::utility::Metadata;
using mspass::utility::MetadataList;
using mspass::utility::ProcessingHistory;
*/

/* For plane wave moveout computations a local cartesian coordinate
system is used wrt to a particular origin.  This approximation is
reasonable until the distance of a station from the origin
becomes a significant fraction of the epicentral distance.
In pseudostation stacking we can minimize the impact of this approximation
by continually translating the origin to the pseudostation point and
computing distances to all stations wrt to that point.  Since large
distance stations receive a low weight, time alignments at large
distances become unimportant.  this function is used to computer
vector distances in this context.

Arguments:
lat0, lon0 - origin to compute dnorth deast from
(IN RADIANS)
lat, lon - lat an lon of point to be converted to dnorth deast
dnorth, deast - output dnorth, deast vector

Author:  G Pavlis
Written:  June 2000
*/
void geographic_to_dne(const double lat0, const double lon0,
const double lat, const double lon,  double *dnorth, double *deast)
{
    double azimuth;
    double d;

    dist(lat0,lon0,lat,lon,&d,&azimuth);
    d *= RADIUS_EARTH;
    *deast = d*sin(azimuth);
    *dnorth = d*cos(azimuth);

}


/* This function computes pseudostation weights for a 2d Gaussian function
to implement pseudostation stacking as described in Neal and Pavlis (1999) grl.

Arguments:
    nsta - number of stations
    dnorth, deast - parallel vectors of length nsta of deast-dnorth coordinates
        relative to a reference geographical position.
    aperture - sigma parameter in Gaussian formula
    cutoff - stations outside this distance (km) are given 0 weight
    w - vector of length nsta (output assumed allocated by caller) to
        hold weights.

The function uses a dnorth-deast approximation.  For every data set I know
of this is not a problem for computing these weights or even for computing
moveout.  The reason is that coherence is lost before this approximation
becomes a problem.

Returns count of nonzero weights.

Author:  G Pavlis
Written:  June 2000
Modified:  Converted to C++ and assimilated into new a working implementation
in March 2003
*/
int compute_pseudostation_weights(const int nsta, const double *dnorth, const double *deast,
const double aperture, const double cutoff, double *w)
{
    int i;
    double distance, denom;
    int count=0;

    denom = 2.0*aperture*aperture;

    for(i=0;i<nsta;++i)
    {
        distance = hypot(dnorth[i],deast[i]);
        if(distance>cutoff)
            w[i] = 0.0;
        else
        {
            w[i] = exp(-distance*distance/denom);
            ++count;
        }
    }
    return(count);
}


/* Computes vector of nsta moveout corrections for slowness vector
 * (ux,uy) using local coordinates stored in parallel deast, dnorth
 * vectors (both assumed length nsta).
 * Assumes moveout has been allocated previusly.
 */
void compute_pwmoveout(const int nsta,const double *deast, const double *dnorth,
const double ux, const double uy, double *moveout)
{
    int i;
    for(i=0;i<nsta;++i)
    {
        moveout[i] = ux*deast[i] + uy*dnorth[i];
    }
}
void dzero(const int n, double *d, const int inc)
{
    int i, ii;
    for(i=0,ii=0;i<n;++i,ii+=inc) d[ii]=0.0;
}


void vscal(const int n, double *w, const int incw, double *y, const int incy)
{
    int i, iw, iy;
    for(i=0,iw=0,iy=0;i<n;++i,iw+=incw,iy+=incy)
        y[iy] = w[iw]*y[iy];
}


void vadd(const int n, double *x,const int  incx, double *y,const int  incy)
{
    int i, ix, iy;
    for(i=0,ix=0,iy=0;i<n;++i,ix+=incx,iy+=incy)
        y[iy] = x[ix] + y[iy];
}



/* Applies a simple ramp taper to any weights component that starts
at zero and then becomes nonzero.  This happens only when the
pseudostation aperture is depth variable.  */
void taper_weights_array(dmatrix& weights, const vector <bool>& use_this_sta,
        const double dt, const double aperture_taper_length)
{
	int nsta=weights.rows();
	int nsamp=weights.columns();
	int i,j,jstart;
	for(i=0;i<nsta;++i)
	{
		if(use_this_sta[i])
		{
			jstart=-1;
			for(j=0;j<nsamp;++j)
			{
				if(weights(i,j)>0.0)
				{
					jstart=j;
					break;
				}
			}
			if(jstart>0)
			{
				double dwpersample;
				dwpersample=dt/aperture_taper_length;
				double scale;
				for(j=jstart,scale=0.0;(j<nsamp) && (scale<1.0);
					++j,scale+=dwpersample)
				{
					weights(i,j)*=scale;
				}
			}
		}
	}
}


/*
Main processing function for this program.  Takes an input
data ensemble and produce a complete suite of plane-wave stacks
defined by the Rectangualr_Slowness_Grid object.
This function uses and enhancement from the original Neal and Pavlis
and Poppeliers and Pavlis papers.  It allows the aperture to
be time variable.  This is generalized, but the normal expectation
is that the apeture would grow wider with time to compensate
somewhat for diffraction.  The aperture widths should increase with lag
or ugly artifacts can appear when traces come in and out of the stack
at variable lag.  This is compensated if a trace appears moving forward
in time lag by a linear taper method.

Arguments
indata - input raw data ensemble (see object definition
the enscapsulates all this requires).  Note the contents are altered by
top mutes.  Also for efficiency the indata ensmble should only be data
within the range defined by the maximum cutoff distance.  Note the
input is passed by value because rotate_to_standard is called on all
members. Passing a reference could have bad side effects.
ugrid - object to define slowness grid for stacking
mute - mute applied to data before stacking
stackmute - mute applied to data after stacking
(This stack is aligned relative to latest
mute time of raw data in stack.)
lat0, lon0 - pseudostation grid point (IN RADIANS)
ux0, uy0 - slowness vector of input data (stacking
is relative to this vector)
tstart:tend - define time period for output stack
The routine pretty much assumes relative timing
so this is normally time wrt the start of
each trace in the ensemble.
aperture - defines variable aperture stack weighting
(see above)
dtcoh - sample interval for coherence calculation
Assumed implicitly to be larger than data dt.
overlap - overlap fraction.  coherence windows at dtcoh
intervals overlpa by this fraction to form an
implicit smoothing.
mdlcopy - defines metadata to copy from the raw data
ensemble to each stack output
mdlout - list of metadata to be saved to output database
(fixed for a given schema)
am - AttributeMap object defining mapping of internal to
external namespace (invariant between calls)
dir and dfile - define file name where output written (dir/dfile)
dbh - DatabaseHandle object for output

Normal return is stack count.  Returns a negative number for different errors.
A zero return is normal, but means there is no data within a the aperture of
the cutoff circle around this pseudostation point.

Changed July 1, 2008
Now returns count of fold for this grid point.  negative or 0
means did nothing.

Throws a MetadataError exception object if there are problems
parsing required metadata from any trace.  Current caller will
abort the program on this condition, but evolution might want
to produce a handler.
Change Jan 15,2004
Used to pass lat0, lon0 by arg list, now passed through the
ensemble metadata.
Change August,2009
Added two parameters to make variable depth aperture work correctly:
centroid_cutoff is used with a centroid calculation.  If the centroid
of the stations retained for a pseudostation fall outside this distance from
the pseudostation point, return -2 and do nothing.  The other is
aperture_taper_length.  A linear taper like a top mute of this length
(in seconds) is applied to weights whenever the first nonzero value is
not at zero lag.

spring 2021:  major changes to adapt the code to MsPASS.  Too many
changes to list
*/
LoggingEnsemble<Seismogram> pwstack_ensemble(LoggingEnsemble<Seismogram>& indata,
  const RectangularSlownessGrid& ugrid,
   TopMute& mute,
    TopMute& stackmute,
     int stack_count_cutoff,
      const double tstart,
       const double tend,
        const DepthDependentAperture& aperture,
         const double aperture_taper_length,
          const double centroid_cutoff,
            const bool save_history,
             const string algid)
{
    const string base_error("pwstack_ensemble:  ");
    const string algname("pwstack");
    // lat0 and lon0 are location of target pseudostation grid point
    //
    double lat0,lon0;
    // Incident wavefield slowness vector components
    double ux0,uy0;
    string gridname;
    int ix1, ix2;
    //int evid;
    string dfile;
    const double WEIGHT_MINIMUM=1.0e-2;
    /* Return immediately doing nothing if the input is already marked dead or
    the ensemble is emtpty.  We return what was received and post no error messages
    for either situation assuming whatever caused the condition was logged then.*/
    if(indata.dead()) return(indata);
    if(indata.member.size()==0) return(indata);

    /* this use of a copy constructor followed by clearing the member
    vector is not efficient, but we expect with the new api using MongoDB
    the input ensembles will not usually be that large. */
    LoggingEnsemble<Seismogram> ensout(indata);
    ensout.member.clear();
    try
    {
        lat0=indata.get_double("lat0");
        lon0=indata.get_double("lon0");
        /* for consistency with other coordinate data lat0 and lon0 
        are posted in degrees.  We need to convert them to radians 
        to be consistent with this function. */
        lat0 = rad(lat0);
        lon0 = rad(lon0);
        ux0=indata.get_double("ux0");
        uy0=indata.get_double("uy0");
        ix1=indata.get_int("ix1");
        ix2=indata.get_int("ix2");
        //evid=indata.get_int("evid");
        gridname=indata.get_string("gridname");
        /* It would be nice if we could handle MongoDB source_id attribute
        here, but unfortunately boost::any doesn't like ObjectIDs - I think.
        We have to do some things in the python wrappers anyway so we handle
        the source_id and related issues on the python side.  This is problematic
        but I see now alternative solution for now. */
    } catch (MsPASSError& serr)
    {
      stringstream ss;
      ss << base_error<<" Error extracting ensemble metadata."<<endl
         << "Error posted:  "<<serr.what()<<endl;
      ensout.elog.log_error(algname,ss.str(),ErrorSeverity::Invalid);
      ensout.kill();
      return ensout;
    }

    int i,j,k;
    vector<Seismogram>::iterator iv;
    int nsta = indata.member.size();
    // This computes the output gather size.  It assumes all data have
    // a common sample rate and we can just grab the first one in the
    // list.  It also quietly assumes a relative time base
    // so all times are computed relative to the start of
    // each trace.  Caller should guarantee this.
    // Current versioh assures this with clean_gather procedure
    //
    double dt=indata.member[0].dt();
    int nsin = indata.member[0].npts();

    Seismogram *stackout;

    int istart = (int)round(tstart/dt);
    int iend = (int)round(tend/dt);
    int nsout = iend-istart+1;

    /* Apply front end mutes to all traces */
    for(iv=indata.member.begin();iv!=indata.member.end();++iv) mute.apply(*iv);

    /* We need dnorth, deast vectors to compute moveout sensibly
    for this program.  Since we use them repeatedly we will
    extract them once from the gather.*/
    vector <double> dnorth;
    vector <double> deast;
    vector <double> elev;
    dnorth.resize(nsta);   deast.resize(nsta);  elev.resize(nsta);
    for(i=0,iv=indata.member.begin();
		iv!=indata.member.end();++iv,++i)
    {
        double lat,lon;
        try
        {
            lat = (*iv).get_double(SEISMICMD_rlat);
            lon = (*iv).get_double(SEISMICMD_rlon);
            // assume metadata store these in degrees so we have to convert
            lat = rad(lat);
            lon = rad(lon);
            geographic_to_dne(lat0, lon0, lat, lon, &(dnorth[i]),&(deast[i]));
            elev[i]=(*iv).get_double(SEISMICMD_relev);
        } catch (MsPASSError& merr)
        {
          stringstream ss;
          ss << base_error
            << "metadata error extracting receiver coordinates from data member="
            << i <<endl
            << "Message posted:"<<endl
            << merr.what()<<endl;
          ensout.elog.log_error(algname,ss.str(),ErrorSeverity::Complaint);
          /* We kill this data member only - later we test for count with
          any kills */
          iv->kill();
        }
    }
    //
    // We want to make sure that all stations are in standard coordinates
    //
    for(i=0,iv=indata.member.begin();iv!=indata.member.end();++i,++iv)
    {
      try{
        (*iv).rotate_to_standard();
      } catch(MsPASSError& merr) {
        stringstream ss;
        ss <<base_error
           << "Failure of rotate_to_standard with member "<<i<<endl
           << "Killing entire ensemble containing this problem datum - message posted follows"
           << merr.what()<<endl;
        ensout.elog.log_error(algname,ss.str(),ErrorSeverity::Complaint);
        /* Here we also kill the member with the problem */
        iv->kill();
      }
    }
    //
    // the weights become a nsta by nsamp matrix to allow
    // variable length apertures.  This algorithm assumes the
    // aperture object input is an nsta by
    // Should use a matrix class for this, but it is easy enough
    // for this simple case to just code inline.

    dmatrix weights(nsta,nsout);
    int stack_count=0;                            //maximum number of traces not zeroed
    vector <double> work(nsta);
    vector <bool> use_this_sta(nsta);
    for(i=0;i<nsta;++i) use_this_sta[i]=false;
    for(i=0;i<nsout;++i)
    {

	// This procedure sets weights to zero outside the cutoff and we depend
	// on this here to build the full weights matrix
  /* Note in mspass conversion the zero weight section should rarely if ever
  be hit.   Retained to maintain that flexibility in the event the geospatial
  query feature of MongoDB has issues I don't anticipate. May want to revise
  this section in the future, but it should have near zero impact on performance.*/
        compute_pseudostation_weights(nsta, &(dnorth[0]),&(deast[0]),
              aperture.get_aperture(tstart+dt*(double)i),
                  aperture.get_cutoff(tstart+dt*(double)i),&(work[0]));
        for(j=0;j<nsta;++j)
        {
            if(work[j]>WEIGHT_MINIMUM)
            {
                weights(j,i)=work[j];
                /* It is a bit inefficient setting this for all lags, but
                this is known to work*/
                if(!use_this_sta[j])use_this_sta[j]=true;
            }
            else
                weights(j,i)=0.0;
        }
    }
    /* Check for kills and force the use_this_sta for those to false.
    This is a retrofit for mating with mspass to match functional model */
    for(i=0,iv=indata.member.begin();iv!=indata.member.end();++i,++iv)
    {
      if(iv->dead()) use_this_sta[i]=false;
    }
    double dncen,decen;  // Centroid of live sta in x,y coord
    double avg_elev,sum_wgt;                      // computed now as weighted sum of station elevations
    dncen=0.0;  decen=0.0;
    for(i=0,stack_count=0,avg_elev=0.0,sum_wgt=0.0;
            i<nsta;++i)
    {
        double w;
        if(use_this_sta[i])
	      {
		      ++stack_count;
          w=weights(i,0);
		      dncen+=w*dnorth[i];
		      decen+=w*deast[i];
        	avg_elev += w*elev[i];
        	sum_wgt += w;
	      }
    }
    ///
    // ERROR RETURN THAT NEEDS TO BE HANDLED GRACEFULLY
    // I don't throw an exception because this should not be viewed as
    // an exception.  It is a case that has to be handled gracefully.
    // it will happen often at the edges of arrays
    //
    if(stack_count<stack_count_cutoff)
    {
      stringstream ss;
      ss << base_error << "stack count="<<stack_count
      <<" was below input threshold="<<stack_count_cutoff;
      ensout.elog.log_error(algname,ss.str(),ErrorSeverity::Invalid);
      ensout.kill();
      return ensout;
    }
    /* Failsafe to avoid possible divide by zero.  If the
   sum_wgt value is too small, just compute avg_elev as a mean.
   Odds are if this is the situation this pseudostation is likely
   to get tossed anyway, but this makes the code more robust.
   Perhaps should post a warning in verbose mode, but odds are small
   this will ever recover if this block is entered.*/
    if(sum_wgt<WEIGHT_MINIMUM)
    {
        dncen=0.0;
        decen=0.0;
        for(i=0,avg_elev=0.0,sum_wgt=0.0;i<nsta;++i)
        {
    	    if(use_this_sta[i])
    	    {
            	avg_elev += elev[i];
              dncen+=dnorth[i];
              decen+=deast[i];
            	sum_wgt += 1.0;
    	    }
	      }
   }
    dncen/=sum_wgt;
    decen/=sum_wgt;
    avg_elev/=sum_wgt;
    if(hypot(dncen,decen)>centroid_cutoff)
    {
      stringstream ss;
      ss << base_error<<"Centroid of pseudostation group has excessive skew and will be discarded"
        << "Computed distance from centroid to center="<<hypot(dncen,decen)
        << " exceed limit imposed ="<<centroid_cutoff;
      ensout.elog.log_error(algname,ss.str(),ErrorSeverity::Invalid);
      ensout.kill();
      return ensout;
    }
    /* New section to handle mspass processing history.   Put here to
    avoid the work for (frequent) situation where the ensemble is dropped. */
    vector<ProcessingHistory*> phvec;
    if(save_history)
    {
      for(iv=indata.member.begin(),i=0;iv!=indata.member.end();++iv,++i)
      {
        if(use_this_sta[i])
          phvec.push_back(dynamic_cast<ProcessingHistory*>(&(indata.member[i])));
      }
    }
    vector <double>moveout(nsta);
    dmatrix stack(3,nsout);
    vector<double>stack_weight(nsout);
    vector<double>twork(nsout);
    // New March 2007:  these matrices hold stack members
    // and associated weights for each sample.
    // They are used for coherence calculations
    vector<dmatrix> gather;
    dmatrix gathwgt(nsout,stack_count);
    gathwgt.zero();
    /* Here we compute a vector of sum of weights at each time step
    and copy only nonzero weight vectors to gathwgt vector.  We do
    this here because this vector and matrix  are common to all slowness vectors
    results we compute here. */
    int icol; // defined here and used in a similar context inside slowness grid loop.
    for(i=0,icol=0;i<nsta;++i)
    {
	    if(use_this_sta[i])
	    {
            vadd(nsout,weights.get_address(i,0),nsta,&(stack_weight[0]),1);
            dcopy(nsout,weights.get_address(i,0),nsta,
                            gathwgt.get_address(0,icol),1);
	          ++icol;
	    }
    }
   /* Find the first sample with a sum of wts above the threshold */
   /* This loop goes with the relic below - it also needs to be deleted
   when confirmed it is not necessary
    int stack_start;
    for(i=0;i<nsout;++i)
    {
	    if(stack_weight[i]>WEIGHT_MINIMUM)
      {
		    stack_start=i;
		    break;
	    }
    }
    */

    /* Return -3 to signal to not use this pseudostation if the
    first sample that satisfied the above test is after the
    stack mute zone.  If running in verbose mode issue a diagnostic
    as a user can get into trouble with variable apertures
    and this necessary detail. */
    /* New TopMute will kill any data with this condition.  For now
    I'll comment this section out */
    /*===============DELETE WHEN clear not needed=================
    if((stack_start*dt)>stackmute.t1)
    {
      stringstream ss;
      ss<<base_error<<"stack start time="<<stack_start*dt
        << " is greater than the specified end of the stack top mute zone="
        << stackmute.t1()<<endl
        <<"Null return as this is known to produce invalid data";
      ensout.elog.log_error(ss.str(),ErrorSeverity::Invalid);
      return ensout;
    }
	return(-3);
  ====================================*/
    /* One last step to clean up the weights array.  If the aperture
    width increases with depth we need to taper the leading edge to
    avoid transients.   This rather inelegant procedure does this.
    It will do nothing if the aperture is constant for all lags */
    taper_weights_array(weights,use_this_sta,dt,aperture_taper_length);

    /* Create the working gather build as 3 matrices, one for each channel,
    with stack_count columns per member.   More logical than a 3d array
    in this case in the authors view. */
    for(i=0;i<3;++i)
    {
        gather.push_back(dmatrix(nsout,stack_count));
    }

    //
    // Loop over slowness grid range storing results in new output ensemble
    //

    int iu,ju,gridid;
    double ux,uy,dux,duy;
    /* This is standard thing done for efficiency on vector containers. */
    ensout.member.reserve(ugrid.ux_size()*ugrid.uy_size());

    /* Note older db version started gridid at 1 instead of 0 like this is */
    for(iu=0,gridid=0;iu<ugrid.ux_size();++iu)
    {
        for(ju=0;ju<ugrid.uy_size();++ju,++gridid)
        {
            int iend_this;
            int ismin;

            stack.zero();
            for(i=0;i<3;++i)gather[i].zero();

            dzero(nsout,&(stack_weight[0]),1);

            /* The input gather is assumed prealigned with the slowness
            vector defined by ux0,uy0.  We use relative moveouts
            from this base moveout for the actual stacks dux,duy*/
            dux = ugrid.ux(iu);
            duy = ugrid.uy(ju);

            ux =  ux0+dux;
            uy =  uy0+duy;
            // moveout computed here and used below assumes the
            // data are aligned on the P arrival
            compute_pwmoveout(nsta,&(deast[0]),&(dnorth[0]),dux,duy,&(moveout[0]));
            for(i=0,icol=0,iv=indata.member.begin(),iend_this=0,ismin=nsout;
                iv!=indata.member.end();++iv,++i)
            {
                int is0,ietest, is, ns_to_copy;
                int j0;
                double lag;
                nsin=iv->npts();
                //
                // Completely drop data for stations not marked with tiny or zero weight
                //
                if(use_this_sta[i])
                {

                    lag = tstart - (iv->t0()) + moveout[i];
                    is0=(int)round(lag/dt);
                    if(is0>=0)
                    {
                        // This block is for positive moveout = negative shift
                        j0=0;
                        is=is0;
                        ismin=min(is,ismin);
                        ietest=is0+nsout-1;
                        if(ietest<=nsin)
                        {
                            ns_to_copy=nsout;
                            iend_this=nsout;
                        }
                        else
                        {
                            ns_to_copy=nsin-is0;
                            iend_this=max(ns_to_copy,iend_this);
                        }
                    }
                    else
                    {
                        // This block is for negative moveout=positive shift
                        j0=-is0;                  // j0 always positive here
                        is=0;
                        ismin=0;
                        ns_to_copy=nsout-j0;
                        if(ns_to_copy>nsin)
                        {
                            ns_to_copy=nsin-j0;
                            iend_this=max(nsin,iend_this);
                        }
                        else if(ns_to_copy<=0)
                        {
                            ns_to_copy=0;
                        }
                        else
                        {
                            iend_this=nsout;
                        }
                    }

                    if(ns_to_copy>0)
                    {
                        for(j=0;j<3;++j)
                        {
                            dzero(nsout,&(twork[0]),1);
                            // This is a slow way to copy these
                            // data to twork, but it is safe.
                            // An earlier version had a problem
                            // with stray indices.
                            int kk,jj;
                            for(k=0,kk=is,jj=j0;
                                (k<ns_to_copy)&&(kk<nsin)&&(jj<nsout);
                                ++k,++jj,++kk)
                            {
                                twork[jj]=iv->u(j,kk);
                                gather[j](jj,icol)=twork[jj];
                            }

                            vscal(nsout,weights.get_address(i,0),nsta,&(twork[0]),1);
                            vadd(nsout,&(twork[0]),1,stack.get_address(j,0),3);
                            //mp.load(stack,string("stack"));
			    /* We accumulate the stack_weight vector only on the first component. */
			    if(j==0) vadd(nsout,weights.get_address(i,0),nsta,&(stack_weight[0]),1);
                        }
                        ++icol;
                    }
                }
            }
            // normalize the stack.
            // not trivial for a variety of reasons
	    // This uses a threshold to avoid divide by zero but
	    // no other complexity.
	    for(i=0;i<nsout;++i)
            {
                if(stack_weight[i]>WEIGHT_MINIMUM)
                {
                    for(j=0;j<3;++j) stack(j,i)/=stack_weight[i];
                }
                else
                {
                    for(j=0;j<3;++j) stack(j,i)/=WEIGHT_MINIMUM;
                }
            }

            // Create the output stack as a 3c trace object and copy
            // metadata from the input into the output object.
            stackout = new Seismogram(nsout);
            /* copy the input ensemble's metadata to stackout.
            We add to it later, but we use this mechanism to pass
            global metadata to from the input ensemble to the stacks
            Note in older versions this was done with a selective copy.
            (copy_selected_metadata function).
            Note this is an obscure trick with inheritance in C++
            */
            stackout->Metadata::operator=(indata);
            /*
            stackout->dt=dt;
            stackout->t0=tstart;
            stackout->live=true;
            stackout->u=stack;
            */
            stackout->set_dt(dt);
            stackout->set_t0(tstart);
            stackout->set_live();
            stackout->set_npts(nsout);
            dcopy(3*nsout,stack.get_address(0,0),1,stackout->u.get_address(0,0),1);
            stackout->set_tref(TimeReferenceType::Relative);
            /*
            copy_selected_metadata(dynamic_cast<Metadata&>(indata),
                dynamic_cast<Metadata&>(*stackout),mdlcopy);
                */
            stackout->put("ix1",ix1);
            stackout->put("ix2",ix2);
            stackout->put("ux",ux);
            stackout->put("uy",uy);
            stackout->put("gridid",gridid);
            stackout->put("dux",dux);
            stackout->put("duy",duy);
            stackout->put("ux0",ux0);
            stackout->put("uy0",uy0);
            /* may want to output a static here, but it is probably better to
            just keep a good estimate of elevation and deal with this in the
            migration algorithm.  Note this is a name change for the key.  In
            the original pwmig package this was just elev.  Changed here to
            mesh with naming in the python wrapper for this function.
            That is fragile to beware. */
            stackout->put("pseudostation_elev",avg_elev);
            stackmute.apply(*stackout);
            if(save_history)
            {
              stackout->new_ensemble_process(algname,algid,
                AtomicType::SEISMOGRAM,phvec,true);
            }
            ensout.member.push_back(*stackout);
            delete stackout;
        }
    }
    /* Essential to do this for new implementation with a live attribute
    for the ensmble.  Could have been done earlier but at the exit of this
    loop is the more obvious place. */
    ensout.set_live();
    return(ensout);
}

} // End namespace pwstack_
} //End namespace pwmig
