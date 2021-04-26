#include <math.h>
#include <cstdio>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "perf.h"
#include "coords.h"
#include "pwstack.h"
#include "seispp.h"

#include "coords.h"
#include "stock.h"

#define RADIUS_EARTH 6378.164
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
void geographic_to_dne(double lat0, double lon0,
double lat, double lon, double *dnorth, double *deast)
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
int compute_pseudostation_weights(int nsta, double *dnorth, double *deast,
double aperture, double cutoff, double *w)
{
    int i;
    double distance, azimuth, denom;
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
void compute_pwmoveout(int nsta,double *deast, double *dnorth,
double ux, double uy, double *moveout)
{
    int i;
    for(i=0;i<nsta;++i)
    {
        moveout[i] = ux*deast[i] + uy*dnorth[i];
    }
}
void dzero(int n, double *d, int inc)
{
    int i, ii;
    for(i=0,ii=0;i<n;++i,ii+=inc) d[ii]=0.0;
}


void vscal(int n, double *w, int incw, double *y, int incy)
{
    int i, iw, iy;
    for(i=0,iw=0,iy=0;i<n;++i,iw+=incw,iy+=incy)
        y[iy] = w[iw]*y[iy];
}


void vadd(int n, double *x,int  incx, double *y,int  incy)
{
    int i, ix, iy;
    for(i=0,ix=0,iy=0;i<n;++i,ix+=incx,iy+=incy)
        y[iy] = x[ix] + y[iy];
}


string virtual_station_name(int ix1, int ix2)
{
    char name[8];
    sprintf(name,"%3.3d%3.3d",ix1,ix2);
    return(string(name));
}


// make station name for pwstack using 3 integer ids
string MakeDfileName(int evid, int x1, int x2)
{
    string dfile;
    ostringstream sbuf(dfile);
    sbuf<<"pwstack_"<<evid<<"_"<<x1<<"_"<<x2;
    return(sbuf.str());
}

/* Applies a simple ramp taper to any weights component that starts
at zero and then becomes nonzero.  This happens only when the 
pseudostation aperture is depth variable.  */
void taper_weights_array(dmatrix& weights, vector <bool>& use_this_sta,
        double dt, double aperture_taper_length)
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

double coherence(double ssqr,double ssqd)
{
    // This is an ambiguous choice to avoid a NaN but
    // I make it 0 as a safer choice.  Avoids 0/0
    if(ssqr<=0.0 || ssqd <= 0.0)
        return(0.0);
    else if(ssqr>ssqd)
        return(0.0);
    else
        return(1.0-sqrt(ssqr)/sqrt(ssqd));
}


/*  Computes stack coherence given an input ensemble sent as
a set of 3 matrices (one for each component) and parallel
weight matrix.  This procedure is linked to pwstack making
assumptions that make it far from general.  Key assumptions
are:  (1)  d and wgt have same t0 as stack
(2)  top mute has times consistent with above.

Uses a moving window computation of coherence over a
specified time gate of length ncoh in stack time sampling
units.
Arguments
d - 3 components stored in matrices with constant start time.
wgt - weight array to form stack.
stack - stack used to compute coherence.
dtcoh - sample interval for output coherence array.
cohwinlen - length (in time units) of coherence averaging window.
mute - top mute that is assumed to have been applied to d.
Coherence is forced to 0 in the zero area of the mute zone.
*/
Coharray compute_stack_coherence(vector<dmatrix> d,
dmatrix& wgt,
ThreeComponentSeismogram& stack,
double dtcoh,
double cohwinlen,
TopMute& mute)
{
    int i,j,k;

    // IMPORTANT assume size of d components and stack are always the same
    // Caller must be sure this is true. This was done for
    // efficiency assuming this code will not be used outside
    // of the pwstack program
    int ns=stack.ns;
    int nd=d[0].columns();
    // This is the number of coherence points for output
    int ncoh=static_cast<int>((stack.endtime()-stack.t0)/dtcoh);
    ++ncoh;                                       // standard interval versus points issue
    Coharray result(stack.t0,dtcoh,ncoh,cohwinlen);
    // Note these arrays are in ThreeComponentSeismogram
    // convension with components in rows of the matrix
    dmatrix r(3,ns);                              // accumulated weighted squared residual array
    // holds accumulated weighted sum of squares for norm d
    dmatrix dwtssq(3,ns);
    dmatrix work(3,ns);                           // simplifies algorithm below
    r.zero();   work.zero();  dwtssq.zero();
    double dtover2=dtcoh/2.0;
    double wtd;
    // compute residual matrix
    for(k=0;k<nd;++k)
    {
        for(i=0;i<3;++i)
        {
            // d is a matrix with seismograms in columns
            double *ptr=d[i].get_address(0,k);
            dcopy(ns,ptr,1,work.get_address(i,0),3);
        }
        for(i=0;i<3;++i)
        {
            for(j=0;j<ns;++j)
            {
                wtd=work(i,j)*wgt(j,k);
                dwtssq(i,j)=dwtssq(i,j)+wtd*wtd;
            }
        }
        work=work-stack.u;
        for(i=0;i<3;++i)
        {
            for(j=0;j<ns;++j)
            {
                wtd=work(i,j)*wgt(j,k);
                r(i,j)=r(i,j)+wtd*wtd;
            }
        }
    }
    // Some useful variables best computed only once for efficiency
    double winhalfwidth=cohwinlen/2.0;
    int nwin=static_cast<int>(cohwinlen/stack.dt) + 1;
    int count;
    double rcompsum[3],dcompsum[3],rsum,dsum;
    // Now we need to fill the output arrays.  This is a simple moving
    // average algorithm.  The oddities here are: (1) the top mute zone
    // is zeroed, and (2) the endtime limit has to deal with the overlap
    // variable.  the time gate is shortened when necessary on both ends
    double t;
    int is;
    for(j=0,t=stack.t0;j<ncoh;++j,t=t+dtcoh)
    {
        TimeWindow gate(t-winhalfwidth,t+winhalfwidth);
        if(t<mute.t0e)
        {
            // Don't need to alter compcoh because creation
            // initializes it to 0s
            result.coh.push_back(0.0);
        }
        else
        {
            for(i=0;i<3;++i)
            {
                rcompsum[i]=0.0;
                dcompsum[i]=0.0;
            }
            rsum=0.0;
            dsum=0.0;
            for(count=0,is=stack.sample_number(t);
                count<nwin && is<ns;
                ++is)
            {
                if(is>=0)
                {
                    ++count;
                    for(i=0;i<3;++i)
                    {
                        rcompsum[i]+=r(i,is);
                        dcompsum[i]+=dwtssq(i,is);
                        rsum+=r(i,is);
                        dsum+=dwtssq(i,is);
                    }
                }
            }
            if(count>0)
            {
                for(i=0;i<3;++i)
                    result.compcoh(i,j)=coherence(rcompsum[i],dcompsum[i]);
                result.coh.push_back(coherence(rsum,dsum));
            }
            else
            {
                for(i=0;i<3;++i)result.compcoh(i,j)=0.0;
                result.coh.push_back(0.0);
            }
        }
    }
    return(result);
}
Coharray::Coharray(const Coharray& parent)
{
	t0=parent.t0;
	dt=parent.dt;
	ns=parent.ns;
	winlen=parent.winlen;
	compcoh=parent.compcoh;
	coh=parent.coh;
}
Coharray& Coharray::operator=(const Coharray& parent)
{
	if(this!=&parent)
	{
		t0=parent.t0;
		dt=parent.dt;
		ns=parent.ns;
		winlen=parent.winlen;
		compcoh=parent.compcoh;
		coh=parent.coh;
	}
	return(*this);
}

/*!  Save coherence data as a set of seismic attribute traces
linked by attribute link table.

Arguments -
    coh - result of coherence calculator used in this program.
    stack - stack from which coh was derived.  It is assumed
    pfhcoh3c - special file handle for pwmig to hold 3c data
    pfhcoh - special file handle for pwmig for scalar coh data
*/

void save_coh(Coharray& coh,
  ThreeComponentSeismogram& stack,
    PwmigFileHandle& pfhcoh3c,
      PwmigFileHandle& pfhcoh)
{
    try
    {
        // First clone the stack.  We'll then modify it's contents to match the
        // coherence data and then save the results using dbsave for ThreeComponentSeismogram
        // and TimeSeries data
        ThreeComponentSeismogram coh3ctrace(stack);
	/* We have to put these into the trace objects as the
	ExtractComponent procedure requires them */
	coh3ctrace.put("samprate",1.0/coh.dt);
	coh3ctrace.put("time",coh.t0);
	coh3ctrace.put("nsamp",coh.ns);
        auto_ptr<TimeSeries>cohtrace(ExtractComponent(coh3ctrace,0));
        coh3ctrace.ns=coh.ns;
        cohtrace->ns=coh.ns;
        coh3ctrace.t0=coh.t0;
        cohtrace->t0=coh.t0;
        coh3ctrace.dt=coh.dt;
        cohtrace->dt=coh.dt;
        cohtrace->s=coh.coh;
        coh3ctrace.u=coh.compcoh;
	pfhcoh3c.save(coh3ctrace);
	pfhcoh.save(*cohtrace);;
    }
    catch (...)
    {
        throw;
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
the enscapsulates all this requires)
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
*/
int pwstack_ensemble(ThreeComponentEnsemble& indata,
RectangularSlownessGrid& ugrid,
TopMute& mute,
TopMute& stackmute,
int stack_count_cutoff,
double tstart,
double tend,
DepthDependentAperture& aperture,
double aperture_taper_length,
double centroid_cutoff,
double dtcoh,
double cohwinlen,
MetadataList& mdlcopy,
PwmigFileHandle& dfh,
PwmigFileHandle& coh3cfh,
PwmigFileHandle& cohfh)
{
    // lat0 and lon0 are location of target pseudostation grid point
    // elev0 is elevation of datum to use for geometric statics
    //
    double lat0,lon0,elev0;
    // Incident wavefield slowness vector components
    double ux0,uy0;
    string gridname;
    int ix1, ix2;
    int evid;
    string dfile;
    const double WEIGHT_MINIMUM=1.0e-2;
    try
    {
        lat0=indata.get_double("lat0");
        lon0=indata.get_double("lon0");
        elev0=indata.get_double("elev0");
        ux0=indata.get_double("ux0");
        uy0=indata.get_double("uy0");
        ix1=indata.get_int("ix1");
        ix2=indata.get_int("ix2");
        evid=indata.get_int("evid");
        gridname=indata.get_string("gridname");
    } catch (SeisppError& serr)
    {
        // The above are all set by main so this is the
        // correct error message.  Could perhaps be dropped
        // after debugging finished.
        serr.log_error();
        cerr << "Coding error requiring a bug fix"<<endl;
        exit(-1);
    }
    // This routine builds a sta name from the index positions.
    // Made a function to allow consistency between programs.
    //
    string sta=virtual_station_name(ix1,ix2);
    int i,j,k;
    vector<ThreeComponentSeismogram>::iterator iv,ov;
    int nsta = indata.member.size();
    // This computes the output gather size.  It assumes all data have
    // a common sample rate and we can just grab the first one in the
    // list.  It also quietly assumes a relative time base
    // so all times are computed relative to the start of
    // each trace.  Caller should guarantee this.
    // Current versioh assures this with clean_gather procedure
    //
    double dt=indata.member[0].dt;
    int nsin = indata.member[0].ns;

    ThreeComponentSeismogram *stackout;

    int istart = SEISPP::nint(tstart/dt);
    int iend = SEISPP::nint(tend/dt);
    int nsout = iend-istart+1;
    int ismute = SEISPP::nint(mute.t1/dt);
    int ismute_this, ie_this;

    /* Apply front end mutes to all traces */
    ApplyTopMute(indata,mute);

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
        int ierr;
        try
        {
            lat = (*iv).get_double("site.lat");
            lon = (*iv).get_double("site.lon");
            // assume metadata store these in degrees so we have to convert
            lat = rad(lat);
            lon = rad(lon);
            geographic_to_dne(lat0, lon0, lat, lon, &(dnorth[i]),&(deast[i]));
            elev[i]=(*iv).get_double("site.elev");
        } catch (MetadataError& merr)
        {
            throw merr;
        }
    }
    //
    // We want to make sure that all stations are in standard coordinates
    //
    for(iv=indata.member.begin();iv!=indata.member.end();++iv)
    {
        (*iv).rotate_to_standard();
    }
    //
    // the weights become are a nsta by nsamp matrix to allow
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
        int nused;

	// This procedure sets weights to zero outside the cutoff and we depend 
	// on this here to build the full weights matrix
        nused=compute_pseudostation_weights(nsta, &(dnorth[0]),&(deast[0]),
            aperture.get_aperture(tstart+dt*(double)i),
            aperture.get_cutoff(tstart+dt*(double)i),&(work[0]));
        for(j=0;j<nsta;++j)
        {
            if(work[j]>WEIGHT_MINIMUM)
            {
                weights(j,i)=work[j];
                if(!use_this_sta[j])use_this_sta[j]=true;
            }
            else
                weights(j,i)=0.0;
        }
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
    if(stack_count<stack_count_cutoff) return(-1);
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
    if(hypot(dncen,decen)>centroid_cutoff) return(-2);
    if(SEISPP_verbose) 
        cout << "Processing data for node ("<<ix1
            <<", "<<ix2<<") with fold="<<stack_count<<endl;
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
    int stack_start;
    for(i=0;i<nsout;++i) 
    {
	if(stack_weight[i]>WEIGHT_MINIMUM)
        {
		stack_start=i;
		break;
	}
    }
    /* Return -3 to signal to not use this pseudostation if the 
    first sample that satisfied the above test is after the 
    stack mute zone.  If running in verbose mode issue a diagnostic
    as a user can get into trouble with variable apertures 
    and this necessary detail. */
    if((stack_start*dt)>stackmute.t1)
	return(-3);
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
    // These are buffers used to store results for each pseudostation point
    // to allow i/o to be concentrated.  This was causing a performance problem
    // this is aimd to address.
    //
    list<ThreeComponentSeismogram> stacklist;
    list <Coharray> coharraylist;
    //
    // Loop over slowness grid range storing results in new output ensemble
    //

    int iu,ju,gridid;
    double ux,uy,dux,duy;
    /* Note older db version started gridid at 1 instead of 0 like this is */
    for(iu=0,gridid=0;iu<ugrid.nux;++iu)
    {
        for(ju=0;ju<ugrid.nuy;++ju,++gridid)
        {
            int iend_this;
            int ismin;
            char buffer[128];

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
            for(i=0,icol=0,iv=indata.member.begin(),ismute_this=0,iend_this=0,ismin=nsout;
                iv!=indata.member.end();++iv,++i)
            {
                int is0,ietest, is, ns_to_copy;
                int j0;
                double lag;
                nsin=iv->ns;
                //
                // Completely drop data for stations not marked with tiny or zero weight
                //
                if(use_this_sta[i])
                {

                    lag = tstart - (iv->t0) + moveout[i];
                    is0=SEISPP::nint(lag/dt);
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
            stackout = new ThreeComponentSeismogram(nsout);
            stackout->dt=dt;
            stackout->t0=tstart;
            stackout->live=true;
            stackout->u=stack;
            copy_selected_metadata(dynamic_cast<Metadata&>(indata),
                dynamic_cast<Metadata&>(*stackout),mdlcopy);
            stackout->put("ix1",ix1);
            stackout->put("ix2",ix2);
            stackout->put("ux",ux);
            stackout->put("uy",uy);
            stackout->put("gridid",gridid);
            stackout->put("dux",dux);
            stackout->put("duy",duy);
            stackout->put("ux0",ux0);
            stackout->put("uy0",uy0);
            // may want to output a static here, but it is probably better to
            // just keep a good estimate of elevation and deal with this in the
            // migration algorithm.
            stackout->put("elev",avg_elev);
            ApplyTopMute(*stackout,stackmute);
            // new March 2007: compute stack coherence
            Coharray coh=compute_stack_coherence(gather,gathwgt,*stackout,
                dtcoh,cohwinlen,stackmute);
            stacklist.push_back(*stackout);
            delete stackout;
            coharraylist.push_back(coh);
        }
    }
    list<ThreeComponentSeismogram>::iterator soutptr;
    list<Coharray>::iterator cohptr;
    /* Loop over the results of above.  Note we assume here that
    coharraylist and stacklist are the same length.*/
    for(soutptr=stacklist.begin(),cohptr=coharraylist.begin();
        soutptr!=stacklist.end();++soutptr,++cohptr)
    {
        try
        {
	    dfh.save(*soutptr);
	    save_coh(*cohptr,*soutptr,coh3cfh,cohfh);
        }
        catch(SeisppError& err)
        {
            err.log_error();
            cerr << "Write failure abort:  cannot continue"<<endl;
            exit(-1);
        }
    }
    return(stack_count);
}
