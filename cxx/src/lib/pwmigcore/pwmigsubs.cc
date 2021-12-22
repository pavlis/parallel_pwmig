#include <vector>
#include <float.h>
#include "misc/blas.h"
#include "pwmig/dsap/coords.h"
#include "pwmig/pwmigcore/pwmig.h"
#include "pwmig/seispp/VelocityModel_1d.h"
#include "pwmig/seispp/interpolator1d.h"
#include "pwmig/seispp/ray1d.h"
#include "pwmig/gclgrid/gclgrid.h"
#include "mspass/utility/dmatrix.h"
#include "mspass/utility/MsPASSError.h"
/* This fortran function is in libgclgrid, but we don't advertise it with an
entry in the gclgrid.h file.  WE thus need this prototype here. */
extern "C" {extern void treex3_(double *,int *,double *,int *, double *);}
namespace pwmig::pwmigcore
{
  using namespace std;
  using namespace pwmig::gclgrid;
  using namespace pwmig::seispp;
  using namespace pwmig::pwmigcore;
  using namespace mspass::seismic;
  using namespace mspass::utility;
  /* New procedure added 2015.  Taken from gclfield2vtk.cc
   * Removes mean for constant x3 slices.  Important for
   * proper display of tomography models showing absolute velocities */
  void remove_mean_x3(GCLscalarfield3d& f)
  {
    double mean;
    int i,j,k;
    double nval=static_cast<double>((f.n1)*(f.n2));
    for(k=0;k<f.n3;++k)
    {
      mean=0.0;
      for(i=0;i<f.n1;++i)
        for(j=0;j<f.n2;++j)
          mean += f.val[i][j][k];
      mean /= nval;
      cout << k << "      "<<mean<<endl;
      for(i=0;i<f.n1;++i)
        for(j=0;j<f.n2;++j)
          f.val[i][j][k] = f.val[i][j][k]-mean;
    }
  }


  /* New procedure added 2015.   Reverses order of travel times stored in
   * a raygrid.   Realized this simplified a lot of logic in the incident
   * ray travel time calculation.
   *
   * f is the travel time field that will be reversed.
   * current_order is a string defining the order f is in on entry.
   * if(current_order=='downward') assume the times are ordered from
   * the top of the grid (n3-1) to the bottom.  Otherwise use 'upward'.
   * Will abort if current_order is anything else.
   */

  void reverse_time_field(GCLscalarfield3d& f, string current_order)
  {
    int i,j,k;
    double deltat,tlast,dt_here;
    if(current_order=="downward")
    {
      for(i=0;i<f.n1;++i)
        for(j=0;j<f.n2;++j)
      {
        deltat=f.val[i][j][0]-f.val[i][j][1];
        tlast=0.0;
        for(k=1;k<f.n3-1;++k)
        {
          dt_here=f.val[i][j][k] - f.val[i][j][k+1];
          f.val[i][j][k]=tlast+deltat;
          deltat=dt_here;
          tlast=f.val[i][j][k];
        }
        /* First and last points are not altered in above loop*/
        f.val[i][j][f.n3-1]=tlast+deltat;
        f.val[i][j][0]=0.0;
      }
    }
    else if(current_order=="upward")
    {
      /* Warning - this block has never been tested */
      for(i=0;i<f.n1;++i)
        for(j=0;j<f.n2;++j)
      {
        deltat=f.val[i][j][f.n3-1]-f.val[i][j][f.n3-2];
        tlast=0.0;
        for(k=f.n3-2;k>0;--k)
        {
          dt_here=f.val[i][j][k-1] - f.val[i][j][k];
          f.val[i][j][k]=tlast+deltat;
          deltat=dt_here;
          tlast=f.val[i][j][k];
        }
        f.val[i][j][f.n3-1]=0.0;
        f.val[i][j][0]=tlast+deltat;
      }
    }
    else
    {
      stringstream ss;
      ss << "reverse_time_field:  coding error.   Illegal current_order parameter="
        <<current_order<<endl;
      throw MsPASSError(ss.str(),ErrorSeverity::Fatal);
    }
  }

  /* small companion to below blindly computes distance for ray path (stored in
   dmatrix ray) to current point (ray(0:2,i0)) from previous point (ray(0:2,i0-1)).
  This is done blindly with pointer arithmetic and assumes i0 is greater than or
  equal to 1.
  */
  double pathdist(dmatrix& ray,int i0)
  {
    double *ptr;                                    // Do this with pointer arithmetic for efficiency
    double dx,r;
    ptr=ray.get_address(0,i0);
    dx=(*ptr)-(*(ptr-3));
    r=dx*dx;
    ++ptr;
    dx=(*ptr)-(*(ptr-3));
    r+=dx*dx;
    ++ptr;
    dx=(*ptr)-(*(ptr-3));
    r+=dx*dx;
    return(sqrt(r));
  }

  /* The procedures below using pathintegral can get in trouble
  if the model is slightly below the reference ellipsoid.  This
  procedure aims to repair that potential thorny problem in a
  systematic way. It uses this algorithm:
  1) scan top surface of model and find the maximum point below
     the reference ellipsoid.
  2) if the point found in 1 is above the datum or 0, do nothing.
  3) otherwise fudge the model upward by the valued computed
    + the input parameter dz.

  Returns the amount the model was displaced.*/

  double fudge_model_upward(GCLscalarfield3d& model,double dz)
  {
    // A GCLgrid always has the top surface at n3-1
    // so we scan that surface for the lowest point relative
    // to r0_ellipse.
    int i,j;
    double lat,r,dr,drmin;
    int ksurface=model.n3 - 1;
    double offset;
    for(i=0;i<model.n1;++i)
    {
      for(j=0;j<model.n2;++j)
      {
        r=model.r(i,j,ksurface);
        lat=model.lat(i,j,ksurface);
        dr=r-r0_ellipse(lat);
        if((i==0) && (j==0))
        {
          drmin=dr;
        }
        else
        {
          drmin=min(dr,drmin);
        }
      }
    }
    if(drmin>=0.0)
      offset=0.0;
    else
    {
      offset=(-drmin)+dz;
      Geographic_point geo;
      Cartesian_point p;
      int k;
      for(i=0;i<model.n1;++i)
        for(j=0;j<model.n2;++j)
          for(k=0;k<model.n3;++k)
          {
            geo.lat=model.lat(i,j,k);
            geo.lon=model.lon(i,j,k);
            geo.r=model.r(i,j,k);
            geo.r+=offset;
            p=model.gtoc(geo);
            model.x1[i][j][k]=p.x1;
            model.x2[i][j][k]=p.x2;
            model.x3[i][j][k]=p.x3;
          }

    }
    return(offset);
  }

  /* computes travel times along a ray path defined by the input dmatrix of Cartesian
  coordinate (path) using a 3D model defined by the GCLsclarfield3d variable U3d.
  The path is assumed to be congruent coordinate system used in the 3d model grid.
  This is true in this code because of code in main that guarantees raygrids and
  the models are congruent.  DO NOT transport this code without dealing with that
  issue.

  The original pwmig code used a grid of total slowness.  This was revised in
  a 2015 revision to require the 3d model be specified as slowness perturbation from
  an (unspecified) background.   Regional tomography models have a finite extent
  so it is highly likely a ray will run outside of the model's bounding box.
  This is complicated by the fact a ray be truncated on either (including possibly
  both) ends.   The previous pwmig version only handled one end member of this family -
  rays that start at the surface but pass outside the model.  This version uses
  a very different algorithm.   Instead of using the GCLgrid pathintegral procedure
  this procedure basically does it's own path integral assuming it should zero
  points that fail the grid lookup procedure.  We assume that indicates a ray passing
  outside the tomography model.   The returned vector is then the accumulated
  time along the specified path.  For efficiency if no points in the path intersect
  the model volume the returned travel time vector will be zero length.  Otherwise it
  is guaranteed to be the same length as the number of points in path.

  Arguments:
  U3d - slowness PERTURBATION field.  It is assume the field stores slowness
  perturbation values with units of s/km.
  path - 3xN matrix of Cartesian points defining the path for integration.

  Returns:
  vector of integral of travel time anomalies along the path.  Direction
  of integration is defined by the path argument with integration done by
  simple summation.   This tacitly assumes that path is more densely
  sampled than U3d, which today is always true of tomography models for
  use with pwmig (i.e. paths must be much more closely sampled than velocities
  for this application).

  */
  vector<double> compute_3Dmodel_time(GCLscalarfield3d& U3d, dmatrix& path)
  {
    double tsum;                                    // holds accumulating sum of times
    int count;                                      // number of points actually accumulated in path integral
    int npts=path.columns();
    vector<double> times;
    times.reserve(npts);
    times.push_back(0.0);                           /* First point is always 0 */
    tsum=0.0;
    count=1;
    /* New code for parallel_lookup.  This is a bit awkward because the std::vector
	   container return doesn't match argument approach I used for holding the
	   last integer position managed by parallel_lookup. */
	  std::vector<int> lookup_origin;
	  lookup_origin = U3d.get_lookup_origin();
	  int ix1_0, ix2_0, ix3_0;
	  ix1_0 = lookup_origin[0];
	  ix2_0 = lookup_origin[1];
  	ix3_0 = lookup_origin[2];
    for(int i=1;i<npts;++i)
    {
      double dx1,dx2,dx3;
      double du;                                    // interpolated slowness perturbation
      /* Only add points when lookup succeeds - returns 0 */
      int iret=U3d.parallel_lookup(path(0,i),path(1,i),path(2,i),
                         ix1_0, ix2_0, ix3_0);
      if(iret==0)
      {
        du=U3d.parallel_interpolate(path(0,i),path(1,i),path(2,i),ix1_0,ix2_0,ix3_0);
        dx1=path(0,i)-path(0,i-1);
        dx2=path(1,i)-path(1,i-1);
        dx3=path(2,i)-path(2,i-1);
        tsum+=sqrt(dx1*dx1+dx2*dx2+dx3*dx3)*du;
        times.push_back(tsum);
        ++count;
        //fprintf(stdout,"i=%d time=%15.5g\n",i,tsum);
      }
    }
    /* clear the contents if there are zero hits.  Caller should test size of
    return for zero length */
    if(count<=1)
      times.clear();
    return times;
  }


  /* applies a cosine taper to the high end of a 3c matrix of data
  from the point mark-taper_length to mark.   This assumes the
  points in d are close to a regular sample interval, which is true
  in pwmig.   This procedure should not be used outside pwmig.
  */
  void cosine_taper_highend(dmatrix& d,int mark, int taper_length)
  {
    /* Shorten the taper length if mark is less than taper_length to
    simplify logic.   */
    if(taper_length>mark) taper_length=mark+1;
    double a=M_PI/static_cast<double>(taper_length-1);
    int i,ii,k;
    for(i=mark-taper_length,ii=0;ii<taper_length;++i,++ii)
    {
      double w;
      w=(cos(a*static_cast<double>(ii)) + 1.0)/2.0;
      for(k=0;k<3;++k) d(k,i)*=w;
    }
    /* Zero all beyond mark */
    for(i=mark;i<d.columns();++i)
      for(k=0;k<3;++k) d(k,i)=0.0;
  }


  /* Returns a 3x nx  matrix of ray path unit vectors pointing in
  direction of tangent to the path defined by the 3 x nx  input matrix
  passed as ray.  The algorithm is a simple forward difference scheme
  with no checking for possible roundoff problems.  Because we use forward
  differences the output vector at point i is computed from the segment
  immediately forward (x[i+1]-x[i]).  Note this means the sign convention
  is determined solely from the path variables.  Each vector is normalized
  to unity.  The last point two points in the output are equal to keep
  sizes consistent and is the only choice with a forward difference scheme.

  Note the dmatrix is created with new here and needs to be cleared by caller.

  Author:  Gary Pavlis
  Written:  August 2003
  */

  dmatrix *ray_path_tangent(dmatrix& ray)
  {
    int nx;
    double dx[3];
    double nrmdx;
    int i,j;

    nx=ray.columns();
    dmatrix *gptr;
    gptr = new dmatrix(3,nx);
    dmatrix& gamma=*gptr;

    //  This shouldn't happen but it is a disaster if it does so we retain
    // a trap that was in the original pwmig code
    if(nx<2 || ray.rows()!=3)
    {
      stringstream ss;
      ss << "ray_path_tangent:  input path size "
        << ray.rows()
        <<" by "
        << ray.columns() <<endl;
      throw MsPASSError(ss.str(),ErrorSeverity::Fatal);
    }

    for(i=0;i<nx-1;++i)
    {
      dx[0]=ray(0,i+1)-ray(0,i);
      dx[1]=ray(1,i+1)-ray(1,i);
      dx[2]=ray(2,i+1)-ray(2,i);
      nrmdx=dnrm2(3,dx,1);
      // Another trap like the one above - shouldn't happen but changed
      // from older exit to throw with fatal condition
      if(nrmdx<=0.0)
      {
        stringstream ss;
        ss<<"ray_path_tangent: coding error two points on path are equal."<<endl;
        throw MsPASSError(ss.str(),ErrorSeverity::Fatal);
      }
      dscal(3,1.0/nrmdx,dx,1);
      for(j=0;j<3;++j) gamma(j,i)=dx[j];
    }
    // copy the last point
    for(j=0;j<3;++j) gamma(j,nx-1)=gamma(j,nx-2);
    return(gptr);
  }


  /* Computes a dmatrix of gradient S vectors for a ray GCLgrid defining the x3
  grid line from x1,x2 grid position ix1,ix2 (convention in this program).
  Velocities used to scale the vectors are derived from the 1D velocity model
  passed as Vs1d.  This is an approximation but one that should make vitually no
  difference in the results for the simplicity and speed it (likely) buys.
  */

  dmatrix compute_gradS(GCLgrid3d& raygrid,int ix1, int ix2, VelocityModel_1d& Vs1d)
  {
    int k,kk,i;
    dmatrix gradS(3,raygrid.n3);
    double vs;
    double nrm;
    // We start at the surface and work down
    // gradS vectors should point upward
    // The raygrid lines are oriented from the bottom up so sign is
    // as below
    for(k=0,kk=raygrid.n3-1;k<raygrid.n3-1;++k,--kk)
    {
      gradS(0,k)=raygrid.x1[ix1][ix2][kk]
        - raygrid.x1[ix1][ix2][kk-1];
      gradS(1,k)=raygrid.x2[ix1][ix2][kk]
        - raygrid.x2[ix1][ix2][kk-1];
      gradS(2,k)=raygrid.x3[ix1][ix2][kk]
        - raygrid.x3[ix1][ix2][kk-1];
      // use the velocity from the upper point
      vs=Vs1d.getv(raygrid.depth(ix1,ix2,kk));
      // gradient is unit vector multiplied by slowness
      nrm=dnrm2(3,gradS.get_address(0,k),1);
      for(i=0;i<3;++i) gradS(i,k)/=(nrm*vs);
    }
    // Final point is just a copy of the second to last as we run out
    // of points in forward difference
    for(i=0;i<3;++i)gradS(i,raygrid.n3-1)=gradS(i,raygrid.n3-2);
    return(gradS);
  }


  vector<double> compute_unit_normal_P(GCLgrid3d& raygrid,
    double x1, double x2, double x3,int& ix1_0, int& ix2_0, int& ix3_0)
  {
    // Ray path tangents do not vary rapidly with position.  We just get the nearest
    // neighbor and compute by finite differences
    vector<double>nu;                               // unit vector here normally pointing generally up and in some direction
    int index[3];
    int err;
    double dx;
    double nrmx;
    double rtest;
    int k;
    Geographic_point geo;
    /* Note this function always alsters x1_0 etc.  and because we
    call it by reference the callers copy is altered */
    err=raygrid.parallel_lookup(x1,x2,x3,ix1_0,ix2_0,ix3_0);
    /* Note:  I think the previous version of this function had a bug in
    how it handled lookup error returns.   It wasn't disasterous it gave
    wrong but not drastically inaccurate results.   Important point for
    validation against previous results. Sept 2021 glp*/
    switch(err)
    {
      case 1:
      case -1:
        // If the lookup failed we start here.  Since we are using a grid
        // constructed of 1d rays we can cheat and jump to the center of the
        // grid and try to use something reasonably approximates the corect
        // ray direction
        raygrid.reset_index();
        raygrid.get_index(index);
        geo=raygrid.ctog(x1,x2,x3);
        rtest=geo.r;
        // hunting upward for the first point on the ray through the
        // origin that is just above the radius defined in geo
        // we set the index position there for the calculation below
        for(k=0;k<raygrid.n3;++k)
        {
          geo=raygrid.geo_coordinates(index[0],index[1],k);
          if(geo.r>rtest)
          {
            index[2]=k;
            break;
          }
        }
        if(k>=raygrid.n3) index[2]= raygrid.n3 - 2; // set to one below surface
        break;
        /* This was the old code.  code 2 was depricated and we don't want
        to set normal return (0) like this any more.  See replacement.
      case 0:
      case 2:
        raygrid.get_index(index);
        */
        case 0:
          index[0] = ix1_0;
          index[1] = ix2_0;
          index[2] = ix3_0;
          break;
        default:
          throw MsPASSError(string("compute_unit_normal_P: ")
             + " parallel_lookup function returned an unexpected error code.\n"
             + "Something is drastically wrong and the code needs a fix",
           ErrorSeverity::Fatal);
      };
      // Standard fix for forward diff.  Back off one if at the top edge
      if(index[2]>=raygrid.n3-1)--index[2];
      dx=raygrid.x1[index[0]][index[1]][index[2]+1];
      dx-=raygrid.x1[index[0]][index[1]][index[2]];
      nu.push_back(dx);
      dx=raygrid.x2[index[0]][index[1]][index[2]+1];
      dx-=raygrid.x2[index[0]][index[1]][index[2]];
      nu.push_back(dx);
      dx=raygrid.x3[index[0]][index[1]][index[2]+1];
      dx-=raygrid.x3[index[0]][index[1]][index[2]];
      nu.push_back(dx);
      nrmx=dnrm2(3,&(nu[0]),1);
      dscal(3,1.0/nrmx,&(nu[0]),1);
    return(nu);
  }


  /* This is a messy algorithm that extends the grid from each edge by
     pad.  There is probably a more compact way to do this, but I did
     it this way to make it more maintainable - i.e comprehendable.

   This routine is not intended to ever be a library routine, so it
   does not handle exceptions that could be thrown by the SlownessVectorMatrix
   class.   In fact, there should not be an exception here unless I made
   a coding error anyway.*/
  SlownessVectorMatrix pad_svm(const SlownessVectorMatrix& svm, const int pad)
  {
    int nr,nc,nrp,ncp;
    nr=svm.rows();
    nc=svm.columns();
    nrp=nr+2*pad;
    ncp=nc+2*pad;
    SlownessVectorMatrix paddedsvm(nrp,ncp);
    SlownessVector uij;
    int i,j,ii,jj;
    /* First fill the unpadded areas */
    for(i=0;i<nr;++i)
      for(j=0;j<nc;++j)
    {
      uij=svm(i,j);
      paddedsvm(i+pad,j+pad)=uij;
    }
    /* Now fill the corners with the vector from each corner */
    /* upper left*/
    uij=svm(0,0);
    for(i=0;i<pad;++i)
      for(j=0;j<pad;++j) paddedsvm(i,j)=uij;
    /* upper right */
    uij=svm(0,nc-1);
    for(i=0;i<pad;++i)
      for(j=pad+nc-1;j<ncp;++j) paddedsvm(i,j)=uij;
    /* lower left*/
    uij=svm(nr-1,0);
    for(i=pad+nr-1;i<nrp;++i)
      for(j=0;j<pad;++j) paddedsvm(i,j)=uij;
    /* lower right */
    uij=svm(nr-1,nc-1);
    for(i=pad+nr-1;i<nrp;++i)
      for(j=pad+nc-1;j<ncp;++j) paddedsvm(i,j)=uij;
    /* Finally do the borders.  Here we extend the value for the adjacent edge.*/
    /* left border */
    for(i=pad,ii=0;i<pad+nr;++i,++ii)
    {
      uij=svm(ii,0);
      for(j=0;j<pad;++j) paddedsvm(i,j)=uij;
    }
    /* Right border */
    for(i=pad,ii=0;ii<nr;++i,++ii)
    {
      uij=svm(ii,nc-1);
      for(j=0;j<pad;++j) paddedsvm(i,j+pad+nc)=uij;
    }
    /* Top border */
    for(j=pad,jj=0;jj<nc;++j,++jj)
    {
      uij=svm(0,jj);
      for(i=0;i<pad;++i)
        paddedsvm(i,j)=uij;
    }
    /* Bottom border */
    for(j=pad,jj=0;jj<nc;++j,++jj)
    {
      uij=svm(nr-1,jj);
      for(i=0;i<pad;++i) paddedsvm(i+pad+nr,j)=uij;
    }
    return(paddedsvm);
  }


  /* This small check routine is a sanity check.   There is a rare possibility
   * that i0 and j0 of the parent pseusostation grid is not consistent with
   * the actual grid contents.   This can only happen if there is an error
   * in the pf file or db that defines the attributes of the GCLgrid.
   * I put this in because the cost is tiny for the potential disaster it
   * could produce. */
  bool grid_mismatched(const GCLgrid& parent, const GCLgrid& padded, const int pad)
  {
    double lat1,lon1,lat2,lon2;
    /* We could test the origin, but a more conservative test
    is to compare the lat,lon of the 0,0 position of the parent*/
    lat1=parent.lat(0,0);
    lon1=parent.lon(0,0);
    lat2=padded.lat(pad,pad);
    lon2=padded.lon(pad,pad);
    const double frac(0.01);
    double test=parent.dx1_nom;
    test*=frac;
    double deltax=hypot(lat1-lat2,lon1-lon2);
    deltax *= EQUATORIAL_EARTH_RADIUS;

    if(deltax>test)
      return true;
    else
      return false;
  }


  /*! \brief Computes incident wave travel times in pwmig
  // Computes a GCLgrid defined field (raygrid) for an incident P wavefield.
  // The result is travel times from a given hypocenter to points inside a GCLgrid3d
  // object.  The times are computed from a 1d refernce model to the base of the model
  // and then using approximate ray tracing from the base up.
  // Actual algorithm computes times to the surface and then subtracts path integral times
  // to get 3d travel times at points in the grid.  The grid geometry is defined by
  // incident wave rays traced with a 1d reference model.  That is, they are like the
  // S wavefield raygrid used in the main program of pwmig.
  //
  // Modified July 2009:  added zdecfac parameter for efficiency.  Found that for higher
  //   sample rate data this created an unnecessarily huge travel time grid.
  //
  // \param pstagrid parent pseudostation grid used to build this grid.
  // \param border_pad is defines additional padding used to extend the pseudostation
  //     grid on which the data are defined.  This is necessary (essential really) because
  //     in pwmig this grid is used for all plane wave components while the S grids are
  //     different for each plane wave component.  This should be made large enough that
  //     it can contain all scattering points connected to S rays from the pseudostationh
  //     grid.
  // \param UP3d P wave slowness defined with a 3d field object
  // \param vp1d P wave 1d reference model.  Used to trace traces
  // \param h hypocenter of source that created wavefield being imaged
  // \param zmax ray trace parameters passed directly to BuildGCLraygrid
  // \param tmax ray trace parameters passed directly to BuildGCLraygrid
  // \param dt ray trace parameters passed directly to BuildGCLraygrid
  // \param zdecfac output ray grid is decimated by this amount in z
  //
  //\returns reference to object containing the travel time field
  //\exception catches all exceptions and passes them along.  Several functions
  //      called in this procedure could throw several possible exceptions I don't
  //      know well which is the reason for the catch all syntax.

  Major change January 2015:   Used to receive a Hypocenter and use
  an absolute travel time reference.  Hypocenter replaced with new concept
  of a SlownessVectorMatrix that defines incident wave direction at each
  pseudostation.   Further, we drop the absolute time reference and use
  the simpler algorithm of lags relative to the arrival time
  reference frame.   Point is this procedure now returns a grid with
  a completely different set of times.  The border padding creates a huge
  complication.   SlownessVectorMatrix is extended on edges by extending
  vectors on edges into pad region.   That is a much better approximation
  than previous version that used constant slowness over the entire grid.
  //
  */
  GCLscalarfield3d* ComputeIncidentWaveRaygrid(GCLgrid& pstagrid,
  int border_pad,
  GCLscalarfield3d& UP3d,
  VelocityModel_1d vp1d,
  SlownessVectorMatrix& svm,
  double zmax,
  double tmax,
  double dt,
  int zdecfac,
  bool use_3d)
  {
    int i,j,k;
    GCLscalarfield3d *Tp=NULL;   // Declared here to avoid memory leak if exceptions thrown
    // Incident wave slowness vector at grid origin
    // Always use 0 elevation datum for this calculation
    // note this isn't always used.
    //SlownessVector u0=h.pslow(pstagrid.lat0,pstagrid.lon0,0.0);
    /* Get u0 from center of the grid - not really all that
       necessary, but will retain it for compatibility */
    //SlownessVector u0=svm(svm.rows()/2,svm.columns()/2);
    // First set up the new expanded grid.  The new grid has border_pad cells
    // added on each side of the parent grid (pstagrid).
    int n1new, n2new;
    int i0new, j0new;
    n1new = pstagrid.n1 + 2*border_pad;
    n2new = pstagrid.n2 + 2*border_pad;
    i0new = pstagrid.i0 + border_pad;
    j0new = pstagrid.j0 + border_pad;
    try
    {
      // this assumes the pstagrid is a "regular" gclgrid created with this same constructor.
      // Otherwise there may be small disconnect, although it should not really matter
      // provided dx1_nome and dx2_nom approximate the parent well.
      GCLgrid ng2d(n1new, n2new, string("ng2d"),
        pstagrid.lat0,pstagrid.lon0,pstagrid.r0,
        pstagrid.azimuth_y, pstagrid.dx1_nom, pstagrid.dx2_nom, i0new, j0new);
      // We MUST make force ng2d to be in the same coordinate
      // system as the pseudostation grid
      remap_grid(ng2d,dynamic_cast<BasicGCLgrid&>(pstagrid));
      if(grid_mismatched(pstagrid,ng2d,border_pad))
        throw GCLgridError("padded GCLgrid geometry mismatch pseudostation.  i0 and/or j0 is probably wrong");
      /* We have to add padding to the the slowness vector matrix.
         This procedure does that.  On creation Tp contains
         1d model travel times*/
      SlownessVectorMatrix svmpadded=pad_svm(svm,border_pad);
      Tp=Build_GCLraygrid(ng2d,svmpadded,vp1d,zmax,tmax,dt);
      /* For the incident wave data we need to reverse the travel
       * time order from that computed by Build_GCLraygrid (top to
       * bottom is returned there).  This is because incident wave
       * is propagating up. */
      reverse_time_field(*Tp,string("downward"));
      /* Now apply a correction for 3D structure if requested.
       * This algorithm uses a path integral in compute_3Dmodel_time.
       * Loop is over each ray path that runs in k order from the
       * base of IncidentRaygrid to the surface. */
      if(use_3d)
      {
        /* This is the accumulated maximum travel time correction for the
        3D model.  This is reported to stdout */
        double dtmax(0.0);
        double dtmin(0.0);
        vector<double>dtP3d;
        for(i=0;i<Tp->n1;++i)
        {
          for(j=0;j<Tp->n2;++j)
          {
            /* This extracts a path from the bottom of the grid to the
             * surface to mesh with 3 coordinate ray paths */
            shared_ptr<dmatrix> path(extract_gridline(*Tp,i,j,0,3,false));
            dtP3d=compute_3Dmodel_time(UP3d,*path);
            /* the above procedure returns a vector of length 1
            if the ray has no intersection with UP3d */
            int dtrange=dtP3d.size();
            //cout << "length dtP3d="<<dtrange<<" ";
            if(dtrange>1)
            {
              /* this min test is not necessary but safe to avoid
              seg faults.  */
              for(k=0;k<min(dtrange,Tp->n3);++k)
              {
                Tp->val[i][j][k] += dtP3d[k];
                if(dtP3d[k]>dtmax) dtmax=dtP3d[k];
                if(dtP3d[k]<dtmin) dtmin=dtP3d[k];
              }
            }
          }
        }
        cout << "Incident Wave Ray Grid:  3d correction range="
          << dtmin<<" to "<<dtmax<<endl;
      }
      if(zdecfac>1)
      {
        GCLscalarfield3d *TPdecimated;
        TPdecimated=decimate(*Tp,1,1,zdecfac);
        delete Tp;
        return TPdecimated;
      }
      else
        return Tp;
    }
    catch (...)
    {
      if(Tp!=NULL) delete Tp;
      throw;
    }
  }


  /* This is a comparable function to Incident Wave calculator immediately
   * above but this compute a vector of S travel times from the surface to
   * positions on ray paths into the surface.   i.e. the travel time vector
   * returned is travel time from the surface to a position in raygrid.
   * Specifically result Stime[k] is the travel time from the surface postion
   * raygrid at [i][j][raygrid.n3-1] to image position
   * raygrid at [i][j][raygrid.n3-k-1].  Confusing indexing is the grid order
   * complication.
   *
   * if use_3d is false only 1d travel times will be returned. */
  vector<double> compute_Stime(GCLscalarfield3d& U3d,
  int i, int j, GCLscalarfield3d& raygrid,bool use_3d)
  {
    /* First we retrieve a vector of the 1D travel times in surface
       to bottom order */
    vector<double> Stime;
    Stime.reserve(raygrid.n3);
    int k,kk;
    for(k=0,kk=raygrid.n3-1;k<raygrid.n3;++k,--kk)
      Stime.push_back(raygrid.val[i][j][kk]);
    if(use_3d)
    {
      /* This retrieves a path from the surface downward from raygrid */
      shared_ptr<dmatrix> path(extract_gridline(raygrid,i,j,raygrid.n3-1,3,true));
      vector<double> dtS3d;
      dtS3d=compute_3Dmodel_time(U3d,*path);
      int dtrange=dtS3d.size();
      if(dtrange>0)
      {
        for(k=0;k<min(dtrange,raygrid.n3);++k)
          Stime[k] += dtS3d[k];
      }
    }
    return Stime;
  }


  /* This function computes the domega term for the inverse Radon transform inversion
  formula.  It works according to this basic algorithm.
  1.  If u0 is the slowness vector for the current ray path, compute four neighboring
  rays at 1/2 grid spacings:  u0+dux1/2, u0-dux1/2, u0+duy1/2, and u0-duy1/2 where the
  1/2 means the 1/2 a grid spacing.  In slowness space this forms a rectangle surrounding
  the u0 point.
  2.  Compute a ray for each of these points using constructors for a RayPathSphere
  object and the routine GCLgrid_Ray_Project.  Note these are computed with the equal
  depth option, not equal time.
  3.  Compute gradS for each point on each ray.
  4.  Interpolate the input gradP vectors to depths (the input is on an irregular
  grid).
  5.  Compute the gradS+gradP sums to build an array of unit vectors for each ray.
  6.  compute solid angles between + and - pairs of finite difference rays using
  a vector cross product and a small angle approximation (if a and b are unit vectors
  axb = |a||b|sin theta = sin theta ~ theta for small angles )  Solid angle is product
  of the two angles.

  Returns a vector of domega terms on the grid defined by the input ray path.

  Arguments:

  u0 - input slowness vector (note passed as SlownessVector object)
  u0 is base slowness for S raygrid passed as g (below)
  dux, duy - base slowness grid increments.
  vmod - 1d velocity model used for ray tracing
  zmax, dz - ray trace parameters.  trace to depth zmax with depth increment dz
  g - gclgrid for S ray geometry
  ix1, ix2  - index position of receiver point in gclgrid g.
  gradP - parallel matrix to path of P time gradient vectors.
  zP - vector of depths parallel to gradP
  (Note algorithm assumes gradP,zP are ordered from surface down)

  Returns an np length STL vector object.

  Author:  Gary Pavlis
  Written:  NOvember 2003
  */

  vector<double> compute_domega_for_path(SlownessVector& u0,double dux, double duy,
  VelocityModel_1d& vmod,
  double zmax,
  double dz,
  GCLgrid3d& g,int ix1, int ix2,
  dmatrix& gradP,
  vector<double>& zP)
  {
    SlownessVector udx1, udx2, udy1, udy2;
    dmatrix *pathdx1, *pathdx2, *pathdy1, *pathdy2;
    dmatrix *gradSdx1, *gradSdx2, *gradSdy1, *gradSdy2;
    int i;
    int npath;                                      // used to define length of valid output

    udx1 = u0;
    udx1.ux -= dux/2.0;
    udx2=u0;
    udx2.ux += dux/2.0;
    udy1 = u0;
    udy1.uy -= duy/2.0;
    udy2=u0;
    udy2.uy += duy/2.0;
    //
    // Now call the constructors for the basic ray paths for these four rays.
    // We set the constructor in depth model and ASSUME zmax is small enough that
    // the ray for each u will not turn within vmod.
    //
    RayPathSphere raydx1(vmod,udx1.mag(),zmax,1.0e99,dz,"z");
    RayPathSphere raydx2(vmod,udx2.mag(),zmax,1.0e99,dz,"z");
    RayPathSphere raydy1(vmod,udy1.mag(),zmax,1.0e99,dz,"z");
    RayPathSphere raydy2(vmod,udy2.mag(),zmax,1.0e99,dz,"z");
    //
    // project these into the GCLgrid coordinate system
    //
    try
    {
      pathdx1 = GCLgrid_Ray_project_down(g,raydx1, udx1.azimuth(),ix1,ix2,
        g.n3-1);
      pathdx2 = GCLgrid_Ray_project_down(g,raydx2, udx2.azimuth(),ix1,ix2,
        g.n3-1);
      pathdy1 = GCLgrid_Ray_project_down(g,raydy1, udy1.azimuth(),ix1,ix2,
        g.n3-1);
      pathdy2 = GCLgrid_Ray_project_down(g,raydy2, udy2.azimuth(),ix1,ix2,
        g.n3-1);
    }  catch (MsPASSError& err)
    {
      /* Rethrow this error to make the job abort.  Errors from these routines
      were originally fatal and should be for a reason */
      throw err;
    }
    //
    // get gradS values by first computing tangents as unit vector and then
    // scaling them by slowness at that depth. Note 1d slowness is used for
    // this calculation, not a 3d value.  This is an approximation that avoid
    // tremendous headaches that could occur otherwise.  It should not be a
    // significant error.
    //
    gradSdx1=ray_path_tangent(*pathdx1);
    gradSdx2=ray_path_tangent(*pathdx2);
    gradSdy1=ray_path_tangent(*pathdy1);
    gradSdy2=ray_path_tangent(*pathdy2);
    delete pathdx1;  delete pathdx2;  delete pathdy1;  delete pathdy2;
    //
    // Check all these matrices have the same size and truncate the
    // output if necessary writing a diagnostic in that situation.
    //
    npath = gradSdx1->columns();
    if(npath>gradSdx2->columns()) npath=gradSdx2->columns();
    if(npath>gradSdy1->columns()) npath=gradSdy1->columns();
    if(npath>gradSdy2->columns()) npath=gradSdy2->columns();
    /* This was in older pwmig code.  I'm commenting it out as in the
    mspass framework without reworking the return there is now way to post
    an error message.  I never saw this error in any of hundreds of runs so
    will assume it won't happen.  Keeping the old code here in case it does.
    if(npath!= gradSdx1->columns())
    {
      cerr << "compute_domega_for_path: irregular path length for bundle of 4 rays used to compute domega" << endl

        << "Slowness grid probably defines turning rays in thie model" << endl;
    }
    */
    for(i=0;i<npath;++i)
    {
      double slow;
      double depth;
      depth = raydx1.depth(i);
      slow = 1.0/vmod.getv(depth);
      dscal(3,slow,gradSdx1->get_address(0,i),1);
      dscal(3,slow,gradSdx2->get_address(0,i),1);
      dscal(3,slow,gradSdy1->get_address(0,i),1);
      dscal(3,slow,gradSdy2->get_address(0,i),1);
    }
    //
    // Now we interpolate the gradS vector arrays onto the gradP set
    // of depths. This produces a domega congruent with gradS and gradP
    // in the sense that the domega values will correspond to the same
    // node depths as those computed for gradS and gradP.  These can
    // be copied on return into the outputs, without any more interpolation
    //  albeit in the reverse order.
    //
    int npath_P;
    double z0=raydx1.depth(0);                      // assume all four rays have same z0
    npath_P = gradP.columns();
    if(npath_P!=zP.size())
    {
      delete gradSdx1;
      delete gradSdx2;
      delete gradSdy1;
      delete gradSdy2;
      stringstream ss;
      ss << "compute_domega_for_path:  vector grid point depths and gradP matrix do not match"
         << endl
         << "columns in gradP matrix="<<npath_P<<" but size of zP vector="<<zP.size()<<endl;
      throw MsPASSError(ss.str(),ErrorSeverity::Fatal);
    }
    dmatrix temp(3,npath_P);                        // Holds interpolated gradS values before replacement
    INTERPOLATOR1D::linear_vector_regular_to_irregular(z0,dz,*gradSdx1,&(zP[0]),temp);
    dmatrix nudx1(gradP+temp);                      // creates vector gradP + gradS
    delete gradSdx1;                                // Dont' need this any longer
    // Same for dx2, dy1, and dy2
    INTERPOLATOR1D::linear_vector_regular_to_irregular(z0,dz,*gradSdx2,&(zP[0]),temp);
    dmatrix nudx2(gradP+temp);
    delete gradSdx2;
    INTERPOLATOR1D::linear_vector_regular_to_irregular(z0,dz,*gradSdy1,&(zP[0]),temp);
    dmatrix nudy1(gradP+temp);
    delete gradSdy1;
    INTERPOLATOR1D::linear_vector_regular_to_irregular(z0,dz,*gradSdy2,&(zP[0]),temp);
    dmatrix nudy2(gradP+temp);
    delete gradSdy2;
    // normalize using antelope function d3norm to make these unit vectors
    for(i=0;i<npath_P;++i)
    {
      dr3norm(nudx1.get_address(0,i));
      dr3norm(nudx2.get_address(0,i));
      dr3norm(nudy1.get_address(0,i));
      dr3norm(nudy2.get_address(0,i));
    }
    //
    // Now get domega using cross products between pairs of unit vectors and
    // small angle approximation (sin theta approx theta when theta is small)
    //
    vector<double>domega;
    for(i=0;i<npath_P;++i)
    {
      double cross[3];
      double dtheta_x, dtheta_y,domega_i;
      dr3cros(nudx1.get_address(0,i),nudx2.get_address(0,i),cross);
      dtheta_x = dr3mag(cross);
      dr3cros(nudy1.get_address(0,i),nudy2.get_address(0,i),cross);
      dtheta_y = dr3mag(cross);
      // abs shouldn't be necessary really, but better
      // safe than sorry given 0 cost
      domega_i=fabs(dtheta_x*dtheta_y);
      domega.push_back(domega_i);
    }
    return(domega);
  }


  /* Much simpler routine to compute weight of this each member of a ray path
   * in generalized radon transform formula.  gradTp and gradTs are gradients
   * for P and S ray directions a defined in the Poppeliers and Pavlis (2003)
   * and Bostock et al. (2001).  Note this function drops the A term that
   * could be computed from geometric spreading because the assumption is
   * the this term is negligible for plane waves
   */
  vector<double> compute_weight_for_path(dmatrix& gradTp,dmatrix& gradTs)
  {
    const double eightpisq=78.95683521;
    double sum[3],nrmsum;
    int np=gradTp.columns();
    int i,j;
    vector<double> weight(np);

    for(j=0;j<np;++j)
    {
      for(i=0;i<3;++i)
      {
        sum[i]=gradTp(i,j)+gradTs(i,j);
      }
      nrmsum=dnrm2(3,sum,1);
      weight[j]=nrmsum*nrmsum*nrmsum/eightpisq;
    }
    return(weight);
  }


  /* Computing a running average of a specific length for a
  series of numbers stored in x.  npts is the number of points
  in the smoother.  npts/2 on the left and right are padded as the
  constants for the first and last npts that can be averaged.
  i.e. there is not endpoint tapering.  Return vector is x
  smoothed as requested.
  */
  vector<double> running_average(vector<double>& x, int ns)
  {
    int nx=x.size();
    int i,ii,j;
    double avg;
    vector<double> result(nx);

    if(nx<=ns)
    {
      for(i=0,avg=0.0;i<nx;++i)
        avg+=x[i];
      avg=avg/static_cast<double>(nx);
      for(i=0;i<nx;++i) result[i]=avg;
    }
    else
    {
      int npo2=(ns+1)/2;
      double avglen=static_cast<double>(ns);
      for(i=npo2-1;i<nx-(ns-npo2);++i)
      {
        for(j=i,ii=0,avg=0.0;ii<ns;++ii,++j)
          avg+=x[j-npo2+1];
        avg/=avglen;
        result[i]=avg;
      }
      // pad front and back
      for(i=npo2-2;i>=0;--i) result[i]=result[i+1];
      for(i=nx-(ns-npo2);i<nx;++i) result[i]=result[i-1];
    }
    return(result);
  }


  /* Overloaded version of same procedure as above for a matrix.  Smoothing
  here is done only along rows.  Used here for coherence grid to handle
  three-component data */
  dmatrix running_average(dmatrix x, int ns)
  {
    int nr=x.rows();
    int nc=x.columns();
    int i,j;
    vector<double> row;
    row.reserve(nc);
    dmatrix result(nr,nc);
    for(i=0;i<nr;++i)
    {
      for(j=0;j<nc;++j) row.push_back(x(i,j));
      row=running_average(row,ns);
      for(j=0;j<nc;++j) result(i,j)=row[j];
      row.clear();
    }
    return(result);
  }



  //
  // We store velocity models externally in the database, but internally
  // slowness is more appropriate.  This converts a velocity field to a
  // slowness field in an efficient way by altering the val array in place.
  //
  void VelocityFieldToSlowness(GCLscalarfield3d& g)
  {
    for(int i=0;i<g.n1;++i)
      for(int j=0;j<g.n2;++j)
        for(int k=0;k<g.n3;++k) g.val[i][j][k]=1.0/g.val[i][j][k];
  }


  /* Returns a 1d velocity model consistent with slowness field
  stored in u.  Algorithm averages node values at each constant x3
  level in u and sets 1d node values to be the same.  It then
  computes gradients to provide linear interpolator between each
  node point in the 1d model.  This makes the result as close to
  the 3d model as reasonably possible in terms of interpolation method. */
  VelocityModel_1d DeriveVM1Dfrom3D(GCLscalarfield3d& u)
  {
    VelocityModel_1d result(u.n3);
    double zavg,uavg;
    int i,j,k;
    int nperlayer=(u.n1)*(u.n2);
    /* We have to count backwards because the grid orientation is
    from the bottom to the top */
    for(k=u.n3-1;k>=0;k--)
    {
      zavg=0.0;
      uavg=0.0;
      for(i=0;i<u.n1;++i)
      {
        for(j=0;j<u.n2;++j)
        {
          zavg+=u.depth(i,j,k);
          uavg+=u.val[i][j][k];
        }
      }
      zavg /= nperlayer;
      uavg /= nperlayer;
      result.z.push_back(zavg);
      result.v.push_back(1.0/uavg);
    }
    /* Now set gradients */
    for(k=0;k<result.nlayers-1;++k)
    {
      double deriv;
      deriv=(result.v[k+1] - result.v[k])
        / (result.z[k+1]-result.z[k]);
      result.grad.push_back(deriv);
    }
    /* Set the gradient for the final point to 0 */
    result.grad.push_back(0.0);
    return(result);
  }


  SlownessVector slowness_average(LoggingEnsemble<Seismogram> *d)
  {
    double ux,uy;
    int n;
    n=d->member.size();
    if(n<=0 || (d->dead())) throw MsPASSError(string("slowness_average:  ")
        +"procedure was passed an empty ensemble ",ErrorSeverity::Fatal);
    ux=0.0;
    uy=0.0;
    for(int i=0;i<n;++i)
    {
      if(d->member[i].live())
      {
        ux+=d->member[i].get_double("ux");
        uy+=d->member[i].get_double("uy");
      }
    }
    SlownessVector result;
    result.ux=ux/static_cast<double>(n);
    result.uy=uy/static_cast<double>(n);
    ux /= static_cast<double>(n);
    uy /= static_cast<double>(n);
    if(hypot(ux,uy)<FLT_EPSILON)
    {
      try
      {
        /* assume n!= 0 or we don't get here */
        ux=d->member[0].get_double("ux0");
        uy=d->member[0].get_double("uy0");
        SlownessVector u0(ux,uy);
        if(u0.mag()<FLT_EPSILON)
          return(SlownessVector(0.0,0.0,0.0));
        else
          return(SlownessVector(0.0,0.0,u0.azimuth()));

      } catch (MsPASSError& mde)
      {
        /* This is bad form in a parallel environment but there is not
        good alternative. */
        cerr << "slowness_average (Warning): slowness_average. "
          << "ux0 and uy0 not defined. using 0"
          <<endl;
        cerr << "Message posted:"<<endl;
        mde.log_error();
        return(SlownessVector(0.0,0.0,0.0));
      }
    }
    else
    {
      return(SlownessVector(ux,uy));
    }
  }


  Geographic_point fetch_TP_point_at_base(GCLscalarfield3d& TP,int i, int j)
  {
    Geographic_point pt;
    pt=TP.geo_coordinates(i,j,0);
    return(pt);
  }


  Geographic_point get_gp_base_TPx(GCLscalarfield3d& TP,Geographic_point xgp)
  {
    /* These matrices hold the Jacobian and inverse Jacobian for cells */
    dmatrix J(3,3),Jinv(3,3);
    /* This is how we get an index to the grid - relic */
    int tpind[3];
    Cartesian_point cp=TP.gtoc(xgp);
    /* WARNING WARNING WARNING:   In current position of this procedure we can be sure
     * that the index pointer is already position and we do not need to waste time on
     * calling the lookup method.  This is a very dangerous assumption made for speed.*/
    TP.get_index(tpind);
    /* Compute unit vectors on 1 and 2 generalized coordinate directions
     * at the image point x.  Use the lookup point */
    int k;
    double x[3],x0[3],x1[3],x2[3],x3[3],dx1[3],dx2[3],dx3[3];
    dvector dx(3),dxunit(3);
    x[0]=cp.x1;  x[1]=cp.x2;  x[2]=cp.x3;
    x0[0]=TP.x1[tpind[0]][tpind[1]][tpind[2]];
    x0[1]=TP.x2[tpind[0]][tpind[1]][tpind[2]];
    x0[2]=TP.x3[tpind[0]][tpind[1]][tpind[2]];
    x1[0]=TP.x1[tpind[0]+1][tpind[1]][tpind[2]];
    x1[1]=TP.x2[tpind[0]+1][tpind[1]][tpind[2]];
    x1[2]=TP.x3[tpind[0]+1][tpind[1]][tpind[2]];
    x2[0]=TP.x1[tpind[0]][tpind[1]+1][tpind[2]];
    x2[1]=TP.x2[tpind[0]][tpind[1]+1][tpind[2]];
    x2[2]=TP.x3[tpind[0]][tpind[1]+1][tpind[2]];
    x3[0]=TP.x1[tpind[0]][tpind[1]][tpind[2]+1];
    x3[1]=TP.x2[tpind[0]][tpind[1]][tpind[2]+1];
    x3[2]=TP.x3[tpind[0]][tpind[1]][tpind[2]+1];
    for(k=0;k<3;++k) dx(k)=x[k]-x0[k];
    for(k=0;k<3;++k) dx1[k]=x1[k]-x0[k];
    for(k=0;k<3;++k) dx2[k]=x2[k]-x0[k];
    for(k=0;k<3;++k) dx3[k]=x3[k]-x0[k];
    /* Here we build the Jacobian and the x point */
    for(k=0;k<3;++k)
    {
      J(k,0)=dx1[k];
      J(k,1)=dx2[k];
      J(k,2)=dx3[k];
    }
    /* This is an analytic FORTRAN routine to compute the inverse of
     * a 3x3 matrix.   In libgclgrid */
    int three(3);
    double det;
    treex3_(J.get_address(0,0),&three,Jinv.get_address(0,0),&three,&det);
    dxunit=Jinv*dx;
    /* Now we do the same Jacobian calculation for the cell at the base of the the TP raygrid.
     * We compute the piecing for the ray linked to the scatter point x by using the Jacobian
     * at this point and a forward an inverse transformation */
    x0[0]=TP.x1[tpind[0]][tpind[1]][0];
    x0[1]=TP.x2[tpind[0]][tpind[1]][0];
    x0[2]=TP.x3[tpind[0]][tpind[1]][0];
    x1[0]=TP.x1[tpind[0]+1][tpind[1]][0];
    x1[1]=TP.x2[tpind[0]+1][tpind[1]][0];
    x1[2]=TP.x3[tpind[0]+1][tpind[1]][0];
    x2[0]=TP.x1[tpind[0]][tpind[1]+1][0];
    x2[1]=TP.x2[tpind[0]][tpind[1]+1][0];
    x2[2]=TP.x3[tpind[0]][tpind[1]+1][0];
    x3[0]=TP.x1[tpind[0]][tpind[1]][1];
    x3[1]=TP.x2[tpind[0]][tpind[1]][1];
    x3[2]=TP.x3[tpind[0]][tpind[1]][1];
    for(k=0;k<3;++k) dx1[k]=x1[k]-x0[k];
    for(k=0;k<3;++k) dx2[k]=x2[k]-x0[k];
    for(k=0;k<3;++k) dx3[k]=x3[k]-x0[k];
    for(k=0;k<3;++k)
    {
      J(k,0)=dx1[k];
      J(k,1)=dx2[k];
      J(k,2)=dx3[k];
    }
    /* We must zero the x3 component of x transformed to unit cell as we
     * want the point at the base */
    dxunit(2)=0.0;
    /* J is the transformation matrix to convert dxunit to physical distance in the
     * cell at the base.*/
    dx=J*dxunit;
    /* x0 is our reference point so we add dx to that.  Type colision of
     * dvector and C array - store in C array */
    for(k=0;k<3;++k) x0[k]+=dx(k);
    Geographic_point gpr0x=TP.ctog(x0[0],x0[1],x0[2]);
    return gpr0x;
  }


  /* New functions added April 2015 for version using relative times.   These
   * are procedures that deal with the p*delta term in the travel time equations.
   * finds an interpolated geographic point at the depth zx
   * along the incident ray defined by the path in TP with top at i,j. Note
   * the issue about border padding makes this not consistent with S ray raygrid
   * positions.*/
  Geographic_point find_TP_at_x_depth(GCLscalarfield3d& TP,int i, int j, double zx)
  {
    /* The grid is always oriented with n3-1 at the earth's surface.  We assume the
     * number of points along a ray here is not huge so we do a simple linear search
     * and use a linear interpolator for lat a lon */
    int k;
    k=TP.n3-1;
    while((zx >= TP.depth(i,j,k)) && (k>0)) --k;
    /* Silently return with the last point coordinate if we hit the bottom (this
     * perhaps should be treated as an exception */
    if(k==0)
      return(TP.geo_coordinates(i,j,0));
    double lat1,lat2,lon1,lon2,z1,z2;
    lat1=TP.lat(i,j,k-1);
    lon1=TP.lon(i,j,k-1);
    z1=TP.depth(i,j,k-1);
    lat2=TP.lat(i,j,k);
    lon2=TP.lon(i,j,k);
    z2=TP.depth(i,j,k);
    Geographic_point gpr0;
    gpr0.lat=INTERPOLATOR1D::linear_scalar(zx,z1,lat1,z2,lat2);
    gpr0.lon=INTERPOLATOR1D::linear_scalar(zx,z1,lon1,z2,lon2);
    gpr0.r=r0_ellipse(gpr0.lat)-zx;
    return gpr0;
  }


  double compute_delta_p_term(Geographic_point r0x, Geographic_point r0,
  SlownessVector u0)
  {
    /* dsap function to computer distance and azimuth.   gcp distance in degrees
     * is depth independent */
    double delta,az;
    dist(r0.lat,r0.lon,r0x.lat,r0x.lon,&delta,&az);
    /*az is the angle from N of gcp from r0 to r0x.  This defines offset component
     * for distance components as gcp angles */
    double delx=delta*sin(az);
    double dely=delta*cos(az);
    //DEBUG
    /*
    cout << "Radius of r0x and r0 points="<<r0x.r<<" "<<r0.r<<" diff="<<r0x.r-r0.r<<endl
        << "Delta (deg)="<<deg(delta)<<" azimuth(deg)="<<deg(az)<<endl
        << "Delx,dely(deg)="<<deg(delx)<<" "<<deg(dely)<<endl
        << "Slowness mag="<<u0.mag()<<" azimuth="<<deg(u0.azimuth())<<endl;
        */

    return(EQUATORIAL_EARTH_RADIUS*(delx*u0.ux + dely*u0.uy));
  }


  /* Extracts du for this enseble.  This procedure should never to extracted
  from this program without major modification.  REason is it make two
  extreme assumptions appropriate here, but anything but general.
  (1)  assume slowness vector attributes extracted from ensemble member
  objects exist and we don't need an error handler for metdata gets
  (2) all members have the same delta u so result can be extracted from
  member 0.
  A less obvious assumption is that the ensemble has data or the request for
  member 0 would throw an exception.
  */
  SlownessVector EnsembleDeltaSlow(ThreeComponentEnsemble *d)
  {
    double ux,uy;
    ux=d->member[0].get_double("ux");
    uy=d->member[0].get_double("uy");
    ux-=d->member[0].get_double("ux0");
    uy-=d->member[0].get_double("uy0");
    return(SlownessVector(ux,uy));
  }
} // End namespace
