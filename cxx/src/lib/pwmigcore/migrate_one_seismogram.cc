//DEBUG - remove this and all cout statements after resolution
#include <iostream>
using namespace std;

#include "mspass/seismic/Seismogram.h"
#include "mspass/utility/dmatrix.h"
#include "pwmig/gclgrid/gclgrid.h"
#include "pwmig/seispp/VelocityModel_1d.h"
#include "pwmig/pwmigcore/PWMIGmigrated_seismogram.h"
#include "pwmig/pwmigcore/pwmig.h"
#include "pwmig/seispp/interpolator1d.h"
#include "pwmig/pwmigcore/Ray_Transformation_Operator.h"
namespace pwmig::pwmigcore {
using namespace mspass::utility;
using namespace mspass::seismic;
using namespace pwmig::pwmigcore;
using namespace pwmig::gclgrid;
using namespace pwmig::seispp;
using INTERPOLATOR1D::linear_vector_regular_to_irregular;

/* The GCL objects can't be const because the the lookup method alters
the object. */
PWMIGmigrated_seismogram migrate_one_seismogram(Seismogram& pwdata,
    GCLgrid& parent,
      GCLscalarfield3d& raygrid,
        GCLscalarfield3d& TPgrid,
           GCLscalarfield3d& Us3d,
             VelocityModel_1d& Vp1d,
               VelocityModel_1d& Vs1d,
                 Metadata& control)
{
  cout << "Entered function migrate_one_seismogram"<<endl;
  int k;  // primary depth index in loops.  Defined outside loops for error logging
  int kk;   // used repeatedly as an index in z running the reverse of index k
  string algname("migrate_one_seismogram");  // used in elog messages as algorithm name
  /* This gets a bit messy to mesh with special class we return that isn't
  quite a Seismogram. */
  if(pwdata.dead())
  {
    cout << "This datum was marked dead - returning deadguy"<<endl;
    PWMIGmigrated_seismogram deadguy;
    deadguy.elog=pwdata.elog;
    if(pwdata.is_defined("ix1"))
    {
      deadguy.ix1=pwdata.get_int("ix1");
    }
    if(pwdata.is_defined("ix2"))
    {
      deadguy.ix2=pwdata.get_int("ix2");
    }
    deadguy.elog=pwdata.elog;
    /* Make sure he is dead*/
    deadguy.live=false;
    return deadguy;
  }
  /* These two metadata gets are closely linked to a map function used
  in the python functoin that calls this one.   The names here must match
  that functions usage or we will surely abort.*/
  double ux,uy;
  ux=pwdata.get_double("ux0");
  uy=pwdata.get_double("uy0");
  SlownessVector u0(ux,uy);
  /* The contents of this function in the original pwmig C++ program were
  the innermost of 3 loops used to compute the GRT migration.   The following
  set of parameters are here extracted from the Metadata container passed
  through the arg with the symbol control.   The contents are not tested for
  validity here for efficiency.  Callers should guarantee that the control
  container has data linked to the keys used here.  That means the input
  should be validated with a python verify function that is planned for
  inclusion in this package.

  I expect to also many of the parameters here also fetching in the top level
  control function for the python script function that drives pwmig and replaces
  the c++ main.   That should further guarantee the following set of gets
  to not throw exceptions */
  bool use_3d_vmodel=control.get_bool("use_3d_velocity_model");
  bool use_grt_weights=control.get_bool("use_grt_weights");
  bool stack_only=control.get_bool("stack_only");
  int border_pad = control.get_int("border_padding");
	double zpad = control.get_double("depth_padding_multiplier");
  double taper_length=control.get_double("taper_length_turning_rays");
  bool rcomp_wt=control.get_bool("recompute_weight_functions");
	int nwtsmooth=control.get_int("weighting_function_smoother_length");
	bool smooth_wt;
	if(nwtsmooth<=0)
		smooth_wt=false;
	else
		smooth_wt=true;
  double dux=control.get_double("slowness_grid_deltau");
  double duy=dux;   // For this implementation slowness increments are the same in x and y
  double dz=control.get_double("ray_trace_depth_increment");
  /* End get calls on control container */

  bool weight_functions_set=false;
  int n3=raygrid.n3;
  double zmaxray;
  vector<double> Stime(n3), SPtime(n3);
  double t0,dt;
  t0=pwdata.t0();
  dt=pwdata.dt();
  /* first decide in which grid cell to place this seismogram
  WARNING:  don't use i, and j as loop indices as is common in the rest of
  this function */
  int i,j,gridid;
  i = pwdata.get_int("ix1");
  j = pwdata.get_int("ix2");
  gridid=pwdata.get_int("gridid");
  cout << "ix1="<<i<<" ix2="<<j<<" gridid="<<gridid<<" n3="<<n3<<endl;
  /* Original pwmig had a static option here.  As a module in mspass it
  would be far better to apply statics as a mspass operator */
  if( (i<0) || (i>raygrid.n1) || (j<0) || (j>raygrid.n2) )
  {
    /* We have this condition throw an exception because it only happens if
    something is corrupted in the data or there is a fatal mismatch of parameters.
    The planned parameter checker should test for this condition*/
    stringstream ss;
    ss << "migrate_one_seismogram:  grid index of data = ("
      << i << ","<< j
      << "is outside range of grid." <<endl
      << "Check grid definitions in database and pf"<<endl;
    throw MsPASSError(ss.str(),ErrorSeverity::Fatal);
  }
  /* We create the output now so we can post error messages to the
  elog of this object that is what we return*/
  PWMIGmigrated_seismogram result(i,j,n3);
  // We need a reference ray for the incident wavefield.
  // This ray is chopped up and transformed to GCLgrid coordinates
  // in the integration loop below.  We compute it here to
  // avoid constantly recomputing it
  zmaxray = raygrid.depth(i,j,0);
  zmaxray *= zpad;                      // Small upward adjustment to avoid rubber banding
  // problems in 3d media

  // Now compute a series of vectorish quantities that
  // are parallel with the ray path associated with the
  // current seismogram (i,j).  First we compute the S wave
  // travel time down this ray path (i.e. ordered from surface down)
  Stime = compute_Stime(Us3d,i,j,raygrid,use_3d_vmodel);
  // Now we compute the gradient in the S ray travel time
  // for each point on the ray.  Could have been done with
  // less memory use in the loop below, but the simplification
  // it provides seems useful to me.  Note this is a 3xn3 matrix
  // like gradTp that is filled inside the loop
  dmatrix gradTs=compute_gradS(raygrid,i,j,Vs1d);
  dmatrix gradTp(3,n3);
  dmatrix nup(3,n3);
  vector<double> zP(n3);
  // Now loop over each point on this ray computing the
  // p wave path and using it to compute the P time and
  // fill the vector gradP matrix
  // Note the double index, one running forward (k) and
  // a second (kk) running backward.  This complication is
  // caused by raygrid having upward oriented curves while
  // the S ray path (pathptr) is oriented downward.
  // Only need to compute P time to surface once for
  // each ray so we compute it before the loop below
  double Tpr=TPgrid.val[i+border_pad][j+border_pad][TPgrid.n3-1];
  double tlag,Tpx,tdelta;
  bool needs_padding;
  int padmark=raygrid.n3-1;
  bool tcompute_problem;                // set true if tt calculation fails away from edges
  /* This is the position of the incident ray at the base of the model.
   * Common for all in next loop and needed for travel time calculation.
   * Note the ugly border pad construct necessary to account for padding
   * efficiently - avoids a lookup */
  Geographic_point rTP_gp=fetch_TP_point_at_base(TPgrid,i+border_pad,j+border_pad);
  TPgrid.reset_index();
  SPtime.clear();
  tcompute_problem=false;
  needs_padding=false;
  /* This is how the slowness vector was set in original pwmig.  Here we pass
  u0 as an arg */
  //SlownessVector u0=svm0(i,j);
  /* This modification was required for parallel processing in revision
  Sept 2021.   The reason is that the lookup method was superceded by
  a parallel_lookup that required the caller to manage the current
  position index of the direction set iteration.   The old library
  maintained that only with internal, private attributes.  That was found to
  NOT be thread safe.  So, we use parallel_lookup.  Note in this loop
  that that happens twice - once in the compute_unit_normal_P function
  and in the direct call to parallel_lookup. */
  std::vector<int> lookup_origin;
	lookup_origin = TPgrid.get_lookup_origin();
	int ix1_0, ix2_0, ix3_0;
	ix1_0 = lookup_origin[0];
	ix2_0 = lookup_origin[1];
	ix3_0 = lookup_origin[2];
  cout << "loopup_origin result:  "<<ix1_0<<", "<<ix2_0<<", "<<ix3_0<<endl;
  for(k=0,kk=raygrid.n3-1;k<raygrid.n3;++k,--kk)
  {
    vector<double>nu;
    double vp;
    /* This is needed below to compute p*delta term */
    Geographic_point x_gp,rxTP_gp0;
    x_gp=raygrid.geo_coordinates(i,j,kk);

    nu = compute_unit_normal_P(TPgrid,raygrid.x1[i][j][kk],
      raygrid.x2[i][j][kk], raygrid.x3[i][j][kk],ix1_0,ix2_0,ix3_0);
    for(auto l=0;l<3;++l) nup(l,kk)=nu[l];
    // zP and gradTp are computed in reverse order
    zP[k]=raygrid.depth(i,j,kk);
    vp=Vp1d.getv(zP[k]);
    for(auto l=0;l<3;++l) gradTp(l,k)=nu[l]/vp;
    /* This section used to be a procedure.  Inlined for
    speed during a cleanup June 2012 */
    int error_lookup;
    cout << "origin passed to parallel_lookup:  "<<ix1_0<<", "<<ix2_0<<", "<<ix3_0<<endl;
    error_lookup=TPgrid.parallel_lookup(raygrid.x1[i][j][kk],
      raygrid.x2[i][j][kk],raygrid.x3[i][j][kk],ix1_0,ix2_0,ix3_0);
    switch(error_lookup)
    {
      case 0:
        cout << "origin passed to parallel_interpolate:  "<<ix1_0<<", "<<ix2_0<<", "<<ix3_0<<endl;
        Tpx=TPgrid.parallel_interpolate(raygrid.x1[i][j][kk],
          raygrid.x2[i][j][kk],raygrid.x3[i][j][kk],
          ix1_0,ix2_0,ix3_0);
        break;
      case 1:
        /* This means the point falls outside the P ray grid.
        Normally that means we should just break the loop.
        Issue a warning only if the kk index is small indicating
        the ray path is very short  */
        needs_padding=true;
        padmark=k;
        break;
      case -1:
      default:
        /* This is a convergence error in lookup and needs to be
        flagged as an error.  NOTE:  may want to consider trying again
        after resetting the origin */
        tcompute_problem=true;
        needs_padding=true;
        padmark=k;
    }
    /* New procedure (Mar 2016) to fix error in approximation found for
     * no antelope version pwmig2.0 */
    try
    {
      rxTP_gp0=get_gp_base_TPx(TPgrid,x_gp);
    }catch(MsPASSError& serr)
    {
      result.elog.log_error(serr);
      result.live=false;
      tcompute_problem=true;
    }
    if(tcompute_problem || needs_padding) break;
    /* Original version had a different form here with 3 terms.
     * We have to add a new one for relative time calculation
     * here to correct for p*delta */
    tdelta=compute_delta_p_term(rxTP_gp0,rTP_gp,u0);
    tlag=Tpx+Stime[k]+tdelta-Tpr;
    SPtime.push_back(tlag);
  }

  // skip this ray if there were travel time computation problems
  if(tcompute_problem)
  {
    cout << "Writing elog entry for tcompute_problem"<<endl;
    stringstream ss;
    ss << "Warning:  slowness gridid "<< gridid
      << ", grid position index ("
      << i << ","
      << j << ").  GCLscalarfield3d lookup failed at ray index="
      << k <<endl
      << "Geographic location: "<<raygrid.lat(i,j,kk)
      <<", "<<deg(raygrid.lon(i,j,kk))
      <<", "<<raygrid.depth(i,j,kk)<<endl
      << "Data from failed ray index downward will be zeroed"<<endl
      << "This may leave ugly holes in the output image."
      <<endl;
    result.elog.log_error(algname,ss.str(),ErrorSeverity::Complaint);
  }
  if (needs_padding)
  {
    cout << "Entered block for needs_padding"<<endl;
    if(padmark==0)
    {
      stringstream ss;
      ss << "Warning:   gridid "<< gridid
        << ", grid position index ("
        << i << ","
        << j << ")." <<endl
        << "P time interpolation failed at earth's surface.  "<<endl
        << "This should not happen and may cause holes in the output."<<endl
        << "This trace will be projected as zeros"<<endl;
      result.elog.log_error(algname,ss.str(),ErrorSeverity::Complaint);
    }
    else
    {
      /* Arbitrarily increment the times by dt intervals.  Will
      zero data beyond pad mark below.  Not the most efficient way
      to do this, but don't expect this to be that common of an issue*/
      double tpadding0=Stime[padmark-1];
      for(k=padmark;k<raygrid.n3;++k)
        Stime.push_back(tpadding0+dt*static_cast<double>(k-padmark+1));
    }
  }
  // We now interpolate the data with tlag values to map
  // the data from time to an absolute location in space
  // Note carefully SPtime and result.migrated_data sense are inverted from
  // the raygrid.  i.e. they are oriented for a ray from the surface
  // to the bottom.  Below we have to reverse this

  linear_vector_regular_to_irregular(t0,dt,pwdata.u,
    &(SPtime[0]),result.migrated_data);
  /* Cosine aper all data that was padded and zero the extension */
  if(needs_padding)
    cosine_taper_highend(result.migrated_data,padmark,taper_length);

  //
  // Compute the transformation vector for each
  // ray point in space.
  //
  cout << "Calling extract gridline with i,j,n3="<<i<<", "<<j<<", "<<n3<<endl;
  dmatrix *pathptr=extract_gridline(raygrid,i,j,raygrid.n3-1,3,true);
  /* Old version used ustack here.  This is more correct as that is an
  approximation.  that was an error in the old version that could create
  a small difference for a test data set comparison */
  Ray_Transformation_Operator trans_operator(parent,*pathptr,u0.azimuth(),nup);
  result.migrated_data=trans_operator.apply(result.migrated_data);
  // done with these now
  delete pathptr;

  // This computes domega for a ray path in constant dz mode
  // We would have to interpolate anyway to mesh with
  // time data choose the faster, coarser dz method here
  //
  if(rcomp_wt || !weight_functions_set)
  {
    /* This is a gross inefficiency in stack_only but since I only expect it to be
    used for CCP stacking equivalent, this should not be a big deal.  Probably should
    do it right some day */
    if(stack_only)
    {
      result.domega.clear();
      for(k=0;k<n3;++k) result.domega.push_back(1.0);
    }
    else
    {
      /* This is the normal block for computing solid angles */
      result.domega=compute_domega_for_path(u0,dux,duy,
        Vs1d, zmaxray,dz,
        raygrid, i, j, gradTp,zP);
    }
    if(use_grt_weights && (!stack_only) )
      result.dweight=compute_weight_for_path(gradTp,gradTs);
    else
      for(k=0;k<n3;++k)result.dweight[k]=1.0;
    if(smooth_wt && (!stack_only))
    {
      result.domega=running_average(result.domega,nwtsmooth);
      // Unnecessary when not using grt weighting
      if(use_grt_weights)
        result.dweight=running_average(result.dweight,nwtsmooth);
    }
    weight_functions_set=true;
  }
  result.live=true;
  cout << "Result object of migrate_one_seismogram"<<endl;
  cout << "ix1="<<result.ix1<<" ix2="<<result.ix2<<endl
       << "dweight and domega vector sizes:"<<result.dweight.size()<<", "<<result.domega.size()<<endl
       << "data matrix size="<<result.migrated_data.columns()<<endl;
  return result;
}

} // End namespace
