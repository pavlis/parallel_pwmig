#ifndef _PWMIG_H_
#define _PWMIG_H_
/* This incldue file is for function prototypes for procedures used in
pwmig.   Most will not be wrapped for python */
#include <vector>
#include "mspass/seismic/Seismogram.h"
#include "mspass/seismic/Ensemble.h"
#include "mspass/seismic/SlownessVector.h"
#include "pwmig/gclgrid/gclgrid.h"
#include "pwmig/gclgrid/PWMIGfielddata.h"
#include "pwmig/pwmigcore/SlownessVectorMatrix.h"
#include "pwmig/pwmigcore/PWMIGmigrated_seismogram.h"
#include "pwmig/seispp/VelocityModel_1d.h"
namespace pwmig::pwmigcore
{
/* These were in the old pwmig.h.   */
mspass::seismic::SlownessVector get_stack_slowness
  (mspass::seismic::LoggingEnsemble<mspass::seismic::Seismogram>& ensemble);
//Hypocenter get_event_hypocenter(ThreeComponentEnsemble& ensemble);
pwmig::gclgrid::GCLscalarfield3d *Build_GCLraygrid(
       pwmig::gclgrid::GCLgrid& parent,
	pwmig::pwmigcore::SlownessVectorMatrix& svm,
	 pwmig::seispp::VelocityModel_1d& vmod,
          double zmax, double tmax, double dt);
mspass::utility::dmatrix *ray_path_tangent(mspass::utility::dmatrix&);
pwmig::gclgrid::GCLscalarfield3d *decimate(pwmig::gclgrid::GCLscalarfield3d& g,
	int dec1, int dec2, int dec3);
/* these functions used to be in the source file for pwmig.cc.  Moved for
mspass to pwmigsubs.cc. These procedures to not use const properly
but since they predate my learning the proper use of const and even clean
compiler support for const.   they are (mostly) not changed because I do not expect
these procedures to be used outside the revised pwmig.*/
void remove_mean_x3(pwmig::gclgrid::GCLscalarfield3d& f);
void reverse_time_field(pwmig::gclgrid::GCLscalarfield3d& f, string current_order);
double pathdist(mspass::utility::dmatrix& ray,int i0);
double fudge_model_upward(pwmig::gclgrid::GCLscalarfield3d& model,double dz);
std::vector<double> compute_3Dmodel_time(pwmig::gclgrid::GCLscalarfield3d& U3d,
	mspass::utility::dmatrix& path);
void cosine_taper_highend(mspass::utility::dmatrix& d,int mark, int taper_length);
mspass::utility::dmatrix compute_gradS(pwmig::gclgrid::GCLgrid3d& raygrid,
	int ix1, int ix2, pwmig::seispp::VelocityModel_1d& Vs1d);
std::vector<double> compute_unit_normal_P(pwmig::gclgrid::GCLgrid3d& raygrid,
	double x1, double x2, double x3,int& ix1_0, int& ix2_0, int& ix3_0);
pwmig::pwmigcore::SlownessVectorMatrix pad_svm(const pwmig::pwmigcore::SlownessVectorMatrix& svm,
	const int pad);
bool grid_mismatched(const pwmig::gclgrid::GCLgrid& parent,
	const pwmig::gclgrid::GCLgrid& padded, const int pad);
/* this one needs some hacking to deal with the auto_ptr */
pwmig::gclgrid::GCLscalarfield3d* ComputeIncidentWaveRaygrid(pwmig::gclgrid::GCLgrid& pstagrid,
  int border_pad,
   pwmig::gclgrid::GCLscalarfield3d& UP3d,
    pwmig::seispp::VelocityModel_1d vp1d,
     pwmig::pwmigcore::SlownessVectorMatrix& svm,
      double zmax,
       double tmax,
        double dt,
         int zdecfac,
          bool use_3d);
std::vector<double> compute_Stime(pwmig::gclgrid::GCLscalarfield3d& U3d,
	  int i, int j, pwmig::gclgrid::GCLscalarfield3d& raygrid,bool use_3d);
std::vector<double> compute_domega_for_path(mspass::seismic::SlownessVector& u0,
	double dux,
	 double duy,
	  pwmig::seispp::VelocityModel_1d& vmod,
		 double zmax,
		  double dz,
		   pwmig::gclgrid::GCLgrid3d& g,
			  int ix1,
				 int ix2,
		      mspass::utility::dmatrix& gradP,
		       std::vector<double>& zP);
std::vector<double> compute_weight_for_path(mspass::utility::dmatrix& gradTp,
	mspass::utility::dmatrix& gradTs);
std::vector<double> running_average(std::vector<double>& x, int ns);
mspass::utility::dmatrix running_average(mspass::utility::dmatrix x, int ns);
mspass::seismic::SlownessVector slowness_average(mspass::seismic::LoggingEnsemble<mspass::seismic::Seismogram> *d);
pwmig::gclgrid::Geographic_point fetch_TP_point_at_base(pwmig::gclgrid::GCLscalarfield3d& TP,
	int i, int j);
pwmig::gclgrid::Geographic_point get_gp_base_TPx(pwmig::gclgrid::GCLscalarfield3d& TP,
	pwmig::gclgrid::Geographic_point xgp);
pwmig::gclgrid::Geographic_point find_TP_at_x_depth(pwmig::gclgrid::GCLscalarfield3d& TP,
	int i, int j, double zx);
double compute_delta_p_term(pwmig::gclgrid::Geographic_point r0x,
	pwmig::gclgrid::Geographic_point r0, mspass::seismic::SlownessVector u0);
mspass::seismic::SlownessVector EnsembleDeltaSlow(mspass::seismic::LoggingEnsemble<mspass::seismic::Seismogram> *d);
void VelocityFieldToSlowness(pwmig::gclgrid::GCLscalarfield3d& g);
pwmig::seispp::VelocityModel_1d DeriveVM1Dfrom3D(pwmig::gclgrid::GCLscalarfield3d& u);
pwmig::pwmigcore::PWMIGmigrated_seismogram migrate_one_seismogram(mspass::seismic::Seismogram& pwdata,
    pwmig::gclgrid::GCLgrid& parent,
      pwmig::gclgrid::GCLscalarfield3d& raygrid,
        pwmig::gclgrid::GCLscalarfield3d& TPgrid,
           pwmig::gclgrid::GCLscalarfield3d& Us3d,
             pwmig::seispp::VelocityModel_1d& Vp1d,
               pwmig::seispp::VelocityModel_1d& Vs1d,
                 mspass::utility::Metadata& control);
/*! Multithreaded function to migrate one plane wave component.

This function takes input assumed to be contained in the input ThreeComponentEnsemble
defining data for one plane wave as output by pwstack.  I uses a ray-based
algorithm to backproject each datum onto a ray path computed from the slowness
vector and the input S wave velocity model.  The output of each such operation
is a Seismogram like object that is then copied into the PWMIGfielddata grid
at the correct index grid position.
*/
pwmig::gclgrid::PWMIGfielddata migrate_component(
  mspass::seismic::ThreeComponentEnsemble& d,
    pwmig::pwmigcore::SlownessVectorMatrix& VPsvm,
      pwmig::gclgrid::GCLgrid& parent,
        pwmig::gclgrid::GCLscalarfield3d& TPgrid,
           pwmig::gclgrid::GCLscalarfield3d& Us3d,
             pwmig::seispp::VelocityModel_1d& Vp1d,
               pwmig::seispp::VelocityModel_1d& Vs1d,
                 mspass::utility::Metadata& control);
} // end namespace
#endif
