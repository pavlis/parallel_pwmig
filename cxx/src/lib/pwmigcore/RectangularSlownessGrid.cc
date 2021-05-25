#include <float.h>
#include "mspass/utility/MsPASSError.h"
#include "mspass/utility/AntelopePf.h"
#include "pwmig/pwmigcore/RectangularSlownessGrid.h"
namespace pwmig::pwmigcore
{
using namespace std;
using namespace pwmig::pwmigcore;
using mspass::utility::AntelopePf;
using mspass::utility::MsPASSError;
using mspass::utility::ErrorSeverity;
// Generic default slowness grid
RectangularSlownessGrid::RectangularSlownessGrid()
{
	name=string("Default_Slowness_Grid");
	uxlow = - 0.5;
	uylow = -0.5;
	nux = 101;
	nuy = 101;
	dux = 0.01;
	duy = 0.01;
}

// brute force constructor
RectangularSlownessGrid::RectangularSlownessGrid(const string nm,
	const double uxl,
		const double uyl,
			const double du1,
				const double du2,
					const int n1,
						const int n2)
{
	name =nm;
	uxlow = uxl;
	uylow = uyl;
	nux = n1;
	nuy = n2;
	dux =du1;
	duy=du2;
}

RectangularSlownessGrid::RectangularSlownessGrid(const AntelopePf& pfsmd,const string tag)
{
    /* Painfully parallel to the pf version above */
	try {
    AntelopePf md=pfsmd.get_branch(tag);
		name=md.get_string("Slowness_Grid_Name");
		uxlow=md.get_double("uxlow");
		uylow=md.get_double("uylow");
		nux = md.get_int("nux");
		nuy = md.get_int("nuy");
		dux = md.get_double("dux");
		duy = md.get_double("duy");
	} catch (...) {throw;}
}
// Copy constructor needed due to string variable (always wise anyway they say)
RectangularSlownessGrid::RectangularSlownessGrid(const RectangularSlownessGrid& rsg)
{
	name=rsg.name;
	uxlow=rsg.uxlow;
	uylow=rsg.uylow;
	nux=rsg.nux;
	nuy=rsg.nuy;
	dux=rsg.dux;
	duy=rsg.duy;
}
// Returns a slowness vector for a grid position i,j
SlownessVector RectangularSlownessGrid::slow(const int i, const int j) const
{
	if(i>=nux || j>=nuy || i<0 || j<0)
	{
		stringstream ss;
		ss << "RectangularSlownessGrid::slow method:  grid index ("<<i<<","<<") is illegal"<<endl
		  <<  "Both must be positive with upper bounds ("<<nux-1<<","<<nuy-1<<")";
		throw MsPASSError(ss.str(),ErrorSeverity::Invalid);
	}
	SlownessVector u;
	u.ux=uxlow+i*dux;
	u.uy=uylow+j*duy;
	return(u);
}
RectangularSlownessGrid& RectangularSlownessGrid::operator=(const RectangularSlownessGrid& parent)
{
	if(this != &parent)
	{
		name=parent.name;
		uxlow=parent.uxlow;
		uylow=parent.uylow;
		nux=parent.nux;
		nuy=parent.nuy;
		dux=parent.dux;
		duy=parent.duy;
	}
	return *this;
}
/* Because of the implementation of implicit value with a regular grid geometry
this particular operator is trivial */
RectangularSlownessGrid& RectangularSlownessGrid::operator+=(const mspass::seismic::SlownessVector& slowvector)
{
	this->uxlow += slowvector.ux;
	this->uylow += slowvector.uy;
	return *this;
}

}  // End namespace SEISPP
