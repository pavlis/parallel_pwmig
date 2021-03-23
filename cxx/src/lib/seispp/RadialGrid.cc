#include <sstream>
#include <vector>
/* We need this for the rad macro and dist calculator */
#include "pwmig/dsap/coords.h"
#include "pwmig/gclgrid/gclgrid.h"
#include "pwmig/seispp/Hypocenter.h"
#include "pwmig/seispp/RadialGrid.h"
#include "mspass/utility/AntelopePf.h"
#include "mspass/utility/MsPASSError.h"
namespace pwmig::seispp{
using namespace std;
using namespace pwmig::seispp;
using namespace pwmig::gclgrid;
using mspass::utility::Metadata;
using mspass::utility::AntelopePf;
using mspass::utility::MsPASSError;
RadialGrid::RadialGrid()
{
	lat0=0.0;
	lon0=0.0;
	naz=0;
	ndelta=0;
	ndelbins=0;
	nazbins=0;
}
RadialGrid::RadialGrid(const AntelopePf& pf)
{
	string base_error("RadialGrid AntelopePf constructor:  ");
	lat0=pf.get_double("origin_latitude");
	lon0=pf.get_double("origin_longitude");
	lat0=rad(lat0);
	lon0=rad(lon0);
	list<string> ptslist;
	/* this will throw an exception if the key is missing.  We let it do
	that for the expected context. */
	ptslist=pf.get_tbl("delta_grid_points");
	for(auto sptr=ptslist.begin();sptr!=ptslist.end();++sptr)
	{
		double valin;
		stringstream ss(*sptr);
		ss >> valin;
		if( (valin<0.0) || (valin>180.0) )
			throw MsPASSError(base_error
				+ string("Illegal distance specification=")
				+ (*sptr));
		delta.push_back(rad(valin));
	}
	ndelta=delta.size();
	ndelbins=ndelta-1;
	ptslist=pf.get_tbl("azimuth_grid_points");
	for(auto sptr=ptslist.begin();sptr!=ptslist.end();++sptr)
	{
		double valin;
		stringstream ss(*sptr);
		ss >> valin;
		if( (valin<-180.0) || (valin>360.0) )
			throw MsPASSError(base_error
				+ string("Illegal azimuth specification=")
				+ (*sptr));
		azimuth.push_back(rad(valin));
	}
	naz=azimuth.size();
	nazbins=naz-1;
}
RadialGrid::RadialGrid(const double azmin, const double azmax, const int nazin,
		const double delmin, const double delmax, const int ndelin,
		const double lat0in, const double lon0in)
{
	const string base_error("RadialGrid parameterized constructor:  ");
	if(azmin>=azmax) throw MsPASSError(base_error
		+ "Illegal azimuth range specified");
	if(delmin>=delmax) throw MsPASSError(base_error
		+ "Illegal distance range specified");
	if( (nazin<=0) || (ndelin<=0) ) throw MsPASSError(base_error
		+ "Number of points in azimuth and delta must be positive");
	naz=nazin;
	ndelta=ndelin;
	nazbins=naz-1;
	ndelbins=ndelta-1;
	lat0=rad(lat0in);
	lon0=rad(lon0in);
	double daz=(azmax-azmin)/static_cast<double>(nazbins);
	double ddel=(delmax-delmin)/static_cast<double>(ndelbins);
	int i;
	for(i=0;i<naz;++i)
	  azimuth.push_back(rad(azmin+daz*static_cast<double>(i)));
	for(i=0;i<ndelta;++i)
	  delta.push_back(rad(delmin+ddel*static_cast<double>(i)));
}
double RadialGrid::lat(const int ir, const int id)
{
    try {
	Geographic_point p=this->grid_point(ir,id);
	return(p.lat);
    } catch (...) {throw;};
}
double RadialGrid::lon(const int ir, const int id)
{
    try {
	Geographic_point p=this->grid_point(ir,id);
	return(p.lon);
    } catch (...) {throw;};
}
/* Note the symbols are very confusing here.  ir is the azimuth and id is
for distance axis.  include file is more rational - I didn't change this
to reduce coding errors*/
Geographic_point RadialGrid::grid_point(const int ir, const int id)
{
	if( (ir<0) || (ir>=nazbins) || (id<0) || (id>=ndelbins) )
	{
		stringstream ss;
		ss << "RadialGrid:  requested index (az,del)=("
			<<ir<<", "<<id<<") out of range"<<endl
			<<"Allowed: distance=(0,"<<ndelbins
			<<") azimuth=(0,"<<nazbins<<")"<<endl;
		throw MsPASSError(ss.str());
	}
	double plat,plon,pr;
	double delcenter,azcenter;
	delcenter=(delta[id]+delta[id+1])/2.0;
	azcenter=(azimuth[ir]+azimuth[ir+1])/2.0;
	latlon(lat0,lon0,delta[id],azimuth[ir],&plat,&plon);
	latlon(lat0,lon0,delcenter,azcenter,&plat,&plon);
	pr=r0_ellipse(plon);
	Geographic_point result;
	result.lat=plat;
	result.lon=plon;
	result.r=pr;
	return(result);
}
Metadata RadialGrid::cell(const int ia, const int id)
{
	try{
		Metadata md;
		/* Note this method returns the grid cell center defined as centroid
		of the bounding points */
		Geographic_point gpt=this->grid_point(ia,id);
		/* store in degrees for metadata  here as the routine is an interface
		to mongodb*/
		md.put("lat",deg(gpt.lat));
		md.put("lon",deg(gpt.lon));
		md.put("delta_minimum",deg(this->delta[id]));
		md.put("delta_maximum",deg(this->delta[id+1]));
		md.put("azimuth_minimum",deg(this->azimuth[ia]));
		md.put("azimuth_maximum",deg(this->azimuth[ia+1]));
		md.put("azimuth_index",ia);
		md.put("distance_index",id);
		return md;
	}catch(...){throw;};
}

SectorTest::SectorTest(RadialGrid& g, const int ia, const int id)
{
	const string base_error("SectorTest constructor:  ");
	if( (ia<0) || (ia>=(g.naz-1) ))
		throw MsPASSError(base_error
		 + "Requested azimuth index outside allowed range");
	if( (id<0) || (id>=(g.ndelta-1)) )
		throw MsPASSError(base_error
		 + "Requested distance index outside allowed range");
	lat0=g.lat0;
	lon0=g.lon0;
	azmin=g.azimuth[ia];
	delmin=g.delta[id];
	azmax=g.azimuth[ia+1];
	delmax=g.delta[id+1];
}
} // end namespace
