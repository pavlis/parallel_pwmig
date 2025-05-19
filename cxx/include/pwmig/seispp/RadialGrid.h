#ifndef _RADIALGRID_H_
#define _RADIALGRID_H_
#include <string>
#include "mspass/utility/AntelopePf.h"
#include "pwmig/gclgrid/gclgrid.h"
/* Include this explicit prototype here normally foudn in coords.h.
that file has some namespace collisions.  To avoid them we just give
the prototype for dist */
extern "C" {
	extern void dist ( double xlat1, double xlong1, double xlat2, double xlong2,
		double *del, double *az );
};
namespace pwmig::seispp
{
/*! \brief Define a global radial grid.

A radial grid on the earth is best thought a set of generalized coordinates
with the an origin that acts like a pole.  Great circle paths to
a specified distance in degrees are defined from that point at uniform
(specified) angles.  This operator maps a port or all of a sphere onto
a circle.   The sectors in the resulting circle define positions within
a range of angles and distances for great circle paths emenating from the
pole.   This grid is used in the telecluster program of pwmig to
assemble sources in a epicentral distance and azimuth bins.  The grid
can be uniform or nonuniform.  Use the parameterized constructor to
create a uniform grid and the AntelopePf constructor to create a
grid with variable distance or azimuth grids.
*/
class RadialGrid
{
public:
	/*! Pole latitude. */
	double lat0;
	/*! Pole longitude. */
	double lon0;
	/*! Number of azimuth grid points. */
	int naz;
	/*! Number of distance grid points. */
	int ndelta;
	/*! Vector of distances defining the grid (in radians).*/
	vector<double> delta;
	/*! Vector of azimuths defining the grid (in radians).  These bin edges
	not grid centers with 0 being the minimum (usually 0) and naz-1 being
	the maximum (usually 360-daz)*/
	vector<double> azimuth;
	/*! Default constuctor.  Initializes all attributes to 0*/
	RadialGrid();
	/*! Construct from an AntelopePf.

	This constructor should be used if you need a nonuniform grid.
	Two Tbl entries with keys  delta_grid_points and azimuth_grid_points
	can be used to define the grid by specific points. In that format
	origin_latitude and origin_longitude are keys used to set lat0 and lon0.
	*/
	RadialGrid(const mspass::utility::AntelopePf& pf);
	/*! Fully parameterized constructor.

	Use this constructor to create a uniform grid in azimuth and distance
	between the specified boundaries.

	\param azmin is the starting point azimuth.
	\param azmax is the largest azimuth to consider (edge not center)
	\param naz is the number of azimuth bins to create (note this is bins not points)
	\param delmin is the starting epicenteral distance (in degrees)
	\paam delmax is the maximum epicentral distance for the grid boundary.
	\param ndel is the number of distance bins (note this is bins not points)
	\param lat0 is the latitude (in degrees) of requested origin.
	\param lon0 is the longitude (in degrees) of requested grid origin.
	*/
	RadialGrid(const double azmin, const double azmax, const int naz,
		const double delmin, const double delmax, const int ndel,const double lat0, const double lon0);
	/*! Return the coordinates of the grid cell center.

	\param iaz is the index position for azimuth of requested cell
	\param idel is the index position in distance of requested cell
	 */
	pwmig::gclgrid::Geographic_point grid_point(const int iaz, const int idel);
	/*! Return latitude in radians of grid index position */
	double lat(const int ir, const int id);
	/*! Return latitude in radians of grid index position */
	double lon(const int ir, const int id);
	/*! Return number of cells in distance (one less than grid points)*/
	int number_distance_bins(){return ndelbins;};
	/*! Return number of cells in azimuth (one less than gride points)*/
	int number_azimuth_bins(){return nazbins;}
	/*! Return a Metadata container with set of attributes fully describing
	the cell - cornerss defined as lat,lon and distanc and azimuth ranges */
	mspass::utility::Metadata cell(const int ia, const int ir);
private:
	/* These were retained out of laziness. They are always one less than
	ndel and naz respectively. */
	double ndelbins;
	double nazbins;
};
/* This function object is used as the predicate for the subset method template in
EventCatalog for a radial grid.  */
class SectorTest
{
public:
	bool operator() (const Hypocenter& h) const
	{
		double del,az;
		dist(lat0,lon0,h.lat,h.lon,&del,&az);
		if( (del<=delmin) || (del>delmax) ) return false;
		if( (az<=azmin) || (az>azmax) ) return false;
		return true;
	};
	SectorTest(RadialGrid& g,const int ia, const int id);
private:
	double delmin,delmax;
	double azmin,azmax;
	double lat0,lon0;
};
}
#endif
