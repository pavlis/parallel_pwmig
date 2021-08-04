#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/operators.h>
#include <pybind11/embed.h>

#include "mspass/utility/Metadata.h"
#include "mspass/utility/AntelopePf.h"
#include "pwmig/gclgrid/gclgrid.h"
#include "pwmig/gclgrid/RegionalCoordinates.h"

/* these are needed to make std::vector containers function propertly.
This was borrowed from mspass pybind11 cc files */
PYBIND11_MAKE_OPAQUE(std::vector<double>);

namespace pwmig {
namespace pwmigpy {

namespace py=pybind11;
using namespace std;
using mspass::utility::Metadata;
using mspass::utility::AntelopePf;
using namespace pwmig::gclgrid;
/* This is what pybind11 calls a trampoline class needed to handle
virtual base classes for gclgrid objects. */
class PyBasicGCLgrid : public BasicGCLgrid
{
public:
  using BasicGCLgrid::BasicGCLgrid;
  void compute_extents() override
  {
    PYBIND11_OVERLOAD_PURE(
      void,
      BasicGCLgrid,
      compute_extents
    );
  }
  void reset_index() override
  {
    PYBIND11_OVERLOAD_PURE(
      void,
      BasicGCLgrid,
      reset_index
    );
  }
  void get_index(int *index) override
  {
    PYBIND11_OVERLOAD_PURE(
      void,
      BasicGCLgrid,
      get_index,
      index
    );
  }
  mspass::utility::Metadata get_attributes() override
  {
    PYBIND11_OVERLOAD_PURE(
      mspass::utility::Metadata,
      BasicGCLgrid,
      get_attributes
    );
  }
  pwmig::gclgrid::Geographic_point ctog(const double x1,
    const double x2, const double x3) const
  {
    PYBIND11_OVERLOAD(
      pwmig::gclgrid::Geographic_point,
      BasicGCLgrid,
      ctog,
      x1,
      x2,
      x3
    );
  }
  pwmig::gclgrid::Geographic_point ctog(const pwmig::gclgrid::Cartesian_point point) const
  {
    PYBIND11_OVERLOAD(
      pwmig::gclgrid::Geographic_point,
      BasicGCLgrid,
      ctog,
      point
    );
  }
  pwmig::gclgrid::Cartesian_point gtoc(const double lat,
    const double lon, const double radius) const
  {
    PYBIND11_OVERLOAD(
      pwmig::gclgrid::Cartesian_point,
      BasicGCLgrid,
      gtoc,
      lat,
      lon,
      radius
    );
  }
  pwmig::gclgrid::Cartesian_point gtoc(const pwmig::gclgrid::Geographic_point point) const
  {
    PYBIND11_OVERLOAD(
      pwmig::gclgrid::Cartesian_point,
      BasicGCLgrid,
      gtoc,
      point
    );
  }
  double depth(const pwmig::gclgrid::Cartesian_point point) const
  {
    PYBIND11_OVERLOAD(
      double,
      BasicGCLgrid,
      depth,
      point
    );
  }
  double depth(const pwmig::gclgrid::Geographic_point point) const
  {
    PYBIND11_OVERLOAD(
      double,
      BasicGCLgrid,
      depth,
      point
    );
  }

};


PYBIND11_MODULE(gclgrid, m) {
/* This is needed to allow vector inputs and outputs */
py::bind_vector<std::vector<double>>(m, "DoubleVector");

py::class_<pwmig::gclgrid::Geographic_point>(m,"Geographic_point","Point on Earth defined in regional cartesian system")
  .def(py::init<>())
  .def(py::init<const Geographic_point&>())
  .def_readwrite("lat",&Geographic_point::lat,"Latitude of a point (radians)")
  .def_readwrite("lon",&Geographic_point::lon,"Longitude of a point (radians)")
  .def_readwrite("r",&Geographic_point::r,"Radial distance from Earth center (km)")
  ;
py::class_<pwmig::gclgrid::Cartesian_point>(m,"Cartesian_point","Point on Earth defined coordinates in radians")
    .def(py::init<>())
    .def(py::init<const Cartesian_point&>())
    .def("coordinates",&Cartesian_point::coordinates,"Return vector of coordinates in standard order")
    .def_readwrite("x1",&Cartesian_point::x1,"x1 coordinate axis value (km)")
    .def_readwrite("x2",&Cartesian_point::x2,"x2 coordinate axis value (km)")
    .def_readwrite("x3",&Cartesian_point::x3,"x3 coordinate axis value (km)")
    ;
py::class_<BasicGCLgrid,PyBasicGCLgrid>(m,"BasicGCLgrid","Base class for family of GCL data objects")
  .def(py::init<>())
  .def("set_transformation_matrix",&BasicGCLgrid::set_transformation_matrix,
      "Sets the tranformation matrix from current values of r0, lat0, lon0, and azimuth_y")
  .def("fetch_transformation_matrix",&BasicGCLgrid::fetch_transformation_matrix,
       "Retrieve the transformation matrix defined for this coordinate system")
  .def("fetch_translation_vector",&BasicGCLgrid::fetch_translation_vector,
       "Retrieve the translation vector of this coordinate system")
  .def("ctog",py::overload_cast<const double, const double, const double>
        (&BasicGCLgrid::ctog,py::const_),"Convert grid Cartesian coordinates to geographic")
  .def("ctog",py::overload_cast<const pwmig::gclgrid::Cartesian_point>
        (&BasicGCLgrid::ctog,py::const_),"Convert grid Cartesian coordinates to geographic")
  .def("gtoc",py::overload_cast<const double, const double, const double>
      (&BasicGCLgrid::gtoc,py::const_),"Convert geographic point to grid cartesian system")
  .def("gtoc",py::overload_cast<const pwmig::gclgrid::Geographic_point>
      (&BasicGCLgrid::gtoc,py::const_),"Convert geographic point to grid cartesian system")
  .def("depth",py::overload_cast<const pwmig::gclgrid::Cartesian_point>(&BasicGCLgrid::depth,py::const_),
      "Return depth from sea level reference ellipsoid of point specified in grid cartesian system")
  .def("depth",py::overload_cast<const pwmig::gclgrid::Geographic_point>(&BasicGCLgrid::depth,py::const_),
      "Return depth from sea level reference ellipsoid of point specified in spherical (geographic) coordinates (not ellipsoid corrected)")
  /* these are pure virtual methods but they still need to be defined here
  to get pybind11 to compile correctly */
  .def("compute_extents",&BasicGCLgrid::compute_extents)
  .def("reset_index",&BasicGCLgrid::reset_index)
  .def("get_index",&BasicGCLgrid::get_index)
  /* public attributes */
  .def_readwrite("name",&BasicGCLgrid::name,"Unique name assigned to this grid object")
  .def_readwrite("lat0",&BasicGCLgrid::lat0)
  .def_readwrite("lon0",&BasicGCLgrid::lon0)
  .def_readwrite("r0",&BasicGCLgrid::r0)
  .def_readwrite("azimuth_y",&BasicGCLgrid::azimuth_y)
  .def_readwrite("dx1_nom",&BasicGCLgrid::dx1_nom)
  .def_readwrite("dx2_nom",&BasicGCLgrid::dx2_nom)
  .def_readwrite("n1",&BasicGCLgrid::n1)
  .def_readwrite("n2",&BasicGCLgrid::n2)
  .def_readwrite("i0",&BasicGCLgrid::i0)
  .def_readwrite("j0",&BasicGCLgrid::j0)
  .def_readwrite("x1low",&BasicGCLgrid::x1low)
  .def_readwrite("x2low",&BasicGCLgrid::x2low)
  .def_readwrite("x3low",&BasicGCLgrid::x3low)
  .def_readwrite("x1high",&BasicGCLgrid::x1high)
  .def_readwrite("x2high",&BasicGCLgrid::x2high)
  .def_readwrite("x3high",&BasicGCLgrid::x3high)
;

py::class_<GCLgrid,BasicGCLgrid>(m,"GCLgrid",py::buffer_protocol(),
                  "Two-dimensional GCL grid object")
  .def(py::init<>())
  .def(py::init<const int, const int>())
  .def(py::init<const int, const int, const string,const double, const double,
    const double, const double, const double, const double,
    const int, const int>())
  .def(py::init<const string, const string>())
  .def(py::init<const GCLgrid&>())
  .def(py::init<const Metadata&>())
  .def("save",&GCLgrid::save,"Save to an external file")
  .def("lookup",&GCLgrid::lookup,"Find point by cartesian coordinates")
  .def("reset_index",&GCLgrid::reset_index,"Initializer for lookup searches - rarely needed")
  .def("get_index",&GCLgrid::get_index,"Return index position found with lookup")
  .def("lat",&GCLgrid::lat,"Get latitude (radians) of a grid point specified by two index ints")
  .def("lon",&GCLgrid::lon,"Get longitude (radians) of a grid point specified by two index ints")
  .def("r",&GCLgrid::r,"Get radial distance from earth center (km) of a grid point specified by two index ints")
  .def("depth",&GCLgrid::depth,"Get depth (km) from 0 reference ellipsoid radius of a grid point specified by two index ints")
  .def("compute_extents",&GCLgrid::compute_extents,"Call after manually building a grid")
  .def("geo_coordinates",&GCLgrid::geo_coordinates,"Fetch grid point defined as geo coordinates")
  .def("get_coordinates",&GCLgrid::get_coordinates,"Fetch grid point defined in Cartesian system")
  .def("get_attributes",&GCLgrid::get_attributes,
     "Fetch all attributes into a Metadata container")
  .def("set_coordinates",py::overload_cast<const Cartesian_point&,
     const int, const int>(&GCLgrid::set_coordinates),
     "Set grid point coordinates using Cartesian system")
  .def("set_coordinates",py::overload_cast<const Cartesian_point&,
      const int, const int>(&GCLgrid::set_coordinates),
      "Set grid point coordinates Geographic point (radians)")
;
py::class_<GCLgrid3d,BasicGCLgrid>(m,"GCLgrid3d",py::buffer_protocol(),
                     "Three-dimensional GCL grid object")
  .def(py::init<>())
  .def(py::init<const int, const int, const int>())
  .def(py::init<const int, const int, const int, const string,
    const double, const double,
    const double, const double, const double, const double, const double,
    const int, const int>())
  .def(py::init<const string, const string,const bool>())
  .def(py::init<const GCLgrid3d&>())
  .def(py::init<const Metadata&>())
  .def("save",&GCLgrid3d::save,"Save to an external file")
  .def("lookup",&GCLgrid3d::lookup,"Find point by cartesian coordinates")
  .def("reset_index",&GCLgrid3d::reset_index,"Initializer for lookup searches - rarely needed")
  .def("get_index",&GCLgrid3d::get_index,"Return index position found with lookup")
  .def("geo_coordinates",&GCLgrid3d::geo_coordinates,"Return geo coordinate struct")
  .def("get_coordinates",&GCLgrid3d::get_coordinates,"Fetch grid point defined in Cartesian system")
  .def("lat",&GCLgrid3d::lat,"Get latitude (radians) of a grid point specified by three index ints")
  .def("lon",&GCLgrid3d::lon,"Get longitude (radians) of a grid point specified by three index ints")
  .def("r",&GCLgrid3d::r,"Get radial distance from earth center (km) of a grid point specified by three index ints")
  .def("depth",&GCLgrid3d::depth,"Get depth (km) from 0 reference ellipsoid radius of a grid point specified by three index ints")
  .def("compute_extents",&GCLgrid3d::compute_extents,"Call after manually building a grid")
  .def("get_attributes",&GCLgrid3d::get_attributes,
     "Fetch all attributes into a Metadata container")
  .def("set_coordinates",py::overload_cast<const Cartesian_point&,
        const int, const int, const int>(&GCLgrid3d::set_coordinates),
        "Set grid point coordinates using Cartesian system")
  .def("set_coordinates",py::overload_cast<const Cartesian_point&,
         const int, const int,const int>(&GCLgrid3d::set_coordinates),
         "Set grid point coordinates Geographic point (radians)")
  .def_readwrite("n3",&GCLgrid3d::n3)
  .def_readwrite("dx3_nom",&GCLgrid3d::dx3_nom)
  .def_readwrite("k0",&GCLgrid3d::k0)
;
py::class_<GCLscalarfield,GCLgrid>(m,"GCLscalarfield","Two-dimensional grid with scalar attributes at each node")
  .def(py::init<>())
  .def(py::init<const int, const int>())
  .def(py::init<const GCLgrid&>())
  .def(py::init<const string, const string, const bool>())
  .def(py::init<const GCLscalarfield&>())
  .def(py::init<const Metadata&>())
  .def("zero",&GCLscalarfield::zero,"Set all field attributes to 0")
  .def("save",&GCLscalarfield::save,"Save contents to a file")
  .def("interpolate",&GCLscalarfield::interpolate,"Interpolate grid to get value at point passed")
  .def("get_attributes",&GCLscalarfield::get_attributes,
     "Fetch all attributes into a Metadata container")
  .def("set_value",&GCLscalarfield::set_value,"Set the value at a specified grid point")
  .def("get_value",&GCLscalarfield::get_value,"Get the value at a specified grid point")
   /* This is normally the right syntax in pybind11 for operator+= but
   not working here for some mysterious reason. Put aside until needed - solvable
   problem just one of those annoying picky pybind11 details*/
  //.def(py::self += py::self)
  .def(py::self *= double())
;

py::class_<GCLvectorfield,GCLgrid>(m,"GCLvectorfield","Two-dimensional grid with vector attributes at each node")
  .def(py::init<>())
  .def(py::init<const int, const int,const int>())
  .def(py::init<const GCLgrid&,const int>())
  .def(py::init<const string, const string, const bool>())
  .def(py::init<const GCLvectorfield&>())
  .def(py::init<const Metadata&>())
  .def("zero",&GCLvectorfield::zero,"Set all field attributes to 0")
  .def("save",&GCLvectorfield::save,"Save contents to a file")
  .def("interpolate",&GCLvectorfield::interpolate,"Interpolate grid to get vector values at point passed")
  .def("get_attributes",&GCLvectorfield::get_attributes,
     "Fetch all attributes into a Metadata container")
  .def("get_value",&GCLvectorfield::get_value,"Get the vector of values at a specified grid point")
  .def("set_value",&GCLvectorfield::set_value,"Set the value at a specified grid point")
  //.def(py::self += py::self)
  .def(py::self *= double())
  .def_readwrite("nv",&GCLvectorfield::nv,"Number of components in each vector")
;
// 3D versions of two class definitions above
py::class_<GCLscalarfield3d,GCLgrid3d>(m,"GCLscalarfield3d","Three-dimensional grid with scalar attributes at each node")
  .def(py::init<>())
  .def(py::init<const int, const int, const int>())
  .def(py::init<const GCLgrid3d&>())
  .def(py::init<const string, const string>())
  .def(py::init<const GCLscalarfield3d&>())
  .def(py::init<const Metadata&>())
  .def("zero",&GCLscalarfield3d::zero,"Set all field attributes to 0")
  .def("save",&GCLscalarfield3d::save,"Save contents to a file")
  .def("interpolate",&GCLscalarfield3d::interpolate,"Interpolate grid to get value at point passed")
  .def("get_attributes",&GCLscalarfield3d::get_attributes,
     "Fetch all attributes into a Metadata container")
  .def("get_value",&GCLscalarfield3d::get_value,"Get the value at a specified grid point")
  .def("set_value",&GCLscalarfield3d::set_value,"Set the value at a specified grid point")
  .def(py::self *= double())
;

py::class_<GCLvectorfield3d,GCLgrid3d>(m,"GCLvectorfield3d","Three-dimensional grid with vector attributes at each node")
  .def(py::init<>())
  .def(py::init<const int, const int,const int,const int>())
  .def(py::init<const GCLgrid3d&,const int>())
  .def(py::init<const string, const string>())
  .def(py::init<const GCLvectorfield3d&>())
  .def(py::init<const Metadata&>())
  .def("zero",&GCLvectorfield3d::zero,"Set all field attributes to 0")
  .def("save",&GCLvectorfield3d::save,"Save contents to a file")
  .def("interpolate",&GCLvectorfield3d::interpolate,"Interpolate grid to get vector values at point passed")
  .def("get_attributes",&GCLvectorfield3d::get_attributes,
     "Fetch all attributes into a Metadata container")
  .def("get_value",&GCLvectorfield3d::get_value,"Get the vector of values at a specified grid point")
  .def("set_value",&GCLvectorfield3d::set_value,"Set the value at a specified grid point")
  //.def(py::self += py::self)
  .def(py::self *= double())
  .def_readwrite("nv",&GCLvectorfield3d::nv,"Number of components in each vector")
;
py::class_<RegionalCoordinates>(m,"RegionalCoordinates",
      "Encapsulates coordinate system used in gclgrid objects")
  .def(py::init<>())
  .def(py::init<const double, const double, const double, const double>())
  .def("cartesian",py::overload_cast<const double, const double, const double>
     (&RegionalCoordinates::cartesian,py::const_),"Return cartesian translation of lat, lon, radius")
  .def("cartesian",py::overload_cast<const Geographic_point>
     (&RegionalCoordinates::cartesian,py::const_),
     "Return cartesian translation of point in geographic struct")
  .def("geographic",py::overload_cast<const double, const double, const double>
    (&RegionalCoordinates::geographic,py::const_),
    "Return geographic points equivalent to three components of a cartesian vector")
  .def("geographic",py::overload_cast<const Cartesian_point>
    (&RegionalCoordinates::geographic,py::const_),
    "Return geographic points equivalent to three components of a cartesian vector")
  .def("aznorth_angle",&RegionalCoordinates::aznorth_angle,"Return internal azimuth angle of x2 axis relative to north")
  .def("origin",&RegionalCoordinates::origin,"Return geographic location of origin of coordinate system")
;

/* gclgrid functions.  Not all are wrapped here - will add them as I need them */
m.def("r0_ellipse",&r0_ellipse,
  "Return the radius of the reference ellipsoid at latitude specified in radians",
  py::return_value_policy::copy,
  py::arg("lat") )
;
m.def("remap_grid",py::overload_cast<GCLgrid&, const BasicGCLgrid&>(&remap_grid),
    "Change coordinate system of a grid to match another",
  py::return_value_policy::copy,
  py::arg("g"),
  py::arg("parent") )
;
m.def("remap_grid",py::overload_cast<GCLgrid3d&, const BasicGCLgrid&>(&remap_grid),
    "Change coordinate system of a grid to match another",
  py::return_value_policy::copy,
  py::arg("g"),
  py::arg("parent") )
;
}
}  // end namespace pwmigpy
}  // end namespace pwmig
