#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/operators.h>
#include <pybind11/embed.h>

#include <boost/archive/basic_archive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/serialization.hpp>

#include "mspass/utility/Metadata.h"
#include "mspass/utility/AntelopePf.h"
#include "pwmig/gclgrid/gclgrid.h"
#include "pwmig/gclgrid/gclgrid_subs.h"
#include "pwmig/gclgrid/RegionalCoordinates.h"
#include "pwmig/gclgrid/PWMIGfielddata.h"

/* these are needed to make std::vector containers function propertly.
This was borrowed from mspass pybind11 cc files */
PYBIND11_MAKE_OPAQUE(std::vector<double>);

namespace pwmig {
namespace pwmigpy {

namespace py=pybind11;
using namespace std;
using namespace mspass::utility;
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
  mspass::utility::Metadata get_attributes() const override
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
  .def(py::pickle(
     [](const Geographic_point &self) {
       return py::make_tuple(self.lat,self.lon,self.r);
     },
     [](py::tuple t) {
       double lat_in = t[0].cast<double>();
       double lon_in = t[1].cast<double>();
       double r_in = t[2].cast<double>();
       return Geographic_point(lat_in, lon_in, r_in);
     }
   )
  )
  ;
py::class_<pwmig::gclgrid::Cartesian_point>(m,"Cartesian_point","Point on Earth defined coordinates in radians")
    .def(py::init<>())
    .def(py::init<const Cartesian_point&>())
    .def("coordinates",&Cartesian_point::coordinates,"Return vector of coordinates in standard order")
    .def_readwrite("x1",&Cartesian_point::x1,"x1 coordinate axis value (km)")
    .def_readwrite("x2",&Cartesian_point::x2,"x2 coordinate axis value (km)")
    .def_readwrite("x3",&Cartesian_point::x3,"x3 coordinate axis value (km)")
    .def(py::pickle(
       [](const Cartesian_point &self) {
         return py::make_tuple(self.x1,self.x2,self.x3);
       },
       [](py::tuple t) {
         double x1_in = t[0].cast<double>();
         double x2_in = t[1].cast<double>();
         double x3_in = t[2].cast<double>();
         return Cartesian_point(x1_in, x2_in, x3_in);
       }
     )
    )
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
  //.def_readwrite("i0",&BasicGCLgrid::i0)
  //.def_readwrite("j0",&BasicGCLgrid::j0)
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
  .def(py::pickle(
    [](const GCLgrid &self) {
      Metadata md;
      md=self.get_attributes();
      pybind11::object sbuf;
      sbuf=serialize_metadata_py(md);
      /* this incantation was borrowed from mspass pickle section
      of pybind11 code for Seismogram*/
      size_t size_arrays=self.n1*self.n2;
      if(size_arrays == 0)
      {
        // We need this for default constructed grids
        py::array_t<double, py::array::f_style> x1arr(size_arrays,NULL);
        py::array_t<double, py::array::f_style> x2arr(size_arrays,NULL);
        py::array_t<double, py::array::f_style> x3arr(size_arrays,NULL);
        return py::make_tuple(sbuf,size_arrays,x1arr,x2arr,x3arr);
      }
      else
      {
        py::array_t<double, py::array::f_style> x1arr(size_arrays,&(self.x1[0][0]));
        py::array_t<double, py::array::f_style> x2arr(size_arrays,&(self.x2[0][0]));
        py::array_t<double, py::array::f_style> x3arr(size_arrays,&(self.x3[0][0]));
        return py::make_tuple(sbuf,size_arrays,x1arr,x2arr,x3arr);
      }
  },
  [](py::tuple t) {
    pybind11::object sbuf=t[0];
    Metadata md=mspass::utility::restore_serialized_metadata_py(sbuf);
    /* Assume these are defined or we are hosed anyway*/
    int n1,n2;
    n1=md.get_int("n1");
    n2=md.get_int("n2");
    size_t array_size_from_md(n1*n2);
    GCLgrid result(n1,n2);
    /* This template function takes all the parameters from md and
    sets the transformation vector and matrix properly.  It does NOT
    alloc the arrays.  Hence for zero size we can just return */
    pwmig::gclgrid::pfload_common_GCL_attributes<GCLgrid>(result,md);

    if(array_size_from_md > 0)
    {
      size_t size_array = t[1].cast<size_t>();
      if(size_array != array_size_from_md)
        throw mspass::utility::MsPASSError("pickle serialization:  metadata grid size does not match buffer size",
           mspass::utility::ErrorSeverity::Fatal);

      py::array_t<double, py::array::f_style> array_buffer;
      array_buffer=t[2].cast<py::array_t<double, py::array::f_style>>();
      py::buffer_info info = array_buffer.request();
      /* x1 is double **.  x[0] is a confusing but standard way to get the
      base pointer.  In gclgrid the base pointer is the first type of
      a contiguous array block of size_array doubles created above*/
      memcpy(result.x1[0],info.ptr,sizeof(double)*size_array);
      /* now pretty much exactly the same for x2 and x3.
      array_buffer is just a pointer we are reassigning so this should work */
      array_buffer=t[3].cast<py::array_t<double, py::array::f_style>>();
      info = array_buffer.request();
      memcpy(result.x2[0],info.ptr,sizeof(double)*size_array);
      array_buffer=t[4].cast<py::array_t<double, py::array::f_style>>();
      info = array_buffer.request();
      memcpy(result.x3[0],info.ptr,sizeof(double)*size_array);
    }
    return result;
  }
  ))
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
  .def("lookup",&GCLgrid3d::lookup,"Find point by cartesian coordinates (depricated)")
  .def("parallel_lookup",py::overload_cast<const double,const double,const double,int&,int&,int&>
      (&GCLgrid3d::parallel_lookup,py::const_),"Thread save lookup method",py::call_guard<py::gil_scoped_release>())
  .def("parallel_lookup",py::overload_cast<const double,const double,const double,std::vector<int>&>
      (&GCLgrid3d::parallel_lookup,py::const_),"Thread save lookup method - vector index overloading",py::call_guard<py::gil_scoped_release>())
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
  .def("get_lookup_origin",&GCLgrid3d::get_lookup_origin,"Fetch the lookup origin - vector of 3 ints")
  .def("set_lookup_origin",py::overload_cast<const int, const int, const int>
      (&GCLgrid3d::set_lookup_origin),"Set lookup origin explicitly (3 int args)")
  .def("set_lookup_origin",py::overload_cast<>(&GCLgrid3d::set_lookup_origin),
      "Set lookup origin to default")
  .def("point_is_inside_grid",&GCLgrid3d::point_is_inside_grid,
      "Test if a point specified by geographic coordinates is inside this grid")
  .def_readwrite("n3",&GCLgrid3d::n3)
  .def_readwrite("dx3_nom",&GCLgrid3d::dx3_nom)
  /* k0 is public, but we should use the new getters and setters that
  are defined above instead of setting this via this mechanism. */
  //.def_readwrite("k0",&GCLgrid3d::k0)
  .def(py::pickle(
    [](const GCLgrid3d &self) {
      Metadata md;
      md=self.get_attributes();
      pybind11::object sbuf;
      sbuf=serialize_metadata_py(md);
      /* this incantation was borrowed from mspass pickle section
      of pybind11 code for Seismogram*/
      size_t size_arrays=self.n1*self.n2*self.n3;
      if(size_arrays == 0)
      {
        // We need this for default constructed grids
        py::array_t<double, py::array::f_style> x1arr(size_arrays,NULL);
        py::array_t<double, py::array::f_style> x2arr(size_arrays,NULL);
        py::array_t<double, py::array::f_style> x3arr(size_arrays,NULL);
        return py::make_tuple(sbuf,size_arrays,x1arr,x2arr,x3arr);
      }
      else
      {
        py::array_t<double, py::array::f_style> x1arr(size_arrays,&(self.x1[0][0][0]));
        py::array_t<double, py::array::f_style> x2arr(size_arrays,&(self.x2[0][0][0]));
        py::array_t<double, py::array::f_style> x3arr(size_arrays,&(self.x3[0][0][0]));
        return py::make_tuple(sbuf,size_arrays,x1arr,x2arr,x3arr);
      }
  },
  [](py::tuple t) {
    pybind11::object sbuf=t[0];
    Metadata md=mspass::utility::restore_serialized_metadata_py(sbuf);
    /* Assume these are defined or we are hosed anyway*/
    int n1,n2,n3;
    n1=md.get_int("n1");
    n2=md.get_int("n2");
    n3=md.get_int("n3");
    size_t array_size_from_md(n1*n2*n3);
    GCLgrid3d result(n1,n2,n3);
    /* This template function takes all the parameters from md and
    sets the transformation vector and matrix properly.  It does NOT
    alloc the arrays.  Hence for zero size we can just return */
    pwmig::gclgrid::pfload_common_GCL_attributes<GCLgrid3d>(result,md);
    /* We need this additional call for 3d grid - a design flaw in gclgrid*/
    pwmig::gclgrid::pfload_3dgrid_attributes<GCLgrid3d>(result,md);
    /* Default constructed result has null arrays so we do nothing with
    them other than load the metadata. When the sizes are nonzero the
    constructor is assumed to have allocated the space for the arrays
    we copy in below*/
    if(array_size_from_md > 0)
    {
      size_t size_array = t[1].cast<size_t>();
      if(size_array != array_size_from_md)
        throw mspass::utility::MsPASSError("pickle serialization:  metadata grid size does not match buffer size",
           mspass::utility::ErrorSeverity::Fatal);;
      py::array_t<double, py::array::f_style> array_buffer;
      array_buffer=t[2].cast<py::array_t<double, py::array::f_style>>();
      py::buffer_info info = array_buffer.request();
      /* x1 is double ***.  x[0][0] is a confusing but standard way to get the
      base pointer.  In gclgrid the base pointer is the first type of
      a contiguous array block of size_array doubles created above*/
      memcpy(result.x1[0][0],info.ptr,sizeof(double)*size_array);
      /* now pretty much exactly the same for x2 and x3.
      array_buffer is just a pointer we are reassigning so this should work */
      array_buffer=t[3].cast<py::array_t<double, py::array::f_style>>();
      info = array_buffer.request();
      memcpy(result.x2[0][0],info.ptr,sizeof(double)*size_array);
      array_buffer=t[4].cast<py::array_t<double, py::array::f_style>>();
      info = array_buffer.request();
      memcpy(result.x3[0][0],info.ptr,sizeof(double)*size_array);
    }
    return result;
  }
  ))
;
py::class_<GCLscalarfield,GCLgrid>(m,"GCLscalarfield","Two-dimensional grid with scalar attributes at each node")
  .def(py::init<>())
  .def(py::init<const int, const int>())
  /* The order here matters.  The copy constructor MUST appear before the
  constructor that uses the GCLgrid base class.   If not the base class
  constructor is called. The reason is that python resolves this kind of
  overloading trying in the order that is defined by these bindings */
  .def(py::init<const GCLscalarfield&>())
  .def(py::init<const GCLgrid&>())
  .def(py::init<const string, const string, const bool>())
  .def(py::init<const Metadata&>())
  .def("zero",&GCLscalarfield::zero,"Set all field attributes to 0")
  .def("save",&GCLscalarfield::save,"Save contents to a file")
  .def("interpolate",&GCLscalarfield::interpolate,"Interpolate grid to get value at point passed")
  .def("get_attributes",&GCLscalarfield::get_attributes,
     "Fetch all attributes into a Metadata container")
  .def("set_value",&GCLscalarfield::set_value,"Set the value at a specified grid point")
  .def("get_value",&GCLscalarfield::get_value,"Get the value at a specified grid point")
   /* This is normally the right syntax in pybind11 for operator+= but
   not working here for some mysterious reason.  After experimentation it is
   clear the problem is that pybind11 demands operator+= have the signature
     T& operator+=(const T& d);
   That won't work with the design of gclgrid because the right hand side is
   not const.   The iteration starting position is held internally.
   Another method of handling sums is needed - adding an accumulate method
   below
  .def(py::self += py::self)
  */
  .def(py::self *= double())
;

py::class_<GCLvectorfield,GCLgrid>(m,"GCLvectorfield","Two-dimensional grid with vector attributes at each node")
  .def(py::init<>())
  .def(py::init<const int, const int,const int>())
  /* The order here matters.  The copy constructor MUST appear before the
  constructor that uses the GCLgrid base class.   If not the base class
  constructor is called. The reason is that python resolves this kind of
  overloading trying in the order that is defined by these bindings */
  .def(py::init<const GCLvectorfield&>())
  .def(py::init<const GCLgrid&,const int>())
  .def(py::init<const string, const string, const bool>())
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
  //.def(py::self += py::self)
  .def_readwrite("nv",&GCLvectorfield::nv,"Number of components in each vector")
;
// 3D versions of two class definitions above
py::class_<GCLscalarfield3d,GCLgrid3d>(m,"GCLscalarfield3d","Three-dimensional grid with scalar attributes at each node")
  .def(py::init<>())
  .def(py::init<const int, const int, const int>())
  /* The order here matters.  The copy constructor MUST appear before the
  constructor that uses the GCLgrid3d base class.   If not the base class
  constructor is called. The reason is that python resolves this kind of
  overloading trying in the order that is defined by these bindings */
  .def(py::init<const GCLscalarfield3d&>())
  .def(py::init<const GCLgrid3d&>())
  .def(py::init<const string, const string>())
  .def(py::init<const Metadata&>())
  .def("zero",&GCLscalarfield3d::zero,"Set all field attributes to 0")
  .def("save",&GCLscalarfield3d::save,"Save contents to a file")
  .def("interpolate",&GCLscalarfield3d::interpolate,"Interpolate grid to get value at point passed")
  .def("get_attributes",&GCLscalarfield3d::get_attributes,
     "Fetch all attributes into a Metadata container")
  .def("get_value",
    py::overload_cast<const int, const int, const int>(&GCLscalarfield3d::get_value,py::const_),
    "Get the value at a specified grid point")
  .def("get_value",
    py::overload_cast<const Geographic_point>(&GCLscalarfield3d::get_value,py::const_),
    "Get the value at a point specified by geographic coordinates")
  .def("set_value",&GCLscalarfield3d::set_value,"Set the value at a specified grid point")
  .def(py::self *= double())
  .def(py::self += py::self)
  .def(py::pickle(
    [](const GCLscalarfield3d &self) {
      Metadata md;
      md=self.get_attributes();
      pybind11::object sbuf;
      sbuf=serialize_metadata_py(md);
      /* this incantation was borrowed from mspass pickle section
      of pybind11 code for Seismogram*/
      size_t size_arrays=self.n1*self.n2*self.n3;
      if(size_arrays == 0)
      {
        // We need this for default constructed grids
        py::array_t<double, py::array::f_style> x1arr(size_arrays,NULL);
        py::array_t<double, py::array::f_style> x2arr(size_arrays,NULL);
        py::array_t<double, py::array::f_style> x3arr(size_arrays,NULL);
        py::array_t<double, py::array::f_style> valarr(size_arrays,NULL);
        return py::make_tuple(sbuf,size_arrays,x1arr,x2arr,x3arr,valarr);
      }
      else
      {
        py::array_t<double, py::array::f_style> x1arr(size_arrays,&(self.x1[0][0][0]));
        py::array_t<double, py::array::f_style> x2arr(size_arrays,&(self.x2[0][0][0]));
        py::array_t<double, py::array::f_style> x3arr(size_arrays,&(self.x3[0][0][0]));
        py::array_t<double, py::array::f_style> valarr(size_arrays,&(self.val[0][0][0]));
        return py::make_tuple(sbuf,size_arrays,x1arr,x2arr,x3arr,valarr);
      }
  },
  [](py::tuple t) {
    pybind11::object sbuf=t[0];
    Metadata md=mspass::utility::restore_serialized_metadata_py(sbuf);
    /* Assume these are defined or we are hosed anyway*/
    int n1,n2,n3;
    n1=md.get_int("n1");
    n2=md.get_int("n2");
    n3=md.get_int("n3");
    size_t array_size_from_md(n1*n2*n3);
    GCLscalarfield3d result(n1,n2,n3);
    /* This template function takes all the parameters from md and
    sets the transformation vector and matrix properly.  It does NOT
    alloc the arrays.  Hence for zero size we can just return */
    pwmig::gclgrid::pfload_common_GCL_attributes<GCLgrid3d>(result,md);
    /* We need this additional call for 3d grid - a design flaw in gclgrid*/
    pwmig::gclgrid::pfload_3dgrid_attributes<GCLgrid3d>(result,md);
    /* If the size is zero the arrays are assumed set NULL - that would
    be behavior for a default constructed object.  Otherwise the constructor
    we used above will create the data arrays we can push the array data to.*/
    if(array_size_from_md>0)
    {
      size_t size_array = t[1].cast<size_t>();
      if(size_array != array_size_from_md)
        throw mspass::utility::MsPASSError("pickle serialization:  metadata grid size does not match buffer size",
           mspass::utility::ErrorSeverity::Fatal);

      py::array_t<double, py::array::f_style> array_buffer;
      array_buffer=t[2].cast<py::array_t<double, py::array::f_style>>();
      py::buffer_info info = array_buffer.request();
      /* x1 is double **.  x[0][0] is a confusing but standard way to get the
      base pointer.  In gclgrid the base pointer is the first type of
      a contiguous array block of size_array doubles created above*/
      memcpy(result.x1[0][0],info.ptr,sizeof(double)*size_array);
      /* now pretty much exactly the same for x2 and x3.
      array_buffer is just a pointer we are reassigning so this should work */
      array_buffer=t[3].cast<py::array_t<double, py::array::f_style>>();
      info = array_buffer.request();
      memcpy(result.x2[0][0],info.ptr,sizeof(double)*size_array);
      array_buffer=t[4].cast<py::array_t<double, py::array::f_style>>();
      info = array_buffer.request();
      memcpy(result.x3[0][0],info.ptr,sizeof(double)*size_array);
      array_buffer=t[5].cast<py::array_t<double, py::array::f_style>>();
      info = array_buffer.request();
      memcpy(result.val[0][0],info.ptr,sizeof(double)*size_array);
    }
    return result;
  }
  ))
;

py::class_<GCLvectorfield3d,GCLgrid3d>(m,"GCLvectorfield3d","Three-dimensional grid with vector attributes at each node")
  .def(py::init<>())
  .def(py::init<const int, const int,const int,const int>())
  /* The order here matters.  The copy constructor MUST appear before the
  constructor that uses the GCLgrid3d base class.   If not the base class
  constructor is called. The reason is that python resolves this kind of
  overloading trying in the order that is defined by these bindings */
  .def(py::init<const GCLvectorfield3d&>())
  .def(py::init<const GCLgrid3d&,const int>())
  .def(py::init<const string, const string>())
  .def(py::init<const Metadata&>())
  .def("zero",&GCLvectorfield3d::zero,"Set all field attributes to 0")
  .def("save",&GCLvectorfield3d::save,"Save contents to a file")
  .def("interpolate",&GCLvectorfield3d::interpolate,"Interpolate grid to get vector values at point passed",py::call_guard<py::gil_scoped_release>())
  .def("get_attributes",&GCLvectorfield3d::get_attributes,
     "Fetch all attributes into a Metadata container")
  .def("get_value",
    py::overload_cast<const int, const int, const int>(&GCLvectorfield3d::get_value,py::const_),
    "Get the vector of values at a specified grid point")
  .def("get_value",
    py::overload_cast<const Geographic_point>(&GCLvectorfield3d::get_value,py::const_),
    "Get the value at a specified point specified by geographic coordinates")
  .def("set_value",&GCLvectorfield3d::set_value,"Set the value at a specified grid point")
  .def(py::self += py::self)
  .def(py::self + py::self)
  .def(py::self *= double())
  .def_readwrite("nv",&GCLvectorfield3d::nv,"Number of components in each vector")
  .def(py::pickle(
    [](const GCLvectorfield3d &self) {
      Metadata md;
      md=self.get_attributes();
      pybind11::object sbuf;
      sbuf=serialize_metadata_py(md);
      /* this incantation was borrowed from mspass pickle section
      of pybind11 code for Seismogram*/
      size_t size_arrays=self.n1*self.n2*self.n3;
      size_t size_val = size_arrays*self.nv;
      if(size_arrays == 0)
      {
        // We need this for default constructed grids
        py::array_t<double, py::array::f_style> x1arr(size_arrays,NULL);
        py::array_t<double, py::array::f_style> x2arr(size_arrays,NULL);
        py::array_t<double, py::array::f_style> x3arr(size_arrays,NULL);
        py::array_t<double, py::array::f_style> valarr(size_val,NULL);
        return py::make_tuple(sbuf,size_arrays,x1arr,x2arr,x3arr,valarr);
      }
      else
      {
        py::array_t<double, py::array::f_style> x1arr(size_arrays,&(self.x1[0][0][0]));
        py::array_t<double, py::array::f_style> x2arr(size_arrays,&(self.x2[0][0][0]));
        py::array_t<double, py::array::f_style> x3arr(size_arrays,&(self.x3[0][0][0]));
        py::array_t<double, py::array::f_style> valarr(size_val,&(self.val[0][0][0][0]));
        return py::make_tuple(sbuf,size_arrays,x1arr,x2arr,x3arr,valarr);
      }
  },
  [](py::tuple t) {
    pybind11::object sbuf=t[0];
    Metadata md=mspass::utility::restore_serialized_metadata_py(sbuf);
    /* Assume these are defined or we are hosed anyway*/
    int n1,n2,n3,nv;
    n1=md.get_int("n1");
    n2=md.get_int("n2");
    n3=md.get_int("n3");
    nv=md.get_int("nv");
    size_t array_size_from_md(n1*n2*n3);
    GCLvectorfield3d result(n1,n2,n3,nv);
    /* This template function takes all the parameters from md and
    sets the transformation vector and matrix properly.  It does NOT
    alloc the arrays.  Hence for zero size we can just return */
    pwmig::gclgrid::pfload_common_GCL_attributes<GCLgrid3d>(result,md);
    /* We need this additional call for 3d grid - a design flaw in gclgrid*/
    pwmig::gclgrid::pfload_3dgrid_attributes<GCLgrid3d>(result,md);
    /* I think we need to set this one specially - the above may not set it */
    result.nv=nv;
    /* If the size is zero the arrays are assumed set NULL - that would
    be behavior for a default constructed object.  Otherwise the constructor
    we used above will create the data arrays we can push the array data to.*/
    if(array_size_from_md>0)
    {
      size_t size_array = t[1].cast<size_t>();
      if(size_array != array_size_from_md)
        throw mspass::utility::MsPASSError("pickle serialization:  metadata grid size does not match buffer size",
           mspass::utility::ErrorSeverity::Fatal);

      py::array_t<double, py::array::f_style> array_buffer;
      array_buffer=t[2].cast<py::array_t<double, py::array::f_style>>();
      py::buffer_info info = array_buffer.request();
      /* x1 is double **.  x[0][0] is a confusing but standard way to get the
      base pointer.  In gclgrid the base pointer is the first type of
      a contiguous array block of size_array doubles created above*/
      memcpy(result.x1[0][0],info.ptr,sizeof(double)*size_array);
      /* now pretty much exactly the same for x2 and x3.
      array_buffer is just a pointer we are reassigning so this should work */
      array_buffer=t[3].cast<py::array_t<double, py::array::f_style>>();
      info = array_buffer.request();
      memcpy(result.x2[0][0],info.ptr,sizeof(double)*size_array);
      array_buffer=t[4].cast<py::array_t<double, py::array::f_style>>();
      info = array_buffer.request();
      memcpy(result.x3[0][0],info.ptr,sizeof(double)*size_array);
      /* We have to compute a new buffer size for the val array because it
      has this nv multiplier.  We also have another level of pointer*/
      size_t size_val=size_array*nv;
      array_buffer=t[5].cast<py::array_t<double, py::array::f_style>>();
      info = array_buffer.request();
      memcpy(result.val[0][0][0],info.ptr,sizeof(double)*size_val);
    }
    return result;
  }
  ))
;
py::class_<PWMIGfielddata,GCLvectorfield3d>(m,"PWMIGfielddata",
    "Extension of generic GCLvectorfield3d for pwmig output - mostly add error logger")
  .def(py::init<>())
  .def(py::init<const pwmig::gclgrid::GCLgrid3d>())
  .def(py::init<const pwmig::gclgrid::PWMIGfielddata>())
  .def("accumulate",&PWMIGfielddata::accumulate,"Sum a single seismogram like migrated output to field",py::call_guard<py::gil_scoped_release>())
  .def_readwrite("elog",&PWMIGfielddata::elog,"ErrorLogger object for handling errors in mspass way")
  /* This is annoyingly parallel to GCLvectorfield3d but I don't see how to do it
  with pybind11 to use inheritance unless we did a double pickle which would be a
  very slow thing to do for this (usually large) object. */
  .def(py::pickle(
    [](const PWMIGfielddata &self) {
      Metadata md;
      md=self.get_attributes();
      pybind11::object sbuf;
      sbuf=serialize_metadata_py(md);
      /* this incantation was borrowed from mspass pickle section
      of pybind11 code for Seismogram*/
      size_t size_arrays=self.n1*self.n2*self.n3;
      size_t size_val = size_arrays*self.nv;
      /* Added for elog part - only difference really from GCLvectorfield3d*/
      stringstream ss;
      boost::archive::text_oarchive ar(ss);
      ar << self.elog;
      string serialized_elog(ss.str());
      if(size_arrays == 0)
      {
        // We need this for default constructed grids
        py::array_t<double, py::array::f_style> x1arr(size_arrays,NULL);
        py::array_t<double, py::array::f_style> x2arr(size_arrays,NULL);
        py::array_t<double, py::array::f_style> x3arr(size_arrays,NULL);
        py::array_t<double, py::array::f_style> valarr(size_val,NULL);
        return py::make_tuple(sbuf,size_arrays,x1arr,x2arr,x3arr,valarr,serialized_elog);
      }
      else
      {
        py::array_t<double, py::array::f_style> x1arr(size_arrays,&(self.x1[0][0][0]));
        py::array_t<double, py::array::f_style> x2arr(size_arrays,&(self.x2[0][0][0]));
        py::array_t<double, py::array::f_style> x3arr(size_arrays,&(self.x3[0][0][0]));
        py::array_t<double, py::array::f_style> valarr(size_val,&(self.val[0][0][0][0]));
        return py::make_tuple(sbuf,size_arrays,x1arr,x2arr,x3arr,valarr,serialized_elog);
      }
  },
  [](py::tuple t) {
    /* this function differs a fair bit from the GCLvectorfield3d version
    due to the need to handle elog */
    pybind11::object sbuf=t[0];
    Metadata md=mspass::utility::restore_serialized_metadata_py(sbuf);
    /* Assume these are defined or we are hosed anyway*/
    int n1,n2,n3,nv;
    n1=md.get_int("n1");
    n2=md.get_int("n2");
    n3=md.get_int("n3");
    nv=md.get_int("nv");
    size_t array_size_from_md(n1*n2*n3);
    /* We have to reconstruct elog first as it is used in the specialized
    constructor called immediately after it is deserialized. */
    ErrorLogger elog_to_clone;
    stringstream ss(t[6].cast<std::string>());
    boost::archive::text_iarchive ar(ss);
    ar>>elog_to_clone;
    PWMIGfielddata result(n1,n2,n3,elog_to_clone);
    /* This template function takes all the parameters from md and
    sets the transformation vector and matrix properly.  It does NOT
    alloc the arrays.  Hence for zero size we can just return */
    pwmig::gclgrid::pfload_common_GCL_attributes<GCLgrid3d>(result,md);
    /* We need this additional call for 3d grid - a design flaw in gclgrid*/
    pwmig::gclgrid::pfload_3dgrid_attributes<GCLgrid3d>(result,md);
    /* I think we need to set this one specially - the above may not set it */
    result.nv=nv;
    /* If the size is zero the arrays are assumed set NULL - that would
    be behavior for a default constructed object.  Otherwise the constructor
    we used above will create the data arrays we can push the array data to.*/
    if(array_size_from_md>0)
    {
      size_t size_array = t[1].cast<size_t>();
      if(size_array != array_size_from_md)
        throw mspass::utility::MsPASSError("pickle serialization:  metadata grid size does not match buffer size",
           mspass::utility::ErrorSeverity::Fatal);

      py::array_t<double, py::array::f_style> array_buffer;
      array_buffer=t[2].cast<py::array_t<double, py::array::f_style>>();
      py::buffer_info info = array_buffer.request();
      /* x1 is double **.  x[0][0] is a confusing but standard way to get the
      base pointer.  In gclgrid the base pointer is the first type of
      a contiguous array block of size_array doubles created above*/
      memcpy(result.x1[0][0],info.ptr,sizeof(double)*size_array);
      /* now pretty much exactly the same for x2 and x3.
      array_buffer is just a pointer we are reassigning so this should work */
      array_buffer=t[3].cast<py::array_t<double, py::array::f_style>>();
      info = array_buffer.request();
      memcpy(result.x2[0][0],info.ptr,sizeof(double)*size_array);
      array_buffer=t[4].cast<py::array_t<double, py::array::f_style>>();
      info = array_buffer.request();
      memcpy(result.x3[0][0],info.ptr,sizeof(double)*size_array);
      /* We have to compute a new buffer size for the val array because it
      has this nv multiplier.  We also have another level of pointer*/
      size_t size_val=size_array*nv;
      array_buffer=t[5].cast<py::array_t<double, py::array::f_style>>();
      info = array_buffer.request();
      memcpy(result.val[0][0][0],info.ptr,sizeof(double)*size_val);
    }
    return result;
  }
  ))
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
m.def("extract_component",[](const GCLvectorfield3d& f,const int icomp) {
  return (std::unique_ptr<GCLscalarfield3d>
     (pwmig::gclgrid::extract_component(f,icomp)));
}
)
;

}
}  // end namespace pwmigpy
}  // end namespace pwmig
