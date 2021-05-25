#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/operators.h>
#include <pybind11/embed.h>

#include "pwmig/pwmigcore/DepthDependentAperture.h"
#include "pwmig/pwmigcore/RectangularSlownessGrid.h"
#include "mspass/utility/Metadata.h"
#include "mspass/utility/AntelopePf.h"


namespace pwmig {
namespace pwmigpy {

namespace py=pybind11;
using namespace std;
using mspass::utility::Metadata;
using mspass::utility::AntelopePf;
using pwmig::pwmigcore::DepthDependentAperture;
using pwmig::pwmigcore::RectangularSlownessGrid;


PYBIND11_MODULE(pwmigcore, m) {
py::class_<DepthDependentAperture>(m,"DepthDependentAperture",
      "Used in plane wave decomposition to allow variable aperture with lag")
  .def(py::init<>())
  .def(py::init<const int>(),"Allocate but do not set number of points passed as arg1")
  .def(py::init<const double, const double,const double, const double,
        const int, const double, const bool>(),"Create with fresnel zone formula - see C++ api documentation")
  .def(py::init<const AntelopePf&,const string>(),"Construct from pf file with tag name given")
  .def(py::init<const DepthDependentAperture&>())
  .def("get_aperture",&DepthDependentAperture::get_aperture,"Return aperture at specified lag time")
  .def("get_cutoff",&DepthDependentAperture::get_cutoff,"Return cutoff distance a specified lag time")
  .def("maximum_cutoff",&DepthDependentAperture::maximum_cutoff,
      "Return largest cutoff distance - used for MongoDB spatial queries")
  .def("maximum_aperture",&DepthDependentAperture::maximum_aperture,
          "Return largest tabulated aperture - used for MongoDB spatial queries")
  .def(py::pickle(
    [](const DepthDependentAperture &self){
      stringstream ss;
      boost::archive::text_oarchive ar(ss);
      ar << self;
      return py::make_tuple(ss.str());
    },
    [](py::tuple t){
      stringstream ss(t[0].cast<std::string>());
      boost::archive::text_iarchive ar(ss);
      DepthDependentAperture aperture;
      ar>>aperture;
      return aperture;
    }
  ))
;
py::class_<RectangularSlownessGrid>(m,"RectangularSlownessGrid",
     "Defines a uniform grid in slowness space for pwstack processing")
  .def(py::init<>())
  .def(py::init<const string, const double, const double, const double, const double,
    const int, const int>(),
      "Fully parameterized constructor - see C++ api doxygen reference")
  .def(py::init<const AntelopePf&,const string>(),
    "Construct from a named branch (Arr block) in an AntelopePf")
  .def(py::init<const RectangularSlownessGrid&>())
  .def("ux",&RectangularSlownessGrid::ux,"Return ux component by grid index")
  .def("uy",&RectangularSlownessGrid::uy,"Return uy component by grid index")
  .def("slow",&RectangularSlownessGrid::slow,
    "Return a slowness vector for point specified by two grid indices")
  .def(py::pickle(
      [](const RectangularSlownessGrid &self){
        stringstream ss;
        boost::archive::text_oarchive ar(ss);
        ar << self;
        return py::make_tuple(ss.str());
    },
    [](py::tuple t){
        stringstream ss(t[0].cast<std::string>());
        boost::archive::text_iarchive ar(ss);
        RectangularSlownessGrid aperture;
        ar>>aperture;
        return aperture;
    }
  ))
;
}
}  // end namespace pwmigpy
}  // end namespace pwmig
