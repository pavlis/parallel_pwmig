#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/operators.h>
#include <pybind11/embed.h>

#include "pwmig/pwmigcore/pwmig.h"
#include "pwmig/pwmigcore/DepthDependentAperture.h"
#include "pwmig/pwmigcore/RectangularSlownessGrid.h"
#include "pwmig/pwmigcore/pwstack.h"
#include "pwmig/pwmigcore/SlownessVectorMatrix.h"
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
using pwmig::pwmigcore::pwstack_ensemble;
using pwmig::pwmigcore::SlownessVectorMatrix;
using pwmig::pwmigcore::Build_GCLraygrid;
using pwmig::pwmigcore::ComputeIncidentWaveRaygrid;
using pwmig::pwmigcore::migrate_one_seismogram;


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
m.def("pwstack_ensemble",&pwstack_ensemble,"Run pwstack algorithm on a SeismogramEnsemble",
   py::return_value_policy::copy,
   py::arg("d"),
   py::arg("ugrid"),
   py::arg("mute"),
   py::arg("stackmute"),
   py::arg("stack_count_cutoff"),
   py::arg("tstart"),
   py::arg("tend"),
   py::arg("aperture"),
   py::arg("aperture_taper_length"),
   py::arg("centroid_cutoff"),
   py::arg("mdlcopy"),
   py::arg("save_history"),
   py::arg("algid")
 );
 py::class_<SlownessVectorMatrix>(m,"SlownessVectorMatrix",
      "Defines a grid parallel to the pseudostation GCLgrid to store slowness vectors at each grid position")
   .def(py::init<const int, const int>(),
       "Create workspace and fill grid with all 0 slowness vectors (allocation and initialization to invalid values)")
   .def(py::init<const SlownessVectorMatrix&>())
   .def("rows",&SlownessVectorMatrix::rows,"Return number of rows (index 1) in the defined grid")
   .def("columns",&SlownessVectorMatrix::columns,"Return number of columns (index 2) in the defined grid")
   .def("set_slowness",&SlownessVectorMatrix::set_slowness,
     "Set slowness vector (arg 0) at grid position i,j (args 1 and 2 )")
   .def("get_slowness",[](SlownessVectorMatrix& self,int i, int j)
     {
       return self(i,j);
     })
  ;
  m.def("Build_GCLraygrid",&Build_GCLraygrid,"Creates a structured grid with lines of constant i,j defined by ray trace geometry",
   py::return_value_policy::copy,
   py::arg("fixed_u"),
   py::arg("parent"),
   py::arg("u"),
   py::arg("svm"),
   py::arg("vmod"),
   py::arg("zmax"),
   py::arg("tmax"),
   py::arg("dt")
 );
 m.def("ComputeIncidentWaveRaygrid",&ComputeIncidentWaveRaygrid,
  "Computes a raygrid for the incident wavefield driven by a grid of slowness values, 1d ray tracing, and 3D model implemented by approximate ray tracing",
  py::return_value_policy::copy,
  py::arg("pstagrid"),
  py::arg("border_pad"),
  py::arg("UP3d"),
  py::arg("vp1d"),
  py::arg("svm"),
  py::arg("zmax"),
  py::arg("tmax"),
  py::arg("dt"),
  py::arg("zdecfac"),
  py::arg("use_3d")
 );
 m.def("migrate_one_seismogram",&migrate_one_seismogram,"pwmig innermost loop function - projects 3C seismogram data along a ray path",
 py::return_value_policy::copy,
 py::arg("pwdata"),
 py::arg("parent"),
 py::arg("raygrid"),
 py::arg("TPgrid"),
 py::arg("Us3d"),
 py::arg("vp1d"),
 py::arg("vs1d"),
 py::arg("control")
);
/* These aren't needed for pwmig python code, but are of broader use so 
I have this binding code for these functions */
 m.def("remove_mean_x3",&pwmig::pwmigcore::remove_mean_x3,"Removes mean for slices in the x3 direction (normally assumed to be a depth variable) from scalar field data",
   py::return_value_policy::copy,
   py::arg("f") 
  );  

}
}  // end namespace pwmigpy
}  // end namespace pwmig
