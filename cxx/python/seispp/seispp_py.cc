#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/operators.h>
#include <pybind11/embed.h>

#include "pwmig/seispp/Hypocenter.h"
#include "pwmig/seispp/RadialGrid.h"
#include "pwmig/seispp/EventCatalog.h"
#include "pwmig/seispp/Stack.h"
#include "pwmig/seispp/VelocityModel_1d.h"
#include "mspass/utility/Metadata.h"
#include "mspass/utility/AntelopePf.h"

/* these are needed to make std::vector containers function propertly.
This was borrowed from mspass pybind11 cc files */
PYBIND11_MAKE_OPAQUE(std::vector<double>);
/* Commented out for reason noted below */
//PYBIND11_MAKE_OPAQUE(std::vector<mspass::seismic::TimeSeries>);

namespace pwmig {
namespace pwmigpy {

namespace py=pybind11;
using namespace std;
using mspass::utility::Metadata;
using mspass::utility::AntelopePf;
using mspass::seismic::TimeSeries;
using mspass::seismic::TimeSeriesEnsemble;
using mspass::algorithms::TimeWindow;
//using namespace pwmig::seispp;
using pwmig::seispp::Hypocenter;
using pwmig::seispp::RadialGrid;
using pwmig::seispp::EventCatalog;
using pwmig::seispp::SectorTest;
using pwmig::seispp::Stack;
using pwmig::seispp::VelocityModel_1d;


PYBIND11_MODULE(seispp, m) {
/* more std::vector requirements borrowed from mspass */
/* Note I had to comment out the bind_vector line below for TimeSeriesVector.  The 
reason seems to be that that symbol is registered global in MsPASS.  DoubleVector, on 
the other hand, uses the default which the documentation for pybind11 says means it 
will be treated as local to only this module.   Hence we have to retain DoubleVector 
but not TimeSeriesVector*/
py::bind_vector<std::vector<double>>(m, "DoubleVector");
//py::bind_vector<std::vector<TimeSeries>>(m, "TimeSeriesVector");


py::class_<Hypocenter>(m,"Hypocenter","Class to hold earthquake space-time coordinates")
  .def(py::init<>())
  .def(py::init<const mspass::utility::Metadata&>())
  .def(py::init<const double, const double, const double, const double>())
  .def("distance",&Hypocenter::distance,
     "Compute epicentral distance in degrees from this source to a point on the earth")
  .def("seaz",&Hypocenter::seaz,
      "Compute station to event azimuth of great circle path from location passed to this source")
  .def("esaz",&Hypocenter::esaz,
      "Compute event to station azimuth of great circle path from this source to specified point")
  .def_readwrite("lat",&Hypocenter::lat,"Source latitude (in radians)")
  .def_readwrite("lon",&Hypocenter::lon,"Source longitude (in radians)")
  .def_readwrite("depth",&Hypocenter::z,"Source depth (in km)")
  .def_readwrite("time",&Hypocenter::time,"Source origin time (Unix epoch time)")
;
py::class_<RadialGrid>(m,"RadialGrid","Defines global radial grid with great circle paths and gcp distances")
  .def(py::init<>())
  .def(py::init<AntelopePf&>())
  .def(py::init<const double, const double, const int,
		const double, const double, const int,const double, const double>())
  .def("lat",&RadialGrid::lat,"Return latitude of lower left corner by grid index value (ir,id)")
  .def("lon",&RadialGrid::lon,"Return longitude of lower left corner by grid index value (ir,id)")
  .def("number_azimuth_bins",&RadialGrid::number_azimuth_bins,"Return number of azimuth bins in grid (one less than points)")
  .def("number_distance_bins",&RadialGrid::number_distance_bins,"Return number of distance bins in grid (one less than points)")
  .def("cell",&RadialGrid::cell,"Return a Metadata container with attributes defining a cell specified by index positions")
  .def_readonly("naz",&RadialGrid::naz,"Number of points in azimuth axis")
  .def_readonly("ndelta",&RadialGrid::ndelta,"Number of points in axis")
  ;
py::class_<SectorTest>(m,"SectorTest","Used for radial grid subsetting - used only to pass to C++ code")
  .def(py::init<RadialGrid&,const int, const int>())
  ;
py::class_<EventCatalog>(m,"EventCatalog","In memory earthquake source data manager")
  .def(py::init<>())
  .def("range",&EventCatalog::range,"Get time span of catalog")
  .def("add",py::overload_cast<const Hypocenter&,const Metadata&>(&EventCatalog::add),
        "Add an event to the catalog already packaged as a hypocenter object (lat lon already converted to radians)" )
  .def("add",py::overload_cast<const Metadata&>(&EventCatalog::replace),
        "Add an event to the catalog using data stored in a Metadata container (angles in deg)" )
  .def("replace",py::overload_cast<const Hypocenter&,const Metadata&>(&EventCatalog::add),
              "Replace an event in the catalog already packaged as a hypocenter object (lat lon already converted to radians)" )
  .def("replace",py::overload_cast<const Metadata&>(&EventCatalog::replace),
              "Replace an event in the catalog using data stored in a Metadata container (angles in deg)" )
  .def("current_aux",&EventCatalog::current_aux,"Return Metadata component of event marked current")
  .def("current",&EventCatalog::current,"Return Hpocenter defined as current (after a find or in iteration)")
  .def("delete_current",&EventCatalog::delete_current,"Delete the event defined as current")
  .def("size",&EventCatalog::size,"Return number of events in the catalog")
  .def("rewind",&EventCatalog::rewind,"Reset current event pointer to top of container")
  /* this won't work because python doesn't have the ++ or -- operators
  If needed will need to implement operator+
  .def(py::self ++ py::self)
  */
  .def("advance",&EventCatalog::advance,"Advance the current event pointer by number of slots requested")
  .def("sector_subset",&EventCatalog::subset<SectorTest>,"Subset using a radial grid")
;
py::class_<Stack>(m,"RobustStack","Ensemble stacker with simple and multiple robust stacking methods")
  .def(py::init<>())
  .def(py::init<TimeSeriesEnsemble&,const TimeWindow>())
  .def(py::init<TimeSeriesEnsemble&,const TimeWindow,const TimeWindow,
      pwmig::seispp::StackType,double>())
  .def_readwrite("stack",&Stack::stack,"Stack produced on construction form input ensemble")
  .def_readwrite("sumwt",&Stack::sumwt,"sum of weights (usuall < fold with robust stacking)")
  .def_readwrite("fold",&Stack::fold,"Number of live TimeSeries in computed stack")
;
py::class_<VelocityModel_1d>(m,"VelocityModel_1d","Data object to hold a layered earth model")
  .def(py::init<>())
  .def(py::init<const std::string,const std::string,const std::string>())
  .def(py::init<const VelocityModel_1d&>())
  .def("getv",&VelocityModel_1d::getv,"Get the velocity at specified depth")
  .def_readwrite("nlayers",&VelocityModel_1d::nlayers,"Number of grid points defining this model")
  .def_readwrite("z",&VelocityModel_1d::z,"Vector of model node depth points (km)")
  .def_readwrite("v",&VelocityModel_1d::v,"Vector of model velocity values (km/s)")
  .def_readwrite("grad",&VelocityModel_1d::grad,"Vector of velocity gradients")
  .def(py::pickle(
    [](const VelocityModel_1d &self) {
      stringstream ss;
      boost::archive::text_oarchive ar(ss);
      ar << self;
      return py::make_tuple(ss.str());
    },
    [](py::tuple t){
      stringstream ss(t[0].cast<std::string>());
      boost::archive::text_iarchive ar(ss);
      VelocityModel_1d vmod;
      ar>>vmod;
      return vmod;
    }
  ))
;
}
}  // end namespace pwmigpy
}  // end namespace pwmig
