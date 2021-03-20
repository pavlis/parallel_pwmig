#include <pybind11/pybind11.h>
//#include <pybind11/numpy.h>
//#include <pybind11/stl.h>
//#include <pybind11/stl_bind.h>
#include <pybind11/operators.h>
//#include <pybind11/embed.h>

#include "pwmig/seispp/Hypocenter.h"
#include "pwmig/seispp/EventCatalog.h"
#include "mspass/utility/Metadata.h"

namespace pwmig {
namespace pwmigpy {

namespace py=pybind11;
using namespace std;
using mspass::utility::Metadata;
using namespace pwmig::seispp;

PYBIND11_MODULE(seispp, m) {
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
py::class_<EventCatalog>(m,"EventCatalog","In memory earthquake source data manager")
  .def(py::init<>())
  .def(py::init<mspass::utility::MetadataList>())
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
  /* this won't work because python doesn't have the ++ or -- operators
  If needed will need to implement operator+
  .def(py::self ++ py::self)
  */
  .def("advance",&EventCatalog::advance,"Advance the current event pointer by number of slots requested")
;
}
}  // end namespace pwmigpy
}  // end namespace pwmig
