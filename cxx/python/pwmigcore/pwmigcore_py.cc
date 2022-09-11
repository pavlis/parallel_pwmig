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
#include "pwmig/pwmigcore/PWMIGmigrated_seismogram.h"
#include "mspass/utility/Metadata.h"
#include "mspass/utility/AntelopePf.h"
#include "mspass/seismic/Seismogram.h"
#include "mspass/seismic/Ensemble.h"


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
using pwmig::pwmigcore::PWMIGmigrated_seismogram;
using pwmig::pwmigcore::migrate_component;

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
     .def(py::pickle(
         [](const SlownessVectorMatrix &self){
           stringstream ss;
           boost::archive::text_oarchive ar(ss);
           ar << self;
           return py::make_tuple(ss.str());
       },
       [](py::tuple t){
           stringstream ss(t[0].cast<std::string>());
           boost::archive::text_iarchive ar(ss);
           SlownessVectorMatrix svm;
           ar>>svm;
           return svm;
       }
     ))
  ;
  py::class_<PWMIGmigrated_seismogram>(m,"PWMIGmigrated_seismogram","Holds time to depth mapped data internal to pwmig.  A bit like Seismogram but not a child")
    .def(py::init<>())
    .def(py::init<const int, const int, const int>(),"Space allocating constructor - ix1, ix2, n (size)")
    .def(py::init<const PWMIGmigrated_seismogram&>(),"copy constructor")
    .def("copy_elog",&PWMIGmigrated_seismogram::copy_elog,"Copy ErrorLogger data from the Seismogram from which this datum was created")
    /* Attributes */
    .def_readwrite("ix1",&PWMIGmigrated_seismogram::ix1,"Surface GCLgrid index location for x1 position")
    .def_readwrite("ix2",&PWMIGmigrated_seismogram::ix1,"Surface GCLgrid index position for x2 position")
    .def_readwrite("live",&PWMIGmigrated_seismogram::live,"Boolean set true of data are valid - do not use if false")
    .def_readwrite("migrated_data",&PWMIGmigrated_seismogram::migrated_data,"dmatrix containing time to depth mapped data")
    .def_readwrite("domega",&PWMIGmigrated_seismogram::domega,"vector of solid angle terms used in pwmig for amplitude normalization using GRT")
    .def_readwrite("dweight",&PWMIGmigrated_seismogram::dweight,"vector of GRT weight terms san solid angle increment")
    .def_readwrite("elog",&PWMIGmigrated_seismogram::elog,"ErrorLogger used to hold errors transferred from Seismogram and any additional logged errors")

    .def(py::pickle(
          [](const PWMIGmigrated_seismogram &self) {
            /* We always have to serialize the error log - for dead data it
            should contain critical log messages */
            stringstream ss_elog;
            boost::archive::text_oarchive ar(ss_elog);
            ar << self.elog;
            /* This handles datum marked dead or empty equally */
            if( self.live && (self.domega.size()>0))
            {
              /* This is borrowed from Seismogram which also contains a dmatrix */
              size_t d_size = self.migrated_data.rows()*self.migrated_data.columns();
              py::array_t<double, py::array::f_style> mdarr(d_size,
                   self.migrated_data.get_address(0,0));
              py::array_t<double, py::array::f_style> domarr(self.domega.size(),
                   &(self.domega[0]));
              py::array_t<double, py::array::f_style> dwarr(self.dweight.size(),
                   &(self.dweight[0]));
              /* We assume domega and dweight have the same size and the size = N
              migrated_data is 3*N.  We store N as tuple component 4 */
              return py::make_tuple(self.ix1,self.ix2,self.live,ss_elog.str(),
                  self.domega.size(),mdarr,domarr,dwarr);
            }
            else
            {
              /* We force component 2 to false to handle case of
              live being true but data vectors empty - should be dead then
              anyway */
              return py::make_tuple(self.ix1,self.ix2,false,ss_elog.str(),NULL,NULL,NULL,NULL);
            }
          },
          [](py::tuple t) {
            int ix1,ix2;
            ix1 = t[0].cast<int>();
            ix2 = t[1].cast<int>();
            bool live=t[2].cast<bool>();
            mspass::utility::ErrorLogger elog;
            stringstream ss_elog(t[3].cast<std::string>());
            boost::archive::text_iarchive arelog(ss_elog);
            arelog >> elog;
            int data_size = t[4].cast<size_t>();
            /* live or dead both use the space allocating constructor that
            sets ix1 and ix2.  When zero depends upon std:;vector property that
            reserve(0) is handled correctly*/
            PWMIGmigrated_seismogram result(ix1,ix2,data_size);
            if(live)
            {
              result.live=true;
              /* The constructor above allocated space for the dmatrix
              holding migrated_data.  Hence, we can use memcpy to copy
              the internal buffer to the object's copy of same */
              size_t matrix_size=3*data_size;
              py::array_t<double, py::array::f_style> darr;
              darr = t[5].cast<py::array_t<double, py::array::f_style>>();
              py::buffer_info info = darr.request();
              /* Perhaps should verify size matches what we would get from
              buffer_info but for planned use we skip that safety and just
              copy */
              memcpy(result.migrated_data.get_address(0,0),
                  info.ptr,sizeof(double)*matrix_size);
              /* Now code similar to TimeSeries to copy the two std::vector
              containers using buffer protocol.  We reuse darr declared above*/
              darr = t[6].cast<py::array_t<double, py::array::f_style>>();
              info = darr.request();
              /* The constructor calls reserve for the two std::vectors.
              In timeSeries we use a resize but here we use a
              more standard loop to copy the vectors from the buffer
              using push_back.   Because reserve was called we can be sure
              there will be no realloc overhead */
              double *ptr=(double *)info.ptr;
              for(auto i=0;i<data_size;++i,++ptr) result.domega.push_back(*ptr);
              /* Same for dweight*/
              darr = t[7].cast<py::array_t<double, py::array::f_style>>();
              info = darr.request();
              ptr=(double *)info.ptr;
              for(auto i=0;i<data_size;++i,++ptr) result.dweight.push_back(*ptr);
            }
            else
            {
              result.elog = elog;
              // Make sure it is marked dead - redundant but clearer
              result.live = false;
            }
            return result;
         }
         ))

  ;

  /* Note the return_value_policy here is ESSENTIAl to avoid memory leaks.
  We normally used copy but that is a bad idea because this function returns
  a pointer to a GCLscalargrid3d object. */
  m.def("Build_GCLraygrid",&Build_GCLraygrid,"Creates a structured grid with lines of constant i,j defined by ray trace geometry",
   py::return_value_policy::take_ownership,
   py::arg("parent"),
   py::arg("svm"),
   py::arg("vmod"),
   py::arg("zmax"),
   py::arg("tmax"),
   py::arg("dt")
 );
 m.def("ComputeIncidentWaveRaygrid",&ComputeIncidentWaveRaygrid,
  "Computes a raygrid for the incident wavefield driven by a grid of slowness values, 1d ray tracing, and 3D model implemented by approximate ray tracing",
  py::return_value_policy::take_ownership,
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
 //m.def("migrate_one_seismogram",&migrate_one_seismogram,"pwmig innermost loop function - projects 3C seismogram data along a ray path",
 m.def("migrate_one_seismogram",static_cast<PWMIGmigrated_seismogram (*)
   (mspass::seismic::Seismogram& pwdata,
       pwmig::gclgrid::GCLgrid&,
             pwmig::gclgrid::GCLscalarfield3d&,
                     pwmig::gclgrid::GCLscalarfield3d&,
                                pwmig::gclgrid::GCLscalarfield3d&,
                                             pwmig::seispp::VelocityModel_1d&,
                                                            pwmig::seispp::VelocityModel_1d&,
                                                                             mspass::utility::Metadata&)>(&migrate_one_seismogram),
   "pwmig innermost loop function - projects 3C seismogram data along a ray path",
 py::return_value_policy::copy,
 py::arg("pwdata"),
 py::arg("parent"),
 py::arg("raygrid"),
 py::arg("TPgrid"),
 py::arg("Us3d"),
 py::arg("Vp1d"),
 py::arg("Vs1d"),
 py::arg("control")
);
/* These aren't needed for pwmig python code, but are of broader use so
I have this binding code for these functions */
 m.def("remove_mean_x3",&pwmig::pwmigcore::remove_mean_x3,"Removes mean for slices in the x3 direction (normally assumed to be a depth variable) from scalar field data",
   py::return_value_policy::copy,
   py::arg("f")
  );

//m.def("migrate_component",&migrate_component,"Migrate one plane wave component with multithreading defined by Metadata argument",

m.def("migrate_component",static_cast<pwmig::gclgrid::PWMIGfielddata (*)
  (mspass::seismic::Ensemble<mspass::seismic::Seismogram>& pwdata,
        pwmig::gclgrid::GCLgrid&,
          pwmig::gclgrid::GCLscalarfield3d&,
            pwmig::pwmigcore::SlownessVectorMatrix&,
              pwmig::gclgrid::GCLscalarfield3d&,
                pwmig::seispp::VelocityModel_1d&,
                  pwmig::seispp::VelocityModel_1d&,
                    mspass::utility::Metadata&)>(&migrate_component),
    "Migrate one plane wave component with multithreading defined by Metadata argument",
  py::return_value_policy::copy,
  py::arg("pwensemble"),
  py::arg("parent"),
  py::arg("TPgrid"),
  py::arg("VPsvm"),
  py::arg("Us3d"),
  py::arg("Vp1d"),
  py::arg("Vs1d"),
  py::arg("control")
);

}
}  // end namespace pwmigpy
}  // end namespace pwmig
