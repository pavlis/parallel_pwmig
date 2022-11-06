# pybind11_INCLUDE and PYTHON_INCLUDE needed because Metadata.h requires it - not a good thing
include_directories(
  ${pybind11_INCLUDE_DIR}
  ${PYTHON_INCLUDE_DIRS}
  ${MSPASSHOME}/include
  #/home/pavlis/mspass/include
  ${PROJECT_BINARY_DIR}/include
  ${PROJECT_SOURCE_DIR}/include)

add_executable (GCL2Dtovtk GCL2Dtovtk.cc vtk_output_GCLgrid.cc)
  target_link_libraries(GCL2Dtovtk PRIVATE pwmig ${MSPASS_LIBRARIES} ${pybind11_LIBRARIES} ${Boost_LIBRARIES})
