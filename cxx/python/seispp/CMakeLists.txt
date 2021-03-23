include_directories(
  ${Boost_INCLUDE_DIRS}
  ${GSL_INCLUDE_DIRS}
  ${pylib11_INCLUDE_DIRS}
  ${PYTHON_INCLUDE_DIRS}
  $ENV{MSPASS_HOME}/include
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_BINARY_DIR}/include)


pybind11_add_module(pyseispp seispp_py.cc)

target_link_libraries(pyseispp PRIVATE seispp dsap gclgrid $ENV{MSPASS_HOME}/lib/libutility.a ${Boost_LIBRARIES})
set_target_properties(pyseispp PROPERTIES OUTPUT_NAME "seispp")
if(DEFINED CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  set_target_properties(pyseispp PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/ccore)
endif()
