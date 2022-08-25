include_directories(
  ${Boost_INCLUDE_DIRS}
  ${GSL_INCLUDE_DIRS}
  ${pylib11_INCLUDE_DIRS}
  ${PYTHON_INCLUDE_DIRS}
  $ENV{MSPASS_HOME}/include
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_BINARY_DIR}/include)


pybind11_add_module(pypwmigcore pwmigcore_py.cc)

target_link_libraries(pypwmigcore PRIVATE pwmig 
  $ENV{MSPASS_HOME}/lib/libmspass.a 
  ${Boost_LIBRARIES})
set_target_properties(pypwmigcore PROPERTIES OUTPUT_NAME "pwmigcore")
if(DEFINED CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  set_target_properties(pypwmigcore PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/ccore)
endif()
