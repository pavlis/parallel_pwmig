FILE(GLOB sources_gclgrid *.cc *.f)
include_directories(
  ${Boost_INCLUDE_DIRS}
  ${pybind11_INCLUDE_DIR}
  ${PYTHON_INCLUDE_DIRS}
  ${MSPASSHOME}/include
  ${PROJECT_BINARY_DIR}/include
  ${PROJECT_SOURCE_DIR}/include)

add_library(gclgrid OBJECT ${sources_gclgrid})
#target_link_libraries(gclgrid PRIVATE dsap ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
#target_link_libraries(gclgrid PRIVATE dsap ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${MSPASS_LIBRARIES})

#install (TARGETS gclgrid DESTINATION lib)
