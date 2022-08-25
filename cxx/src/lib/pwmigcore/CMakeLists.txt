FILE(GLOB sources_pwmigcore *.cc)
include_directories(
  ${Boost_INCLUDE_DIRS}
  ${pybind11_INCLUDE_DIR}
  ${PYTHON_INCLUDE_DIRS}
  ${MSPASSHOME}/include
  ${PROJECT_BINARY_DIR}/include
  ${PROJECT_SOURCE_DIR}/include)

add_library(pwmigcore OBJECT ${sources_pwmigcore})
#target_link_libraries(pwmigcore PRIVATE dsap seispp gclgrid ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} )
#target_link_libraries(pwmigcore PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${MSPASS_LIBRARIES})

#install (TARGETS pwmigcore DESTINATION lib)
