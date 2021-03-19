FILE(GLOB sources_pwmigcore *.cc)
include_directories(
  ${Boost_INCLUDE_DIRS}
  ${pybind11_INCLUDE_DIR}
  ${PYTHON_INCLUDE_DIRS}
  ${MSPASSHOME}/include
  ${PROJECT_BINARY_DIR}/include
  ${PROJECT_SOURCE_DIR}/include)

add_library(pwmigcore STATIC ${sources_pwmigcore})
target_link_libraries(pwmigcore PRIVATE pwmigcore ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

install (TARGETS pwmigcore DESTINATION lib)
