FILE(GLOB sources_dsap *.c)
include_directories(
  ${Boost_INCLUDE_DIRS}
  ${pybind11_INCLUDE_DIR}
  ${PYTHON_INCLUDE_DIRS}
  ${MSPASSHOME}/include
  ${PROJECT_BINARY_DIR}/include
  ${PROJECT_SOURCE_DIR}/include)

add_library(dsap STATIC ${sources_dsap})
target_link_libraries(dsap PRIVATE utility ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

install (TARGETS dsap DESTINATION lib)
