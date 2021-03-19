FILE(GLOB sources_seispp *.cc)
include_directories(
  ${Boost_INCLUDE_DIRS}
  ${pybind11_INCLUDE_DIR}
  ${PYTHON_INCLUDE_DIRS}
  ${MSPASSHOME}/include
  ${PROJECT_BINARY_DIR}/include
  ${PROJECT_SOURCE_DIR}/include)

add_library(seispp STATIC ${sources_seispp})
target_link_libraries(seispp PRIVATE dsap ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

install (TARGETS seispp DESTINATION lib)
