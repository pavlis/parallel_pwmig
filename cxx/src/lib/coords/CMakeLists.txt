FILE(GLOB sources_coords *.c)
include_directories(
  ${PROJECT_SOURCE_DIR}/include/pwmig/utility)

add_library(coords STATIC ${sources_coords})
target_link_libraries(coords PRIVATE utility ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

install (TARGETS coords DESTINATION lib)
