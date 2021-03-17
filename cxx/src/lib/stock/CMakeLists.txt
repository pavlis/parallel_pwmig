FILE(GLOB sources_stock *.c)
include_directories(
  ${PROJECT_SOURCE_DIR}/include/pwmig/utility)

add_library(stock STATIC ${sources_stock})
target_link_libraries(stock PRIVATE utility ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

install (TARGETS stock DESTINATION lib)
