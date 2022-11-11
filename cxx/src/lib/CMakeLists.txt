add_subdirectory(dsap)
add_subdirectory(gclgrid)
add_subdirectory(pwmigcore)
add_subdirectory(seispp)

add_library(pwmig $<TARGET_OBJECTS:dsap> $<TARGET_OBJECTS:gclgrid> $<TARGET_OBJECTS:pwmigcore> $<TARGET_OBJECTS:seispp> )
target_link_libraries(pwmig PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}  ${PYTHON_LIBRARIES} ${GSL_LIBRARIES} OpenMP::OpenMP_CXX ${MSPASS_LIBRARIES})

install (TARGETS pwmig DESTINATION lib)

