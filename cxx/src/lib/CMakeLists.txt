add_subdirectory(dsap)
add_subdirectory(gclgrid)
# Temporarily disable to get other libraries functional first
add_subdirectory(pwmigcore)
add_subdirectory(seispp)


file(WRITE ${CMAKE_BINARY_DIR}/src/dummy.cc "")
add_library(pwmig ${CMAKE_BINARY_DIR}/src/dummy.cc)
target_link_libraries(pwmig PUBLIC dsap gclgrid pwmigcore seispp)
#
#install (TARGETS pwmig DESTINATION lib)
