add_subdirectory(coords)
#add_subdirectory(pwmigcore)
add_subdirectory(stock)
add_subdirectory(utility)


#file(WRITE ${CMAKE_BINARY_DIR}/src/dummy.cc "")
#add_library(mspass ${CMAKE_BINARY_DIR}/src/dummy.cc)
#target_link_libraries(mspass PUBLIC seismic utility deconvolution algorithms)
#
#install (TARGETS mspass DESTINATION lib)
