add_subdirectory(gclgrid)
add_subdirectory(seispp)

file(WRITE ${CMAKE_BINARY_DIR}/src/dummy.cc "")
add_library(pwmigutil ${CMAKE_BINARY_DIR}/src/dummy.cc)
target_link_libraries(pwmigutil PUBLIC utility)

install (TARGETS pwmigutil DESTINATION lib)

