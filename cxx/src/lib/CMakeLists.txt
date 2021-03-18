add_subdirectory(coords)
#add_subdirectory(pwmigcore)
add_subdirectory(stock)
add_subdirectory(utility)


file(WRITE ${CMAKE_BINARY_DIR}/src/dummy.cc "")
add_library(pwmig ${CMAKE_BINARY_DIR}/src/dummy.cc)
target_link_libraries(pwmig PUBLIC pwmigutil)
#
#install (TARGETS pwmig DESTINATION lib)
