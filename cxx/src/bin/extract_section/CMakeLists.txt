# pybind11_INCLUDE and PYTHON_INCLUDE needed because Metadata.h requires it - not a good thing
include_directories(
  ${pybind11_INCLUDE_DIR}
  ${PYTHON_INCLUDE_DIRS}
  #${MSPASSHOME}/include}
  /home/pavlis/mspass/include
  ${PROJECT_BINARY_DIR}/include
  ${PROJECT_SOURCE_DIR}/include)

add_executable (extract_section extract_section.cc)
target_link_libraries(extract_section PRIVATE pwmig ${MSPASS_LIBRARIES} ${pybind11_LIBRARIES})
