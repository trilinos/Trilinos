# Install an executable and dependent files to the same directory as
# ${CMAKE_CURRENT_SOURCE_DIR}, but under ${CMAKE_INSTALL_PREFIX}.
#
# For example, if called in packages/muelu/example/basic/CMakeLists.txt, then we will install to
#   ${CMAKE_INSTALL_PREFIX}/packages/muelu/example/basic/
#
# Usage::
#
#   muelu_install( <exe_root_name>
#     [ADDITIONAL_FILES <file0> <file1> ...]
#     )
#
# For example, we might have added an executable and a couple of input files like so:
#
# tribits_add_executable(
#      Driver
#      SOURCES driver_source.cpp
# )
#
# set(XML_INPUTS
#    mediocre_input.xml amazing_input.xml
# )
#
# tribits_copy_files_to_binary_dir(
#     inputs_cp
#     SOURCE_FILES
#     ${XML_INPUTS}
# )
#
# In order to add these to the install, we do
#
# muelu_install(
#    Driver
#    ADDITIONAL_FILES
#    ${XML_INPUTS}
#  )
#
function(muelu_install TARGET)

  cmake_parse_arguments(
    PARSE
    ""
    ""
    "ADDITIONAL_FILES"
    ${ARGN}
  )

  # set the install directory to the current source folder, but relative to install prefix
  set(INSTALL_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  file(RELATIVE_PATH INSTALL_DIRECTORY ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
  set(INSTALL_DIRECTORY ${CMAKE_INSTALL_PREFIX}/${INSTALL_DIRECTORY})

  # install executable
  INSTALL(TARGETS "${PACKAGE_NAME}_${TARGET}" COMPONENT ${PACKAGE_NAME} RUNTIME DESTINATION ${INSTALL_DIRECTORY})

  # install additional files
  INSTALL(FILES ${PARSE_ADDITIONAL_FILES} DESTINATION ${INSTALL_DIRECTORY})

endfunction()
