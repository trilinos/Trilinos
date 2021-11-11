# Source: https://devblogs.microsoft.com/cppblog/clear-functional-c-documentation-with-sphinx-breathe-doxygen-cmake/
# Author: Evan Harvey <eharvey@sandia.gov>
find_package(Doxygen REQUIRED)
find_package(Sphinx REQUIRED)

set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/docs/sphinx)
set(KOKKOS_INCLUDE_DIR ${Kokkos_DIR}/../../../include)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/conf.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/index.rst DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Find all public headers in KokkosKernels and Kokkos::kokkos
file(GLOB_RECURSE ${PROJECT_NAME}_PUBLIC_HEADERS ${PROJECT_SOURCE_DIR}/src/*.hpp ${KOKKOS_INCLUDE_DIR}/*.hpp)

set(DOXYGEN_KOKKOSKERNELS_INPUT_DIR ${PROJECT_SOURCE_DIR}/src)
set(DOXYGEN_KOKKOS_INPUT_DIR ${KOKKOS_INCLUDE_DIR})
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/doxygen/)
set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/xml/index.xml)
set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

#Replace variables inside @@ with the current values
configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT})

file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR}) #Doxygen won't create this for us

add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
        DEPENDS ${${PROJECT_NAME}_PUBLIC_HEADERS}
        COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        DEPENDS Doxyfile
        COMMENT "Generating docs")

add_custom_target(Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})


add_custom_target(Sphinx ALL
        COMMAND ${SPHINX_EXECUTABLE} -b html
        # Tell Breathe where to find the Doxygen output
        -Dbreathe_projects.${PROJECT_NAME}=${DOXYGEN_OUTPUT_DIR}/xml
        ${SPHINX_SOURCE} ${SPHINX_BUILD}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        DEPENDS Doxygen
        COMMENT "Generating documentation with Sphinx")