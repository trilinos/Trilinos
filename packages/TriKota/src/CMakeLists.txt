
INCLUDE(PackageLibraryMacros)

#
# Define the configure Dakota libraries
#

EXECUTE_PROCESS(
  COMMAND echo "-->TriKota: Launching Dakota configure script -- takes several minutes."
  WORKING_DIRECTORY ${PACKAGE_SOURCE_DIR}/Dakota
)

EXECUTE_PROCESS(
  WORKING_DIRECTORY ${PACKAGE_SOURCE_DIR}/Dakota
  COMMAND ./configure --without-graphics
    CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} MPICXX=${CMAKE_CXX_COMPILER}
    FC=${CMAKE_Fortran_COMPILER} F77=${CMAKE_Fortran_COMPILER}
    --prefix=${PACKAGE_BINARY_DIR}/Dakota
    --with-teuchos-include=${PACKAGE_SOURCE_DIR}/../teuchos/src
    --with-teuchos-lib=${PACKAGE_BINARY_DIR}/../teuchos/src
    OUTPUT_FILE ${PACKAGE_SOURCE_DIR}/Dakota/configure.out
  )

## ToDo: Add *all* compiler options.  See Makefile.export.in and
## PackageWriteExportMakefile.cmake for how to grab *all* the compiler options
## that you need to be consistent with building against Trilinos.
#  ...try these...
#    CC=mpicc CXX=mpicxx F77=mpif90 FC=mpif90

EXECUTE_PROCESS(
  COMMAND echo "<--TriKota: Finished Dakota configure script."
  WORKING_DIRECTORY ${PACKAGE_SOURCE_DIR}/Dakota
)

# Customize the build process, by adding a custom make
# target, which makes Dakota (in source tree) and then
# uses Dakota's PREFIX to install it in the build tree.
ADD_CUSTOM_COMMAND(
  OUTPUT ${PACKAGE_SOURCE_DIR}/Dakota/src/dakota
  COMMAND make -j 4 && make -j 4 install
  WORKING_DIRECTORY ${PACKAGE_SOURCE_DIR}/Dakota
  )

ADD_CUSTOM_TARGET( ${PACKAGE_NAME}_Dakota_libs ALL
  DEPENDS ${PACKAGE_SOURCE_DIR}/Dakota/src/dakota
  )

#
# Build the TriKota library
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../Dakota/src
  ${CMAKE_CURRENT_SOURCE_DIR}/../Dakota
  )

APPEND_SET(HEADERS
    TriKota_ConfigDefs.hpp
    TriKota_DirectApplicInterface.hpp
    TriKota_Driver.hpp
    TriKota_Version.hpp
  ) 

APPEND_SET(SOURCES
    TriKota_DirectApplicInterface.cpp
    TriKota_Driver.cpp
    TriKota_Version.cpp
  )

# Add Dakota libraries to the list of TriKota libraries
#SET(DAKOTA_LIBS
#    -ldakota -lpecos -lfftw3 -llhs -levidence -lsurfpack 
#    -lconmin -lddace -lfsudace -ljega -lnlpql -lnpsol -lopt
#    -lpsuade -lnewmat -lncsuopt -lquadrature -lcoliny -lcolin
#    -lpebbl -lutilib -l3po -lnappspack -lappspack -lconveyor
#    -lshared -lcdd -lcport -lamplsolver -llhs
#  )

SET(DAKOTA_LIBS
    dakota pecos fftw3 lhs evidence surfpack 
    conmin ddace fsudace jega nlpql npsol opt
    psuade newmat ncsuopt quadrature coliny colin
    pebbl utilib 3po nappspack appspack conveyor
    shared cdd cport amplsolver lhs
  )
LINK_DIRECTORIES(${PACKAGE_BINARY_DIR}/Dakota/lib)

PACKAGE_ADD_LIBRARY(
  trikota
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS ${DAKOTA_LIBS}
  )

# Custom "make install" command to copy contents of Dakota lib
FILE(GLOB_RECURSE dakota_lib_files ${PACKAGE_BINARY_DIR}/Dakota/lib/*)
INSTALL(FILES ${dakota_lib_files} DESTINATION lib)

# Custom "make install" command for dakota includes,
# preserving directory structure
SET(DAKINC  ${PACKAGE_BINARY_DIR}/Dakota/include)

# Q: need *hpp.inl  as well?
FILE(GLOB dakota_header_files ${DAKINC}/*.h  ${DAKINC}/*.H  ${DAKINC}/*.hpp)
INSTALL(FILES ${dakota_header_files} DESTINATION include)

FILE(GLOB tpo_header_files ${DAKINC}/3po/*.h)
INSTALL(FILES ${tpo_header_files} DESTINATION include/3po)

# Q: need ampl/*.hd as well?
FILE(GLOB ampl_header_files ${DAKINC}/ampl/*.h)
INSTALL(FILES ${ampl_header_files} DESTINATION include/ampl)

# Q: need boost subdirectories as well?
FILE(GLOB boost_header_files ${DAKINC}/boost/*.hpp)
INSTALL(FILES ${boost_header_files} DESTINATION include/boost)

FILE(GLOB colin_header_files ${DAKINC}/colin/*.h)
INSTALL(FILES ${colin_header_files} DESTINATION include/colin)

FILE(GLOB coliny_header_files ${DAKINC}/coliny/*.h)
INSTALL(FILES ${coliny_header_files} DESTINATION include/coliny)

FILE(GLOB coliny_header_files ${DAKINC}/coliny/*.h)
INSTALL(FILES ${coliny_header_files} DESTINATION include/coliny)

FILE(GLOB pebbl_header_files ${DAKINC}/pebbl/*.h)
INSTALL(FILES ${pebbl_header_files} DESTINATION include/pebbl)

FILE(GLOB pecos_header_files ${DAKINC}/pecos/*.h ${DAKINC}/pecos/*.hpp)
INSTALL(FILES ${pecos_header_files} DESTINATION include/pecos)

FILE(GLOB utilib_header_files ${DAKINC}/utilib/*.h)
INSTALL(FILES ${utilib_header_files} DESTINATION include/utilib)
