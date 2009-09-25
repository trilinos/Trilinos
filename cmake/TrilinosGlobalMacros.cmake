

#
# This file contains global-level macros that are specific to Trilinos
#



#
# Macro that defines Trilinos testing support
#

MACRO(TRILINOS_SETUP_TESTING_SUPPORT)
  
  IF (WIN32)
    SET(Trilinos_ENABLE_FRAMEWORK_UNIT_TESTS_DEFAULT OFF)
  ELSE()
    SET(Trilinos_ENABLE_FRAMEWORK_UNIT_TESTS_DEFAULT ${${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE})
  ENDIF()
  
  # 2008/10/17: rabartl: Above, I can not turn these tests on by default
  # with cygwin because the custom script target is not working for some
  # reason.
  
  ADVANCED_OPTION(Trilinos_ENABLE_FRAMEWORK_UNIT_TESTS
    "Enable all Trilinos framework unit tests by default."
    ${Trilinos_ENABLE_FRAMEWORK_UNIT_TESTS_DEFAULT}
    )
  
  ADVANCED_OPTION(Trilinos_ENABLE_DEPENCENCY_UNIT_TESTS
    "Enable Trilinos Framework dependency unit tests."
    ${Trilinos_ENABLE_FRAMEWORK_UNIT_TESTS}
    )
  
  ADVANCED_OPTION(Trilinos_ENABLE_TESTING_UNIT_TESTS
    "Enable Trilinos CTest testing support unit tests."
    ${Trilinos_ENABLE_FRAMEWORK_UNIT_TESTS}
    )
  
  ADVANCED_OPTION(Trilinos_ENABLE_PYTHON_UNIT_TESTS
    "Enable Trilinos python script unit tests."
    ${Trilinos_ENABLE_FRAMEWORK_UNIT_TESTS}
    )

  # Add the directory for the unit tests
  ADD_SUBDIRECTORY(cmake)

  CONFIGURE_FILE(
    ${Trilinos_SOURCE_DIR}/cmake/ctest/CTestCustom.ctest.in
    ${Trilinos_BINARY_DIR}/CTestCustom.ctest
    )

ENDMACRO()


#
#  Function for helping set up exclude files only for the packages
#  that will not be supporting autotools.
#  Returns a list of the given file name with a path for packages
#  that are not supporting autotools anymore.
#
#  example: PACKAGE_APPLY_TO_NO_AUTOTOOLS_PACKAGES("configure.ac" list)
#    assuming that the packages epetra and teuchos are not supporting 
#    autotools anymore then the return value would be:
#    "epetra/configure.ac;teuchos/configure.ac"
#
#

FUNCTION(APPLY_TO_NO_AUTOTOOLS_PACKAGES FILE_NAME LIST_RETURN)
  SET(NON_AUTOTOOLS_PACKAGES
    /packages/amesos
    /packages/anasazi
    /packages/aztecoo
    /packages/belos
    /packages/didasko
    /packages/epetra
    /packages/epetraext
    /packages/fei
    /packages/galeri
    /packages/ifpack
    /packages/intrepid
    /packages/isorropia
    /packages/kokkos
    /packages/komplex
    /packages/meros
    /packages/ml
    /packages/moertel
    /packages/moocho
    /packages/nox
    /packages/pamgen
    /packages/phalanx
    /packages/phdmesh
    /packages/pliris
    /packages/PyTrilinos
    /packages/rtop
    /packages/rythmos
    /packages/sacado
    /packages/shards
    /packages/stratimikos
    /packages/Sundance
    /packages/teuchos
    /packages/ThreadPool
    /packages/thyra
    /packages/tpetra
    /packages/trilinoscouplings
    /packages/triutils
  )
  
  FOREACH(PACKAGE ${NON_AUTOTOOLS_PACKAGES})
    SET(LIST_RETURN_TMP ${LIST_RETURN_TMP} ${PACKAGE}/${FILE_NAME} ${PACKAGE}/\(.*/\)*${FILE_NAME})
  ENDFOREACH()
  
  SET(${LIST_RETURN} ${LIST_RETURN_TMP} PARENT_SCOPE)
ENDFUNCTION()

#
# Macro that defines Trilinos packaging options:
#

MACRO(TRILINOS_DEFINE_PACKAGING)
  APPLY_TO_NO_AUTOTOOLS_PACKAGES("configure.ac" CONFIGURE_AC_LIST)
  APPLY_TO_NO_AUTOTOOLS_PACKAGES("configure"    CONFIGURE_LIST)
  APPLY_TO_NO_AUTOTOOLS_PACKAGES("Makefile.am"  MAKEFILE_AM_LIST)
  APPLY_TO_NO_AUTOTOOLS_PACKAGES("Makefile.in"  MAKEFILE_AC_LIST)
  APPLY_TO_NO_AUTOTOOLS_PACKAGES(".*.m4"        M4_LIST)
  APPLY_TO_NO_AUTOTOOLS_PACKAGES("bootstrap"    BOOTSTRAP_LIST)
  APPLY_TO_NO_AUTOTOOLS_PACKAGES("config/"      CONFIG_LIST)

    
  SET(CPACK_SOURCE_IGNORE_FILES
    /CVS/
    ".cvsignore"
    classicMakefile
    /packages/CTrilinos
    /packages/ForTrilinos
    /packages/ITAPS
    /packages/globipack
    /packages/mesquite
    /packages/optika
    /packages/optipack
    /packages/stokhos
    /packages/tifpack
    /packages/TriKota
    /packages/aristos
    /packages/claps
    /packages/external
    /packages/jpetra
    /packages/new_package
    /packages/rbgen
    /packages/WebTrilinos
    ${CONFIGURE_AC_LIST}
    ${CONFIGURE_LIST}
    ${MAKEFILE_AM_LIST}
    ${MAKEFILE_AC_LIST}
    ${M4_LIST}
    ${BOOTSTRAP_LIST}
    ${CONFIG_LIST}
    /packages/configure.ac
    /packages/configure
    /packages/Makefile.am
    /packages/Makefile.in
    Trilinos/configure.ac
    Trilinos/configure
    Trilinos/Makefile.am
    Trilinos/Makefile.in
    Trilinos/bootstrap
    Trilinos/config
    Trilinos/doc/[^b]
    ".*.pyc"
    /SIERRA/
    /commonTools/test
    /commonTools/scripts
    /packages/PyTrilinos/Notes.txt
    /packages/PyTrilinos/aclocal.m4
    /packages/PyTrilinos/bootstrap
    /packages/PyTrilinos/config
    /packages/PyTrilinos/lib
    /packages/PyTrilinos/macdist
    /packages/PyTrilinos/shared
    /packages/PyTrilinos/src/NOX
    /packages/PyTrilinos/src/PyTrilinos_config.h.in
    /packages/PyTrilinos/src/depend
    /packages/PyTrilinos/src/setup.py
    /packages/PyTrilinos/src-boost
    /packages/zoltan/test/ch_brack2_3
    /packages/zoltan/test/ch_bug
    /packages/zoltan/test/ch_degenerate
    /packages/zoltan/test/ch_degenerateAA
    /packages/zoltan/test/ch_drake
    /packages/zoltan/test/ch_ewgt
    /packages/zoltan/test/ch_grid20x19
    /packages/zoltan/test/ch_hammond
    /packages/zoltan/test/ch_hammond2
    /packages/zoltan/test/ch_nograph
    /packages/zoltan/test/ch_onedbug
    /packages/zoltan/test/ch_random
    /packages/zoltan/test/ch_serial
    /packages/zoltan/test/ch_slac
    /packages/zoltan/test/ch_vwgt
    /packages/zoltan/test/ch_vwgt2
    /packages/zoltan/test/hg_cage10
    /packages/zoltan/test/hg_diag500_4
    /packages/zoltan/test/hg_ewgt
    /packages/zoltan/test/hg_felix
    /packages/zoltan/test/hg_ibm03
    /packages/zoltan/test/hg_ml27
    /packages/zoltan/test/hg_nograph
    /packages/zoltan/test/hg_vwgt
    /packages/zoltan/test/nem_ti_20k
    /packages/zoltan/test/nem_ti_4k
    /packages/zoltan/test/misc_siefert
    /packages/zoltan/test/th
    /packages/zoltan/test/bin
    /packages/zoltan/doc/Zoltan_html/tu_html
    /packages/zoltan/src/ZoltanComponent
    /packages/zoltan/src/driver_old
    /packages/zoltan/src/fdriver_old
    /packages/amesos/doc/AmesosOverview
    /packages/amesos/doc/PARA06
    /packages/anasazi/doc/TOMS
    /packages/anasazi/doc/OrthoStudy
    /packages/anasazi/doc/ThyraPerf
    /packages/aztecoo/doc/AZ_capture_matrix_howto.txt
    /packages/aztecoo/doc/Aztec2.0
    /packages/aztecoo/doc/Aztec2.1
    /packages/aztecoo/doc/Managing_conditioning_howto.txt
    /packages/aztecoo/doc/UserGuide
    /packages/aztecoo/doc/azteclogo.gif
    /packages/aztecoo/doc/read_captured_matrix.c
    /packages/aztecoo/example/AztecOO_RecursiveCall
    /packages/aztecoo/example/Epetra_MsrMatrix_AztecOO
    /packages/aztecoo/example/Epetra_MsrMatrix_PowerMethod
    /packages/aztecoo/example/IfpackIctAztecOO
    /packages/aztecoo/example/IfpackAztecOO
    /packages/aztecoo/example/IfpackVbrAztecOO
    /packages/aztecoo/example/MLAztecOO
    /packages/aztecoo/example/azoo_iterate_hb
    /packages/aztecoo/example/aztec_app
    /packages/aztecoo/example/aztec_hb
    /packages/galeri/src-pfem
    /packages/galeri/example-pfem
    /packages/tpetra/doc/CodingGuidelines
    /packages/tpetra/doc/TpetraDesign
    /packages/kokkos/doc
  )
  
  IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("Exclude files when building source packages")
    FOREACH(item ${CPACK_SOURCE_IGNORE_FILES})
      MESSAGE(${item})
    ENDFOREACH()
  ENDIF()
  
  SET(CPACK_PACKAGE_DESCRIPTION "Trilinos provides algorithms and technologies for the solution of large-scale, complex multi-physics engineering and scientific problems.")
  SET(CPACK_PACKAGE_FILE_NAME "trilinos-setup-${Trilinos_VERSION}")
  SET(CPACK_PACKAGE_INSTALL_DIRECTORY "Trilinos ${Trilinos_VERSION}")
  SET(CPACK_PACKAGE_REGISTRY_KEY "Trilinos ${Trilinos_VERSION}")
  SET(CPACK_PACKAGE_NAME "trilinos")
  SET(CPACK_PACKAGE_VENDOR "Sandia National Laboratories")
  SET(CPACK_PACKAGE_VERSION "${Trilinos_VERSION}")
  SET(CPACK_RESOURCE_FILE_README "${Trilinos_SOURCE_DIR}/README")
  SET(CPACK_RESOURCE_FILE_LICENSE "${Trilinos_SOURCE_DIR}/README")
  SET(CPACK_SOURCE_GENERATOR "TGZ;TBZ2")
  SET(CPACK_SOURCE_FILE_NAME "trilinos-source-${Trilinos_VERSION}")
  SET(CPACK_COMPONENTS_ALL ${Trilinos_PACKAGES})
  
  PACKAGE_ARCH_GET_ENABLED_LIST( Trilinos_PACKAGES Trilinos ON
    FALSE ENABLED_PACKAGES NUM_ENABLED)
  string(REPLACE " " ";" ENABLED_PACKAGES "${ENABLED_PACKAGES}")
  
  #message("ENABLED PACKAGES: ${ENABLED_PACKAGES} ${NUM_ENABLED}")
  FOREACH(PKG ${ENABLED_PACKAGES})
    IF(NOT "${${PKG}_LIB_REQUIRED_DEP_PACKAGES}" STREQUAL "")
        string(TOUPPER ${PKG} UPPER_PKG)
        #message("${UPPER_PKG} depends on : ${${PKG}_LIB_REQUIRED_DEP_PACKAGES}")
        SET(CPACK_COMPONENT_${UPPER_PKG}_DEPENDS ${${PKG}_LIB_REQUIRED_DEP_PACKAGES})
    ENDIF()
    #message("${PKG} depends on : ${${PKG}_LIB_REQUIRED_DEP_PACKAGES}")
  ENDFOREACH()

  
  IF(WIN32)
    SET(CPACK_GENERATOR "NSIS")
    SET(CPACK_NSIS_MODIFY_PATH OFF)
  ENDIF()
  
  INCLUDE(CPack)

ENDMACRO()
