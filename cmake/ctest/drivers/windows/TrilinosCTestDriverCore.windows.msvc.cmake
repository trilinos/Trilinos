# Get some variables from the environment.
FILE(TO_CMAKE_PATH "$ENV{GIT_EXE}" GIT_EXE)

IF(DEFINED ENV{SEMS_DIR})
  FILE(TO_CMAKE_PATH "$ENV{SEMS_DIR}" SEMS_DIR)
ELSE()
  SET(SEMS_DIR "C:/projects/sems/install/win-x86_64")
ENDIF()


INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../TrilinosCTestDriverCore.cmake")

MACRO(TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER)

  # Base of Trilinos/cmake/ctest then BUILD_DIR_NAME
  SET( CTEST_DASHBOARD_ROOT "${TRILINOS_CMAKE_DIR}/../../${BUILD_DIR_NAME}" )
  
  SET( CTEST_NOTES_FILES "${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}" )
  
  SET( CTEST_CMAKE_GENERATOR "Ninja" )
  SET( COMPILER_VERSION "MSVC14")
  
  SET(Trilinos_PACKAGES AztecOO Epetra ML MueLu Pamgen ROL Xpetra Zoltan
                        Triutils Teuchos Belos EpetraExt Amesos Ifpack TrilinosSS)

  # Blacklist packages					  
  SET(Trilinos_EXCLUDE_PACKAGES Galeri Kokkos Sacado Intrepid Thyra)
  FOREACH(PKG ${Trilinos_EXCLUDE_PACKAGES})
    LIST(APPEND EXTRA_SYSTEM_CONFIGURE_OPTIONS "-DTrilinos_ENABLE_${PKG}:BOOL=OFF")
  ENDFOREACH()
  
  # Setup other configure options
  LIST(APPEND EXTRA_SYSTEM_CONFIGURE_OPTIONS
       "-DCMAKE_BUILD_TYPE:STRING=${BUILD_TYPE}"
       "-DBUILD_SHARED_LIBS:BOOL=OFF"
       "-DTrilinos_ENABLE_TESTS:BOOL=ON"
       "-DTrilinos_ENABLE_EXAMPLES:BOOL=OFF"
       "-DTrilinos_ENABLE_FORTRAN:BOOL=OFF"
       "-DTrilinos_ENABLE_DEBUG:BOOL=OFF"
       "-DTrilinos_ENABLE_CXX11:BOOL=ON"
       "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
       "-DTrilinos_GENERATE_REPO_VERSION_FILE:BOOL=OFF"
       "-DCMAKE_OBJECT_PATH_MAX=500"
       "-DTPL_ENABLE_MPI:BOOL=ON"
       "-DMPI_BASE_DIR:PATH=${SEMS_DIR}/compiler/Microsoft MPI/8.1.12438.1084"
       "-DBLAS_LIBRARY_NAMES:STRING=blas\\\;libf2c"
       "-DBLAS_LIBRARY_DIRS:STRING=${SEMS_DIR}/tpl/clapack/3.2.1/VisualStudio/14.0/base/lib"
       "-DLAPACK_LIBRARY_NAMES:STRING=lapack"
       "-DLAPACK_LIBRARY_DIRS:STRING=${SEMS_DIR}/tpl/clapack/3.2.1/VisualStudio/14.0/base/lib"
       "-DPERL_EXECUTABLE:FILEPATH=${SEMS_DIR}/compiler/strawberry_perl/5.24.1.1/perl/bin/perl.exe"
       "-DTPL_ENABLE_DLlib:BOOL=OFF"
       "-DTPL_ENABLE_Pthread:BOOL=OFF")
  
  TRILINOS_CTEST_DRIVER()
ENDMACRO()
