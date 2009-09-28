
IF(Trilinos_ENABLE_Fortran AND BUILD_ZOLTAN_F90_INTERFACE)

INCLUDE(PackageAddExecutable)

SET(ZFDRIVE_SOURCES
    fdr_migrate.f90 
    farg_typical.f
    mpistring.c
    lampmpi.f90 
    fdr_sort.f90 
    fdr_main.f90 
    fdr_loadbal.f90 
    fdr_input.f90 
    fdr_param_file.f90 
    fdr_chaco_io.f90 
    fdr_mm_io.f90 
    fdr_const.f90 
    mmio.f
    mpi_h.f
    )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

PACKAGE_ADD_EXECUTABLE(
  zfdrive
  NOEXEPREFIX
  SOURCES ${ZFDRIVE_SOURCES}
  COMM serial mpi
  )

#  What is a good location for executables zdrive and zCPPdrive?
#  When should they be moved there?  During install?  Or build?
#INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/zfdrive.exe DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

ENDIF()
