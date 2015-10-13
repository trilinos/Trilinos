if(${CMAKE_PROJECT_NAME} STREQUAL "SEACAS")
  SET(SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
    Ioss        src               PT  REQUIRED
    Ioex        src/exodus        PT  OPTIONAL
    Iofx        src/exo_fpp       PT  OPTIONAL
    Iopx        src/exo_par       PT  OPTIONAL
    Ioexo_fac   src/exo_fac       PT  OPTIONAL
    Iogn        src/generated     PT  REQUIRED
    Iopg        src/pamgen        PT  OPTIONAL
    Iohb        src/heartbeat     PT  REQUIRED
    Iotr        src/transform     PT  REQUIRED
    Ionit       src/init          PT  REQUIRED
    IoMain	src/main          PT  REQUIRED
    IoUtest	src/utest         PT  REQUIRED
    )

  SET(LIB_OPTIONAL_DEP_PACKAGES Exodus)
  SET(LIB_OPTIONAL_DEP_TPLS XDMF HDF5 Pamgen Zoltan MPI)

else()
  SET(LIB_OPTIONAL_DEP_PACKAGES SEACASExodus Pamgen Zoltan)
  SET(LIB_OPTIONAL_DEP_TPLS XDMF HDF5)
endif()

SET(LIB_REQUIRED_DEP_PACKAGES)
SET(TEST_REQUIRED_DEP_PACKAGES)
SET(TEST_OPTIONAL_DEP_PACKAGES)
SET(LIB_REQUIRED_DEP_TPLS)
SET(TEST_REQUIRED_DEP_TPLS)
SET(TEST_OPTIONAL_DEP_TPLS)

