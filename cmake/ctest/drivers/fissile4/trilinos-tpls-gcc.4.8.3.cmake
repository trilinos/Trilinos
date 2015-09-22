#
# Trilinos TPLs
#

# Give a default TPL_INSTALL_DIR if not set!
IF (NOT TPL_INSTALL_DIR)
  SET(TPL_INSTALL_DIR_DEFAULT_STATIC
    "/projects/vera/gcc-4.8.3/tpls/opt_static")
  IF ( (NOT BUILD_SHARED_LIBS) AND (EXISTS "${TPL_INSTALL_DIR_DEFAULT_STATIC}") )
    SET(TPL_INSTALL_DIR_DEFAULT "${TPL_INSTALL_DIR_DEFAULT_STATIC}")
  ELSE()
    SET(TPL_INSTALL_DIR_DEFAULT "/projects/vera/gcc-4.8.3/tpls/opt")
  ENDIF()
  MESSAGE("-- " "TPL_INSTALL_DIR not set so setting by default to"
    " ${TPL_INSTALL_DIR_DEFAULT}!")
  SET(TPL_INSTALL_DIR  "${TPL_INSTALL_DIR_DEFAULT}"  CACHE  PATH
    "Set by default in vera-tpls-std.cmake")
ENDIF()

# NOTE: Above we turn on all the TPLs, reguardless of packages are enabled.
# If a VERA repo is not enabled that defines one of these TPLs, the enables
# are no-ops.  This ensures that the same XXX_config.h files get generated.
# For exmaple, Teuchos defines an optional dependence on QT.  However, if you
# don't enable the SCALE packages, QT (which is required for SCALE) does not
# get enabled and Teuchos_config.h does not define HAVE_TEUCHOS_QT.  This
# causes the Teuchos_config.h file to keep changing with different package
# enables.  If you don't want to enable one of these TPLs (because you know
# you don't need it) you can just pass in -DTPL_ENABLE_<TPLNAME>=OFF and it
# will override the above cache enable.

SET(BLAS_LIBRARY_DIRS "${TPL_INSTALL_DIR}/lapack-3.3.1/lib" CACHE FILEPATH "")
SET(LAPACK_LIBRARY_DIRS "${TPL_INSTALL_DIR}/lapack-3.3.1/lib" CACHE FILEPATH "")

SET(Boost_DIR "${TPL_INSTALL_DIR}/boost-1.55.0")
SET(Boost_INCLUDE_DIRS "${Boost_DIR}/include" CACHE FILEPATH "")

SET(BoostLib_LIBRARY_NAMES "boost_system;boost_thread;boost_program_options" CACHE STRING  "")
SET(BoostLib_INCLUDE_DIRS ${Boost_INCLUDE_DIRS}  CACHE FILEPATH "")
SET(BoostLib_LIBRARY_DIRS "${Boost_DIR}/lib" CACHE FILEPATH "")

SET(Zlib_INCLUDE_DIRS "${TPL_INSTALL_DIR}/zlib-1.2.7/include" CACHE FILEPATH "" )
SET(Zlib_LIBRARY_DIRS "${TPL_INSTALL_DIR}/zlib-1.2.7/lib" CACHE FILEPATH "" )

SET(HDF5_DIR "${TPL_INSTALL_DIR}/hdf5-1.8.10")
SET(HDF5_LIB_DIR "${HDF5_DIR}/lib")
SET(HDF5_LIBRARY_NAMES "hdf5_hl;hdf5_cpp;hdf5_fortran;hdf5;z" CACHE STRING "")
SET(HDF5_INCLUDE_DIRS
  "${HDF5_DIR}/include;${TPL_COMMON_DIR}/zlib-1.2.7/include"
  CACHE FILEPATH "" )
SET(HDF5_LIBRARY_DIRS
  "${HDF5_LIB_DIR};${TPL_INSTALL_DIR}/zlib-1.2.7/lib"
  CACHE FILEPATH "" )

SET(QT_REQUIRED_VERSION 4.7.1 CACHE STRING "" )
SET(QT_QMAKE_EXECUTABLE ${TPL_INSTALL_DIR}/qt-4.8.2/bin/qmake CACHE FILEPATH "" )

# TriKota Hack!
# Must force this varible in the Dakota configure process
SET(BOOST_INCLUDEDIR ${Boost_INCLUDE_DIRS}  CACHE PATH  "" FORCE)

# Disable TPLs we don't have on this system by default
SET(TPL_ENABLE_Netcdf  OFF  CACHE  BOOL
  "Set in trilinos-tpls-gcc.4.8.3.cmake" )
