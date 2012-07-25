#
# Define all of the TPL header paths, library dirs, and libraries for all CASL
# VRI software that does *not* depend on the C++ (GCC or Intel) compiler used.
#

SET(TPL_ENABLE_BinUtils  ON  CACHE BOOL "")
SET(TPL_ENABLE_JDK  ON  CACHE BOOL "")
SET(TPL_ENABLE_OpenSSL ON  CACHE BOOL "")
SET(TPL_ENABLE_TCL  ON  CACHE BOOL "")
SET(TPL_ENABLE_Zlib  ON  CACHE BOOL "")
SET(TPL_ENABLE_Netcdf ON CACHE BOOL "")
SET(TPL_ENABLE_SPRNG  OFF CACHE BOOL "")
#SET(TPL_ENABLE_MOOSE_HYPRE ON CACHE BOOL "")
SET(TPL_ENABLE_Libmesh ON CACHE BOOL "")
SET(TPL_ENABLE_Libmesh_Contrib ON CACHE BOOL "")

# We don't have the Matio TPL for SEACAS
SET(TPL_ENABLE_Matio OFF CACHE BOOL "")

SET(ANT_PATH  /usr/bin  CACHE FILEPATH "")
SET(TCLSH_PATH  /usr/bin  CACHE FILEPATH "")
SET(Netcdf_INCLUDE_DIRS /opt/gcc-4.5.1/tpls/netcdf-4.1.1/include CACHE FILEPATH "")
SET(Netcdf_LIBRARY_DIRS /opt/gcc-4.5.1/tpls/netcdf-4.1.1/lib CACHE FILEPATH "")
SET(TCL_LIBRARY_NAMES "tcl8.5"  CACHE STRING  "")
SET(JDK_INCLUDE_DIRS "/usr/lib/jvm/java/include;/usr/lib/jvm/java/include/linux"  CACHE FILEPATH "")
SET(JDK_LIBRARY_DIRS /usr/lib/jvm/java/jre/lib/amd64/server  CACHE FILEPATH "")
SET(SILO_INCLUDE_DIRS  /projects/gcc-4.6.1/tpls/silo-4.7.2/serial/include CACHE FILEPATH "")
SET(SILO_LIBRARY_DIRS  /projects/gcc-4.6.1/tpls/silo-4.7.2/serial/lib  CACHE FILEPATH "")
SET(LIBXML2_INCLUDE_DIRS  /usr/include/libxml2  CACHE FILEPATH "")
SET(LIBXML2_LIBRARY_DIRS  /usr/lib64  CACHE FILEPATH "")
SET(SCALE_INSTALL_PATH "/opt/scale/Linux_x86_64"         CACHE FILEPATH "")
SET(SCALE_INC_DIRS  AmpxLib       fulcrum              Resource   Cell               Serialize      Input        ModuleFactory   XSProc   Module          TextFileWriter   reporter    CACHE FILEPATH "")
SET(SCALE_INC_FILES AmpxLibrary.h SimpleFulcrumInput.h Resource.h CellLibRegistrar.h Serializable.h Input.h      ModuleFactory.h XSProc.h ControlModule.h textfilewriter.h reporter.h  CACHE FILEPATH "")
SET(SCALE_LIB_NAMES AmpxLib       fulcrumlib           Resource   Cell               Serialize      InputLibrary ModuleFactory   XSProc   Module          TextFileWriter   reporter    CACHE FILEPATH "")
SET(HYDRATPL_INCLUDE_DIRS /projects/TPLs/hydrath/gcc461_jul19/include CACHE FILEPATH "")
SET(HYDRATPL_LIBRARY_DIRS /projects/TPLs/hydrath/gcc461_jul19/lib     CACHE FILEPATH "")
SET(HYDRATPL_BINARY_DIRS  /projects/TPLs/hydrath/gcc461_jul19/bin     CACHE FILEPATH "")
SET(Libmesh_INCLUDE_DIRS  /projects/TPLs/moose_libs/libmesh/include/flat_headers CACHE FILEPATH "")
SET(Libmesh_LIBRARY_DIRS  /projects/TPLs/moose_libs/libmesh/lib/x86_64-unknown-linux-gnu_opt CACHE FILEPATH "")
SET(Libmesh_Contrib_INCLUDE_DIRS  /projects/TPLs/moose_libs/libmesh/contrib/flat_headers CACHE FILEPATH "")
SET(Libmesh_Contrib_LIBRARY_DIRS  /projects/TPLs/moose_libs/libmesh/contrib/lib/x86_64-unknown-linux-gnu_opt CACHE FILEPATH "")
