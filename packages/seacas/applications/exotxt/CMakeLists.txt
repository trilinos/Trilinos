INCLUDE(TribitsSubPackageMacros)
INCLUDE(TribitsAddExecutable)
INCLUDE(FortranSettings)

TRIBITS_SUBPACKAGE(Exotxt)

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES 
   et_dbinam.f  et_dbivtt.f  et_wress.f   et_wrname.f  et_wrstep.f
   et_dbielb.f  et_dbiqa.f   et_dbpini.f  et_wrprop.f  et_wrinit.f
   et_wrnps.f   et_wrxyz.f   et_dbign.f   et_dbiste.f  et_dbvini.f
   et_wrelb.f   et_wrmap.f   et_wrqa.f)

TRIBITS_ADD_LIBRARY(exotxt_app_lib SOURCES ${SOURCES} DEPLIBS supes suplib exodus_for)
TRIBITS_ADD_EXECUTABLE(exotxt NOEXEPREFIX NOEXESUFFIX LINKER_LANGUAGE Fortran SOURCES exotxt.f DEPLIBS exotxt_app_lib supes suplib exodus_for COMM serial mpi)
install_executable(exotxt)

TRIBITS_SUBPACKAGE_POSTPROCESS()

