INCLUDE(SubPackageMacros)
INCLUDE(PackageAddExecutable)
INCLUDE(FortranSettings)

SUBPACKAGE(Exotxt)

SET(CMAKE_EXECUTABLE_SUFFIX "")

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES 
   dbicon.f  dbinam.f  dbivtt.f  wress.f   wrname.f  wrstep.f
   dbielb.f  dbiqa.f   dbpini.f  rwpval.f  wrinit.f  wrnps.f   wrxyz.f
   dbign.f   dbiste.f  dbvini.f  wrelb.f   wrmap.f   wrqa.f)

PACKAGE_ADD_LIBRARY(exotxt_app_lib SOURCES ${SOURCES} DEPLIBS supes suplib exodus_for)
PACKAGE_ADD_EXECUTABLE(exotxt NOEXEPREFIX NOEXESUFFIX LINKER_LANGUAGE Fortran SOURCES exotxt.f DEPLIBS exotxt_app_lib supes suplib exodus_for)
install_executable(exotxt)

SUBPACKAGE_POSTPROCESS()

