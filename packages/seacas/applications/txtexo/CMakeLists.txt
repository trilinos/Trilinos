INCLUDE(TribitsSubPackageMacros)
INCLUDE(TribitsAddExecutable)
INCLUDE(FortranSettings)

TRIBITS_SUBPACKAGE(Txtexo)

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES 
   dboelb.f  rdelb.f   rdinit.f  rdnm1.f   rdxyz.f   rwpval.f  rwqa.f
   rdeb1.f   rdess.f   rdmap.f   rdnps.f   rwname.f  rwpval1.f rwstep.f)

TRIBITS_ADD_LIBRARY(txtexo_app_lib SOURCES ${SOURCES} DEPLIBS supes suplib exodus_for exodus)
TRIBITS_ADD_EXECUTABLE(txtexo NOEXEPREFIX NOEXESUFFIX LINKER_LANGUAGE Fortran SOURCES txtexo.f DEPLIBS txtexo_app_lib COMM serial mpi)
install_executable(txtexo)

TRIBITS_SUBPACKAGE_POSTPROCESS()

