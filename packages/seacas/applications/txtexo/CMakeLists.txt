INCLUDE(SubPackageMacros)
INCLUDE(PackageAddExecutable)
INCLUDE(FortranSettings)

SUBPACKAGE(Txtexo)

SET(CMAKE_EXECUTABLE_SUFFIX "")

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES 
   dboelb.f  rdelb.f   rdinit.f  rdnm1.f   rdxyz.f   rwpval.f  rwqa.f
   rdeb1.f   rdess.f   rdmap.f   rdnps.f   rwname.f  rwpval1.f rwstep.f)

PACKAGE_ADD_LIBRARY(txtexo_app_lib SOURCES ${SOURCES} DEPLIBS supes suplib exodus_for exodus)
PACKAGE_ADD_EXECUTABLE(txtexo NOEXEPREFIX NOEXESUFFIX HEADERS ${HEADERS} SOURCES txtexo.f DEPLIBS txtexo_app_lib COMM serial)
install_executable(txtexo)

SUBPACKAGE_POSTPROCESS()

