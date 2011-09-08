INCLUDE(SubPackageMacros)
INCLUDE(PackageAddExecutable)
INCLUDE(FortranSettings)

SUBPACKAGE(Gjoin)

SET(CMAKE_EXECUTABLE_SUFFIX "")

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES prtelb.F prtess.F prtnps.F
   chgint.f dbmirr.f expxyz.f initin.f matxyz.f movrea.f munnps.f ordix.f  rdgen.f  rotxyz.f wrgen.f  zmmap.f
   comand.f dbpini.f fixmap.f indexx.f intcnt.f mirss.f  munelb.f munxyz.f ordstr.f renelb.f setitl.f zmess.f  zmnps.f
   dbmir1.f dorot.f  getall.f inigen.f irennp.f movint.f muness.f offset.f outlog.f renix.f  setsta.f zmfixd.f zmxyz.f)

PACKAGE_ADD_LIBRARY(gjoin_app_lib SOURCES ${SOURCES} DEPLIBS supes suplib exodus_for)
PACKAGE_ADD_EXECUTABLE(gjoin NOEXEPREFIX NOEXESUFFIX LINKER_LANGUAGE Fortran SOURCES gjoin.f DEPLIBS gjoin_app_lib COMM serial)

install_executable(gjoin)

SUBPACKAGE_POSTPROCESS()

