INCLUDE(SubPackageMacros)
INCLUDE(PackageAddExecutable)
INCLUDE(FortranSettings)

SUBPACKAGE(Grope)

SET(CMAKE_EXECUTABLE_SUFFIX "")

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES 
   chkcnt.f ckess.f  dbsbel.f indexx.f nxmmax.f prelem.f prinit.f prnode.f prssv.f  prxyz.f  rdess.f  rdname.f rdstep.f rixwrd.f setprc.f
   chkrng.f ckmap.f  dbselb.f invcon.f preb1.f  press.f  prmap.f  prnps.f  prstep.f rdcord.f rdmap.f  rdnmap.f rdtims.f selblk.f tostep.f
   ckeb1.f  cknps.f  getprc.f limits.f prebn.f  prfrm.f  prmmax.f prnsv.f  prtims.f rdeb1.f  rdmmax.f rdnps.f  rixid.f  selmap.f wdberr.f
   ckelb.f  comand.f numstr.f prelb.f  prglob.f prname.f prqa.f   prvers.f rdelb.f  rdnam1.f rdqa.f   rixint.f selset.f)

PACKAGE_ADD_LIBRARY(grope_app_lib SOURCES ${SOURCES} DEPLIBS supes suplib exodus_for)
PACKAGE_ADD_EXECUTABLE(grope NOEXEPREFIX NOEXESUFFIX HEADERS ${HEADERS} SOURCES grope.f DEPLIBS grope_app_lib COMM serial)

install_executable(grope)

SUBPACKAGE_POSTPROCESS()

