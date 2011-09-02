INCLUDE(SubPackageMacros)
INCLUDE(PackageAddExecutable)
INCLUDE(FortranSettings)

SUBPACKAGE(Mapvar)

SET(CMAKE_EXECUTABLE_SUFFIX "")

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES 
   avg.f     clsfil.f  error.f   frge.f    ininod.f  invcon.f  mkei.f    qadsrc.f  rdb2.f    shlsrc.f  srchh.f   tetsrc.f  wrtc.f
   bldpte.f  cntr.f    exth.f    getbnd.f  initlz.f  jacobn.f  mklstv.f  rank.f    rdinpt.f  sintpe.f  srchq.f   tranab.f
   bldptn.f  elgrad.f  extq.f    hexsrc.f  intrp3.f  lenstr.f  mkrnk.f   rda1.f    seton0.f  sintpn.f  srchs.f   trutbl.f
   bldsrf.f  elton0.f  exts.f    indexx.f  intrpe.f  node.f    rda2.f    seton1.f  srchge.f  srcht.f   version.f
   bs.f      elton1.f  flgrad.f  inielt.f  intrpn.f  mke.f     opnfil.f  rdb1.f    shapef.f  srchgt.f  stran.f   vol.f)

PACKAGE_ADD_LIBRARY(mapvar_app_lib SOURCES ${SOURCES} DEPLIBS supes suplib exodus_for)
PACKAGE_ADD_EXECUTABLE(mapvar NOEXEPREFIX NOEXESUFFIX HEADERS ${HEADERS} SOURCES mapvar.f DEPLIBS mapvar_app_lib COMM serial)
install_executable(mapvar)

SUBPACKAGE_POSTPROCESS()

