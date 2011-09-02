INCLUDE(SubPackageMacros)
INCLUDE(PackageAddExecutable)
INCLUDE(FortranSettings)

SUBPACKAGE(Mapvar-kd)

SET(CMAKE_EXECUTABLE_SUFFIX "")

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES 
   bldsrf.f    ebbyeb.blk  exts.f      inival.blk  mke.f       optkd.h     schdat.blk  srchq.f     tranab.f
   bmesh.blk   elgrad.f    flgrad.f    intrp3.f    mkei.f      qadsrc.f    seton0.f    srchs.f     trutbl.f
   aexds1.blk  bs.f        elton0.f    frge.f      intrpe.f    naivekd.c   rda1.f      seton1.f    srcht.f     varept.blk
   aexds2.blk  clsfil.f    elton1.f    header.blk  intrpn.f    naivekd.h   rda2.f      shapef.f    steps.blk   varnpt.blk
   amesh.blk   cntr.f      error.f     hexsrc.f    invcon.f    node.f      rdb1.f      shlsrc.f    stran.f     version.f
   avg.f       contrl.blk  ex2tp.blk   inielt.f    jacobn.f    ntpdat.blk  rdb2.f      sintpe.f    tapes.blk   vol.f
   bldpte.f    debg.blk    exth.f      ininod.f    lenstr.f    opnfil.f    rdinpt.f    sintpn.f    tetsrc.f    wrtc.f
   bldptn.f    debug.f     extq.f      initlz.f    optkd.c     rundat.blk  srchh.f     toldat.blk)

PACKAGE_ADD_LIBRARY(mapvar-kd_app_lib SOURCES ${SOURCES} DEPLIBS supes suplib exodus_for)
PACKAGE_ADD_EXECUTABLE(mapvar-kd NOEXEPREFIX NOEXESUFFIX HEADERS ${HEADERS} SOURCES mapvar-kd.f DEPLIBS mapvar-kd_app_lib COMM serial)
install_executable(mapvar-kd)

SUBPACKAGE_POSTPROCESS()

