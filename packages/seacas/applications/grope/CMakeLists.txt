INCLUDE(TribitsSubPackageMacros)
INCLUDE(TribitsAddExecutable)
INCLUDE(FortranSettings)

TRIBITS_SUBPACKAGE(Grope)

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES 
   gr_chkcnt.f gr_ckess.f  gr_dbsbel.f gr_nxmmax.f gr_prelem.f gr_prinit.f gr_prnode.f
   gr_prssv.f  gr_prxyz.f  gr_rdess.f  gr_rdname.f gr_rdstep.f gr_rixwrd.f gr_setprc.f
   gr_chkrng.f gr_ckmap.f  gr_dbselb.f gr_invcon.f gr_preb1.f  gr_press.f  gr_prmap.f
   gr_prnps.f  gr_prstep.f gr_rdcord.f gr_rdmap.f  gr_rdnmap.f gr_rdtims.f gr_selblk.f
   gr_tostep.f gr_ckeb1.f  gr_cknps.f  gr_getprc.f gr_limits.f gr_prebn.f  gr_prfrm.f
   gr_prmmax.f gr_prnsv.f  gr_prtims.f gr_rdeb1.f  gr_rdmmax.f gr_rdnps.f  gr_rixid.f
   gr_selmap.f gr_wdberr.f gr_ckelb.f  gr_comand.f gr_prelb.f  gr_prglob.f gr_selset.f
   gr_prname.f gr_prqa.f   gr_prvers.f gr_rdelb.f  gr_rdnam1.f gr_rdqa.f   gr_rixint.f)

TRIBITS_ADD_LIBRARY(grope_app_lib SOURCES ${SOURCES} DEPLIBS supes suplib exodus_for)
TRIBITS_ADD_EXECUTABLE(grope NOEXEPREFIX NOEXESUFFIX LINKER_LANGUAGE Fortran SOURCES grope.f DEPLIBS grope_app_lib COMM serial mpi)

install_executable(grope)

TRIBITS_SUBPACKAGE_POSTPROCESS()

