INCLUDE(TribitsSubPackageMacros)
INCLUDE(TribitsAddExecutable)
INCLUDE(FortranSettings)

TRIBITS_SUBPACKAGE(Gjoin)

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES gj_prtelb.F gj_prtess.F gj_prtnps.F
	    gj_chgint.f gj_dbmirr.f gj_expxyz.f gj_initin.f gj_matxyz.f gj_movrea.f
	    gj_munnps.f gj_ordix.f  gj_rdgen.f  gj_rotxyz.f gj_wrgen.f  gj_zmmap.f
	    gj_comand.f gj_dbpini.f gj_fixmap.f gj_intcnt.f gj_mirss.f  gj_munelb.f
	    gj_munxyz.f gj_ordstr.f gj_renelb.f gj_setitl.f gj_zmess.f  gj_zmnps.f
	    gj_dbmir1.f gj_dorot.f  gj_getall.f gj_inigen.f gj_irennp.f gj_movint.f
	    gj_muness.f gj_offset.f gj_outlog.f gj_renix.f  gj_setsta.f gj_zmfixd.f gj_zmxyz.f)

TRIBITS_ADD_LIBRARY(gjoin_app_lib SOURCES ${SOURCES} DEPLIBS supes suplib exodus_for)
TRIBITS_ADD_EXECUTABLE(gjoin NOEXEPREFIX NOEXESUFFIX LINKER_LANGUAGE Fortran SOURCES gjoin.f DEPLIBS gjoin_app_lib COMM serial mpi)

install_executable(gjoin)

TRIBITS_SUBPACKAGE_POSTPROCESS()

