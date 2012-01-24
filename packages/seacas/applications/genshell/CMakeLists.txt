INCLUDE(TribitsSubPackageMacros)
INCLUDE(TribitsAddExecutable)
INCLUDE(FortranSettings)

TRIBITS_SUBPACKAGE(Genshell)

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES gs_arcxyz.f gs_comand.f gs_dbieb1.f gs_dbielb.f gs_dbmir1.f gs_dbmirr.f gs_dboelb.f gs_dboqa.f
	    gs_dbpini.f gs_exparc.f gs_felcen.f gs_fnpcen.f gs_getint.f gs_getpro.f
	    gs_getspl.f gs_getspt.f gs_idblnk.f gs_inigrd.f gs_ixlink.f gs_makcol.f gs_makrow.f gs_mirss.f
	    gs_mrkcen.f gs_newel1.f gs_newess.f gs_newid.f  gs_newini.f gs_newmap.f gs_newnps.f gs_newxyz.f
	    gs_proxyz.f gs_rdelb.f  gs_rdspln.f gs_relblk.f gs_renumb.f gs_rotxyz.f gs_show.f gs_shwint.f
	    gs_spline.f gs_splxyz.f gs_sptxyz.f gs_trnxyz.f gs_twixyz.f gs_usblk.f gs_usids.f gs_wrelb.f
	    gs_wress.f  gs_wrnps.f gs_wrpxyz.f)

TRIBITS_ADD_LIBRARY(genshell_app_lib SOURCES ${SOURCES} DEPLIBS supes suplib exodus_for)
TRIBITS_ADD_EXECUTABLE(genshell NOEXEPREFIX NOEXESUFFIX LINKER_LANGUAGE Fortran SOURCES genshell.f DEPLIBS genshell_app_lib COMM serial mpi)

install_executable(genshell)

TRIBITS_SUBPACKAGE_POSTPROCESS()


