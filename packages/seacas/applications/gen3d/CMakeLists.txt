INCLUDE(TribitsSubPackageMacros)
INCLUDE(TribitsAddExecutable)
INCLUDE(FortranSettings)

TRIBITS_SUBPACKAGE(Gen3D)

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")

SET(SOURCES g3_arcxyz.f   g3_arcyxz.f   g3_comand.f   g3_dbmir1.f   g3_dbmirr.f   g3_dboqa.f
	    g3_dbpini.f   g3_exparc.f   g3_felcen.f   g3_fnpcen.f   g3_getint.f  g3_getpro.f   g3_getspl.f
	    g3_getspt.f   g3_idblnk.f   g3_inigrd.f   g3_ixlink.f   g3_makcol.f   g3_makrow.f  g3_mirss.f    g3_mrkcen.f
	    g3_newel1.f   g3_newess.f   g3_newid.f    g3_newini.f   g3_newmap.f   g3_newnps.f  g3_newxyz.f   g3_proxyz.f
	    g3_rdelb.f    g3_rdspln.f   g3_relblk.f   g3_renumb.f   g3_rotxyz.f   g3_show.f    g3_shwint.f   g3_spline.f
	    g3_splxyz.f   g3_sptxyz.f   g3_trnxyz.f   g3_twixyz.f   g3_usblk.f    g3_usids.f   g3_wrelb.f    g3_wress.f
	    g3_wrnps.f    g3_wrpxyz.f)

TRIBITS_ADD_LIBRARY(gen3d_app_lib SOURCES ${SOURCES} DEPLIBS supes suplib exodus_for)
TRIBITS_ADD_EXECUTABLE(gen3d NOEXEPREFIX NOEXESUFFIX LINKER_LANGUAGE Fortran SOURCES gen3d.f DEPLIBS gen3d_app_lib COMM serial mpi)

install_executable(gen3d)

TRIBITS_SUBPACKAGE_POSTPROCESS()

