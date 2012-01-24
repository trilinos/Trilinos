INCLUDE(TribitsSubPackageMacros)
INCLUDE(TribitsAddExecutable)
INCLUDE(FortranSettings)

TRIBITS_SUBPACKAGE(Grepos)

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")

SET(SOURCES gp_centroid.f gp_chkss.f  gp_comand.f gp_comid.f
	    gp_dbieb1.f   gp_dbielb.f gp_dbinam.f gp_dbinm1.f gp_dbist1.f
	    gp_dbist2.f   gp_dbiste.f gp_dbmir1.f gp_dbmirr.f gp_dboelb.f
	    gp_dbonam.f   gp_dboqa.f  gp_dboste.f gp_dbpelb.f gp_dbpini.f gp_deform.f
	    gp_delid.f    gp_dobnd.f  gp_donrm2.f gp_donrm3.f gp_elementize.f
	    gp_exfcon.f   gp_expcon.f gp_expxyz.f gp_expxyzn.f gp_getall.f
	    gp_getattnam.f gp_getspl.f gp_getss.f gp_hunt.f   gp_incid.f  gp_intcnt.f
	    gp_limits.f   gp_mapev.f  gp_mapvar.f gp_matxyz.f gp_mirss.f  gp_move.f
	    gp_movint.f   gp_movnd2.f gp_movnod.f gp_movrea.f gp_munelb.f
	    gp_muness.f   gp_munnps.f gp_muntt.f  gp_newatr.f gp_newid.f  gp_newnam.f
	    gp_newxyz.f   gp_ordix.f  gp_ordstr.f gp_putss.f  gp_rdspln.f gp_refnd1.f
	    gp_refnod.f   gp_remap.f  gp_renelb.f gp_renix.f  gp_rotatr.f gp_rotxyz.f
	    gp_sclatr.f   gp_setnor.f gp_setnr2.f gp_show.f   gp_showfl.f gp_smogs.f
	    gp_smogs2.f   gp_smogs3.f gp_snap.f   gp_snpnd2.f gp_snpnod.f gp_spline.f
	    gp_spoff.f    gp_swpess.f gp_swpid.f  gp_vrfyfl.f gp_wrpxyz.f gp_xyzmod.f
	    gp_zess.f     gp_zlnk.f   gp_zmap.f   gp_zmess.f  gp_zmfixd.f gp_zmmap.f
	    gp_zmnps.f    gp_zmxyz.f  gp_znps.f   gp_zxyz.f )

TRIBITS_ADD_LIBRARY(grepos_app_lib SOURCES ${SOURCES} DEPLIBS supes suplib exodus_for)
TRIBITS_ADD_EXECUTABLE(grepos NOEXEPREFIX NOEXESUFFIX LINKER_LANGUAGE Fortran SOURCES grepos.f DEPLIBS grepos_app_lib COMM serial mpi)

install_executable(grepos)

TRIBITS_SUBPACKAGE_POSTPROCESS()

