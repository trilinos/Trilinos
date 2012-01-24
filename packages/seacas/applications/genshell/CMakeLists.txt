INCLUDE(TribitsSubPackageMacros)
INCLUDE(TribitsAddExecutable)
INCLUDE(FortranSettings)

TRIBITS_SUBPACKAGE(Genshell)

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")

SET(SOURCES gs_dbielb.f gs_dboelb.f gs_getspl.f gs_newid.f gs_rotxyz.f
	    gs_spline.f gs_usids.f gs_wrpxyz.f gs_comand.f gs_dbmir1.f gs_dboqa.f
	    gs_mirss.f gs_newini.f gs_show.f gs_splxyz.f gs_wress.f gs_dbieb1.f
	    gs_dbmirr.f gs_dbpini.f gs_newess.f gs_newxyz.f gs_shwint.f
	    gs_trnxyz.f gs_wrnps.f)


TRIBITS_ADD_LIBRARY(genshell_app_lib SOURCES ${SOURCES} DEPLIBS supes suplib exodus_for)
TRIBITS_ADD_EXECUTABLE(genshell NOEXEPREFIX NOEXESUFFIX LINKER_LANGUAGE Fortran SOURCES genshell.f DEPLIBS genshell_app_lib COMM serial mpi)

install_executable(genshell)

TRIBITS_SUBPACKAGE_POSTPROCESS()


