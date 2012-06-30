INCLUDE(TribitsSubPackageMacros)
INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsLibraryMacros)

TRIBITS_SUBPACKAGE(Nemspread)

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES  el_elm_info.C            el_exoII_io.C            pe_check_inp.C
	     pe_input.C               pe_load_lb_info.C
	     pe_str_util.C            pe_util.C                pe_write_parExo_info.C
	     ps_restart.C             rf_allo.C
	     rf_mesh_param.C          rf_pre_proc.C            rf_util.C
	     sort_utils.C	      md_timer_getrusage.C md_wrap_generic_c.C)

TRIBITS_ADD_LIBRARY(nem_spread_app_lib SOURCES ${SOURCES} DEPLIBS exodus nemesis suplib)
TRIBITS_ADD_EXECUTABLE(nem_spread NOEXEPREFIX NOEXESUFFIX SOURCES nem_spread.C DEPLIBS nem_spread_app_lib COMM serial mpi)
install_executable(nem_spread)

TRIBITS_SUBPACKAGE_POSTPROCESS()

