INCLUDE(TribitsSubPackageMacros)
INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsLibraryMacros)

TRIBITS_SUBPACKAGE(Nemspread)

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES el_check_monot.c      el_elm_info.c   el_exoII_io.c
	    el_exoII_io_util.c    pe_check_inp.c
	    pe_init_info.c        pe_input.c      pe_load_lb_info.c
	    pe_str_util.c         pe_util.c       pe_write_parExo_info.c
	    ps_restart.c          rf_allo.c       rf_comm.c
	    rf_mesh_param.c       rf_pre_proc.c   rf_util.c
	    sort_utils.c          md_timer_getrusage.c  md_wrap_generic_c.c)

TRIBITS_ADD_LIBRARY(nem_spread_app_lib SOURCES ${SOURCES} DEPLIBS exodus nemesis suplib)
TRIBITS_ADD_EXECUTABLE(nem_spread NOEXEPREFIX NOEXESUFFIX SOURCES pe_exoII_2par_exoII.c DEPLIBS nem_spread_app_lib COMM serial mpi)
install_executable(nem_spread)

TRIBITS_SUBPACKAGE_POSTPROCESS()

