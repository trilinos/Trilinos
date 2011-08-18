INCLUDE(SubPackageMacros)
INCLUDE(PackageAddExecutable)
INCLUDE(PackageLibraryMacros)

SUBPACKAGE(Nemspread)

SET(CMAKE_EXECUTABLE_SUFFIX "")

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
	    sort_utils.c          add_to_log.c
            md_timer_getrusage.c  md_wrap_generic_c.c)

PACKAGE_ADD_LIBRARY(nem_spread_app_lib SOURCES ${SOURCES} DEPLIBS exodus nemesis)
PACKAGE_ADD_EXECUTABLE(nem_spread NOEXEPREFIX NOEXESUFFIX SOURCES pe_exoII_2par_exoII.c DEPLIBS nem_spread_app_lib COMM serial)
install_executable(nem_spread)

SUBPACKAGE_POSTPROCESS()

