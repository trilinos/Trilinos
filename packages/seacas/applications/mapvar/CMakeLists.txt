INCLUDE(TribitsSubPackageMacros)
INCLUDE(TribitsAddExecutable)
INCLUDE(FortranSettings)

TRIBITS_SUBPACKAGE(Mapvar)

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES getbnd.f  mkrnk.f   srchge.f  srchh.f   srchs.f   version.f
	    mklstv.f  rank.f    srchgt.f  srchq.f   srcht.f)

TRIBITS_ADD_LIBRARY(mapvar_app_lib SOURCES ${SOURCES} DEPLIBS mapvarlib supes suplib exodus_for)
TRIBITS_ADD_EXECUTABLE(mapvar NOEXEPREFIX NOEXESUFFIX LINKER_LANGUAGE Fortran SOURCES mapvar.f DEPLIBS mapvar_app_lib COMM serial)
install_executable(mapvar)

TRIBITS_SUBPACKAGE_POSTPROCESS()

