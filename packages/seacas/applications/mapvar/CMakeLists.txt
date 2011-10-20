INCLUDE(SubPackageMacros)
INCLUDE(PackageAddExecutable)
INCLUDE(FortranSettings)

SUBPACKAGE(Mapvar)

SET(CMAKE_EXECUTABLE_SUFFIX "")

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES 
	    getbnd.f  ininod.f  mkrnk.f   rdinpt.f  srchge.f  srchh.f   srchs.f   version.f
	    inielt.f  initlz.f  mklstv.f  rank.f    sintpe.f  srchgt.f  srchq.f   srcht.f)

PACKAGE_ADD_LIBRARY(mapvar_app_lib SOURCES ${SOURCES} DEPLIBS mapvarlib supes suplib exodus_for)
PACKAGE_ADD_EXECUTABLE(mapvar NOEXEPREFIX NOEXESUFFIX LINKER_LANGUAGE Fortran SOURCES mapvar.f DEPLIBS mapvar_app_lib COMM serial)
install_executable(mapvar)

SUBPACKAGE_POSTPROCESS()

