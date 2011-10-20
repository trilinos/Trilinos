INCLUDE(SubPackageMacros)
INCLUDE(PackageAddExecutable)
INCLUDE(FortranSettings)

SUBPACKAGE(Mapvar-kd)

SET(CMAKE_EXECUTABLE_SUFFIX "")

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES 
	    debug.f     ininod.f    sintpe.f    srchq.f     srcht.f     optkd.c
	    inielt.f    initlz.f    rdinpt.f    srchh.f     srchs.f     version.f)

PACKAGE_ADD_LIBRARY(mapvar-kd_app_lib SOURCES ${SOURCES} DEPLIBS mapvarlib supes suplib exodus_for)
PACKAGE_ADD_EXECUTABLE(mapvar-kd NOEXEPREFIX NOEXESUFFIX LINKER_LANGUAGE Fortran SOURCES mapvar-kd.f DEPLIBS mapvar-kd_app_lib COMM serial)
install_executable(mapvar-kd)

SUBPACKAGE_POSTPROCESS()

