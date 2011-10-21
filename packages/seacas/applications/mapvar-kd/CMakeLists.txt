INCLUDE(SubPackageMacros)
INCLUDE(PackageAddExecutable)
INCLUDE(FortranSettings)

SUBPACKAGE(Mapvar-kd)

SET(CMAKE_EXECUTABLE_SUFFIX "")

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES srchq-kd.f  srcht-kd.f   optkd.c  srchh-kd.f   srchs-kd.f   version-kd.f)

PACKAGE_ADD_LIBRARY(mapvar-kd_app_lib SOURCES ${SOURCES} DEPLIBS mapvarlib supes suplib exodus_for)
PACKAGE_ADD_EXECUTABLE(mapvar-kd NOEXEPREFIX NOEXESUFFIX LINKER_LANGUAGE Fortran SOURCES mapvar-kd.f DEPLIBS mapvar-kd_app_lib COMM serial)
install_executable(mapvar-kd)

SUBPACKAGE_POSTPROCESS()

