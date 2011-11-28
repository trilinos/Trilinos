INCLUDE(TribitsSubPackageMacros)
INCLUDE(TribitsAddExecutable)
INCLUDE(FortranSettings)

TRIBITS_SUBPACKAGE(Mapvar-kd)

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES srchq-kd.f  srcht-kd.f   optkd.c  srchh-kd.f   srchs-kd.f   version-kd.f)

TRIBITS_ADD_LIBRARY(mapvar-kd_app_lib SOURCES ${SOURCES} DEPLIBS mapvarlib supes suplib exodus_for)
TRIBITS_ADD_EXECUTABLE(mapvar-kd NOEXEPREFIX NOEXESUFFIX LINKER_LANGUAGE Fortran SOURCES mapvar-kd.f DEPLIBS mapvar-kd_app_lib COMM serial mpi)
install_executable(mapvar-kd)

TRIBITS_SUBPACKAGE_POSTPROCESS()

