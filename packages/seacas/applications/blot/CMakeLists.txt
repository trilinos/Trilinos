INCLUDE(TribitsSubPackageMacros)
INCLUDE(TribitsAddExecutable)
INCLUDE(FortranSettings)

TRIBITS_SUBPACKAGE(Blot)

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")

APPEND_GLOB(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.f ${CMAKE_CURRENT_SOURCE_DIR}/*.F)
list(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/blot.f")

TRIBITS_ADD_LIBRARY(
	blotlib
	HEADERS	${HEADERS}
	SOURCES ${SOURCES}
	DEPLIBS supes suplib exodus_for plt
)

TRIBITS_ADD_EXECUTABLE(blot_x11 NOEXEPREFIX NOEXESUFFIX
				LINKER_LANGUAGE Fortran
			        SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/blot.f
			        DEPLIBS blotlib svdi_cgi_x11 
			        COMM serial mpi)
install_executable(blot_x11)

TRIBITS_ADD_EXECUTABLE(blot_xcps NOEXEPREFIX NOEXESUFFIX
				LINKER_LANGUAGE Fortran
			        SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/blot.f
			        DEPLIBS blotlib vdx11cps vdicps_dual 
			        COMM serial mpi)
install_executable(blot_xcps)

TRIBITS_ADD_EXECUTABLE(blot_cps NOEXEPREFIX NOEXESUFFIX
				LINKER_LANGUAGE Fortran
			        SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/blot.f
			        DEPLIBS blotlib svdi_cgi_cps
			        COMM serial mpi)
install_executable(blot_cps)

TRIBITS_SUBPACKAGE_POSTPROCESS()

