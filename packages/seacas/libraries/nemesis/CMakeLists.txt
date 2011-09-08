INCLUDE(SubPackageMacros)
INCLUDE(PackageLibraryMacros)
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(FortranSettings)

SUBPACKAGE(Nemesis)

SET(NEMESIS_VERSION_MAJOR "3")
SET(NEMESIS_VERSION_MINOR "09")
SET(NEMESIS_VERSION_PATCH "0")
SET(NEMESIS_VERSION "${NEMESIS_VERSION_MAJOR}.${NEMESIS_VERSION_MINOR}")
SET(NEMESIS_VERSION_FULL "${NEMESIS_VERSION}.${NEMESIS_VERSION_PATCH}")

SET(HEADERS "ne_nemesisI.h")

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/nemesis_cfg.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/nemesis_cfg.h
  @ONLY IMMEDIATE)

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_BINARY_DIR}"
)

SET(SOURCES
ne_gcmp.c
ne_gebig.c
ne_gelcm.c
ne_gelm.c
ne_gelt.c
ne_gennm.c
ne_genv.c
ne_giinf.c
ne_ginig.c
ne_glbp.c
ne_gncm.c
ne_gncor.c
ne_gneat.c
ne_gnelc.c
ne_gnm.c
ne_gnnnm.c
ne_gnns.c
ne_gnnsd.c
ne_gnnv.c
ne_gnspg.c
ne_gnss.c
ne_gnssd.c
ne_gsspg.c
ne_pcmp.c
ne_pcmpc.c
ne_pebig.c
ne_pelcm.c
ne_pelm.c
ne_pennm.c
ne_pevs.c
ne_piinf.c
ne_pinig.c
ne_plbp.c
ne_plbpc.c
ne_pncm.c
ne_pncor.c
ne_pneat.c
ne_pnelc.c
ne_pnm.c
ne_pnnnm.c
ne_pnns.c
ne_pnnsd.c
ne_pnspg.c
ne_pnss.c
ne_pnssd.c
ne_pnvs.c
ne_psspg.c
ne_util.c
)

PACKAGE_ADD_LIBRARY(
	nemesis
	HEADERS	${HEADERS}
	SOURCES ${SOURCES}
	DEPLIBS exodus
)

PACKAGE_ADD_TEST(
	ne_test
	NAME ne_test
	NUM_MPI_PROCS 1
)

PACKAGE_ADD_EXECUTABLE(
 ne_test
 SOURCES ne_test.c
)

SUBPACKAGE_POSTPROCESS()

