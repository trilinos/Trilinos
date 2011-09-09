
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

# 
#   Exception added to deal with fixed form Fortran code on XL
#
if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "XL|VisualAge")
     set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qfixed=72 -qxlines")
   endif()

#
# A) Define the package
#

PACKAGE(AztecOO)

#
# B) Set up package-specific options
#
ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_AZLU
  HAVE_AZLU
  "Enables the azlu preconditioner.  Requires the y12m TPL."
  "${${PACKAGE_NAME}_ENABLE_y12m}")

IF (NOT ${PACKAGE_NAME}_ENABLE_y12m AND ${PACKAGE_NAME}_ENABLE_AZLU)
  MESSAGE(FATAL_ERROR "You have enabled the AZLU support in AztecOO with the flag \"AztecOO_ENABLE_AZLU\", but this requires the y12m TPL which is currently disabled.  Please enable this TPL or disable AZLU support.")
ENDIF()

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_TEUCHOS_TIME_MONITOR
  AZ_ENABLE_TIMEMONITOR
  "Add timers to AztecOO GMRES solver."
  NO
  )

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_TEST_DIRECTORIES(test)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

#
# Exclude files for source package.
#

PACKAGE_ARCH_EXCLUDE_AUTOTOOLS_FILES()

PACKAGE_ARCH_EXCLUDE_FILES(
  doc/AZ_capture_matrix_howto.txt
  doc/Aztec2.0
  doc/Aztec2.1
  doc/Managing_conditioning_howto.txt
  doc/UserGuide
  doc/azteclogo.gif
  doc/read_captured_matrix.c
  example/AztecOO_RecursiveCall
  example/Epetra_MsrMatrix_AztecOO
  example/Epetra_MsrMatrix_PowerMethod
  example/IfpackIctAztecOO
  example/IfpackAztecOO
  example/IfpackVbrAztecOO
  example/MLAztecOO
  example/azoo_iterate_hb
  example/aztec_app
  example/aztec_hb
  example/AztecOO/adapt_main.mk
  example/AztecOO/vbr_main.mk
  example/AztecOO_MatlabInput/A.dat
  example/AztecOO_MatlabInput/Ainv.dat
  src/AztecOO_string_maps.txt
  src/AztecOO_string_maps_builder.pl
  src/az_comm_.*
  src/md_timer_intel.c
  src/md_timer_ncube.c
  src/md_timer_sol.c
  src/md_timer_sp2.c
  src/md_timer_sun.c
  src/md_timer_win2000.c
  src/md_wrap_intel_c.c
  src/md_wrap_ncube_c.c
  src/md_wrap_puma_c.c
  src/md_wrap_sp2_c.c
  src/stamp-h.in
  test/scripts/daily/serial/Ex_AztecOO_UserOpUserMat
  )

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
