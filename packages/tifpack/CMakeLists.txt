
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE( Tifpack )

#
# B) Set up package-specific options
#

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_Trilinos
  HAVE_${PACKAGE_NAME}_Trilinos
  "Building as a Trilinos package"
  ON
  )

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_DEBUG
  HAVE_TIFPACK_DEBUG
  "Enable runtime debug checking."
  ${Trilinos_ENABLE_DEBUG}
  )

IF (TPL_ENABLE_QD)
ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_QD
  HAVE_TIFPACK_QD
  "Enable use of the QD TPL."
  ${TPL_ENABLE_QD}
  )
ENDIF()

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_MPI
  HAVE_TIFPACK_MPI
  "Enable use of the MPI TPL."
  ${Trilinos_ENABLE_MPI}
  )

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(test)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()


