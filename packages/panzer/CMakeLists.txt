INCLUDE(CombinedOption)

#
# A) Define the package
#

TRIBITS_PACKAGE_DECL(Panzer)

#
# B) Set up package-specific options
#

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEBUG
  PANZER_DEBUG
  "Build with debug code enabled."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_FADTYPE
  PANZER_FADTYPE
  "Choose the Sacado automatic differentiation scalart type (default is DFad)."
  "Sacado::Fad::DFad<RealType>" )

TRIBITS_ADD_ENABLE_TEUCHOS_TIME_MONITOR_OPTION()

IF(Teuchos_ENABLE_LONG_LONG_INT AND  Tpetra_INST_INT_LONG_LONG)
   GLOBAL_SET(PANZER_HAVE_LONG_LONG_INT ON)
ELSE()
   GLOBAL_SET(PANZER_HAVE_LONG_LONG_INT OFF)
ENDIF()

# Explicit Template Instantiation

# ALWAYS FORCE ETI ON FOR PANZER!  ETI in Panzer is for evaluation
# types.  For now we set flags, until we can remove all the ETI code
# for evaluation types.
# TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()
# Copy of above function with hard coded default of ON. The above
# function is hard coded to the Project level default.
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION
  HAVE_${PACKAGE_NAME_UC}_EXPLICIT_INSTANTIATION
  "Enable the use of explicit template instantiation."
  ON
  )
SET(${PARENT_PACKAGE_NAME}_TEMPLATE_DECL_EXT "_decl.hpp")
SET(${PARENT_PACKAGE_NAME}_TEMPLATE_DEF_EXT "_impl.hpp")

TRIBITS_PROCESS_SUBPACKAGES()

TRIBITS_PACKAGE_DEF()

TRIBITS_PACKAGE_POSTPROCESS()
