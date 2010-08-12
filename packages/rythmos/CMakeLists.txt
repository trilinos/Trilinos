
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)
INCLUDE(CombinedOption)

#
# A) Define the package
#

PACKAGE(Rythmos ENABLE_SHADOWING_WARNINGS CLEANED )

#
# B) Set up package-specific options
#

ADD_OPTION_AND_DEFINE( ${PACKAGE_NAME}_ENABLE_Experimental
    HAVE_RYTHMOS_EXPERIMENTAL
    "Enables Rythmos experimental features"
    OFF
    )

ADD_OPTION_AND_DEFINE( ${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION
    HAVE_RYTHMOS_EXPLICIT_INSTANTIATION
    "Turn on explicit instantiation for Rythmos."
    ${Trilinos_ENABLE_EXPLICIT_INSTANTIATION}
    )

ADD_OPTION_AND_DEFINE( ${PACKAGE_NAME}_ENABLE_DEBUG
    RYTHMOS_DEBUG
    "Turn on Rythmos debugging code."
    ${Trilinos_ENABLE_DEBUG}
    )

#
# GAASP options:
#

# 05/04/09 tscoffe:  You must ask for GAASP support.  If you also provide a
# separate GaaspTpl, then it will be used, otherwise the internal GAASP TPL
# will be used.

OPTION( ${PACKAGE_NAME}_ENABLE_GAASP
    "Enables GAASP adapters in Rythmos for global error estimation and control, requires Boost and its support in Teuchos"
    OFF
    )
#PRINT_VAR( ${PACKAGE_NAME}_ENABLE_GAASP )

# Report error if Boost is not turned on
IF( ${PACKAGE_NAME}_ENABLE_GAASP AND NOT Teuchos_ENABLE_BOOST ) 
  MESSAGE( FATAL_ERROR
      "Error: you can not enable option ${PACKAGE_NAME}_ENABLE_GAASP unless"
      " you also enable the option Teuchos_ENABLE_BOOST.")
ENDIF()

OPTION( ${PACKAGE_NAME}_ENABLE_GAASP_EXTERNAL_TPL
    "Enables Rythmos use of an externally provided GAASP TPL"
    ${TPL_ENABLE_GaaspTpl}
    )
#PRINT_VAR( ${PACKAGE_NAME}_ENABLE_GAASP_EXTERNAL_TPL )

# Report error if Rythmos external GAASP TPL is turned on, but there is no TPL provided
IF( ${PACKAGE_NAME}_ENABLE_GAASP_EXTERNAL_TPL AND NOT TPL_ENABLE_GaaspTpl )
  MESSAGE( FATAL_ERROR
      "Error: you can not enable option ${PACKAGE_NAME}_ENABLE_GAASP_EXTERNAL_TPL unless"
      " you also enable the option TPL_ENABLE_GaaspTpl.")
ENDIF()


SET(${PACKAGE_NAME}_ENABLE_GAASP_INTERNAL_TPL_DEFAULT OFF)
IF( ${PACKAGE_NAME}_ENABLE_GAASP AND NOT ${PACKAGE_NAME}_ENABLE_GAASP_EXTERNAL_TPL )  
  SET(${PACKAGE_NAME}_ENABLE_GAASP_INTERNAL_TPL_DEFAULT ON)
ENDIF()
#PRINT_VAR( ${PACKAGE_NAME}_ENABLE_GAASP_INTERNAL_TPL_DEFAULT )

OPTION( ${PACKAGE_NAME}_ENABLE_GAASP_INTERNAL_TPL
    "Enables Rythmos use of internal GAASP TPL"
    ${${PACKAGE_NAME}_ENABLE_GAASP_INTERNAL_TPL_DEFAULT}
    )
#PRINT_VAR( ${PACKAGE_NAME}_ENABLE_GAASP_INTERNAL_TPL )

# Report error if Boost is not turned on
IF( ${PACKAGE_NAME}_ENABLE_GAASP_INTERNAL_TPL AND NOT Teuchos_ENABLE_BOOST ) 
  MESSAGE( FATAL_ERROR
      "Error: you can not enable option ${PACKAGE_NAME}_ENABLE_GAASP_INTERNAL_TPL unless"
      " you also enable the option Teuchos_ENABLE_BOOST.")
ENDIF()

# Report error if both internal and external GAASP TPLs are turned on
IF( ${PACKAGE_NAME}_ENABLE_GAASP_INTERNAL_TPL AND ${PACKAGE_NAME}_ENABLE_GAASP_EXTERNAL_TPL )
  MESSAGE( FATAL_ERROR
      "Error: you can not enable option ${PACKAGE_NAME}_ENABLE_GAASP_INTERNAL_TPL and"
      " also enable the option ${PACKAGE_NAME}_ENABLE_GAASP_EXTERNAL_TPL.")
ENDIF()
# Report error if no GAASP TPLs are turned on
IF( ${PACKAGE_NAME}_ENABLE_GAASP AND NOT ${PACKAGE_NAME}_ENABLE_GAASP_INTERNAL_TPL AND NOT ${PACKAGE_NAME}_ENABLE_GAASP_EXTERNAL_TPL)
  MESSAGE( FATAL_ERROR
      "Error: you can not enable option ${PACKAGE_NAME}_ENABLE_GAASP unless"
      " you also enable one of ${PACKAGE_NAME}_ENABLE_GAASP_INTERNAL_TPL"
      " or ${PACKAGE_NAME}_ENABLE_GAASP_EXTERNAL_TPL.")
ENDIF()

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

PACKAGE_ADD_TEST_DIRECTORIES(test)

ADD_SUBDIRECTORY(adapters)

#
# Exclude files for source package.
#

PACKAGE_ARCH_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
