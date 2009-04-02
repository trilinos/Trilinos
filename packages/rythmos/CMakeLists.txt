
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)
INCLUDE(CombinedOption)

#
# A) Define the package
#

PACKAGE(Rythmos)

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
    HAVE_RYTHMOS_DEBUG
    "Turn on Rythmos debugging code."
    ${Trilinos_ENABLE_DEBUG}
    )

#
# GAASP options:
#

#COMBINED_OPTION( ${PACKAGE_NAME}_ENABLE_GAASP
#  DEP_OPTIONS_NAMES Teuchos_ENABLE_BOOST
#  DOCSTR "Enables GAASP adapters in Rythmos for global error estimation and control, requires Boost and its support in Teuchos"
#  OFF
#  )

#COMBINED_OPTION( ${PACKAGE_NAME}_ENABLE_GAASP_INTERNAL
#  DEP_OPTIONS_NAMES ${PACKAGE_NAME}_ENABLE_GAASP
#  DOCSTR "Enables internal TPL libraries for GAASP"
#  )


#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

PACKAGE_ADD_TEST_DIRECTORIES(test)

ADD_SUBDIRECTORY(adapters)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
