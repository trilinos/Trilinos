
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE( Anasazi ENABLE_SHADOWING_WARNINGS )

# 2009/05/01: rabartl: I have removed 'CLEANED' from above from now because it
# fails to compile with gcc 3.4.6 (see bug 4395).

#
# B) Set up package-specific options
#

# Anasazi includes tests that may require complex arithmetic in order
# to build.  These include especially (but not exclusively) the
# Tpetra-specific tests, since Tpetra objects are templated on the
# Scalar data type.  We don't want to build tests that depend on
# complex arithmetic if the rest of Trilinos (esp. Teuchos, since we
# need e.g., Teuchos::ScalarTraits< Scalar >) was not built with
# complex arithmetic support.  Thus, we define an appropriate package
# configuration variable ${PACKAGE_NAME}_ENABLE_COMPLEX, which we can
# use to decide whether or not to build particular tests.
#
# The default value of this config variable should be the same as
# Teuchos_ENABLE_COMPLEX, because most of Anasazi's complex-arithmetic
# tests depend heavily on Teuchos' complex arithmetic support.
SET(${PACKAGE_NAME}_ENABLE_COMPLEX_DEFAULT ${Teuchos_ENABLE_COMPLEX})

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_COMPLEX
  HAVE_ANASAZI_COMPLEX
  "Enable Anasazi features for complex data types."
  ${${PACKAGE_NAME}_ENABLE_COMPLEX_DEFAULT} )

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)


IF (Anasazi_ENABLE_Epetra)
  ADD_SUBDIRECTORY(epetra)
ENDIF()

PACKAGE_ADD_TEST_DIRECTORIES(test)

IF (Anasazi_ENABLE_Thyra)
  ADD_SUBDIRECTORY(thyra)
ENDIF()


IF (Anasazi_ENABLE_Tpetra)
  ADD_SUBDIRECTORY(tpetra)
ENDIF()

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
