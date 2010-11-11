
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

# The "KrylovSupport" package will contain essential interfaces and
# utilities for the Krylov subspace methods implemented in Anasazi
# (for eigenvalue problems) and Belos (for linear systems).  It will
# supersede MultiVecTraits, OperatorTraits, and the orthogonalization
# managers which are currently duplicated code in both packages.  It
# will also contain new kernel interfaces (such as the matrix powers
# kernel) and specializations of those interfaces for linear algebra
# libraries such as Tpetra.  Furthermore, it will contain additional
# support code that performs numerically sensitive computations useful
# for both Anasazi and Belos, relating to the formation of Krylov
# subspace bases (for example, the computation of a Leja ordering of a
# sequence of Newton-basis shifts in the complex plane).
PACKAGE(KrylovSupport)

#
# B) Set up package-specific options
#

#ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_SOMETHING
#  HAVE_${PACKAGE_NAME}_SOMETHING
#  "Enable some functionality."
#  ON )

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_TEST_DIRECTORIES(test)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
