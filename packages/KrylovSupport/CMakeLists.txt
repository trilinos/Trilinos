
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

# The "KrylovSupport" package will contain essential interfaces and
# utilities for the Krylov subspace methods implemented in Anasazi
# (for eigenvalue problems) and Belos (for linear systems).  It will
# supersede MultiVecTraits, OperatorTraits, and the OrthoManager and
# MatOrthoManager classes and subclasses, which are currently
# duplicated code in both packages.  It will also contain new kernel
# interfaces (such as the matrix powers kernel) and specializations of
# those interfaces for linear algebra libraries such as Tpetra.
# Furthermore, it will contain additional support code that performs
# numerically sensitive computations useful for both Anasazi and
# Belos, relating to the formation of Krylov subspace bases (for
# example, the computation of a Leja ordering of a sequence of
# Newton-basis shifts in the complex plane).
PACKAGE(KrylovSupport)

#
# B) Set up package-specific options
#

# Krylov support includes traits classes and wrappers for various
# Epetra objects.  These are enabled by default when building Trilinos
# with Epetra.
ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_Epetra
  HAVE_KRYLOV_EPETRA
  "Enable Epetra support in the KrylovSupport package."
  "${Trilinos_ENABLE_Epetra}"
  )

# Krylov support includes traits classes and wrappers for various
# Tpetra objects.  These are enabled by default when building Trilinos
# with Tpetra.
ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_Tpetra
  HAVE_KRYLOV_TPETRA
  "Enable Tpetra support in the KrylovSupport package."
  "${Trilinos_ENABLE_Tpetra}"
  )

# Krylov support includes traits classes and wrappers for various
# Thyra objects.  These are enabled by default when building Trilinos
# with Thyra.
ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_Thyra
  HAVE_KRYLOV_TPETRA
  "Enable Thyra support in the KrylovSupport package."
  "${Trilinos_ENABLE_Thyra}"
  )

# The EpetraExt package has some utilities that are helpful for
# testing things like OrthoManager subclasses, but this package can
# still build without EpetraExt support.
ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_EpetraExt
  HAVE_KRYLOV_EPETRAEXT
  "Enable EpetraExt support in the KrylovSupport package."
  "${Trilinos_ENABLE_EpetraExt}"
  )

# The Triutils package has some utilities that are helpful for
# testing, but this package can still build without Triutils support.
ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_Triutils
  HAVE_KRYLOV_TRIUTILS
  "Enable Triutils support in the KrylovSupport package."
  "${Trilinos_ENABLE_Triutils}"
  )

# Kokkos includes an implementation of the Tall Skinny QR (TSQR)
# factorization.  TSQR is used to implement an OrthoManager
# (Krylov::TsqrOrthoManager) and a MatOrthoManager (which only uses
# TSQR in the Euclidean-inner-product case).  TSQR is enabled by
# default if enabled in Kokkos.  If Teuchos_ENABLE_COMPLEX is enabled,
# TSQR is built with complex arithmetic support.
ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_TSQR
  HAVE_KRYLOV_TSQR
  "Enable TSQR (Tall Skinny QR factorization) support in KrylovSupport."
  "${Kokkos_ENABLE_TSQR}" 
  )

# Enable timers in the various computational routines.  We use
# Teuchos::TimeMonitor to wrap the raw timers.
ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_TIMERS
  HAVE_KRYLOV_TIMERS
  "Build with timers enabled."
  NO
  )

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_TEST_DIRECTORIES(test)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
