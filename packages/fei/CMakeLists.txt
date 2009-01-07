
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(FEI)

#
# B) Set up package-specific options
#

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Petsc
  HAVE_FEI_PETSC
  "Enable Petsc functionality."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Prometheus
  HAVE_FEI_PROMETHEUS
  "Enable Prometheus functionality."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Hypre
  HAVE_FEI_HYPRE
  "Enable Prometheus functionality."
  OFF )

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(base)
ADD_SUBDIRECTORY(support-Trilinos)
ADD_SUBDIRECTORY(test_utils)

IF(FEI_ENABLE_Petsc AND TPL_ENABLE_MPI)
  MESSAGE("PETSc is enabbled")
  ADD_SUBDIRECTORY(support-PETSc)
ENDIF()


PACKAGE_ADD_TEST_DIRECTORIES(test)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
