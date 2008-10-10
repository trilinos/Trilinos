
INCLUDE(TrilinosPackageHelpers)
INCLUDE(Trilinos_Add_Option)

#
# A) Define the package
#

TRILINOS_PACKAGE(EpetraExt)

#
# B) Set up package-specific options
#

TRILINOS_ADD_OPTION(${PACKAGE_NAME}_BUILD_TRANSFORM
  HAVE_TRANSFORM
  "Enable transform functionality."
  ON )

TRILINOS_ADD_OPTION(${PACKAGE_NAME}_BUILD_GRAPH_REORDERINGS
  HAVE_GRAPH_REORDERINGS
  "Enable experimental RCM and AMD graph reorderings."
   OFF )
   
TRILINOS_ADD_OPTION(${PACKAGE_NAME}_BUILD_EXPERIMENTAL
  HAVE_EXPERIMENTAL
  "Enable experimental features"
  OFF )

TRILINOS_ADD_OPTION(${PACKAGE_NAME}_BUILD_BTF
  HAVE_BTF
  "Enable BTF functionality"
  OFF )

TRILINOS_ADD_OPTION(${PACKAGE_NAME}_BUILD_INOUT
  HAVE_INOUT
  "Enable input/output functionality"
  ON )

TRILINOS_ADD_OPTION(${PACKAGE_NAME}_BUILD_Zoltan
  HAVE_ZOLTAN
  "Enable zoltan functionality"
  OFF )

TRILINOS_ADD_OPTION(${PACKAGE_NAME}_BUILD_COLORING
  HAVE_COLORING
  "Enable coloring functionality"
  ON )

TRILINOS_ADD_OPTION(${PACKAGE_NAME}_BUILD_MATLAB
  HAVE_MATLAB
  "Enable Matlab interface"
  OFF )

TRILINOS_ADD_OPTION(${PACKAGE_NAME}_BUILD_MODEL_EVALUATOR
  HAVE_MODEL_EVALUATOR
  "Enable the EpetraExt::ModelEvaluator class and supporting software"
  ON )

TRILINOS_ADD_OPTION(${PACKAGE_NAME}_BUILD_BLOCK
  HAVE_BLOCK
  "Enable support for blocked Epetra objects"
  ON )

TRILINOS_ADD_OPTION(${PACKAGE_NAME}_BUILD_RESTRICT
  HAVE_RESTRICT
  "Enable support for creating restricted Epetra objects"
  ON )

TRILINOS_ADD_OPTION(${PACKAGE_NAME}_USING_PETSC
  HAVE_PETSC
  "Enable Epetra/PETSc adapters"
  OFF )

TRILINOS_ADD_OPTION(${PACKAGE_NAME}_USING_HTF5
  HAVE_EPETRAEXT_HDF5
  "Enable EpetraExt interface support for HDF5. This interface requires an already installed HDF5 library\; the include, the library path and the library name must be specified as well."
  OFF )

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

TRILINOS_PACKAGE_ADD_TEST_DIRECTORIES(test)

#
# D) Do standard postprocessing
#

TRILINOS_PACKAGE_POSTPROCESS()
