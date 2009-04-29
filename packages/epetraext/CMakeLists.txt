
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(EpetraExt)

#
# B) Set up package-specific options
#

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_BUILD_TRANSFORM
  HAVE_TRANSFORM
  "Enable transform functionality."
  ON )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_BUILD_GRAPH_REORDERINGS
  HAVE_GRAPH_REORDERINGS
  "Enable experimental RCM and AMD graph reorderings."
   OFF )
   
ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_BUILD_EXPERIMENTAL
  HAVE_EXPERIMENTAL
  "Enable experimental features"
  OFF )

# Experimental and Graph reordering require amd package
IF(${PACKAGE_NAME}_BUILD_EXPERIMENTAL)
  IF(${PACKAGE_NAME}_BUILD_GRAPH_REORDERINGS)
    IF(NOT ${PACKAGE_NAME}_ENABLE_AMD AND NOT ${PACKAGE_NAME}_ENABLE_UMFPACK) 
      MESSAGE(FATAL_ERROR "You have enabled ${PACKAGE_NAME}_BUILD_GRAPH_REORDERINGS and ${PACKAGE_NAME}_BUILD_EXPERIMENTAL options in EpetraExt.  This requires the amd library.  You must enable either the UMFPACK TPL or the AMD TPL to bring in the required library.")
    ENDIF()
  ENDIF()
ENDIF()
   

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_BUILD_BTF
  HAVE_BTF
  "Enable BTF functionality"
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_BUILD_INOUT
  HAVE_INOUT
  "Enable input/output functionality"
  ON )

# Requires the TPL ZoltanTpl or Zoltan package 
ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_BUILD_Zoltan
  HAVE_ZOLTAN
  "Enable zoltan functionality"
  OFF )
IF(${PACKAGE_NAME}_BUILD_Zoltan)
  IF(NOT ${PACKAGE_NAME}_ENABLE_Zoltan AND NOT ${PACKAGE_NAME}_ENABLE_ZoltanTpl)
    MESSAGE(FATAL_ERROR "Zoltan support is being built since the flag ${PACKAGE_NAME}_BUILD_Zoltan was enabled, but neither the TPl ZoltanTpl or the Trilinos package Zoltan was enabled.  Please enable the Zoltan package or the TPL ZoltanTpl.")
  ENDIF()
ENDIF()

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_BUILD_COLORING
  HAVE_COLORING
  "Enable coloring functionality"
  ON )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_BUILD_MATLAB
  HAVE_MATLAB
  "Enable Matlab interface"
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_BUILD_MODEL_EVALUATOR
  HAVE_MODEL_EVALUATOR
  "Enable the EpetraExt::ModelEvaluator class and supporting software"
  ON )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_BUILD_BLOCK
  HAVE_BLOCK
  "Enable support for blocked Epetra objects"
  ON )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_BUILD_RESTRICT
  HAVE_RESTRICT
  "Enable support for creating restricted Epetra objects"
  ON )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_USING_PETSC
  HAVE_PETSC
  "Enable Epetra/PETSc adapters"
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_USING_HDF5
  HAVE_EPETRAEXT_HDF5
  "Enable EpetraExt interface support for HDF5. This interface requires an already installed HDF5 library\; the include, the library path and the library name must be specified as well."
  OFF )

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_TEST_DIRECTORIES(test)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
