
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(EpetraExt)

#
# B) Set up package-specific options
#

PACKAGE_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

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

# Requires the Zoltan package 
ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_BUILD_Zoltan
  HAVE_ZOLTAN
  "Enable zoltan functionality"
  OFF )

IF(${PACKAGE_NAME}_ENABLE_Zoltan)
  SET(${PACKAGE_NAME}_BUILD_Zoltan ON)
  SET(HAVE_ZOLTAN ON)
ENDIF()

IF(${PACKAGE_NAME}_BUILD_Zoltan)
  IF(NOT ${PACKAGE_NAME}_ENABLE_Zoltan)
    MESSAGE(FATAL_ERROR "Zoltan support is being built since the flag ${PACKAGE}_BUILD_Zoltan was enabled, but the Trilinos package Zoltan was not enabled.  Please enable the Zoltan package.")
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

#Adding a define for when the TPL HYPRE is enabled.
#This replaces an option that was reusing an existing name.
SET(HAVE_HYPRE "${${PACKAGE_NAME}_ENABLE_HYPRE}")

IF(${PACKAGE_NAME}_ENABLE_HDF5)
  SET(${PACKAGE_NAME}_ENABLE_HDF5_DEFAULT ON)
ELSE()
  SET(${PACKAGE_NAME}_ENABLE_HDF5_DEFAULT OFF)
ENDIF()

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_USING_HDF5
  HAVE_EPETRAEXT_HDF5
  "Enable EpetraExt interface support for HDF5. This interface requires an already installed HDF5 library\; the include, the library path and the library name must be specified as well."
  ${${PACKAGE_NAME}_ENABLE_HDF5_DEFAULT} )

#Adding a define for when PyTrilinos is enabled.
#This replaces an option that was reusing an existing name.
SET(HAVE_PYTRILINOS "${Trilinos_ENABLE_PyTrilinos}")

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_TEST_DIRECTORIES(test)

#
# Exclude files for source package.
#

PACKAGE_ARCH_EXCLUDE_AUTOTOOLS_FILES()

PACKAGE_ARCH_EXCLUDE_FILES(
  doc/UserGuide
  doc/inout
  doc/matlab.README
  example/MapColoring/sample_map
  example/MapColoring/sample_matrix
  example/inout/build
  example/model_evaluator/GLpApp/Parallel2DMeshGeneratorFormat.pdf
  example/model_evaluator/GLpApp/README
  example/model_evaluator/GLpApp/data/square/generate-serial-meshes-1-2.sh
  example/model_evaluator/GLpApp/data/square/generate-serial-meshes-1-2.sh.out
  example/model_evaluator/GLpApp/data/square/square.1.2
  example/model_evaluator/GLpApp/data/square/square.1.2.000
  example/model_evaluator/GLpApp/data/square/square.1.2.001
  example/model_evaluator/GLpApp/data/square/square.1.2.edge
  example/model_evaluator/GLpApp/data/square/square.1.2.ele
  example/model_evaluator/GLpApp/data/square/square.1.2.epart.2
  example/model_evaluator/GLpApp/data/square/square.1.2.node
  example/model_evaluator/GLpApp/data/square/square.1.2.npart.2
  example/model_evaluator/GLpApp/data/square/square.1.2.poly
  example/model_evaluator/GLpApp/data/square/square.1.edge
  example/model_evaluator/GLpApp/data/square/square.1.ele
  example/model_evaluator/GLpApp/data/square/square.1.node
  example/model_evaluator/GLpApp/data/square/square.1.poly
  example/model_evaluator/GLpApp/data/square/square.2.edge
  example/model_evaluator/GLpApp/data/square/square.2.ele
  example/model_evaluator/GLpApp/data/square/square.2.node
  example/model_evaluator/GLpApp/data/square/square.2.poly
  example/model_evaluator/GLpApp/data/square/square.poly
  example/model_evaluator/GLpApp/from-triangle-to-serial-input-mesh.pl
  src/btf/pothen/btf_notes
  src/btf/pothen/main.f
  src/distdir
  src/transform/EpetraExt_Dirichlet_.*
  src/transform/EpetraExt_StaticCondensation_LinearProblem..*
  src/transform/EpetraExt_SubCopy_CrsMatrix..*
  src/zoltan/EpetraExt_ZoltanMpi.*
  test/Copy
  test/Makefile.template
  test/Zoltan/Dummy
  test/inout/build
  test/testAll.*
  )

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
