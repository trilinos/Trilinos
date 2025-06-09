INCLUDE(TrilinosCreateClientTemplateHeaders)
INCLUDE(ETI_functions)


#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.h)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)
IF (NOT Tpetra_ENABLE_DEPRECATED_CODE OR NOT TpetraCore_ENABLE_Epetra)
  LIST(REMOVE_ITEM SOURCES ${DIR}/Epetra_TsqrMessenger.cpp)
  LIST(REMOVE_ITEM SOURCES ${DIR}/Tpetra_EpetraRowMatrix.cpp)
ENDIF()
IF (NOT HAVE_TPETRACORE_MPI)
  LIST(REMOVE_ITEM HEADERS ${DIR}/Tpetra_Details_Ialltofewv.hpp)
  LIST(REMOVE_ITEM SOURCES ${DIR}/Tpetra_Details_Ialltofewv.cpp)
ENDIF()
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

# Pull in the Kokkos refactor code.
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/kokkos_refactor)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

# Must glob the binary dir last to get all of the auto-generated headers
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_SET(HEADERS ${DIR}/${PACKAGE_NAME}_config.h )
APPEND_SET(HEADERS ${DIR}/${PACKAGE_NAME}_ETIHelperMacros.h)

IF (${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION)
  SET (MultiVector_ETI_SCALARS ${TpetraCore_ETI_SCALARS})
  IF (NOT Tpetra_INST_INT_INT)
    # GO = int is disabled, so the Scalar = GlobalOrdinal case doesn't
    # cover Scalar = int.  We need Scalar = int for communication of
    # process ranks.
    LIST(APPEND MultiVector_ETI_SCALARS "int")
  ENDIF()

  # Generate ETI .cpp files for Tpetra::MultiVector.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(MULTIVECTOR_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "MultiVector" "MULTIVECTOR" "${MultiVector_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${MULTIVECTOR_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::Vector.
  # Use the same set of Scalar types as with MultiVector.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(VECTOR_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "Vector" "VECTOR" "${MultiVector_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${VECTOR_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::FEMultiVector.
  # Use the same set of Scalar types as with MultiVector.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(FEMULTIVECTOR_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "FEMultiVector" "FEMULTIVECTOR" "${MultiVector_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${FEMULTIVECTOR_OUTPUT_FILES})

  # Tpetra::CrsMatrix ETI Support
  ASSERT_DEFINED(Tpetra_ENABLE_Reduced_ETI)
  IF(Tpetra_ENABLE_Reduced_ETI)
    # For Reduced_ETI we don't include int
    MESSAGE(STATUS "Tpetra: Using Reduced ETI")
    SET (CrsMatrix_ETI_SCALARS  ${TpetraCore_ETI_SCALARS})
  ELSE()
    # Zoltan2 wants Scalar = int for CrsMatrix (Bug 6298), and thus
    # indirectly for RowMatrix as well.
    MESSAGE(STATUS "Tpetra: Using Full ETI")
    SET (CrsMatrix_ETI_SCALARS ${MultiVector_ETI_SCALARS})
  ENDIF()

  # Generate ETI .cpp files for Tpetra::CrsMatrix.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(CRSMATRIX_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "CrsMatrix" "CRSMATRIX" "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${CRSMATRIX_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::FECrsMatrix.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(FECRSMATRIX_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "FECrsMatrix" "FECRSMATRIX" "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${FECRSMATRIX_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::RowMatrix.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(ROWMATRIX_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "RowMatrix" "ROWMATRIX" "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${ROWMATRIX_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::RowMatrixTransposer.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(ROWMATRIXTRANSPOSER_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "RowMatrixTransposer" "ROWMATRIXTRANSPOSER" "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${ROWMATRIXTRANSPOSER_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::CrsGraphTransposer.
  TPETRA_PROCESS_ALL_LGN_TEMPLATES(CRSGRAPHTRANSPOSER_OUTPUT_FILES "Tpetra_ETI_LO_GO_NT.tmpl" "CrsGraphTransposer" "CRSGRAPHTRANSPOSER" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}")
  LIST(APPEND SOURCES ${CRSGRAPHTRANSPOSER_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::BlockCrsMatrix.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(BLOCKCRSMATRIX_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "BlockCrsMatrix" "BLOCKCRSMATRIX" "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${BLOCKCRSMATRIX_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::BlockMultiVector.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(BLOCKMULTIVECTOR_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "BlockMultiVector" "BLOCKMULTIVECTOR" "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE )
  LIST(APPEND SOURCES ${BLOCKMULTIVECTOR_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::LinearProblem.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(LINEARPROBLEM_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "LinearProblem" "LINEARPROBLEM" "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${LINEARPROBLEM_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::Reindex_LinearProblem.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(REINDEXLINEARPROBLEM_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "Reindex_LinearProblem" "REINDEXLINEARPROBLEM" "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${REINDEXLINEARPROBLEM_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::Reindex_CrsMatrix.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(REINDEXCRSMATRIX_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "Reindex_CrsMatrix" "REINDEXCRSMATRIX" "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${REINDEXCRSMATRIX_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::Reindex_MultiVector.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(REINDEXMULTIVECTOR_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "Reindex_MultiVector" "REINDEXMULTIVECTOR" "${MultiVector_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${REINDEXMULTIVECTOR_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::SolverMap_LinearProblem.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(SOLVERMAPLINEARPROBLEM_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "SolverMap_LinearProblem" "SOLVERMAPLINEARPROBLEM" "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${SOLVERMAPLINEARPROBLEM_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::SolverMap_CrsMatrix.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(SOLVERMAPCRSMATRIX_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "SolverMap_CrsMatrix" "SOLVERMAPCRSMATRIX" "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${SOLVERMAPCRSMATRIX_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::CrsSingletonFilter_LinearProblem.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(CRSSINGLETONFILTER_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "CrsSingletonFilter_LinearProblem" "CRSSINGLETONFILTER" "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${CRSSINGLETONFILTER_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::BlockVector.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(BLOCKVECTOR_OUTPUT_FILES
    "Tpetra_ETI_SC_LO_GO_NT.tmpl" "BlockVector" "BLOCKVECTOR"
    "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}"
    "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${BLOCKVECTOR_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::BlockCrsMatrix helpers.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(BLOCKCRSMATRIX_HELPER_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl"
  "BlockCrsMatrix_Helpers" "BLOCKCRSMATRIX_HELPERS" "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${BLOCKCRSMATRIX_HELPER_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::Details::Random
  TPETRA_PROCESS_ALL_N_TEMPLATES(RANDOM_OUTPUT_FILES "Tpetra_ETI_NT.tmpl" "Details_Random" "DETAILS_RANDOM" "${TpetraCore_ETI_NODES}")
  LIST(APPEND SOURCES ${RANDOM_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::Details::getDiagCopyWithoutOffsets.
  # Do so for the same set of Scalar types as CrsMatrix (see above),
  # because this function is an implementation detail of CrsMatrix.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(GETDIAGCOPYWITHOUTOFFSETS_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "Details_getDiagCopyWithoutOffsets" "DETAILS_GETDIAGCOPYWITHOUTOFFSETS" "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${GETDIAGCOPYWITHOUTOFFSETS_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::Details::packCrsMatrix and
  # Tpetra::Details::packCrsMatrixWithOwningPIDs.  Do so for the same
  # set of Scalar types as CrsMatrix (see above), because this
  # function is an implementation detail of CrsMatrix.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(PACKCRSMATRIX_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "Details_packCrsMatrix" "DETAILS_PACKCRSMATRIX" "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${PACKCRSMATRIX_OUTPUT_FILES})

  # Generate ETI .cpp files for
  # Tpetra::Details::unpackCrsMatrixAndCombine,
  # Tpetra::Details::unpackAndCombineIntoCrsArrays, and
  # Tpetra::Details::unpackAndCombineWithOwningPIDsCount.  Do so for the same
  # set of Scalar types as CrsMatrix (see above), because this
  # function is an implementation detail of CrsMatrix.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(UNPACKCRSMATRIXANDCOMBINE_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "Details_unpackCrsMatrixAndCombine" "DETAILS_UNPACKCRSMATRIXANDCOMBINE" "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${UNPACKCRSMATRIXANDCOMBINE_OUTPUT_FILES})

  # DistObject requires the same Scalar type instantiations as
  # MultiVector, plus Scalar = char (for CrsMatrix).
  SET (DistObject_ETI_SCALARS ${MultiVector_ETI_SCALARS})
  LIST(APPEND DistObject_ETI_SCALARS "char")

  # Generate ETI .cpp files for Tpetra::DistObject.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(DISTOBJECT_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "DistObject" "DISTOBJECT" "${DistObject_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${DISTOBJECT_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::leftAndOrRightScaleCrsMatrix.
  # Do this only for non-integer Scalar types, since we really only
  # need this function for linear solvers.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(LEFTANDORRIGHTSCALECRSMATRIX_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "leftAndOrRightScaleCrsMatrix" "LEFTANDORRIGHTSCALECRSMATRIX" "${TpetraCore_ETI_SCALARS_NO_ORDS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" FALSE)
  LIST(APPEND SOURCES ${LEFTANDORRIGHTSCALECRSMATRIX_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::computeRowAndColumnOneNorms.
  # Do this only for non-integer Scalar types, since we really only
  # need this function for linear solvers.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(COMPUTEROWANDCOLUMNONENORMS_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "computeRowAndColumnOneNorms" "COMPUTEROWANDCOLUMNONENORMS" "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${COMPUTEROWANDCOLUMNONENORMS_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::replaceDiagonalCrsMatrix.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(REPLACEDIAGONALCRSMATRIX_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "replaceDiagonalCrsMatrix" "REPLACEDIAGONALCRSMATRIX" "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" FALSE)
  LIST(APPEND SOURCES ${REPLACEDIAGONALCRSMATRIX_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::CrsMatrix::convert.
  TPETRA_PROCESS_ALL_CONVERT_TEMPLATES(CRSMATRIX_CONVERT_OUTPUT_FILES "Tpetra_ETI_SOUT_SIN_LO_GO_NT.tmpl" "CrsMatrix" "CRSMATRIX_CONVERT" "${TpetraCore_ETI_SCALARS_NO_ORDS}" "${TpetraCore_ETI_SCALARS_NO_ORDS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" FALSE)
  LIST(APPEND SOURCES ${CRSMATRIX_CONVERT_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::MultiVector::convert.
  TPETRA_PROCESS_ALL_CONVERT_TEMPLATES(MULTIVECTOR_CONVERT_OUTPUT_FILES "Tpetra_ETI_SOUT_SIN_LO_GO_NT.tmpl" "MultiVector" "MULTIVECTOR_CONVERT" "${TpetraCore_ETI_SCALARS_NO_ORDS}" "${TpetraCore_ETI_SCALARS_NO_ORDS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" FALSE)
  LIST(APPEND SOURCES ${MULTIVECTOR_CONVERT_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::LocalCrsMatrixOperator.
  TPETRA_PROCESS_ALL_SN_TEMPLATES(LOCALCRSMATRIXOPERATOR_OUTPUT_FILES
    "Tpetra_ETI_SC_NT.tmpl" "LocalCrsMatrixOperator"
    "LOCALCRSMATRIXOPERATOR" "${CrsMatrix_ETI_SCALARS}"
    "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE TRUE)
  LIST(APPEND SOURCES ${LOCALCRSMATRIXOPERATOR_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::LocalCrsMatrixOperator with mixed Scalars.
  ASSERT_DEFINED (Tpetra_INST_DOUBLE)
  ASSERT_DEFINED (Tpetra_INST_FLOAT)
  IF (Tpetra_INST_DOUBLE AND Tpetra_INST_FLOAT)
    TPETRA_PROCESS_ALL_CONVERT_TEMPLATES(LOCALCRSMATRIXOPERATOR_MIXED_OUTPUT_FILES
      "Tpetra_ETI_SOUT_SIN_LO_GO_NT.tmpl" "LocalCrsMatrixOperator"
      "LOCALCRSMATRIXOPERATOR_MIXED"
      "double" "float"
      "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" True True)
    LIST(APPEND SOURCES ${LOCALCRSMATRIXOPERATOR_MIXED_OUTPUT_FILES})
  ENDIF()
  ASSERT_DEFINED (Tpetra_INST_COMPLEX_DOUBLE)
  ASSERT_DEFINED (Tpetra_INST_COMPLEX_FLOAT)
  IF (Tpetra_INST_COMPLEX_DOUBLE AND Tpetra_INST_COMPLEX_FLOAT)
    TPETRA_PROCESS_ALL_CONVERT_TEMPLATES(LOCALCRSMATRIXOPERATOR_MIXED_OUTPUT_FILES
      "Tpetra_ETI_SOUT_SIN_LO_GO_NT.tmpl" "LocalCrsMatrixOperator"
      "LOCALCRSMATRIXOPERATOR_MIXED"
      "std::complex<double>" "std::complex<float>"
      "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" True True)
    LIST(APPEND SOURCES ${LOCALCRSMATRIXOPERATOR_MIXED_OUTPUT_FILES})
  ENDIF()

  # Generate ETI .cpp files for Tpetra::CrsGraph.
  TPETRA_PROCESS_ALL_LGN_TEMPLATES(CRSGRAPH_OUTPUT_FILES
    "Tpetra_ETI_LO_GO_NT.tmpl"
    "CrsGraph"
    "CRSGRAPH"
    "${TpetraCore_ETI_LORDS}"
    "${TpetraCore_ETI_GORDS}"
    "${TpetraCore_ETI_NODES}")
  LIST(APPEND SOURCES ${CRSGRAPH_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::FECrsGraph.
  TPETRA_PROCESS_ALL_LGN_TEMPLATES(FECRSGRAPH_OUTPUT_FILES
    "Tpetra_ETI_LO_GO_NT.tmpl"
    "FECrsGraph"
    "FECRSGRAPH"
    "${TpetraCore_ETI_LORDS}"
    "${TpetraCore_ETI_GORDS}"
    "${TpetraCore_ETI_NODES}")
  LIST(APPEND SOURCES ${FECRSGRAPH_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::Directory.
  TPETRA_PROCESS_ALL_LGN_TEMPLATES(DIRECTORY_OUTPUT_FILES
    "Tpetra_ETI_LO_GO_NT.tmpl"
    "Directory"
    "DIRECTORY"
    "${TpetraCore_ETI_LORDS}"
    "${TpetraCore_ETI_GORDS}"
    "${TpetraCore_ETI_NODES}")
  LIST(APPEND SOURCES ${DIRECTORY_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::Directory's implementation classes.
  TPETRA_PROCESS_ALL_LGN_TEMPLATES(DIRECTORYIMPL_OUTPUT_FILES
    "Tpetra_ETI_LO_GO_NT.tmpl"
    "DirectoryImpl"
    "DIRECTORYIMPL"
    "${TpetraCore_ETI_LORDS}"
    "${TpetraCore_ETI_GORDS}"
    "${TpetraCore_ETI_NODES}")
  LIST(APPEND SOURCES ${DIRECTORYIMPL_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::Export.
  TPETRA_PROCESS_ALL_LGN_TEMPLATES(EXPORT_OUTPUT_FILES
    "Tpetra_ETI_LO_GO_NT.tmpl"
    "Export"
    "EXPORT"
    "${TpetraCore_ETI_LORDS}"
    "${TpetraCore_ETI_GORDS}"
    "${TpetraCore_ETI_NODES}")
  LIST(APPEND SOURCES ${EXPORT_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::Import.
  TPETRA_PROCESS_ALL_LGN_TEMPLATES(IMPORT_OUTPUT_FILES
    "Tpetra_ETI_LO_GO_NT.tmpl"
    "Import"
    "IMPORT"
    "${TpetraCore_ETI_LORDS}"
    "${TpetraCore_ETI_GORDS}"
    "${TpetraCore_ETI_NODES}")
  LIST(APPEND SOURCES ${IMPORT_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::ImportExportData.
  TPETRA_PROCESS_ALL_LGN_TEMPLATES(IMPORTEXPORTDATA_OUTPUT_FILES
    "Tpetra_ETI_LO_GO_NT.tmpl"
    "ImportExportData"
    "IMPORTEXPORTDATA"
    "${TpetraCore_ETI_LORDS}"
    "${TpetraCore_ETI_GORDS}"
    "${TpetraCore_ETI_NODES}")
  LIST(APPEND SOURCES ${IMPORTEXPORTDATA_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::Details::makeColMap.
  TPETRA_PROCESS_ALL_LGN_TEMPLATES(MAKECOLMAP_OUTPUT_FILES
    "Tpetra_ETI_LO_GO_NT.tmpl"
    "Details_makeColMap"
    "DETAILS_MAKECOLMAP"
    "${TpetraCore_ETI_LORDS}"
    "${TpetraCore_ETI_GORDS}"
    "${TpetraCore_ETI_NODES}")
  LIST(APPEND SOURCES ${MAKECOLMAP_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::Details::packCrsGraph.
  TPETRA_PROCESS_ALL_LGN_TEMPLATES(PACKCRSGRAPH_OUTPUT_FILES
    "Tpetra_ETI_LO_GO_NT.tmpl"
    "Details_packCrsGraph"
    "DETAILS_PACKCRSGRAPH"
    "${TpetraCore_ETI_LORDS}"
    "${TpetraCore_ETI_GORDS}"
    "${TpetraCore_ETI_NODES}")
  LIST(APPEND SOURCES ${PACKCRSGRAPH_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::Details::unpackCrsGraphAndCombine.
  TPETRA_PROCESS_ALL_LGN_TEMPLATES(UNPACKCRSGRAPHANDCOMBINE_OUTPUT_FILES
    "Tpetra_ETI_LO_GO_NT.tmpl"
    "Details_unpackCrsGraphAndCombine"
    "DETAILS_UNPACKCRSGRAPHANDCOMBINE"
    "${TpetraCore_ETI_LORDS}"
    "${TpetraCore_ETI_GORDS}"
    "${TpetraCore_ETI_NODES}")
  LIST(APPEND SOURCES ${UNPACKCRSGRAPHANDCOMBINE_OUTPUT_FILES})

  # Tpetra::Details::FixedHashTable ETI Support
  # This needs GO=int always.
  SET (FixedHashTable_ETI_GORDS ${TpetraCore_ETI_GORDS})
  IF (NOT Tpetra_INST_INT_INT)
    LIST(APPEND FixedHashTable_ETI_GORDS "int")
  ENDIF()

  TPETRA_PROCESS_ALL_LGN_TEMPLATES(FIXEDHASHTABLE_OUTPUT_FILES "Tpetra_ETI_LO_GO_NT.tmpl" "Details_FixedHashTable" "DETAILS_FIXEDHASHTABLE" "${TpetraCore_ETI_LORDS}" "${FixedHashTable_ETI_GORDS}" "${TpetraCore_ETI_NODES}")
  LIST(APPEND SOURCES ${FIXEDHASHTABLE_OUTPUT_FILES})

  #MESSAGE(STATUS "SOURCES = ${SOURCES}")
ENDIF()

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  tpetra
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  ADDED_LIB_TARGET_NAME_OUT TPETRA_LIBNAME
  )

# We need to set the linker language explicitly here for CUDA builds.
SET_PROPERTY(
  TARGET ${TPETRA_LIBNAME}
  APPEND PROPERTY LINKER_LANGUAGE CXX
  )

#
# Make a trivial change to this comment if you add / remove a file either to
# / from this directory, or to / from the 'impl' subdirectory.  That ensures
# that running "make" will also rerun CMake in order to regenerate Makefiles.
#
# Her eis another change.
