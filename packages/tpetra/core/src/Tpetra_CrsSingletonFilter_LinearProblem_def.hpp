// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_CRSSINGLETONFILTER_LINEARPROBLEM_DEF_HPP
#define TPETRA_CRSSINGLETONFILTER_LINEARPROBLEM_DEF_HPP

/// \file Tpetra_CrsSingletonFilter_LinearProblem_def.hpp
/// \brief Definition of the Tpetra::CrsSingletonFilter_LinearProblem class

#include "Tpetra_CrsSingletonFilter_LinearProblem_decl.hpp"

namespace Tpetra {

//==============================================================================

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    CrsSingletonFilter_LinearProblem(bool verbose)
  : globalNumSingletonRows_(Teuchos::as<local_ordinal_type>(0))
  , globalNumSingletonCols_(Teuchos::as<local_ordinal_type>(0))
  , RatioOfDimensions_(0.0)
  , RatioOfNonzeros_(0.0)
  , HaveReducedProblem_(false)
  , AnalysisDone_(false)
  , SymmetricElimination_(true)
  , localMaxNumRowEntries_(Teuchos::as<local_ordinal_type>(0))
  , FullMatrixIsCrsMatrix_(false)
  , verbose_(verbose) {
}

//==============================================================================

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::NewType
CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
operator()(const CrsSingletonFilter_LinearProblem::OriginalType& originalLinearProblem) {
  analyze(originalLinearProblem);
  return construct();
}

//==============================================================================

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    analyze(const CrsSingletonFilter_LinearProblem::OriginalType& originalLinearProblem) {
  this->origObj_ = originalLinearProblem;

  FullMatrix_ = originalLinearProblem->getMatrix();

  Analyze(FullMatrix());

  if (verbose_ && FullMatrix()->getComm()->getRank() == 0) {
    std::cout << "\nAnalyzed Singleton Problem:\n";
    std::cout << "---------------------------\n";
  }
  if (SingletonsDetected()) {
    if (verbose_ && FullMatrix()->getComm()->getRank() == 0) {
      std::cout << "Singletons Detected!" << std::endl;
      std::cout << "Num Singletons:      " << NumSingletons() << std::endl;
    }
  } else {
    if (verbose_ && FullMatrix()->getComm()->getRank() == 0)
      std::cout << "No Singletons Detected!" << std::endl;
  }

  if (verbose_ && FullMatrix()->getComm()->getRank() == 0)
    std::cout << "---------------------------\n\n";

  return;
}

//==============================================================================

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::NewType
CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    construct() {
  const char tfecfFuncName[] = "construct: ";

  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(this->origObj_ == Teuchos::null, std::runtime_error,
                                        "Linear problem (orgiObj_) have not been set.  Call analyze() first.");

  ConstructReducedProblem(this->origObj_);

  this->newObj_ = ReducedProblem();

  if (verbose_ && SingletonsDetected() && FullMatrix()->getComm()->getRank() == 0) {
    std::cout << "\nConstructedSingleton Problem:\n";
    std::cout << "---------------------------\n";
    std::cout << "RatioOfDimensions:   " << RatioOfDimensions() << std::endl;
    std::cout << "RatioOfNonzeros:     " << RatioOfNonzeros() << std::endl;
    std::cout << "---------------------------\n\n";
  }

  return this->newObj_;
}

//==============================================================================

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    fwd() {
  UpdateReducedProblem(FullProblem());
}

//==============================================================================

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    rvs() {
  // GOT TO HERE!!
  // ComputeFullSolution();
}

//==============================================================================

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Analyze(const Teuchos::RCP<row_matrix_type>& fullMatrix) {
  const char tfecfFuncName[] = "Analyze: ";

  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(AnalysisDone_, std::runtime_error,
                                        "Analyze() already done once.  Cannot do it again.");

  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(fullMatrix == Teuchos::null,
                                        std::runtime_error, "Input matrix is Teuchos::null.");

  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(fullMatrix->getGlobalNumRows() == 0,
                                        std::runtime_error, "Full matrix has zero dimension.");

  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(fullMatrix->getGlobalNumEntries() == 0,
                                        std::runtime_error, "Full matrix has no nonzero terms.");

  FullMatrix_ = fullMatrix;

  // First check for columns with single entries and find columns with singleton rows
  vector_type_int ColProfiles(FullMatrixColMap());
  ColProfiles.putScalar(0);
  vector_type_int ColHasRowWithSingleton(FullMatrixColMap());
  ColHasRowWithSingleton.putScalar(0);

  // localRowIDofSingletonCol[j] will contain the local row ID associated with the jth column,
  // if the jth col has a single entry
  vector_type_LO localRowIDofSingletonCol(FullMatrixColMap());
  localRowIDofSingletonCol.putScalar(-1);

  // Define MapColoring objects
  RowMapColors_ = Teuchos::rcp(new vector_type_int(FullMatrixRowMap()));  // Initial colors are all 0
  ColMapColors_ = Teuchos::rcp(new vector_type_int(FullMatrixColMap()));

  local_ordinal_type localNumRows = FullMatrix()->getLocalNumRows();
  local_ordinal_type localNumCols = FullMatrix()->getLocalNumCols();

  // Set up for accessing full matrix.  Will do so row-by-row.
  InitFullMatrixAccess();

  // Scan matrix for singleton rows, build up column profiles
  size_t NumIndices = 1;
  // int * localIndices;
  // nonconst_local_inds_host_view_type localIndices;
  Teuchos::Array<local_ordinal_type> localIndices;
  localNumSingletonRows_ = 0;

  auto ColProfilesData              = ColProfiles.getLocalViewHost(Tpetra::Access::ReadWrite);
  auto localRowIDofSingletonColData = localRowIDofSingletonCol.getLocalViewHost(Tpetra::Access::ReadWrite);
  auto ColHasRowWithSingletonData   = ColHasRowWithSingleton.getLocalViewHost(Tpetra::Access::ReadWrite);
  auto RowMapColors_Data            = RowMapColors_->getLocalViewHost(Tpetra::Access::ReadWrite);
  auto ColMapColors_Data            = ColMapColors_->getLocalViewHost(Tpetra::Access::ReadWrite);

  for (int i = 0; i < localNumRows; i++) {
    // Get ith row
    GetRow(i, NumIndices, localIndices);
    for (size_t j = 0; j < NumIndices; j++) {
      local_ordinal_type ColumnIndex = localIndices[j];

      // Bounds check for ColumnIndex
      if (static_cast<size_t>(ColumnIndex) >= ColProfilesData.extent(0)) {
        std::cout << "Error: ColumnIndex out of bounds: " << ColumnIndex << std::endl;
        std::abort();
      }

      ColProfilesData(ColumnIndex, 0)++;  // Increment column count

      if (static_cast<size_t>(ColumnIndex) >= localRowIDofSingletonColData.extent(0)) {
        std::cout << "Error: ColumnIndex out of bounds for localRowIDofSingletonColData: "
                  << ColumnIndex << std::endl;
        std::abort();
      }

      // Record local row ID for current column
      // will use to identify row to eliminate if column is a singleton
      localRowIDofSingletonColData(ColumnIndex, 0) = i;
    }
    // If row has single entry, color it and associated column with color=1
    if (NumIndices == 1) {
      int j2 = localIndices[0];
      ColHasRowWithSingletonData(j2, 0)++;
      RowMapColors_Data(i, 0)  = 1;
      ColMapColors_Data(j2, 0) = 1;
      localNumSingletonRows_++;
    }
  }

  // 1) The vector ColProfiles has column nonzero counts for each processor's contribution
  // Combine these to get total column profile information and then redistribute to processors
  // so each can determine if it is the owner of the row associated with the singleton column
  // 2) The vector ColHasRowWithSingleton[i] contain count of singleton rows  that are associated with
  // the ith column on this processor.  Must tell other processors that they should also eliminate
  // these columns.

  // Make a copy of ColProfiles for later use when detecting columns that disappear locally
  vector_type_int NewColProfiles(ColProfiles.getMap());
  NewColProfiles.update(1.0, ColProfiles, 0.0);
  auto NewColProfilesData = NewColProfiles.getLocalViewHost(Tpetra::Access::ReadWrite);

  // If importer is non-trivial, we need to perform a gather/scatter to accumulate results
  auto importer = FullCrsMatrix()->getCrsGraph()->getImporter();
  if (importer != Teuchos::null) {
    vector_type_int tmpVec(FullMatrixDomainMap());  // Use for gather/scatter of column vectors
    tmpVec.putScalar(0);
    tmpVec.doExport(ColProfiles, *importer, ADD);
    ColProfiles.doImport(tmpVec, *importer, INSERT);

    tmpVec.putScalar(0);
    tmpVec.doExport(ColHasRowWithSingleton, *importer, ADD);
    ColHasRowWithSingleton.doImport(tmpVec, *importer, INSERT);
  }

  for (size_t j = 0; j < NumIndices; j++) {
    local_ordinal_type ColumnIndex = localIndices[j];

    // Bounds check for ColumnIndex
    if (static_cast<size_t>(ColumnIndex) >= ColProfilesData.extent(0)) {
      std::cout << "Error: ColumnIndex out of bounds: " << ColumnIndex << std::endl;
      std::abort();
    }
  }
  for (size_t j = 0; j < NumIndices; j++) {
    local_ordinal_type ColumnIndex = localIndices[j];

    // Bounds check for ColumnIndex
    if (static_cast<size_t>(ColumnIndex) >= ColHasRowWithSingletonData.extent(0)) {
      std::cout << "Error: ColumnIndex out of bounds for ColHasRowWithSingletonData: "
                << ColumnIndex << std::endl;
      std::abort();
    }
  }

  // ColProfiles now contains the nonzero column entry count for all columns that have
  // an entry on this processor.
  // ColHasRowWithSingleton now contains a count of singleton rows associated with the corresponding
  // local column.  Next we check to make sure no column is associated with more than one singleton row.

  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(ColHasRowWithSingleton.normInf() > 1, std::runtime_error,
                                        "At least one column is associated with two singleton rows, can't handle it.");

  vector_type_int RowHasColWithSingleton(FullMatrix()->getRowMap());  // Use to check for errors
  RowHasColWithSingleton.putScalar(0);
  auto RowHasColWithSingletonData = RowHasColWithSingleton.getLocalViewHost(Tpetra::Access::ReadWrite);

  localNumSingletonCols_ = 0;
  // Count singleton columns (that were not already counted as singleton rows)
  for (local_ordinal_type j = 0; j < localNumCols; j++) {
    local_ordinal_type i2 = localRowIDofSingletonColData(j, 0);
    // Check if column is a singleton
    if (ColProfilesData(j, 0) == 1) {
      // Check to see if this column already eliminated by the row check above
      if (RowMapColors_Data(i2, 0) != 1) {
        RowHasColWithSingletonData(i2, 0)++;  // Increment col singleton counter for ith row
        RowMapColors_Data(i2, 0) = 2;         // Use 2 for now, to distinguish between row eliminated directly or via column singletons
        ColMapColors_Data(j, 0)  = 1;
        localNumSingletonCols_++;
        // If we delete a row, we need to keep track of associated column entries that were also deleted
        // in case all entries in a column are eventually deleted, in which case the column should
        // also be deleted.
        GetRow(i2, NumIndices, localIndices);
        for (size_t jj = 0; jj < NumIndices; jj++) {
          NewColProfilesData(localIndices[jj], 0)--;
        }
      }
    }
    // Check if some other processor eliminated this column
    else if (ColHasRowWithSingletonData(j, 0) == 1 && RowMapColors_Data(i2, 0) != 1) {
      ColMapColors_Data(j, 0) = 1;
    }
  }

  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(RowHasColWithSingleton.normInf() > 1, std::runtime_error,
                                        "At least one row is associated with two singleton columns, can't handle it.");

  // Generate arrays that keep track of column singleton row, col and pivot info needed for post-solve phase
  CreatePostSolveArrays(localRowIDofSingletonCol, ColProfiles, NewColProfiles, ColHasRowWithSingleton);

  for (local_ordinal_type i = 0; i < localNumRows; i++) {
    if (RowMapColors_Data(i, 0) == 2) RowMapColors_Data(i, 0) = 1;  // Convert all eliminated rows to same color
  }

  const Teuchos::Ptr<local_ordinal_type> gRowsPtr(&globalNumSingletonRows_);
  const Teuchos::Ptr<local_ordinal_type> gColsPtr(&globalNumSingletonCols_);
  auto rowComm = FullMatrix()->getRowMap()->getComm();
  Teuchos::reduceAll<local_ordinal_type, local_ordinal_type>(*rowComm, Teuchos::REDUCE_SUM, localNumSingletonRows_, gRowsPtr);
  Teuchos::reduceAll<local_ordinal_type, local_ordinal_type>(*rowComm, Teuchos::REDUCE_SUM, localNumSingletonCols_, gColsPtr);

  AnalysisDone_ = true;
  return;
}

//==============================================================================

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    GenerateReducedMap(const Teuchos::RCP<const map_type>& originalMap,
                       const Teuchos::RCP<vector_type_int>& mapColors, int color) {
  // Vector to hold the reduced global indices
  std::vector<GlobalOrdinal> myReducedGlobalIndices;

  // Get the global IDs owned by the calling process
  auto myGlobalIndices = originalMap->getMyGlobalIndices();

  // Iterate through the global indices to find the ones that match the color
  for (size_t i = 0; i < myGlobalIndices.size(); i++) {
    // Access the value in mapColors using the appropriate method
    int colorValue = mapColors->getData()[i];  // Use getData() to access the vector data
    if (colorValue == color) {
      myReducedGlobalIndices.push_back(myGlobalIndices[i]);
    }
  }

  // Create the reduced map using the collected indices
  Teuchos::RCP<const map_type> reducedMap = createNonContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(
      Teuchos::ArrayView<const GlobalOrdinal>(myReducedGlobalIndices.data(), myReducedGlobalIndices.size()),
      originalMap->getComm());

  return reducedMap;
}

//==============================================================================

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ConstructReducedProblem(const Teuchos::RCP<linear_problem_type>& Problem) {
  const char tfecfFuncName[] = "ConstructReducedProblem: ";

  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(HaveReducedProblem_, std::runtime_error,
                                        "Already have a reduced problem.  Cannot do it again.");

  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(Problem == Teuchos::null,
                                        std::runtime_error, "Problem is Teuchos::null.");

  FullProblem_ = Problem;
  FullMatrix_  = Teuchos::rcp_dynamic_cast<row_matrix_type>(Problem->getMatrix());
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(FullMatrix() == Teuchos::null,
                                        std::runtime_error, "Need a RowMatrix.");
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(Problem->getRHS() == Teuchos::null,
                                        std::runtime_error, "Need a RHS.");
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(Problem->getLHS() == Teuchos::null,
                                        std::runtime_error, "Need a LHS.");

  // Generate reduced row and column maps
  if (SingletonsDetected()) {
    ReducedMatrixRowMap_ = GenerateReducedMap(FullMatrixRowMap(), RowMapColors_, 0);
    ReducedMatrixColMap_ = GenerateReducedMap(FullMatrixColMap(), ColMapColors_, 0);

    // Create domain and range map colorings by exporting map coloring of column and row maps
    auto importer = FullCrsMatrix()->getCrsGraph()->getImporter();
    if (importer != Teuchos::null) {
      Teuchos::RCP<vector_type_int> DomainMapColors = Teuchos::rcp(new vector_type_int(FullMatrixDomainMap()));
      DomainMapColors->doExport(*ColMapColors_, *importer, Tpetra::ABSMAX);
      OrigReducedMatrixDomainMap_ = GenerateReducedMap(FullMatrixDomainMap(), DomainMapColors, 0);
    } else {
      OrigReducedMatrixDomainMap_ = ReducedMatrixColMap();
    }

    auto exporter = FullCrsMatrix()->getCrsGraph()->getImporter();
    if (FullMatrixIsCrsMatrix_) {
      if (exporter != Teuchos::null) {  // Non-trivial exporter
        Teuchos::RCP<vector_type_int> RangeMapColors = Teuchos::rcp(new vector_type_int(FullMatrixRangeMap()));
        RangeMapColors->doExport(*ColMapColors_, *exporter, Tpetra::ABSMAX);
        ReducedMatrixRangeMap_ = GenerateReducedMap(FullMatrixRangeMap(), RangeMapColors, 0);
      } else
        ReducedMatrixRangeMap_ = ReducedMatrixRowMap();
    } else
      ReducedMatrixRangeMap_ = ReducedMatrixRowMap();

    // Check to see if the reduced system domain and range maps are the same.
    // If not, we need to remap entries of the LHS multivector so that they are distributed
    // conformally with the rows of the reduced matrix and the RHS multivector
    SymmetricElimination_ = ReducedMatrixRangeMap()->isSameAs(*OrigReducedMatrixDomainMap_);
    if (!SymmetricElimination_) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(SymmetricElimination_ != 0,
                                            std::runtime_error, "Symmetric Elimination has not been tested or transitioned from Epetra!");
      // ConstructRedistributeExporter(OrigReducedMatrixDomainMap_, ReducedMatrixRangeMap_,
      //                               RedistributeDomainExporter_, ReducedMatrixDomainMap_);
    } else {
      ReducedMatrixDomainMap_     = OrigReducedMatrixDomainMap_;
      OrigReducedMatrixDomainMap_ = Teuchos::null;
      // RedistributeDomainExporter_ = 0;
    }

    // Create pointer to Full RHS, LHS
    Teuchos::RCP<multivector_type> FullRHS = FullProblem()->getRHS();
    Teuchos::RCP<multivector_type> FullLHS = FullProblem()->getLHS();
    int NumVectors                         = FullLHS->getNumVectors();

    // Create importers
    Full2ReducedLHSImporter_ = Teuchos::rcp(new import_type(FullMatrixDomainMap(), ReducedMatrixDomainMap()));
    Full2ReducedRHSImporter_ = Teuchos::rcp(new import_type(FullRHS->getMap(), ReducedMatrixRowMap()));

    // Construct Reduced Matrix
    LocalOrdinal maxNumEntries = FullCrsMatrix()->getLocalMaxNumRowEntries();
    ReducedMatrix_             = Teuchos::rcp(new crs_matrix_type(ReducedMatrixRowMap(), ReducedMatrixColMap(), maxNumEntries));

    // Create storage for temporary X values due to explicit elimination of rows
    tempExportX_ = Teuchos::rcp(new multivector_type(FullMatrixColMap(), NumVectors));

    size_t NumEntries = 0;
    Teuchos::ArrayView<const Scalar> Values;
    Teuchos::Array<GlobalOrdinal> Indices;
    LocalOrdinal localNumRows        = FullMatrix()->getLocalNumRows();
    LocalOrdinal ColSingletonCounter = 0;

    for (LocalOrdinal i = 0; i < localNumRows; i++) {
      GlobalOrdinal curGRID = FullMatrixRowMap()->getGlobalElement(i);
      if (ReducedMatrixRowMap()->isNodeGlobalElement(curGRID)) {  // Check if this row should go into reduced matrix
        GetRowGCIDs(i, NumEntries, Values, Indices);              // Get current row (Indices are global)

        // The ReducedMatrix has the row and column indices already removed so cannot "insert" them.
        // Need to filter them to only insert remaining.
        // Filter indices and values
        Teuchos::Array<GlobalOrdinal> filteredIndices;
        Teuchos::Array<Scalar> filteredValues;

        for (typename Teuchos::Array<GlobalOrdinal>::size_type j = 0; j < Indices.size(); ++j) {
          if (ReducedMatrixColMap()->isNodeGlobalElement(Indices[j])) {
            filteredIndices.push_back(Indices[j]);
            filteredValues.push_back(Values[j]);
          }
        }

        // Insert filtered values into the matrix
        if (!filteredIndices.empty()) {
          ReducedMatrix()->insertGlobalValues(curGRID, filteredIndices(), filteredValues());
        }

      } else {
        Teuchos::ArrayView<const LocalOrdinal> localIndices;
        GetRow(i, NumEntries, Values, localIndices);  // Get current row
        if (NumEntries == 1) {
          Scalar pivot = Values[0];
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(pivot == 0.0, std::runtime_error,
                                                "Encountered zero row, unable to continue.");  // Should improve this comparison to zero.
          LocalOrdinal indX = localIndices[0];
          for (LocalOrdinal j = 0; j < NumVectors; j++) {
            auto rhsData     = FullRHS->getData(j);               // Get the underlying data for vector j
            auto exportData  = tempExportX_->getDataNonConst(j);  // Get the underlying data for vector j (non-const)
            exportData[indX] = rhsData[i] / pivot;
          }
        } else {  // Singleton column
          LocalOrdinal targetCol = ColSingletonColLIDs_[ColSingletonCounter];
          for (size_t j = 0; j < NumEntries; j++) {
            if (localIndices[j] == targetCol) {
              Scalar pivot = Values[j];
              TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(pivot == 0.0, std::runtime_error,
                                                    "Encountered zero column, unable to continue.");  // Should improve this comparison to zero.
              ColSingletonPivotLIDs_[ColSingletonCounter] = j;
              ColSingletonPivots_[ColSingletonCounter]    = pivot;
              ColSingletonCounter++;
              break;
            }
          }
        }
      }
    }

    // Now convert to local indexing.  We have constructed things so that the domain and range of the
    // matrix will have the same map.  If the reduced matrix domain and range maps were not the same, the
    // differences were addressed in the ConstructRedistributeExporter() method
    ReducedMatrix()->fillComplete(ReducedMatrixDomainMap(), ReducedMatrixRangeMap());

    // 1) The vector ColProfiles has column nonzero counts for each processor's contribution
    // Construct Reduced LHS (Puts any initial guess values into reduced system)

    // Create the ReducedLHS_ MultiVector
    ReducedLHS_ = Teuchos::rcp(new multivector_type(ReducedMatrixDomainMap(), NumVectors));
    ReducedLHS_->doImport(*FullLHS, *Full2ReducedLHSImporter_, Tpetra::INSERT);

    FullLHS->putScalar(0.0);

    // Construct Reduced RHS

    // First compute influence of already-known values of X on RHS
    tempX_ = Teuchos::rcp(new multivector_type(FullMatrixDomainMap(), NumVectors));
    tempB_ = Teuchos::rcp(new multivector_type(FullRHS->getMap(), NumVectors));

    // Inject known X values into tempX for purpose of computing tempB = FullMatrix*tempX
    //  Also inject into full X since we already know the solution
    if (importer != Teuchos::null) {
      tempX_->doExport(*tempExportX_, *importer, Tpetra::ADD);
      FullLHS->doExport(*tempExportX_, *importer, Tpetra::ADD);
    } else {
      tempX_->update(1.0, *tempExportX_, 0.0);   // tempX_ = 1.0 * tempExportX_ + 0.0 * tempX_
      FullLHS->update(1.0, *tempExportX_, 0.0);  // FullLHS = 1.0 * tempExportX_ + 0.0 * FullLHS
    }

    FullMatrix()->apply(*tempX_, *tempB_);
    tempB_->update(1.0, *FullRHS, -1.0);

    ReducedRHS_ = Teuchos::rcp(new multivector_type(ReducedMatrixRowMap(), FullRHS->getNumVectors()));
    ReducedRHS_->doImport(*tempB_, *Full2ReducedRHSImporter_, Tpetra::INSERT);

    // Finally construct Reduced Linear Problem
    ReducedProblem_ = Teuchos::rcp(new Tpetra::LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
        ReducedMatrix(), ReducedLHS_, ReducedRHS_));
  } else {
    // There are no singletons, so don't bother building a reduced problem.
    ReducedProblem_ = Problem;
    ReducedMatrix_  = Teuchos::rcp(dynamic_cast<crs_matrix_type*>(Problem->getMatrix().getRawPtr()), false);
  }

  double fn   = (double)FullMatrix()->getGlobalNumRows();
  double fnnz = (double)FullMatrix()->getGlobalNumEntries();

  double rn   = (double)ReducedMatrix()->getGlobalNumRows();
  double rnnz = (double)ReducedMatrix()->getGlobalNumEntries();

  RatioOfDimensions_  = rn / fn;
  RatioOfNonzeros_    = rnnz / fnnz;
  HaveReducedProblem_ = true;

  return;
}

//==============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    UpdateReducedProblem(const Teuchos::RCP<linear_problem_type>& Problem) {
  const char tfecfFuncName[] = "UpdateReducedProblem: ";

  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(!HaveReducedProblem_, std::runtime_error,
                                        "Must have a reduced problem.");
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(Problem == Teuchos::null,
                                        std::runtime_error, "Problem is Teuchos::null.");
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(FullMatrix() == Teuchos::null,
                                        std::runtime_error, "Need a RowMatrix.");
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(Problem->getRHS() == Teuchos::null,
                                        std::runtime_error, "Need a RHS.");
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(Problem->getLHS() == Teuchos::null,
                                        std::runtime_error, "Need a LHS.");

  if (SingletonsDetected()) {
    // Create pointer to Full RHS, LHS
    Teuchos::RCP<multivector_type> FullRHS = FullProblem()->getRHS();
    Teuchos::RCP<multivector_type> FullLHS = FullProblem()->getLHS();

    int NumVectors = FullLHS->getNumVectors();
    tempExportX_->putScalar(0.0);

    size_t NumEntries = 0;
    Teuchos::ArrayView<const Scalar> Values;
    Teuchos::Array<GlobalOrdinal> Indices;
    LocalOrdinal localNumRows        = FullMatrix()->getLocalNumRows();
    LocalOrdinal ColSingletonCounter = 0;

    for (LocalOrdinal i = 0; i < localNumRows; i++) {
      GlobalOrdinal curGRID = FullMatrixRowMap()->getGlobalElement(i);
      if (ReducedMatrixRowMap()->isNodeGlobalElement(curGRID)) {  // Check if this row should go into reduced matrix
        GetRowGCIDs(i, NumEntries, Values, Indices);

        // Filter indices and values
        Teuchos::Array<GlobalOrdinal> filteredIndices;
        Teuchos::Array<Scalar> filteredValues;

        for (typename Teuchos::Array<GlobalOrdinal>::size_type j = 0; j < Indices.size(); ++j) {
          if (ReducedMatrixColMap()->isNodeGlobalElement(Indices[j])) {
            filteredIndices.push_back(Indices[j]);
            filteredValues.push_back(Values[j]);
          }
        }

        // Insert filtered values into the matrix
        if (!filteredIndices.empty()) {
          ReducedMatrix()->replaceGlobalValues(curGRID, filteredIndices(), filteredValues());
        }

      }
      // Otherwise if singleton row we explicitly eliminate this row and solve for corresponding X value
      else {
        Teuchos::ArrayView<const LocalOrdinal> localIndices;
        GetRow(i, NumEntries, Values, localIndices);  // Get current row
        if (NumEntries == 1) {
          Scalar pivot = Values[0];
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(pivot == 0.0, std::runtime_error,
                                                "Encountered zero row, unable to continue.");  // Should improve this comparison to zero.
          LocalOrdinal indX = localIndices[0];
          for (LocalOrdinal j = 0; j < NumVectors; j++) {
            auto rhsData     = FullRHS->getData(j);               // Get the underlying data for vector j
            auto exportData  = tempExportX_->getDataNonConst(j);  // Get the underlying data for vector j (non-const)
            exportData[indX] = rhsData[i] / pivot;
          }
        } else {  // Singleton column
          LocalOrdinal j = ColSingletonPivotLIDs_[ColSingletonCounter];

          Scalar pivot = Values[j];
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(pivot == 0.0, std::runtime_error,
                                                "Encountered zero column, unable to continue.");  // Should improve this comparison to zero.
          ColSingletonPivots_[ColSingletonCounter] = pivot;
          ColSingletonCounter++;
        }
      }
    }

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(!(ColSingletonCounter == localNumSingletonCols_),
                                          std::runtime_error, "Sanity Check.");

    // Update Reduced LHS (Puts any initial guess values into reduced system)

    ReducedLHS_->putScalar(0.0);  // zero out Reduced LHS
    ReducedRHS_->doImport(*FullLHS, *Full2ReducedRHSImporter_, Tpetra::INSERT);

    FullLHS->putScalar(0.0);  // zero out Full LHS since we will inject values as we get them

    // Construct Reduced RHS

    // Zero out temp space
    tempX_->putScalar(0.0);
    tempB_->putScalar(0.0);

    // Inject known X values into tempX for purpose of computing tempB = FullMatrix*tempX
    // Also inject into full X since we already know the solution

    auto importer = FullCrsMatrix()->getCrsGraph()->getImporter();
    if (importer != Teuchos::null) {
      tempX_->doExport(*tempExportX_, *importer, Tpetra::ADD);
      FullLHS->doExport(*tempExportX_, *importer, Tpetra::ADD);
    } else {
      tempX_->update(1.0, *tempExportX_, 0.0);   // tempX_ = 1.0 * tempExportX_ + 0.0 * tempX_
      FullLHS->update(1.0, *tempExportX_, 0.0);  // FullLHS = 1.0 * tempExportX_ + 0.0 * FullLHS
    }

    FullMatrix()->apply(*tempX_, *tempB_);
    tempB_->update(1.0, *FullRHS, -1.0);

    ReducedRHS_->putScalar(0.0);
    ReducedRHS_->doImport(*tempB_, *Full2ReducedRHSImporter_, Tpetra::INSERT);

  } else {
    // There are no singletons, so don't bother building a reduced problem.
    ReducedProblem_ = Problem;
    ReducedMatrix_  = Teuchos::rcp(dynamic_cast<crs_matrix_type*>(Problem->getMatrix().getRawPtr()), false);
  }

  return;
}

////==============================================================================
// template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
// void CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
// ConstructRedistributeExporter(Epetra_Map * SourceMap, Epetra_Map * TargetMap,
//                               Epetra_Export * & RedistributeExporter,
//                               Epetra_Map * & RedistributeMap) {
//
//   int_type IndexBase = (int_type) SourceMap->IndexBase64();
//   if (IndexBase!=(int_type) TargetMap->IndexBase64()) EPETRA_CHK_ERR(-1);
//
//   const Epetra_Comm & Comm = TargetMap->Comm();
//
//   int TargetNumMyElements = TargetMap->NumMyElements();
//   int SourceNumMyElements = SourceMap->NumMyElements();
//
//   // ContiguousTargetMap has same number of elements per PE as TargetMap, but uses contigious indexing
//   Epetra_Map ContiguousTargetMap((int_type) -1, TargetNumMyElements, IndexBase,Comm);
//
//   // Same for ContiguousSourceMap
//   Epetra_Map ContiguousSourceMap((int_type) -1, SourceNumMyElements, IndexBase, Comm);
//
//   assert(ContiguousSourceMap.NumGlobalElements64()==ContiguousTargetMap.NumGlobalElements64());
//
//   // Now create a vector that contains the global indices of the Source Epetra_MultiVector
//   int_type* SourceMapMyGlobalElements = 0;
//   SourceMap->MyGlobalElementsPtr(SourceMapMyGlobalElements);
//   typename Epetra_GIDTypeVector<int_type>::impl SourceIndices(View, ContiguousSourceMap, SourceMapMyGlobalElements);
//
//   // Create an exporter to send the SourceMap global IDs to the target distribution
//   Epetra_Export Exporter(ContiguousSourceMap, ContiguousTargetMap);
//
//   // Create a vector to catch the global IDs in the target distribution
//   typename Epetra_GIDTypeVector<int_type>::impl TargetIndices(ContiguousTargetMap);
//   TargetIndices.Export(SourceIndices, Exporter, Insert);
//
//   // Create a new map that describes how the Source MultiVector should be laid out so that it has
//   // the same number of elements on each processor as the TargetMap
//   RedistributeMap = new Epetra_Map((int_type) -1, TargetNumMyElements, TargetIndices.Values(), IndexBase, Comm);
//
//   // This exporter will finally redistribute the Source MultiVector to the same layout as the TargetMap
//   RedistributeExporter = new Epetra_Export(*SourceMap, *RedistributeMap);
//   return(0);
// }
//
////==============================================================================
// int CrsSingletonFilter_LinearProblem::ComputeFullSolution() {
//
//   if ( SingletonsDetected() ) {
//     int jj, k;
//
//     Epetra_MultiVector * FullLHS = FullProblem()->GetLHS();
//     Epetra_MultiVector * FullRHS = FullProblem()->GetRHS();
//
//     tempX_->PutScalar(0.0); tempExportX_->PutScalar(0.0);
//     // Inject values that the user computed for the reduced problem into the full solution vector
//     EPETRA_CHK_ERR(tempX_->Export(*ReducedLHS_, *Full2ReducedLHSImporter_, Add));
//
//     FullLHS->Update(1.0, *tempX_, 1.0);
//
//     // Next we will use our full solution vector which is populated with pre-filter solution
//     // values and reduced system solution values to compute the sum of the row contributions
//     // that must be subtracted to get the post-filter solution values
//
//     EPETRA_CHK_ERR(FullMatrix()->Multiply(false, *FullLHS, *tempB_));
//
//     // Finally we loop through the local rows that were associated with column singletons and compute the
//     // solution for these equations.
//
//     int NumVectors = tempB_->NumVectors();
//     for (k=0; k<localNumSingletonCols_; k++) {
//       int i = ColSingletonRowLIDs_[k];
//       int j = ColSingletonColLIDs_[k];
//       double pivot = ColSingletonPivots_[k];
//       for (jj=0; jj<NumVectors; jj++)
//         (*tempExportX_)[jj][j]= ((*FullRHS)[jj][i] - (*tempB_)[jj][i])/pivot;
//     }
//
//     // Finally, insert values from post-solve step and we are done!!!!
//
//     if (FullMatrix()->RowMatrixImporter()!=0) {
//       EPETRA_CHK_ERR(tempX_->Export(*tempExportX_, *FullMatrix()->RowMatrixImporter(), Add));
//     }
//     else {
//       tempX_->Update(1.0, *tempExportX_, 0.0);
//     }
//
//     FullLHS->Update(1.0, *tempX_, 1.0);
//   }
//
//   return(0);
// }
//==============================================================================

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    InitFullMatrixAccess() {
  localMaxNumRowEntries_ = FullMatrix()->getLocalMaxNumRowEntries();

  // Cast to CrsMatrix, if possible.  Can save some work.
  FullCrsMatrix_ = Teuchos::rcp_dynamic_cast<crs_matrix_type>(FullMatrix());
  FullMatrixIsCrsMatrix_ = (FullCrsMatrix() != Teuchos::null);  // Pointer is non-null if cast worked
  Indices_               = Teuchos::arcp(new global_ordinal_type[localMaxNumRowEntries_], 0, localMaxNumRowEntries_, true);
  for (auto i = 0; i < Indices_.size(); ++i) {
    Indices_[i] = 0;
  }
  Values_ = Teuchos::arcp(new scalar_type[localMaxNumRowEntries_], 0, localMaxNumRowEntries_, true);
  for (auto i = 0; i < Values_.size(); ++i) {
    Values_[i] = 0;
  }

  return;
}

//==============================================================================

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    GetRow(local_ordinal_type localRow, size_t& NumIndices,
           Teuchos::Array<local_ordinal_type>& localIndices) {
  // Must get the values, but we ignore them.
  size_t maxNumEntries = FullCrsMatrix()->getLocalMaxNumRowEntries();
  nonconst_local_inds_host_view_type Indices("indices", maxNumEntries);
  nonconst_values_host_view_type Values("values", maxNumEntries);
  if (FullMatrixIsCrsMatrix_) {  // View of current row
    FullCrsMatrix()->getLocalRowCopy(localRow, Indices, Values, NumIndices);
  } else {  // Copy of current row (we must get the values, but we ignore them)
    FullMatrix()->getLocalRowCopy(localRow, Indices, Values, NumIndices);
  }
  // If LocalRow does not belong to the calling process, then the method sets
  // NumIndices to Teuchos::OrdinalTraits<size_t>::invalid(), and does not
  // modify localIndices or Values.

  // Copy the indices to the output array
  localIndices.resize(NumIndices);
  for (size_t i = 0; i < NumIndices; ++i) {
    localIndices[i] = Indices(i);
  }
  return;
}

//==============================================================================

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    GetRow(LocalOrdinal Row, size_t& NumIndices, Teuchos::ArrayView<const Scalar>& Values,
           Teuchos::ArrayView<const LocalOrdinal>& Indices) {
  typename Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_inds_host_view_type localIndices;
  typename Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::values_host_view_type rowValues;

  if (FullMatrixIsCrsMatrix_) {  // View of current row
    FullCrsMatrix()->getLocalRowView(Row, localIndices, rowValues);

    NumIndices = localIndices.size();
    Values     = Teuchos::ArrayView<const Scalar>(rowValues.data(), rowValues.size());
    Indices    = Teuchos::ArrayView<const LocalOrdinal>(localIndices.data(), localIndices.size());

  } else {  // Copy of current row
    typename Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::nonconst_local_inds_host_view_type localIndicesCopy("localIndicesCopy", FullMatrix()->getLocalMaxNumRowEntries());
    typename Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::nonconst_values_host_view_type rowValuesCopy("rowValuesCopy", FullMatrix()->getLocalMaxNumRowEntries());

    FullMatrix()->getLocalRowCopy(Row, localIndicesCopy, rowValuesCopy, NumIndices);

    Values  = Teuchos::ArrayView<const Scalar>(rowValuesCopy.data(), NumIndices);
    Indices = Teuchos::ArrayView<const LocalOrdinal>(localIndicesCopy.data(), NumIndices);
  }
  return;
}

//==============================================================================

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    GetRowGCIDs(
        LocalOrdinal localRow,
        size_t& NumIndices,
        Teuchos::ArrayView<const Scalar>& Values,
        Teuchos::Array<GlobalOrdinal>& GlobalIndices) {
  // Extract the row data (local indices and values)
  typename Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_inds_host_view_type LocalIndices;
  typename Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::values_host_view_type RowValues;

  FullMatrix()->getLocalRowView(localRow, LocalIndices, RowValues);

  // Convert local indices to global indices
  NumIndices = LocalIndices.size();
  GlobalIndices.resize(NumIndices);  // Resize the array to hold global indices
  for (size_t j = 0; j < NumIndices; ++j) {
    GlobalIndices[j] = FullMatrixColMap()->getGlobalElement(LocalIndices[j]);
  }
  // Copy values into the provided ArrayView
  Values = Teuchos::ArrayView<const Scalar>(RowValues.data(), RowValues.size());

  return;
}

//==============================================================================

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    CreatePostSolveArrays(vector_type_LO localRowIDofSingletonCol,
                          vector_type_LO ColProfiles,
                          vector_type_LO NewColProfiles,
                          vector_type_LO ColHasRowWithSingleton) {
  const char tfecfFuncName[] = "CreatePostSolveArrays: ";

  if (localNumSingletonCols_ == 0) return;  // Nothing to do

  // We will need these arrays for the post-solve phase
  ColSingletonRowLIDs_   = Teuchos::arcp(new local_ordinal_type[localNumSingletonCols_], 0, localNumSingletonCols_, true);
  ColSingletonColLIDs_   = Teuchos::arcp(new local_ordinal_type[localNumSingletonCols_], 0, localNumSingletonCols_, true);
  ColSingletonPivotLIDs_ = Teuchos::arcp(new local_ordinal_type[localNumSingletonCols_], 0, localNumSingletonCols_, true);
  ColSingletonPivots_    = Teuchos::arcp(new scalar_type[localNumSingletonCols_], 0, localNumSingletonCols_, true);

  // Register singleton columns (that were not already counted as singleton rows)
  // Check to see if any columns disappeared because all associated rows were eliminated
  int NumMyColSingletonstmp         = 0;
  local_ordinal_type localNumCols   = FullMatrix()->getLocalNumCols();
  auto localRowIDofSingletonColData = localRowIDofSingletonCol.getLocalViewHost(Tpetra::Access::ReadWrite);
  auto ColProfilesData              = ColProfiles.getLocalViewHost(Tpetra::Access::ReadWrite);
  auto RowMapColors_Data            = RowMapColors_->getLocalViewHost(Tpetra::Access::ReadWrite);
  auto NewColProfilesData           = NewColProfiles.getLocalViewHost(Tpetra::Access::ReadWrite);
  auto ColHasRowWithSingletonData   = ColHasRowWithSingleton.getLocalViewHost(Tpetra::Access::ReadWrite);
  auto ColMapColors_Data            = ColMapColors_->getLocalViewHost(Tpetra::Access::ReadWrite);

  for (int j = 0; j < localNumCols; j++) {
    int i = localRowIDofSingletonColData(j, 0);
    if (ColProfilesData(j, 0) == 1 && RowMapColors_Data(i, 0) != 1) {
      ColSingletonRowLIDs_[NumMyColSingletonstmp] = i;
      ColSingletonColLIDs_[NumMyColSingletonstmp] = j;
      NumMyColSingletonstmp++;
    }
    // Also check for columns that were eliminated implicitly by
    // having all associated row eliminated
    else if (NewColProfilesData(j, 0) == 0 && ColHasRowWithSingletonData(j, 0) != 1 && RowMapColors_Data(i, 0) == 0) {
      ColMapColors_Data(j, 0) = 1;
    }
  }

  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(NumMyColSingletonstmp != localNumSingletonCols_,
                                        std::runtime_error, "Sanity check.");

  Tpetra::sort2(ColSingletonRowLIDs_.begin(), ColSingletonRowLIDs_.end(), ColSingletonColLIDs_.begin());

  return;
}

}  // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_CRSSINGLETONFILTER_INSTANT(SCALAR, LO, GO, NODE) \
  template class CrsSingletonFilter_LinearProblem<SCALAR, LO, GO, NODE>;

#endif  //  TPETRA_CRSSINGLETONFILTER_LINEARPROBLEM_DEF_HPP
