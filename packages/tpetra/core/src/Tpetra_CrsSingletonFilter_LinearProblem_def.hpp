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

// Functor to create post-solve arrays
template <class LocalLO_type, class ColIDView_type, class ColMapColor_type, class RowMapColor_type, class Counter_type>
struct CreatePostSolveArraysFunctor {
  LocalLO_type LocalRowIDofSingletonColData;
  ColIDView_type ColSingletonRowLIDs;
  ColIDView_type ColSingletonColLIDs;

  LocalLO_type NewColProfilesData;
  LocalLO_type ColProfilesData;
  LocalLO_type ColHasRowWithSingletonData;
  ColMapColor_type ColMapColors_Data;
  RowMapColor_type RowMapColors_Data;
  Counter_type numSingletonCols;

  CreatePostSolveArraysFunctor(LocalLO_type localRowIDofSingletonColData,
                               ColIDView_type colSingletonRowLIDs, ColIDView_type colSingletonColLIDs,
                               LocalLO_type newColProfilesData, LocalLO_type colProfilesData, LocalLO_type colHasRowWithSingletonData,
                               ColMapColor_type colMapColors_Data, RowMapColor_type rowMapColors_Data, Counter_type in_numSingletonCols)
    : LocalRowIDofSingletonColData(localRowIDofSingletonColData)
    , ColSingletonRowLIDs(colSingletonRowLIDs)
    , ColSingletonColLIDs(colSingletonColLIDs)
    , NewColProfilesData(newColProfilesData)
    , ColProfilesData(colProfilesData)
    , ColHasRowWithSingletonData(colHasRowWithSingletonData)
    , ColMapColors_Data(colMapColors_Data)
    , RowMapColors_Data(rowMapColors_Data)
    , numSingletonCols(in_numSingletonCols) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t j) const {
    int i = LocalRowIDofSingletonColData(j, 0);
    if (ColProfilesData(j, 0) == 1 && RowMapColors_Data(i, 0) != 1) {
      // These will be sorted by rowLIDs, so no need to be in any particular order
      auto lclNumSingletonCols                 = Kokkos::atomic_fetch_add(&numSingletonCols(0), 1);
      ColSingletonRowLIDs(lclNumSingletonCols) = i;
      ColSingletonColLIDs(lclNumSingletonCols) = j;
    }
    // Also check for columns that were eliminated implicitly by
    // having all associated row eliminated
    else if (NewColProfilesData(j, 0) == 0 && ColHasRowWithSingletonData(j, 0) != 1 && RowMapColors_Data(i, 0) == 0) {
      ColMapColors_Data(j, 0) = 1;
    }
  }
};

// Functors to construct or update reduced problem
template <class LocalOrdinal, class GlobalOrdinal, class LocalMap_type, class LocalMatrix_type, class row_map_type, class entries_type, class values_type>
struct ConstructReducedProblemFunctor {
  LocalMap_type lclFullRowMap;
  LocalMap_type lclFullColMap;
  LocalMap_type lclReducedRowMap;
  LocalMap_type lclReducedColMap;

  LocalMatrix_type lclFullMatrix;

  row_map_type lclReducedRowPtr;
  entries_type lclReducedColInd;
  values_type lclReducedValues;

  // Constructor
  ConstructReducedProblemFunctor(LocalMap_type lclFullRowMap_, LocalMap_type lclFullColMap_, LocalMap_type lclReducedRowMap_, LocalMap_type lclReducedColMap_,
                                 LocalMatrix_type lclFullMatrix_, row_map_type lclReducedRowPtr_, entries_type lclReducedColInd_, values_type lclReducedValues_)
    : lclFullRowMap(lclFullRowMap_)
    , lclFullColMap(lclFullColMap_)
    , lclReducedRowMap(lclReducedRowMap_)
    , lclReducedColMap(lclReducedColMap_)
    , lclFullMatrix(lclFullMatrix_)
    , lclReducedRowPtr(lclReducedRowPtr_)
    , lclReducedColInd(lclReducedColInd_)
    , lclReducedValues(lclReducedValues_) {}

  // Functor to count nonzero entries
  KOKKOS_INLINE_FUNCTION
  void operator()(const LocalOrdinal i, size_t& nnz) const {
    const LocalOrdinal INVALID = Tpetra::Details::OrdinalTraits<LocalOrdinal>::invalid();

    auto lclFullRowPtr = lclFullMatrix.graph.row_map;
    auto lclFullColInd = lclFullMatrix.graph.entries;
    auto lclFullValues = lclFullMatrix.values;

    GlobalOrdinal glbRowID = lclFullRowMap.getGlobalElement(i);
    LocalOrdinal lclRowID  = lclReducedRowMap.getLocalElement(glbRowID);

    if (lclRowID != INVALID) {
      for (size_t j = lclFullRowPtr(i); j < lclFullRowPtr(i + 1); j++) {
        GlobalOrdinal glbColID = lclFullColMap.getGlobalElement(lclFullColInd[j]);
        LocalOrdinal lclColID  = lclReducedColMap.getLocalElement(glbColID);
        if (lclColID != INVALID) {
          lclReducedRowPtr(1 + lclRowID)++;
        }
      }
    }
    nnz += lclReducedRowPtr(1 + lclRowID);
    ;
  }

  // Functor to insert nonzero entries
  KOKKOS_INLINE_FUNCTION
  void operator()(const LocalOrdinal i) const {
    const LocalOrdinal INVALID = Tpetra::Details::OrdinalTraits<LocalOrdinal>::invalid();

    auto lclFullRowPtr = lclFullMatrix.graph.row_map;
    auto lclFullColInd = lclFullMatrix.graph.entries;
    auto lclFullValues = lclFullMatrix.values;

    GlobalOrdinal glbRowID = lclFullRowMap.getGlobalElement(i);
    LocalOrdinal lclRowID  = lclReducedRowMap.getLocalElement(glbRowID);

    if (lclRowID != INVALID) {
      size_t k = lclReducedRowPtr(lclRowID);
      for (size_t j = lclFullRowPtr(i); j < lclFullRowPtr(i + 1); j++) {
        GlobalOrdinal glbColID = lclFullColMap.getGlobalElement(lclFullColInd[j]);
        LocalOrdinal lclColID  = lclReducedColMap.getLocalElement(glbColID);
        if (lclColID != INVALID) {
          lclReducedColInd[k] = lclColID;
          lclReducedValues[k] = lclFullValues[j];
          k++;
        }
      }
    }
  }
};

template <class LocalOrdinal, class GlobalOrdinal, class LocalMap_type, class LocalMatrix_type>
struct UpdateReducedProblemFunctor {
  LocalMap_type lclFullRowMap;
  LocalMap_type lclFullColMap;
  LocalMap_type lclReducedRowMap;
  LocalMap_type lclReducedColMap;

  LocalMatrix_type lclFullMatrix;
  LocalMatrix_type lclReducedMatrix;

  // Constructor
  UpdateReducedProblemFunctor(LocalMap_type lclFullRowMap_, LocalMap_type lclFullColMap_,
                              LocalMap_type lclReducedRowMap_, LocalMap_type lclReducedColMap_,
                              LocalMatrix_type lclFullMatrix_, LocalMatrix_type lclReducedMatrix_)
    : lclFullRowMap(lclFullRowMap_)
    , lclFullColMap(lclFullColMap_)
    , lclReducedRowMap(lclReducedRowMap_)
    , lclReducedColMap(lclReducedColMap_)
    , lclFullMatrix(lclFullMatrix_)
    , lclReducedMatrix(lclReducedMatrix_) {}

  // Functor to update nonzero entries
  KOKKOS_INLINE_FUNCTION
  void operator()(const LocalOrdinal i) const {
    const LocalOrdinal INVALID = Tpetra::Details::OrdinalTraits<LocalOrdinal>::invalid();

    auto lclFullRowPtr = lclFullMatrix.graph.row_map;
    auto lclFullColInd = lclFullMatrix.graph.entries;
    auto lclFullValues = lclFullMatrix.values;

    auto lclReducedRowPtr = lclReducedMatrix.graph.row_map;
    auto lclReducedColInd = lclReducedMatrix.graph.entries;
    auto lclReducedValues = lclReducedMatrix.values;

    GlobalOrdinal glbRowID = lclFullRowMap.getGlobalElement(i);
    LocalOrdinal lclRowID  = lclReducedRowMap.getLocalElement(glbRowID);

    if (lclRowID != INVALID) {
      for (size_t j = lclFullRowPtr(i); j < lclFullRowPtr(i + 1); j++) {
        GlobalOrdinal glbColID = lclFullColMap.getGlobalElement(lclFullColInd(j));
        LocalOrdinal lclColID  = lclReducedColMap.getLocalElement(glbColID);
        if (lclColID != INVALID) {
          for (size_t k = lclReducedRowPtr(lclRowID); k < lclReducedRowPtr(lclRowID + 1); k++) {
            if (lclReducedColInd[k] == lclColID) {
              lclReducedValues[k] = lclFullValues[j];
            }
          }
        }
      }
    }
  }
};

// Functor to solve singleton problem
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node,
          class LocalMap_type, class LocalMatrix_type, class LocalX_type, class LocalB_type,
          class View_type_int, class View_type_scalar, class Err_type>
struct SolveSingletonProblemFunctor {
  using local_ordinal_type = LocalOrdinal;
  using impl_scalar_type   = typename Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::impl_scalar_type;

  local_ordinal_type NumVectors;
  local_ordinal_type localNumSingletonCols;

  Err_type error_code;
  LocalMap_type lclFullRowMap;
  LocalMap_type lclFullColMap;
  LocalMap_type lclReducedRowMap;

  LocalMatrix_type lclFullMatrix;

  LocalX_type localExportX;
  LocalB_type localRHS;

  View_type_int ColSingletonColLIDs;
  View_type_int ColSingletonRowLIDs;
  View_type_int ColSingletonPivotLIDs;
  View_type_scalar ColSingletonPivots;

  // Constructor
  SolveSingletonProblemFunctor(local_ordinal_type NumVectors_, local_ordinal_type localNumSingletonCols_, Err_type error_code_,
                               LocalMap_type lclFullRowMap_, LocalMap_type lclFullColMap_, LocalMap_type lclReducedRowMap_,
                               LocalMatrix_type lclFullMatrix_, LocalX_type localExportX_, LocalB_type localRHS_,
                               View_type_int colSingletonColLIDs_, View_type_int colSingletonRowLIDs_,
                               View_type_int colSingletonPivotLIDs_, View_type_scalar colSingletonPivots_)
    : NumVectors(NumVectors_)
    , localNumSingletonCols(localNumSingletonCols_)
    , error_code(error_code_)
    , lclFullRowMap(lclFullRowMap_)
    , lclFullColMap(lclFullColMap_)
    , lclReducedRowMap(lclReducedRowMap_)
    , lclFullMatrix(lclFullMatrix_)
    , localExportX(localExportX_)
    , localRHS(localRHS_)
    , ColSingletonColLIDs(colSingletonColLIDs_)
    , ColSingletonRowLIDs(colSingletonRowLIDs_)
    , ColSingletonPivotLIDs(colSingletonPivotLIDs_)
    , ColSingletonPivots(colSingletonPivots_) {}

  // Functor used in ConstructReducedProblem with parallel-for
  KOKKOS_INLINE_FUNCTION
  void operator()(const local_ordinal_type i) const {
    const LocalOrdinal INVALID = Tpetra::Details::OrdinalTraits<LocalOrdinal>::invalid();
    const impl_scalar_type zero(0);

    auto lclFullRowPtr = lclFullMatrix.graph.row_map;
    auto lclFullColInd = lclFullMatrix.graph.entries;
    auto lclFullValues = lclFullMatrix.values;

    GlobalOrdinal glbID = lclFullRowMap.getGlobalElement(i);
    LocalOrdinal lclID  = lclReducedRowMap.getLocalElement(glbID);

    if (lclID == INVALID) {
      if (lclFullRowPtr(i + 1) == lclFullRowPtr(i) + 1) {
        // Singleton row
        LocalOrdinal indX      = lclFullColInd(lclFullRowPtr(i));
        impl_scalar_type pivot = lclFullValues(lclFullRowPtr(i));
        if (pivot == zero) {
          Kokkos::atomic_exchange(&error_code(0), 1);
        } else {
          for (LocalOrdinal j = 0; j < NumVectors; j++) {
            localExportX(indX, j) = localRHS(i, j) / pivot;
          }
        }
      } else {
        // Not singleton row, but need to be removed == Singleton column
        //  look for matching row ID
        int myNum = -1;
        for (local_ordinal_type j = 0; j < localNumSingletonCols; j++) {
          if (ColSingletonRowLIDs[j] == i) {
            myNum = j;
          }
        }
        if (myNum == -1) {
          Kokkos::atomic_exchange(&error_code(0), 3);
        } else {
          //  look for matching col ID
          LocalOrdinal targetCol = ColSingletonColLIDs[myNum];
          for (size_t j = lclFullRowPtr(i); j < lclFullRowPtr(i + 1); j++) {
            if (lclFullColMap.getGlobalElement(lclFullColInd[j]) == targetCol) {
              impl_scalar_type pivot = lclFullValues[j];
              if (pivot == zero) {
                Kokkos::atomic_exchange(&error_code(0), 2);
              } else {
                ColSingletonPivotLIDs[myNum] = j;
                ColSingletonPivots[myNum]    = pivot;
              }
              break;
            }
          }
        }
      }
    }
  }

  // Functor used in UpdateReducedProblem with parallel-reduce
  KOKKOS_INLINE_FUNCTION
  void operator()(const local_ordinal_type i, local_ordinal_type& lclNumSingletonCols) const {
    const LocalOrdinal INVALID = Tpetra::Details::OrdinalTraits<LocalOrdinal>::invalid();
    const impl_scalar_type zero(0);

    auto lclFullRowPtr = lclFullMatrix.graph.row_map;
    auto lclFullColInd = lclFullMatrix.graph.entries;
    auto lclFullValues = lclFullMatrix.values;

    GlobalOrdinal glbID = lclFullRowMap.getGlobalElement(i);
    LocalOrdinal lclID  = lclReducedRowMap.getLocalElement(glbID);

    if (lclID == INVALID) {
      if (lclFullRowPtr(i + 1) == lclFullRowPtr(i) + 1) {
        // Singleton row
        LocalOrdinal indX      = lclFullColInd(lclFullRowPtr(i));
        impl_scalar_type pivot = lclFullValues(lclFullRowPtr(i));
        if (pivot == zero) {
          Kokkos::atomic_exchange(&error_code(0), 1);
        } else {
          for (LocalOrdinal j = 0; j < NumVectors; j++) {
            localExportX(indX, j) = localRHS(i, j) / pivot;
          }
        }
      } else {
        // Not singleton row, but need to be removed == Singleton column
        //  look for matching row ID
        int myNum = -1;
        for (local_ordinal_type j = 0; j < localNumSingletonCols; j++) {
          if (ColSingletonRowLIDs[j] == i) {
            myNum = j;
          }
        }
        if (myNum == -1) {
          Kokkos::atomic_exchange(&error_code(0), 3);
        } else {
          //  look for matching col ID
          LocalOrdinal targetCol = ColSingletonColLIDs[myNum];
          for (size_t j = lclFullRowPtr(i); j < lclFullRowPtr(i + 1); j++) {
            if (lclFullColMap.getGlobalElement(lclFullColInd[j]) == targetCol) {
              impl_scalar_type pivot = lclFullValues[j];
              if (pivot == zero) {
                Kokkos::atomic_exchange(&error_code(0), 2);
              } else {
                ColSingletonPivots[myNum] = pivot;
              }
              break;
            }
          }
          lclNumSingletonCols++;
        }
      }
    }
  }
};

// Functor to compute full singleton problem
template <class View_type_int, class View_type_scalar, class SOL_type, class RHS_type>
struct ComputeFullSolutionFunctor {
  using scalar_type = typename View_type_scalar::non_const_value_type;

  View_type_int ColSingletonRowLIDs;
  View_type_int ColSingletonColLIDs;
  View_type_scalar ColSingletonPivots;

  SOL_type localX;
  RHS_type localRHS;
  RHS_type localB;

  ComputeFullSolutionFunctor(View_type_int rowLIDs, View_type_int colLIDs, View_type_scalar pivots,
                             SOL_type X, RHS_type RHS, RHS_type B)
    : ColSingletonRowLIDs(rowLIDs)
    , ColSingletonColLIDs(colLIDs)
    , ColSingletonPivots(pivots)
    , localX(X)
    , localRHS(RHS)
    , localB(B) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t k) const {
    int i             = ColSingletonRowLIDs[k];
    int j             = ColSingletonColLIDs[k];
    scalar_type pivot = ColSingletonPivots[k];

    size_t NumVectors = localX.extent(1);
    for (size_t jj = 0; jj < NumVectors; jj++) {
      localX(j, jj) = (localRHS(i, jj) - localB(i, jj)) / pivot;
    }
  }
};

template <class Map_color_type, class LocalMap_type, class IndexList_type>
struct GenerateReducedMapFunctor {
  int color;
  Map_color_type mapColors;
  LocalMap_type lclMap;
  IndexList_type indexList;

  // .. internal view used to count match ..
  using device_type = typename Map_color_type::device_type;
  Kokkos::View<size_t*, device_type> numView;

  // constructor for counting
  GenerateReducedMapFunctor(int color_, Map_color_type mapColors_)
    : color(color_)
    , mapColors(mapColors_) {}

  // constructor for inserting
  GenerateReducedMapFunctor(int color_, Map_color_type mapColors_, LocalMap_type lclMap_, IndexList_type indexList_)
    : color(color_)
    , mapColors(mapColors_)
    , lclMap(lclMap_)
    , indexList(indexList_) {
    Kokkos::resize(numView, 1);
    Kokkos::deep_copy(numView, 0);
  }

  // Functor to count number of elements with matching color (parallel-reduce)
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t i, size_t& num) const {
    if (mapColors(i, 0) == color) {
      num++;
    }
  }

  // Functor to insert elements with matching color (parallel-for)
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t i) const {
    if (mapColors(i, 0) == color) {
      auto num       = Kokkos::atomic_fetch_add(&numView(0), 1);  // ..atomic to count..
      indexList(num) = lclMap.getGlobalElement(i);                // to global
    }
  }
};

//==============================================================================

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    CrsSingletonFilter_LinearProblem(bool run_on_host, bool verbose)
  : globalNumSingletonRows_(Teuchos::as<local_ordinal_type>(0))
  , globalNumSingletonCols_(Teuchos::as<local_ordinal_type>(0))
  , RatioOfDimensions_(0.0)
  , RatioOfNonzeros_(0.0)
  , HaveReducedProblem_(false)
  , AnalysisDone_(false)
  , SymmetricElimination_(true)
  , localMaxNumRowEntries_(Teuchos::as<local_ordinal_type>(0))
  , FullMatrixIsCrsMatrix_(false)
  , run_on_host_(run_on_host)
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
  ComputeFullSolution();
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
  localRowIDofSingletonCol.putScalar(Teuchos::OrdinalTraits<LocalOrdinal>::invalid());

  // Define MapColoring objects
  RowMapColors_ = Teuchos::rcp(new vector_type_int(FullMatrixRowMap()));  // Initial colors are all 0
  ColMapColors_ = Teuchos::rcp(new vector_type_int(FullMatrixColMap()));

  local_ordinal_type localNumRows = FullMatrix()->getLocalNumRows();
  local_ordinal_type localNumCols = FullMatrix()->getLocalNumCols();

  // Set up for accessing full matrix.  Will do so row-by-row.
  InitFullMatrixAccess();

  // Scan matrix for singleton rows, build up column profiles
  localNumSingletonRows_ = 0;
  if (run_on_host_) {
    size_t NumIndices = 1;
    // int * localIndices;
    // nonconst_local_inds_host_view_type localIndices;
    Teuchos::Array<local_ordinal_type> localIndices;
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
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(!FullMatrixIsCrsMatrix_, std::runtime_error,
                                          "Error: FullMatrix is not CrsMatrix");

    auto ColProfilesData              = ColProfiles.getLocalViewDevice(Tpetra::Access::ReadWrite);
    auto localRowIDofSingletonColData = localRowIDofSingletonCol.getLocalViewDevice(Tpetra::Access::ReadWrite);
    auto ColHasRowWithSingletonData   = ColHasRowWithSingleton.getLocalViewDevice(Tpetra::Access::ReadWrite);
    auto RowMapColors_Data            = RowMapColors_->getLocalViewDevice(Tpetra::Access::ReadWrite);
    auto ColMapColors_Data            = ColMapColors_->getLocalViewDevice(Tpetra::Access::ReadWrite);

    Kokkos::View<int*, execution_space, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::Atomic>>
        ColProfilesAtomic(ColProfilesData.data(), ColProfilesData.extent(0));
    Kokkos::View<int*, execution_space, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::Atomic>>
        ColHasRowWithSingletonAtomic(ColHasRowWithSingletonData.data(), ColHasRowWithSingletonData.extent(0));
    Kokkos::View<int*, execution_space, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::Atomic>>
        ColMapColors_Atomic(ColMapColors_Data.data(), ColMapColors_Data.extent(0));
    Kokkos::View<int*, execution_space> error_code("singleton-bound-check", 2);
    Kokkos::deep_copy(error_code, 0);

    auto lclRowPtr = FullCrsMatrix_->getLocalRowPtrsDevice();
    auto lclColInd = FullCrsMatrix_->getLocalIndicesDevice();

    Kokkos::parallel_reduce(
        "CrsSingletonFilter_LinearProblem:Analyze(find-singleton-row)", range_policy(0, localNumRows),
        KOKKOS_LAMBDA(const size_t i, local_ordinal_type& lclNumSingletonRows) {
          // Get ith row
          for (size_t k = lclRowPtr(i); k < lclRowPtr(i + 1); k++) {
            local_ordinal_type ColumnIndex = lclColInd(k);

            // Bounds check for ColumnIndex
            if (static_cast<size_t>(ColumnIndex) >= ColProfilesData.extent(0)) {
              Kokkos::atomic_exchange(&error_code(0), ColumnIndex + 1);
            } else {
              ColProfilesAtomic(ColumnIndex)++;  // Increment column count
            }

            // Bounds check for ColumnIndex
            if (static_cast<size_t>(ColumnIndex) >= localRowIDofSingletonColData.extent(0)) {
              Kokkos::atomic_exchange(&error_code(1), ColumnIndex + 1);
            } else {
              // Record local row ID for current column
              // will use to identify row to eliminate if column is a singleton
              // It will store the last local row index where this column has non-zero entry
              //  (it may not be singleton row?)
              Kokkos::atomic_max(&localRowIDofSingletonColData(ColumnIndex, 0), i);
            }
          }
          // If row has single entry, color it and associated column with color=1
          if (lclRowPtr(i + 1) == lclRowPtr(i) + 1) {
            local_ordinal_type ColumnIndex = lclColInd(lclRowPtr(i));
            RowMapColors_Data(i, 0)        = 1;

            // this column will be eliminated by this singleton row
            ColHasRowWithSingletonAtomic(ColumnIndex)++;
            ColMapColors_Atomic(ColumnIndex) = 1;
            lclNumSingletonRows++;
          }
        },
        localNumSingletonRows_);

    // Check error code for bound check
    // NOTE: Should I broadcast the error-code?
    auto h_error = Kokkos::create_mirror_view(error_code);
    Kokkos::deep_copy(h_error, error_code);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(h_error(0) > 0, std::runtime_error,
                                          "Error: ColumnIndex out of bounds: " + std::to_string(h_error(0)));
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(h_error(1) > 0, std::runtime_error,
                                          "Error: ColumnIndex out of bounds for localRowIDofSingletonColData " + std::to_string(h_error(1)));
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

  // ColProfiles now contains the nonzero column entry count for all columns that have
  // an entry on this processor.
  // ColHasRowWithSingleton now contains a count of singleton rows associated with the corresponding
  // local column.

  // Next we check to make sure no column is associated with more than one singleton row.
  // If multiple singleton rows have nonzero entries at the same column, the matrix is singular?
  // [0 0 0 0 0 x 0 0]
  // [0 0 0 0 0 x 0 0]
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(ColHasRowWithSingleton.normInf() > 1, std::runtime_error,
                                        "At least one column is associated with two singleton rows, can't handle it.");

  localNumSingletonCols_ = 0;
  vector_type_int RowHasColWithSingleton(FullMatrix()->getRowMap());  // Use to check for errors
  RowHasColWithSingleton.putScalar(0);
  if (run_on_host_) {
    size_t NumIndices = 1;
    // int * localIndices;
    // nonconst_local_inds_host_view_type localIndices;
    Teuchos::Array<local_ordinal_type> localIndices;

    auto ColProfilesData              = ColProfiles.getLocalViewHost(Tpetra::Access::ReadOnly);
    auto localRowIDofSingletonColData = localRowIDofSingletonCol.getLocalViewHost(Tpetra::Access::ReadOnly);
    auto ColHasRowWithSingletonData   = ColHasRowWithSingleton.getLocalViewHost(Tpetra::Access::ReadOnly);

    auto RowMapColors_Data = RowMapColors_->getLocalViewHost(Tpetra::Access::ReadWrite);
    auto ColMapColors_Data = ColMapColors_->getLocalViewHost(Tpetra::Access::ReadWrite);

    auto NewColProfilesData         = NewColProfiles.getLocalViewHost(Tpetra::Access::ReadWrite);
    auto RowHasColWithSingletonData = RowHasColWithSingleton.getLocalViewHost(Tpetra::Access::ReadWrite);

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
  } else {
    auto ColProfilesData              = ColProfiles.getLocalViewDevice(Tpetra::Access::ReadOnly);
    auto localRowIDofSingletonColData = localRowIDofSingletonCol.getLocalViewDevice(Tpetra::Access::ReadOnly);
    auto ColHasRowWithSingletonData   = ColHasRowWithSingleton.getLocalViewDevice(Tpetra::Access::ReadOnly);

    auto RowMapColors_Data = RowMapColors_->getLocalViewDevice(Tpetra::Access::ReadWrite);
    auto ColMapColors_Data = ColMapColors_->getLocalViewDevice(Tpetra::Access::ReadWrite);

    auto NewColProfilesData         = NewColProfiles.getLocalViewDevice(Tpetra::Access::ReadWrite);
    auto RowHasColWithSingletonData = RowHasColWithSingleton.getLocalViewDevice(Tpetra::Access::ReadWrite);

    auto lclRowPtr = FullCrsMatrix_->getLocalRowPtrsDevice();
    auto lclColInd = FullCrsMatrix_->getLocalIndicesDevice();

    // Count singleton columns (that were not already counted as singleton rows)
    Kokkos::parallel_reduce(
        "CrsSingletonFilter_LinearProblem:Analyze(find-singleton-column)", range_policy(0, localNumCols),
        KOKKOS_LAMBDA(const size_t j, local_ordinal_type& lclNumSingletonCols) {
          // Check if column is a singleton
          if (ColProfilesData(j, 0) == 1) {
            // i = id of "last" local row with non-zero in this col
            // since this is singletone column , "i" is the row id of the single nonzero entry in this column
            local_ordinal_type i = localRowIDofSingletonColData(j, 0);

            // Check to see if this column already eliminated by the row check above
            //  RowMapColors(i,0) : 0 = ith row was not singleton, 1 = was singleton, 2 = was not singleton, but processed
            //  Multiple singleton columns cannot have nonzero entry in the same row, i
            //   Otherwise, the matrix will be singular
            //  So, rowMapColors(i,0) should be checked by one singleton column (note 1, also see check 2)
            //   not atomic needed
            if (RowMapColors_Data(i, 0) != 1) {  // that row is not singleton, and hence this col has not been removed
              RowMapColors_Data(i, 0) = 2;
              ColMapColors_Data(j, 0) = 1;
              lclNumSingletonCols++;

              // Increment col singleton counter for ith row (Only one j should update (i) because of "note 1")
              RowHasColWithSingletonData(i, 0)++;

              // If we delete a row, we need to keep track of associated column entries that were also deleted
              // in case all entries in a column are eventually deleted, in which case the column should
              // also be deleted.
              // Only one j should execute this loop for a specific "i" (because of "note 1")
              for (size_t k = lclRowPtr(i); k < lclRowPtr(i + 1); k++) {
                NewColProfilesData(lclColInd(k), 0)--;
              }
            }
          }
        },
        localNumSingletonCols_);
    Kokkos::parallel_for(
        "CrsSingletonFilter_LinearProblem:Analyze(find-non-singleton-column)", range_policy(0, localNumCols),
        KOKKOS_LAMBDA(const size_t j) {
          // Check if column is *NOT* a singleton
          if (ColProfilesData(j, 0) != 1) {
            // Since this is not singletone col,
            //  "i" is the "last" local row with non-zero in this col and may not corresponds to the singleton row?
            local_ordinal_type i = localRowIDofSingletonColData(j, 0);

            // Check if some other processor has "previously" eliminated this column
            // (note: multiple singleton columns cannot have nonzero in the same column, colhasrowwithsingleton(j) is 0 or 1)
            //  This column has "one" row, which is singleton
            //   that singleton row has the single nonzero entry at jth col
            //  At jth iteration, if RowMapColors_Data(i, 0) == 1, then it stays as 1
            //             otherwise it has been processed into 2 by previous column, or still 0
            if (ColHasRowWithSingletonData(j, 0) == 1 && RowMapColors_Data(i, 0) != 1) {
              // RowMapColors_Data(i, 0) was 0, or processed by a previous column
              // TODO: check the second and third conditions
              ColMapColors_Data(j, 0) = 1;
            }
          }
        });
  }

  // Next we check to make sure no row is associated with more than one singleton col (check 2)
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(RowHasColWithSingleton.normInf() > 1, std::runtime_error,
                                        "At least one row is associated with two singleton columns, can't handle it.");

  // Generate arrays that keep track of column singleton row, col and pivot info needed for post-solve phase
  CreatePostSolveArrays(localRowIDofSingletonCol, ColProfiles, NewColProfiles, ColHasRowWithSingleton);
  if (run_on_host_) {
    auto RowMapColors_Data = RowMapColors_->getLocalViewHost(Tpetra::Access::ReadWrite);
    for (local_ordinal_type i = 0; i < localNumRows; i++) {
      if (RowMapColors_Data(i, 0) == 2) RowMapColors_Data(i, 0) = 1;  // Convert all eliminated rows to same color
    }
  } else {
    auto RowMapColors_Data = RowMapColors_->getLocalViewDevice(Tpetra::Access::ReadWrite);
    Kokkos::parallel_for(
        "CrsSingletonFilter_LinearProblem:Analyze(convert-rowmap-color)", range_policy(0, localNumRows),
        KOKKOS_LAMBDA(const size_t i) {
          if (RowMapColors_Data(i, 0) == 2) RowMapColors_Data(i, 0) = 1;  // Convert all eliminated rows to same color
        });
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
                       const Teuchos::RCP<vector_type_int>& mapColors, int color,
                       bool locally_sort_gids) {
  Teuchos::RCP<const map_type> reducedMap;
  bool canRunOnHost = std::is_same_v<typename device_type::memory_space, Kokkos::HostSpace>;
  if (run_on_host_ && canRunOnHost) {
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
    reducedMap = createNonContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(
        Teuchos::ArrayView<const GlobalOrdinal>(myReducedGlobalIndices.data(), myReducedGlobalIndices.size()),
        originalMap->getComm());
  } else {
    // Iterate through the global indices to find the ones that match the color
    using IndexListType = Kokkos::View<GlobalOrdinal*, device_type>;
    using functor_type  = GenerateReducedMapFunctor<local_vector_int_type, local_map_type, IndexListType>;
    auto mapColors_Data = mapColors->getLocalViewDevice(Tpetra::Access::ReadWrite);

    // count number of elements with the matching color
    size_t lclNumReducedElements = 0;
    size_t lclNumElements        = originalMap->getLocalNumElements();
    {
      functor_type functor(color, mapColors_Data);
      Kokkos::parallel_reduce(
          "CrsSingletonFilter_LinearProblem:GenerateReduceMap(count)", range_policy(0, lclNumElements),
          functor, lclNumReducedElements);
    }
    // create the list of elements with the matching color
    IndexListType indexList;
    {
      auto lclMap = originalMap()->getLocalMap();
      Kokkos::resize(indexList, lclNumReducedElements);
      functor_type functor(color, mapColors_Data, lclMap, indexList);
      Kokkos::parallel_for(
          "CrsSingletonFilter_LinearProblem:GenerateReduceMap(insert)", range_policy(0, lclNumElements),
          functor);
    }
    if (locally_sort_gids) {
      // locally-sort the matched GIDs
      Kokkos::sort(execution_space(), indexList);
    }

    // Create the reduced map using the collected indices
    const Tpetra::global_size_t INVALID = Tpetra::Details::OrdinalTraits<Tpetra::global_size_t>::invalid();
    reducedMap                          = Teuchos::rcp(new map_type(INVALID, indexList,
                                                                    originalMap->getIndexBase(),
                                                                    originalMap->getComm()));
  }

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

  Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();

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
      ConstructRedistributeExporter(OrigReducedMatrixDomainMap_, ReducedMatrixRangeMap_,
                                    RedistributeDomainExporter_, ReducedMatrixDomainMap_);
    } else {
      ReducedMatrixDomainMap_     = OrigReducedMatrixDomainMap_;
      OrigReducedMatrixDomainMap_ = Teuchos::null;
      RedistributeDomainExporter_ = Teuchos::null;
    }

    // Create pointer to Full RHS, LHS
    Teuchos::RCP<multivector_type> FullRHS = FullProblem()->getRHS();
    Teuchos::RCP<multivector_type> FullLHS = FullProblem()->getLHS();
    local_ordinal_type NumVectors          = FullLHS->getNumVectors();

    // Create importers
    Full2ReducedLHSImporter_ = Teuchos::rcp(new import_type(FullMatrixDomainMap(), ReducedMatrixDomainMap()));
    Full2ReducedRHSImporter_ = Teuchos::rcp(new import_type(FullRHS->getMap(), ReducedMatrixRowMap()));

    // Create storage for temporary X values due to explicit elimination of rows
    tempExportX_ = Teuchos::rcp(new multivector_type(FullMatrixColMap(), NumVectors));

    size_t NumEntries = 0;
    Teuchos::ArrayView<const Scalar> Values;
    Teuchos::Array<GlobalOrdinal> Indices;
    LocalOrdinal localNumRows        = FullMatrix()->getLocalNumRows();
    LocalOrdinal ColSingletonCounter = 0;

    // Check if ColSingletonPivotLIDs_, ColSingletonPivot_ can be accessed on host
    bool canRunOnHost = std::is_same_v<typename device_type::memory_space, Kokkos::HostSpace>;
    if (run_on_host_ && canRunOnHost) {
      // Construct Reduced Matrix
      LocalOrdinal maxNumEntries = FullCrsMatrix()->getLocalMaxNumRowEntries();
      ReducedMatrix_             = Teuchos::rcp(new crs_matrix_type(ReducedMatrixRowMap(), ReducedMatrixColMap(), maxNumEntries));

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
            TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(pivot == zero, std::runtime_error,
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
                TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(pivot == zero, std::runtime_error,
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
    } else {
      // The ReducedMatrix has the row and column indices already removed so cannot "insert" them.
      // Need to filter them to only insert remaining.
      // Filter indices and values

      // Construct the reduced matrix
      // * Only the maps have been setup, and hence the rowptr, colind, nzvals need to be still allocated.
      {
        using graph_type   = typename local_matrix_type::StaticCrsGraphType;
        using row_map_type = typename graph_type::row_map_type::non_const_type;
        using entries_type = typename graph_type::entries_type::non_const_type;
        using values_type  = typename local_matrix_type::values_type;
        using functor_type = ConstructReducedProblemFunctor<local_ordinal_type, global_ordinal_type, local_map_type, local_matrix_type, row_map_type, entries_type, values_type>;

        size_t nnzA           = 1;
        auto lclFullRowMap    = FullMatrixRowMap()->getLocalMap();
        auto lclFullColMap    = FullMatrixColMap()->getLocalMap();
        auto lclReducedRowMap = ReducedMatrixRowMap()->getLocalMap();
        auto lclReducedColMap = ReducedMatrixColMap()->getLocalMap();

        auto lclFullMatrix = FullCrsMatrix_->getLocalMatrixDevice();

        LocalOrdinal localNumReducedRows = ReducedMatrixRowMap()->getLocalNumElements();
        row_map_type rowmap_view("rowmap_view", localNumReducedRows + 1);
        entries_type column_view("colind_view", nnzA);
        values_type values_view("values_view", nnzA);

        // launch functor to count nnz per row
        {
          Kokkos::deep_copy(rowmap_view, 0);

          functor_type functor(lclFullRowMap, lclFullColMap, lclReducedRowMap, lclReducedColMap,
                               lclFullMatrix, rowmap_view, column_view, values_view);
          Kokkos::parallel_reduce(
              "CrsSingletonFilter_LinearProblem:CountReducedProblem", range_policy(0, localNumRows),
              functor, nnzA);
        }

        // convert to rowptrs
        KokkosKernels::Impl::kk_inclusive_parallel_prefix_sum<execution_space>(1 + localNumReducedRows, rowmap_view);

        // insert non-zero entries
        {
          Kokkos::resize(column_view, nnzA);
          Kokkos::resize(values_view, nnzA);

          // functor with new view sizes
          functor_type functor(lclFullRowMap, lclFullColMap, lclReducedRowMap, lclReducedColMap,
                               lclFullMatrix, rowmap_view, column_view, values_view);
          Kokkos::parallel_for(
              "CrsSingletonFilter_LinearProblem:InsertReducedProblem", range_policy(0, localNumRows),
              functor);
        }

        // construct local crsmatrix
        graph_type static_graph(column_view, rowmap_view);
        local_matrix_type crsmat("CrsMatrix", localNumReducedRows, values_view, static_graph);

        // Construct Reduced Matrix with localmatrix + rowmap, colmap, domainmap, & rangemap
        // We have constructed things so that the domain and range of the
        // matrix will have the same map.  If the reduced matrix domain and range maps were not the same, the
        // differences were addressed in the ConstructRedistributeExporter() method
        ReducedMatrix_ = Teuchos::rcp(new crs_matrix_type(crsmat, ReducedMatrixRowMap(), ReducedMatrixColMap(), ReducedMatrixDomainMap(), ReducedMatrixRangeMap()));
      }

      // Solve with singleton part
      {
        using error_code_type = typename Kokkos::View<int*, execution_space>;
        using functor_type    = SolveSingletonProblemFunctor<scalar_type, local_ordinal_type, global_ordinal_type, Node, local_map_type, local_matrix_type,
                                                          local_multivector_type, const_local_multivector_type,
                                                          vector_view_type_int, vector_view_type_scalar, error_code_type>;
        auto lclFullRowMap    = FullMatrixRowMap()->getLocalMap();
        auto lclFullColMap    = FullMatrixColMap()->getLocalMap();
        auto lclReducedRowMap = ReducedMatrixRowMap()->getLocalMap();

        auto lclFullMatrix = FullCrsMatrix_->getLocalMatrixDevice();

        auto localExportX = tempExportX_->getLocalViewDevice(Tpetra::Access::ReadWrite);
        auto localRHS     = FullRHS->getLocalViewDevice(Tpetra::Access::ReadOnly);

        error_code_type error_code("singleton-bound-check", 2);
        Kokkos::deep_copy(error_code, 0);

        // launch functor
        functor_type functor(NumVectors, localNumSingletonCols_, error_code,
                             lclFullRowMap, lclFullColMap, lclReducedRowMap,
                             lclFullMatrix, localExportX, localRHS,
                             ColSingletonColLIDs_, ColSingletonRowLIDs_, ColSingletonPivotLIDs_, ColSingletonPivots_);
        Kokkos::parallel_for(
            "CrsSingletonFilter_LinearProblem:SolveSingletonProblem", range_policy(0, localNumRows),
            functor);
        auto h_error = Kokkos::create_mirror_view(error_code);
        Kokkos::deep_copy(h_error, error_code);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(h_error(0) == 1, std::runtime_error,
                                              "ConstructReducedProblem: Encountered zero row, unable to continue.");  // Should improve this comparison to zero.
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(h_error(0) == 2, std::runtime_error,
                                              "ConstructReducedProblem: Encountered zero column, unable to continue.");  // Should improve this comparison to zero.
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(h_error(0) == 3, std::runtime_error,
                                              "ConstructReducedProblem: Unable to find my column.");
      }
    }

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

  Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();

  if (SingletonsDetected()) {
    // Create pointer to Full RHS, LHS
    Teuchos::RCP<multivector_type> FullRHS = FullProblem()->getRHS();
    Teuchos::RCP<multivector_type> FullLHS = FullProblem()->getLHS();

    local_ordinal_type NumVectors = FullLHS->getNumVectors();
    tempExportX_->putScalar(0.0);

    size_t NumEntries = 0;
    Teuchos::ArrayView<const Scalar> Values;
    Teuchos::Array<GlobalOrdinal> Indices;
    LocalOrdinal localNumRows        = FullMatrix()->getLocalNumRows();
    LocalOrdinal ColSingletonCounter = 0;

    // Check if ColSingletonPivotLIDs_, ColSingletonPivot_ can be accessed on host
    bool canRunOnHost = std::is_same_v<typename device_type::memory_space, Kokkos::HostSpace>;
    if (run_on_host_ && canRunOnHost) {
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
            TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(pivot == zero, std::runtime_error,
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
            TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(pivot == zero, std::runtime_error,
                                                  "Encountered zero column, unable to continue.");  // Should improve this comparison to zero.
            ColSingletonPivots_[ColSingletonCounter] = pivot;
            ColSingletonCounter++;
          }
        }
      }
    } else {
      // Not part of the reduced matrix
      {
        using error_code_type = typename Kokkos::View<int*, execution_space>;
        using functor_type    = SolveSingletonProblemFunctor<scalar_type, local_ordinal_type, global_ordinal_type, Node, local_map_type, local_matrix_type,
                                                          local_multivector_type, const_local_multivector_type,
                                                          vector_view_type_int, vector_view_type_scalar, error_code_type>;

        auto lclFullRowMap    = FullMatrixRowMap()->getLocalMap();
        auto lclFullColMap    = FullMatrixColMap()->getLocalMap();
        auto lclReducedRowMap = ReducedMatrixRowMap()->getLocalMap();

        auto lclFullMatrix = FullCrsMatrix_->getLocalMatrixDevice();

        auto localRHS     = FullRHS->getLocalViewDevice(Tpetra::Access::ReadOnly);
        auto localExportX = tempExportX_->getLocalViewDevice(Tpetra::Access::ReadWrite);

        error_code_type error_code("singleton-bound-check", 2);
        Kokkos::deep_copy(error_code, 0);

        // launch functor
        functor_type functor(NumVectors, localNumSingletonCols_, error_code,
                             lclFullRowMap, lclFullColMap, lclReducedRowMap,
                             lclFullMatrix, localExportX, localRHS,
                             ColSingletonColLIDs_, ColSingletonRowLIDs_, ColSingletonPivotLIDs_, ColSingletonPivots_);
        Kokkos::parallel_reduce(
            "CrsSingletonFilter_LinearProblem:SolveSingletonProblem", range_policy(0, localNumRows),
            functor, ColSingletonCounter);

        // get error
        auto h_error = Kokkos::create_mirror_view(error_code);
        Kokkos::deep_copy(h_error, error_code);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(h_error(0) == 1, std::runtime_error,
                                              "UpdateReducedProblem: Encountered zero row, unable to continue.");  // Should improve this comparison to zero.
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(h_error(0) == 2, std::runtime_error,
                                              "UpdateReducedProblem: Encountered zero column, unable to continue.");  // Should improve this comparison to zero.
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(h_error(0) == 3, std::runtime_error,
                                              "UpdateReducedProblem: Unable to find local column.");
      }

      // Part of the reduced matrix
      {
        using functor_type    = UpdateReducedProblemFunctor<local_ordinal_type, global_ordinal_type, local_map_type, local_matrix_type>;
        auto lclFullRowMap    = FullMatrixRowMap()->getLocalMap();
        auto lclFullColMap    = FullMatrixColMap()->getLocalMap();
        auto lclReducedRowMap = ReducedMatrixRowMap()->getLocalMap();
        auto lclReducedColMap = ReducedMatrixColMap()->getLocalMap();

        auto lclReducedMatrix = ReducedMatrix_->getLocalMatrixDevice();
        auto lclFullMatrix    = FullCrsMatrix_->getLocalMatrixDevice();

        // launch functor
        functor_type functor(lclFullRowMap, lclFullColMap, lclReducedRowMap, lclReducedColMap,
                             lclFullMatrix, lclReducedMatrix);
        Kokkos::parallel_for(
            "CrsSingletonFilter_LinearProblem:UpdateReducedProblem", range_policy(0, localNumRows),
            functor);
      }
    }

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(!(ColSingletonCounter == localNumSingletonCols_),
                                          std::runtime_error, "Sanity Check " + std::to_string(ColSingletonCounter) + " vs " + std::to_string(localNumSingletonCols_) + " (UpdateReducedProblem).");

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

//==============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ConstructRedistributeExporter(Teuchos::RCP<const map_type> SourceMap, Teuchos::RCP<const map_type> TargetMap,
                                  Teuchos::RCP<export_type>& RedistributeExporter,
                                  Teuchos::RCP<const map_type>& RedistributeMap) {
  const char tfecfFuncName[] = "ConstructRedistributeExporter: ";

  auto IndexBase = SourceMap->getIndexBase();
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(IndexBase != TargetMap->getIndexBase(), std::runtime_error,
                                        "Base indices do not match.");

  auto Comm = TargetMap->getComm();

  size_t TargetNumMyElements = TargetMap->getLocalNumElements();
  size_t SourceNumMyElements = SourceMap->getLocalNumElements();

  // ContiguousTargetMap has the same number of elements per PE as TargetMap, but uses contiguous indexing
  Teuchos::RCP<const map_type> ContiguousTargetMap =
      Teuchos::rcp(new const map_type(
          Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
          TargetNumMyElements, IndexBase, Comm));

  // Same for ContiguousSourceMap
  Teuchos::RCP<const map_type> ContiguousSourceMap =
      Teuchos::rcp(new const map_type(
          Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
          SourceNumMyElements, IndexBase, Comm));

  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ContiguousSourceMap->getGlobalNumElements() != ContiguousTargetMap->getGlobalNumElements(),
      std::runtime_error, "Global number of elements do not match.");

  // Retrieve the global indices owned by the calling process as a Kokkos::View
  auto SourceMapMyGlobalIndices = SourceMap->getMyGlobalIndices();

  // Create a Tpetra::Vector to hold the global indices
  Teuchos::RCP<Tpetra::MultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>> SourceIndices =
      Teuchos::rcp(new Tpetra::MultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>(ContiguousSourceMap, 1));

  // Fill the Tpetra::Vector with the global indices
  for (size_t i = 0; i < static_cast<size_t>(SourceMapMyGlobalIndices.size()); ++i) {
    SourceIndices->replaceLocalValue(i, 0, static_cast<GlobalOrdinal>(SourceMapMyGlobalIndices[i]));
  }

  // Get the map associated with the vector
  auto map = SourceIndices->getMap();

  Teuchos::RCP<Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node>> Exporter =
      Teuchos::rcp(new Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node>(ContiguousSourceMap, ContiguousTargetMap));

  Teuchos::RCP<Tpetra::MultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>> TargetIndices =
      Teuchos::rcp(new Tpetra::MultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>(ContiguousTargetMap, 1));

  TargetIndices->doExport(*SourceIndices, *Exporter, Tpetra::INSERT);

  // Create a new map that describes how the Source MultiVector should be laid out so that it has
  // the same number of elements on each processor as the TargetMap
  RedistributeMap = Tpetra::createNonContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(
      Teuchos::ArrayView<const GlobalOrdinal>(TargetIndices->getData(0).getRawPtr(), TargetIndices->getLocalLength()),
      Comm);

  RedistributeExporter =
      Teuchos::rcp(new Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node>(SourceMap, RedistributeMap));

  return;
}

//==============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ComputeFullSolution() {
  if (SingletonsDetected()) {
    Teuchos::RCP<multivector_type> FullLHS = FullProblem()->getLHS();
    Teuchos::RCP<multivector_type> FullRHS = FullProblem()->getRHS();

    tempX_->putScalar(Scalar(0.0));
    tempExportX_->putScalar(Scalar(0.0));

    // Inject values that the user computed for the reduced problem into the full solution vector
    tempX_->doExport(*ReducedLHS_, *Full2ReducedLHSImporter_, Tpetra::ADD);
    FullLHS->update(Scalar(1.0), *tempX_, Scalar(1.0));  // FullLHS = 1.0 * tempX_ + 1.0 * FullLHS

    // Next we will use our full solution vector which is populated with pre-filter solution
    // values and reduced system solution values to compute the sum of the row contributions
    // that must be subtracted to get the post-filter solution values
    FullMatrix()->apply(*FullLHS, *tempB_);

    // Finally we loop through the local rows that were associated with column singletons and compute the
    // solution for these equations.
    size_t NumVectors = tempB_->getNumVectors();
    bool canRunOnHost = std::is_same_v<typename device_type::memory_space, Kokkos::HostSpace>;
    if (run_on_host_ && canRunOnHost) {
      for (int k = 0; k < localNumSingletonCols_; k++) {
        LocalOrdinal i = ColSingletonRowLIDs_[k];
        LocalOrdinal j = ColSingletonColLIDs_[k];
        Scalar pivot   = ColSingletonPivots_[k];
        for (size_t jj = 0; jj < NumVectors; jj++) {
          auto tempExportXData = tempExportX_->getDataNonConst(jj);
          auto FullRHSData     = FullRHS->getData(jj);
          auto tempBData       = tempB_->getData(jj);
          tempExportXData[j]   = (FullRHSData[i] - tempBData[i]) / pivot;
        }
      }
    } else {
      using functor_type = ComputeFullSolutionFunctor<vector_view_type_int, vector_view_type_scalar,
                                                      local_multivector_type, const_local_multivector_type>;
      auto localB        = tempB_->getLocalViewDevice(Tpetra::Access::ReadOnly);
      auto localRHS      = FullRHS->getLocalViewDevice(Tpetra::Access::ReadOnly);
      auto localX        = tempExportX_->getLocalViewDevice(Tpetra::Access::ReadWrite);

      functor_type functor(ColSingletonRowLIDs_, ColSingletonColLIDs_, ColSingletonPivots_, localX, localRHS, localB);
      Kokkos::parallel_for(
          "CrsSingletonFilter_LinearProblem:ComputeFullSolution", range_policy(0, localNumSingletonCols_),
          functor);
    }

    // Finally, insert values from post-solve step and we are done!!!!
    auto importer = FullMatrix()->getGraph()->getImporter();
    if (importer != Teuchos::null) {
      tempX_->doExport(*tempExportX_, *importer, Tpetra::ADD);
    } else {
      tempX_->update(Scalar(1.0), *tempExportX_, Scalar(0.0));
    }

    FullLHS->update(Scalar(1.0), *tempX_, Scalar(1.0));
  }

  return;
}
//==============================================================================

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    InitFullMatrixAccess() {
  localMaxNumRowEntries_ = FullMatrix()->getLocalMaxNumRowEntries();

  // Cast to CrsMatrix, if possible.  Can save some work.
  FullCrsMatrix_         = Teuchos::rcp_dynamic_cast<crs_matrix_type>(FullMatrix());
  FullMatrixIsCrsMatrix_ = (FullCrsMatrix() != Teuchos::null);  // Pointer is non-null if cast worked

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
    Values     = Teuchos::ArrayView<const Scalar>(reinterpret_cast<const Scalar*>(rowValues.data()), rowValues.size());
    Indices    = Teuchos::ArrayView<const LocalOrdinal>(localIndices.data(), localIndices.size());

  } else {  // Copy of current row
    typename Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::nonconst_local_inds_host_view_type localIndicesCopy("localIndicesCopy", FullMatrix()->getLocalMaxNumRowEntries());
    typename Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::nonconst_values_host_view_type rowValuesCopy("rowValuesCopy", FullMatrix()->getLocalMaxNumRowEntries());

    FullMatrix()->getLocalRowCopy(Row, localIndicesCopy, rowValuesCopy, NumIndices);

    Values  = Teuchos::ArrayView<const Scalar>(reinterpret_cast<const Scalar*>(rowValuesCopy.data()), NumIndices);
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
  Values = Teuchos::ArrayView<const Scalar>(reinterpret_cast<const Scalar*>(RowValues.data()), RowValues.size());

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
  Kokkos::resize(ColSingletonRowLIDs_, localNumSingletonCols_);
  Kokkos::resize(ColSingletonColLIDs_, localNumSingletonCols_);
  Kokkos::resize(ColSingletonPivotLIDs_, localNumSingletonCols_);
  Kokkos::resize(ColSingletonPivots_, localNumSingletonCols_);

  int NumMyColSingletonstmp = 0;

  // Check if ColSingletonRowLIDs_, ColSingletonColLIDs_ can be accessed on host
  bool canRunOnHost = std::is_same_v<typename device_type::memory_space, Kokkos::HostSpace>;
  if (run_on_host_ && canRunOnHost) {
    // Register singleton columns (that were not already counted as singleton rows)
    // Check to see if any columns disappeared because all associated rows were eliminated
    local_ordinal_type localNumCols   = FullMatrix()->getLocalNumCols();
    auto localRowIDofSingletonColData = localRowIDofSingletonCol.getLocalViewHost(Tpetra::Access::ReadOnly);
    auto ColProfilesData              = ColProfiles.getLocalViewHost(Tpetra::Access::ReadOnly);
    auto RowMapColors_Data            = RowMapColors_->getLocalViewHost(Tpetra::Access::ReadOnly);
    auto NewColProfilesData           = NewColProfiles.getLocalViewHost(Tpetra::Access::ReadOnly);
    auto ColHasRowWithSingletonData   = ColHasRowWithSingleton.getLocalViewHost(Tpetra::Access::ReadOnly);
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

  } else {
    {
      // Register singleton columns (that were not already counted as singleton rows)
      // Check to see if any columns disappeared because all associated rows were eliminated
      local_ordinal_type localNumCols   = FullMatrix()->getLocalNumCols();
      auto localRowIDofSingletonColData = localRowIDofSingletonCol.getLocalViewDevice(Tpetra::Access::ReadOnly);
      auto ColProfilesData              = ColProfiles.getLocalViewDevice(Tpetra::Access::ReadOnly);
      auto RowMapColors_Data            = RowMapColors_->getLocalViewDevice(Tpetra::Access::ReadOnly);
      auto NewColProfilesData           = NewColProfiles.getLocalViewDevice(Tpetra::Access::ReadOnly);
      auto ColHasRowWithSingletonData   = ColHasRowWithSingleton.getLocalViewDevice(Tpetra::Access::ReadOnly);
      auto ColMapColors_Data            = ColMapColors_->getLocalViewDevice(Tpetra::Access::ReadWrite);

      // counter to be atomic-incremented
      Kokkos::View<int*, execution_space> numSingletonCols("numSingletonCols", 1);
      Kokkos::deep_copy(numSingletonCols, 0);

      // launch functor
      CreatePostSolveArraysFunctor functor(localRowIDofSingletonColData,
                                           ColSingletonRowLIDs_, ColSingletonColLIDs_,
                                           NewColProfilesData, ColProfilesData, ColHasRowWithSingletonData,
                                           ColMapColors_Data, RowMapColors_Data, numSingletonCols);
      Kokkos::parallel_for(
          "CrsSingletonFilter_LinearProblem:CreatePostSolveArrays", range_policy(0, localNumCols),
          functor);

      // get the counter
      auto h_numSingletonCols = Kokkos::create_mirror_view(numSingletonCols);
      Kokkos::deep_copy(h_numSingletonCols, numSingletonCols);
      NumMyColSingletonstmp = h_numSingletonCols(0);
    }
  }

  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(NumMyColSingletonstmp != localNumSingletonCols_,
                                        std::runtime_error, "Sanity Check " + std::to_string(NumMyColSingletonstmp) + " vs " + std::to_string(localNumSingletonCols_) + " (CreatePostSolveArrays).");

  Kokkos::Experimental::sort_by_key(execution_space(), ColSingletonRowLIDs_, ColSingletonColLIDs_);

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
