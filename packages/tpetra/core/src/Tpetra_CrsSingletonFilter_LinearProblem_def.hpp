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

//#include "Epetra_ConfigDefs.h"
//#include "Epetra_Map.h"
//#include "Epetra_Util.h"
//#include "Epetra_Export.h"
//#include "Epetra_Import.h"
//#include "Epetra_MultiVector.h"
//#include "Epetra_Vector.h"
//#include "Epetra_GIDTypeVector.h"
//#include "Epetra_Comm.h"
//#include "Epetra_LinearProblem.h"
//#include "Epetra_MapColoring.h"
//#include "EpetraExt_CrsSingletonFilter_LinearProblem.h"

#include "Tpetra_CrsSingletonFilter_LinearProblem_decl.hpp"


namespace Tpetra {

  #define PRINT(outputRank_, message) \
    do { \
        if (outputRank_) { \
            std::cout << "petra : " << message << std::endl; \
        } \
        Tpetra::getDefaultComm()->barrier(); \
    } while (0)

  #define PRINT_NB(outputRank_, message) \
    do { \
        if (outputRank_) { \
            std::cout << "petra : " << message << std::endl; \
        } \
    } while (0)

  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
  void printMultiVector(
    const std::string& label,
    const Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& multiVector)
  {
    // Print the label
    std::cout << "petra : " << label << std::endl;
  
    // Loop over columns (vectors) in the MultiVector
    for (size_t col = 0; col < multiVector->getNumVectors(); ++col) {
      auto data = multiVector->getData(col); // Get the data for column `col`
      std::cout << "Column " << col << ": ";
      
      // Loop over rows (entries) in the column
      for (size_t row = 0; row < multiVector->getLocalLength(); ++row) {
          std::cout << "(" << row << ", " << data[row] << ") ";
      }
      std::cout << std::endl;
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  PrintTpetraRowMatrix(const Teuchos::RCP<const row_matrix_type>& matrix) {
    auto comm = matrix->getComm();
    auto numRows = matrix->getLocalNumRows();
    auto maxNumEntries = matrix->getLocalMaxNumRowEntries();

    nonconst_local_inds_host_view_type indices("indices", maxNumEntries);
    nonconst_values_host_view_type values("values", maxNumEntries);

    for (size_t i = 0; i < numRows; ++i) {
        size_t numEntries;
        matrix->getLocalRowCopy(i, indices, values, numEntries);

        std::cout << "Row " << i << ": ";
        for (size_t j = 0; j < numEntries; ++j) {
            std::cout << "(" << indices(j) << ", " << values(j) << ") ";
        }
        std::cout << std::endl;
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  PrintTpetraCrsMatrix(const Teuchos::RCP<const crs_matrix_type>& matrix) {
    auto comm = matrix->getComm();
    auto numRows = matrix->getLocalNumRows();

    for (size_t i = 0; i < numRows; ++i) {

        size_t numEntries = matrix->getNumEntriesInLocalRow(i); // Get the actual number of entries in the row

        nonconst_local_inds_host_view_type indices("indices", numEntries);
        nonconst_values_host_view_type values("values", numEntries);

        matrix->getLocalRowCopy(i, indices, values, numEntries);

        std::cout << "Row " << i << ": ";
        for (size_t j = 0; j < numEntries; ++j) {
            std::cout << "(" << indices(j) << ", " << values(j) << ") ";
        }
        std::cout << std::endl;
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  printMap(const Teuchos::RCP<const map_type> & map) const {
        Teuchos::RCP<Teuchos::FancyOStream> fancyOut = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
        map->describe(*fancyOut, Teuchos::VERB_HIGH);
        std::cout << "         MyPID    " << "       Local Index  " << "      Global Index  " << std::endl;
        for (LocalOrdinal ii = 0; ii < static_cast<int>(map->getLocalNumElements()); ++ii) {
            GlobalOrdinal globalIndex = map->getGlobalElement(ii);
            std::cout << std::setw(14) << std::right << FullMatrix_->getComm()->getRank() << "    "
                      << std::setw(14) << std::right << ii << "    "
                      << std::setw(14) << std::right << globalIndex << "    "
                      << std::endl;
        }
    }

//==============================================================================

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  CrsSingletonFilter_LinearProblem( bool verbose )
    : globalNumSingletonRows_(Teuchos::as<local_ordinal_type>(0)),
      globalNumSingletonCols_(Teuchos::as<local_ordinal_type>(0)),
      RatioOfDimensions_(0.0),
      RatioOfNonzeros_(0.0),
      HaveReducedProblem_(false),
      AnalysisDone_(false),
      SymmetricElimination_(true),
      localMaxNumRowEntries_(Teuchos::as<local_ordinal_type>(0)),
      FullMatrixIsCrsMatrix_(false),
      verbose_(verbose)
  {
  }

//==============================================================================

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  ~CrsSingletonFilter_LinearProblem()
  {
//    if (ReducedLHS_!=0) delete ReducedLHS_;
//    if (ReducedRHS_!=0) delete ReducedRHS_;
//    if (ReducedMatrixDomainMap_!=ReducedMatrixColMap_) delete ReducedMatrixDomainMap_;
//    if (OrigReducedMatrixDomainMap_!=ReducedMatrixColMap_ &&
//        OrigReducedMatrixDomainMap_!=0) delete OrigReducedMatrixDomainMap_;
//    if (ReducedMatrixRangeMap_!=ReducedMatrixRowMap_) delete ReducedMatrixRangeMap_;
//    if (ReducedMatrixRowMap_!=0) delete ReducedMatrixRowMap_;
//    if (ReducedMatrixColMap_!=0) delete ReducedMatrixColMap_;
//    if (Full2ReducedRHSImporter_!=0) delete Full2ReducedRHSImporter_;
//    if (Full2ReducedLHSImporter_!=0) delete Full2ReducedLHSImporter_;
//    if (RedistributeDomainExporter_!=0) delete RedistributeDomainExporter_;
//    if (RowMapColors_!=0) delete RowMapColors_;
//    if (ColMapColors_!=0) delete ColMapColors_;
//  
//    if (ColSingletonRowLIDs_ != 0) delete [] ColSingletonRowLIDs_;
//    if (ColSingletonColLIDs_ != 0) delete [] ColSingletonColLIDs_;
//    if (ColSingletonPivotLIDs_ != 0) delete [] ColSingletonPivotLIDs_;
//    if (ColSingletonPivots_ != 0) delete [] ColSingletonPivots_;
//    if (tempExportX_ != 0) delete tempExportX_;
//    if (Indices_int_ != 0) delete [] Indices_int_;
//  #ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
//    if (Indices_LL_ != 0) delete [] Indices_LL_;
//  #endif
//    if (tempX_ != 0) delete tempX_;
//    if (tempB_ != 0) delete tempB_;
 }

//==============================================================================

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::NewType
  CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  operator()(const CrsSingletonFilter_LinearProblem::OriginalType & originalLinearProblem)
  {
    const size_t outRank = 0;
    auto Comm = Tpetra::getDefaultComm ();
    const size_t myRank = Comm->getRank();
    outputRank_ = outRank == myRank;

    analyze( originalLinearProblem );
    //return originalLinearProblem;
    return construct();
  }

//==============================================================================

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  analyze(const CrsSingletonFilter_LinearProblem::OriginalType & originalLinearProblem )
{
  this->origObj_ = originalLinearProblem;

  FullMatrix_ = originalLinearProblem->getMatrix();

  Analyze( FullMatrix_ );

  if ( verbose_ && FullMatrix_->getComm()->getRank()==0 ) {
    std::cout << "\nAnalyzed Singleton Problem:\n";
    std::cout << "---------------------------\n";
  }
  if ( SingletonsDetected() ) {
    if ( verbose_ && FullMatrix_->getComm()->getRank()==0 ) {
      std::cout << "Singletons Detected!" << std::endl;;
      std::cout << "Num Singletons:      " << NumSingletons() << std::endl;
    }
  }
  else {
    if ( verbose_ && FullMatrix_->getComm()->getRank()==0 )
        std::cout << "No Singletons Detected!" << std::endl;
  }
/*
    std::cout << "List of Row Singletons:\n";
    int * slist = RowMapColors_->ColorLIDList(1);
    for( int i = 0; i < RowMapColors_->NumElementsWithColor(1); ++i )
      std::cout << slist[i] << " ";
    std::cout << "\n";
    std::cout << "List of Col Singletons:\n";
    slist = ColMapColors_->ColorLIDList(1);
    for( int i = 0; i < ColMapColors_->NumElementsWithColor(1); ++i )
      std::cout << slist[i] << " ";
    std::cout << "\n";
*/
  if ( verbose_ && FullMatrix_->getComm()->getRank()==0 )
    std::cout << "---------------------------\n\n";

  return;
}

//==============================================================================

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::NewType
  CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  construct()
  {
    const char tfecfFuncName[] = "construct: ";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(this->origObj_==Teuchos::null, std::runtime_error,
      "Linear problem (orgiObj_) have not been set.  Call analyze() first.");

    ConstructReducedProblem( this->origObj_ );
  
    this->newObj_ = ReducedProblem();
  
    if( verbose_ && SingletonsDetected() && FullMatrix_->getComm()->getRank()==0 )
    {
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
  fwd()
  {
    UpdateReducedProblem( FullProblem_ );
  }

//==============================================================================

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  rvs()
  {
// GOT TO HERE!!
    //ComputeFullSolution();
  }

//==============================================================================

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  Analyze(const Teuchos::RCP<row_matrix_type> & fullMatrix)
  {
    PRINT(outputRank_, "Analyze() A");
    const char tfecfFuncName[] = "Analyze: ";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(AnalysisDone_, std::runtime_error,
      "Analysis() already done once.  Cannot do it again.");
  
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(fullMatrix == Teuchos::null,
      std::runtime_error, "Input matrix is Teuchos::null.");
  
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(fullMatrix->getGlobalNumRows()==0,
      std::runtime_error, "Full matrix has zero dimension.");
  
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(fullMatrix->getGlobalNumEntries()==0,
      std::runtime_error, "Full matrix has no nonzero terms.");
  
    FullMatrix_ = fullMatrix;

    PRINT(outputRank_, "Analyze() B");
    //std::cout << "FullMatrix_ = " << FullMatrix_ << std::endl;
    //std::cout << "FullMatrix_->getColMap() = " << FullMatrix_->getColMap() << std::endl;

    // First check for columns with single entries and find columns with singleton rows
    vector_type_int ColProfiles(FullMatrixColMap()); ColProfiles.putScalar(0);
    vector_type_int ColHasRowWithSingleton(FullMatrixColMap()); ColHasRowWithSingleton.putScalar(0);
  
    // localRowIDofSingletonCol[j] will contain the local row ID associated with the jth column,
    // if the jth col has a single entry
    vector_type_LO localRowIDofSingletonCol(FullMatrixColMap());
    localRowIDofSingletonCol.putScalar(-1);
  
    PRINT(outputRank_, "Analyze() C");

    // Define MapColoring objects
    //RowMapColors_ = new Epetra_MapColoring(FullMatrixRowMap());  // Initial colors are all 0
    //ColMapColors_ = new Epetra_MapColoring(FullMatrixColMap());
    RowMapColors_ = Teuchos::rcp(new vector_type_int(FullMatrixRowMap()));  // Initial colors are all 0
    ColMapColors_ = Teuchos::rcp(new vector_type_int(FullMatrixColMap()));
    //Epetra_MapColoring & rowMapColors = *RowMapColors_;
    //Epetra_MapColoring & colMapColors = *ColMapColors_;
  
    local_ordinal_type localNumRows = FullMatrix_->getLocalNumRows();
    local_ordinal_type localNumCols = FullMatrix_->getLocalNumCols();
  
    // Set up for accessing full matrix.  Will do so row-by-row.
    InitFullMatrixAccess();
  
    // Scan matrix for singleton rows, build up column profiles
    size_t NumIndices = 1;
    PRINT(outputRank_, "Analyze() Initialize NumIndices = " << NumIndices);
    //int * localIndices;
    //nonconst_local_inds_host_view_type localIndices;
    Teuchos::Array<local_ordinal_type> localIndices;
    localNumSingletonRows_ = 0;

    PRINT(outputRank_, "Analyze() D");

    auto ColProfilesData = ColProfiles.getLocalViewHost(Tpetra::Access::ReadWrite);
    auto localRowIDofSingletonColData = localRowIDofSingletonCol.getLocalViewHost(Tpetra::Access::ReadWrite);
    auto ColHasRowWithSingletonData = ColHasRowWithSingleton.getLocalViewHost(Tpetra::Access::ReadWrite);
    auto RowMapColors_Data = RowMapColors_->getLocalViewHost(Tpetra::Access::ReadWrite);
    auto ColMapColors_Data = ColMapColors_->getLocalViewHost(Tpetra::Access::ReadWrite);

    PRINT(outputRank_, "Analyze() E");

    for (int i=0; i<localNumRows; i++) {
      PRINT(outputRank_, "Analyze() E  i = " << i);
      PRINT(outputRank_, "Analyze() E  NumIndices = " << NumIndices);
      // Get ith row
      GetRow(i, NumIndices, localIndices);
      PRINT(outputRank_, "Analyze() E  NumIndices = " << NumIndices);
      if (outputRank_) std::cout << "petra : Analyze() E      localIndices = ";
      for (int ii = 0; ii < localIndices.size(); ++ii) {
        std::cout << localIndices[ii] << " ";
      } std::cout << std::endl;
      PRINT(outputRank_, "Analyze() F");
      for (size_t j=0; j<NumIndices; j++) {
        PRINT(outputRank_, "Analyze() F    j = " << j);
        local_ordinal_type ColumnIndex = localIndices[j];
        PRINT(outputRank_, "Analyze() F    ColumnIndex = " << ColumnIndex);
        PRINT(outputRank_, "Analyze() F    ColProfiles[ColumnIndex] = " << ColProfilesData(ColumnIndex,0));
        ColProfilesData(ColumnIndex,0)++; // Increment column count
        PRINT(outputRank_, "Analyze() F    ColProfiles[ColumnIndex] = " << ColProfilesData(ColumnIndex,0));

        // Record local row ID for current column
        // will use to identify row to eliminate if column is a singleton
        localRowIDofSingletonColData(ColumnIndex,0) = i;
        PRINT(outputRank_, "Analyze() F    RowIDs[ColumnIndex] = " << localRowIDofSingletonColData(ColumnIndex,0));
      }
      PRINT(outputRank_, "Analyze() G");
      // If row has single entry, color it and associated column with color=1
      if (NumIndices == 1) {
        int j2 = localIndices[0];
        ColHasRowWithSingletonData(j2,0)++;
        RowMapColors_Data(i,0) = 1;
        ColMapColors_Data(j2,0) = 1;
        localNumSingletonRows_++;
      }
      PRINT(outputRank_, "Analyze() H");
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
/// auto importer = FullMatrix_->getImporter();
    auto importer = FullCrsMatrix_->getCrsGraph()->getImporter();
    if (importer != Teuchos::null) {
      vector_type_int tmpVec(FullMatrixDomainMap()); // Use for gather/scatter of column vectors
      tmpVec.putScalar(0);
      tmpVec.doExport(ColProfiles, *importer, ADD);
      ColProfiles.doImport(tmpVec, *importer, INSERT);
  
      tmpVec.putScalar(0);
      tmpVec.doExport(ColHasRowWithSingleton, *importer, ADD);
      ColHasRowWithSingleton.doImport(tmpVec, *importer, INSERT);
    }


    PRINT(outputRank_, "Analyze() I    ColProfiles");
    PRINT(outputRank_, "Analyze() I    NumIndices = " << NumIndices);
    for (size_t j=0; j<NumIndices; j++) {
      PRINT(outputRank_, "Analyze() I    j = " << j);
      local_ordinal_type ColumnIndex = localIndices[j];
      PRINT(outputRank_, "Analyze() I    ColumnIndex = " << ColumnIndex);
      PRINT(outputRank_, "Analyze() I    ColProfiles[ColumnIndex] = " << ColProfilesData(ColumnIndex,0));
    }
    PRINT(outputRank_, "Analyze() I    ColHasRowWithSingleton");
    PRINT(outputRank_, "Analyze() I    NumIndices = " << NumIndices);
    for (size_t j=0; j<NumIndices; j++) {
      PRINT(outputRank_, "Analyze() I    j = " << j);
      local_ordinal_type ColumnIndex = localIndices[j];
      PRINT(outputRank_, "Analyze() I    ColumnIndex = " << ColumnIndex);
      PRINT(outputRank_, "Analyze() I    ColHasRowWithSingleton[ColumnIndex] = " << ColProfilesData(ColumnIndex,0));
    }

    // ColProfiles now contains the nonzero column entry count for all columns that have
    // an entry on this processor.
    // ColHasRowWithSingleton now contains a count of singleton rows associated with the corresponding
    // local column.  Next we check to make sure no column is associated with more than one singleton row.
  
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(ColHasRowWithSingleton.normInf() > 1, std::runtime_error,
      "At least one column is associated with two singleton rows, can't handle it.");
    PRINT(outputRank_, "Analyze() J    ColHasRowWithSingleton MaxValue = " << ColHasRowWithSingleton.normInf());

    vector_type_int RowHasColWithSingleton(FullMatrix_->getRowMap()); // Use to check for errors
    RowHasColWithSingleton.putScalar(0);
    auto RowHasColWithSingletonData = RowHasColWithSingleton.getLocalViewHost(Tpetra::Access::ReadWrite);
  
    localNumSingletonCols_ = 0;
    // Count singleton columns (that were not already counted as singleton rows)
    for (local_ordinal_type j = 0; j<localNumCols; j++) {
      local_ordinal_type i2 = localRowIDofSingletonColData(j,0);
      PRINT(outputRank_, "Analyze() K    i2 = " << i2);
      // Check if column is a singleton
      if (ColProfilesData(j, 0) == 1) {
        // Check to make sure RowID is not invalid
  //      assert(i!=-1);
        // Check to see if this column already eliminated by the row check above
        if (RowMapColors_Data(i2, 0) != 1) {
          RowHasColWithSingletonData(i2, 0)++; // Increment col singleton counter for ith row
          RowMapColors_Data(i2, 0) = 2; // Use 2 for now, to distinguish between row eliminated directly or via column singletons
          ColMapColors_Data(j, 0) = 1;
          localNumSingletonCols_++;
          // If we delete a row, we need to keep track of associated column entries that were also deleted
          // in case all entries in a column are eventually deleted, in which case the column should
          // also be deleted.
          GetRow(i2, NumIndices, localIndices);
          for (size_t jj=0; jj<NumIndices; jj++) {
            NewColProfilesData(localIndices[jj], 0)--;
          }
        }
      }
      // Check if some other processor eliminated this column
      else if (ColHasRowWithSingletonData(j, 0)==1 && RowMapColors_Data(i2, 0)!=1) {
          ColMapColors_Data(j, 0) = 1;
      }
    }
    for (int j=0; j<localNumCols; j++) {
      PRINT(outputRank_, "Analyze() DEBUG    ColProfiles["<<j<<"] = " << ColProfilesData(j,0));
    }
  
     
    PRINT(outputRank_, "Analyze() K    NumMyColSingletons_ = " << localNumSingletonCols_);

    auto size1 = RowHasColWithSingletonData.size();
    for (int j=0; j<int(size1); j++) {
      if (RowHasColWithSingletonData(j,0) == 0) continue;
      PRINT_NB(outputRank_, "Analyze() K    RowHasColWithSingleton["<<j<<"] = " << RowHasColWithSingletonData(j,0));
    }
    size1 = RowMapColors_Data.size();
    for (int j=0; j<int(size1); j++) {
      if (RowMapColors_Data(j,0) != 0) {
        PRINT_NB(outputRank_, "Analyze() K    rowMapColors["<<j<<"] = " << RowMapColors_Data(j,0));
      }
      if (ColMapColors_Data(j,0) != 0) {
        PRINT_NB(outputRank_, "Analyze() K    colMapColors["<<j<<"] = " << ColMapColors_Data(j,0));
      }
    }
    size1 = RowHasColWithSingletonData.size();
    for (int j=0; j<int(size1); j++) {
      PRINT(outputRank_, "Analyze() K    NewColProfiles["<<j<<"] = " << NewColProfilesData(j,0));
    }


    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(RowHasColWithSingleton.normInf() > 1, std::runtime_error,
      "At least one row is associated with two singleton columns, can't handle it.");
  
    // Generate arrays that keep track of column singleton row, col and pivot info needed for post-solve phase
    CreatePostSolveArrays(localRowIDofSingletonCol, ColProfiles, NewColProfiles, ColHasRowWithSingleton);
  
    for (local_ordinal_type j = 0; j<localNumCols; j++) {
      PRINT(outputRank_, "Analyze() L    RowIDs["<<j<<"] = " << localRowIDofSingletonColData(j,0));
      PRINT(outputRank_, "Analyze() L    ColProfiles["<<j<<"] = " << ColProfilesData(j,0));
    }
    size1 = RowHasColWithSingletonData.size();
    for (int j=0; j<int(size1); j++) {
      PRINT(outputRank_, "Analyze() L    NewColProfiles["<<j<<"] = " << NewColProfilesData(j,0));
    }
    size1 = ColHasRowWithSingletonData.size();
    for (int j=0; j<int(size1); j++) {
      PRINT(outputRank_, "Analyze() L    ColHasRowWithSingleton["<<j<<"] = " << ColHasRowWithSingletonData(j,0));
    }


    for (local_ordinal_type i = 0; i<localNumRows; i++) {
      if (RowMapColors_Data(i, 0) == 2) RowMapColors_Data(i, 0) = 1; // Convert all eliminated rows to same color
    }
    size1 = RowMapColors_Data.size();
    for (int j=0; j<int(size1); j++) {
      if (RowMapColors_Data(j,0) != 0) {
        PRINT_NB(outputRank_, "Analyze() L    rowMapColors["<<j<<"] = " << RowMapColors_Data(j,0));
      }
    }
  
    const Teuchos::Ptr<local_ordinal_type> gRowsPtr(&globalNumSingletonRows_);
    const Teuchos::Ptr<local_ordinal_type> gColsPtr(&globalNumSingletonCols_);
    auto rowComm = FullMatrix_->getRowMap()->getComm();
    Teuchos::reduceAll<local_ordinal_type,local_ordinal_type>(*rowComm, Teuchos::REDUCE_SUM, localNumSingletonRows_, gRowsPtr);
    Teuchos::reduceAll<local_ordinal_type,local_ordinal_type>(*rowComm, Teuchos::REDUCE_SUM, localNumSingletonCols_, gColsPtr);

    PRINT(outputRank_, "Analyze() M    NumGlobalRowSingletons_ = " << globalNumSingletonRows_);
    PRINT(outputRank_, "Analyze() M    NumGlobalColSingletons_ = " << globalNumSingletonCols_);

    AnalysisDone_ = true;
    return;
  }

//==============================================================================

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> 
CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
GenerateReducedMap(const Teuchos::RCP<const map_type>& originalMap, 
                   const Teuchos::RCP<vector_type_int>& mapColors, int color)
{
    // Vector to hold the reduced global indices
    std::vector<GlobalOrdinal> myReducedGlobalIndices;

    // Get the global IDs owned by the calling process
    auto myGlobalIndices = originalMap->getMyGlobalIndices();

    // Iterate through the global indices to find the ones that match the color
    for (size_t i = 0; i < myGlobalIndices.size(); i++) {
        // Access the value in mapColors using the appropriate method
        int colorValue = mapColors->getData()[i]; // Use getData() to access the vector data
        if (colorValue == color) {
            myReducedGlobalIndices.push_back(myGlobalIndices[i]);
        }
    }

    // Create the reduced map using the collected indices
    Teuchos::RCP<const map_type> reducedMap = createNonContigMap<LocalOrdinal, GlobalOrdinal>(
        Teuchos::ArrayView<const GlobalOrdinal>(myReducedGlobalIndices.data(), myReducedGlobalIndices.size()),
        originalMap->getComm()
    );

    return reducedMap;

//    Teuchos::ArrayView<const global_ordinal_type> myReducedGlobalIndices;
//  
//    // Get the global IDs owned by the calling process
//    auto myGlobalIndices = originalMap->getMyGlobalIndices();
//    for (unsigned i = 0; i < myGlobalIndices->size(); i++) {
//      if (mapColors(i) == color) myReducedGlobalIndices.push_back(myGlobalIndices[i]);
//    }
//    Teuchos::RCP<const map_type> reducedMap = createNonContigMap(myReducedGlobalIndices, originalMap->getComm();
//    return reducedMap;
  }

//==============================================================================

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  ConstructReducedProblem(const Teuchos::RCP<linear_problem_type> & Problem) {

    PRINT(outputRank_, "ConstructReducedProblem() A");
    const char tfecfFuncName[] = "ConstructReducedProblem: ";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(HaveReducedProblem_, std::runtime_error,
      "Already have a reduced problem.  Cannot do it again.");
  
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(Problem == Teuchos::null,
      std::runtime_error, "Problem is Teuchos::null.");
  
    FullProblem_ = Problem;
    FullMatrix_ = Teuchos::rcp_dynamic_cast<row_matrix_type>(Problem->getMatrix());
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(FullMatrix_ == Teuchos::null,
      std::runtime_error, "Need a RowMatrix.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(Problem->getRHS() == Teuchos::null,
      std::runtime_error, "Need a RHS.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(Problem->getLHS() == Teuchos::null,
      std::runtime_error, "Need a LHS.");

    PRINT(outputRank_, "ConstructReducedProblem() A    SingletonsDetected() = " << SingletonsDetected());
    // Generate reduced row and column maps
    if ( SingletonsDetected() ) {
  
      //Epetra_MapColoring & RowMapColors_ = *RowMapColors_;
      //Epetra_MapColoring & ColMapColors_ = *ColMapColors_;
  
      //ReducedMatrixRowMap_ = RowMapColors_.GenerateMap(0);
      //ReducedMatrixColMap_ = ColMapColors_.GenerateMap(0);
  
      auto RowMapColors_Data = RowMapColors_->getLocalViewHost(Tpetra::Access::ReadWrite);
      auto ColMapColors_Data = ColMapColors_->getLocalViewHost(Tpetra::Access::ReadWrite);
      auto size1 = RowMapColors_Data.size();
      for (int j=0; j<int(size1); j++) {
        if (RowMapColors_Data(j,0) != 0) {
          PRINT_NB(outputRank_, "ConstructReducedProblem() A    rowMapColors["<<j<<"] = " << RowMapColors_Data(j,0));
        }
        if (ColMapColors_Data(j,0) != 0) {
          PRINT_NB(outputRank_, "ConstructReducedProblem() A    colMapColors["<<j<<"] = " << ColMapColors_Data(j,0));
        }
      }

      ReducedMatrixRowMap_ = GenerateReducedMap(FullMatrixRowMap(), RowMapColors_, 0);
      printMap(ReducedMatrixRowMap_);
  
      ReducedMatrixColMap_ = GenerateReducedMap(FullMatrixColMap(), ColMapColors_, 0);
      printMap(ReducedMatrixColMap_);

      // Create domain and range map colorings by exporting map coloring of column and row maps
  
      auto importer = FullCrsMatrix_->getCrsGraph()->getImporter();
      if (importer != Teuchos::null) {
        Teuchos::RCP<vector_type_int> DomainMapColors = Teuchos::rcp(new vector_type_int(FullMatrixDomainMap()));
        DomainMapColors->doImport(*ColMapColors_, *importer, Tpetra::ABSMAX);
        OrigReducedMatrixDomainMap_ = GenerateReducedMap(FullMatrixDomainMap(), DomainMapColors, 0);
      } else
        OrigReducedMatrixDomainMap_ = ReducedMatrixColMap_;
      printMap(OrigReducedMatrixDomainMap_);
  
      auto exporter = FullCrsMatrix_->getCrsGraph()->getImporter();
      if (FullMatrixIsCrsMatrix_) {
        if (exporter != Teuchos::null) { // Non-trivial exporter
          Teuchos::RCP<vector_type_int>RangeMapColors = Teuchos::rcp(new vector_type_int(FullMatrixRangeMap()));
          RangeMapColors->doExport(*ColMapColors_, *exporter, Tpetra::ABSMAX);
          ReducedMatrixRangeMap_ = GenerateReducedMap(FullMatrixRangeMap(), RangeMapColors, 0);
        }
        else
          ReducedMatrixRangeMap_ = ReducedMatrixRowMap_;
      }
      else
        ReducedMatrixRangeMap_ = ReducedMatrixRowMap_;
      printMap(ReducedMatrixRangeMap_);


      PRINT(outputRank_, "ConstructReducedProblem() C    SymmetricElimination_ = " << SymmetricElimination_);

      // Check to see if the reduced system domain and range maps are the same.
      // If not, we need to remap entries of the LHS multivector so that they are distributed
      // conformally with the rows of the reduced matrix and the RHS multivector
      SymmetricElimination_ = ReducedMatrixRangeMap_->isSameAs(*OrigReducedMatrixDomainMap_);
      if (!SymmetricElimination_) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(SymmetricElimination_!=0,
          std::runtime_error, "Symmetric Elimination has not been tested or transitioned from Epetra!");
        //ConstructRedistributeExporter(OrigReducedMatrixDomainMap_, ReducedMatrixRangeMap_,
        //                              RedistributeDomainExporter_, ReducedMatrixDomainMap_);
      } else {
        ReducedMatrixDomainMap_ = OrigReducedMatrixDomainMap_;
        OrigReducedMatrixDomainMap_ = Teuchos::null;
        //RedistributeDomainExporter_ = 0;
      }

      // Create pointer to Full RHS, LHS
      Teuchos::RCP<multivector_type> FullRHS = FullProblem()->getRHS();
      Teuchos::RCP<multivector_type> FullLHS = FullProblem()->getLHS();
      int NumVectors = FullLHS->getNumVectors();

      PRINT(outputRank_, "ConstructReducedProblem() D    NumVectors = " << NumVectors);

      PRINT(outputRank_, "ConstructReducedProblem() E    ReducedMatrixDomainMap() = " << ReducedMatrixDomainMap());
      //PRINT(outputRank_, "ConstructReducedProblem() E    FullMatrixDomainMap() = " << FullMatrixDomainMap());
      PRINT(outputRank_, "ConstructReducedProblem() E    ReducedMatrixRowMap() = " << ReducedMatrixRowMap());
      //PRINT(outputRank_, "ConstructReducedProblem() E    FullRHS->Map() = " << FullRHS->getMap());

      // Create importers
      Full2ReducedLHSImporter_ = Teuchos::rcp(new import_type(ReducedMatrixDomainMap(), FullMatrixDomainMap()));
      printMap(ReducedMatrixRowMap());
      printMap(FullRHS->getMap());
      Full2ReducedRHSImporter_ = Teuchos::rcp(new import_type(FullRHS->getMap(), ReducedMatrixRowMap()));
      if (Full2ReducedRHSImporter_ != Teuchos::null) {
          // Create a FancyOStream for output
          Teuchos::RCP<Teuchos::FancyOStream> fancyOut = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

          // Print detailed information about the importer
          Full2ReducedRHSImporter_->describe(*fancyOut, Teuchos::VERB_HIGH);
      } else {
          std::cout << "Full2ReducedRHSImporter_ is null." << std::endl;
      }


      // Construct Reduced Matrix
      LocalOrdinal maxNumEntries = FullCrsMatrix_->getLocalMaxNumRowEntries();
      ReducedMatrix_ = Teuchos::rcp(new crs_matrix_type(ReducedMatrixRowMap(), ReducedMatrixColMap(), maxNumEntries));

      if (outputRank_) PrintTpetraCrsMatrix(ReducedMatrix_);

      // Create storage for temporary X values due to explicit elimination of rows
      tempExportX_ = Teuchos::rcp(new multivector_type(FullMatrixColMap(), NumVectors));

      size_t NumEntries = 0;
      Teuchos::ArrayView<const Scalar> Values;
      Teuchos::Array<GlobalOrdinal> Indices;
      LocalOrdinal localNumRows = FullMatrix()->getLocalNumRows();
      LocalOrdinal ColSingletonCounter = 0;

      for (LocalOrdinal i=0; i<localNumRows; i++) {
        GlobalOrdinal curGRID = FullMatrixRowMap()->getGlobalElement(i);
        if (ReducedMatrixRowMap()->isNodeGlobalElement(curGRID)) { // Check if this row should go into reduced matrix
          GetRowGCIDs(i, NumEntries, Values, Indices); // Get current row (Indices are global)
          

          //// Debug: Print row and column indices being inserted
          //std::cout << "Inserting into row " << curGRID << " with column indices: ";
          //for (typename Teuchos::Array<GlobalOrdinal>::size_type j = 0; j < Indices.size(); ++j) {
          //    std::cout << Indices[j] << " ";
          //}
          //std::cout << std::endl;

          //// Debug: Print the column map on this process
          //Teuchos::ArrayView<const GlobalOrdinal> columnMapIndices = ReducedMatrixColMap()->getLocalElementList();
          //std::cout << "Column map on this process: ";
          //for (typename Teuchos::Array<GlobalOrdinal>::size_type j = 0; j < columnMapIndices.size(); ++j) {
          //    std::cout << columnMapIndices[j] << " ";
          //}
          //std::cout << std::endl;

          //// Debug: Check for invalid column indices
          //for (typename Teuchos::Array<GlobalOrdinal>::size_type j = 0; j < Indices.size(); ++j) {
          //    if (!ReducedMatrixColMap()->isNodeGlobalElement(Indices[j])) {
          //        std::cout << "Error: Column index " << Indices[j] << " is not in the column map on this process!" << std::endl;
          //    }
          //}


          //ReducedMatrix()->insertGlobalValues(curGRID, Indices, Values);

          // Filter indices and values
          Teuchos::Array<GlobalOrdinal> filteredIndices;
          Teuchos::Array<Scalar> filteredValues;

          for (typename Teuchos::Array<GlobalOrdinal>::size_type j = 0; j < Indices.size(); ++j) {
              if (ReducedMatrixColMap()->isNodeGlobalElement(Indices[j])) {
                  filteredIndices.push_back(Indices[j]);
                  filteredValues.push_back(Values[j]);
              //} else {
              //    std::cout << "Excluding column index " << Indices[j] << " as it is not in the column map on this process!" << std::endl;
              }
          }

          // Insert filtered values into the matrix
          if (!filteredIndices.empty()) {
              ReducedMatrix()->insertGlobalValues(curGRID, filteredIndices(), filteredValues());
          //} else {
          //    std::cout << "Skipping insertion for row " << curGRID << " as no valid column indices remain." << std::endl;
          }

          // Insert into reduce matrix
//        // Positive errors will occur because we are submitting col entries that are not part of
//        // reduced system.  However, because we specified a column map to the ReducedMatrix constructor
//        // these extra column entries will be ignored and we will be politely reminded by a positive
//        // error code
//        if (ierr<0) EPETRA_CHK_ERR(ierr);
        }
        else {
          Teuchos::ArrayView<const LocalOrdinal> localIndices;
          PRINT(outputRank_, "ConstructReducedProblem() F    i, NumEntries = " << i << ", " << NumEntries);
          GetRow(i, NumEntries, Values, localIndices); // Get current row
          if (NumEntries==1) {
            Scalar pivot = Values[0];
            TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(pivot == 0.0, std::runtime_error,
              "Encountered zero row, unable to continue.");  // Should improve this comparison to zero.
            LocalOrdinal indX = localIndices[0];
            for (LocalOrdinal j=0; j<NumVectors; j++) {
              //(*tempExportX_)(j, indX) = (*FullRHS)(j, i) / pivot;
              //tempExportX_->replaceLocalValue(indX, j, FullRHS->getData(i, j) / pivot);
              auto rhsData = FullRHS->getData(j); // Get the underlying data for vector j
              auto exportData = tempExportX_->getDataNonConst(j); // Get the underlying data for vector j (non-const)
              exportData[indX] = rhsData[i] / pivot;
            }
          } else { // Singleton column
            LocalOrdinal targetCol = ColSingletonColLIDs_[ColSingletonCounter];
            for (size_t j=0; j<NumEntries; j++) {
              if (localIndices[j]==targetCol) {
                Scalar pivot = Values[j];
                TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(pivot == 0.0, std::runtime_error,
                  "Encountered zero column, unable to continue.");  // Should improve this comparison to zero.
                ColSingletonPivotLIDs_[ColSingletonCounter] = j;
                ColSingletonPivots_[ColSingletonCounter] = pivot;
                ColSingletonCounter++;
                break;
              }
            }
          }
        }
      }

      PRINT(outputRank_, "ConstructReducedProblem() G");

      // Now convert to local indexing.  We have constructed things so that the domain and range of the
      // matrix will have the same map.  If the reduced matrix domain and range maps were not the same, the
      // differences were addressed in the ConstructRedistributeExporter() method
      ReducedMatrix()->fillComplete(ReducedMatrixDomainMap(), ReducedMatrixRangeMap());

      // 1) The vector ColProfiles has column nonzero counts for each processor's contribution
      // Construct Reduced LHS (Puts any initial guess values into reduced system)
  
      // Create the ReducedLHS_ MultiVector
      ReducedLHS_ = Teuchos::rcp(new multivector_type(ReducedMatrixDomainMap(), NumVectors));
      // Import data from FullLHS to ReducedLHS_
      ReducedLHS_->doImport(*FullLHS, *Full2ReducedLHSImporter_, Tpetra::INSERT);
      if (outputRank_) printMultiVector("ConstructReducedProblem() G    ReducedLHS_ = ", ReducedLHS_);
      FullLHS->putScalar(0.0);


      
  
      // Construct Reduced RHS
  
      // First compute influence of already-known values of X on RHS
//    tempX_ = new Epetra_MultiVector(FullMatrixDomainMap(), NumVectors);
//    tempB_ = new Epetra_MultiVector(FullRHS->Map(), NumVectors);
      tempX_ = Teuchos::rcp(new multivector_type(FullMatrixDomainMap(), NumVectors));
      tempB_ = Teuchos::rcp(new multivector_type(FullRHS->getMap(), NumVectors));

      //Inject known X values into tempX for purpose of computing tempB = FullMatrix*tempX
      // Also inject into full X since we already know the solution

      //auto importer = FullCrsMatrix_->getCrsGraph()->getImporter();
      if (importer != Teuchos::null) {
//      EPETRA_CHK_ERR(tempX_->Export(*tempExportX_, *FullMatrix()->RowMatrixImporter(), Add));
//      EPETRA_CHK_ERR(FullLHS->Export(*tempExportX_, *FullMatrix()->RowMatrixImporter(), Add));
        tempX_->doExport(*tempExportX_, *importer, Tpetra::ADD);
        FullLHS->doExport(*tempExportX_, *importer, Tpetra::ADD);

      }
      else {
//      tempX_->Update(1.0, *tempExportX_, 0.0);
//      FullLHS->Update(1.0, *tempExportX_, 0.0);
        tempX_->update(1.0, *tempExportX_, 0.0); // tempX_ = 1.0 * tempExportX_ + 0.0 * tempX_
        FullLHS->update(1.0, *tempExportX_, 0.0); // FullLHS = 1.0 * tempExportX_ + 0.0 * FullLHS

        PRINT(outputRank_, "ConstructReducedProblem() L");
        if (outputRank_) printMultiVector("ConstructReducedProblem() L    tempX_ = ", tempX_);
        if (outputRank_) printMultiVector("ConstructReducedProblem() L    FullLHS = ", FullLHS);

      }

//    EPETRA_CHK_ERR(FullMatrix()->Multiply(false, *tempX_, *tempB_));
//
//    EPETRA_CHK_ERR(tempB_->Update(1.0, *FullRHS, -1.0)); // tempB now has influence of already-known X values

      FullMatrix()->apply(*tempX_, *tempB_);
      tempB_->update(1.0, *FullRHS, -1.0);

      PRINT(outputRank_, "ConstructReducedProblem() M");
      if (outputRank_) printMultiVector("ConstructReducedProblem() M    tempB_ = ", tempB_);

//    ReducedRHS_ = new Epetra_MultiVector(*ReducedMatrixRowMap(), FullRHS->NumVectors());
//    EPETRA_CHK_ERR(ReducedRHS_->Import(*tempB_, *Full2ReducedRHSImporter_, Insert));
      ReducedRHS_ = Teuchos::rcp(new multivector_type(ReducedMatrixRowMap(), FullRHS->getNumVectors()));
      if (outputRank_) printMultiVector("ConstructReducedProblem() S    ReducedRHS_ = ", ReducedRHS_);
      ReducedRHS_->doImport(*tempB_, *Full2ReducedRHSImporter_, Tpetra::INSERT);

      if (Full2ReducedRHSImporter_ != Teuchos::null) {
          // Create a FancyOStream for output
          Teuchos::RCP<Teuchos::FancyOStream> fancyOut = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

          // Print detailed information about the importer
          Full2ReducedRHSImporter_->describe(*fancyOut, Teuchos::VERB_HIGH);
      } else {
          std::cout << "Full2ReducedRHSImporter_ is null." << std::endl;
      }

      PRINT(outputRank_, "ConstructReducedProblem() S");
      if (outputRank_) printMultiVector("ConstructReducedProblem() S    ReducedLHS_ = ", ReducedLHS_);
      if (outputRank_) printMultiVector("ConstructReducedProblem() S    ReducedRHS_ = ", ReducedRHS_);

//    // Finally construct Reduced Linear Problem
      ReducedProblem_ = Teuchos::rcp(new Tpetra::LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
        ReducedMatrix_, ReducedLHS_, ReducedRHS_));
   }
   else {
  
     // There are no singletons, so don't bother building a reduced problem.
     ReducedProblem_ = Problem;
     ReducedMatrix_ = Teuchos::rcp(dynamic_cast<crs_matrix_type *>(Problem->getMatrix().getRawPtr()), false);

   }

    PRINT(outputRank_, "ConstructReducedProblem() Z");
    if (outputRank_) PrintTpetraCrsMatrix(ReducedMatrix_);

 
    double fn = (double) FullMatrix()->getGlobalNumRows();
    double fnnz = (double) FullMatrix()->getGlobalNumEntries();
    PRINT(outputRank_, "ConstructReducedProblem() Z    fn = " << fn);
    PRINT(outputRank_, "ConstructReducedProblem() Z    fnnz = " << fnnz);

    double rn = (double) ReducedMatrix()->getGlobalNumRows();
    double rnnz = (double) ReducedMatrix()->getGlobalNumEntries();
    PRINT(outputRank_, "ConstructReducedProblem() Z    rn = " << rn);
    PRINT(outputRank_, "ConstructReducedProblem() Z    rnnz = " << rnnz);
 
    RatioOfDimensions_ = rn/fn;
    RatioOfNonzeros_ = rnnz/fnnz;
    HaveReducedProblem_ = true;
 
    return;
  }


//==============================================================================
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  UpdateReducedProblem(const Teuchos::RCP<linear_problem_type> & Problem) {
  
    PRINT(outputRank_, "UpdateReducedProblem() A");
    const char tfecfFuncName[] = "UpdateReducedProblem: ";
    
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(!HaveReducedProblem_, std::runtime_error,
      "Must have a reduced problem.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(Problem == Teuchos::null,
      std::runtime_error, "Problem is Teuchos::null.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(FullMatrix_ == Teuchos::null,
      std::runtime_error, "Need a RowMatrix.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(Problem->getRHS() == Teuchos::null,
      std::runtime_error, "Need a RHS.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(Problem->getLHS() == Teuchos::null,
      std::runtime_error, "Need a LHS.");
 
  
    PRINT(outputRank_, "UpdateReducedProblem() A    SingletonsDetected() = " << SingletonsDetected());
    if ( SingletonsDetected() ) {

      // Create pointer to Full RHS, LHS
      Teuchos::RCP<multivector_type> FullRHS = FullProblem()->getRHS();
      Teuchos::RCP<multivector_type> FullLHS = FullProblem()->getLHS();
      int NumVectors = FullLHS->getNumVectors();

      PRINT(outputRank_, "UpdateReducedProblem() B    NumVectors = " << NumVectors);

      tempExportX_->putScalar(0.0);

      size_t NumEntries = 0;
      Teuchos::ArrayView<const Scalar> Values;
      Teuchos::Array<GlobalOrdinal> Indices;
      LocalOrdinal localNumRows = FullMatrix()->getLocalNumRows();
      LocalOrdinal ColSingletonCounter = 0;

      for (LocalOrdinal i=0; i<localNumRows; i++) {
        GlobalOrdinal curGRID = FullMatrixRowMap()->getGlobalElement(i);
        if (ReducedMatrixRowMap()->isNodeGlobalElement(curGRID)) { // Check if this row should go into reduced matrix
          GetRowGCIDs(i, NumEntries, Values, Indices);

          // Filter indices and values
          Teuchos::Array<GlobalOrdinal> filteredIndices;
          Teuchos::Array<Scalar> filteredValues;

          for (typename Teuchos::Array<GlobalOrdinal>::size_type j = 0; j < Indices.size(); ++j) {
              if (ReducedMatrixColMap()->isNodeGlobalElement(Indices[j])) {
                  filteredIndices.push_back(Indices[j]);
                  filteredValues.push_back(Values[j]);
              //} else {
              //    std::cout << "Excluding column index " << Indices[j] << " as it is not in the column map on this process!" << std::endl;
              }
          }

          // Insert filtered values into the matrix
          if (!filteredIndices.empty()) {
              ReducedMatrix()->replaceGlobalValues(curGRID, filteredIndices(), filteredValues());
          //} else {
          //    std::cout << "Skipping insertion for row " << curGRID << " as no valid column indices remain." << std::endl;
          }


        }
        // Otherwise if singleton row we explicitly eliminate this row and solve for corresponding X value
        else {

          Teuchos::ArrayView<const LocalOrdinal> localIndices;
          GetRow(i, NumEntries, Values, localIndices); // Get current row
          if (NumEntries==1) {
            Scalar pivot = Values[0];
            TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(pivot == 0.0, std::runtime_error,
              "Encountered zero row, unable to continue.");  // Should improve this comparison to zero.
            LocalOrdinal indX = localIndices[0];
            for (LocalOrdinal j=0; j<NumVectors; j++) {
              auto rhsData = FullRHS->getData(j); // Get the underlying data for vector j
              auto exportData = tempExportX_->getDataNonConst(j); // Get the underlying data for vector j (non-const)
              exportData[indX] = rhsData[i] / pivot;
            }
          } else { // Singleton column
            LocalOrdinal j = ColSingletonColLIDs_[ColSingletonCounter];
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

      ReducedLHS_->putScalar(0.0); // zero out Reduced LHS
      ReducedRHS_->doImport(*FullLHS, *Full2ReducedRHSImporter_, Tpetra::INSERT);
      
      FullLHS->putScalar(0.0); // zero out Full LHS since we will inject values as we get them

      // Construct Reduced RHS

      // Zero out temp space
      tempX_->putScalar(0.0);
      tempB_->putScalar(0.0);
    
  //  //Inject known X values into tempX for purpose of computing tempB = FullMatrix*tempX
  //  // Also inject into full X since we already know the solution

      auto importer = FullCrsMatrix_->getCrsGraph()->getImporter();
      if (importer != Teuchos::null) {
        tempX_->doExport(*tempExportX_, *importer, Tpetra::ADD);
        FullLHS->doExport(*tempExportX_, *importer, Tpetra::ADD);
      }
      else {
        tempX_->update(1.0, *tempExportX_, 0.0); // tempX_ = 1.0 * tempExportX_ + 0.0 * tempX_
        FullLHS->update(1.0, *tempExportX_, 0.0); // FullLHS = 1.0 * tempExportX_ + 0.0 * FullLHS
      }

      FullMatrix()->apply(*tempX_, *tempB_);
      tempB_->update(1.0, *FullRHS, -1.0);

      ReducedRHS_->putScalar(0.0);
      ReducedRHS_->doImport(*tempB_, *Full2ReducedRHSImporter_, Tpetra::INSERT);

      if (outputRank_) printMultiVector("UpdateReducedProblem() Z    FullLHS = ", FullLHS);
      if (outputRank_) printMultiVector("UpdateReducedProblem() Z    FullRHS = ", FullRHS);
    }
    else {
    
      // There are no singletons, so don't bother building a reduced problem.
      ReducedProblem_ = Problem;
      ReducedMatrix_ = Teuchos::rcp(dynamic_cast<crs_matrix_type *>(Problem->getMatrix().getRawPtr()), false);
    }

    if (outputRank_) PrintTpetraRowMatrix(FullMatrix());
    if (outputRank_) printMultiVector("UpdateReducedProblem() Z    tempX_ = ", tempX_);
    if (outputRank_) printMultiVector("UpdateReducedProblem() Z    tempB_ = ", tempB_);
    if (outputRank_) printMultiVector("UpdateReducedProblem() Z    tempExportX_ = ", tempExportX_);
    if (outputRank_) PrintTpetraRowMatrix(ReducedMatrix());
    if (outputRank_) printMultiVector("UpdateReducedProblem() Z    ReducedLHS_ = ", ReducedLHS_);
    if (outputRank_) printMultiVector("UpdateReducedProblem() Z    ReducedRHS_ = ", ReducedRHS_);
  
    return;
  }

////==============================================================================
//template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
//void CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
//ConstructRedistributeExporter(Epetra_Map * SourceMap, Epetra_Map * TargetMap,
//                              Epetra_Export * & RedistributeExporter,
//                              Epetra_Map * & RedistributeMap) {
//
//  int_type IndexBase = (int_type) SourceMap->IndexBase64();
//  if (IndexBase!=(int_type) TargetMap->IndexBase64()) EPETRA_CHK_ERR(-1);
//
//  const Epetra_Comm & Comm = TargetMap->Comm();
//
//  int TargetNumMyElements = TargetMap->NumMyElements();
//  int SourceNumMyElements = SourceMap->NumMyElements();
//
//  // ContiguousTargetMap has same number of elements per PE as TargetMap, but uses contigious indexing
//  Epetra_Map ContiguousTargetMap((int_type) -1, TargetNumMyElements, IndexBase,Comm);
//
//  // Same for ContiguousSourceMap
//  Epetra_Map ContiguousSourceMap((int_type) -1, SourceNumMyElements, IndexBase, Comm);
//
//  assert(ContiguousSourceMap.NumGlobalElements64()==ContiguousTargetMap.NumGlobalElements64());
//
//  // Now create a vector that contains the global indices of the Source Epetra_MultiVector
//  int_type* SourceMapMyGlobalElements = 0;
//  SourceMap->MyGlobalElementsPtr(SourceMapMyGlobalElements);
//  typename Epetra_GIDTypeVector<int_type>::impl SourceIndices(View, ContiguousSourceMap, SourceMapMyGlobalElements);
//
//  // Create an exporter to send the SourceMap global IDs to the target distribution
//  Epetra_Export Exporter(ContiguousSourceMap, ContiguousTargetMap);
//
//  // Create a vector to catch the global IDs in the target distribution
//  typename Epetra_GIDTypeVector<int_type>::impl TargetIndices(ContiguousTargetMap);
//  TargetIndices.Export(SourceIndices, Exporter, Insert);
//
//  // Create a new map that describes how the Source MultiVector should be laid out so that it has
//  // the same number of elements on each processor as the TargetMap
//  RedistributeMap = new Epetra_Map((int_type) -1, TargetNumMyElements, TargetIndices.Values(), IndexBase, Comm);
//
//  // This exporter will finally redistribute the Source MultiVector to the same layout as the TargetMap
//  RedistributeExporter = new Epetra_Export(*SourceMap, *RedistributeMap);
//  return(0);
//}
//
////==============================================================================
//int CrsSingletonFilter_LinearProblem::ComputeFullSolution() {
//
//  if ( SingletonsDetected() ) {
//    int jj, k;
//
//    Epetra_MultiVector * FullLHS = FullProblem()->GetLHS();
//    Epetra_MultiVector * FullRHS = FullProblem()->GetRHS();
//
//    tempX_->PutScalar(0.0); tempExportX_->PutScalar(0.0);
//    // Inject values that the user computed for the reduced problem into the full solution vector
//    EPETRA_CHK_ERR(tempX_->Export(*ReducedLHS_, *Full2ReducedLHSImporter_, Add));
//
//    FullLHS->Update(1.0, *tempX_, 1.0);
//
//    // Next we will use our full solution vector which is populated with pre-filter solution
//    // values and reduced system solution values to compute the sum of the row contributions
//    // that must be subtracted to get the post-filter solution values
//
//    EPETRA_CHK_ERR(FullMatrix()->Multiply(false, *FullLHS, *tempB_));
//
//    // Finally we loop through the local rows that were associated with column singletons and compute the
//    // solution for these equations.
//
//    int NumVectors = tempB_->NumVectors();
//    for (k=0; k<localNumSingletonCols_; k++) {
//      int i = ColSingletonRowLIDs_[k];
//      int j = ColSingletonColLIDs_[k];
//      double pivot = ColSingletonPivots_[k];
//      for (jj=0; jj<NumVectors; jj++)
//        (*tempExportX_)[jj][j]= ((*FullRHS)[jj][i] - (*tempB_)[jj][i])/pivot;
//    }
//
//    // Finally, insert values from post-solve step and we are done!!!!
//
//    if (FullMatrix()->RowMatrixImporter()!=0) {
//      EPETRA_CHK_ERR(tempX_->Export(*tempExportX_, *FullMatrix()->RowMatrixImporter(), Add));
//    }
//    else {
//      tempX_->Update(1.0, *tempExportX_, 0.0);
//    }
//
//    FullLHS->Update(1.0, *tempX_, 1.0);
//  }
//
//  return(0);
//}
//==============================================================================

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
CrsSingletonFilter_LinearProblem::InitFullMatrixAccess() {

  PRINT(outputRank_, "InitFullMatrixAccess() A");
  PRINT(outputRank_, "InitFullMatrixAccess() A    FullMatrix() = " << FullMatrix());
  if (outputRank_) PrintTpetraRowMatrix(FullMatrix());
  Tpetra::getDefaultComm()->barrier();

  //MaxNumMyEntries_ = FullMatrix()->MaxNumEntries();
  localMaxNumRowEntries_ = FullMatrix()->getLocalMaxNumRowEntries();
  PRINT(outputRank_, "InitFullMatrixAccess() B    localMaxNumRowEntries_ = " << localMaxNumRowEntries_);

  // Cast to CrsMatrix, if possible.  Can save some work.
  FullCrsMatrix_ = Teuchos::rcp_dynamic_cast<crs_matrix_type>(FullMatrix());
  PRINT(outputRank_, "InitFullMatrixAccess() C    FullCrsMatrix_ = " << FullCrsMatrix_);
  if (outputRank_) PrintTpetraCrsMatrix(FullCrsMatrix_);
  Tpetra::getDefaultComm()->barrier();
  FullMatrixIsCrsMatrix_ = (FullCrsMatrix_ != Teuchos::null); // Pointer is non-null if cast worked
  PRINT(outputRank_, "InitFullMatrixAccess() D    FullMatrixIsCrsMatrix_ = " << FullMatrixIsCrsMatrix_);
  Indices_ = Teuchos::arcp(new global_ordinal_type[localMaxNumRowEntries_],0,localMaxNumRowEntries_,true);
  PRINT(outputRank_, "InitFullMatrixAccess() E    Indices_ = " << Indices_);
  Values_  = Teuchos::arcp(new scalar_type[localMaxNumRowEntries_],0,localMaxNumRowEntries_,true);
  PRINT(outputRank_, "InitFullMatrixAccess() F    Values_ = " << Values_);

  return;
}
//==============================================================================

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
GetRow(local_ordinal_type localRow, size_t & NumIndices,
       Teuchos::Array<local_ordinal_type> & localIndices) {

  PRINT(outputRank_, "GetRow() A      NumIndices = " << NumIndices);

  // Must get the values, but we ignore them.
  size_t maxNumEntries = FullCrsMatrix_->getLocalMaxNumRowEntries();
  nonconst_local_inds_host_view_type Indices("indices", maxNumEntries);
  nonconst_values_host_view_type Values("values", maxNumEntries);
  PRINT(outputRank_, "GetRow() B      NumIndices = " << NumIndices);
  PRINT(outputRank_, "GetRow() B      FullMatrixIsCrsMatrix_ = " << FullMatrixIsCrsMatrix_);
  if (FullMatrixIsCrsMatrix_) { // View of current row
    // EPETRA_CHK_ERR(FullCrsMatrix()->Graph().ExtractMyRowView(localRow, NumIndices, localIndices));
    //FullCrsMatrix_->getLocalRowCopy(localRow, localIndices, Values, NumIndices);
    PRINT(outputRank_, "GetRow() C      NumIndices = " << NumIndices);
    PRINT(outputRank_, "GetRow() C      localRow = " << localRow);
    FullCrsMatrix_->getLocalRowCopy(localRow, Indices, Values, NumIndices);
  }
  else { // Copy of current row (we must get the values, but we ignore them)
    // EPETRA_CHK_ERR(FullMatrix()->ExtractMyRowCopy(localRow, localMaxNumRowEntries_, NumIndices, Values_.Values(), Indices_int_));
    //FullMatrix_->getLocalRowCopy(localRow, localIndices, Values, NumIndices);
    PRINT(outputRank_, "GetRow() D      NumIndices = " << NumIndices);
    PRINT(outputRank_, "GetRow() D      localRow = " << localRow);
    FullMatrix_->getLocalRowCopy(localRow, Indices, Values, NumIndices);
  }
  // If LocalRow does not belong to the calling process, then the method sets 
  // NumIndices to Teuchos::OrdinalTraits<size_t>::invalid(), and does not
  // modify localIndices or Values. 

  // Copy the indices to the output array
  localIndices.resize(NumIndices);
  for (size_t i = 0; i < NumIndices; ++i) {
     localIndices[i] = Indices(i);
  }
  PRINT(outputRank_, "GetRow() E      NumIndices = " << NumIndices);
  return;
}
////==============================================================================

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
GetRow(LocalOrdinal Row, size_t &NumIndices, Teuchos::ArrayView<const Scalar> &Values, 
       Teuchos::ArrayView<const LocalOrdinal> &Indices) 
{
    typename Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_inds_host_view_type localIndices;
    typename Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::values_host_view_type rowValues;

    PRINT(outputRank_, "GetRow(,,,) A   Row, NumIndices = " << Row << ", " << NumIndices);

    if (FullMatrixIsCrsMatrix_) { // View of current row
        FullCrsMatrix_->getLocalRowView(Row, localIndices, rowValues);

        NumIndices = localIndices.size();
        Values = Teuchos::ArrayView<const Scalar>(rowValues.data(), rowValues.size());
        Indices = Teuchos::ArrayView<const LocalOrdinal>(localIndices.data(), localIndices.size());

        PRINT(outputRank_, "GetRow(,,,) B   Row, NumIndices = " << Row << ", " << NumIndices);
        PRINT(outputRank_, "GetRow(,,,) B   Values = " << Values);
        PRINT(outputRank_, "GetRow(,,,) B   Indices = " << Indices);

    } else { // Copy of current row
        typename Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::nonconst_local_inds_host_view_type localIndicesCopy("localIndicesCopy", FullMatrix()->getLocalMaxNumRowEntries());
        typename Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::nonconst_values_host_view_type rowValuesCopy("rowValuesCopy", FullMatrix()->getLocalMaxNumRowEntries());

        FullMatrix()->getLocalRowCopy(Row, localIndicesCopy, rowValuesCopy, NumIndices);

        Values = Teuchos::ArrayView<const Scalar>(rowValuesCopy.data(), NumIndices);
        Indices = Teuchos::ArrayView<const LocalOrdinal>(localIndicesCopy.data(), NumIndices);

        PRINT(outputRank_, "GetRow(,,,) C   Row, NumIndices = " << Row << ", " << NumIndices);
        PRINT(outputRank_, "GetRow(,,,) C   Values = " << Values);
        PRINT(outputRank_, "GetRow(,,,) C   Indices = " << Indices);
    }
    return;
}

//==============================================================================

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
GetRowGCIDs(
    LocalOrdinal localRow, 
    size_t & NumIndices, 
    Teuchos::ArrayView<const Scalar> & Values, 
    Teuchos::Array<GlobalOrdinal> & GlobalIndices) 
{

    PRINT(outputRank_, "GetRowGCIDs() A      Row = " << localRow);

    // Extract the row data (local indices and values)
    typename Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_inds_host_view_type LocalIndices;
    typename Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::values_host_view_type RowValues;

    FullMatrix()->getLocalRowView(localRow, LocalIndices, RowValues);

    // Convert local indices to global indices
    NumIndices = LocalIndices.size();
    PRINT(outputRank_, "GetRowGCIDs() A      NumIndices = " << NumIndices);
    GlobalIndices.resize(NumIndices); // Resize the array to hold global indices
    for (size_t j = 0; j < NumIndices; ++j) {
        GlobalIndices[j] = FullMatrixColMap()->getGlobalElement(LocalIndices[j]);
    }
    // Copy values into the provided ArrayView
    Values = Teuchos::ArrayView<const Scalar>(RowValues.data(), RowValues.size());

    // Print the Indices array
    std::cout << "petra : GetRowGCIDs() A      Row " << localRow << " Indices: ";
    for (LocalOrdinal j = 0; j < GlobalIndices.size(); ++j) {
        std::cout << GlobalIndices[j] << " ";
    }
    std::cout << std::endl;
    // Print the Values array
    std::cout << "petra : GetRowGCIDs() A      Row " << localRow << " Values: ";
    for (LocalOrdinal j = 0; j < Values.size(); ++j) {
        std::cout << Values[j] << " ";
    }
    std::cout << std::endl;

    return;
}
////template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
////void CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
////GetRowGCIDs(local_ordinal_type localRow, size_t & NumIndices, 
////nonconst_values_host_view_type & Values, Teuchos::Array<local_ordinal_type> & GlobalIndices) {
////
////  
////  PRINT(outputRank_, "GetRowGCIDs() A      Row = " << localRow);
////  PRINT(outputRank_, "GetRowGCIDs() A      NumIndices = " << NumIndices);
////
////  //EPETRA_CHK_ERR(FullMatrix()->ExtractMyRowCopy(Row, localMaxNumRowEntries_, NumIndices,
////  //                                              Values_.Values(), Indices_int_));
////  FullMatrix_->getLocalRowCopy(localRow, GlobalIndices, Values, NumIndices);
////
////  for (int j=0; j<NumIndices; j++) Indices_int_[j] = FullMatrixColMap().GID(Indices_int_[j]);
////  Values = Values_.Values();
////  GlobalIndices = Indices_int_;
////
////  return;
////}

//#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
//int CrsSingletonFilter_LinearProblem::GetRowGCIDs(int Row, int & NumIndices,
//                                           double * & Values, int * & GlobalIndices) {
//
//    EPETRA_CHK_ERR(FullMatrix()->ExtractMyRowCopy(Row, localMaxNumRowEntries_, NumIndices,
//                                                  Values_.Values(), Indices_int_));
//    for (int j=0; j<NumIndices; j++) Indices_int_[j] = FullMatrixColMap().GID(Indices_int_[j]);
//    Values = Values_.Values();
//    GlobalIndices = Indices_int_;
//  return(0);
//}
//#endif
//
//#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
//int CrsSingletonFilter_LinearProblem::GetRowGCIDs(int Row, int & NumIndices,
//                                           double * & Values, long long * & GlobalIndices) {
//    EPETRA_CHK_ERR(FullMatrix()->ExtractMyRowCopy(Row, localMaxNumRowEntries_, NumIndices,
//                                                  Values_.Values(), Indices_int_));
//    for (int j=0; j<NumIndices; j++) Indices_LL_[j] = FullMatrixColMap().GID64(Indices_int_[j]);
//    Values = Values_.Values();
//    GlobalIndices = Indices_LL_;
//  return(0);
//}
//#endif

//==============================================================================

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CrsSingletonFilter_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
CreatePostSolveArrays(vector_type_LO localRowIDofSingletonCol,
                      vector_type_LO ColProfiles,
                      vector_type_LO NewColProfiles,
                      vector_type_LO ColHasRowWithSingleton) {

   const char tfecfFuncName[] = "CreatePostSolveArrays: ";

   if (localNumSingletonCols_ == 0) return; // Nothing to do

  // We will need these arrays for the post-solve phase

  // ColSingletonRowLIDs_ = new int[localNumSingletonCols_];
  // ColSingletonColLIDs_ = new int[localNumSingletonCols_];
  // ColSingletonPivotLIDs_ = new int[localNumSingletonCols_];
  // ColSingletonPivots_ = new double[localNumSingletonCols_];

  ColSingletonRowLIDs_ = Teuchos::arcp(new local_ordinal_type[localNumSingletonCols_],0,localNumSingletonCols_,true);
  ColSingletonColLIDs_ = Teuchos::arcp(new local_ordinal_type[localNumSingletonCols_],0,localNumSingletonCols_,true);
  ColSingletonPivotLIDs_ = Teuchos::arcp(new local_ordinal_type[localNumSingletonCols_],0,localNumSingletonCols_,true);
  ColSingletonPivots_ = Teuchos::arcp(new scalar_type[localNumSingletonCols_],0,localNumSingletonCols_,true);

  // Register singleton columns (that were not already counted as singleton rows)
  // Check to see if any columns disappeared because all associated rows were eliminated
  int NumMyColSingletonstmp = 0;
  local_ordinal_type localNumCols = FullMatrix_->getLocalNumCols();
  PRINT(outputRank_, "CreatePostSolveArrays() A    NumMyCols = " << localNumCols);
  auto localRowIDofSingletonColData = localRowIDofSingletonCol.getLocalViewHost(Tpetra::Access::ReadWrite);
  auto ColProfilesData = ColProfiles.getLocalViewHost(Tpetra::Access::ReadWrite);
  auto RowMapColors_Data = RowMapColors_->getLocalViewHost(Tpetra::Access::ReadWrite);
  auto NewColProfilesData = NewColProfiles.getLocalViewHost(Tpetra::Access::ReadWrite);
  auto ColHasRowWithSingletonData = ColHasRowWithSingleton.getLocalViewHost(Tpetra::Access::ReadWrite);
  auto ColMapColors_Data = ColMapColors_->getLocalViewHost(Tpetra::Access::ReadWrite);

  for (int j=0; j<localNumCols; j++) {
    int i = localRowIDofSingletonColData(j,0);
    PRINT(outputRank_, "CreatePostSolveArrays() B    j, i = " <<j<< ", " <<i);
    PRINT(outputRank_, "CreatePostSolveArrays() B    ColProfiles["<<j<<"] = " << ColProfilesData(j,0));
    PRINT(outputRank_, "CreatePostSolveArrays() B    rowMapColors["<<i<<"] = " << RowMapColors_Data(i,0));
    if ( ColProfilesData(j,0)==1 && RowMapColors_Data(i,0)!=1 ) {
      ColSingletonRowLIDs_[NumMyColSingletonstmp] = i;
      ColSingletonColLIDs_[NumMyColSingletonstmp] = j;
      NumMyColSingletonstmp++;
    }
    // Also check for columns that were eliminated implicitly by
    // having all associated row eliminated
    else if (NewColProfilesData(j,0)==0 && ColHasRowWithSingletonData(j,0)!=1 && RowMapColors_Data(i,0)==0) {
      ColMapColors_Data(j,0) = 1;
    }
  }

  for (int j=0; j<localNumSingletonCols_; j++) {
    PRINT(outputRank_, "CreatePostSolveArrays() C    ColSingletonRowLIDs_["<<j<<"] = " << ColSingletonRowLIDs_[j]); 
    PRINT(outputRank_, "CreatePostSolveArrays() C    ColSingletonColLIDs_["<<j<<"] = " << ColSingletonColLIDs_[j]);
  }

  auto size1 = RowMapColors_Data.size();
  for (int j=0; j<int(size1); j++) {
    if (ColMapColors_Data(j,0) != 0) {
      PRINT_NB(outputRank_, "CreatePostSolveArrays() C    colMapColors["<<j<<"] = " << ColMapColors_Data(j,0));
    }
    Tpetra::getDefaultComm()->barrier();
  }

  PRINT(outputRank_, "CreatePostSolveArrays() C    NumMyColSingletonstmp = " << NumMyColSingletonstmp);
  PRINT(outputRank_, "CreatePostSolveArrays() C    NumMyColSingletons_ = " << localNumSingletonCols_);

  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(NumMyColSingletonstmp != localNumSingletonCols_,
    std::runtime_error, "Sanity check.");
  
  Tpetra::sort2(ColSingletonRowLIDs_.begin(), ColSingletonRowLIDs_.end(), ColSingletonColLIDs_.begin());

  for (int j=0; j<localNumSingletonCols_; j++) {
    PRINT(outputRank_, "CreatePostSolveArrays() D    ColSingletonRowLIDs_["<<j<<"] = " << ColSingletonRowLIDs_[j]);
    PRINT(outputRank_, "CreatePostSolveArrays() D    ColSingletonColLIDs_["<<j<<"] = " << ColSingletonColLIDs_[j]);
  }

  return;
}

} //namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_CRSSINGLETONFILTER_INSTANT(SCALAR,LO,GO,NODE) \
  template class CrsSingletonFilter_LinearProblem< SCALAR , LO , GO , NODE >;


#endif //  TPETRA_CRSSINGLETONFILTER_LINEARPROBLEM_DEF_HPP
