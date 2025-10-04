
//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER

#include "Epetra_ConfigDefs.h"
#include "Epetra_Map.h"
#include "Epetra_Util.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_GIDTypeVector.h"
#include "Epetra_Comm.h"
#include "Epetra_MpiComm.h"                                //  REMOVE TEMP!!
#include "Epetra_LinearProblem.h"
#include "Epetra_MapColoring.h"
#include "EpetraExt_CrsSingletonFilter_LinearProblem.h"


void printEpetraMultiVector(
    const std::string& label,
    const Epetra_MultiVector* multiVector)
{
    // Print the label
    std::cout << "petra : " << label << std::endl;

    // Loop over columns (vectors) in the MultiVector
    for (int col = 0; col < multiVector->NumVectors(); ++col) {
        double* data;
        int stride;

        // Extract the view of the data for the current column
        multiVector->ExtractView(&data, &stride);

        std::cout << "Column " << col << ": ";

        // Loop over rows (entries) in the column
        for (int row = 0; row < multiVector->MyLength(); ++row) {
            std::cout << "(" << row << ", " << data[row + col * stride] << ") ";
        }
        std::cout << std::endl;
    }
}


void PrintEpetraRowMatrix(const Epetra_RowMatrix& matrix) {
    int numRows = matrix.NumMyRows();
    int maxNumEntries = matrix.MaxNumEntries();
    
    std::vector<int> indices(maxNumEntries);
    std::vector<double> values(maxNumEntries);

    for (int i = 0; i < numRows; ++i) {
        int numEntries;
        matrix.ExtractMyRowCopy(i, maxNumEntries, numEntries, values.data(), indices.data());

        std::cout << "Row " << i << ": ";
        for (int j = 0; j < numEntries; ++j) {
            std::cout << "(" << indices[j] << ", " << values[j] << ") ";
        }
        std::cout << std::endl;
    }
}

void PrintEpetraCrsMatrix(const Epetra_CrsMatrix& matrix) {
    int numRows = matrix.NumMyRows();
    int maxNumEntries = matrix.MaxNumEntries(); // Maximum number of entries across all rows

    std::vector<int> indices(maxNumEntries);
    std::vector<double> values(maxNumEntries);

    for (int i = 0; i < numRows; ++i) {
        int numEntries;
        matrix.ExtractMyRowCopy(i, maxNumEntries, numEntries, values.data(), indices.data());

        std::cout << "Row " << i << ": ";
        for (int j = 0; j < numEntries; ++j) {
            std::cout << "(" << indices[j] << ", " << values[j] << ") ";
        }
        std::cout << std::endl;
    }
}

void PrintEpetraCrsMatrixGlobal(const Epetra_CrsMatrix& matrix) {
    int numRows = matrix.NumGlobalRows(); // Get the total number of global rows

    for (int i = 0; i < numRows; ++i) {

        int maxNumEntries = matrix.NumMyEntries(i);
        //std::cout << "maxNumEntries = " << maxNumEntries << std::endl;

        std::vector<int> indices(maxNumEntries);
        std::vector<double> values(maxNumEntries);

        int numEntries = 0;
        int LenofIndices = matrix.NumMyEntries(i);
        //std::cout << "LenofIndices = " << LenofIndices << std::endl;
        
        std::cout << "Row " << i << ": ";
        if (LenofIndices != 0) {

          matrix.ExtractGlobalRowCopy(i, LenofIndices, numEntries, values.data(), indices.data());
          //int ierr = matrix.ExtractGlobalRowCopy(i, LenofIndices, numEntries, values.data(), indices.data());
          //if (ierr != 0) {
          //    std::cerr << "Error extracting row " << i << ": Epetra error code " << ierr << std::endl;
          //    continue;
          //}

          for (int j = 0; j < numEntries; ++j) {
              std::cout << "(" << indices[j] << ", " << values[j] << ") ";
          }
        }
        std::cout << std::endl;
    }
}


namespace EpetraExt {

  #define PRINT(outputRank_, message) \
    do { \
        if (outputRank_) { \
            std::cout << "petra : " << message << std::endl; \
        } \
        MPI_Barrier(MPI_COMM_WORLD); \
    } while (0)

  #define PRINT_NB(outputRank_, message) \
    do { \
        if (outputRank_) { \
            std::cout << "petra : " << message << std::endl; \
        } \
    } while (0)

//==============================================================================
LinearProblem_CrsSingletonFilter::
LinearProblem_CrsSingletonFilter( bool verbose )
: verbose_(verbose)
{
  InitializeDefaults();
}
//==============================================================================
LinearProblem_CrsSingletonFilter::
~LinearProblem_CrsSingletonFilter()
{
  if (ReducedLHS_!=0) delete ReducedLHS_;
  if (ReducedRHS_!=0) delete ReducedRHS_;
  if (ReducedMatrixDomainMap_!=ReducedMatrixColMap_) delete ReducedMatrixDomainMap_;
  if (OrigReducedMatrixDomainMap_!=ReducedMatrixColMap_ &&
      OrigReducedMatrixDomainMap_!=0) delete OrigReducedMatrixDomainMap_;
  if (ReducedMatrixRangeMap_!=ReducedMatrixRowMap_) delete ReducedMatrixRangeMap_;
  if (ReducedMatrixRowMap_!=0) delete ReducedMatrixRowMap_;
  if (ReducedMatrixColMap_!=0) delete ReducedMatrixColMap_;
  if (Full2ReducedRHSImporter_!=0) delete Full2ReducedRHSImporter_;
  if (Full2ReducedLHSImporter_!=0) delete Full2ReducedLHSImporter_;
  if (RedistributeDomainExporter_!=0) delete RedistributeDomainExporter_;
  if (RowMapColors_!=0) delete RowMapColors_;
  if (ColMapColors_!=0) delete ColMapColors_;

  if (ColSingletonRowLIDs_ != 0) delete [] ColSingletonRowLIDs_;
  if (ColSingletonColLIDs_ != 0) delete [] ColSingletonColLIDs_;
  if (ColSingletonPivotLIDs_ != 0) delete [] ColSingletonPivotLIDs_;
  if (ColSingletonPivots_ != 0) delete [] ColSingletonPivots_;
  if (tempExportX_ != 0) delete tempExportX_;
  if (Indices_int_ != 0) delete [] Indices_int_;
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if (Indices_LL_ != 0) delete [] Indices_LL_;
#endif
  if (tempX_ != 0) delete tempX_;
  if (tempB_ != 0) delete tempB_;

}

//==============================================================================
LinearProblem_CrsSingletonFilter::NewTypeRef
LinearProblem_CrsSingletonFilter::
operator()( LinearProblem_CrsSingletonFilter::OriginalTypeRef orig )
{
  int outRank = 0;
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  int myRank = Comm.MyPID();
  outputRank_ = outRank == myRank;

  analyze( orig );
  return construct();
}

//==============================================================================
bool
LinearProblem_CrsSingletonFilter::
analyze( LinearProblem_CrsSingletonFilter::OriginalTypeRef orig )
{
  origObj_ = &orig;

  FullMatrix_ = orig.GetMatrix();


#ifdef NDEBUG
  (void) Analyze( FullMatrix_ );
#else
  // assert() statements go away if NDEBUG is defined.  Don't declare
  // the 'flag' variable if it never gets used.
  int flag = Analyze( FullMatrix_ );
  assert( flag >= 0 );
#endif // NDEBUG

  if ( verbose_ && FullMatrix_->Comm().MyPID()==0 ) {
    std::cout << "\nAnalyzed Singleton Problem:\n";
    std::cout << "---------------------------\n";
  }
  if ( SingletonsDetected() ) {
    if ( verbose_ && FullMatrix_->Comm().MyPID()==0 ) {
      std::cout << "Singletons Detected!" << std::endl;;
      std::cout << "Num Singletons:      " << NumSingletons() << std::endl;
    }
  }
  else {
    if ( verbose_ && FullMatrix_->Comm().MyPID()==0 )
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
  if ( verbose_ && FullMatrix_->Comm().MyPID()==0 )
    std::cout << "---------------------------\n\n";

  return true;
}

//==============================================================================
LinearProblem_CrsSingletonFilter::NewTypeRef
LinearProblem_CrsSingletonFilter::
construct()
{
  if( !origObj_ ) abort();

#ifdef NDEBUG
  (void) ConstructReducedProblem( origObj_ );
#else
  // assert() statements go away if NDEBUG is defined.  Don't declare
  // the 'flag' variable if it never gets used.
  int flag = ConstructReducedProblem( origObj_ );
  assert( flag >= 0 );
#endif // NDEBUG

  newObj_ = ReducedProblem();

  if( verbose_ && SingletonsDetected() && FullMatrix_->Comm().MyPID()==0 )
  {
    std::cout << "\nConstructedSingleton Problem:\n";
    std::cout << "---------------------------\n";
    std::cout << "RatioOfDimensions:   " << RatioOfDimensions() << std::endl;
    std::cout << "RatioOfNonzeros:     " << RatioOfNonzeros() << std::endl;
    std::cout << "---------------------------\n\n";
  }

  return *newObj_;
}

//==============================================================================
bool LinearProblem_CrsSingletonFilter::fwd()
{
  int ierr = UpdateReducedProblem( FullProblem_ );
  if( ierr ) std::cout << "EDT_LinearProblem_CrsSingletonFilter::UpdateReducedProblem FAILED!\n";

  return (ierr==0);
}

//==============================================================================
bool LinearProblem_CrsSingletonFilter::rvs()
{
  int ierr = ComputeFullSolution();
  if( ierr ) std::cout << "EDT_LinearProblem_CrsSingletonFilter::ComputeFullSolution FAILED!\n";

  return (ierr==0);
}

//==============================================================================
void LinearProblem_CrsSingletonFilter::InitializeDefaults() {

// Initialize all attributes that have trivial default values

  FullProblem_ = 0;
  FullMatrix_ = 0;
  ReducedRHS_ = 0;
  ReducedLHS_ = 0;
  ReducedMatrixRowMap_ = 0;
  ReducedMatrixColMap_ = 0;
  ReducedMatrixDomainMap_ = 0;
  ReducedMatrixRangeMap_ = 0;
  OrigReducedMatrixDomainMap_ = 0;
  Full2ReducedRHSImporter_ = 0;
  Full2ReducedLHSImporter_ = 0;
  RedistributeDomainExporter_ = 0;

  ColSingletonRowLIDs_ = 0;
  ColSingletonColLIDs_ = 0;
  ColSingletonPivotLIDs_ = 0;
  ColSingletonPivots_ = 0;

  AbsoluteThreshold_ = 0;
  RelativeThreshold_ = 0;

  NumMyRowSingletons_ = -1;
  NumMyColSingletons_ = -1;
  NumGlobalRowSingletons_ = -1;
  NumGlobalColSingletons_ = -1;
  RatioOfDimensions_ = -1.0;
  RatioOfNonzeros_ = -1.0;

  HaveReducedProblem_ = false;
  UserDefinedEliminateMaps_ = false;
  AnalysisDone_ = false;
  SymmetricElimination_ = true;

  tempExportX_ = 0;
  tempX_ = 0;
  tempB_ = 0;

  Indices_int_ = 0;
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  Indices_LL_ = 0;
#endif

  RowMapColors_ = 0;
  ColMapColors_ = 0;

  FullMatrixIsCrsMatrix_ = false;
  MaxNumMyEntries_ = 0;
  return;
}
//==============================================================================
int LinearProblem_CrsSingletonFilter::Analyze(Epetra_RowMatrix * fullMatrix) {
  PRINT(outputRank_, "Analyze() A");
  FullMatrix_ = fullMatrix;

  if (AnalysisDone_) EPETRA_CHK_ERR(-1); // Analysis already done once.  Cannot do it again
  if (fullMatrix==0) EPETRA_CHK_ERR(-2); // Input matrix pointer is zero
  if (fullMatrix->NumGlobalRows64()==0) EPETRA_CHK_ERR(-3); // Full matrix has zero dimension.
  if (fullMatrix->NumGlobalNonzeros64()==0) EPETRA_CHK_ERR(-4); // Full matrix has no nonzero terms.

  PRINT(outputRank_, "Analyze() B");
  MPI_Barrier(MPI_COMM_WORLD);
  // First check for columns with single entries and find columns with singleton rows
  Epetra_IntVector ColProfiles(FullMatrixColMap()); ColProfiles.PutValue(0);
  Epetra_IntVector ColHasRowWithSingleton(FullMatrixColMap()); ColHasRowWithSingleton.PutValue(0);

  // RowIDs[j] will contain the local row ID associated with the jth column,
  // if the jth col has a single entry
  Epetra_IntVector RowIDs(FullMatrixColMap()); RowIDs.PutValue(-1);

  PRINT(outputRank_, "Analyze() C");
  MPI_Barrier(MPI_COMM_WORLD);
  // Define MapColoring objects
  RowMapColors_ = new Epetra_MapColoring(FullMatrixRowMap());  // Initial colors are all 0
  ColMapColors_ = new Epetra_MapColoring(FullMatrixColMap());
  Epetra_MapColoring & rowMapColors = *RowMapColors_;
  Epetra_MapColoring & colMapColors = *ColMapColors_;

  int NumMyRows = fullMatrix->NumMyRows();
  int NumMyCols = fullMatrix->NumMyCols();

  // Set up for accessing full matrix.  Will do so row-by-row.
  EPETRA_CHK_ERR(InitFullMatrixAccess());

  // Scan matrix for singleton rows, build up column profiles
  int NumIndices = 1;
  PRINT(outputRank_, "Analyze() Initialize NumIndices = " << NumIndices);
  int * Indices;
  NumMyRowSingletons_ = 0;
  PRINT(outputRank_, "Analyze() D");
  PRINT(outputRank_, "Analyze() E");
  for (int i=0; i<NumMyRows; i++) {
    PRINT(outputRank_, "Analyze() E  i = " << i);
    PRINT(outputRank_, "Analyze() E  NumIndices = " << NumIndices);
    // Get ith row
    EPETRA_CHK_ERR(GetRow(i, NumIndices, Indices));
    PRINT(outputRank_, "Analyze() E  NumIndices = " << NumIndices);
    if (outputRank_) std::cout << "petra : Analyze() E      localIndices = ";
  MPI_Barrier(MPI_COMM_WORLD);
    for (int ii = 0; ii  < NumIndices; ++ii) {
        std::cout << Indices[ii] << " ";
    } std::cout << std::endl;
    PRINT(outputRank_, "Analyze() F");
    for (int j=0; j<NumIndices; j++) {
      PRINT(outputRank_, "Analyze() F    j = " << j);
      int ColumnIndex = Indices[j];
      PRINT(outputRank_, "Analyze() F    ColumnIndex = " << ColumnIndex);
      PRINT(outputRank_, "Analyze() F    ColProfiles[ColumnIndex] = " << ColProfiles[ColumnIndex]);
      ColProfiles[ColumnIndex]++; // Increment column count
      PRINT(outputRank_, "Analyze() F    ColProfiles[ColumnIndex] = " << ColProfiles[ColumnIndex]);
      // Record local row ID for current column
      // will use to identify row to eliminate if column is a singleton
      RowIDs[ColumnIndex] = i;
      PRINT(outputRank_, "Analyze() F    RowIDs[ColumnIndex] = " << RowIDs[ColumnIndex]);
    }
    PRINT(outputRank_, "Analyze() G");
    // If row has single entry, color it and associated column with color=1
    if (NumIndices==1) {
      int j2 = Indices[0];
      ColHasRowWithSingleton[j2]++;
      rowMapColors[i] = 1;
      colMapColors[j2] = 1;
      NumMyRowSingletons_++;
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

  Epetra_IntVector NewColProfiles(ColProfiles);

  // If RowMatrixImporter is non-trivial, we need to perform a gather/scatter to accumulate results

  if (fullMatrix->RowMatrixImporter()!=0) {
    Epetra_IntVector tmpVec(FullMatrixDomainMap()); // Use for gather/scatter of column vectors
    EPETRA_CHK_ERR(tmpVec.PutValue(0));
    EPETRA_CHK_ERR(tmpVec.Export(ColProfiles, *fullMatrix->RowMatrixImporter(), Add));
    EPETRA_CHK_ERR(ColProfiles.Import(tmpVec, *fullMatrix->RowMatrixImporter(), Insert));

    EPETRA_CHK_ERR(tmpVec.PutValue(0));
    EPETRA_CHK_ERR(tmpVec.Export(ColHasRowWithSingleton, *fullMatrix->RowMatrixImporter(), Add));
    EPETRA_CHK_ERR(ColHasRowWithSingleton.Import(tmpVec, *fullMatrix->RowMatrixImporter(), Insert));
  }

  PRINT(outputRank_, "Analyze() I    ColProfiles");
  PRINT(outputRank_, "Analyze() I    NumIndices = " << NumIndices);
  MPI_Barrier(MPI_COMM_WORLD);
  for (int j=0; j<NumIndices; j++) {
    PRINT(outputRank_, "Analyze() I    j = " << j);
    int ColumnIndex = Indices[j];
    PRINT(outputRank_, "Analyze() I    ColumnIndex = " << ColumnIndex);
    PRINT(outputRank_, "Analyze() I    ColProfiles[ColumnIndex] = " << ColProfiles[ColumnIndex]);
  }
  PRINT(outputRank_, "Analyze() I    ColHasRowWithSingleton");
  PRINT(outputRank_, "Analyze() I    NumIndices = " << NumIndices);
  for (int j=0; j<NumIndices; j++) {
    PRINT(outputRank_, "Analyze() I    j = " << j);
    int ColumnIndex = Indices[j];
    PRINT(outputRank_, "Analyze() I    ColumnIndex = " << ColumnIndex);
    PRINT(outputRank_, "Analyze() I    ColHasRowWithSingleton[ColumnIndex] = " << ColProfiles[ColumnIndex]);
  }

  // ColProfiles now contains the nonzero column entry count for all columns that have
  // an entry on this processor.
  // ColHasRowWithSingleton now contains a count of singleton rows associated with the corresponding
  // local column.  Next we check to make sure no column is associated with more than one singleton row.

  if (ColHasRowWithSingleton.MaxValue()>1) {
    EPETRA_CHK_ERR(-2); // At least one col is associated with two singleton rows, can't handle it.
  }
  PRINT(outputRank_, "Analyze() J    ColHasRowWithSingleton MaxValue = " << ColHasRowWithSingleton.MaxValue());

  Epetra_IntVector RowHasColWithSingleton(fullMatrix->RowMatrixRowMap()); // Use to check for errors
  RowHasColWithSingleton.PutValue(0);

  NumMyColSingletons_ = 0;
  // Count singleton columns (that were not already counted as singleton rows)
  for (int j=0; j<NumMyCols; j++) {
    int i2 = RowIDs[j];
    PRINT(outputRank_, "Analyze() K    i2 = " << i2);
    // Check if column is a singleton
    if (ColProfiles[j]==1 ) {
      // Check to make sure RowID is not invalid
//      assert(i!=-1);
      // Check to see if this column already eliminated by the row check above
      if (rowMapColors[i2]!=1) {
        RowHasColWithSingleton[i2]++; // Increment col singleton counter for ith row
        rowMapColors[i2] = 2; // Use 2 for now, to distinguish between row eliminated directly or via column singletons
        colMapColors[j] = 1;
        NumMyColSingletons_++;
        // If we delete a row, we need to keep track of associated column entries that were also deleted
        // in case all entries in a column are eventually deleted, in which case the column should
        // also be deleted.
        EPETRA_CHK_ERR(GetRow(i2, NumIndices, Indices));
        for (int jj=0; jj<NumIndices; jj++) {
          NewColProfiles[Indices[jj]]--;
        }
      }
    }
    // Check if some other processor eliminated this column
    else if (ColHasRowWithSingleton[j]==1 && rowMapColors[i2]!=1) {
        colMapColors[j] = 1;
    }
  }
  for (int j=0; j<NumMyCols; j++) {
    PRINT(outputRank_, "Analyze() DEBUG    ColProfiles["<<j<<"] = " << ColProfiles[j]);
  }


  PRINT(outputRank_, "Analyze() K    NumMyColSingletons_ = " << NumMyColSingletons_);
  auto size1 = RowHasColWithSingleton.MyLength();
  for (int j=0; j<size1; j++) {
    if (RowHasColWithSingleton[j] == 0) continue;
    PRINT(outputRank_, "Analyze() K    RowHasColWithSingleton["<<j<<"] = " << RowHasColWithSingleton[j]);
  }
  size1 = NumMyRows;
  for (int j=0; j<size1; j++) {
    if (rowMapColors[j] != 0) {
      PRINT_NB(outputRank_, "Analyze() K    rowMapColors["<<j<<"] = " << rowMapColors[j]);
    }
    if (colMapColors[j] != 0) {
      PRINT_NB(outputRank_, "Analyze() K    colMapColors["<<j<<"] = " << colMapColors[j]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  size1 = NewColProfiles.MyLength();
  for (int j=0; j<size1; j++) {
    PRINT(outputRank_, "Analyze() K    NewColProfiles["<<j<<"] = " << NewColProfiles[j]);
  }



  if (RowHasColWithSingleton.MaxValue()>1) {
    EPETRA_CHK_ERR(-3); // At lease one row is associated with two singleton cols, can't handle it.
  }

 // Generate arrays that keep track of column singleton row, col and pivot info needed for post-solve phase
  EPETRA_CHK_ERR(CreatePostSolveArrays(RowIDs, rowMapColors, ColProfiles, NewColProfiles,
                                       ColHasRowWithSingleton));

  for (int j=0; j<NumMyCols; j++) {
    PRINT(outputRank_, "Analyze() L    RowIDs["<<j<<"] = " << RowIDs[j]);
    PRINT(outputRank_, "Analyze() L    ColProfiles["<<j<<"] = " << ColProfiles[j]);
  }
  size1 = NewColProfiles.MyLength();
  for (int j=0; j<size1; j++) {
    PRINT(outputRank_, "Analyze() L    NewColProfiles["<<j<<"] = " << NewColProfiles[j]);
  }
  size1 = ColHasRowWithSingleton.MyLength();
  for (int j=0; j<size1; j++) {
    PRINT(outputRank_, "Analyze() L    ColHasRowWithSingleton["<<j<<"] = " << ColHasRowWithSingleton[j]);
  }


  for (int i=0; i<NumMyRows; i++) {
    if (rowMapColors[i]==2) rowMapColors[i] = 1; // Convert all eliminated rows to same color
  }
  size1 = NumMyRows;
  for (int j=0; j<size1; j++) {
    if (rowMapColors[j] != 0) {
      PRINT_NB(outputRank_, "Analyze() L    rowMapColors["<<j<<"] = " << rowMapColors[j]);
    }
  }

  fullMatrix->RowMatrixRowMap().Comm().SumAll(&NumMyRowSingletons_, &NumGlobalRowSingletons_, 1);
  fullMatrix->RowMatrixRowMap().Comm().SumAll(&NumMyColSingletons_, &NumGlobalColSingletons_, 1);
  PRINT(outputRank_, "Analyze() M    NumGlobalRowSingletons_ = " << NumGlobalRowSingletons_);
  PRINT(outputRank_, "Analyze() M    NumGlobalColSingletons_ = " << NumGlobalColSingletons_);
  AnalysisDone_ = true;
  return(0);
}
//==============================================================================
template<typename int_type>
int LinearProblem_CrsSingletonFilter::TConstructReducedProblem(Epetra_LinearProblem * Problem) {

  PRINT(outputRank_, "ConstructReducedProblem() A");
  
  int i, j;
  if (HaveReducedProblem_) EPETRA_CHK_ERR(-1); // Setup already done once.  Cannot do it again
  if (Problem==0) EPETRA_CHK_ERR(-2); // Null problem pointer

  FullProblem_ = Problem;
  FullMatrix_ = dynamic_cast<Epetra_RowMatrix *>(Problem->GetMatrix());
  if (FullMatrix_==0) EPETRA_CHK_ERR(-3); // Need a RowMatrix
  if (Problem->GetRHS()==0) EPETRA_CHK_ERR(-4); // Need a RHS
  if (Problem->GetLHS()==0) EPETRA_CHK_ERR(-5); // Need a LHS
  // Generate reduced row and column maps

  PRINT(outputRank_, "ConstructReducedProblem() A    SingletonsDetected() = " << SingletonsDetected());
  if ( SingletonsDetected() ) {

    Epetra_MapColoring & rowMapColors = *RowMapColors_;
    Epetra_MapColoring & colMapColors = *ColMapColors_;

    auto size1 = FullMatrix()->NumMyRows();
    for (int jj=0; jj<size1; jj++) {
      if (rowMapColors[jj] != 0) {
        PRINT_NB(outputRank_, "ConstructReducedProblem() A    rowMapColors["<<jj<<"] = " << rowMapColors[jj]);
      }
      if (colMapColors[jj] != 0) {
        PRINT_NB(outputRank_, "ConstructReducedProblem() A    colMapColors["<<jj<<"] = " << colMapColors[jj]);
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }

    ReducedMatrixRowMap_ = rowMapColors.GenerateMap(0);
    ReducedMatrixColMap_ = colMapColors.GenerateMap(0);

    ReducedMatrixRowMap_->Print(std::cout);
    ReducedMatrixColMap_->Print(std::cout);

    // Create domain and range map colorings by exporting map coloring of column and row maps

    if (FullMatrix()->RowMatrixImporter()!=0) {
      Epetra_MapColoring DomainMapColors(FullMatrixDomainMap());
      EPETRA_CHK_ERR(DomainMapColors.Export(*ColMapColors_, *FullMatrix()->RowMatrixImporter(), AbsMax));
      OrigReducedMatrixDomainMap_ = DomainMapColors.GenerateMap(0);
    }
    else
      OrigReducedMatrixDomainMap_ = ReducedMatrixColMap_;
    OrigReducedMatrixDomainMap_->Print(std::cout);

    if (FullMatrixIsCrsMatrix_) {
      if (FullCrsMatrix()->Exporter()!=0) { // Non-trivial exporter
        Epetra_MapColoring RangeMapColors(FullMatrixRangeMap());
        EPETRA_CHK_ERR(RangeMapColors.Export(*RowMapColors_, *FullCrsMatrix()->Exporter(),
                                           AbsMax));
        ReducedMatrixRangeMap_ = RangeMapColors.GenerateMap(0);
      }
      else
        ReducedMatrixRangeMap_ = ReducedMatrixRowMap_;
    }
    else
      ReducedMatrixRangeMap_ = ReducedMatrixRowMap_;
    ReducedMatrixRangeMap_->Print(std::cout);


    PRINT(outputRank_, "ConstructReducedProblem() C    SymmetricElimination_ = " << SymmetricElimination_);

    // Check to see if the reduced system domain and range maps are the same.
    // If not, we need to remap entries of the LHS multivector so that they are distributed
    // conformally with the rows of the reduced matrix and the RHS multivector
    SymmetricElimination_ = ReducedMatrixRangeMap_->SameAs(*OrigReducedMatrixDomainMap_);
    if (!SymmetricElimination_)
      ConstructRedistributeExporter(OrigReducedMatrixDomainMap_, ReducedMatrixRangeMap_,
                                    RedistributeDomainExporter_, ReducedMatrixDomainMap_);
    else {
      ReducedMatrixDomainMap_ = OrigReducedMatrixDomainMap_;
      OrigReducedMatrixDomainMap_ = 0;
      RedistributeDomainExporter_ = 0;
    }

    // Create pointer to Full RHS, LHS
    Epetra_MultiVector * FullRHS = FullProblem()->GetRHS();
    Epetra_MultiVector * FullLHS = FullProblem()->GetLHS();
    int NumVectors = FullLHS->NumVectors();

    PRINT(outputRank_, "ConstructReducedProblem() D    NumVectors = " << NumVectors);

    PRINT(outputRank_, "ConstructReducedProblem() E    ReducedMatrixDomainMap() = " << ReducedMatrixDomainMap());
    //PRINT(outputRank_, "ConstructReducedProblem() E    FullMatrixDomainMap() = " << FullMatrixDomainMap());
    PRINT(outputRank_, "ConstructReducedProblem() E    ReducedMatrixRowMap() = " << ReducedMatrixRowMap());
    //PRINT(outputRank_, "ConstructReducedProblem() E    FullRHS->Map() = " << FullRHS->Map());
  
    // Create importers
    Full2ReducedLHSImporter_ = new Epetra_Import(*ReducedMatrixDomainMap(), FullMatrixDomainMap());
    ReducedMatrixRowMap_->Print(std::cout);
    FullRHS->Map().Print(std::cout);
    Full2ReducedRHSImporter_ = new Epetra_Import(*ReducedMatrixRowMap(), FullRHS->Map());
    if (Full2ReducedRHSImporter_ != nullptr) {
        // Print detailed information about the importer
        Full2ReducedRHSImporter_->Print(std::cout);
    } else {
        std::cout << "Full2ReducedRHSImporter_ is null." << std::endl;
    }


    // Construct Reduced Matrix
    ReducedMatrix_ = Teuchos::rcp( new Epetra_CrsMatrix(Copy, *ReducedMatrixRowMap(), *ReducedMatrixColMap(), 0) );

    if (outputRank_) PrintEpetraCrsMatrix(*ReducedMatrix_);

    // Create storage for temporary X values due to explicit elimination of rows
    tempExportX_ = new Epetra_MultiVector(FullMatrixColMap(), NumVectors);

    int NumEntries;
    double * Values;
    int NumMyRows = FullMatrix()->NumMyRows();
    int ColSingletonCounter = 0;
    for (i=0; i<NumMyRows; i++) {
      int_type curGRID = (int_type) FullMatrixRowMap().GID64(i);
      if (ReducedMatrixRowMap()->MyGID(curGRID)) { // Check if this row should go into reduced matrix
        int_type * Indices;
        EPETRA_CHK_ERR(GetRowGCIDs(i, NumEntries, Values, Indices)); // Get current row (Indices are global)


        //// Debug: Print row and column indices being inserted
        //std::cout << "Inserting into row " << curGRID << " with column indices: ";
        //for (int j = 0; j < NumEntries; ++j) {
        //    std::cout << Indices[j] << " ";
        //}
        //std::cout << std::endl;
        //
        //// Debug: Print the column map on this process
        //std::cout << "Column map on this process: ";
        //for (int j = 0; j < ReducedMatrixColMap()->NumMyElements(); ++j) {
        //    std::cout << ReducedMatrixColMap()->GID64(j) << " "; // Get global column indices
        //}
        //std::cout << std::endl;
        //
        //// Debug: Check for invalid column indices
        //for (int j = 0; j < NumEntries; ++j) {
        //    if (!ReducedMatrixColMap()->MyGID(Indices[j])) { // Check if column index is in the column map
        //        std::cout << "Error: Column index " << Indices[j] << " is not in the column map on this process!" << std::endl;
        //    }
        //}


        int ierr = ReducedMatrix()->InsertGlobalValues(curGRID, NumEntries,
                                                     Values, Indices); // Insert into reduce matrix
        // Positive errors will occur because we are submitting col entries that are not part of
        // reduced system.  However, because we specified a column map to the ReducedMatrix constructor
        // these extra column entries will be ignored and we will be politely reminded by a positive
        // error code
        if (ierr<0) EPETRA_CHK_ERR(ierr);
      }
      else {
        int * Indices;
        PRINT(outputRank_, "ConstructReducedProblem() F    i, NumEntries = " << i << ", " << NumEntries);
        EPETRA_CHK_ERR(GetRow(i, NumEntries, Values, Indices)); // Get current row
        if (NumEntries==1) {
          double pivot = Values[0];
          if (pivot==0.0) EPETRA_CHK_ERR(-1); // Encountered zero row, unable to continue
          int indX = Indices[0];
          for (j=0; j<NumVectors; j++)
            (*tempExportX_)[j][indX] = (*FullRHS)[j][i]/pivot;
        }
        // Otherwise, this is a singleton column and we will scan for the pivot element needed
        // for post-solve equations
        else {
          int targetCol = ColSingletonColLIDs_[ColSingletonCounter];
          for (j=0; j<NumEntries; j++) {
            if (Indices[j]==targetCol) {
              double pivot = Values[j];
              if (pivot==0.0) EPETRA_CHK_ERR(-2); // Encountered zero column, unable to continue
              ColSingletonPivotLIDs_[ColSingletonCounter] = j; // Save for later use
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
    EPETRA_CHK_ERR(ReducedMatrix()->FillComplete(*ReducedMatrixDomainMap(), *ReducedMatrixRangeMap()));

    // 1) The vector ColProfiles has column nonzero counts for each processor's contribution
    // Construct Reduced LHS (Puts any initial guess values into reduced system)

    ReducedLHS_ = new Epetra_MultiVector(*ReducedMatrixDomainMap(), NumVectors);
    EPETRA_CHK_ERR(ReducedLHS_->Import(*FullLHS, *Full2ReducedLHSImporter_, Insert));
    if (outputRank_) printEpetraMultiVector("ConstructReducedProblem() G    ReducedLHS_ = ", ReducedLHS_);
    FullLHS->PutScalar(0.0); // zero out Full LHS since we will inject values as we get them

    // Construct Reduced RHS

    // First compute influence of already-known values of X on RHS
    tempX_ = new Epetra_MultiVector(FullMatrixDomainMap(), NumVectors);
    tempB_ = new Epetra_MultiVector(FullRHS->Map(), NumVectors);

    //Inject known X values into tempX for purpose of computing tempB = FullMatrix*tempX
    // Also inject into full X since we already know the solution

    if (FullMatrix()->RowMatrixImporter()!=0) {
      EPETRA_CHK_ERR(tempX_->Export(*tempExportX_, *FullMatrix()->RowMatrixImporter(), Add));
      EPETRA_CHK_ERR(FullLHS->Export(*tempExportX_, *FullMatrix()->RowMatrixImporter(), Add));
    }
    else {
      tempX_->Update(1.0, *tempExportX_, 0.0);
      FullLHS->Update(1.0, *tempExportX_, 0.0);

      PRINT(outputRank_, "ConstructReducedProblem() L");
      if (outputRank_) printEpetraMultiVector("ConstructReducedProblem() L    tempX_ = ", tempX_);
      if (outputRank_) printEpetraMultiVector("ConstructReducedProblem() L    FullLHS = ", FullLHS);

    }

    EPETRA_CHK_ERR(FullMatrix()->Multiply(false, *tempX_, *tempB_));

    EPETRA_CHK_ERR(tempB_->Update(1.0, *FullRHS, -1.0)); // tempB now has influence of already-known X values

    PRINT(outputRank_, "ConstructReducedProblem() M");
    if (outputRank_) printEpetraMultiVector("ConstructReducedProblem() M    tempB_ = ", tempB_);

    ReducedRHS_ = new Epetra_MultiVector(*ReducedMatrixRowMap(), FullRHS->NumVectors());
    if (outputRank_) printEpetraMultiVector("ConstructReducedProblem() S    ReducedRHS_ = ", ReducedRHS_);
    EPETRA_CHK_ERR(ReducedRHS_->Import(*tempB_, *Full2ReducedRHSImporter_, Insert));

    if (Full2ReducedRHSImporter_ != nullptr) {
        // Print detailed information about the importer
        Full2ReducedRHSImporter_->Print(std::cout);
    } else {
        std::cout << "Full2ReducedRHSImporter_ is null." << std::endl;
    }

    PRINT(outputRank_, "ConstructReducedProblem() S");
    if (outputRank_) printEpetraMultiVector("ConstructReducedProblem() S    ReducedLHS_ = ", ReducedLHS_);
    if (outputRank_) printEpetraMultiVector("ConstructReducedProblem() S    ReducedRHS_ = ", ReducedRHS_);

    // Finally construct Reduced Linear Problem
    ReducedProblem_ = Teuchos::rcp( new Epetra_LinearProblem(ReducedMatrix_.get(), ReducedLHS_, ReducedRHS_) );
  }
  else {

    // There are no singletons, so don't bother building a reduced problem.
    ReducedProblem_ = Teuchos::rcp( Problem, false );
    ReducedMatrix_ = Teuchos::rcp( dynamic_cast<Epetra_CrsMatrix *>(Problem->GetMatrix()), false );
  }


  PRINT(outputRank_, "ConstructReducedProblem() Z");
  if (outputRank_) PrintEpetraCrsMatrix(*ReducedMatrix_);


  double fn = (double) FullMatrix()->NumGlobalRows64();
  double fnnz = (double) FullMatrix()->NumGlobalNonzeros64();
  PRINT(outputRank_, "ConstructReducedProblem() Z    fn = " << fn);
  PRINT(outputRank_, "ConstructReducedProblem() Z    fnnz = " << fnnz);

  double rn = (double) ReducedMatrix()->NumGlobalRows64();
  double rnnz = (double) ReducedMatrix()->NumGlobalNonzeros64();
  PRINT(outputRank_, "ConstructReducedProblem() Z    rn = " << rn);
  PRINT(outputRank_, "ConstructReducedProblem() Z    rnnz = " << rnnz);

  RatioOfDimensions_ = rn/fn;
  RatioOfNonzeros_ = rnnz/fnnz;
  HaveReducedProblem_ = true;

  return(0);
}

int LinearProblem_CrsSingletonFilter::ConstructReducedProblem(Epetra_LinearProblem * Problem) {
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(Problem->GetMatrix()->RowMatrixRowMap().GlobalIndicesInt()) {
    return TConstructReducedProblem<int>(Problem);
  }
  else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(Problem->GetMatrix()->RowMatrixRowMap().GlobalIndicesLongLong()) {
    return TConstructReducedProblem<long long>(Problem);
  }
  else
#endif
    throw "LinearProblem_CrsSingletonFilter::ConstructReducedProblem: ERROR, GlobalIndices type unknown.";
}

//==============================================================================
template<typename int_type>
int LinearProblem_CrsSingletonFilter::TUpdateReducedProblem(Epetra_LinearProblem * Problem) {

  PRINT(outputRank_, "UpdateReducedProblem() A");

  int i, j;

  if (Problem==0) EPETRA_CHK_ERR(-1); // Null problem pointer

  FullProblem_ = Problem;
  FullMatrix_ = dynamic_cast<Epetra_RowMatrix *>(Problem->GetMatrix());
  if (FullMatrix_==0) EPETRA_CHK_ERR(-2); // Need a RowMatrix
  if (Problem->GetRHS()==0) EPETRA_CHK_ERR(-3); // Need a RHS
  if (Problem->GetLHS()==0) EPETRA_CHK_ERR(-4); // Need a LHS
  if (!HaveReducedProblem_) EPETRA_CHK_ERR(-5); // Must have set up reduced problem

  PRINT(outputRank_, "UpdateReducedProblem() A    SingletonsDetected() = " << SingletonsDetected());
  if ( SingletonsDetected() ) {

    // Create pointer to Full RHS, LHS
    Epetra_MultiVector * FullRHS = FullProblem()->GetRHS();
    Epetra_MultiVector * FullLHS = FullProblem()->GetLHS();

    int NumVectors = FullLHS->NumVectors();
    PRINT(outputRank_, "UpdateReducedProblem() B    NumVectors = " << NumVectors);
    tempExportX_->PutScalar(0.0);

    int NumEntries;
    double * Values;
    int NumMyRows = FullMatrix()->NumMyRows();
    int ColSingletonCounter = 0;
    for (i=0; i<NumMyRows; i++) {
      int_type curGRID = (int_type) FullMatrixRowMap().GID64(i);
      if (ReducedMatrixRowMap()->MyGID(curGRID)) { // Check if this row should go into reduced matrix
        int_type * Indices;
        EPETRA_CHK_ERR(GetRowGCIDs(i, NumEntries, Values, Indices)); // Get current row (indices global)
        int ierr = ReducedMatrix()->ReplaceGlobalValues(curGRID, NumEntries,
                                                      Values, Indices);
        // Positive errors will occur because we are submitting col entries that are not part of
        // reduced system.  However, because we specified a column map to the ReducedMatrix constructor
        // these extra column entries will be ignored and we will be politely reminded by a positive
        // error code
        if (ierr<0) EPETRA_CHK_ERR(ierr);
      }
      // Otherwise if singleton row we explicitly eliminate this row and solve for corresponding X value
      else {
        int * Indices;
        EPETRA_CHK_ERR(GetRow(i, NumEntries, Values, Indices)); // Get current row
        if (NumEntries==1) {
          double pivot = Values[0];
          if (pivot==0.0) EPETRA_CHK_ERR(-1); // Encountered zero row, unable to continue
          int indX = Indices[0];
          for (j=0; j<NumVectors; j++)
            (*tempExportX_)[j][indX] = (*FullRHS)[j][i]/pivot;
        }
        // Otherwise, this is a singleton column and we will scan for the pivot element needed
        // for post-solve equations
        else {
          j = ColSingletonPivotLIDs_[ColSingletonCounter];
          double pivot = Values[j];
          if (pivot==0.0) EPETRA_CHK_ERR(-2); // Encountered zero column, unable to continue
          ColSingletonPivots_[ColSingletonCounter] = pivot;
          ColSingletonCounter++;
        }
      }
    }

    assert(ColSingletonCounter==NumMyColSingletons_); // Sanity test

    // Update Reduced LHS (Puts any initial guess values into reduced system)

    ReducedLHS_->PutScalar(0.0); // zero out Reduced LHS
    EPETRA_CHK_ERR(ReducedLHS_->Import(*FullLHS, *Full2ReducedLHSImporter_, Insert));
    FullLHS->PutScalar(0.0); // zero out Full LHS since we will inject values as we get them

    // Construct Reduced RHS

    // Zero out temp space
    tempX_->PutScalar(0.0);
    tempB_->PutScalar(0.0);

    //Inject known X values into tempX for purpose of computing tempB = FullMatrix*tempX
    // Also inject into full X since we already know the solution

    if (FullMatrix()->RowMatrixImporter()!=0) {
      EPETRA_CHK_ERR(tempX_->Export(*tempExportX_, *FullMatrix()->RowMatrixImporter(), Add));
      EPETRA_CHK_ERR(FullLHS->Export(*tempExportX_, *FullMatrix()->RowMatrixImporter(), Add));
    }
    else {
      tempX_->Update(1.0, *tempExportX_, 0.0);
      FullLHS->Update(1.0, *tempExportX_, 0.0);
    }

    EPETRA_CHK_ERR(FullMatrix()->Multiply(false, *tempX_, *tempB_));

    EPETRA_CHK_ERR(tempB_->Update(1.0, *FullRHS, -1.0)); // tempB now has influence of already-known X values

    ReducedRHS_->PutScalar(0.0);
    EPETRA_CHK_ERR(ReducedRHS_->Import(*tempB_, *Full2ReducedRHSImporter_, Insert));

    if (outputRank_) printEpetraMultiVector("UpdateReducedProblem() Z    FullLHS = ", FullLHS);
    if (outputRank_) printEpetraMultiVector("UpdateReducedProblem() Z    FullRHS = ", FullRHS);
  }
  else {

    // There are no singletons, so don't bother building a reduced problem.
    ReducedProblem_ = Teuchos::rcp( Problem, false );
    ReducedMatrix_ = Teuchos::rcp( dynamic_cast<Epetra_CrsMatrix *>(Problem->GetMatrix()), false );
  }

  if (outputRank_) PrintEpetraRowMatrix(*FullMatrix());
  if (outputRank_) printEpetraMultiVector("UpdateReducedProblem() Z    tempX_ = ", tempX_);
  if (outputRank_) printEpetraMultiVector("UpdateReducedProblem() Z    tempB_ = ", tempB_);
  if (outputRank_) printEpetraMultiVector("UpdateReducedProblem() Z    tempExportX_ = ", tempExportX_);
  if (outputRank_) PrintEpetraRowMatrix(*ReducedMatrix());
  if (outputRank_) printEpetraMultiVector("UpdateReducedProblem() Z    ReducedLHS_ = ", ReducedLHS_);
  if (outputRank_) printEpetraMultiVector("UpdateReducedProblem() Z    ReducedRHS_ = ", ReducedRHS_);

  return(0);
}

int LinearProblem_CrsSingletonFilter::UpdateReducedProblem(Epetra_LinearProblem * Problem) {
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(Problem->GetMatrix()->RowMatrixRowMap().GlobalIndicesInt()) {
    return TUpdateReducedProblem<int>(Problem);
  }
  else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(Problem->GetMatrix()->RowMatrixRowMap().GlobalIndicesLongLong()) {
    return TUpdateReducedProblem<long long>(Problem);
  }
  else
#endif
    throw "LinearProblem_CrsSingletonFilter::UpdateReducedProblem: ERROR, GlobalIndices type unknown.";
}
//==============================================================================
template<typename int_type>
int LinearProblem_CrsSingletonFilter::TConstructRedistributeExporter(Epetra_Map * SourceMap, Epetra_Map * TargetMap,
                                                             Epetra_Export * & RedistributeExporter,
                                                             Epetra_Map * & RedistributeMap) {

  int_type IndexBase = (int_type) SourceMap->IndexBase64();
  if (IndexBase!=(int_type) TargetMap->IndexBase64()) EPETRA_CHK_ERR(-1);

  const Epetra_Comm & Comm = TargetMap->Comm();

  int TargetNumMyElements = TargetMap->NumMyElements();
  int SourceNumMyElements = SourceMap->NumMyElements();

  // ContiguousTargetMap has same number of elements per PE as TargetMap, but uses contigious indexing
  Epetra_Map ContiguousTargetMap((int_type) -1, TargetNumMyElements, IndexBase,Comm);

  // Same for ContiguousSourceMap
  Epetra_Map ContiguousSourceMap((int_type) -1, SourceNumMyElements, IndexBase, Comm);

  assert(ContiguousSourceMap.NumGlobalElements64()==ContiguousTargetMap.NumGlobalElements64());

  // Now create a vector that contains the global indices of the Source Epetra_MultiVector
  int_type* SourceMapMyGlobalElements = 0;
  SourceMap->MyGlobalElementsPtr(SourceMapMyGlobalElements);
  typename Epetra_GIDTypeVector<int_type>::impl SourceIndices(View, ContiguousSourceMap, SourceMapMyGlobalElements);

  // Create an exporter to send the SourceMap global IDs to the target distribution
  Epetra_Export Exporter(ContiguousSourceMap, ContiguousTargetMap);

  // Create a vector to catch the global IDs in the target distribution
  typename Epetra_GIDTypeVector<int_type>::impl TargetIndices(ContiguousTargetMap);
  TargetIndices.Export(SourceIndices, Exporter, Insert);

  // Create a new map that describes how the Source MultiVector should be laid out so that it has
  // the same number of elements on each processor as the TargetMap
  RedistributeMap = new Epetra_Map((int_type) -1, TargetNumMyElements, TargetIndices.Values(), IndexBase, Comm);

  // This exporter will finally redistribute the Source MultiVector to the same layout as the TargetMap
  RedistributeExporter = new Epetra_Export(*SourceMap, *RedistributeMap);
  return(0);
}

int LinearProblem_CrsSingletonFilter::ConstructRedistributeExporter(Epetra_Map * SourceMap, Epetra_Map * TargetMap,
                                                             Epetra_Export * & RedistributeExporter,
                                                             Epetra_Map * & RedistributeMap) {
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(SourceMap->GlobalIndicesInt() && TargetMap->GlobalIndicesInt()) {
    return TConstructRedistributeExporter<int>(SourceMap, TargetMap, RedistributeExporter, RedistributeMap);
  }
  else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(SourceMap->GlobalIndicesLongLong() && TargetMap->GlobalIndicesLongLong()) {
    return TConstructRedistributeExporter<long long>(SourceMap, TargetMap, RedistributeExporter, RedistributeMap);
  }
  else
#endif
    throw "LinearProblem_CrsSingletonFilter::ConstructRedistributeExporter: ERROR, GlobalIndices type unknown.";
}
//==============================================================================
int LinearProblem_CrsSingletonFilter::ComputeFullSolution() {

  if ( SingletonsDetected() ) {
    int jj, k;

    Epetra_MultiVector * FullLHS = FullProblem()->GetLHS();
    Epetra_MultiVector * FullRHS = FullProblem()->GetRHS();

    tempX_->PutScalar(0.0); tempExportX_->PutScalar(0.0);
    // Inject values that the user computed for the reduced problem into the full solution vector
    EPETRA_CHK_ERR(tempX_->Export(*ReducedLHS_, *Full2ReducedLHSImporter_, Add));

    FullLHS->Update(1.0, *tempX_, 1.0);

    // Next we will use our full solution vector which is populated with pre-filter solution
    // values and reduced system solution values to compute the sum of the row contributions
    // that must be subtracted to get the post-filter solution values

    EPETRA_CHK_ERR(FullMatrix()->Multiply(false, *FullLHS, *tempB_));

    // Finally we loop through the local rows that were associated with column singletons and compute the
    // solution for these equations.

    int NumVectors = tempB_->NumVectors();
    for (k=0; k<NumMyColSingletons_; k++) {
      int i = ColSingletonRowLIDs_[k];
      int j = ColSingletonColLIDs_[k];
      double pivot = ColSingletonPivots_[k];
      for (jj=0; jj<NumVectors; jj++)
        (*tempExportX_)[jj][j]= ((*FullRHS)[jj][i] - (*tempB_)[jj][i])/pivot;
    }

    // Finally, insert values from post-solve step and we are done!!!!

    if (FullMatrix()->RowMatrixImporter()!=0) {
      EPETRA_CHK_ERR(tempX_->Export(*tempExportX_, *FullMatrix()->RowMatrixImporter(), Add));
    }
    else {
      tempX_->Update(1.0, *tempExportX_, 0.0);
    }

    FullLHS->Update(1.0, *tempX_, 1.0);
  }

  return(0);
}
//==============================================================================
int LinearProblem_CrsSingletonFilter::InitFullMatrixAccess() {

  PRINT(outputRank_, "InitFullMatrixAccess() A");
  PRINT(outputRank_, "InitFullMatrixAccess() A    FullMatrix() = " << FullMatrix());
  if (outputRank_) PrintEpetraRowMatrix(*FullMatrix());
  MPI_Barrier(MPI_COMM_WORLD);

  MaxNumMyEntries_ = FullMatrix()->MaxNumEntries();
  PRINT(outputRank_, "InitFullMatrixAccess() B    localMaxNumRowEntries_ = " << MaxNumMyEntries_);

  // Cast to CrsMatrix, if possible.  Can save some work.
  FullCrsMatrix_ = dynamic_cast<Epetra_CrsMatrix *>(FullMatrix());
  PRINT(outputRank_, "InitFullMatrixAccess() C    FullCrsMatrix_ = " << FullCrsMatrix_);
  if (outputRank_) PrintEpetraCrsMatrix(*FullCrsMatrix_);
  MPI_Barrier(MPI_COMM_WORLD);
  FullMatrixIsCrsMatrix_ = (FullCrsMatrix_!=0); // Pointer is non-zero if cast worked
  PRINT(outputRank_, "InitFullMatrixAccess() D    FullMatrixIsCrsMatrix_ = " << FullMatrixIsCrsMatrix_);
  Indices_int_ = new int[MaxNumMyEntries_];
  PRINT(outputRank_, "InitFullMatrixAccess() E    Indices_ = " << Indices_int_);
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(FullMatrix()->RowMatrixRowMap().GlobalIndicesLongLong())
    Indices_LL_ = new long long[MaxNumMyEntries_];
#endif
  Values_.Size(MaxNumMyEntries_);
  PRINT(outputRank_, "InitFullMatrixAccess() F    Values_ = " << Values_);

  return(0);
}
//==============================================================================
int LinearProblem_CrsSingletonFilter::GetRow(int Row, int & NumIndices, int * & Indices) {

  PRINT(outputRank_, "GetRow() A      NumIndices = " << NumIndices);
  PRINT(outputRank_, "GetRow() B      NumIndices = " << NumIndices);
  PRINT(outputRank_, "GetRow() B      FullMatrixIsCrsMatrix_ = " << FullMatrixIsCrsMatrix_);

  if (FullMatrixIsCrsMatrix_) { // View of current row
    PRINT(outputRank_, "GetRow() C      NumIndices = " << NumIndices);
    PRINT(outputRank_, "GetRow() C      localRow = " << Row);
    EPETRA_CHK_ERR(FullCrsMatrix()->Graph().ExtractMyRowView(Row, NumIndices, Indices));
  }
  else { // Copy of current row (we must get the values, but we ignore them)
    PRINT(outputRank_, "GetRow() D      NumIndices = " << NumIndices);
    PRINT(outputRank_, "GetRow() D      localRow = " << Row);
    EPETRA_CHK_ERR(FullMatrix()->ExtractMyRowCopy(Row, MaxNumMyEntries_, NumIndices,
                                                  Values_.Values(), Indices_int_));
    Indices = Indices_int_;
  }
  PRINT(outputRank_, "GetRow() E      NumIndices = " << NumIndices);
  MPI_Barrier(MPI_COMM_WORLD);
  return(0);
}
//==============================================================================
int LinearProblem_CrsSingletonFilter::GetRow(int Row, int & NumIndices,
                                      double * & Values, int * & Indices) {

  PRINT(outputRank_, "GetRow(,,,) A   Row, NumIndices = " << Row << ", " << NumIndices);

  if (FullMatrixIsCrsMatrix_) { // View of current row
    EPETRA_CHK_ERR(FullCrsMatrix_->ExtractMyRowView(Row, NumIndices, Values, Indices));

        // Format Values and Indices for output
        std::ostringstream valuesStream;
        valuesStream << "{";
        for (int i = 0; i < NumIndices; ++i) {
            valuesStream << std::scientific << std::setprecision(17) << Values[i];
            if (i < NumIndices - 1) {
                valuesStream << ", ";
            }
        }
        valuesStream << "}";

        std::ostringstream indicesStream;
        indicesStream << "{";
        for (int i = 0; i < NumIndices; ++i) {
            indicesStream << Indices[i];
            if (i < NumIndices - 1) {
                indicesStream << ", ";
            }
        }
        indicesStream << "}";

        PRINT(outputRank_, "GetRow(,,,) B   Row, NumIndices = " << Row << ", " << NumIndices);
        PRINT(outputRank_, "GetRow(,,,) B   Values = " << valuesStream.str());
        PRINT(outputRank_, "GetRow(,,,) B   Indices = " << indicesStream.str());
  }
  else { // Copy of current row (we must get the values, but we ignore them)
    EPETRA_CHK_ERR(FullMatrix()->ExtractMyRowCopy(Row, MaxNumMyEntries_, NumIndices,
                                                  Values_.Values(), Indices_int_));
    Values = Values_.Values();
    Indices = Indices_int_;

        // Format Values and Indices for output
        std::ostringstream valuesStream;
        valuesStream << "{";
        for (int i = 0; i < NumIndices; ++i) {
            valuesStream << std::scientific << std::setprecision(17) << Values[i];
            if (i < NumIndices - 1) {
                valuesStream << ", ";
            }
        }
        valuesStream << "}";

        std::ostringstream indicesStream;
        indicesStream << "{";
        for (int i = 0; i < NumIndices; ++i) {
            indicesStream << Indices[i];
            if (i < NumIndices - 1) {
                indicesStream << ", ";
            }
        }
        indicesStream << "}";

        PRINT(outputRank_, "GetRow(,,,) B   Row, NumIndices = " << Row << ", " << NumIndices);
        PRINT(outputRank_, "GetRow(,,,) B   Values = " << valuesStream.str());
        PRINT(outputRank_, "GetRow(,,,) B   Indices = " << indicesStream.str());
  }
  return(0);
}
//==============================================================================
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int LinearProblem_CrsSingletonFilter::GetRowGCIDs(int Row, int & NumIndices,
                                           double * & Values, int * & GlobalIndices) {

  PRINT(outputRank_, "GetRowGCIDs() A      Row = " << Row);

    EPETRA_CHK_ERR(FullMatrix()->ExtractMyRowCopy(Row, MaxNumMyEntries_, NumIndices,
                                                  Values_.Values(), Indices_int_));
  PRINT(outputRank_, "GetRowGCIDs() A      NumIndices = " << NumIndices);
    for (int j=0; j<NumIndices; j++) Indices_int_[j] = FullMatrixColMap().GID(Indices_int_[j]);
    Values = Values_.Values();
    GlobalIndices = Indices_int_;

    // Print the Indices array
    std::cout << "petra : GetRowGCIDs() A      Row " << Row << " Indices: ";
    for (int j = 0; j < NumIndices; ++j) {
        std::cout << GlobalIndices[j] << " ";
    }
    std::cout << std::endl;

    // Print the Values array
    std::cout << "petra : GetRowGCIDs() A      Row " << Row << " Values: ";
    for (int j = 0; j < NumIndices; ++j) {
        std::cout << Values[j] << " ";
    }
    std::cout << std::endl;

  return(0);
}
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int LinearProblem_CrsSingletonFilter::GetRowGCIDs(int Row, int & NumIndices,
                                           double * & Values, long long * & GlobalIndices) {
    EPETRA_CHK_ERR(FullMatrix()->ExtractMyRowCopy(Row, MaxNumMyEntries_, NumIndices,
                                                  Values_.Values(), Indices_int_));
    for (int j=0; j<NumIndices; j++) Indices_LL_[j] = FullMatrixColMap().GID64(Indices_int_[j]);
    Values = Values_.Values();
    GlobalIndices = Indices_LL_;
  return(0);
}
#endif
//==============================================================================
int LinearProblem_CrsSingletonFilter::CreatePostSolveArrays(const Epetra_IntVector & RowIDs,
                                                     const Epetra_MapColoring & rowMapColors,
                                                     const Epetra_IntVector & ColProfiles,
                                                     const Epetra_IntVector & NewColProfiles,
                                                     const Epetra_IntVector & ColHasRowWithSingleton) {

  int j;

  if (NumMyColSingletons_==0) return(0); // Nothing to do

  Epetra_MapColoring & colMapColors = *ColMapColors_;

  int NumMyCols = FullMatrix()->NumMyCols();
  PRINT(outputRank_, "CreatePostSolveArrays() A    NumMyCols = " << NumMyCols);

  // We will need these arrays for the post-solve phase
  ColSingletonRowLIDs_ = new int[NumMyColSingletons_];
  ColSingletonColLIDs_ = new int[NumMyColSingletons_];
  ColSingletonPivotLIDs_ = new int[NumMyColSingletons_];
  ColSingletonPivots_ = new double[NumMyColSingletons_];

  // Register singleton columns (that were not already counted as singleton rows)
  // Check to see if any columns disappeared because all associated rows were eliminated
  int NumMyColSingletonstmp = 0;
  for (j=0; j<NumMyCols; j++) {
    int i = RowIDs[j];
    PRINT(outputRank_, "CreatePostSolveArrays() B    j, i = " <<j<< ", " <<i);
    PRINT(outputRank_, "CreatePostSolveArrays() B    ColProfiles["<<j<<"] = " << ColProfiles[j]);
    PRINT(outputRank_, "CreatePostSolveArrays() B    rowMapColors["<<i<<"] = " << rowMapColors[i]);

    if ( ColProfiles[j]==1 && rowMapColors[i]!=1 ) {
      ColSingletonRowLIDs_[NumMyColSingletonstmp] = i;
      ColSingletonColLIDs_[NumMyColSingletonstmp] = j;
      NumMyColSingletonstmp++;
    }
    // Also check for columns that were eliminated implicitly by
    // having all associated row eliminated
    else if (NewColProfiles[j]==0 && ColHasRowWithSingleton[j]!=1 && rowMapColors[i]==0) {
          colMapColors[j] = 1;
    }
  }

  for (j=0; j<NumMyColSingletons_; j++) {
    PRINT(outputRank_, "CreatePostSolveArrays() C    ColSingletonRowLIDs_["<<j<<"] = " << ColSingletonRowLIDs_[j]);
    PRINT(outputRank_, "CreatePostSolveArrays() C    ColSingletonColLIDs_["<<j<<"] = " << ColSingletonColLIDs_[j]);
  }
  auto size1 = FullMatrix()->NumMyRows();
  for (j=0; j<size1; j++) {
    if (colMapColors[j] != 0) {
      PRINT_NB(outputRank_, "CreatePostSolveArrays() C    colMapColors["<<j<<"] = " << colMapColors[j]);
    }
  }

  PRINT(outputRank_, "CreatePostSolveArrays() C    NumMyColSingletonstmp = " << NumMyColSingletonstmp);
  PRINT(outputRank_, "CreatePostSolveArrays() C    NumMyColSingletons_ = " << NumMyColSingletons_);

  assert(NumMyColSingletonstmp==NumMyColSingletons_); //Sanity check
  Epetra_Util sorter;
  sorter.Sort(true, NumMyColSingletons_, ColSingletonRowLIDs_, 0, 0, 1, &ColSingletonColLIDs_);

  for (j=0; j<NumMyColSingletons_; j++) {
    PRINT(outputRank_, "CreatePostSolveArrays() D    ColSingletonRowLIDs_["<<j<<"] = " << ColSingletonRowLIDs_[j]);
    PRINT(outputRank_, "CreatePostSolveArrays() D    ColSingletonColLIDs_["<<j<<"] = " << ColSingletonColLIDs_[j]);
  }

  return(0);
}

} //namespace EpetraExt

