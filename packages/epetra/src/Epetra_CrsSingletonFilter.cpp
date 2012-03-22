/*
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
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
// ************************************************************************
//@HEADER
*/

  // TODO this file needs to be changed for long long

#include "Epetra_Map.h"
#include "Epetra_Util.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Comm.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_MapColoring.h"
#include "Epetra_CrsSingletonFilter.h"

//==============================================================================
Epetra_CrsSingletonFilter::Epetra_CrsSingletonFilter() {
  InitializeDefaults();
}
//==============================================================================
Epetra_CrsSingletonFilter::~Epetra_CrsSingletonFilter(){


  if (ReducedProblem_!=0) delete ReducedProblem_;
  if (ReducedMatrix_!=0) delete ReducedMatrix_;
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
  if (RowMapColors_ != 0) delete RowMapColors_;
  if (ColMapColors_ != 0) delete ColMapColors_;

  if (ColSingletonRowLIDs_ != 0) delete [] ColSingletonRowLIDs_;
  if (ColSingletonColLIDs_ != 0) delete [] ColSingletonColLIDs_;
  if (ColSingletonPivotLIDs_ != 0) delete [] ColSingletonPivotLIDs_;
  if (ColSingletonPivots_ != 0) delete [] ColSingletonPivots_;
  if (tempExportX_ != 0) delete tempExportX_;
  if (Indices_ != 0) delete [] Indices_;
  if (tempX_ != 0) delete tempX_;
  if (tempB_ != 0) delete tempB_;

}
//==============================================================================
void Epetra_CrsSingletonFilter::InitializeDefaults() { 

// Initialize all attributes that have trivial default values

  FullProblem_ = 0;
  ReducedProblem_ = 0;
  FullMatrix_ = 0;
  ReducedMatrix_ = 0;
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

  Indices_ = 0;

  RowMapColors_ = 0;
  ColMapColors_ = 0;

  FullMatrixIsCrsMatrix_ = false;
  MaxNumMyEntries_ = 0;
  return;
}
//==============================================================================
int Epetra_CrsSingletonFilter::Analyze(Epetra_RowMatrix * fullMatrix) {

  int i, j, jj;

  FullMatrix_ = fullMatrix; 

  if (AnalysisDone_) EPETRA_CHK_ERR(-1); // Analysis already done once.  Cannot do it again
  if (fullMatrix==0) EPETRA_CHK_ERR(-2); // Input matrix pointer is zero
  if (fullMatrix->NumGlobalRows()==0) EPETRA_CHK_ERR(-3); // Full matrix has zero dimension.
  if (fullMatrix->NumGlobalNonzeros()==0) EPETRA_CHK_ERR(-4); // Full matrix has no nonzero terms.


  // First check for columns with single entries and find columns with singleton rows
  Epetra_IntVector ColProfiles(FullMatrixColMap()); ColProfiles.PutValue(0);
  Epetra_IntVector ColHasRowWithSingleton(FullMatrixColMap()); ColHasRowWithSingleton.PutValue(0);

  // RowIDs[j] will contain the local row ID associated with the jth column, 
  // if the jth col has a single entry
  Epetra_IntVector RowIDs(FullMatrixColMap()); RowIDs.PutValue(-1);

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
  int NumIndices;
  int * Indices;
  NumMyRowSingletons_ = 0;
  for (i=0; i<NumMyRows; i++) {
    // Get ith row
    EPETRA_CHK_ERR(GetRow(i, NumIndices, Indices));
    for (j=0; j<NumIndices; j++) {
      int ColumnIndex = Indices[j];
      ColProfiles[ColumnIndex]++; // Increment column count
      // Record local row ID for current column
      // will use to identify row to eliminate if column is a singleton
      RowIDs[ColumnIndex] = i;
    }
    // If row has single entry, color it and associated column with color=1
    if (NumIndices==1) {
      j = Indices[0];
      ColHasRowWithSingleton[j]++;
      rowMapColors[i] = 1;
      colMapColors[j] = 1;
      NumMyRowSingletons_++;
    }
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
    EPETRA_CHK_ERR(tmpVec.Export(ColProfiles, *fullMatrix->RowMatrixImporter(), Add));
    EPETRA_CHK_ERR(ColProfiles.Import(tmpVec, *fullMatrix->RowMatrixImporter(), Insert));
    
    EPETRA_CHK_ERR(tmpVec.PutValue(0));
    EPETRA_CHK_ERR(tmpVec.Export(ColHasRowWithSingleton, *fullMatrix->RowMatrixImporter(), Add));
    EPETRA_CHK_ERR(ColHasRowWithSingleton.Import(tmpVec, *fullMatrix->RowMatrixImporter(), Insert));
  }
  // ColProfiles now contains the nonzero column entry count for all columns that have
  // an entry on this processor.
  // ColHasRowWithSingleton now contains a count of singleton rows associated with the corresponding
  // local column.  Next we check to make sure no column is associated with more than one singleton row.

  if (ColHasRowWithSingleton.MaxValue()>1) {
    EPETRA_CHK_ERR(-2); // At least one col is associated with two singleton rows, can't handle it.
  }


  Epetra_IntVector RowHasColWithSingleton(fullMatrix->RowMatrixRowMap()); // Use to check for errors
  RowHasColWithSingleton.PutValue(0);
 
  NumMyColSingletons_ = 0;
  // Count singleton columns (that were not already counted as singleton rows)
  for (j=0; j<NumMyCols; j++) {
    i = RowIDs[j];
    // Check if column is a singleton
    if (ColProfiles[j]==1) {
      // Check to see if this column already eliminated by the row check above
      if (rowMapColors[i]!=1) {
	RowHasColWithSingleton[i]++; // Increment col singleton counter for ith row
	rowMapColors[i] = 2; // Use 2 for now, to distinguish between row eliminated directly or via column singletons
	colMapColors[j] = 1;
	NumMyColSingletons_++;
	// If we delete a row, we need to keep track of associated column entries that were also deleted 
	// in case all entries in a column are eventually deleted, in which case the column should
	// also be deleted.
	EPETRA_CHK_ERR(GetRow(i, NumIndices, Indices));
	for (jj=0; jj<NumIndices; jj++) NewColProfiles[Indices[jj]]--;
	
      }
    }
    // Check if some other processor eliminated this column    
    else if (ColHasRowWithSingleton[j]==1 && rowMapColors[i]!=1) { 
	colMapColors[j] = 1;
    }
  }
  if (RowHasColWithSingleton.MaxValue()>1) {
    EPETRA_CHK_ERR(-3); // At lease one row is associated with two singleton cols, can't handle it.
  }


 // Generate arrays that keep track of column singleton row, col and pivot info needed for post-solve phase
  EPETRA_CHK_ERR(CreatePostSolveArrays(RowIDs, rowMapColors, ColProfiles, NewColProfiles,
				       ColHasRowWithSingleton));

  for (i=0; i<NumMyRows; i++) if (rowMapColors[i]==2) rowMapColors[i] = 1; // Convert all eliminated rows to same color

  fullMatrix->RowMatrixRowMap().Comm().SumAll(&NumMyRowSingletons_, &NumGlobalRowSingletons_, 1);
  fullMatrix->RowMatrixRowMap().Comm().SumAll(&NumMyColSingletons_, &NumGlobalColSingletons_, 1);
  AnalysisDone_ = true;
  return(0);
}
//==============================================================================
int Epetra_CrsSingletonFilter::ConstructReducedProblem(Epetra_LinearProblem * Problem) {

  int i, j;
  if (HaveReducedProblem_) EPETRA_CHK_ERR(-1); // Setup already done once.  Cannot do it again
  if (Problem==0) EPETRA_CHK_ERR(-2); // Null problem pointer

  FullProblem_ = Problem;
  FullMatrix_ = dynamic_cast<Epetra_RowMatrix *>(Problem->GetMatrix());
  if (FullMatrix_==0) EPETRA_CHK_ERR(-3); // Need a RowMatrix
  if (Problem->GetRHS()==0) EPETRA_CHK_ERR(-4); // Need a RHS
  if (Problem->GetLHS()==0) EPETRA_CHK_ERR(-5); // Need a LHS
  // Generate reduced row and column maps

  Epetra_MapColoring & rowMapColors = *RowMapColors_;
  Epetra_MapColoring & colMapColors = *ColMapColors_;

  ReducedMatrixRowMap_ = rowMapColors.GenerateMap(0);
  ReducedMatrixColMap_ = colMapColors.GenerateMap(0);

  // Create domain and range map colorings by exporting map coloring of column and row maps

  if (FullMatrix()->RowMatrixImporter()!=0) {
    Epetra_MapColoring DomainMapColors(FullMatrixDomainMap());
    EPETRA_CHK_ERR(DomainMapColors.Export(*ColMapColors_, *FullMatrix()->RowMatrixImporter(), AbsMax));
    OrigReducedMatrixDomainMap_ = DomainMapColors.GenerateMap(0);
  }
  else
    OrigReducedMatrixDomainMap_ = ReducedMatrixColMap_;

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

  // Create importers
  Full2ReducedLHSImporter_ = new Epetra_Import(*ReducedMatrixDomainMap(), FullMatrixDomainMap());
  Full2ReducedRHSImporter_ = new Epetra_Import(*ReducedMatrixRowMap(), FullRHS->Map());

  // Construct Reduced Matrix
  ReducedMatrix_ = new Epetra_CrsMatrix(Copy, *ReducedMatrixRowMap(), *ReducedMatrixColMap(), 0);

  // Create storage for temporary X values due to explicit elimination of rows
  tempExportX_ = new Epetra_MultiVector(FullMatrixColMap(), NumVectors);

  int NumEntries;
  int * Indices;
  double * Values;
  int NumMyRows = FullMatrix()->NumMyRows();
  int ColSingletonCounter = 0;
  for (i=0; i<NumMyRows; i++) {
    int curGRID = FullMatrixRowMap().GID(i);
    if (ReducedMatrixRowMap()->MyGID(curGRID)) { // Check if this row should go into reduced matrix

      EPETRA_CHK_ERR(GetRowGCIDs(i, NumEntries, Values, Indices)); // Get current row (Indices are global)
      
      int ierr = ReducedMatrix()->InsertGlobalValues(curGRID, NumEntries, 
						     Values, Indices); // Insert into reduce matrix
      // Positive errors will occur because we are submitting col entries that are not part of
      // reduced system.  However, because we specified a column map to the ReducedMatrix constructor
      // these extra column entries will be ignored and we will be politely reminded by a positive
      // error code
      if (ierr<0) EPETRA_CHK_ERR(ierr); 
    }
    else {
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

  // Now convert to local indexing.  We have constructed things so that the domain and range of the
  // matrix will have the same map.  If the reduced matrix domain and range maps were not the same, the
  // differences were addressed in the ConstructRedistributeExporter() method
  EPETRA_CHK_ERR(ReducedMatrix()->FillComplete(*ReducedMatrixDomainMap(), *ReducedMatrixRangeMap()));

  // Construct Reduced LHS (Puts any initial guess values into reduced system)

  ReducedLHS_ = new Epetra_MultiVector(*ReducedMatrixDomainMap(), NumVectors);
  EPETRA_CHK_ERR(ReducedLHS_->Import(*FullLHS, *Full2ReducedLHSImporter_, Insert));
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
  }


  EPETRA_CHK_ERR(FullMatrix()->Multiply(false, *tempX_, *tempB_));

  EPETRA_CHK_ERR(tempB_->Update(1.0, *FullRHS, -1.0)); // tempB now has influence of already-known X values

  ReducedRHS_ = new Epetra_MultiVector(*ReducedMatrixRowMap(), FullRHS->NumVectors());
  EPETRA_CHK_ERR(ReducedRHS_->Import(*tempB_, *Full2ReducedRHSImporter_, Insert));

  // Finally construct Reduced Linear Problem

  ReducedProblem_ = new Epetra_LinearProblem(ReducedMatrix_, ReducedLHS_, ReducedRHS_);

  double fn = FullMatrix()->NumGlobalRows();
  double fnnz = FullMatrix()->NumGlobalNonzeros();
  double rn = ReducedMatrix()->NumGlobalRows();
  double rnnz = ReducedMatrix()->NumGlobalNonzeros();

  RatioOfDimensions_ = (fn-rn)/fn;
  RatioOfNonzeros_ = (fnnz-rnnz)/fnnz;
  HaveReducedProblem_ = true;
  
  return(0);
}
//==============================================================================
int Epetra_CrsSingletonFilter::UpdateReducedProblem(Epetra_LinearProblem * Problem) {

  int i, j;

  if (Problem==0) EPETRA_CHK_ERR(-1); // Null problem pointer

  FullProblem_ = Problem;
  FullMatrix_ = dynamic_cast<Epetra_RowMatrix *>(Problem->GetMatrix());
  if (FullMatrix_==0) EPETRA_CHK_ERR(-2); // Need a RowMatrix
  if (Problem->GetRHS()==0) EPETRA_CHK_ERR(-3); // Need a RHS
  if (Problem->GetLHS()==0) EPETRA_CHK_ERR(-4); // Need a LHS
  if (!HaveReducedProblem_) EPETRA_CHK_ERR(-5); // Must have set up reduced problem

  // Create pointer to Full RHS, LHS
  Epetra_MultiVector * FullRHS = FullProblem()->GetRHS();
  Epetra_MultiVector * FullLHS = FullProblem()->GetLHS();
  int NumVectors = FullLHS->NumVectors();

  int NumEntries;
  int * Indices;
  double * Values;
  int NumMyRows = FullMatrix()->NumMyRows();
  int ColSingletonCounter = 0;
  for (i=0; i<NumMyRows; i++) {
    int curGRID = FullMatrixRowMap().GID(i);
    if (ReducedMatrixRowMap()->MyGID(curGRID)) { // Check if this row should go into reduced matrix
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
    return(0);
}
//==============================================================================
int Epetra_CrsSingletonFilter::ConstructRedistributeExporter(Epetra_Map * SourceMap, Epetra_Map * TargetMap,
							     Epetra_Export * & RedistributeExporter,
							     Epetra_Map * & RedistributeMap) {

  int IndexBase = SourceMap->IndexBase();
  if (IndexBase!=TargetMap->IndexBase()) EPETRA_CHK_ERR(-1);

  const Epetra_Comm & Comm = TargetMap->Comm();

  int TargetNumMyElements = TargetMap->NumMyElements();
  int SourceNumMyElements = SourceMap->NumMyElements();

  // ContiguousTargetMap has same number of elements per PE as TargetMap, but uses contigious indexing 
  Epetra_Map ContiguousTargetMap(-1, TargetNumMyElements, IndexBase,Comm);

  // Same for ContiguousSourceMap
  Epetra_Map ContiguousSourceMap(-1, SourceNumMyElements, IndexBase, Comm);

  assert(ContiguousSourceMap.NumGlobalElements()==ContiguousTargetMap.NumGlobalElements());

  // Now create a vector that contains the global indices of the Source Epetra_MultiVector
  Epetra_IntVector *SourceIndices;
  Epetra_LongLongVector *SourceIndices_LL;
  if(SourceMap->GlobalIndicesInt()) 
    SourceIndices = new Epetra_IntVector(View, ContiguousSourceMap, SourceMap->MyGlobalElements());
  else
    SourceIndices_LL = new Epetra_LongLongVector(View, ContiguousSourceMap, SourceMap->MyGlobalElements());

  // Create an exporter to send the SourceMap global IDs to the target distribution
  Epetra_Export Exporter(ContiguousSourceMap, ContiguousTargetMap);
  
  // Create a vector to catch the global IDs in the target distribution
  Epetra_IntVector *TargetIndices;
  Epetra_LongLongVector *TargetIndices_LL;  
  if(TargetMap->GlobalIndicesInt()) {
    TargetIndices = new Epetra_IntVector(ContiguousTargetMap);
    TargetIndices->Export(*SourceIndices, Exporter, Insert);
  } else {
    TargetIndices_LL = new Epetra_IntVector(ContiguousTargetMap);
    TargetIndices_LL->Export(*SourceIndices_LL, Exporter, Insert);
  }

  // Create a new map that describes how the Source MultiVector should be laid out so that it has
  // the same number of elements on each processor as the TargetMap
  if(TargetMap->GlobalIndicesInt()) 
    RedistributeMap = new Epetra_Map(-1, TargetNumMyElements, TargetIndices->Values(), IndexBase, Comm);
  else
    RedistributeMap = new Epetra_Map(-1, TargetNumMyElements, TargetIndices_LL->Values(), IndexBase, Comm);

  // This exporter will finally redistribute the Source MultiVector to the same layout as the TargetMap
  RedistributeExporter = new Epetra_Export(*SourceMap, *RedistributeMap);
  return(0);
}
//==============================================================================
int Epetra_CrsSingletonFilter::ComputeFullSolution() {

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
    
  return(0);
}
//==============================================================================
int Epetra_CrsSingletonFilter::InitFullMatrixAccess() {

  MaxNumMyEntries_ = FullMatrix()->MaxNumEntries(); 

  // Cast to CrsMatrix, if possible.  Can save some work.
  FullCrsMatrix_ = dynamic_cast<Epetra_CrsMatrix *>(FullMatrix());
  FullMatrixIsCrsMatrix_ = (FullCrsMatrix_!=0); // Pointer is non-zero if cast worked
  Indices_ = new int[MaxNumMyEntries_];
  Values_.Size(MaxNumMyEntries_);

  return(0);
}
//==============================================================================
int Epetra_CrsSingletonFilter::GetRow(int Row, int & NumIndices, int * & Indices) {

  if (FullMatrixIsCrsMatrix_) { // View of current row
    EPETRA_CHK_ERR(FullCrsMatrix()->Graph().ExtractMyRowView(Row, NumIndices, Indices)); 
  }
  else { // Copy of current row (we must get the values, but we ignore them)
    EPETRA_CHK_ERR(FullMatrix()->ExtractMyRowCopy(Row, MaxNumMyEntries_, NumIndices, 
						  Values_.Values(), Indices_));
    Indices = Indices_;
  } 
  return(0);
}
//==============================================================================
int Epetra_CrsSingletonFilter::GetRow(int Row, int & NumIndices, 
				      double * & Values, int * & Indices) {

  if (FullMatrixIsCrsMatrix_) { // View of current row
    EPETRA_CHK_ERR(FullCrsMatrix_->ExtractMyRowView(Row, NumIndices, Values, Indices)); 
  }
  else { // Copy of current row (we must get the values, but we ignore them)
    EPETRA_CHK_ERR(FullMatrix()->ExtractMyRowCopy(Row, MaxNumMyEntries_, NumIndices, 
						  Values_.Values(), Indices_));
    Values = Values_.Values();
    Indices = Indices_;
  } 
  return(0);
}
//==============================================================================
int Epetra_CrsSingletonFilter::GetRowGCIDs(int Row, int & NumIndices, 
					   double * & Values, int * & GlobalIndices) {

    EPETRA_CHK_ERR(FullMatrix()->ExtractMyRowCopy(Row, MaxNumMyEntries_, NumIndices, 
						  Values_.Values(), Indices_));
    for (int j=0; j<NumIndices; j++) Indices_[j] = FullMatrixColMap().GID(Indices_[j]);
    Values = Values_.Values();
    GlobalIndices = Indices_;
  return(0);
}
//==============================================================================
int Epetra_CrsSingletonFilter::CreatePostSolveArrays(const Epetra_IntVector & RowIDs,
						     const Epetra_MapColoring & rowMapColors,
						     const Epetra_IntVector & ColProfiles,
						     const Epetra_IntVector & NewColProfiles,
						     const Epetra_IntVector & ColHasRowWithSingleton) {

  int j;

  if (NumMyColSingletons_==0) return(0); // Nothing to do

  Epetra_MapColoring & colMapColors = *ColMapColors_;

  int NumMyCols = FullMatrix()->NumMyCols();

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
    if ( ColProfiles[j]==1 && rowMapColors[i]!=1) {
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

  assert(NumMyColSingletonstmp==NumMyColSingletons_); //Sanity check
  Epetra_Util sorter;
  sorter.Sort(true, NumMyColSingletons_, ColSingletonRowLIDs_, 0, 0, 1, &ColSingletonColLIDs_, 0, 0);
    
  return(0);
}
