
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

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
#include "Komplex_LinearProblem.h"

//==============================================================================
Komplex_LinearProblem::Komplex_LinearProblem() {
  InitializeDefaults();
}
//==============================================================================
Komplex_LinearProblem::~Komplex_LinearProblem(){


  if (KomplexProblem_!=0) delete KomplexProblem_;
  if (KomplexMatrix_!=0) delete KomplexMatrix_;
  if (KomplexLHS_!=0) delete KomplexLHS_;
  if (KomplexRHS_!=0) delete KomplexRHS_;
  if (KomplexMatrixDomainMap_!=KomplexMatrixColMap_) delete KomplexMatrixDomainMap_;
  if (KomplexMatrixRangeMap_!=KomplexMatrixRowMap_) delete KomplexMatrixRangeMap_;
  if (KomplexMatrixRowMap_!=0) delete KomplexMatrixRowMap_;
  if (KomplexMatrixColMap_!=0) delete KomplexMatrixColMap_;

}
//==============================================================================
void Komplex_LinearProblem::InitializeDefaults() { 

// Initialize all attributes that have trivial default values

  KomplexProblem_ = 0;
  FullMatrix_ = 0;
  KomplexMatrix_ = 0;
  KomplexRHS_ = 0;
  KomplexLHS_ = 0;
  KomplexMatrixRowMap_ = 0;
  KomplexMatrixColMap_ = 0;
  KomplexMatrixDomainMap_ = 0;
  KomplexMatrixRangeMap_ = 0;

  HaveKomplexProblem_ = false;
  AnalysisDone_ = false;

  Values_ = 0;
  Indices_ = 0;


  UserMatrixIsCrsMatrix_ = false;
  UserMatrixIsVbrMatrix_ = false;
  MaxNumMyEntries_ = 0;
  return;
}
//==============================================================================
int Komplex_LinearProblem::Analyze(Epetra_RowMatrix * FullMatrix) {

  int i, j, jj;

  FullMatrix_ = FullMatrix; 

  if (AnalysisDone_) EPETRA_CHK_ERR(-1); // Analysis already done once.  Cannot do it again

  // First check for columns with single entries and find columns with singleton rows
  Epetra_IntVector ColProfiles(FullMatrixColMap()); ColProfiles.PutValue(0);
  Epetra_IntVector ColHasRowWithSingleton(FullMatrixColMap()); ColHasRowWithSingleton.PutValue(0);

  // Define MapColoring objects
  RowMapColors_ = new Epetra_MapColoring(FullMatrixRowMap());  // Initial colors are all 0
  ColMapColors_ = new Epetra_MapColoring(FullMatrixColMap());
  Epetra_MapColoring & RowMapColors = *RowMapColors_;
  Epetra_MapColoring & ColMapColors = *ColMapColors_;


  int NumMyRows = FullMatrix->NumMyRows();
  int NumMyCols = FullMatrix->NumMyCols();

  // RowIDs[j] will contain the local row ID associated with the jth column, 
  // if the jth col has a single entry
  int * RowIDs = 0;
  if (NumMyCols>0) RowIDs = new int[NumMyCols];
  for (i=0; i<NumMyCols; i++) RowIDs[i] = -1;

  // Set up for accessing full matrix.  Will do so row-by-row.
  EPETRA_CHK_ERR(InitFullMatrixAccess());

  // Scan matrix for singleton rows, build up column profiles
  int NumIndices;
  int * Indices;
  NumRowSingletons_ = 0;
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
      int j = Indices[0];
      //cout << "i, j = " << i << "  " << j << endl;
      ColHasRowWithSingleton[j]++;
      RowMapColors[i] = 1;
      ColMapColors[j] = 1;
      NumRowSingletons_++;
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

  if (FullMatrix->RowMatrixImporter()!=0) {
    Epetra_IntVector tmpVec(FullMatrixDomainMap()); // Use for gather/scatter of column vectors
    EPETRA_CHK_ERR(tmpVec.Export(ColProfiles, *FullMatrix->RowMatrixImporter(), Add));
    EPETRA_CHK_ERR(ColProfiles.Import(tmpVec, *FullMatrix->RowMatrixImporter(), Insert));
    
    EPETRA_CHK_ERR(tmpVec.PutValue(0));
    EPETRA_CHK_ERR(tmpVec.Export(ColHasRowWithSingleton, *FullMatrix->RowMatrixImporter(), Add));
    EPETRA_CHK_ERR(ColHasRowWithSingleton.Import(tmpVec, *FullMatrix->RowMatrixImporter(), Insert));
  }
  // ColProfiles now contains the nonzero column entry count for all columns that have
  // an entry on this processor.
  // ColHasRowWithSingleton now contains a count of singleton rows associated with the corresponding
  // local column.  Next we check to make sure no column is associated with more than one singleton row.

  if (ColHasRowWithSingleton.MaxValue()>1) {
    EPETRA_CHK_ERR(-2); // At least one col is associated with two singleton rows, can't handle it.
  }


  Epetra_IntVector RowHasColWithSingleton(FullMatrix->RowMatrixRowMap()); // Use to check for errors
  RowHasColWithSingleton.PutValue(0);
 
  NumColSingletons_ = 0;
  // Count singleton columns (that were not already counted as singleton rows)
  for (j=0; j<NumMyCols; j++) {
    int i = RowIDs[j];
    // Check if column is a singleton
    if (ColProfiles[j]==1) {
      // Check to see if this column already eliminated by the row check above
      if (RowMapColors[i]!=1) {
	RowHasColWithSingleton[i]++; // Increment col singleton counter for ith row
	RowMapColors[i] = 2; // Use 2 for now, to distinguish between row eliminated directly or via column singletons
	ColMapColors[j] = 1;
	NumColSingletons_++;
	// If we delete a row, we need to keep track of associated column entries that were also deleted 
	// in case all entries in a column are eventually deleted, in which case the column should
	// also be deleted.
	EPETRA_CHK_ERR(GetRow(i, NumIndices, Indices));
	for (jj=0; jj<NumIndices; jj++) NewColProfiles[Indices[jj]]--;
	
      }
    }
    // Check if some other processor eliminated this column    
    else if (ColHasRowWithSingleton[j]==1 && RowMapColors[i]!=1) { 
	ColMapColors[j] = 1;
    }
  }
  if (RowHasColWithSingleton.MaxValue()>1) {
    EPETRA_CHK_ERR(-3); // At lease one row is associated with two singleton cols, can't handle it.
  }

 // Generate arrays that keep track of column singleton row, col and pivot info needed for post-solve phase
  EPETRA_CHK_ERR(CreatePostSolveArrays(RowIDs, RowMapColors, ColProfiles, NewColProfiles,
				       ColHasRowWithSingleton));

  for (i=0; i<NumMyRows; i++) if (RowMapColors[i]==2) RowMapColors[i] = 1; // Convert all eliminated rows to same color

  if (RowIDs!=0) delete [] RowIDs;

  AnalysisDone_ = true;
  return(0);
}
//==============================================================================
int Komplex_LinearProblem::ConstructKomplexProblem(Epetra_LinearProblem * Problem) {

  int i, j;
  if (HaveKomplexProblem_) EPETRA_CHK_ERR(-1); // Setup already done once.  Cannot do it again
  if (Problem==0) EPETRA_CHK_ERR(-2); // Null problem pointer

  FullProblem_ = Problem;
  FullMatrix_ = dynamic_cast<Epetra_RowMatrix *>(Problem->GetMatrix());
  if (FullMatrix_==0) EPETRA_CHK_ERR(-3); // Need a RowMatrix
  if (Problem->GetRHS()==0) EPETRA_CHK_ERR(-4); // Need a RHS
  if (Problem->GetLHS()==0) EPETRA_CHK_ERR(-5); // Need a LHS
  // Generate reduced row and column maps

  Epetra_MapColoring & RowMapColors = *RowMapColors_;
  Epetra_MapColoring & ColMapColors = *ColMapColors_;

  KomplexMatrixRowMap_ = RowMapColors.GenerateMap(0);
  KomplexMatrixColMap_ = ColMapColors.GenerateMap(0);

  // Create domain and range map colorings by exporting map coloring of column and row maps

  if (FullMatrix()->RowMatrixImporter()!=0) {
    Epetra_MapColoring DomainMapColors(FullMatrixDomainMap());
    EPETRA_CHK_ERR(DomainMapColors.Export(*ColMapColors_, *FullMatrix()->RowMatrixImporter(), AbsMax));
    OrigKomplexMatrixDomainMap_ = DomainMapColors.GenerateMap(0);
  }
  else
    OrigKomplexMatrixDomainMap_ = KomplexMatrixColMap_;

  if (FullMatrixIsCrsMatrix_) {
    if (FullCrsMatrix()->Exporter()!=0) { // Non-trivial exporter
      Epetra_MapColoring RangeMapColors(FullMatrixRangeMap());
      EPETRA_CHK_ERR(RangeMapColors.Export(*RowMapColors_, *FullCrsMatrix()->Exporter(), 
					   AbsMax));
      KomplexMatrixRangeMap_ = RangeMapColors.GenerateMap(0);
    }
    else
      KomplexMatrixRangeMap_ = KomplexMatrixRowMap_;
  }
  else
    KomplexMatrixRangeMap_ = KomplexMatrixRowMap_;

  // Check to see if the reduced system domain and range maps are the same.
  // If not, we need to remap entries of the LHS multivector so that they are distributed
  // conformally with the rows of the reduced matrix and the RHS multivector
  SymmetricElimination_ = KomplexMatrixRangeMap_->SameAs(*OrigKomplexMatrixDomainMap_);
  if (!SymmetricElimination_) 
    ConstructRedistributeExporter(KomplexMatrixRangeMap_, OrigKomplexMatrixDomainMap_, 
				  RedistributeDomainExporter_, KomplexMatrixDomainMap_);
  else {
    KomplexMatrixDomainMap_ = OrigKomplexMatrixDomainMap_;
    OrigKomplexMatrixDomainMap_ = 0;
  }
  
  // Create pointer to Full RHS, LHS
  Epetra_MultiVector * FullRHS = FullProblem()->GetRHS();
  Epetra_MultiVector * FullLHS = FullProblem()->GetLHS();
  int NumVectors = FullLHS->NumVectors();

  // Create importers
  Full2KomplexLHSImporter_ = new Epetra_Import(*KomplexMatrixDomainMap(), FullMatrixDomainMap());
  Full2KomplexRHSImporter_ = new Epetra_Import(*KomplexMatrixRowMap(), FullRHS->Map());

  // Construct Komplex Matrix
  KomplexMatrix_ = new Epetra_CrsMatrix(Copy, *KomplexMatrixRowMap(), *KomplexMatrixColMap(), 0);

  // Create storage for temporary X values due to explicit elimination of rows
  tempExportX_ = new Epetra_MultiVector(FullMatrixColMap(), NumVectors);

  int NumEntries;
  int * Indices;
  double * Values;
  int NumMyRows = FullMatrix()->NumMyRows();
  int ColSingletonCounter = 0;
  for (i=0; i<NumMyRows; i++) {
    int curGRID = FullMatrixRowMap().GID(i);
    if (KomplexMatrixRowMap()->MyGID(curGRID)) { // Check if this row should go into reduced matrix

      EPETRA_CHK_ERR(GetRowGCIDs(i, NumEntries, Values, Indices)); // Get current row (Indices are global)
      
      int ierr = KomplexMatrix()->InsertGlobalValues(curGRID, NumEntries, 
						     Values, Indices); // Insert into reduce matrix
      // Positive errors will occur because we are submitting col entries that are not part of
      // reduced system.  However, because we specified a column map to the KomplexMatrix constructor
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
	//if (i!=ColSingletonRowLIDs_[ColSingletonCounter]) cout << "i = "<<i<<" ColSingletonRowLIDs_["<<ColSingletonCounter<<"] = "
	//						       << ColSingletonRowLIDs_[ColSingletonCounter]<< endl;
	//assert(i==ColSingletonRowLIDs_[ColSingletonCounter]);  // Sanity test
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

  // Now convert to local indexing.
  EPETRA_CHK_ERR(KomplexMatrix()->TransformToLocal(KomplexMatrixDomainMap(), KomplexMatrixRangeMap()));

  // Construct Komplex LHS (Puts any initial guess values into reduced system)

  KomplexLHS_ = new Epetra_MultiVector(*KomplexMatrixDomainMap(), NumVectors);
  EPETRA_CHK_ERR(KomplexLHS_->Import(*FullLHS, *Full2KomplexLHSImporter_, Insert));
  FullLHS->PutScalar(0.0); // zero out Full LHS since we will inject values as we get them

  // Construct Komplex RHS

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

  //cout << "tempExportX = " << endl << *tempExportX << endl;
  //cout << "tempX = " << endl << *tempX << endl;
  //cout << "tempB = " << endl << *tempB << endl;

  EPETRA_CHK_ERR(FullMatrix()->Multiply(false, *tempX_, *tempB_));

  EPETRA_CHK_ERR(tempB_->Update(1.0, *FullRHS, -1.0)); // tempB now has influence of already-known X values

  KomplexRHS_ = new Epetra_MultiVector(*KomplexMatrixRowMap(), FullRHS->NumVectors());
  EPETRA_CHK_ERR(KomplexRHS_->Import(*tempB_, *Full2KomplexRHSImporter_, Insert));

  // Finally construct Komplex Linear Problem

  KomplexProblem_ = new Epetra_LinearProblem(KomplexMatrix_, KomplexLHS_, KomplexRHS_);

  HaveKomplexProblem_ = true;
  
  return(0);
}
//==============================================================================
int Komplex_LinearProblem::UpdateKomplexProblem(Epetra_LinearProblem * Problem) {

  int i, j;

  if (Problem==0) EPETRA_CHK_ERR(-1); // Null problem pointer

  FullProblem_ = Problem;
  FullMatrix_ = dynamic_cast<Epetra_RowMatrix *>(Problem->GetMatrix());
  if (FullMatrix_==0) EPETRA_CHK_ERR(-2); // Need a RowMatrix
  if (Problem->GetRHS()==0) EPETRA_CHK_ERR(-3); // Need a RHS
  if (Problem->GetLHS()==0) EPETRA_CHK_ERR(-4); // Need a LHS
  if (!HaveKomplexProblem_) EPETRA_CHK_ERR(-5); // Must have set up reduced problem

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
    if (KomplexMatrixRowMap()->MyGID(curGRID)) { // Check if this row should go into reduced matrix
      EPETRA_CHK_ERR(GetRowGCIDs(i, NumEntries, Values, Indices)); // Get current row (indices global)
      int ierr = KomplexMatrix()->ReplaceGlobalValues(curGRID, NumEntries, 
						      Values, Indices);
      // Positive errors will occur because we are submitting col entries that are not part of
      // reduced system.  However, because we specified a column map to the KomplexMatrix constructor
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
	break;
      }
    }
  }

  assert(ColSingletonCounter==NumColSingletons_); // Sanity test

  // Update Komplex LHS (Puts any initial guess values into reduced system)

  KomplexLHS_->PutScalar(0.0); // zero out Komplex LHS
  EPETRA_CHK_ERR(KomplexLHS_->Import(*FullLHS, *Full2KomplexLHSImporter_, Insert));
  FullLHS->PutScalar(0.0); // zero out Full LHS since we will inject values as we get them

  // Construct Komplex RHS

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

  //cout << "tempExportX = " << endl << *tempExportX << endl;
  //cout << "tempX = " << endl << *tempX << endl;
  //cout << "tempB = " << endl << *tempB << endl;

  EPETRA_CHK_ERR(FullMatrix()->Multiply(false, *tempX_, *tempB_));

  EPETRA_CHK_ERR(tempB_->Update(1.0, *FullRHS, -1.0)); // tempB now has influence of already-known X values

  KomplexRHS_->PutScalar(0.0);
  EPETRA_CHK_ERR(KomplexRHS_->Import(*tempB_, *Full2KomplexRHSImporter_, Insert));
    return(0);
}
//==============================================================================
int Komplex_LinearProblem::Statistics() const {

  if (!HaveKomplexProblem_) {
    cout << "ConstructKomplexProblem method must be called first." << endl;
    return(0);
  }

  double fn = FullMatrix()->NumGlobalRows();
  double fnnz = FullMatrix()->NumGlobalNonzeros();
  double rn = KomplexMatrix()->NumGlobalRows();
  double rnnz = KomplexMatrix()->NumGlobalNonzeros();
  if (fn==0.0 || fnnz==0.0) {
    cout << "Full problem has zero size." << endl;
    return(0);
  }

  cout << "Full System characteristics:" << endl << endl
       << "  Dimension                             = " << FullMatrix()->NumGlobalRows() << endl
       << "  Number of nonzeros                    = " << FullMatrix()->NumGlobalNonzeros() << endl
       << "  Maximum Number of Row Entries         = " << MaxNumMyEntries_ << endl << endl
       << "Komplex System characteristics:" << endl << endl
       << "  Dimension                             = " << KomplexMatrix()->NumGlobalRows() << endl
       << "  Number of nonzeros                    = " << KomplexMatrix()->NumGlobalNonzeros() << endl
       << "  Maximum Number of Row Entries         = " << KomplexMatrix()->GlobalMaxNumEntries() << endl << endl
       << "Singleton information: " << endl
       << "  Number of rows with single entries    = " << NumRowSingletons_ << endl
       << "  Number of columns with single entries " << endl
       << "    (that were not already counted as " << endl
       << "     row singletons)                    = " << NumColSingletons_ << endl << endl
       << "Ratios: " << endl
       << "  Percent reduction in dimension        = " << (fn-rn)/fn*100.0 << endl
       << "  Percent reduction in nonzero count    = " << (fnnz-rnnz)/fnnz*100.0 << endl << endl;

  return(0);

}
//==============================================================================
int Komplex_LinearProblem::ConstructRedistributeExporter(Epetra_Map * SourceMap, Epetra_Map * TargetMap,
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

  // Now create a vector that contains the global indices of the Source Epetra_MultiVector
  Epetra_IntVector SourceIndices(View, ContiguousSourceMap, SourceMap->MyGlobalElements());

  cout << "Contiguous Source Map " << ContiguousSourceMap << endl
       << "Contiguous Target Map " << ContiguousTargetMap << endl;
  // Create an exporter to send the SourceMap global IDs to the target distribution
  Epetra_Export Exporter(ContiguousSourceMap, ContiguousTargetMap);
  
  // Create a vector to catch the global IDs in the target distribution
  Epetra_IntVector TargetIndices(ContiguousTargetMap);
  TargetIndices.Export(SourceIndices, Exporter, Insert);

  // Create a new map that describes how the Source MultiVector should be laid out so that it has
  // the same number of elements on each processor as the TargetMap
  RedistributeMap = new Epetra_Map(-1, TargetNumMyElements, TargetIndices.Values(), IndexBase, Comm);

  // This exporter will finally redistribute the Source MultiVector to the same layout as the TargetMap
  RedistributeExporter = new Epetra_Export(*SourceMap, *RedistributeMap);
  return(0);
}
//==============================================================================
int Komplex_LinearProblem::Analyze(int AbsoluteThreshold, double RelativeThreshold) {
  EPETRA_CHK_ERR(-1); // Not implemented
  return(0);
}
//==============================================================================
int Komplex_LinearProblem::ComputeFullSolution() {

  int jj, k;

  Epetra_MultiVector * FullLHS = FullProblem()->GetLHS(); 
  Epetra_MultiVector * FullRHS = FullProblem()->GetRHS(); 

  tempX_->PutScalar(0.0); tempExportX_->PutScalar(0.0);
  // Inject values that the user computed for the reduced problem into the full solution vector
  EPETRA_CHK_ERR(tempX_->Export(*KomplexLHS_, *Full2KomplexLHSImporter_, Add));
  FullLHS->Update(1.0, *tempX_, 1.0);

  // Next we will use our full solution vector which is populated with pre-filter solution
  // values and reduced system solution values to compute the sum of the row contributions
  // that must be subtracted to get the post-filter solution values

  EPETRA_CHK_ERR(FullMatrix()->Multiply(false, *FullLHS, *tempB_));



  // Finally we loop through the local rows that were associated with column singletons and compute the
  // solution for these equations.

  int NumVectors = tempB_->NumVectors();
  for (k=0; k<NumColSingletons_; k++) {
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
int Komplex_LinearProblem::InitFullMatrixAccess() {

  MaxNumMyEntries_ = FullMatrix()->MaxNumEntries(); 

  // Cast to CrsMatrix, if possible.  Can save some work.
  FullCrsMatrix_ = dynamic_cast<Epetra_CrsMatrix *>(FullMatrix());
  FullMatrixIsCrsMatrix_ = (FullCrsMatrix_!=0); // Pointer is non-zero if cast worked
  Indices_ = new int[MaxNumMyEntries_];
  Values_ = new double[MaxNumMyEntries_];

  return(0);
}
//==============================================================================
int Komplex_LinearProblem::GetRow(int Row, int & NumIndices, int * & Indices) {

  if (FullMatrixIsCrsMatrix_) { // View of current row
    EPETRA_CHK_ERR(FullCrsMatrix()->Graph().ExtractMyRowView(Row, NumIndices, Indices)); 
  }
  else { // Copy of current row (we must get the values, but we ignore them)
    EPETRA_CHK_ERR(FullMatrix()->ExtractMyRowCopy(Row, MaxNumMyEntries_, NumIndices, 
						  Values_, Indices_));
    Indices = Indices_;
  } 
  return(0);
}
//==============================================================================
int Komplex_LinearProblem::GetRow(int Row, int & NumIndices, 
				      double * & Values, int * & Indices) {

  if (FullMatrixIsCrsMatrix_) { // View of current row
    EPETRA_CHK_ERR(FullCrsMatrix_->ExtractMyRowView(Row, NumIndices, Values, Indices)); 
  }
  else { // Copy of current row (we must get the values, but we ignore them)
    EPETRA_CHK_ERR(FullMatrix()->ExtractMyRowCopy(Row, MaxNumMyEntries_, NumIndices, 
						  Values_, Indices_));
    Values = Values_;
    Indices = Indices_;
  } 
  return(0);
}
//==============================================================================
int Komplex_LinearProblem::GetRowGCIDs(int Row, int & NumIndices, 
					   double * & Values, int * & GlobalIndices) {

    EPETRA_CHK_ERR(FullMatrix()->ExtractMyRowCopy(Row, MaxNumMyEntries_, NumIndices, 
						  Values_, Indices_));
    for (int j=0; j<NumIndices; j++) Indices_[j] = FullMatrixColMap().GID(Indices_[j]);
    Values = Values_;
    GlobalIndices = Indices_;
  return(0);
}
//==============================================================================
int Komplex_LinearProblem::CreatePostSolveArrays(int * RowIDs,
						     const Epetra_MapColoring & RowMapColors,
						     const Epetra_IntVector & ColProfiles,
						     const Epetra_IntVector & NewColProfiles,
						     const Epetra_IntVector & ColHasRowWithSingleton) {

  int j;

  if (NumColSingletons_==0) return(0); // Nothing to do

  Epetra_MapColoring & ColMapColors = *ColMapColors_;

  int NumMyCols = FullMatrix()->NumMyCols();

  // We will need these arrays for the post-solve phase
  ColSingletonRowLIDs_ = new int[NumColSingletons_];
  ColSingletonColLIDs_ = new int[NumColSingletons_];
  ColSingletonPivotLIDs_ = new int[NumColSingletons_];
  ColSingletonPivots_ = new double[NumColSingletons_];
  
  // Register singleton columns (that were not already counted as singleton rows)
  // Check to see if any columns disappeared because all associated rows were eliminated
  int NumColSingletonstmp = 0;
  for (j=0; j<NumMyCols; j++) {
    int i = RowIDs[j];
    if ( ColProfiles[j]==1 && RowMapColors[i]!=1) {
      ColSingletonRowLIDs_[NumColSingletonstmp] = i;
      ColSingletonColLIDs_[NumColSingletonstmp] = j;
      NumColSingletonstmp++;
    }
    // Also check for columns that were eliminated implicitly by 
    // having all associated row eliminated
    else if (NewColProfiles[j]==0 && ColHasRowWithSingleton[j]!=1 && RowMapColors[i]!=1) {
	  ColMapColors[j] = 1;
    }
  }

  assert(NumColSingletonstmp==NumColSingletons_); //Sanity check
  return(0);
}
