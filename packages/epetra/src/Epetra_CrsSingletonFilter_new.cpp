
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
#include "Epetra_CrsSingletonFilter_new.h"

//==============================================================================
Epetra_CrsSingletonFilter::Epetra_CrsSingletonFilter(Epetra_LinearProblem * Problem) {
  InitializeDefaults();
  int ierr = Setup(Problem);
  if (ierr!=0) throw Problem->GetRHS()->ReportError("Error in CrsSingletonFilter Setup routine",ierr);
}
//==============================================================================
Epetra_CrsSingletonFilter::~Epetra_CrsSingletonFilter(){


  if (ReducedProblem_!=0) delete ReducedProblem_;
  if (ReducedMatrix_!=0) delete ReducedMatrix_;
  if (ReducedLHS_!=0) delete ReducedLHS_;
  if (ReducedRHS_!=0) delete ReducedRHS_;
  if (RowEliminateMap_!=0) delete RowEliminateMap_;
  if (ColEliminateMap_!=0) delete ColEliminateMap_;
  if (ReducedMatrixDomainMap_!=ReducedMatrixRowMap_) delete ReducedMatrixDomainMap_;
  if (ReducedMatrixRowMap_!=0) delete ReducedMatrixRowMap_;
  if (Full2ReducedRHSImporter_!=0) delete Full2ReducedRHSImporter_;
  if (Full2ReducedLHSImporter_!=0) delete Full2ReducedLHSImporter_;
  if (RedistributeDomainExporter_!=0) delete RedistributeDomainExporter_;

  if (ColSingletonRowLIDs_ != 0) delete ColSingletonRowLIDs_;
  if (ColSingletonColLIDs_ != 0) delete ColSingletonColLIDs_;
  if (ColSingletonPivots_ != 0) delete ColSingletonPivots_;
  if (tempExportX_ != 0) delete tempExportX_;
  if (tempX_ != 0) delete tempX_;
  if (tempB_ != 0) delete tempB_;
  if (ReducedIndices_ != 0) delete [] ReducedIndices_;
  if (ReducedValues_ != 0) delete [] ReducedValues_;

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
  RowEliminateMap_ = 0;
  ColEliminateMap_ = 0;
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
  ColSingletonPivots_ = 0;

  AbsoluteThreshold_ = 0;
  RelativeThreshold_ = 0;

  HaveReducedProblem_ = false;
  UserDefinedEliminateMaps_ = false;
  AnalysisDone_ = false;
  SymmetricElimination_ = true;

  tempExportX_ = 0;
  tempX_ = 0;
  tempB_ = 0;
  ReducedIndices_ = 0;
  ReducedValues_ = 0;

  Values_ = 0;
  Indices_ = 0;

  RowMapColors_ = 0;
  ColMapColors_ = 0;

  FullMatrixIsCrsMatrix_ = false;
  MaxNumMyEntries_ = 0;
  return;
}
//==============================================================================
int  Epetra_CrsSingletonFilter::Setup(Epetra_LinearProblem * Problem) {

  if (Problem==0) return(-1); // Null problem pointer

  FullProblem_ = Problem;
  FullMatrix_ = dynamic_cast<Epetra_RowMatrix *>(Problem->GetMatrix());
  if (FullMatrix_==0) EPETRA_CHK_ERR(-2); // Need a RowMatrix
  if (Problem->GetRHS()==0) EPETRA_CHK_ERR(-3); // Need a RHS
  if (Problem->GetLHS()==0) EPETRA_CHK_ERR(-4); // Need a LHS
  return(0);
}

//==============================================================================
int Epetra_CrsSingletonFilter::Analyze() {

  int i, j, jj;

  if (AnalysisDone_) EPETRA_CHK_ERR(-1); // Analysis already done once.  Cannot do it again

  // First check for columns with single entries and find columns with singleton rows
  Epetra_IntVector ColProfiles(FullMatrixColMap());
  Epetra_IntVector ColHasRowWithSingleton(FullMatrixColMap());

  // Define MapColoring objects
  RowMapColors_ = new Epetra_MapColoring(FullMatrixRowMap());  // Initial colors are all 0
  ColMapColors_ = new Epetra_MapColoring(FullMatrixColMap());
  Epetra_MapColoring & RowMapColors = *RowMapColors_;
  Epetra_MapColoring & ColMapColors = *ColMapColors_;


  int NumMyRows = FullMatrix()->NumMyRows();
  int NumMyCols = FullMatrix()->NumMyCols();

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

  Epetra_IntVector tmpVec(FullMatrixDomainMap()); // Use for gather/scatter of column vectors
  EPETRA_CHK_ERR(tmpVec.Export(ColProfiles, *FullMatrix()->RowMatrixImporter(), Add));
  EPETRA_CHK_ERR(ColProfiles.Import(tmpVec, *FullMatrix()->RowMatrixImporter(), Insert));
  
  EPETRA_CHK_ERR(tmpVec.PutValue(0));
  EPETRA_CHK_ERR(tmpVec.Export(ColHasRowWithSingleton, *FullMatrix()->RowMatrixImporter(), Add));
  EPETRA_CHK_ERR(ColHasRowWithSingleton.Import(tmpVec, *FullMatrix()->RowMatrixImporter(), Insert));

  // ColProfiles now contains the nonzero column entry count for all columns that have
  // an entry on this processor.
  // ColHasRowWithSingleton now contains a count of singleton rows associated with the corresponding
  // local column.  Next we check to make sure no column is associated with more than one singleton row.

  if (ColHasRowWithSingleton.MaxValue()>1) {
    EPETRA_CHK_ERR(-2); // At least one col is associated with two singleton rows, can't handle it.
  }


  Epetra_IntVector RowHasColWithSingleton(FullMatrix()->RowMatrixRowMap()); // Use to check for errors
 
  NumColSingletons_ = 0;
  // Count singleton columns (that were not already counted as singleton rows)
  for (j=0; j<NumMyCols; j++) {
    int i = RowIDs[j];
    // Check if column is a singleton
    if (ColProfiles[j]==1) {
      // Check to see if this column already eliminated by the row check above
      if (RowMapColors[i]!=1) {
	RowHasColWithSingleton[i]++; // Increment col singleton counter for ith row
	RowMapColors[i] = 1;
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

  RowEliminateMap_ = RowMapColors.GenerateMap(1); // Generate map associated with all eliminated elements.
  ColEliminateMap_ = ColMapColors.GenerateMap(1);

  if (RowIDs!=0) delete [] RowIDs;

  AnalysisDone_ = true;
  return(0);
}
//==============================================================================
int Epetra_CrsSingletonFilter::ConstructReducedProblem() {

  int i, j;
  // Generate reduced row and column maps

  Epetra_MapColoring & RowMapColors = *RowMapColors_;
  Epetra_MapColoring & ColMapColors = *ColMapColors_;

  ReducedMatrixRowMap_ = RowMapColors.GenerateMap(0);
  ReducedMatrixColMap_ = ColMapColors.GenerateMap(0);

  // Create domain and range map colorings by exporting map coloring of column and row maps

  Epetra_MapColoring DomainMapColors(FullMatrixDomainMap());
  EPETRA_CHK_ERR(DomainMapColors.Export(*ColMapColors_, *FullMatrix()->RowMatrixImporter(), AbsMax));
  OrigReducedMatrixDomainMap_ = DomainMapColors.GenerateMap(0);

  if (FullMatrixIsCrsMatrix_) {
    if (FullCrsMatrix()->Exporter()!=0) { // Non-trivial exporter
      Epetra_MapColoring RangeMapColors(FullMatrixRangeMap());
      EPETRA_CHK_ERR(RangeMapColors.Export(*RowMapColors_, *FullCrsMatrix()->Exporter(), 
					   AbsMax));
      ReducedMatrixRangeMap_ = RangeMapColors.GenerateMap(0);
    }
    else
      ReducedMatrixRangeMap_ =ReducedMatrixRowMap_;
  }
  else
    ReducedMatrixRangeMap_ =ReducedMatrixRowMap_;

  // Check to see if the reduced system domain and range maps are the same.
  // If not, we need to remap entries of the LHS multivector so that they are distributed
  // conformally with the rows of the reduced matrix and the RHS multivector
  SymmetricElimination_ = ReducedMatrixRangeMap_->SameAs(*OrigReducedMatrixDomainMap_);
  if (!SymmetricElimination_) 
    ConstructRedistributeExporter(ReducedMatrixRangeMap_, OrigReducedMatrixDomainMap_, 
				  RedistributeDomainExporter_, ReducedMatrixDomainMap_);
  else {
    ReducedMatrixDomainMap_ = OrigReducedMatrixDomainMap_;
    OrigReducedMatrixDomainMap_ = 0;
  }
  
  // Create pointer to Full RHS, LHS
  Epetra_MultiVector * FullRHS = FullProblem()->GetRHS();
  Epetra_MultiVector * FullLHS = FullProblem()->GetLHS();
  int NumVectors = FullLHS->NumVectors();

  // Create importers
  Full2ReducedLHSImporter_ = new Epetra_Import(*ReducedMatrixDomainMap(), FullMatrixDomainMap());
  Full2ReducedRHSImporter_ = new Epetra_Import(*ReducedMatrixRowMap(), FullRHS->Map());

  // Construct Reduced Matrix

  if (ReducedMatrix_!=0) delete ReducedMatrix_; // Get rid of any existing Reduced Matrix
  ReducedMatrix_ = new Epetra_CrsMatrix(Copy, *ReducedMatrixRowMap(), *ReducedMatrixColMap(), 0);

  // Create storage for temporary X values due to explicit elimination of rows
  tempExportX_ = new Epetra_MultiVector(FullMatrixColMap(), NumVectors);

  //cout << "FullLHS = " << endl << *FullLHS << endl;
  //cout << "FullRHS = " << endl << *FullRHS << endl;
  //cout << "tempExportX = " << endl << *tempExportX << endl;

  if (MaxNumMyEntries_>0) {
    ReducedIndices_ = new int[MaxNumMyEntries_];
    ReducedValues_ = new double[MaxNumMyEntries_];
  }
  int NumEntries;
  int * Indices;
  double * Values;
  int NumMyRows = FullMatrix()->NumMyRows();
  int ColSingletonCounter = 0;
  for (i=0; i<NumMyRows; i++) {
    int curGRID = FullMatrixRowMap().GID(i);
    EPETRA_CHK_ERR(GetRow(i, NumEntries, Values, Indices)); // Get current row
    if (ReducedMatrixRowMap()->MyGID(curGRID)) { // Check if this row should go into reduced matrix
      int ReducedNumEntries = 0;
      for (j=0; j<NumEntries; j++) {
	int ColumnIndex = FullMatrixColMap().GID(Indices[j]);
	if (ReducedMatrixColMap()->MyGID(ColumnIndex)) {
	  ReducedIndices_[ReducedNumEntries] = ColumnIndex;
	  ReducedValues_[ReducedNumEntries++] = Values[j];
	}
      }
      EPETRA_CHK_ERR(ReducedMatrix()->InsertGlobalValues(curGRID, ReducedNumEntries, 
							 ReducedValues_, ReducedIndices_));
    }
    // Otherwise if singleton row we explicitly eliminate this row and solve for corresponding X value
    else if (NumEntries==1) {
      double pivot = Values[0];
      if (pivot==0.0) EPETRA_CHK_ERR(-1); // Encountered zero row, unable to continue
      int indX = Indices[0];
      for (j=0; j<NumVectors; j++)
	(*tempExportX_)[j][indX] = (*FullRHS)[j][i]/pivot;
    }
    // Otherwise, this is a singleton column and we will scan for the pivot element needed 
    // for post-solve equations
    else {
      assert(i==ColSingletonRowLIDs_[ColSingletonCounter]);  // Sanity test
      int targetCol = ColSingletonColLIDs_[ColSingletonCounter];
      for (j=0; j<NumEntries; j++) {
	if (Indices[j]==targetCol) {
	  double pivot = Values[j];
	  if (pivot==0.0) EPETRA_CHK_ERR(-2); // Encountered zero column, unable to continue
	  ColSingletonPivots_[ColSingletonCounter] = pivot;
	  ColSingletonCounter++;
	  break;
	}
      }
    }
  }

  assert(ColSingletonCounter==NumColSingletons_); // Sanity test

  // Now convert to local indexing.
  EPETRA_CHK_ERR(ReducedMatrix()->TransformToLocal(ReducedMatrixDomainMap(), ReducedMatrixRangeMap()));

  // Construct Reduced LHS (Puts any initial guess values into reduced system)

  if (ReducedLHS_!=0) delete ReducedLHS_; // Get rid of any existing Reduced LHS
  ReducedLHS_ = new Epetra_MultiVector(*ReducedMatrixDomainMap(), NumVectors);
  EPETRA_CHK_ERR(ReducedLHS_->Import(*FullLHS, *Full2ReducedLHSImporter_, Insert));
  FullLHS->PutScalar(0.0); // zero out Full LHS since we will inject values as we get them

  // Construct Reduced RHS

  // First compute influence of already-known values of X on RHS
  tempX_ = new Epetra_MultiVector(FullMatrixDomainMap(), NumVectors);
  tempB_ = new Epetra_MultiVector(FullRHS->Map(), NumVectors);
  
  //Inject known X values into tempX for purpose of computing tempB = FullMatrix*tempX
  // Also inject into full X since we already know the solution

  EPETRA_CHK_ERR(tempX_->Export(*tempExportX_, *FullMatrix()->RowMatrixImporter(), Add));
  EPETRA_CHK_ERR(FullLHS->Export(*tempExportX_, *FullMatrix()->RowMatrixImporter(), Add));

  //cout << "tempExportX = " << endl << *tempExportX << endl;
  //cout << "tempX = " << endl << *tempX << endl;
  //cout << "tempB = " << endl << *tempB << endl;

  EPETRA_CHK_ERR(FullMatrix()->Multiply(false, *tempX_, *tempB_));

  EPETRA_CHK_ERR(tempB_->Update(1.0, *FullRHS, -1.0)); // tempB now has influence of already-known X values

  if (ReducedRHS_!=0) delete ReducedRHS_; // Get rid of any existing Reduced RHS
  ReducedRHS_ = new Epetra_MultiVector(*ReducedMatrixRowMap(), FullRHS->NumVectors());
  EPETRA_CHK_ERR(ReducedRHS_->Import(*tempB_, *Full2ReducedRHSImporter_, Insert));

  // Finally construct Reduced Linear Problem

  ReducedProblem_ = new Epetra_LinearProblem(ReducedMatrix_, ReducedLHS_, ReducedRHS_);

  HaveReducedProblem_ = true;

  
  return(0);
}
//==============================================================================
int Epetra_CrsSingletonFilter::UpdateReducedProblem() {

  if (!HaveReducedProblem_) EPETRA_CHK_ERR(-4); // Must have an existing reduced problem


    return(0);
}
//==============================================================================
int Epetra_CrsSingletonFilter::Statistics() const {

  if (!HaveReducedProblem_) {
    cout << "ConstructReducedProblem method must be called first." << endl;
    return(0);
  }

  double fn = FullMatrix()->NumGlobalRows();
  double fnnz = FullMatrix()->NumGlobalNonzeros();
  double rn = ReducedMatrix()->NumGlobalRows();
  double rnnz = ReducedMatrix()->NumGlobalNonzeros();
  if (fn==0.0 || fnnz==0.0) {
    cout << "Full problem has zero size." << endl;
    return(0);
  }

  cout << "Full System characteristics:" << endl << endl
       << "  Dimension                             = " << FullMatrix()->NumGlobalRows() << endl
       << "  Number of nonzeros                    = " << FullMatrix()->NumGlobalNonzeros() << endl
       << "  Maximum Number of Row Entries         = " << MaxNumMyEntries_ << endl << endl
       << "Reduced System characteristics:" << endl << endl
       << "  Dimension                             = " << ReducedMatrix()->NumGlobalRows() << endl
       << "  Number of nonzeros                    = " << ReducedMatrix()->NumGlobalNonzeros() << endl
       << "  Maximum Number of Row Entries         = " << ReducedMatrix()->GlobalMaxNumEntries() << endl << endl
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

  // Now create a vector that contains the global indices of the Source Epetra_MultiVector
  Epetra_IntVector SourceIndices(View, ContiguousSourceMap, SourceMap->MyGlobalElements());

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
int Epetra_CrsSingletonFilter::Analyze(int AbsoluteThreshold, double RelativeThreshold) {
  return(0);
}
//==============================================================================
int Epetra_CrsSingletonFilter::ComputeFullSolution() {

  int i, jj, k;

  Epetra_MultiVector * FullLHS = FullProblem()->GetLHS(); 
  Epetra_MultiVector * FullRHS = FullProblem()->GetRHS(); 
  Epetra_MultiVector tempX(FullLHS->Map(),FullLHS->NumVectors());

  // Inject values that the user computed for the reduced problem into the full solution vector
  EPETRA_CHK_ERR(tempX.Export(*ReducedLHS_, *Full2ReducedLHSImporter_, Add));
  FullLHS->Update(1.0, tempX, 1.0);

  // Next we will use our full solution vector which is populated with pre-filter solution
  // values and reduced system solution values to compute the sum of the row contributions
  // that must be subtracted to get the post-filter solution values

  Epetra_MultiVector tempAx(*FullRHS);
  EPETRA_CHK_ERR(FullMatrix()->Multiply(false, *FullLHS, tempAx));


  // Create storage for temporary X since we may be updating X values that belong to other processors
  Epetra_MultiVector tempImportX(FullMatrixColMap(), FullLHS->NumVectors());

  // Finally we loop through the local rows that were associated with column singletons and compute the
  // solution for these equations.

  for (k=0; k<NumColSingletons_; k++) {
    int i = ColSingletonRowLIDs_[k];
    int j = ColSingletonColLIDs_[k];
    double pivot = ColSingletonPivots_[k];
    int NumVectors = tempAx.NumVectors();
    for (jj=0; jj<NumVectors; jj++)
      tempImportX[jj][j]= ((*FullRHS)[jj][i] - tempAx[jj][i])/pivot;
  }

  // Finally, insert values from post-solve step and we are done!!!!

  
  Epetra_MultiVector tempX2(FullLHS->Map(),FullLHS->NumVectors());
  EPETRA_CHK_ERR(tempX.Export(tempImportX, *FullMatrix()->RowMatrixImporter(), Add));

  FullLHS->Update(1.0, tempX, 1.0);
    
  return(0);
}
//==============================================================================
int Epetra_CrsSingletonFilter::InitFullMatrixAccess() {

  int i;
  // Cast to CrsMatrix, if possible.  Can save some work.
  FullCrsMatrix_ = dynamic_cast<Epetra_CrsMatrix *>(FullMatrix());
  FullMatrixIsCrsMatrix_ = (FullCrsMatrix_!=0); // Pointer is non-zero if cast worked


  // Set up for scan of matrix graph, depends on whether matrix is CrsMatrix or general RowMatrix
  if (FullMatrixIsCrsMatrix_)
    MaxNumMyEntries_ = FullCrsMatrix_->MaxNumEntries(); 
  else {
    int NumMyRows = FullMatrix()->NumMyRows();
    int NumEntries;
    MaxNumMyEntries_ = 0;
    for (i=0; i<NumMyRows; i++) { 
      FullMatrix()->NumMyRowEntries(i, NumEntries);
      MaxNumMyEntries_ = EPETRA_MAX(MaxNumMyEntries_, NumEntries);
    }
    Indices_ = new int[MaxNumMyEntries_];
    Values_ = new double[MaxNumMyEntries_];
  }
  
  return(0);
}
//==============================================================================
int Epetra_CrsSingletonFilter::GetRow(int Row, int & NumIndices, int * & Indices) {
  int i;
  if (FullMatrixIsCrsMatrix_) { // View of current row
    EPETRA_CHK_ERR(FullCrsMatrix()->Graph().ExtractMyRowView(i, NumIndices, Indices)); 
  }
  else { // Copy of current row (we must get the values, but we ignore them)
    EPETRA_CHK_ERR(FullMatrix()->ExtractMyRowCopy(i, MaxNumMyEntries_, NumIndices, 
						  Values_, Indices_));
    Indices = Indices_;
  } 
  return(0);
}
//==============================================================================
int Epetra_CrsSingletonFilter::GetRow(int Row, int & NumIndices, 
				      double * & Values, int * & Indices) {

  int i;
  if (FullMatrixIsCrsMatrix_) { // View of current row
    EPETRA_CHK_ERR(FullCrsMatrix_->ExtractMyRowView(i, NumIndices, Values, Indices)); 
  }
  else { // Copy of current row (we must get the values, but we ignore them)
    EPETRA_CHK_ERR(FullMatrix()->ExtractMyRowCopy(i, MaxNumMyEntries_, NumIndices, 
						  Values_, Indices_));
    Values = Values_;
    Indices = Indices_;
  } 
  return(0);
}
//==============================================================================
int Epetra_CrsSingletonFilter::CreatePostSolveArrays(int * RowIDs,
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
