
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
#include "Epetra_Import.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_CrsSingletonFilter.h"

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
  if (ReducedMatrixRowMap_!=0) delete ReducedMatrixRowMap_;
  if (ReducedMatrixDomainMap_!=0) delete ReducedMatrixDomainMap_;
  if (Full2ReducedRHSImporter_!=0) delete Full2ReducedRHSImporter_;
  if (Full2ReducedLHSImporter_!=0) delete Full2ReducedLHSImporter_;

  if (ColSingletonRowLIDs_ != 0) delete ColSingletonRowLIDs_;
  if (ColSingletonColLIDs_ != 0) delete ColSingletonColLIDs_;
  if (ColSingletonPivots_ != 0) delete ColSingletonPivots_;

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
  ReducedMatrixDomainMap_ = 0;
  Full2ReducedRHSImporter_ = 0;
  Full2ReducedLHSImporter_ = 0;
    
  ColSingletonRowLIDs_ = 0;
  ColSingletonColLIDs_ = 0;
  ColSingletonPivots_ = 0;

  AbsoluteThreshold_ = 0;
  RelativeThreshold_ = 0;

  HaveReducedProblem_ = false;
  UserDefinedEliminateMaps_ = false;
  AnalysisDone_ = false;
  return;
}
//==============================================================================
int  Epetra_CrsSingletonFilter::Setup(Epetra_LinearProblem * Problem) {

  if (Problem==0) return(-1); // Null problem pointer

  FullProblem_ = Problem;
  FullMatrix_ = dynamic_cast<Epetra_CrsMatrix *>(Problem->GetMatrix());
  if (FullMatrix_==0) EPETRA_CHK_ERR(-2); // Matrix is not an Epetra_CrsMatrix (required for this conversion).
  if (Problem->GetRHS()==0) EPETRA_CHK_ERR(-3); // Need a RHS
  if (Problem->GetLHS()==0) EPETRA_CHK_ERR(-4); // Need a LHS
  return(0);
}

//==============================================================================
int Epetra_CrsSingletonFilter::Analyze() {

  if (AnalysisDone_) EPETRA_CHK_ERR(-1); // Analysis already done once.  Cannot do it again

  // First check for columns with single entries and find columns with singleton rows
  Epetra_Vector LocalColProfiles(FullMatrixImportMap());
  Epetra_Vector ColProfiles(FullMatrixDomainMap());

  Epetra_Vector LocalColWithSingletonRows(FullMatrixImportMap());
  Epetra_Vector ColWithSingletonRows(FullMatrixDomainMap());

  int NumMyRows = FullMatrix()->NumMyRows();
  int NumMyCols = FullMatrix()->NumMyCols();
  int IndexBase = FullMatrixRowMap().IndexBase();

  int * LocalRowIDs = 0;
  // LocalRowIDs[j] will contain the local row ID associated with the jth column, if the jth col has a single entry
  if (NumMyCols>0) LocalRowIDs = new int[NumMyCols];
  for (int i=0; i<NumMyCols; i++) LocalRowIDs[i] = -1;

  // The follow two arrays will keep track of any rows and columns that should be eliminated.
  int * EliminatedRowIDs = 0;
  int * EliminatedColIDs = 0;
  if (NumMyRows>0) EliminatedRowIDs = new int[NumMyRows];
  if (NumMyCols>0) EliminatedColIDs = new int[NumMyCols];

  int NumEliminatedRows = 0;
  int NumEliminatedCols = 0;

  int NumIndices;
  int * Indices;
  {for (int i=0; i<NumMyRows; i++) {
    EPETRA_CHK_ERR(FullMatrix()->Graph().ExtractMyRowView(i, NumIndices, Indices)); // View of current row
    for (int j=0; j<NumIndices; j++) {
      int ColumnIndex = Indices[j];
      LocalColProfiles[ColumnIndex] += 1.0; // Increment column count
      // Record local row ID for current column
      // will use to identify row to eliminate if column is a singleton
      LocalRowIDs[ColumnIndex] = i;
    }
    // Now check for rows with single entries and
    // Add a 1 in columns eliminated by singleton row
    if (NumIndices==1) {
      int j = Indices[0];
      //cout << "i, j = " << i << "  " << j << endl;
      LocalColWithSingletonRows[j] += 1.0;
      EliminatedRowIDs[NumEliminatedRows++] = FullMatrixRowMap().GID(i);
      EliminatedColIDs[NumEliminatedCols++] = FullMatrixImportMap().GID(j);
      NumRowSingletons_++;
    }
  }}

  // 1) The vector LocalColProfiles has column nonzero counts for each processor's contribution
  // Combine these to get total column profile information and then redistribute to processors
  // so each can determine if it is the owner of the row associated with the singleton column
  // 2) The vector LocalColWithSingletonRows[i] contain count of singleton rows  that are associated with 
  // the ith column on this processor.  Must tell other processors that they should also eliminate these columns

  // Make a copy of LocalColProfiles for later use when detecting columns that disappear locally
  Epetra_Vector LocalColProfilesCopy(LocalColProfiles);

  if (LocalColProfiles.Map().DistributedGlobal()) {
    EPETRA_CHK_ERR(ColProfiles.Export(LocalColProfiles, *FullMatrix()->Importer(), Add));
    EPETRA_CHK_ERR(LocalColProfiles.Import(ColProfiles, *FullMatrix()->Importer(), Insert));
    EPETRA_CHK_ERR(ColWithSingletonRows.Export(LocalColWithSingletonRows, *FullMatrix()->Importer(), Add));
    EPETRA_CHK_ERR(LocalColWithSingletonRows.Import(ColWithSingletonRows, *FullMatrix()->Importer(), Insert));
  }

  // LocalColProfiles now contains the nonzero column entry count for all columns that have
  // an entry on this processor.
  // LocalColWithSingletonRows now contains a count of singleton rows associated with the corresponding
  // local column.


  // If the ith element of this array is true, then ith row is associated with a column with a singleton
  bool * RowsWithColSingletons = new bool[NumMyRows];
  {for (int i=0; i<NumMyRows; i++) RowsWithColSingletons[i] = false;} // Init to false

  NumColSingletons_ = 0;

  // Count singleton columns (that were not already counted as singleton rows)
  for (int j=0; j<NumMyCols; j++) {
    if (LocalColProfiles[j]==1.0) {
      int i = LocalRowIDs[j];
      if (i==-1) { // If column has no row entry on this processor, just eliminate the column
	EliminatedColIDs[NumEliminatedCols++] = FullMatrixImportMap().GID(j);
      }
      else { // Check to see if this column already eliminated by the row check above
	int k = FullMatrix()->NumMyEntries(i);
	if (k!=1) {
	  if (RowsWithColSingletons[i]) EPETRA_CHK_ERR(-2); // This row has two column singletons, can't handle it.
	  RowsWithColSingletons[i] = true;
	  EliminatedRowIDs[NumEliminatedRows++] = FullMatrixRowMap().GID(i);
	  EliminatedColIDs[NumEliminatedCols++] = FullMatrixImportMap().GID(j);
	  NumColSingletons_++;
	  // If we delete a row, we need to keep track of associated column entries that were also deleted 
	  // in case all entries in a column are eventually deleted, in which case the column should
	  // also be deleted.
	  EPETRA_CHK_ERR(FullMatrix()->Graph().ExtractMyRowView(i, NumIndices, Indices)); // View of current row
	  for (int jj=0; jj<NumIndices; jj++) LocalColProfilesCopy[Indices[jj]] -= 1.0;

	}
      }
    }
    else if (LocalColWithSingletonRows[j]==1.0) { // Check if some other processor eliminated this column
      if (LocalRowIDs[j]!=-1) { // If column has no row entry on this processor, just eliminate the column
	int k = FullMatrix()->NumMyEntries(LocalRowIDs[j]);
	if (k!=1 && !FullMatrixRowMap().MyLID(j))
	  EliminatedColIDs[NumEliminatedCols++] = FullMatrixImportMap().GID(j);
      }
    }
    else if (LocalColWithSingletonRows[j]>1.0) {
      EPETRA_CHK_ERR(-3); // This column has two row singletons, can't handle it.
    }
  }

  delete [] RowsWithColSingletons;


  // We will need these arrays for the post-solve phase
  if (NumColSingletons_>0) {
    ColSingletonRowLIDs_ = new int[NumColSingletons_];
    ColSingletonColLIDs_ = new int[NumColSingletons_];
    ColSingletonPivots_ = new double[NumColSingletons_];
  }

  // Register singleton columns (that were not already counted as singleton rows)
  // Check to see if any columns disappeared because all associated rows were eliminated
  int NumColSingletonstmp = 0;
  {for (int j=0; j<NumMyCols; j++) {
    if ( LocalColProfiles[j]==1.0) {
      int i = LocalRowIDs[j];
      if (i!=-1) { // Make sure column has an entry on this PE
        int k = FullMatrix()->NumMyEntries(i);
        if (k!=1) {
          ColSingletonRowLIDs_[NumColSingletonstmp] = i;
          ColSingletonColLIDs_[NumColSingletonstmp] = j;
	        NumColSingletonstmp++;
        }
      }
    }
    else if (LocalColProfilesCopy[j]==0.0 && LocalColWithSingletonRows[j]!=1.0) {
      if (LocalRowIDs[j]!=-1) { // If column has no row entry on this processor, just eliminate the column
	int k = FullMatrix()->NumMyEntries(LocalRowIDs[j]);
	if (k!=1)      
	  EliminatedColIDs[NumEliminatedCols++] = FullMatrixImportMap().GID(j);
      }
    }
  }}

  assert(NumColSingletonstmp==NumColSingletons_); //Sanity check

  // Sort ColSingleton LID arrays so that Row LIDs are in ascending order
  Epetra_Util util;
  util.Sort(true, NumColSingletons_, ColSingletonRowLIDs_, 0, 0, 1, &ColSingletonColLIDs_);

  RowEliminateMap_ = new Epetra_Map(-1, NumEliminatedRows, EliminatedRowIDs, IndexBase, FullMatrixImportMap().Comm());
  ColEliminateMap_ = new Epetra_Map(-1, NumEliminatedCols, EliminatedColIDs, IndexBase, FullMatrixImportMap().Comm());

  if (LocalRowIDs!=0) delete [] LocalRowIDs;
  if (EliminatedRowIDs!=0) delete [] EliminatedRowIDs;
  if (EliminatedColIDs!=0) delete [] EliminatedColIDs;

  int * RowReducedList = 0;
  int * ColReducedList = 0;

  int NumReducedRows = NumMyRows - NumEliminatedRows;
  int NumReducedCols = NumMyCols - NumEliminatedCols;


  if (NumReducedRows>0) RowReducedList = new int[NumReducedRows];
  if (NumReducedCols>0) ColReducedList = new int[NumReducedCols];

  int NumReducedRowstmp = 0;
  int NumReducedColstmp = 0;

  // Register reduced rows by determining if each row of full matrix is part of RowEliminateMap
  {for (int i=0; i<NumMyRows; i++) {
    int curGID = FullMatrixRowMap().GID(i); 
    if (!RowEliminateMap()->MyGID(curGID)) RowReducedList[NumReducedRowstmp++] = curGID;
  }}

  // Register reduced domain by determining if each column of full matrix is part of ColEliminateMap
  {for (int i=0; i<NumMyCols; i++) {
    int curGID = FullMatrixImportMap().GID(i); 
    if (!ColEliminateMap()->MyGID(curGID)) ColReducedList[NumReducedColstmp++] = curGID;
  }}
  //cout << "ColumnEliminateMap = "<< *ColEliminateMap_ << endl;
  //cout << "Row   EliminateMap = "<< *RowEliminateMap_ << endl;
  assert(NumReducedRowstmp==NumReducedRows); //Sanity check
  assert(NumReducedColstmp==NumReducedCols); //Sanity check

  {for (int i=0; i<NumEliminatedRows; i++) {
    int curGID = RowEliminateMap()->GID(i); 
    if (!ColEliminateMap()->MyGID(curGID)) EPETRA_CHK_ERR(-4); // Cannot handle case of nonsymmetric elimination (yet)
  }}
  

  // Construct Reduced matrix maps
  ReducedMatrixRowMap_ = new Epetra_Map(-1, NumReducedRows, RowReducedList, IndexBase,
					FullMatrixImportMap().Comm());
  ReducedMatrixDomainMap_ = new Epetra_Map(-1, NumReducedRows, RowReducedList, IndexBase,
					   FullMatrixImportMap().Comm());
  // ReducedMatrixDomainMap_ = new Epetra_Map(-1, NumReducedCols, ColReducedList, IndexBase,
  //					   FullMatrixImportMap().Comm());


  // Delete local arrays
  if (RowReducedList!=0) delete [] RowReducedList;
  if (ColReducedList!=0) delete [] ColReducedList;

  // Create importers
  Full2ReducedLHSImporter_ = new Epetra_Import(*ReducedMatrixDomainMap(), FullMatrixDomainMap());
  Full2ReducedRHSImporter_ = new Epetra_Import(*ReducedMatrixRowMap(), FullProblem()->GetRHS()->Map());
  
  AnalysisDone_ = true;
  return(0);
}
//==============================================================================
int Epetra_CrsSingletonFilter::Statistics() {

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
       << "  Maximum Number of Row Entries         = " << FullMatrix()->GlobalMaxNumEntries() << endl << endl
       << "Reduced System characteristics:" << endl << endl
       << "  Dimension                             = " << ReducedMatrix()->NumGlobalRows() << endl
       << "  Number of nonzeros                    = " << ReducedMatrix()->NumGlobalNonzeros() << endl
       << "  Maximum Number of Row Entries         = " << ReducedMatrix()->GlobalMaxNumEntries() << endl << endl
       << "Singleton information: " << endl
       << "  Number of rows with single entries    = " << NumRowSingletons_ << endl
       << "  Number of columns with single entries = " << NumColSingletons_ << endl << endl
       << "Ratios: " << endl
       << "  Percent reduction in dimension        = " << (fn-rn)/fn*100.0 << endl
       << "  Percent reduction in nonzero count    = " << (fnnz-rnnz)/fnnz*100.0 << endl << endl;

  return(0);

}
//==============================================================================
int Epetra_CrsSingletonFilter::ConstructReducedProblem() {

  // Create pointer to Full RHS, LHS
  Epetra_MultiVector * FullRHS = FullProblem()->GetRHS();
  Epetra_MultiVector * FullLHS = FullProblem()->GetLHS();

  // Construct Reduced Matrix

  if (ReducedMatrix_!=0) delete ReducedMatrix_; // Get rid of any existing Reduced Matrix
  ReducedMatrix_ = new Epetra_CrsMatrix(Copy, *ReducedMatrixRowMap(), 0);

  // Create storage for temporary X values due to explicit elimination of rows
  Epetra_MultiVector tempExportX(FullMatrixImportMap(), FullLHS->NumVectors());

  //cout << "FullLHS = " << endl << *FullLHS << endl;
  //cout << "FullRHS = " << endl << *FullRHS << endl;
  //cout << "tempExportX = " << endl << tempExportX << endl;

  int NumMyRows = FullMatrix()->NumMyRows();

  int * ReducedIndices;
  double * ReducedValues;
  int MaxNumEntries = FullMatrix()->MaxNumEntries();
  if (MaxNumEntries>0) {
    ReducedIndices = new int[MaxNumEntries];
    ReducedValues = new double[MaxNumEntries];
  }
  int NumEntries;
  int * Indices;
  double * Values;
  int ColSingletonCounter = 0;
  for (int i=0; i<NumMyRows; i++) {
    int curGRID = FullMatrixRowMap().GID(i);
    EPETRA_CHK_ERR(FullMatrix()->ExtractMyRowView(i, NumEntries, Values, Indices)); // View of current row
    if (ReducedMatrixRowMap()->MyGID(curGRID)) { // Check if this row should go into reduced matrix
      int ReducedNumEntries = 0;
      for (int j=0; j<NumEntries; j++) {
	int ColumnIndex = FullMatrixImportMap().GID(Indices[j]);
	if (!ColEliminateMap()->MyGID(ColumnIndex)) {
	  ReducedIndices[ReducedNumEntries] = ColumnIndex;
	  ReducedValues[ReducedNumEntries++] = Values[j];
	}
      }
      EPETRA_CHK_ERR(ReducedMatrix()->InsertGlobalValues(curGRID, 
							 ReducedNumEntries, ReducedValues, ReducedIndices));
    }
    // Otherwise if singleton row we explicitly eliminate this row and solve for corresponding X value
    else if (NumEntries==1) {
      double pivot = Values[0];
      if (pivot==0.0) EPETRA_CHK_ERR(-1); // Encountered zero row, unable to continue
      int indX = Indices[0];
      int NumVectors = tempExportX.NumVectors();
      for (int j=0; j<NumVectors; j++)
	tempExportX[j][indX] = (*FullRHS)[j][i]/pivot;
    }
    // Otherwise, this is a singleton column and we will scan for the pivot element needed for post-solve equations
    else {
      assert(i==ColSingletonRowLIDs_[ColSingletonCounter]);  // Sanity test
      int targetCol = ColSingletonColLIDs_[ColSingletonCounter];
      for (int j=0; j<NumEntries; j++) {
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

  EPETRA_CHK_ERR(ReducedMatrix()->TransformToLocal(ReducedMatrixDomainMap(), ReducedMatrixRowMap()));

  if (MaxNumEntries>0) {
    delete [] ReducedIndices;
    delete [] ReducedValues;
  }

  // Construct Reduced LHS (Puts any initial guess values into reduced system)

  if (ReducedLHS_!=0) delete ReducedLHS_; // Get rid of any existing Reduced LHS
  ReducedLHS_ = new Epetra_MultiVector(*ReducedMatrixDomainMap(), FullLHS->NumVectors());
  EPETRA_CHK_ERR(ReducedLHS_->Import(*FullLHS, *Full2ReducedLHSImporter_, Insert));
  FullLHS->PutScalar(0.0); // zero out Full LHS since we will inject values as we get them

  // Construct Reduced RHS

  // First compute influence of already-known values of X on RHS
  Epetra_MultiVector tempX(FullMatrixDomainMap(), FullLHS->NumVectors());
  Epetra_MultiVector tempB(FullRHS->Map(), FullRHS->NumVectors());
  
  //Inject known X values into tempX for purpose of computing tempB = FullMatrix*tempX
  // Also inject into full X since we already know the solution

  EPETRA_CHK_ERR(tempX.Export(tempExportX, *FullMatrix()->Importer(), Add));
  EPETRA_CHK_ERR(FullLHS->Export(tempExportX, *FullMatrix()->Importer(), Add));

  //cout << "tempExportX = " << endl << tempExportX << endl;
  //cout << "tempX = " << endl << tempX << endl;
  //cout << "tempB = " << endl << tempB << endl;

  EPETRA_CHK_ERR(FullMatrix()->Multiply(false, tempX, tempB));

  EPETRA_CHK_ERR(tempB.Update(1.0, *FullRHS, -1.0)); // tempB now has influence of already-known X values

  if (ReducedRHS_!=0) delete ReducedRHS_; // Get rid of any existing Reduced RHS
  ReducedRHS_ = new Epetra_MultiVector(*ReducedMatrixRowMap(), FullRHS->NumVectors());
  EPETRA_CHK_ERR(ReducedRHS_->Import(tempB, *Full2ReducedRHSImporter_, Insert));

  // Finally construct Reduced Linear Problem

  ReducedProblem_ = new Epetra_LinearProblem(ReducedMatrix_, ReducedLHS_, ReducedRHS_);

  HaveReducedProblem_ = true;


    return(0);
}
//==============================================================================
int Epetra_CrsSingletonFilter::Analyze(int AbsoluteThreshold, double RelativeThreshold) {
  return(0);
}
//==============================================================================
int Epetra_CrsSingletonFilter::ComputeFullSolution() {

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
  Epetra_MultiVector tempImportX(FullMatrixImportMap(), FullLHS->NumVectors());

  // Finally we loop through the local rows that were associated with column singletons and compute the
  // solution for these equations.

  for (int k=0; k<NumColSingletons_; k++) {
    int i = ColSingletonRowLIDs_[k];
    int j = ColSingletonColLIDs_[k];
    double pivot = ColSingletonPivots_[k];
    int NumVectors = tempAx.NumVectors();
    for (int jj=0; jj<NumVectors; jj++)
      tempImportX[jj][j]= ((*FullRHS)[jj][i] - tempAx[jj][i])/pivot;
  }

  // Finally, insert values from post-solve step and we are done!!!!

  
  Epetra_MultiVector tempX2(FullLHS->Map(),FullLHS->NumVectors());
  EPETRA_CHK_ERR(tempX.Export(tempImportX, *FullMatrix()->Importer(), Add));

  FullLHS->Update(1.0, tempX, 1.0);
    
  return(0);
}


