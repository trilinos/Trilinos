
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
#include "Epetra_Import.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_ReducedLinearProblem.h"

//==============================================================================
Epetra_ReducedLinearProblem::Epetra_ReducedLinearProblem(Epetra_LinearProblem * Problem) {
  InitializeDefaults();
  int ierr = Setup(Problem);
  if (ierr!=0) throw Problem->GetRHS()->ReportError("Error in ReducedLinearProblem Setup routine",ierr);
}
//==============================================================================
Epetra_ReducedLinearProblem::~Epetra_ReducedLinearProblem(){

  if (!FullMatrixRowMap().SameAs(FullProblem()->GetRHS()->Map()))
    if (Full2ReducedRHSImporter_!=0) delete Full2ReducedRHSImporter_;

  if (ReducedMatrix_!=0) delete ReducedMatrix_;
  if (ReducedLHS_!=0) delete ReducedLHS_;
  if (ReducedRHS_!=0) delete ReducedRHS_;
  if (ReducedProblem_!=0) delete ReducedProblem_;
  if (RowEliminateMap_!=0) delete RowEliminateMap_;
  if (ColEliminateMap_!=0) delete ColEliminateMap_;
  if (ReducedMatrixRowMap_!=0) delete ReducedMatrixRowMap_;
  if (ReducedMatrixDomainMap_!=0) delete ReducedMatrixDomainMap_;
  if (Full2ReducedMatrixImporter_!=0) delete Full2ReducedMatrixImporter_;
  if (Full2ReducedRHSImporter_!=0) delete Full2ReducedRHSImporter_;
  if (Full2ReducedLHSImporter_!=0) delete Full2ReducedLHSImporter_;

}
//==============================================================================
void Epetra_ReducedLinearProblem::InitializeDefaults() { 

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
  Full2ReducedMatrixImporter_ = 0;
  Full2ReducedRHSImporter_ = 0;
  Full2ReducedLHSImporter_ = 0;
    
  AbsoluteThreshold_ = 0;
  RelativeThreshold_ = 0;

  HaveReducedProblem_ = false;
  UserDefinedEliminateMaps_ = false;
  AnalysisDone_ = false;
  return;
}
//==============================================================================
int  Epetra_ReducedLinearProblem::Setup(Epetra_LinearProblem * Problem) {

  if (Problem==0) return(-1); // Null problem pointer

  FullMatrix_ = dynamic_cast<Epetra_CrsMatrix *>(Problem->GetOperator());
  if (FullMatrix_==0) EPETRA_CHK_ERR(-2); // Matrix is not an Epetra_CrsMatrix (required for this conversion).
  if (Problem->GetRHS()==0) EPETRA_CHK_ERR(-3); // Need a RHS
  if (Problem->GetLHS()==0) EPETRA_CHK_ERR(-4); // Need a LHS
  return(0);
}

//==============================================================================
int Epetra_ReducedLinearProblem::Analyze() {

  if (AnalysisDone_) EPETRA_CHK_ERR(-1); // Analysis already done once.  Cannot do it again

  // First check for columns with single entries
  Epetra_Vector LocalColProfiles(FullMatrixImportMap());
  Epetra_Vector ColProfiles(FullMatrixDomainMap());

  int NumMyRows = FullMatrix().NumMyRows();
  int NumMyCols = FullMatrix().NumMyCols();
  int * LocalRowIDs = 0;
  if (NumMyCols>0) LocalRowIDs = new int[NumMyCols];
  int NumIndices;
  int * Indices;
  for (int i=0; i<NumMyRows; i++) {
    FullMatrix().Graph().ExtractMyRowView(i, NumIndices, Indices); // View of current row
    for (int j=0; j<NumIndices; j++) {
      int ColumnIndex = Indices[i];
      LocalColProfiles[ColumnIndex] += 1.0; // Increment column count
      // Record local row ID for current column
      // will use to identify row to eliminate if column is a singleton
      LocalRowIDs[ColumnIndex] = i;
    }
  }

  // The vector LocalColProfiles has column nonzero counts for each processor's contribution
  // Combine these to get total column profile information and then redistribute to processors
  // so each can determine if it is the owner of the row associated with the singleton column
  ColProfiles.Export(LocalColProfiles, *FullMatrix().Importer(), Add);
  LocalColProfiles.Import(ColProfiles, *FullMatrix().Importer(), Insert);

  // LocalColProfiles now contains the nonzero column entry count for all columns that have
  // an entry on this processor.

  
  // Now check for rows with single entries
  int NumSingletons = 0;
  // Count singleton rows
  for (int i=0; i<NumMyRows; i++)
    if (FullMatrix().NumMyEntries(i)==1) NumSingletons++;
  // Count singleton columns (that were not already counted as singleton rows)
  for (int i=0; i<NumMyCols; i++)
    if (ColProfiles[i]==1.0 && FullMatrix().NumMyEntries(LocalRowIDs[i])!=1) NumSingletons++;

  int * RowSingletonList = 0;
  int * ColSingletonList = 0;

  if (NumSingletons>0) {
    RowSingletonList = new int[NumSingletons];
    ColSingletonList = new int[NumSingletons];
  }
  NumRowSingletons_ = 0;
  NumColSingletons_ = 0;
  int ColumnIndex;
  int * ColumnIndexP = &ColumnIndex;
  // Register singleton rows
  for (int i=0; i<NumMyRows; i++)
    if (FullMatrix().NumMyEntries(i)==1) {
      FullMatrix().Graph().ExtractMyRowView(i, NumIndices, ColumnIndexP); // View of current row
      assert(NumIndices==1); // Sanity test
      RowSingletonList[NumRowSingletons_] = FullMatrixRowMap().GID(i);
      ColSingletonList[NumRowSingletons_] = FullMatrixImportMap().GID(ColumnIndex);
      NumRowSingletons_++;
    }


  // Register singleton columns (that were not already counted as singleton rows)
  for (int i=0; i<NumMyCols; i++)
    if (ColProfiles[i]==1.0 && FullMatrix().NumMyEntries(LocalRowIDs[i])!=1) {
      assert(NumIndices==1); // Sanity test
      RowSingletonList[NumRowSingletons_+NumColSingletons_] = FullMatrixRowMap().GID(LocalRowIDs[i]);
      ColSingletonList[NumRowSingletons_+NumColSingletons_] = FullMatrixImportMap().GID(i);
      NumColSingletons_++;
    }

  RowEliminateMap_ = new Epetra_Map(-1, NumSingletons, RowSingletonList, 0, FullMatrixImportMap().Comm());
  ColEliminateMap_ = new Epetra_Map(-1, NumSingletons, ColSingletonList, 0, FullMatrixImportMap().Comm());
  if (LocalRowIDs!=0) delete [] LocalRowIDs;
  if (RowSingletonList!=0) delete [] RowSingletonList;
  if (ColSingletonList!=0) delete [] ColSingletonList;

  int * RowReducedList = 0;
  int * ColReducedList = 0;

  if (NumMyRows - NumSingletons>0) RowReducedList = new int[NumMyRows - NumSingletons];
  if (NumMyCols - NumSingletons>0) ColReducedList = new int[NumMyCols - NumSingletons];

  int NumReducedRows = 0;
  int NumReducedCols = 0;

  // Register reduced rows by determining if each row of full matrix is part of RowEliminateMap
  for (int i=0; i<NumMyRows; i++) {
    int curGID = FullMatrixRowMap().GID(i); 
    if (!RowEliminateMap()->MyGID(curGID)) RowReducedList[NumReducedRows++] = curGID;
  }

  // Register reduced domain by determining if each column of full matrix is part of ColEliminateMap
  for (int i=0; i<NumMyCols; i++) {
    int curGID = FullMatrixDomainMap().GID(i); 
    if (!ColEliminateMap()->MyGID(curGID)) ColReducedList[NumReducedCols++] = curGID;
  }

  // Sanity tests
  if (NumReducedRows!=NumMyRows - NumSingletons) EPETRA_CHK_ERR(-1);
  if (NumReducedCols!=NumMyCols - NumSingletons) EPETRA_CHK_ERR(-2);

  // Construct Reduced matrix maps
  ReducedMatrixRowMap_ = new Epetra_Map(-1, NumReducedRows, RowReducedList, 0, FullMatrixImportMap().Comm());
  ReducedMatrixDomainMap_ = new Epetra_Map(-1, NumReducedCols, ColReducedList, 0, FullMatrixImportMap().Comm());

  // Delete local arrays
  if (RowReducedList!=0) delete [] RowReducedList;
  if (ColReducedList!=0) delete [] ColReducedList;

  // Create importers
  Full2ReducedMatrixImporter_ = new Epetra_Import(*ReducedMatrixRowMap(), FullMatrixRowMap());
  Full2ReducedLHSImporter_ = new Epetra_Import(*ReducedMatrixDomainMap(), FullMatrixDomainMap());

  // If full matrix and full RHS have same layout, don't need separate RHS Importer
  if (FullMatrixRowMap().SameAs(FullProblem()->GetRHS()->Map()))
    Full2ReducedRHSImporter_ = Full2ReducedMatrixImporter_;
  else
    Full2ReducedRHSImporter_ = new Epetra_Import(*ReducedMatrixDomainMap(), FullProblem()->GetRHS()->Map());
  
  AnalysisDone_ = true;
  return(0);
}
//==============================================================================
int Epetra_ReducedLinearProblem::ConstructReducedProblem() {

  // Construct Reduced Matrix

  if (ReducedMatrix_!=0) delete ReducedMatrix_; // Get rid of any existing Reduced Matrix
  ReducedMatrix_ = new Epetra_CrsMatrix(Copy, *ReducedMatrixRowMap(), 0);
  ReducedMatrix()->Import(FullMatrix(), *Full2ReducedMatrixImporter_, Insert);
  ReducedMatrix()->TransformToLocal(ReducedMatrixRowMap(), ReducedMatrixDomainMap());

  // Construct Reduced RHS

  if (ReducedRHS_!=0) delete ReducedRHS_; // Get rid of any existing Reduced RHS
  ReducedRHS_ = new Epetra_MultiVector(*ReducedMatrixRowMap(), FullProblem()->GetRHS()->NumVectors());
  ReducedRHS_->Import(*FullProblem()->GetRHS(), *Full2ReducedRHSImporter_, Insert);
  
  // Construct Reduced LHS

  if (ReducedLHS_!=0) delete ReducedLHS_; // Get rid of any existing Reduced LHS
  ReducedLHS_ = new Epetra_MultiVector(*ReducedMatrixDomainMap(), FullProblem()->GetLHS()->NumVectors());
  ReducedLHS_->Import(*FullProblem()->GetLHS(), *Full2ReducedLHSImporter_, Insert);

  HaveReducedProblem_ = true;

    return(0);
}
//==============================================================================
int Epetra_ReducedLinearProblem::Analyze(int AbsoluteThreshold, double RelativeThreshold) 
{
  return(0);
}

