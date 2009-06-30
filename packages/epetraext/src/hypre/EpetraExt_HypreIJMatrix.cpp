//@HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

#include "EpetraExt_ConfigDefs.h"
#ifdef HAVE_HYPRE

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "EpetraExt_HypreIJMatrix.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Array.hpp"
#include <set>

//=======================================================
EpetraExt_HypreIJMatrix::EpetraExt_HypreIJMatrix(HYPRE_IJMatrix matrix)
#ifdef HAVE_MPI
  : Epetra_BasicRowMatrix(Epetra_MpiComm(hypre_IJMatrixComm(matrix))),
#else
  : Epetra_BasicRowMatrix(Epetra_SerialComm()),
#endif
    Matrix_(matrix),
    ParMatrix_(0),
    NumMyRows_(-1),
    NumGlobalRows_(-1),
    NumGlobalCols_(-1),
    MyRowStart_(-1),
    MyRowEnd_(-1),
    MatType_(-1), 
    SolverCreated_(false),
    PrecondCreated_(false),
    TransposeSolve_(false),
    SolveOrPrec_(Solver)
{
  // Initialize default values for global variables
  int ierr = 0;
  ierr += InitializeDefaults();
  TEST_FOR_EXCEPTION(ierr != 0, std::logic_error, "Couldn't initialize default values.");
  
  // Create array of global row ids
  Teuchos::Array<int> GlobalRowIDs;  GlobalRowIDs.resize(NumMyRows_);
  
  for (int i = MyRowStart_; i <= MyRowEnd_; i++) {
    GlobalRowIDs[i-MyRowStart_] = i;
  }
  
  // Create array of global column ids
  int new_value = 0; int entries = 0;
  std::set<int> Columns;
  int num_entries;
  double *values;
  int *indices;
  for(int i = 0; i < NumMyRows_; i++){
    ierr += HYPRE_ParCSRMatrixGetRow(ParMatrix_, i+MyRowStart_, &num_entries, &indices, &values);
    ierr += HYPRE_ParCSRMatrixRestoreRow(ParMatrix_, i+MyRowStart_,&num_entries,&indices,&values);
    TEST_FOR_EXCEPTION(ierr != 0, std::logic_error, "Couldn't get row of matrix.");
    entries = num_entries;
    for(int j = 0; j < num_entries; j++){
      // Insert column ids from this row into set
      new_value = indices[j];
      Columns.insert(new_value);
    }
  }
  int NumMyCols = Columns.size(); 
  Teuchos::Array<int> GlobalColIDs; GlobalColIDs.resize(NumMyCols);
  
  std::set<int>::iterator it;
  int counter = 0;
  for (it = Columns.begin(); it != Columns.end(); it++) {
    // Get column ids in order
    GlobalColIDs[counter] = *it;
    counter = counter + 1;
  }
  //printf("Proc[%d] Rows from %d to %d, num = %d\n", Comm().MyPID(), MyRowStart_,MyRowEnd_, NumMyRows_);
  
  Epetra_Map RowMap(-1, NumMyRows_, &GlobalRowIDs[0], 0, Comm());
  Epetra_Map ColMap(-1, NumMyCols, &GlobalColIDs[0], 0, Comm());
  
  //Need to call SetMaps()
  SetMaps(RowMap, ColMap);
 
  // Get an MPI_Comm to create vectors.
  // The vectors will be reused in Multiply(), so that they aren't recreated every time.   
  MPI_Comm comm;
  ierr += HYPRE_ParCSRMatrixGetComm(ParMatrix_, &comm);
  TEST_FOR_EXCEPTION(ierr != 0, std::logic_error, "Couldn't get communicator from Hypre Matrix.");
  
  ierr += HYPRE_ParCSRMatrixGetComm(ParMatrix_, &comm);
  ierr += HYPRE_IJVectorCreate(comm, MyRowStart_, MyRowEnd_, &X_hypre);
  ierr += HYPRE_IJVectorSetObjectType(X_hypre, HYPRE_PARCSR);
  ierr += HYPRE_IJVectorInitialize(X_hypre);
  ierr += HYPRE_IJVectorAssemble(X_hypre);
  ierr += HYPRE_IJVectorGetObject(X_hypre, (void**) &par_x);
  TEST_FOR_EXCEPTION(ierr != 0, std::logic_error, "Couldn't create Hypre X vector.");

  ierr += HYPRE_IJVectorCreate(comm, MyRowStart_, MyRowEnd_, &Y_hypre);
  ierr += HYPRE_IJVectorSetObjectType(Y_hypre, HYPRE_PARCSR);
  ierr += HYPRE_IJVectorInitialize(Y_hypre);
  ierr += HYPRE_IJVectorAssemble(Y_hypre);
  ierr += HYPRE_IJVectorGetObject(Y_hypre, (void**) &par_y);
  TEST_FOR_EXCEPTION(ierr != 0, std::logic_error, "Couldn't create Hypre Y vector.");

  x_vec = (hypre_ParVector *) hypre_IJVectorObject(((hypre_IJVector *) X_hypre));
  x_local = hypre_ParVectorLocalVector(x_vec);

  y_vec = (hypre_ParVector *) hypre_IJVectorObject(((hypre_IJVector *) Y_hypre));
  y_local = hypre_ParVectorLocalVector(y_vec);


} //EpetraExt_HYPREIJMatrix(Hypre_IJMatrix) Constructor

//=======================================================

EpetraExt_HypreIJMatrix::~EpetraExt_HypreIJMatrix(){
  HYPRE_IJVectorDestroy(X_hypre);
  HYPRE_IJVectorDestroy(Y_hypre);

  /* Destroy solver and preconditioner */
  if(SolverCreated_){
    SolverDestroyPtr_(Solver_);
  } 
  if(PrecondCreated_){
    PrecondDestroyPtr_(Preconditioner_);
  } 
} //EpetraExt_HypreIJMatrix destructor

//=======================================================

int EpetraExt_HypreIJMatrix::ExtractMyRowCopy(int Row, int Length, int & NumEntries, 
                     double * Values, int * Indices) const 
{
  int ierr = 0;
  // Get values and indices of ith row of matrix
  int *indices;
  double *values;
  int num_entries;
  ierr += HYPRE_ParCSRMatrixGetRow(ParMatrix_, Row+MyRowStart_, &num_entries, &indices, &values);
  TEST_FOR_EXCEPTION(ierr != 0, std::logic_error, "Couldn't get row from Hypre Matrix.");
  ierr += HYPRE_ParCSRMatrixRestoreRow(ParMatrix_, Row+MyRowStart_, &num_entries, &indices, &values);
  TEST_FOR_EXCEPTION(ierr != 0, std::logic_error, "Couldn't restore row from Hypre Matrix.");
  
  NumEntries = num_entries;
  
  if(Length < NumEntries){
    printf("The arrays passed in are not large enough. Allocate more space.\n");
    return -2;
  }

  for(int i = 0; i < NumEntries; i++){
    Values[i] = values[i];
    Indices[i] = RowMatrixColMap().LID(indices[i]);
  }

  return ierr;
}

//=======================================================

int EpetraExt_HypreIJMatrix::NumMyRowEntries(int Row, int & NumEntries) const 
{
  int ierr = 0;
  int my_row[1];
  my_row[0] = Row+MyRowStart_;
  int nentries[1];
  ierr += HYPRE_IJMatrixGetRowCounts(Matrix_, 1, my_row, nentries);
  TEST_FOR_EXCEPTION(ierr != 0, std::logic_error, "Couldn't get row count from Hypre Matrix.");
  NumEntries = nentries[0];
  return ierr;
}
//=======================================================

 int EpetraExt_HypreIJMatrix::ExtractMyEntryView(int CurEntry, double *&Value, int &RowIndex, int &ColIndex)
{/* 
  This gives a copy, not a view of the values, so is not implemented.
  if(CurEntry >= NumMyNonzeros() || CurEntry < 0){
    return -1;
  }
  int counter = 0;
  for(int i = 0; i < NumMyRows(); i++){
    int entries;
    NumMyRowEntries(i, entries);
    counter += entries;
    if(counter > CurEntry) {
      counter -= entries;
      RowIndex = i;
      break;
    }
  }
  int entries;
  NumMyRowEntries(RowIndex, entries);
  int *indices;
  double *values;
  int num_entries;
  HYPRE_ParCSRMatrixGetRow(ParMatrix_, RowIndex+MyRowStart_, &num_entries, &indices, &values);
  HYPRE_ParCSRMatrixRestoreRow(ParMatrix_, RowIndex+MyRowStart_, &num_entries, &indices, &values);
  ColIndex = RowMatrixColMap().LID(indices[CurEntry-counter]);
  Value = &values[CurEntry-counter];
  return 0;*/return -1;
}
 
 //=======================================================
 int EpetraExt_HypreIJMatrix::ExtractMyEntryView(int CurEntry, double const * & Value, int & RowIndex, int & ColIndex) const 
{/*
   if(CurEntry >= NumMyNonzeros() || CurEntry < 0){
     return -1;
   }
   int counter = 0;
   for(int i = 0; i < NumMyRows(); i++){
     int entries;
     NumMyRowEntries(i, entries);
     counter += entries;
     if(counter > CurEntry) {
       counter -= entries;
       RowIndex = i;
       break;
     }
   }
   int entries;
   NumMyRowEntries(RowIndex, entries);
  int *indices;
  double *values;
  int num_entries;
  HYPRE_ParCSRMatrixGetRow(ParMatrix_, RowIndex+MyRowStart_, &num_entries, &indices, &values);
  HYPRE_ParCSRMatrixRestoreRow(ParMatrix_, RowIndex+MyRowStart_, &num_entries, &indices, &values);
  ColIndex = RowMatrixColMap().LID(indices[CurEntry-counter]);
  Value = &values[CurEntry-counter];
  return 0;*/ return -1;
}

//=======================================================
int EpetraExt_HypreIJMatrix::Multiply(bool TransA,
                               const Epetra_MultiVector& X,
                               Epetra_MultiVector& Y) const
{
  
  //printf("Proc[%d], Row start: %d, Row End: %d\n", Comm().MyPID(), MyRowStart_, MyRowEnd_);
  int ierr = 0;
  bool SameVectors = false; 
  int NumVectors = X.NumVectors();
  if (NumVectors != Y.NumVectors()) return -1;  // X and Y must have same number of vectors
  if(X.Pointers() == Y.Pointers()){
    SameVectors = true;
  }
  for(int VecNum = 0; VecNum < NumVectors; VecNum++) {
    //Get values for current vector in multivector.
    double * x_values;
    double * y_values;
    ierr += (*X(VecNum)).ExtractView(&x_values);
    double *x_temp = x_local->data; 
    double *y_temp = y_local->data;
    if(!SameVectors){
      ierr += (*Y(VecNum)).ExtractView(&y_values);
    } else {
      y_values = new double[X.MyLength()];
    }
    y_local->data = y_values;
    HYPRE_ParVectorSetConstantValues(par_y,0.0);
    // Temporarily make a pointer to data in Hypre for end
    // Replace data in Hypre vectors with epetra values
    x_local->data = x_values;
    // Do actual computation.
    if(TransA) {
      // Use transpose of A in multiply
      ierr += HYPRE_ParCSRMatrixMatvecT(1.0, ParMatrix_, par_x, 1.0, par_y);
      TEST_FOR_EXCEPTION(ierr != 0, std::logic_error, "Couldn't do Transpose MatVec in Hypre.");
    } else {
      ierr += HYPRE_ParCSRMatrixMatvec(1.0, ParMatrix_, par_x, 1.0, par_y);
      TEST_FOR_EXCEPTION(ierr != 0, std::logic_error, "Couldn't do MatVec in Hypre.");
    }
    if(SameVectors){
      /*double *invalid_ptr;
      (*Y(VecNum)).ExtractView(&invalid_ptr);
      (*Y(VecNum)).ResetView(y_values);
      delete[] invalid_ptr;*/
      
      int NumEntries = Y.MyLength();
      std::vector<double> new_values; new_values.resize(NumEntries);
      std::vector<int> new_indices; new_indices.resize(NumEntries);
      for(int i = 0; i < NumEntries; i++){
        new_values[i] = y_values[i];
        new_indices[i] = i;
      }
      (*Y(VecNum)).ReplaceMyValues(NumEntries, &new_values[0], &new_indices[0]);
      delete[] y_values;
    }
    x_local->data = x_temp;
    y_local->data = y_temp;
  }

  double flops = (double) NumVectors * (double) NumGlobalNonzeros();
  UpdateFlops(flops);
  return(ierr);
} //Multiply() 

int EpetraExt_HypreIJMatrix::SetParameter(Hypre_Chooser chooser, int parameter, int (*pt2Func)(HYPRE_Solver, int)){
  int result;
  if(chooser){
    result = pt2Func(Preconditioner_, parameter);
  }  else {
    result = pt2Func(Solver_, parameter);
  }
  return result;
}

int EpetraExt_HypreIJMatrix::SetParameter(Hypre_Chooser chooser, double parameter, int (*pt2Func)(HYPRE_Solver, double)){
  int result;
  if(chooser){
    result = pt2Func(Preconditioner_, parameter);
  } else {
    result = pt2Func(Solver_, parameter);
  }
  return result;
}

int EpetraExt_HypreIJMatrix::SetParameter(Hypre_Chooser chooser, double parameter1, int parameter2, int (*pt2Func)(HYPRE_Solver, double, int)){
  int result;
  if(chooser){
    result = pt2Func(Preconditioner_, parameter1, parameter2);
  } else {
    result = pt2Func(Solver_, parameter1, parameter2);
  }
  return result;
}

int EpetraExt_HypreIJMatrix::SetParameter(Hypre_Chooser chooser, int parameter1, int parameter2, int (*pt2Func)(HYPRE_Solver, int, int)){
  int result;
  if(chooser){
    result = pt2Func(Preconditioner_, parameter1, parameter2);
  } else {
    result = pt2Func(Solver_, parameter1, parameter2);
  }
  return result;
}

int EpetraExt_HypreIJMatrix::SetParameter(Hypre_Chooser chooser, double* parameter, int (*pt2Func)(HYPRE_Solver, double*)){
  int result;
  if(chooser){
    result = pt2Func(Preconditioner_, parameter);
  } else {
    result = pt2Func(Solver_, parameter);
  }
  return result;
}

int EpetraExt_HypreIJMatrix::SetParameter(Hypre_Chooser chooser, int* parameter, int (*pt2Func)(HYPRE_Solver, int*)){
  int result;
  if(chooser){
    result = pt2Func(Preconditioner_, parameter);
  } else {
    result = pt2Func(Solver_, parameter);
  }
  return result;
}

int EpetraExt_HypreIJMatrix::SetSolverType(Hypre_Solver Solver, bool transpose){
  if(transpose && Solver != BoomerAMG){
    return -1;
  }
  switch(Solver) {
    case BoomerAMG:
      SolverCreatePtr_ = &EpetraExt_HypreIJMatrix::Hypre_BoomerAMGCreate;
      SolverDestroyPtr_ = &HYPRE_BoomerAMGDestroy;
      SolverSetupPtr_ = &HYPRE_BoomerAMGSetup;
      SolverPrecondPtr_ = NULL;
      if(transpose){
        TransposeSolve_ = true;
        SolverSolvePtr_ = &HYPRE_BoomerAMGSolveT;
      } else {
        SolverSolvePtr_ = &HYPRE_BoomerAMGSolve;
      }
      break;
    case AMS:
      SolverCreatePtr_ = &EpetraExt_HypreIJMatrix::Hypre_AMSCreate;
      SolverDestroyPtr_ = &HYPRE_AMSDestroy;
      SolverSetupPtr_ = &HYPRE_AMSSetup;
      SolverSolvePtr_ = &HYPRE_AMSSolve;
      SolverPrecondPtr_ = NULL;
      break;
    case Hybrid:
      SolverCreatePtr_ = &EpetraExt_HypreIJMatrix::Hypre_ParCSRHybridCreate;
      SolverDestroyPtr_ = &HYPRE_ParCSRHybridDestroy;
      SolverSetupPtr_ = &HYPRE_ParCSRHybridSetup;
      SolverSolvePtr_ = &HYPRE_ParCSRHybridSolve;
      SolverPrecondPtr_ = &HYPRE_ParCSRHybridSetPrecond;
      break;
    case PCG:
      SolverCreatePtr_ = &EpetraExt_HypreIJMatrix::Hypre_ParCSRPCGCreate;
      SolverDestroyPtr_ = &HYPRE_ParCSRPCGDestroy;
      SolverSetupPtr_ = &HYPRE_ParCSRPCGSetup;
      SolverSolvePtr_ = &HYPRE_ParCSRPCGSolve;
      SolverPrecondPtr_ = &HYPRE_ParCSRPCGSetPrecond;
      break;
    case GMRES:
      SolverCreatePtr_ = &EpetraExt_HypreIJMatrix::Hypre_ParCSRGMRESCreate;
      SolverDestroyPtr_ = &HYPRE_ParCSRGMRESDestroy;
      SolverSetupPtr_ = &HYPRE_ParCSRGMRESSetup;
      SolverSolvePtr_ = &HYPRE_ParCSRGMRESSolve;
      SolverPrecondPtr_ = &HYPRE_ParCSRGMRESSetPrecond;
      break;
    case FlexGMRES:
      SolverCreatePtr_ = &EpetraExt_HypreIJMatrix::Hypre_ParCSRFlexGMRESCreate;
      SolverDestroyPtr_ = &HYPRE_ParCSRFlexGMRESDestroy;
      SolverSetupPtr_ = &HYPRE_ParCSRFlexGMRESSetup;
      SolverSolvePtr_ = &HYPRE_ParCSRFlexGMRESSolve;
      SolverPrecondPtr_ = &HYPRE_ParCSRFlexGMRESSetPrecond;
      break;
    case LGMRES:
      SolverCreatePtr_ = &EpetraExt_HypreIJMatrix::Hypre_ParCSRLGMRESCreate;
      SolverDestroyPtr_ = &HYPRE_ParCSRLGMRESDestroy;
      SolverSetupPtr_ = &HYPRE_ParCSRLGMRESSetup;
      SolverSolvePtr_ = &HYPRE_ParCSRLGMRESSolve;
      SolverPrecondPtr_ = &HYPRE_ParCSRLGMRESSetPrecond;
      break;
    case BiCGSTAB:
      SolverCreatePtr_ = &EpetraExt_HypreIJMatrix::Hypre_ParCSRBiCGSTABCreate;
      SolverDestroyPtr_ = &HYPRE_ParCSRBiCGSTABDestroy;
      SolverSetupPtr_ = &HYPRE_ParCSRBiCGSTABSetup;
      SolverSolvePtr_ = &HYPRE_ParCSRBiCGSTABSolve;
      SolverPrecondPtr_ = &HYPRE_ParCSRBiCGSTABSetPrecond;
      break;
    default:
      return -1;
    }
  CreateSolver();
  return 0;
}

int EpetraExt_HypreIJMatrix::SetPrecondType(Hypre_Solver Precond){
  switch(Precond) {
    case BoomerAMG:
      PrecondCreatePtr_ = &EpetraExt_HypreIJMatrix::Hypre_BoomerAMGCreate;
      PrecondDestroyPtr_ = &HYPRE_BoomerAMGDestroy;
      PrecondSetupPtr_ = &HYPRE_BoomerAMGSetup;
      PrecondSolvePtr_ = &HYPRE_BoomerAMGSolve;
      break;
    case ParaSails:
      PrecondCreatePtr_ = &EpetraExt_HypreIJMatrix::Hypre_ParaSailsCreate;
      PrecondDestroyPtr_ = &HYPRE_ParaSailsDestroy;
      PrecondSetupPtr_ = &HYPRE_ParaSailsSetup;
      PrecondSolvePtr_ = &HYPRE_ParaSailsSolve;
      break;
    case Euclid:
      PrecondCreatePtr_ = &EpetraExt_HypreIJMatrix::Hypre_EuclidCreate;
      PrecondDestroyPtr_ = &HYPRE_EuclidDestroy;
      PrecondSetupPtr_ = &HYPRE_EuclidSetup;
      PrecondSolvePtr_ = &HYPRE_EuclidSolve;
      break;
    case Pilut:
      PrecondCreatePtr_ = &EpetraExt_HypreIJMatrix::Hypre_ParCSRPilutCreate;
      PrecondDestroyPtr_ = &HYPRE_ParCSRPilutDestroy;
      PrecondSetupPtr_ = &HYPRE_ParCSRPilutSetup;
      PrecondSolvePtr_ = &HYPRE_ParCSRPilutSolve;
      break;
    case AMS:
      PrecondCreatePtr_ = &EpetraExt_HypreIJMatrix::Hypre_AMSCreate;
      PrecondDestroyPtr_ = &HYPRE_AMSDestroy;
      PrecondSetupPtr_ = &HYPRE_AMSSetup;
      PrecondSolvePtr_ = &HYPRE_AMSSolve;
      break;
    default:
      return -1;
    }
  CreatePrecond();
  return 0;

}

int EpetraExt_HypreIJMatrix::CreateSolver(){
  SolverCreated_ = true;
  MPI_Comm comm;
  HYPRE_ParCSRMatrixGetComm(ParMatrix_, &comm);
  return (this->*SolverCreatePtr_)(comm, &Solver_);
}

int EpetraExt_HypreIJMatrix::CreatePrecond(){
  PrecondCreated_ = true;
  MPI_Comm comm;
  HYPRE_ParCSRMatrixGetComm(ParMatrix_, &comm);
  return (this->*PrecondCreatePtr_)(comm, &Preconditioner_);
}

int EpetraExt_HypreIJMatrix::SetPreconditioner(){
  if(SolverPrecondPtr_ == NULL){
    return -1;
  } else {
    SolverPrecondPtr_(Solver_, PrecondSolvePtr_, PrecondSetupPtr_, Preconditioner_);
    return 0;
  }
}

int EpetraExt_HypreIJMatrix::SolveOrPrecondition(Hypre_Chooser answer){
  switch(answer) {
    case Solver:
      if(SolverCreated_){
        SolveOrPrec_ = answer;
        return 0;
      } else {
        return -1;
      }
    case Preconditioner:
      if(PrecondCreated_){
        SolveOrPrec_ = answer;
        return 0;
      } else {
        return -1;
      }
  }
  return -1;
}

int EpetraExt_HypreIJMatrix::Solve(bool Upper, bool transpose, bool UnitDiagonal, const Epetra_MultiVector & X, Epetra_MultiVector & Y) const {
  int ierr = 0;
  bool SameVectors = false;
  int NumVectors = X.NumVectors();
  if (NumVectors != Y.NumVectors()) return -1;  // X and Y must have same number of vectors
  if(X.Pointers() == Y.Pointers()){
    SameVectors = true;
  }
     
  for(int VecNum = 0; VecNum < NumVectors; VecNum++) {
    //Get values for current vector in multivector.
    double * x_values;
    ierr += (*X(VecNum)).ExtractView(&x_values);
    double * y_values;
    if(!SameVectors){
      ierr += (*Y(VecNum)).ExtractView(&y_values);
    } else {
      y_values = new double[X.MyLength()];
    }
    // Temporarily make a pointer to data in Hypre for end
    double *x_temp = x_local->data; 
    // Replace data in Hypre vectors with epetra values
    x_local->data = x_values;
    double *y_temp = y_local->data;
    y_local->data = y_values;
    
    HYPRE_ParVectorSetConstantValues(par_y, 0.0);
    if(!SolverCreated_ && !PrecondCreated_){
      // Need to create a solver or preconditioner
      return -1;
    }
    if(transpose && !TransposeSolve_){
      // User requested a transpose solve, but the solver selected doesn't provide one
      return -1;
    }
    if(SolveOrPrec_ == Solver){
      // Use the solver methods
      ierr += SolverSetupPtr_(Solver_, ParMatrix_, par_x, par_y);
      ierr += SolverSolvePtr_(Solver_, ParMatrix_, par_x, par_y);
    } else {
      // Apply the preconditioner
      ierr += PrecondSetupPtr_(Preconditioner_, ParMatrix_, par_x, par_y);
      ierr += PrecondSolvePtr_(Preconditioner_, ParMatrix_, par_x, par_y);
    }
    TEST_FOR_EXCEPTION(ierr != 0, std::logic_error, "Hypre solve failed.");
    
    if(SameVectors){
      int NumEntries = Y.MyLength();
      std::vector<double> new_values; new_values.resize(NumEntries);
      std::vector<int> new_indices; new_indices.resize(NumEntries);
      for(int i = 0; i < NumEntries; i++){
        new_values[i] = y_values[i];
        new_indices[i] = i;
      }
      (*Y(VecNum)).ReplaceMyValues(NumEntries, &new_values[0], &new_indices[0]);
      delete[] y_values;
    }
    x_local->data = x_temp;
    y_local->data = y_temp;
  }

  double flops = (double) NumVectors * (double) NumGlobalNonzeros();
  UpdateFlops(flops);
  return(ierr);
}

//=======================================================
int EpetraExt_HypreIJMatrix::LeftScale(const Epetra_Vector& X) {
  int ierr = 0;
  
  for(int i = 0; i < NumMyRows_; i++){
    //Vector-scalar mult on ith row
    int num_entries;
    int *indices;
    double *values;
    ierr += HYPRE_ParCSRMatrixGetRow(ParMatrix_,i+MyRowStart_, &num_entries, &indices, &values);
    ierr += HYPRE_ParCSRMatrixRestoreRow(ParMatrix_,i+MyRowStart_, &num_entries, &indices, &values);
    Teuchos::Array<double> new_values; new_values.resize(num_entries);
    Teuchos::Array<int> new_indices; new_indices.resize(num_entries);
    for(int j = 0; j < num_entries; j++){
      // Scale this row with the appropriate values from the vector
      new_values[j] = X[i]*values[j];
      new_indices[j] = indices[j];
    }
    int rows[1];
    rows[0] = i+MyRowStart_;
    ierr += HYPRE_IJMatrixSetValues(Matrix_, 1, &num_entries, rows, &new_indices[0], &new_values[0]);
    // Finally set values of the Matrix for this row
    
    TEST_FOR_EXCEPTION(ierr != 0, std::logic_error, "Couldn't set values in Hypre Matrix.");
  }
  
  HaveNumericConstants_ = false;
  UpdateFlops(NumGlobalNonzeros());
  return(ierr);
} //LeftScale()

//=======================================================
int EpetraExt_HypreIJMatrix::RightScale(const Epetra_Vector& X) {
  int ierr = 0;
  // First we need to import off-processor values of the vector
  Epetra_Import Importer(RowMatrixColMap(), RowMatrixRowMap());
  Epetra_Vector Import_Vector(RowMatrixColMap(), true);
  ierr += Import_Vector.Import(X, Importer, Insert, 0);
  TEST_FOR_EXCEPTION(ierr != 0, std::logic_error, "Couldn't import values ierr = " << ierr << ".");
  
  for(int i = 0; i < NumMyRows_; i++){
    //Vector-scalar mult on ith column
    int num_entries;
    double *values;
    int *indices;
    // Get values and indices of ith row of matrix
    ierr += HYPRE_ParCSRMatrixGetRow(ParMatrix_,i+MyRowStart_, &num_entries, &indices, &values);
    ierr += HYPRE_ParCSRMatrixRestoreRow(ParMatrix_,i+MyRowStart_,&num_entries,&indices,&values);
    TEST_FOR_EXCEPTION(ierr != 0, std::logic_error, "Couldn't get row.");
    Teuchos::Array<int> new_indices; new_indices.resize(num_entries);
    Teuchos::Array<double> new_values; new_values.resize(num_entries);
    for(int j = 0; j < num_entries; j++){
      // Multiply column j with jth element
      int index = RowMatrixColMap().LID(indices[j]);
      TEST_FOR_EXCEPTION(index < 0, std::logic_error, "Index is negtive.");
      new_values[j] = values[j] * Import_Vector[index];
      new_indices[j] = indices[j];
    }
    // Finally set values of the Matrix for this row
    int rows[1];
    rows[0] = i+MyRowStart_;
    ierr += HYPRE_IJMatrixSetValues(Matrix_, 1, &num_entries, rows, &new_indices[0], &new_values[0]);
    TEST_FOR_EXCEPTION(ierr != 0, std::logic_error, "Couldn't set values in Hypre Matrix.");
  }
  
  HaveNumericConstants_ = false;
  UpdateFlops(NumGlobalNonzeros());
  return(ierr);
} //RightScale()

int EpetraExt_HypreIJMatrix::Hypre_BoomerAMGCreate(MPI_Comm comm, HYPRE_Solver *solver){
  return HYPRE_BoomerAMGCreate(solver);
}

int EpetraExt_HypreIJMatrix::Hypre_ParaSailsCreate(MPI_Comm comm, HYPRE_Solver *solver){
  return HYPRE_ParaSailsCreate(comm, solver);
}

int EpetraExt_HypreIJMatrix::Hypre_EuclidCreate(MPI_Comm comm, HYPRE_Solver *solver){
  return HYPRE_EuclidCreate(comm, solver);
}

int EpetraExt_HypreIJMatrix::Hypre_ParCSRPilutCreate(MPI_Comm comm, HYPRE_Solver *solver){
  return HYPRE_ParCSRPilutCreate(comm, solver);
}

int EpetraExt_HypreIJMatrix::Hypre_AMSCreate(MPI_Comm comm, HYPRE_Solver *solver){
  return HYPRE_AMSCreate(solver);
}

int EpetraExt_HypreIJMatrix::Hypre_ParCSRHybridCreate(MPI_Comm comm, HYPRE_Solver *solver){
  return HYPRE_ParCSRHybridCreate(solver);
}

int EpetraExt_HypreIJMatrix::Hypre_ParCSRPCGCreate(MPI_Comm comm, HYPRE_Solver *solver){
  return HYPRE_ParCSRPCGCreate(comm, solver);
}

int EpetraExt_HypreIJMatrix::Hypre_ParCSRGMRESCreate(MPI_Comm comm, HYPRE_Solver *solver){
  return HYPRE_ParCSRGMRESCreate(comm, solver);
}

int EpetraExt_HypreIJMatrix::Hypre_ParCSRFlexGMRESCreate(MPI_Comm comm, HYPRE_Solver *solver){
  return HYPRE_ParCSRFlexGMRESCreate(comm, solver);
}

int EpetraExt_HypreIJMatrix::Hypre_ParCSRLGMRESCreate(MPI_Comm comm, HYPRE_Solver *solver){
  return HYPRE_ParCSRLGMRESCreate(comm, solver);
}

int EpetraExt_HypreIJMatrix::Hypre_ParCSRBiCGSTABCreate(MPI_Comm comm, HYPRE_Solver *solver){
  return HYPRE_ParCSRBiCGSTABCreate(comm, solver);
}

int EpetraExt_HypreIJMatrix::InitializeDefaults(){

  int ierr = 0;
  
  int my_type;
  // Get type of Hypre IJ Matrix
  ierr += HYPRE_IJMatrixGetObjectType(Matrix_, &my_type);
  TEST_FOR_EXCEPTION(ierr != 0, std::logic_error, "Couldn't get object type of Hypre Matrix.");
  MatType_ = my_type;
  // Currently only ParCSR is supported
  TEST_FOR_EXCEPTION(MatType_ != HYPRE_PARCSR, std::logic_error, "Object is not type PARCSR");

  // Get the actual ParCSR object from the IJ matrix
  ierr += HYPRE_IJMatrixGetObject(Matrix_, (void**) &ParMatrix_);
  TEST_FOR_EXCEPTION(ierr != 0, std::logic_error, "Couldn't get Hypre_ParCSR Matrix object.");
  
  int numRows, numCols;
  
  // Get dimensions of the matrix and store as global variables
  ierr +=HYPRE_ParCSRMatrixGetDims(ParMatrix_, &numRows, &numCols);
  TEST_FOR_EXCEPTION(ierr != 0, std::logic_error, "Couldn't get dimensions of Hypre Matrix.");
  
  NumGlobalRows_ = numRows;
  NumGlobalCols_ = numCols;
  
  // Get the local dimensions of the matrix
  int ColStart, ColEnd;
  ierr += HYPRE_ParCSRMatrixGetLocalRange(ParMatrix_, &MyRowStart_, &MyRowEnd_, &ColStart, &ColEnd);
  TEST_FOR_EXCEPTION(ierr != 0, std::logic_error, "Couldn't get local range of Hypre Matrix.");
  
  // Determine number of local rows
  NumMyRows_ = MyRowEnd_ - MyRowStart_+1;
  
  return ierr;
}  //InitializeDefaults() 

#endif /*HAVE_HYPRE*/
