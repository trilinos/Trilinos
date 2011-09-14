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
  : Epetra_BasicRowMatrix(Epetra_MpiComm(hypre_IJMatrixComm(matrix))),
    Matrix_(matrix),
    ParMatrix_(0),
    NumMyRows_(-1),
    NumGlobalRows_(-1),
    NumGlobalCols_(-1),
    MyRowStart_(-1),
    MyRowEnd_(-1),
    MatType_(-1), 
    TransposeSolve_(false),
    SolveOrPrec_(Solver)
{
  IsSolverSetup_ = new bool[1];
  IsPrecondSetup_ = new bool[1];
  IsSolverSetup_[0] = false;
  IsPrecondSetup_[0] = false;
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

  SolverCreatePtr_ = &EpetraExt_HypreIJMatrix::Hypre_ParCSRPCGCreate;
  SolverDestroyPtr_ = &HYPRE_ParCSRPCGDestroy;
  SolverSetupPtr_ = &HYPRE_ParCSRPCGSetup;
  SolverSolvePtr_ = &HYPRE_ParCSRPCGSolve;
  SolverPrecondPtr_ = &HYPRE_ParCSRPCGSetPrecond;
  CreateSolver();

  PrecondCreatePtr_ = &EpetraExt_HypreIJMatrix::Hypre_EuclidCreate;
  PrecondDestroyPtr_ = &HYPRE_EuclidDestroy;
  PrecondSetupPtr_ = &HYPRE_EuclidSetup;
  PrecondSolvePtr_ = &HYPRE_EuclidSolve;
  CreatePrecond();
  ComputeNumericConstants();
  ComputeStructureConstants();
} //EpetraExt_HYPREIJMatrix(Hypre_IJMatrix) Constructor

//=======================================================
EpetraExt_HypreIJMatrix::~EpetraExt_HypreIJMatrix(){
  int ierr = 0;
  ierr += HYPRE_IJVectorDestroy(X_hypre);
  TEST_FOR_EXCEPTION(ierr != 0, std::logic_error, "Couldn't destroy the X Vector.");
  ierr += HYPRE_IJVectorDestroy(Y_hypre);
  TEST_FOR_EXCEPTION(ierr != 0, std::logic_error, "Couldn't destroy the Y Vector.");

  /* Destroy solver and preconditioner */
  if(IsSolverSetup_[0] == true){
    ierr += SolverDestroyPtr_(Solver_);
    TEST_FOR_EXCEPTION(ierr != 0, std::logic_error, "Couldn't destroy the Solver.");
  }
  if(IsPrecondSetup_[0] == true){
    ierr += PrecondDestroyPtr_(Preconditioner_);
    TEST_FOR_EXCEPTION(ierr != 0, std::logic_error, "Couldn't destroy the Preconditioner.");
  }
  delete[] IsSolverSetup_;
  delete[] IsPrecondSetup_;
} //EpetraExt_HypreIJMatrix destructor

//=======================================================
int EpetraExt_HypreIJMatrix::ExtractMyRowCopy(int Row, int Length, int & NumEntries, 
                     double * Values, int * Indices) const 
{
  // Get values and indices of ith row of matrix
  int *indices;
  double *values;
  int num_entries;
  EPETRA_CHK_ERR(HYPRE_ParCSRMatrixGetRow(ParMatrix_, Row+RowMatrixRowMap().MinMyGID(), &num_entries, &indices, &values));
  EPETRA_CHK_ERR(HYPRE_ParCSRMatrixRestoreRow(ParMatrix_, Row+RowMatrixRowMap().MinMyGID(), &num_entries, &indices, &values));
  
  NumEntries = num_entries;
  
  if(Length < NumEntries){
    printf("The arrays passed in are not large enough. Allocate more space.\n");
    return -2;
  }

  for(int i = 0; i < NumEntries; i++){
    Values[i] = values[i];
    Indices[i] = RowMatrixColMap().LID(indices[i]);
  }

  return 0;
} //ExtractMyRowCopy()

//=======================================================
int EpetraExt_HypreIJMatrix::NumMyRowEntries(int Row, int & NumEntries) const 
{
  int my_row[1];
  my_row[0] = Row+RowMatrixRowMap().MinMyGID();
  int nentries[1];
  EPETRA_CHK_ERR(HYPRE_IJMatrixGetRowCounts(Matrix_, 1, my_row, nentries));
  NumEntries = nentries[0];
  return 0;
} //NumMyRowEntries()

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
} //ExtractMyEntryView() - not implemented
 
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
} //ExtractMyEntryView() - const version, not implemented

//=======================================================
int EpetraExt_HypreIJMatrix::Multiply(bool TransA,
                               const Epetra_MultiVector& X,
                               Epetra_MultiVector& Y) const
{
  
  //printf("Proc[%d], Row start: %d, Row End: %d\n", Comm().MyPID(), MyRowStart_, MyRowEnd_);
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
    EPETRA_CHK_ERR((*X(VecNum)).ExtractView(&x_values));
    double *x_temp = x_local->data; 
    double *y_temp = y_local->data;
    if(!SameVectors){
      EPETRA_CHK_ERR((*Y(VecNum)).ExtractView(&y_values));
    } else {
      y_values = new double[X.MyLength()];
    }
    y_local->data = y_values;
    EPETRA_CHK_ERR(HYPRE_ParVectorSetConstantValues(par_y,0.0));
    // Temporarily make a pointer to data in Hypre for end
    // Replace data in Hypre vectors with epetra values
    x_local->data = x_values;
    // Do actual computation.
    if(TransA) {
      // Use transpose of A in multiply
      EPETRA_CHK_ERR(HYPRE_ParCSRMatrixMatvecT(1.0, ParMatrix_, par_x, 1.0, par_y));
    } else {
      EPETRA_CHK_ERR(HYPRE_ParCSRMatrixMatvec(1.0, ParMatrix_, par_x, 1.0, par_y));
    }
    if(SameVectors){
      int NumEntries = Y.MyLength();
      std::vector<double> new_values; new_values.resize(NumEntries);
      std::vector<int> new_indices; new_indices.resize(NumEntries);
      for(int i = 0; i < NumEntries; i++){
        new_values[i] = y_values[i];
        new_indices[i] = i;
      }
      EPETRA_CHK_ERR((*Y(VecNum)).ReplaceMyValues(NumEntries, &new_values[0], &new_indices[0]));
      delete[] y_values;
    }
    x_local->data = x_temp;
    y_local->data = y_temp;
  }
  double flops = (double) NumVectors * (double) NumGlobalNonzeros();
  UpdateFlops(flops);
  return 0;
} //Multiply() 

//=======================================================
int EpetraExt_HypreIJMatrix::SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, int), int parameter){
  if(chooser == Preconditioner){
    EPETRA_CHK_ERR(pt2Func(Preconditioner_, parameter));
  }  else {
    EPETRA_CHK_ERR(pt2Func(Solver_, parameter));
  }
  return 0;
} //SetParameter() - int function pointer

//=======================================================
int EpetraExt_HypreIJMatrix::SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, double), double parameter){
  if(chooser == Preconditioner){
    EPETRA_CHK_ERR(pt2Func(Preconditioner_, parameter));
  } else {
    EPETRA_CHK_ERR(pt2Func(Solver_, parameter));
  }
  return 0;
} //SetParamater() - double function pointer

//=======================================================
int EpetraExt_HypreIJMatrix::SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, double, int), double parameter1, int parameter2){
  if(chooser == Preconditioner){
    EPETRA_CHK_ERR(pt2Func(Preconditioner_, parameter1, parameter2));
  } else {
    EPETRA_CHK_ERR(pt2Func(Solver_, parameter1, parameter2));
  }
  return 0;
} //SetParameter() - double,int function pointer

//=======================================================
int EpetraExt_HypreIJMatrix::SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, int, int), int parameter1, int parameter2){
  if(chooser == Preconditioner){
    EPETRA_CHK_ERR(pt2Func(Preconditioner_, parameter1, parameter2));
  } else {
    EPETRA_CHK_ERR(pt2Func(Solver_, parameter1, parameter2));
  }
  return 0;
} //SetParameter() - int,int function pointer

//=======================================================
int EpetraExt_HypreIJMatrix::SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, double*), double* parameter){
  if(chooser == Preconditioner){
    EPETRA_CHK_ERR(pt2Func(Preconditioner_, parameter));
  } else {
    EPETRA_CHK_ERR(pt2Func(Solver_, parameter));
  }
  return 0;
} //SetParameter() - double* function pointer

//=======================================================
int EpetraExt_HypreIJMatrix::SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, int*), int* parameter){
  if(chooser == Preconditioner){
    EPETRA_CHK_ERR(pt2Func(Preconditioner_, parameter));
  } else {
    EPETRA_CHK_ERR(pt2Func(Solver_, parameter));
  }
  return 0;
} //SetParameter() - int* function pointer

//=======================================================
int EpetraExt_HypreIJMatrix::SetParameter(Hypre_Chooser chooser, Hypre_Solver solver, bool transpose){
  if(chooser == Solver){
  if(transpose && solver != BoomerAMG){
    EPETRA_CHK_ERR(-1);
  }
  switch(solver) {
    case BoomerAMG:
      if(IsSolverSetup_[0]){
        SolverDestroyPtr_(Solver_);
        IsSolverSetup_[0] = false;
      }
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
      if(IsSolverSetup_[0]){
        SolverDestroyPtr_(Solver_);
        IsSolverSetup_[0] = false;
      }
      SolverCreatePtr_ = &EpetraExt_HypreIJMatrix::Hypre_AMSCreate;
      SolverDestroyPtr_ = &HYPRE_AMSDestroy;
      SolverSetupPtr_ = &HYPRE_AMSSetup;
      SolverSolvePtr_ = &HYPRE_AMSSolve;
      SolverPrecondPtr_ = NULL;
      break;
    case Hybrid:
      if(IsSolverSetup_[0]){
        SolverDestroyPtr_(Solver_);
        IsSolverSetup_[0] = false;
      }
      SolverCreatePtr_ = &EpetraExt_HypreIJMatrix::Hypre_ParCSRHybridCreate;
      SolverDestroyPtr_ = &HYPRE_ParCSRHybridDestroy;
      SolverSetupPtr_ = &HYPRE_ParCSRHybridSetup;
      SolverSolvePtr_ = &HYPRE_ParCSRHybridSolve;
      SolverPrecondPtr_ = &HYPRE_ParCSRHybridSetPrecond;
      break;
    case PCG:
      if(IsSolverSetup_[0]){
        SolverDestroyPtr_(Solver_);
        IsSolverSetup_[0] = false;
      }
      SolverCreatePtr_ = &EpetraExt_HypreIJMatrix::Hypre_ParCSRPCGCreate;
      SolverDestroyPtr_ = &HYPRE_ParCSRPCGDestroy;
      SolverSetupPtr_ = &HYPRE_ParCSRPCGSetup;
      SolverSolvePtr_ = &HYPRE_ParCSRPCGSolve;
      SolverPrecondPtr_ = &HYPRE_ParCSRPCGSetPrecond;
      break;
    case GMRES:
      if(IsSolverSetup_[0]){
        SolverDestroyPtr_(Solver_);
        IsSolverSetup_[0] = false;
      }
      SolverCreatePtr_ = &EpetraExt_HypreIJMatrix::Hypre_ParCSRGMRESCreate;
      SolverDestroyPtr_ = &HYPRE_ParCSRGMRESDestroy;
      SolverSetupPtr_ = &HYPRE_ParCSRGMRESSetup;
      SolverSolvePtr_ = &HYPRE_ParCSRGMRESSolve;
      SolverPrecondPtr_ = &HYPRE_ParCSRGMRESSetPrecond;
      break;
    case FlexGMRES:
      if(IsSolverSetup_[0]){
        SolverDestroyPtr_(Solver_);
        IsSolverSetup_[0] = false;
      }
      SolverCreatePtr_ = &EpetraExt_HypreIJMatrix::Hypre_ParCSRFlexGMRESCreate;
      SolverDestroyPtr_ = &HYPRE_ParCSRFlexGMRESDestroy;
      SolverSetupPtr_ = &HYPRE_ParCSRFlexGMRESSetup;
      SolverSolvePtr_ = &HYPRE_ParCSRFlexGMRESSolve;
      SolverPrecondPtr_ = &HYPRE_ParCSRFlexGMRESSetPrecond;
      break;
    case LGMRES:
      if(IsSolverSetup_[0]){
        SolverDestroyPtr_(Solver_);
        IsSolverSetup_[0] = false;
      }
      SolverCreatePtr_ = &EpetraExt_HypreIJMatrix::Hypre_ParCSRLGMRESCreate;
      SolverDestroyPtr_ = &HYPRE_ParCSRLGMRESDestroy;
      SolverSetupPtr_ = &HYPRE_ParCSRLGMRESSetup;
      SolverSolvePtr_ = &HYPRE_ParCSRLGMRESSolve;
      SolverPrecondPtr_ = &HYPRE_ParCSRLGMRESSetPrecond;
      break;
    case BiCGSTAB:
      if(IsSolverSetup_[0]){
        SolverDestroyPtr_(Solver_);
        IsSolverSetup_[0] = false;
      }
      SolverCreatePtr_ = &EpetraExt_HypreIJMatrix::Hypre_ParCSRBiCGSTABCreate;
      SolverDestroyPtr_ = &HYPRE_ParCSRBiCGSTABDestroy;
      SolverSetupPtr_ = &HYPRE_ParCSRBiCGSTABSetup;
      SolverSolvePtr_ = &HYPRE_ParCSRBiCGSTABSolve;
      SolverPrecondPtr_ = &HYPRE_ParCSRBiCGSTABSetPrecond;
      break;
    default:
      EPETRA_CHK_ERR(-1);
    }
  CreateSolver();
  } else {
  // Preconditioner
  switch(solver) {
    case BoomerAMG:
      if(IsPrecondSetup_[0]){
        PrecondDestroyPtr_(Preconditioner_);
        IsPrecondSetup_[0] = false;
      }
      PrecondCreatePtr_ = &EpetraExt_HypreIJMatrix::Hypre_BoomerAMGCreate;
      PrecondDestroyPtr_ = &HYPRE_BoomerAMGDestroy;
      PrecondSetupPtr_ = &HYPRE_BoomerAMGSetup;
      PrecondSolvePtr_ = &HYPRE_BoomerAMGSolve;
      break;
    case ParaSails:
      if(IsPrecondSetup_[0]){
        PrecondDestroyPtr_(Preconditioner_);
        IsPrecondSetup_[0] = false;
      }
      PrecondCreatePtr_ = &EpetraExt_HypreIJMatrix::Hypre_ParaSailsCreate;
      PrecondDestroyPtr_ = &HYPRE_ParaSailsDestroy;
      PrecondSetupPtr_ = &HYPRE_ParaSailsSetup;
      PrecondSolvePtr_ = &HYPRE_ParaSailsSolve;
      break;
    case Euclid:
      if(IsPrecondSetup_[0]){
        PrecondDestroyPtr_(Preconditioner_);
        IsPrecondSetup_[0] = false;
      }
      PrecondCreatePtr_ = &EpetraExt_HypreIJMatrix::Hypre_EuclidCreate;
      PrecondDestroyPtr_ = &HYPRE_EuclidDestroy;
      PrecondSetupPtr_ = &HYPRE_EuclidSetup;
      PrecondSolvePtr_ = &HYPRE_EuclidSolve;
      break;
    case AMS:
      if(IsPrecondSetup_[0]){
        PrecondDestroyPtr_(Preconditioner_);
        IsPrecondSetup_[0] = false;
      }
      PrecondCreatePtr_ = &EpetraExt_HypreIJMatrix::Hypre_AMSCreate;
      PrecondDestroyPtr_ = &HYPRE_AMSDestroy;
      PrecondSetupPtr_ = &HYPRE_AMSSetup;
      PrecondSolvePtr_ = &HYPRE_AMSSolve;
      break;
    default:
      EPETRA_CHK_ERR(-1);
    }
  CreatePrecond();

  }
  return 0;
} //SetParameter() - Choose solver or preconditioner type

//=======================================================
int EpetraExt_HypreIJMatrix::CreateSolver(){
  MPI_Comm comm;
  HYPRE_ParCSRMatrixGetComm(ParMatrix_, &comm);
  return (this->*SolverCreatePtr_)(comm, &Solver_);
} //CreateSolver()

//=======================================================
int EpetraExt_HypreIJMatrix::CreatePrecond(){
  MPI_Comm comm;
  HYPRE_ParCSRMatrixGetComm(ParMatrix_, &comm);
  return (this->*PrecondCreatePtr_)(comm, &Preconditioner_);
} //CreatePrecond()

//=======================================================
int EpetraExt_HypreIJMatrix::SetParameter(bool UsePreconditioner){
  if(UsePreconditioner == false){
    return 0;
  } else {
    if(SolverPrecondPtr_ == NULL){
      return -1;
    } else {
      SolverPrecondPtr_(Solver_, PrecondSolvePtr_, PrecondSetupPtr_, Preconditioner_);
      return 0;
    }
  }
} //SetParameter() - Set the precondioner pointer for solver

//=======================================================
int EpetraExt_HypreIJMatrix::SetParameter(Hypre_Chooser answer){
  SolveOrPrec_ = answer;
  return 0;
} //SetParameter() - Choose to solve or precondition

//=======================================================
int EpetraExt_HypreIJMatrix::SetupSolver() const{
  SolverSetupPtr_(Solver_, ParMatrix_, par_x, par_y);
  IsSolverSetup_[0] = true;
  return 0;
} //SetupSolver()

//=======================================================
int EpetraExt_HypreIJMatrix::SetupPrecond() const{
  PrecondSetupPtr_(Preconditioner_, ParMatrix_, par_x, par_y);
  IsPrecondSetup_[0] = true;
  return 0;
} //SetupPrecond()

//=======================================================
int EpetraExt_HypreIJMatrix::Solve(bool Upper, bool transpose, bool UnitDiagonal, const Epetra_MultiVector & X, Epetra_MultiVector & Y) const {
  bool SameVectors = false;
  int NumVectors = X.NumVectors();
  if (NumVectors != Y.NumVectors()) return -1;  // X and Y must have same number of vectors
  if(X.Pointers() == Y.Pointers()){
    SameVectors = true;
  }
  if(SolveOrPrec_ == Solver){
    if(IsSolverSetup_[0] == false){
      SetupSolver();
    }
  } else {
    if(IsPrecondSetup_[0] == false){
      SetupPrecond();
    }
  }
  for(int VecNum = 0; VecNum < NumVectors; VecNum++) {
    //Get values for current vector in multivector.
    double * x_values;
    EPETRA_CHK_ERR((*X(VecNum)).ExtractView(&x_values));
    double * y_values;
    if(!SameVectors){
      EPETRA_CHK_ERR((*Y(VecNum)).ExtractView(&y_values));
    } else {
      y_values = new double[X.MyLength()];
    }
    // Temporarily make a pointer to data in Hypre for end
    double *x_temp = x_local->data; 
    // Replace data in Hypre vectors with epetra values
    x_local->data = x_values;
    double *y_temp = y_local->data;
    y_local->data = y_values;
    
    EPETRA_CHK_ERR(HYPRE_ParVectorSetConstantValues(par_y, 0.0));
    if(transpose && !TransposeSolve_){
      // User requested a transpose solve, but the solver selected doesn't provide one
      EPETRA_CHK_ERR(-1);
    }
    if(SolveOrPrec_ == Solver){
      // Use the solver methods
      EPETRA_CHK_ERR(SolverSolvePtr_(Solver_, ParMatrix_, par_x, par_y));
    } else {
      // Apply the preconditioner
      EPETRA_CHK_ERR(PrecondSolvePtr_(Preconditioner_, ParMatrix_, par_x, par_y));
    }
    if(SameVectors){
      int NumEntries = Y.MyLength();
      std::vector<double> new_values; new_values.resize(NumEntries);
      std::vector<int> new_indices; new_indices.resize(NumEntries);
      for(int i = 0; i < NumEntries; i++){
        new_values[i] = y_values[i];
        new_indices[i] = i;
      }
      EPETRA_CHK_ERR((*Y(VecNum)).ReplaceMyValues(NumEntries, &new_values[0], &new_indices[0]));
      delete[] y_values;
    }
    x_local->data = x_temp;
    y_local->data = y_temp;
  }
  double flops = (double) NumVectors * (double) NumGlobalNonzeros();
  UpdateFlops(flops);
  return 0;
} //Solve()

//=======================================================
int EpetraExt_HypreIJMatrix::LeftScale(const Epetra_Vector& X) {
  for(int i = 0; i < NumMyRows_; i++){
    //Vector-scalar mult on ith row
    int num_entries;
    int *indices;
    double *values;
    EPETRA_CHK_ERR(HYPRE_ParCSRMatrixGetRow(ParMatrix_,i+MyRowStart_, &num_entries, &indices, &values));
    EPETRA_CHK_ERR(HYPRE_ParCSRMatrixRestoreRow(ParMatrix_,i+MyRowStart_, &num_entries, &indices, &values));
    Teuchos::Array<double> new_values; new_values.resize(num_entries);
    Teuchos::Array<int> new_indices; new_indices.resize(num_entries);
    for(int j = 0; j < num_entries; j++){
      // Scale this row with the appropriate values from the vector
      new_values[j] = X[i]*values[j];
      new_indices[j] = indices[j];
    }
    int rows[1];
    rows[0] = i+MyRowStart_;
    EPETRA_CHK_ERR(HYPRE_IJMatrixSetValues(Matrix_, 1, &num_entries, rows, &new_indices[0], &new_values[0]));
    // Finally set values of the Matrix for this row
  }
  HaveNumericConstants_ = false;
  UpdateFlops(NumGlobalNonzeros());
  return 0;
} //LeftScale()

//=======================================================
int EpetraExt_HypreIJMatrix::RightScale(const Epetra_Vector& X) {
  // First we need to import off-processor values of the vector
  Epetra_Import Importer(RowMatrixColMap(), RowMatrixRowMap());
  Epetra_Vector Import_Vector(RowMatrixColMap(), true);
  EPETRA_CHK_ERR(Import_Vector.Import(X, Importer, Insert, 0));
  
  for(int i = 0; i < NumMyRows_; i++){
    //Vector-scalar mult on ith column
    int num_entries;
    double *values;
    int *indices;
    // Get values and indices of ith row of matrix
    EPETRA_CHK_ERR(HYPRE_ParCSRMatrixGetRow(ParMatrix_,i+MyRowStart_, &num_entries, &indices, &values));
    EPETRA_CHK_ERR(HYPRE_ParCSRMatrixRestoreRow(ParMatrix_,i+MyRowStart_,&num_entries,&indices,&values));
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
    EPETRA_CHK_ERR(HYPRE_IJMatrixSetValues(Matrix_, 1, &num_entries, rows, &new_indices[0], &new_values[0]));
  }
  
  HaveNumericConstants_ = false;
  UpdateFlops(NumGlobalNonzeros());
  return 0;
} //RightScale()

//=======================================================
int EpetraExt_HypreIJMatrix::InitializeDefaults(){
  int my_type;
  // Get type of Hypre IJ Matrix
  EPETRA_CHK_ERR(HYPRE_IJMatrixGetObjectType(Matrix_, &my_type));
  MatType_ = my_type;
  // Currently only ParCSR is supported
  TEST_FOR_EXCEPTION(MatType_ != HYPRE_PARCSR, std::logic_error, "Object is not type PARCSR");

  // Get the actual ParCSR object from the IJ matrix
  EPETRA_CHK_ERR(HYPRE_IJMatrixGetObject(Matrix_, (void**) &ParMatrix_));
  
  int numRows, numCols;
  
  // Get dimensions of the matrix and store as global variables
  EPETRA_CHK_ERR(HYPRE_ParCSRMatrixGetDims(ParMatrix_, &numRows, &numCols));
  
  NumGlobalRows_ = numRows;
  NumGlobalCols_ = numCols;
  
  // Get the local dimensions of the matrix
  int ColStart, ColEnd;
  EPETRA_CHK_ERR(HYPRE_ParCSRMatrixGetLocalRange(ParMatrix_, &MyRowStart_, &MyRowEnd_, &ColStart, &ColEnd));
  
  // Determine number of local rows
  NumMyRows_ = MyRowEnd_ - MyRowStart_+1;
  
  return 0;
}  //InitializeDefaults() 

#endif /*HAVE_HYPRE*/
