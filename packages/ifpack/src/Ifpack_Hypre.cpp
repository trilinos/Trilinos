/*@HEADER
// ***********************************************************************
// 
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
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
*/
#include "Ifpack_Hypre.h"
#if defined(HAVE_HYPRE) && defined(HAVE_MPI)


#include "Ifpack_Utils.h"
//#include "io.h"
#include "Epetra_MpiComm.h"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#ifdef IFPACK_NODE_AWARE_CODE
extern int ML_NODE_ID;
#endif

using Teuchos::RCP;
using Teuchos::rcp;



Ifpack_Hypre::Ifpack_Hypre(Epetra_RowMatrix* A):
  A_(rcp(A,false)),
  UseTranspose_(false),
  IsParallel_(false),
  IsInitialized_(false),
  IsComputed_(false),
  Label_(),
  NumInitialize_(0),
  NumCompute_(0),
  NumApplyInverse_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyInverseTime_(0.0),
  ComputeFlops_(0.0),
  ApplyInverseFlops_(0.0),
  Time_(A_->Comm()),
  SolveOrPrec_(Solver),
  NumFunsToCall_(0),
  SolverType_(PCG),
  PrecondType_(Euclid),
  UsePreconditioner_(false),
  NiceRowMap_(true)
{
  IsSolverSetup_ = new bool[1];
  IsPrecondSetup_ = new bool[1];
  IsSolverSetup_[0] = false;
  IsPrecondSetup_[0] = false;
  if(Comm().NumProc() != 1) IsParallel_ = true;
  else IsParallel_ = false;
  MPI_Comm comm = GetMpiComm();
  // First make a copy of the Epetra Matrix as a Hypre Matrix
  int ilower = A_->RowMatrixRowMap().MinMyGID();
  int iupper = A_->RowMatrixRowMap().MaxMyGID();
  // Need to check if the RowMap is the way Hypre expects (if not more difficult)
  std::vector<int> ilowers; ilowers.resize(Comm().NumProc());
  std::vector<int> iuppers; iuppers.resize(Comm().NumProc());
  int myLower[1]; myLower[0] = ilower;
  int myUpper[1]; myUpper[0] = iupper;
  Comm().GatherAll(myLower, &ilowers[0], 1);
  Comm().GatherAll(myUpper, &iuppers[0], 1);
  for(int i = 0; i < Comm().NumProc()-1; i++){
    NiceRowMap_ = (NiceRowMap_ && iuppers[i]+1 == ilowers[i+1]);
  }
  if(!NiceRowMap_){
    ilower = (A_->NumGlobalRows() / Comm().NumProc())*Comm().MyPID();
    iupper = (A_->NumGlobalRows() / Comm().NumProc())*(Comm().MyPID()+1)-1;
    if(Comm().MyPID() == Comm().NumProc()-1){
      iupper = A_-> NumGlobalRows()-1;
    }
  }

  // Next create vectors that will be used when ApplyInverse() is called
  IFPACK_CHK_ERRV(HYPRE_IJVectorCreate(comm, ilower, iupper, &XHypre_));
  IFPACK_CHK_ERRV(HYPRE_IJVectorSetObjectType(XHypre_, HYPRE_PARCSR));
  IFPACK_CHK_ERRV(HYPRE_IJVectorInitialize(XHypre_));
  IFPACK_CHK_ERRV(HYPRE_IJVectorAssemble(XHypre_));
  IFPACK_CHK_ERRV(HYPRE_IJVectorGetObject(XHypre_, (void**) &ParX_));

  IFPACK_CHK_ERRV(HYPRE_IJVectorCreate(comm, ilower, iupper, &YHypre_));
  IFPACK_CHK_ERRV(HYPRE_IJVectorSetObjectType(YHypre_, HYPRE_PARCSR));
  IFPACK_CHK_ERRV(HYPRE_IJVectorInitialize(YHypre_));
  IFPACK_CHK_ERRV(HYPRE_IJVectorAssemble(YHypre_));
  IFPACK_CHK_ERRV(HYPRE_IJVectorGetObject(YHypre_, (void**) &ParY_));

  XVec_ = (hypre_ParVector *) hypre_IJVectorObject(((hypre_IJVector *) XHypre_));
  XLocal_ = hypre_ParVectorLocalVector(XVec_);

  YVec_ = (hypre_ParVector *) hypre_IJVectorObject(((hypre_IJVector *) YHypre_));
  YLocal_ = hypre_ParVectorLocalVector(YVec_);
  std::vector<int> rows; rows.resize(iupper - ilower +1);
  for(int i = ilower; i <= iupper; i++){
    rows[i-ilower] = i;
  }
  MySimpleMap_ = rcp(new Epetra_Map(-1, iupper-ilower+1, &rows[0], 0, Comm()));
  Aptr_ = dynamic_cast<Epetra_CrsMatrix*>(A_.get());
}

void Ifpack_Hypre::Destroy(){
  if(IsInitialized()){
    IFPACK_CHK_ERRV(HYPRE_IJMatrixDestroy(HypreA_));
  } 
  IFPACK_CHK_ERRV(HYPRE_IJVectorDestroy(XHypre_));
  IFPACK_CHK_ERRV(HYPRE_IJVectorDestroy(YHypre_));
  if(IsSolverSetup_[0]){
    IFPACK_CHK_ERRV(SolverDestroyPtr_(Solver_));
  }
  if(IsPrecondSetup_[0]){
    IFPACK_CHK_ERRV(PrecondDestroyPtr_(Preconditioner_));
  }
  delete[] IsSolverSetup_;
  delete[] IsPrecondSetup_;
}

int Ifpack_Hypre::Initialize(){
  Time_.ResetStartTime();
  MPI_Comm comm = GetMpiComm();
  int ilower = MySimpleMap_->MinMyGID();
  int iupper = MySimpleMap_->MaxMyGID();
  IFPACK_CHK_ERR(HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, &HypreA_));
  IFPACK_CHK_ERR(HYPRE_IJMatrixSetObjectType(HypreA_, HYPRE_PARCSR));
  IFPACK_CHK_ERR(HYPRE_IJMatrixInitialize(HypreA_));
  for(int i = 0; i < A_->NumMyRows(); i++){
    int numElements;
    IFPACK_CHK_ERR(A_->NumMyRowEntries(i,numElements));
    std::vector<int> indices; indices.resize(numElements);
    std::vector<double> values; values.resize(numElements);
    int numEntries;
    IFPACK_CHK_ERR(A_->ExtractMyRowCopy(i, numElements, numEntries, &values[0], &indices[0]));
    for(int j = 0; j < numEntries; j++){
      indices[j] = A_->RowMatrixColMap().GID(indices[j]);
    }
    int GlobalRow[1];
    GlobalRow[0] = A_->RowMatrixRowMap().GID(i);
    IFPACK_CHK_ERR(HYPRE_IJMatrixSetValues(HypreA_, 1, &numEntries, GlobalRow, &indices[0], &values[0]));
  }
  IFPACK_CHK_ERR(HYPRE_IJMatrixAssemble(HypreA_));
  IFPACK_CHK_ERR(HYPRE_IJMatrixGetObject(HypreA_, (void**)&ParMatrix_));
  IsInitialized_=true;
  NumInitialize_ = NumInitialize_ + 1;
  InitializeTime_ = InitializeTime_ + Time_.ElapsedTime();
  return 0;
}

int Ifpack_Hypre::SetParameters(Teuchos::ParameterList& list){
  List_ = list;
  Hypre_Solver solType = list.get("Solver", PCG);
  SolverType_ = solType;
  Hypre_Solver precType = list.get("Preconditioner", Euclid);
  PrecondType_ = precType;
  Hypre_Chooser chooser = list.get("SolveOrPrecondition", Solver);
  SolveOrPrec_ = chooser;
  bool SetPrecond = list.get("SetPreconditioner", false);
  IFPACK_CHK_ERR(SetParameter(SetPrecond));
  int NumFunctions = list.get("NumFunctions", 0);
  FunsToCall_.clear();
  NumFunsToCall_ = 0;
  if(NumFunctions > 0){
    RCP<FunctionParameter>* params = list.get<RCP<FunctionParameter>*>("Functions");
    for(int i = 0; i < NumFunctions; i++){
      IFPACK_CHK_ERR(AddFunToList(params[i]));
    }
  }
  return 0;
}

int Ifpack_Hypre::AddFunToList(RCP<FunctionParameter> NewFun){
  NumFunsToCall_ = NumFunsToCall_+1;
  FunsToCall_.resize(NumFunsToCall_);
  FunsToCall_[NumFunsToCall_-1] = NewFun;
  return 0;
}

int Ifpack_Hypre::SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, int), int parameter){
  RCP<FunctionParameter> temp = rcp(new FunctionParameter(chooser, pt2Func, parameter));
  IFPACK_CHK_ERR(AddFunToList(temp));
  return 0;
}

int Ifpack_Hypre::SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, double), double parameter){
  RCP<FunctionParameter> temp = rcp(new FunctionParameter(chooser, pt2Func, parameter));
  IFPACK_CHK_ERR(AddFunToList(temp));
  return 0;
}

int Ifpack_Hypre::SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, double, int), double parameter1, int parameter2){
  RCP<FunctionParameter> temp = rcp(new FunctionParameter(chooser, pt2Func, parameter1, parameter2));
  IFPACK_CHK_ERR(AddFunToList(temp));
  return 0;
}

int Ifpack_Hypre::SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, int, int), int parameter1, int parameter2){
  RCP<FunctionParameter> temp = rcp(new FunctionParameter(chooser, pt2Func, parameter1, parameter2));
  IFPACK_CHK_ERR(AddFunToList(temp));
  return 0;
}

int Ifpack_Hypre::SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, double*), double* parameter){
  RCP<FunctionParameter> temp = rcp(new FunctionParameter(chooser, pt2Func, parameter));
  IFPACK_CHK_ERR(AddFunToList(temp));
  return 0;
}

int Ifpack_Hypre::SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, int*), int* parameter){
  RCP<FunctionParameter> temp = rcp(new FunctionParameter(chooser, pt2Func, parameter));
  IFPACK_CHK_ERR(AddFunToList(temp));
  return 0;
}

int Ifpack_Hypre::SetParameter(Hypre_Chooser chooser, Hypre_Solver solver){
  if(chooser == Solver){
    SolverType_ = solver;
  } else {
    PrecondType_ = solver;
  }
  return 0;
}


int Ifpack_Hypre::SetParameter(Hypre_Chooser answer){
  SolveOrPrec_ = answer;
  return 0;
}

int Ifpack_Hypre::SetParameter(bool UsePreconditioner){
  UsePreconditioner_ = UsePreconditioner;
  return 0;
}

int Ifpack_Hypre::Compute(){
  if(IsInitialized() == false){
    IFPACK_CHK_ERR(Initialize());
  }
  Time_.ResetStartTime();
  IFPACK_CHK_ERR(SetSolverType(SolverType_));
  IFPACK_CHK_ERR(SetPrecondType(PrecondType_));
  for(int i = 0; i < NumFunsToCall_; i++){
    IFPACK_CHK_ERR(FunsToCall_[i]->CallFunction(Solver_, Preconditioner_));
  }
  if(UsePreconditioner_){
    if(SolverPrecondPtr_ != NULL){
      IFPACK_CHK_ERR(SolverPrecondPtr_(Solver_, PrecondSolvePtr_, PrecondSetupPtr_, Preconditioner_));
    }
  }
  if(SolveOrPrec_ == Solver){
    IFPACK_CHK_ERR(SolverSetupPtr_(Solver_, ParMatrix_, ParX_, ParY_));
    IsSolverSetup_[0] = true;
  } else {
    IFPACK_CHK_ERR(PrecondSetupPtr_(Preconditioner_, ParMatrix_, ParX_, ParY_));
    IsPrecondSetup_[0] = true;
  }
  IsComputed_ = true;
  NumCompute_ = NumCompute_ + 1;
  ComputeTime_ = ComputeTime_ + Time_.ElapsedTime();
  return 0;
}

int Ifpack_Hypre::CallFunctions() const{
  for(int i = 0; i < NumFunsToCall_; i++){
    IFPACK_CHK_ERR(FunsToCall_[i]->CallFunction(Solver_, Preconditioner_));
  }
  
  return 0;
}

const Epetra_Map& Ifpack_Hypre::OperatorDomainMap() const{
  return *MySimpleMap_;
}

const Epetra_Map& Ifpack_Hypre::OperatorRangeMap() const{
  return *MySimpleMap_;
}


int Ifpack_Hypre::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const{
  if(IsComputed() == false){
    IFPACK_CHK_ERR(-1);
  }
  Time_.ResetStartTime();
  bool SameVectors = false;
  int NumVectors = X.NumVectors();
  if (NumVectors != Y.NumVectors()) IFPACK_CHK_ERR(-1);  // X and Y must have same number of vectors
  if(X.Pointers() == Y.Pointers()){
    SameVectors = true;
  }
  for(int VecNum = 0; VecNum < NumVectors; VecNum++) {
    //Get values for current vector in multivector.
    double * XValues;
    IFPACK_CHK_ERR((*X(VecNum)).ExtractView(&XValues));
    double * YValues;
    if(!SameVectors){
      IFPACK_CHK_ERR((*Y(VecNum)).ExtractView(&YValues));
    } else {
      YValues = new double[X.MyLength()];
    }
    // Temporarily make a pointer to data in Hypre for end
    double *XTemp = XLocal_->data;
    // Replace data in Hypre vectors with epetra values
    XLocal_->data = XValues;
    double *YTemp = YLocal_->data;
    YLocal_->data = YValues;

    IFPACK_CHK_ERR(HYPRE_ParVectorSetConstantValues(ParY_, 0.0));
    if(SolveOrPrec_ == Solver){
      // Use the solver methods
      IFPACK_CHK_ERR(SolverSolvePtr_(Solver_, ParMatrix_, ParX_, ParY_));
    } else {
      // Apply the preconditioner
      IFPACK_CHK_ERR(PrecondSolvePtr_(Preconditioner_, ParMatrix_, ParX_, ParY_));
    }
    if(SameVectors){
      int NumEntries = Y.MyLength();
      std::vector<double> new_values; new_values.resize(NumEntries);
      std::vector<int> new_indices; new_indices.resize(NumEntries);
      for(int i = 0; i < NumEntries; i++){
        new_values[i] = YValues[i];
        new_indices[i] = i;
      }
      IFPACK_CHK_ERR((*Y(VecNum)).ReplaceMyValues(NumEntries, &new_values[0], &new_indices[0]));
      delete[] YValues;
    }
    XLocal_->data = XTemp;
    YLocal_->data = YTemp;
  }
  NumApplyInverse_ = NumApplyInverse_ + 1;
  ApplyInverseTime_ = ApplyInverseTime_ + Time_.ElapsedTime();
  return 0;
}

int Ifpack_Hypre::Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const{
  if(IsInitialized() == false){
    IFPACK_CHK_ERR(-1);
  }
  bool SameVectors = false;
  int NumVectors = X.NumVectors();
  if (NumVectors != Y.NumVectors()) IFPACK_CHK_ERR(-1);  // X and Y must have same number of vectors
  if(X.Pointers() == Y.Pointers()){
    SameVectors = true;
  }
  for(int VecNum = 0; VecNum < NumVectors; VecNum++) {
    //Get values for current vector in multivector.
    double * XValues;
    double * YValues;
    IFPACK_CHK_ERR((*X(VecNum)).ExtractView(&XValues));
    double *XTemp = XLocal_->data;
    double *YTemp = YLocal_->data;
    if(!SameVectors){
      IFPACK_CHK_ERR((*Y(VecNum)).ExtractView(&YValues));
    } else {
      YValues = new double[X.MyLength()];
    }
    YLocal_->data = YValues;
    IFPACK_CHK_ERR(HYPRE_ParVectorSetConstantValues(ParY_,0.0));
    // Temporarily make a pointer to data in Hypre for end
    // Replace data in Hypre vectors with epetra values
    XLocal_->data = XValues;
    // Do actual computation.
    if(TransA) {
      // Use transpose of A in multiply
      IFPACK_CHK_ERR(HYPRE_ParCSRMatrixMatvecT(1.0, ParMatrix_, ParX_, 1.0, ParY_));
    } else {
      IFPACK_CHK_ERR(HYPRE_ParCSRMatrixMatvec(1.0, ParMatrix_, ParX_, 1.0, ParY_));
    }
    if(SameVectors){
      int NumEntries = Y.MyLength();
      std::vector<double> new_values; new_values.resize(NumEntries);
      std::vector<int> new_indices; new_indices.resize(NumEntries);
      for(int i = 0; i < NumEntries; i++){
        new_values[i] = YValues[i];
        new_indices[i] = i;
      }
      IFPACK_CHK_ERR((*Y(VecNum)).ReplaceMyValues(NumEntries, &new_values[0], &new_indices[0]));
      delete[] YValues;
    }
    XLocal_->data = XTemp;
    YLocal_->data = YTemp;
  }
  return 0;
}

ostream& Ifpack_Hypre::Print(ostream& os) const{
  if (!Comm().MyPID()) {
    os << endl;
    os << "================================================================================" << endl;
    os << "Ifpack_Hypre: " << Label () << endl << endl;
    os << "Using " << Comm().NumProc() << " processors." << endl;
    os << "Global number of rows            = " << A_->NumGlobalRows() << endl;
    os << "Global number of nonzeros        = " << A_->NumGlobalNonzeros() << endl;
    os << "Condition number estimate = " << Condest() << endl;
    os << endl;
    os << "Phase           # calls   Total Time (s)       Total MFlops     MFlops/s" << endl;
    os << "-----           -------   --------------       ------------     --------" << endl;
    os << "Initialize()    "   << std::setw(5) << NumInitialize_
       << "  " << std::setw(15) << InitializeTime_
       << "              0.0              0.0" << endl;
    os << "Compute()       "   << std::setw(5) << NumCompute_
       << "  " << std::setw(15) << ComputeTime_
       << "  " << std::setw(15) << 1.0e-6 * ComputeFlops_;
    if (ComputeTime_ != 0.0)
      os << "  " << std::setw(15) << 1.0e-6 * ComputeFlops_ / ComputeTime_ << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "ApplyInverse()  "   << std::setw(5) << NumApplyInverse_
       << "  " << std::setw(15) << ApplyInverseTime_
       << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops_;
    if (ApplyInverseTime_ != 0.0)
      os << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops_ / ApplyInverseTime_ << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "================================================================================" << endl;
    os << endl;
  }
  return os;
}


double Ifpack_Hypre::Condest(const Ifpack_CondestType CT, 
                             const int MaxIters,
                             const double Tol,
                             Epetra_RowMatrix* Matrix_in){
  if (!IsComputed()) // cannot compute right now
    return(-1.0);

  Condest_ = Ifpack_Condest(*this, CT, MaxIters, Tol, Matrix_in);

  return(Condest_);
}

int Ifpack_Hypre::SetSolverType(Hypre_Solver Solver){
  switch(Solver) {
    case BoomerAMG:
      if(IsSolverSetup_[0]){
        SolverDestroyPtr_(Solver_);
        IsSolverSetup_[0] = false;
      }
      SolverCreatePtr_ = &Ifpack_Hypre::Hypre_BoomerAMGCreate;
      SolverDestroyPtr_ = &HYPRE_BoomerAMGDestroy;
      SolverSetupPtr_ = &HYPRE_BoomerAMGSetup;
      SolverPrecondPtr_ = NULL;
      SolverSolvePtr_ = &HYPRE_BoomerAMGSolve;
      break;
    case AMS:
      if(IsSolverSetup_[0]){
        SolverDestroyPtr_(Solver_);
        IsSolverSetup_[0] = false;
      }
      SolverCreatePtr_ = &Ifpack_Hypre::Hypre_AMSCreate;
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
      SolverCreatePtr_ = &Ifpack_Hypre::Hypre_ParCSRHybridCreate;
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
      SolverCreatePtr_ = &Ifpack_Hypre::Hypre_ParCSRPCGCreate;
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
      SolverCreatePtr_ = &Ifpack_Hypre::Hypre_ParCSRGMRESCreate;
      SolverDestroyPtr_ = &HYPRE_ParCSRGMRESDestroy;
      SolverSetupPtr_ = &HYPRE_ParCSRGMRESSetup;
    
      SolverPrecondPtr_ = &HYPRE_ParCSRGMRESSetPrecond;
      break;
    case FlexGMRES:
      if(IsSolverSetup_[0]){
        SolverDestroyPtr_(Solver_);
        IsSolverSetup_[0] = false;
      }
      SolverCreatePtr_ = &Ifpack_Hypre::Hypre_ParCSRFlexGMRESCreate;
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
      SolverCreatePtr_ = &Ifpack_Hypre::Hypre_ParCSRLGMRESCreate;
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
      SolverCreatePtr_ = &Ifpack_Hypre::Hypre_ParCSRBiCGSTABCreate;
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

int Ifpack_Hypre::SetPrecondType(Hypre_Solver Precond){
  switch(Precond) {
    case BoomerAMG:
      if(IsPrecondSetup_[0]){
        PrecondDestroyPtr_(Preconditioner_);
        IsPrecondSetup_[0] = false;
      }
      PrecondCreatePtr_ = &Ifpack_Hypre::Hypre_BoomerAMGCreate;
      PrecondDestroyPtr_ = &HYPRE_BoomerAMGDestroy;
      PrecondSetupPtr_ = &HYPRE_BoomerAMGSetup;
      PrecondSolvePtr_ = &HYPRE_BoomerAMGSolve;
      break;
    case ParaSails:
      if(IsPrecondSetup_[0]){
        PrecondDestroyPtr_(Preconditioner_);
        IsPrecondSetup_[0] = false;
      }
      PrecondCreatePtr_ = &Ifpack_Hypre::Hypre_ParaSailsCreate;
      PrecondDestroyPtr_ = &HYPRE_ParaSailsDestroy;
      PrecondSetupPtr_ = &HYPRE_ParaSailsSetup;
      PrecondSolvePtr_ = &HYPRE_ParaSailsSolve;
      break;
    case Euclid:
      if(IsPrecondSetup_[0]){
        PrecondDestroyPtr_(Preconditioner_);
        IsPrecondSetup_[0] = false;
      }
      PrecondCreatePtr_ = &Ifpack_Hypre::Hypre_EuclidCreate;
      PrecondDestroyPtr_ = &HYPRE_EuclidDestroy;
      PrecondSetupPtr_ = &HYPRE_EuclidSetup;
      PrecondSolvePtr_ = &HYPRE_EuclidSolve;
      break;
    case AMS:
      if(IsPrecondSetup_[0]){
        PrecondDestroyPtr_(Preconditioner_);
        IsPrecondSetup_[0] = false;
      }
      PrecondCreatePtr_ = &Ifpack_Hypre::Hypre_AMSCreate;
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

int Ifpack_Hypre::CreateSolver(){
  MPI_Comm comm;
  HYPRE_ParCSRMatrixGetComm(ParMatrix_, &comm);
  return (this->*SolverCreatePtr_)(comm, &Solver_);
}

int Ifpack_Hypre::CreatePrecond(){
  MPI_Comm comm;
  HYPRE_ParCSRMatrixGetComm(ParMatrix_, &comm);
  return (this->*PrecondCreatePtr_)(comm, &Preconditioner_);
}

int Ifpack_Hypre::Hypre_BoomerAMGCreate(MPI_Comm comm, HYPRE_Solver *solver){
  return HYPRE_BoomerAMGCreate(solver);
}

int Ifpack_Hypre::Hypre_ParaSailsCreate(MPI_Comm comm, HYPRE_Solver *solver){
  return HYPRE_ParaSailsCreate(comm, solver);
}

int Ifpack_Hypre::Hypre_EuclidCreate(MPI_Comm comm, HYPRE_Solver *solver){
  return HYPRE_EuclidCreate(comm, solver);
}

int Ifpack_Hypre::Hypre_AMSCreate(MPI_Comm comm, HYPRE_Solver *solver){
  return HYPRE_AMSCreate(solver);
}

int Ifpack_Hypre::Hypre_ParCSRHybridCreate(MPI_Comm comm, HYPRE_Solver *solver){
  return HYPRE_ParCSRHybridCreate(solver);
}

int Ifpack_Hypre::Hypre_ParCSRPCGCreate(MPI_Comm comm, HYPRE_Solver *solver){
  return HYPRE_ParCSRPCGCreate(comm, solver);
}

int Ifpack_Hypre::Hypre_ParCSRGMRESCreate(MPI_Comm comm, HYPRE_Solver *solver){
  return HYPRE_ParCSRGMRESCreate(comm, solver);
}

int Ifpack_Hypre::Hypre_ParCSRFlexGMRESCreate(MPI_Comm comm, HYPRE_Solver *solver){
  return HYPRE_ParCSRFlexGMRESCreate(comm, solver);
}

int Ifpack_Hypre::Hypre_ParCSRLGMRESCreate(MPI_Comm comm, HYPRE_Solver *solver){
  return HYPRE_ParCSRLGMRESCreate(comm, solver);
}

int Ifpack_Hypre::Hypre_ParCSRBiCGSTABCreate(MPI_Comm comm, HYPRE_Solver *solver){
  return HYPRE_ParCSRBiCGSTABCreate(comm, solver);
}

#endif // HAVE_HYPRE && HAVE_MPI
