/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#include "ml_common.h"
#ifdef HAVE_ML_MLAPI
#include "ml_RowMatrix.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_PointRelaxation.h"
#include "Ifpack_IC.h"
#include "Ifpack_ICT.h"
#include "Ifpack_ILU.h"
#include "Ifpack_ILUT.h"
#include "Ifpack_AdditiveSchwarz.h"
#include "Ifpack_Amesos.h"
#include "MLAPI_Error.h"
#include "MLAPI_Space.h"
#include "MLAPI_BaseOperator.h"
#include "MLAPI_CompObject.h"
#include "MLAPI_TimeObject.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_MultiVector_Utils.h"
#include "MLAPI_Workspace.h"
#include "MLAPI_LoadBalanceInverseOperator.h"

namespace MLAPI
{


LoadBalanceInverseOperator::LoadBalanceInverseOperator(const LoadBalanceInverseOperator& RHS)
: InverseOperator(RHS)
{
  Op_           = RHS.GetOperator();
  RCPRowMatrix_ = RHS.RCPRowMatrix();
  RCPData_      = RHS.GetRCPData();
  SetLabel(RHS.GetLabel());
}

LoadBalanceInverseOperator& LoadBalanceInverseOperator::operator=(const LoadBalanceInverseOperator& RHS)
{
  if (this == &RHS)
    return(*this);

  Op_           = RHS.GetOperator();
  RCPRowMatrix_ = RHS.RCPRowMatrix();
  RCPData_      = RHS.GetRCPData();

  SetLabel(RHS.GetLabel());
  return(*this);
}

void LoadBalanceInverseOperator::Reshape()
{
  Destroy();
}


void LoadBalanceInverseOperator::Reshape(Ifpack_Preconditioner* prec,
                                         const LoadBalanceOperator& Op,
                                         const bool ownership)
{
  ResetTimer();
  StackPush();

  Op_ = Op;

  if (GetParticipation()) RCPRowMatrix_ = Op.GetRCPRowMatrix();
  else                    RCPRowMatrix_ = Teuchos::null;

  if (GetParticipation()) RCPData_ = Teuchos::rcp(prec,ownership);
  else                    RCPData_ = Teuchos::null;

  StackPop();
  UpdateTime();
}


const Space LoadBalanceInverseOperator::GetOperatorRangeSpace() const {
  return(Op_.GetRangeSpace());
}

const Space LoadBalanceInverseOperator::GetOperatorDomainSpace() const {
  return(Op_.GetDomainSpace());
}

const Space LoadBalanceInverseOperator::GetRangeSpace() const {
  return(Op_.GetRangeSpace());
}

const Space LoadBalanceInverseOperator::GetDomainSpace() const {
  return(Op_.GetDomainSpace());
}

const Teuchos::RCP<Epetra_RowMatrix> LoadBalanceInverseOperator::RCPRowMatrix() const
{
  return(RCPRowMatrix_);
}

Epetra_RowMatrix* LoadBalanceInverseOperator::RowMatrix() const
{
  return(RCPRowMatrix_.get());
}

const LoadBalanceOperator& LoadBalanceInverseOperator::GetOperator() const
{
  return(Op_);
}

Teuchos::RCP<Ifpack_Preconditioner>& LoadBalanceInverseOperator::GetRCPData()
{
  return(RCPData_);
}

const Teuchos::RCP<Ifpack_Preconditioner>& LoadBalanceInverseOperator::GetRCPData() const
{
  return(RCPData_);
}

int LoadBalanceInverseOperator::Apply(const MultiVector& x, MultiVector& y) const
{
  ResetTimer();
  StackPush();

  if (GetDomainSpace() != x.GetVectorSpace())
    ML_THROW("DomainSpace and x.GetVectorSpace() differ", -1);

  if (GetRangeSpace() != y.GetVectorSpace())
    ML_THROW("RangeSpace and y.GetVectorSpace() differ", -1);

  int x_nv = x.GetNumVectors();
  int y_nv = y.GetNumVectors();
  double FL = 0.0;
  if (RCPData_ != Teuchos::null)
    FL = RCPData_->ComputeFlops();

  if (x_nv != y_nv)
    ML_THROW("Number of vectors of x and y differ (" +
             GetString(x_nv) + " vs. " + GetString(x_nv), -1);

  for (int v = 0 ; v < x_nv ; ++v) {
    if (GetParticipation()) // some procs have no internal RCPData_
    {
      Epetra_Vector x_Epetra(View,RowMatrix()->OperatorDomainMap(),
                             (double*)&(x(0,v)));
      Epetra_Vector y_Epetra(View,RowMatrix()->OperatorRangeMap(),
                             (double*)&(y(0,v)));

      if (RCPData_ != Teuchos::null)
        RCPData_->ApplyInverse(x_Epetra,y_Epetra);
      else
        ML_THROW("Ifpack smoother is not properly set up", -1);
    }
  }

  StackPop();
  if (RCPData_ != Teuchos::null)
    UpdateFlops(RCPData_->ComputeFlops() - FL);
  UpdateTime();

  return(0);
}

MultiVector LoadBalanceInverseOperator::operator()(const MultiVector& LHS)
{
  StackPush();

  MultiVector RHS(LHS.GetVectorSpace());
  RHS = 0.0;
  Apply(LHS,RHS);

  StackPop();

  return(RHS);
}

MultiVector LoadBalanceInverseOperator::operator()(const MultiVector& LHS,
                       const MultiVector& RHS)
{
  MultiVector RHS2 = Duplicate(RHS);
  Apply(LHS,RHS2);
  return(RHS2);
}

std::ostream& LoadBalanceInverseOperator::Print(std::ostream& os, const bool /* verbose */) const
{

  StackPush();

  if (GetMyPID() == 0) {
    os << "***MLAPI::InverseOperator" << std::endl;
    os << "Label             = " << GetLabel() << std::endl;
    os << "Number of rows    = " << GetRangeSpace().GetNumGlobalElements() << std::endl;
    os << "Number of columns = " << GetRangeSpace().GetNumGlobalElements() << std::endl;
    os << "Flop count        = " << GetFlops() << std::endl;
    os << "Cumulative time   = " << GetTime() << std::endl;
    if (GetTime() != 0.0)
      os << "MFlops rate       = " << 1.0e-6 * GetFlops() / GetTime() << std::endl;
    else
      os << "MFlops rate       = 0.0" << std::endl;
    os << std::endl;
  }

  StackPop();

  return(os);

}

void LoadBalanceInverseOperator::Destroy()
{
  Op_.Reshape();
  RCPRowMatrix_ = Teuchos::null;
  RCPData_      = Teuchos::null;
}

} // namespace MLAPI
#endif
