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
#include "MLAPI_InverseOperator.h"

namespace MLAPI
{

InverseOperator::InverseOperator(const Operator& Op, const string Type)
{
  Reshape(Op, Type);
}

InverseOperator::InverseOperator(const Operator& Op, const string Type,
                Teuchos::ParameterList& List)
{
  Reshape(Op, Type, List);
}

InverseOperator::InverseOperator(const InverseOperator& RHS)
{
  Op_           = RHS.GetOperator();
  RCPRowMatrix_ = RHS.RCPRowMatrix();
  RCPData_      = RHS.GetRCPData();
  RCPMLPrec_    = RHS.GetRCPMLPrec();
  SetLabel(RHS.GetLabel());
}

InverseOperator& InverseOperator::operator=(const InverseOperator& RHS)
{
  if (this == &RHS)
    return(*this);

  Op_           = RHS.GetOperator();
  RCPRowMatrix_ = RHS.RCPRowMatrix();
  RCPData_      = RHS.GetRCPData();
  RCPMLPrec_    = RHS.GetRCPMLPrec();

  SetLabel(RHS.GetLabel());
  return(*this);
}

void InverseOperator::Reshape()
{
  Destroy();
}

void InverseOperator::Reshape(const Operator& Op, const string Type)
{
  Teuchos::ParameterList List;
  Reshape(Op, Type, List);
}

void InverseOperator::Reshape(Ifpack_Preconditioner* prec, const Operator& Op, 
                              const bool ownership)
{
  ResetTimer();
  StackPush();

  Op_ = Op;

  RCPRowMatrix_ = Op.GetRCPRowMatrix();

  RCPData_ = Teuchos::rcp(prec,ownership);

  StackPop();
  UpdateTime();
}

void InverseOperator::Reshape(const Operator& Op, const string Type,
             Teuchos::ParameterList& List, Teuchos::ParameterList* pushlist)
{
  ResetTimer();
  StackPush();

  Op_ = Op;

  RCPRowMatrix_ = Op.GetRCPRowMatrix();

  // FIXME: to add overlap and level-of-fill
  int NumSweeps   = List.get("smoother: sweeps", 1);
  double Damping  = List.get("smoother: damping factor", 0.67); 
  int LOF_ilu     = List.get("smoother: ilu fill", 0);
  double LOF_ict  = List.get("smoother: ilut fill", 1.0);
  double LOF_ilut = List.get("smoother: ict fill", 1.0);
  string reorder  = List.get("schwarz: reordering type","rcm");

  Teuchos::ParameterList IFPACKList;

  // any parameters from the main list List are overwritten by the pushlist
  IFPACKList.set("relaxation: sweeps", NumSweeps);
  IFPACKList.set("relaxation: damping factor", Damping);
  IFPACKList.set("fact: level-of-fill", LOF_ilu);
  IFPACKList.set("fact: ict level-of-fill", LOF_ict);
  IFPACKList.set("fact: ilut level-of-fill", LOF_ilut);
  IFPACKList.set("relaxation: zero starting solution", false);
  IFPACKList.set("schwarz: reordering type",reorder);
  IFPACKList.set("fact: relative threshold",1.0);

  // if present, the pushlist is assumed to be a preconstructed ifpack list
  // that is copied straight to the ifpack list here
  // entries in pushlist overwrite previous list entries
  if (pushlist)
    IFPACKList.setParameters(*pushlist);


  bool verbose = false; //(GetMyPID() == 0 && GetPrintLevel() > 5);

  // the ML smoother
  RCPMLPrec_ = Teuchos::null;
  // The Ifpack smoother
  Ifpack_Preconditioner* Prec = NULL;

  if (Type == "Jacobi") {
    if (verbose) {
      cout << "Damping factor = " << Damping 
        << ", sweeps = " << NumSweeps << endl;
      cout << endl;
    }
    IFPACKList.set("relaxation: type", "Jacobi");
    Prec = new Ifpack_PointRelaxation(RowMatrix());
  }
  else if (Type == "Gauss-Seidel") {
    if (verbose) {
      cout << "Damping factor = " << Damping 
        << ", sweeps = " << NumSweeps << endl;
      cout << endl;
    }
    IFPACKList.set("relaxation: type", "Gauss-Seidel");
    Prec = new Ifpack_PointRelaxation(RowMatrix());
  }
  else if (Type == "symmetric Gauss-Seidel") {
    if (verbose) {
      cout << "Damping factor = " << Damping 
        << ", sweeps = " << NumSweeps << endl;
      cout << endl;
    }
    IFPACKList.set("relaxation: type", "symmetric Gauss-Seidel");

    Prec = new Ifpack_PointRelaxation(RowMatrix());
  }
  else if (Type == "ILU") {
    if (verbose) {
      cout << "ILU factorization, ov = 0, no reordering, LOF = "
        << LOF_ilu << endl;
      cout << endl;
    }
    
    // use the Additive Schwarz class because it does reordering
    Prec = new Ifpack_AdditiveSchwarz<Ifpack_ILU>(RowMatrix());
  }
  else if (Type == "ILUT") {
    if (verbose) {
      cout << "ILUT factorization, ov = 0, no reordering, LOF = "
        << LOF_ilu << endl;
      cout << endl;
    }
    Prec = new Ifpack_ILUT(RowMatrix());
  }
  else if (Type == "IC") {
    if (verbose) {
      cout << "IC factorization, ov = 0, no reordering, LOF = "
        << LOF_ilu << endl;
      cout << endl;
    }
    Prec = new Ifpack_IC(RowMatrix());
  }
  else if (Type == "ICT") {
    if (verbose) {
      cout << "ICT factorization, ov = 0, no reordering, LOF = "
        << LOF_ilu << endl;
      cout << endl;
    }
    Prec = new Ifpack_ICT(RowMatrix());
  }
  else if (Type == "LU") {
    if (verbose) {
      cout << "LU factorization, ov = 0, local solver = KLU" << endl;
      cout << endl;
    }
    Prec = new Ifpack_AdditiveSchwarz<Ifpack_Amesos>(RowMatrix());
  }
  else if (Type == "Amesos" || Type == "Amesos-KLU")  {
    if (verbose) {
      cout << "Amesos-KLU direct solver" << endl;
      cout << endl;
    }
    Prec = new Ifpack_Amesos(RowMatrix());
  }
  else if (Type == "MLS" || Type == "ML MLS" || 
           Type == "ML symmetric Gauss-Seidel" ||
           Type == "ML Gauss-Seidel") 
  {
    if (verbose) {
      cout << "ML's MLS smoother" << endl;
      cout << endl;
      }
    Teuchos::ParameterList mlparams;
    ML_Epetra::SetDefaults("SA",mlparams);
    int output = List.get("ML output",-47);
    if (output == -47) output = List.get("output",-1);
    if (output != -1) mlparams.set("ML output",output);
    mlparams.set("max levels",1);
    int sweeps = List.get("smoother: sweeps",1);
    mlparams.set("coarse: sweeps",sweeps);
    double damp = List.get("smoother: damping factor",0.67);
    mlparams.set("coarse: damping factor",damp);
    mlparams.set("zero starting solution", false);
    if (Type == "MLS" || Type == "ML MLS")
    {
      mlparams.set("coarse: type","MLS"); // MLS symmetric Gauss-Seidel Amesos-KLU
      int poly = List.get("smoother: MLS polynomial order",3);
      mlparams.set("coarse: MLS polynomial order",poly);
    }
    else if (Type == "ML symmetric Gauss-Seidel")
      mlparams.set("coarse: type","symmetric Gauss-Seidel"); // MLS symmetric Gauss-Seidel Amesos-KLU
    else if (Type == "ML Gauss-Seidel")
      mlparams.set("coarse: type","Gauss-Seidel");
    else if (Type == "ML Jacobi")
      mlparams.set("coarse: type","Jacobi");
    else
      ML_THROW("Requested type (" + Type + ") not recognized", -1);
    RCPMLPrec_ = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*RowMatrix(),mlparams,true));
  }
  else
    ML_THROW("Requested type (" + Type + ") not recognized", -1);

  if (Prec)
  {
    RCPData_ = Teuchos::rcp(Prec);
    RCPData_->SetParameters(IFPACKList);
    RCPData_->Initialize();
    RCPData_->Compute();
    UpdateFlops(RCPData_->InitializeFlops());
    UpdateFlops(RCPData_->ComputeFlops());
  }
  else 
    RCPData_ = Teuchos::null;
    
  StackPop();
  UpdateTime();

}

const Space InverseOperator::GetOperatorRangeSpace() const {
  return(Op_.GetRangeSpace());
}

const Space InverseOperator::GetOperatorDomainSpace() const {
  return(Op_.GetDomainSpace());
}

const Space InverseOperator::GetRangeSpace() const {
  return(Op_.GetRangeSpace());
}

const Space InverseOperator::GetDomainSpace() const {
  return(Op_.GetDomainSpace());
}

const Teuchos::RefCountPtr<Epetra_RowMatrix> InverseOperator::RCPRowMatrix() const
{
  return(RCPRowMatrix_);
}

Epetra_RowMatrix* InverseOperator::RowMatrix() const
{
  return(RCPRowMatrix_.get());
}

const Operator& InverseOperator::GetOperator() const
{
  return(Op_);
}

Teuchos::RefCountPtr<Ifpack_Preconditioner>& InverseOperator::GetRCPData() 
{
  return(RCPData_);
}

const Teuchos::RefCountPtr<Ifpack_Preconditioner>& 
InverseOperator::GetRCPData() const
{
  return(RCPData_);
}

Teuchos::RefCountPtr<ML_Epetra::MultiLevelPreconditioner>& InverseOperator::GetRCPMLPrec() 
{
  return(RCPMLPrec_);
}

const Teuchos::RefCountPtr<ML_Epetra::MultiLevelPreconditioner>& 
InverseOperator::GetRCPMLPrec() const
{
  return(RCPMLPrec_);
}

int 
InverseOperator::Apply(const MultiVector& x, MultiVector& y) const
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

    Epetra_Vector x_Epetra(View,RowMatrix()->OperatorDomainMap(),
                           (double*)&(x(0,v)));
    Epetra_Vector y_Epetra(View,RowMatrix()->OperatorRangeMap(),
                           (double*)&(y(0,v)));

    if (RCPData_ != Teuchos::null)
      RCPData_->ApplyInverse(x_Epetra,y_Epetra);
    else if (RCPMLPrec_ != Teuchos::null)
      RCPMLPrec_->ApplyInverse(x_Epetra,y_Epetra);
    else
      ML_THROW("Neither Ifpack nor ML smoother is properly set up", -1);
  }

  StackPop();
  if (RCPData_ != Teuchos::null)
    UpdateFlops(RCPData_->ComputeFlops() - FL);
  UpdateTime();

  return(0);
}

MultiVector 
InverseOperator::operator()(const MultiVector& LHS)
{
  StackPush();

  MultiVector RHS(LHS.GetVectorSpace());
  RHS = 0.0;
  Apply(LHS,RHS);

  StackPop();

  return(RHS);
}

MultiVector 
InverseOperator::operator()(const MultiVector& LHS,
                       const MultiVector& RHS)
{
  MultiVector RHS2 = Duplicate(RHS);
  Apply(LHS,RHS2);
  return(RHS2);
}

ostream& 
InverseOperator::Print(std::ostream& os, const bool verbose) const
{

  StackPush();

  if (GetMyPID() == 0) {
    os << "***MLAPI::InverseOperator" << endl;
    os << "Label             = " << GetLabel() << endl;
    os << "Number of rows    = " << GetRangeSpace().GetNumGlobalElements() << endl;
    os << "Number of columns = " << GetRangeSpace().GetNumGlobalElements() << endl;
    os << "Flop count        = " << GetFlops() << endl;
    os << "Cumulative time   = " << GetTime() << endl;
    if (GetTime() != 0.0)
      os << "MFlops rate       = " << 1.0e-6 * GetFlops() / GetTime() << endl;
    else
      os << "MFlops rate       = 0.0" << endl;
    os << endl;
  }

  StackPop();

  return(os);

}

void InverseOperator::Destroy()
{
  Op_.Reshape();
  RCPRowMatrix_ = Teuchos::null;
  RCPData_      = Teuchos::null;
  RCPMLPrec_    = Teuchos::null;
}

} // namespace MLAPI
#endif
