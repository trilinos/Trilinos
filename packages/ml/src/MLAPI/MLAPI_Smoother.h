#ifndef ML_SMOOTHER_H
#define ML_SMOOTHER_H
#include "ml_include.h"
#include "ml_RowMatrix.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_PointRelaxation.h"
#include "Ifpack_Amesos.h"
#include "MLAPI_Space.h"
#include "MLAPI_DoubleVector.h"

using namespace std;

namespace MLAPI {

class Smoother {

public:

  Smoother()
  {}
  
  Smoother(const Smoother& RHS)
  {
    Op_ = RHS.GetOperator();
    ApplyInverseTemp_.Reshape(Op_.RangeSpace());

    RowMatrix_ = RHS.RowMatrix();
    Data_ = RHS.Data();
  }

  void Reshape(const Operator& Op, const string Type,
               Teuchos::ParameterList& List)
  {
    Op_ = Op;
    ApplyInverseTemp_.Reshape(Op.RangeSpace());

    RowMatrix_ = Teuchos::rcp(new ML_Epetra::RowMatrix(Op.GetOperator(),
                                                       &GetEpetraComm()));
    Ifpack_Preconditioner* Prec;
    if (Type == "SGS") 
      Prec = new Ifpack_PointRelaxation(RowMatrix_.get());
    else if (Type == "Amesos") 
      Prec = new Ifpack_Amesos(RowMatrix_.get());

    Data_ = Teuchos::rcp(Prec);

    Data_->SetParameters(List);
    Data_->Initialize();
    Data_->Compute();
  }

  ~Smoother()
  {}

  const Operator& GetOperator() const
  {
    return(Op_);
  }

  Smoother& operator=(const Smoother& RHS)
  {
    if (this == &RHS)
      return(*this);

    Op_ = RHS.GetOperator();
    ApplyInverseTemp_.Reshape(Op_.RangeSpace());

    RowMatrix_ = RHS.RowMatrix();
    Data_ = RHS.Data();
  }

  int ApplyInverse(const DoubleVector& lhs, DoubleVector& rhs) const
  {
    Epetra_Vector elhs(View,RowMatrix_->OperatorDomainMap(),
                       (double*)&(lhs(0)));
    Epetra_Vector erhs(View,RowMatrix_->OperatorRangeMap(),
                       (double*)&(rhs(0)));

    Data_->ApplyInverse(elhs,erhs);
  }

  DoubleVector& ApplyInverse(const DoubleVector& lhs) const
  {
    ApplyInverseTemp_= lhs;
    Epetra_Vector elhs(View,RowMatrix_->OperatorDomainMap(),
                       (double*)&(lhs(0)));
    Epetra_Vector erhs(View,RowMatrix_->OperatorRangeMap(),
                       (double*)&(ApplyInverseTemp_)(0));

    Data_->ApplyInverse(elhs,erhs);
    return(ApplyInverseTemp_);
  }

  const Space& RangeSpace() const {
    return(Op_.RangeSpace());
  }

  const Space& DomainSpace() const {
    return(Op_.DomainSpace());
  }

  const Teuchos::RefCountPtr<ML_Epetra::RowMatrix> RowMatrix() const
  {
    return(RowMatrix_);
  }

  const Teuchos::RefCountPtr<Ifpack_Preconditioner> Data() const
  {
    return(Data_);
  }

private:
  Operator Op_;
  mutable DoubleVector ApplyInverseTemp_;
  Teuchos::RefCountPtr<ML_Epetra::RowMatrix> RowMatrix_;
  Teuchos::RefCountPtr<Ifpack_Preconditioner> Data_;

}; // Smoother

} // namespace MLAPI
#endif // ML_SMOOTHER_H
