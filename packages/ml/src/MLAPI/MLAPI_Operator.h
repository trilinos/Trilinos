#ifndef ML_OPERATOR_H
#define ML_OPERATOR_H
#include "ml_config.h"
#include <iostream>
#include "ml_operator.h"
#include "ml_epetra.h"
#include "ml_amesos.h"
#include "ml_epetra_utils.h"
#include "ml_amesos_wrap.h"
#include "MLAPI_Space.h"
#include "MLAPI_Vector.h"
#include "MLAPI_Utils.h"

using namespace std;

namespace MLAPI {

class Operator {

public:
  Operator(const Space& DomainSpace, const Space& RangeSpace,
           ML_Operator* Op, bool Ownership = true) :
    Operator_(Op),
    Ownership_(Ownership),
    ApplyTemp_(0)
  {
    RangeSize_ = RangeSpace.NumMyElements();
    RangeSpace_ = new Space(RangeSpace);
    DomainSize_ = DomainSpace.NumMyElements();
    DomainSpace_ = new Space(DomainSpace);
    ApplyTemp_ = new Vector(*RangeSpace_);
  }

  Operator(const Space& DomainSpace, const Space& RangeSpace,
           Epetra_RowMatrix& Matrix) :
    Operator_(0),
    Ownership_(true),
    ApplyTemp_(0)
  {
    RangeSize_ = RangeSpace.NumMyElements();
    RangeSpace_ = new Space(RangeSpace);
    DomainSize_ = DomainSpace.NumMyElements();
    DomainSpace_ = new Space(DomainSpace);
    ApplyTemp_ = new Vector(*RangeSpace_);

    Operator_ = ML_Operator_Create(MLAPI::GetMLComm());
    Epetra2MLMatrix(&Matrix, Operator_);
  }

  ~Operator()
  {
    if (Ownership_) {
      if (Operator_)
        ML_Operator_Destroy(&Operator_);
    }
    if (ApplyTemp_)
      delete ApplyTemp_;
  }

  int Apply(const Vector& lhs, Vector& rhs) const
  {
    if (Operator_ == 0)
      throw("Operator not set");

    int (*func)(ML_Operator*,int,double*,int,double*) = Operator_->matvec->func_ptr;

    (*func)(Operator_,(int)DomainSize_,(double*)&lhs[0],(int)RangeSize_,(double*)&rhs[0]);
    return(0);
  }

  Vector& Apply(const Vector& lhs) const
  {
    Apply(lhs,*ApplyTemp_);
    return(*ApplyTemp_);
  }

  Space& RangeSpace() const {
    return(*RangeSpace_);
  }

  Space& DomainSpace() const {
    return(*DomainSpace_);
  }

  ML_Operator* GetOperator() const
  {
    return(Operator_);
  }

private:
  int RangeSize_;
  int DomainSize_;
  Space* DomainSpace_;
  Space* RangeSpace_;
  ML_Operator* Operator_;
  mutable Vector* ApplyTemp_;
  bool Ownership_;
}; // Operator

Operator* RAP(const Operator& R, const Operator& A, 
              const Operator& P, int matrix_type = ML_CSR_MATRIX)
{
  ML_Operator* Rmat = R.GetOperator();
  ML_Operator* Amat = A.GetOperator();
  ML_Operator* Pmat = P.GetOperator();
  ML_Operator* result = 0;

  result = ML_Operator_Create (Rmat->comm);

  ML_rap(Rmat, Amat, Pmat, result, matrix_type);

  Operator* op = new Operator(P.DomainSpace(),P.DomainSpace(),
                              result);
  return(op);

}

#include "ml_aggregate.h"
#include "ml_agg_METIS.h"

Operator* BuildP(const Operator& A, Teuchos::ParameterList& List)
{
  ML_Aggregate* agg_object;
  ML_Aggregate_Create(&agg_object);
  ML_Aggregate_Set_MaxLevels(agg_object,2);
  ML_Aggregate_Set_StartLevel(agg_object,0);
  
  ML_Operator* ML_Ptent = 0;
  ML_Ptent = ML_Operator_Create(GetMLComm());
  Operator* Ptent = 0;
  ML_Aggregate_CoarsenUncoupled(agg_object, A.GetOperator(),
                                &ML_Ptent, GetMLComm());
  
  // FIXME: set null space
  //double* OldNullspace = List.get("nullspace", (double*)0);

  int NumMyElements = ML_Ptent->invec_leng;
  Space CoarseSpace(NumMyElements,A.DomainSpace().Comm());
  Ptent = new Operator(CoarseSpace,A.DomainSpace(),ML_Ptent,true);

  return(Ptent);

}

Operator* BuildTranspose(const Operator& A) 
{
  ML_Operator* ML_transp;
  ML_transp = ML_Operator_Create(GetMLComm());
  ML_Operator_Transpose_byrow(A.GetOperator(),ML_transp);

  Operator* transp = new Operator(A.RangeSpace(),A.DomainSpace(),
                                  ML_transp,true);
  return(transp);
}

} // namespace MLAPI
#endif // ML_OPERATOR_H

