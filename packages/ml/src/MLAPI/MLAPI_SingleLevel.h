#ifndef MLAPI_SINGLELEVEL_H
#define MLAPI_SINGLELEVEL_H

#include "ml_include.h"
#include "Teuchos_ParameterList.hpp"
#include "MLAPI_Operator.h"
#include "MLAPI_DoubleVector.h"
#include "MLAPI_Smoother.h"
#include "MLAPI_JacobiSmoother.h"
#include "MLAPI_SGSSmoother.h"
#include "MLAPI_AmesosSmoother.h"
#include "MLAPI_Preconditioner.h"
#include "MLAPI_Workspace.h"

namespace MLAPI {

class SingleLevel : public Preconditioner {

public:

  SingleLevel(const Operator* FineMatrix,
              Teuchos::ParameterList& MLList) :
    FineMatrix_(*FineMatrix),
    Smoother_(0)
  {

    Smoother_ = new SGSSmoother(*FineMatrix,MLList);
  }

  ~SingleLevel()
  {
    if (Smoother_)
      delete Smoother_;
  }

  int Solve(const DoubleVector& b_f, DoubleVector& x_f) const
  {
    Smoother& S = *Smoother_;
    x_f = S / b_f;
    return(0);
  }

  const Space& DomainSpace() const {
    return(FineMatrix_.DomainSpace());
  }

  const Space& RangeSpace() const {
    return(FineMatrix_.RangeSpace());
  }

private:
  const Operator& FineMatrix_;
  Smoother* Smoother_;

};
} // namespace MLAPI

#endif
