#ifndef MLAPI_MULTILEVEL_H
#define MLAPI_MULTILEVEL_H

#include "ml_config.h"
#include "MLAPI_Operator.h"
#include "MLAPI_Vector.h"
#include "MLAPI_Smoother.h"
#include "MLAPI_JacobiSmoother.h"
#include "MLAPI_SGSSmoother.h"
#include "MLAPI_AmesosSmoother.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_Preconditioner.h"
#include "MLAPI_SingleLevel.h"
#include "MLAPI_Utils.h"

namespace MLAPI {

class MultiLevel : public Preconditioner {

public:

  MultiLevel(const Operator* FineMatrix,
             Teuchos::ParameterList& MLList) :
    FineMatrix_(*FineMatrix)
  {

    MaxLevels_ = MLList.get("max levels",2);
    Hierarchy_ = new SingleLevel[10];

    Hierarchy_[0].SetA(FineMatrix);

    int i;
    for (i = 0 ; i < MaxLevels_ - 1 ; ++i) {
      const Operator* Amat = Hierarchy_[i].A();
      Operator* Pmat = BuildP(*Amat,MLList);
      Operator* Rmat = BuildTranspose(*Pmat);
      Operator* Cmat = RAP(*Rmat,*Amat,*Pmat);
      // build smoothers
      Smoother* PreSmoother = new SGSSmoother(*Amat,MLList);
      Smoother* PostSmoother = new SGSSmoother(*Amat,MLList);
      // put operators in hierarchy
      Hierarchy_[i].SetR(Rmat);
      Hierarchy_[i].SetP(Pmat);
      Hierarchy_[i+1].SetA(Cmat);
      // put smoothers in hierarchy
      Hierarchy_[i].SetPreSmoother(PreSmoother);
      Hierarchy_[i].SetPostSmoother(PostSmoother);
      // break if coarse matrix is below specified tolerance
      if (Cmat->DomainSpace().NumGlobalElements() < 128)
        break;
    }
    // set coarse solver
    Smoother* PostSmoother = new AmesosSmoother(*(Hierarchy_[i].A()),MLList);
    Hierarchy_[i].SetPostSmoother(PostSmoother);
    MaxLevels_ = i + 1;
  }

  ~MultiLevel()
  {
    delete Hierarchy_;
  }

  int Solve(const Vector& b_f, Vector& x_f) const
  {
    SolveMultiLevel(b_f,x_f,0);
    return(0);
  }

  int SolveMultiLevel(const Vector& b_f,Vector& x_f, int level) const 
  {
    const Smoother* Spre = Hierarchy_[level].PreSmoother();
    const Smoother* Spost = Hierarchy_[level].PostSmoother();

    if (level == MaxLevels_ - 1) {
      x_f = (*Spost) / b_f;
#ifdef CHECK
      const Operator* A_f = Hierarchy_[level].A();
      Vector xxx(A_f->DomainSpace());
      xxx = b_f - (*A_f) * x_f;
      cout << xxx * xxx << endl;
#endif
      return(0);
    }

    const Operator* A = Hierarchy_[level].A();
    const Operator* R = Hierarchy_[level].R();
    const Operator* P = Hierarchy_[level].P();

    Vector r_f(P->RangeSpace());
    Vector r_c(P->DomainSpace());
    Vector z_c(P->DomainSpace());

    // apply post-smoother
    if (Spre != 0)
      x_f = (*Spre) / b_f;
    else
      x_f = b_f;
    // new residual
    r_f = b_f - (*A) * x_f;
    // restrict to coarse
    r_c = (*R) * r_f;
    // solve coarse problem
    SolveMultiLevel(r_c,z_c,level + 1);
    // prolongate back and add to solution
    x_f = x_f + (*P) * z_c;
    // apply post-smoother
    if (Spost != 0)
      Spost->ApplyInverse(b_f,x_f);

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
  SingleLevel* Hierarchy_;
  int MaxLevels_;

};
} // namespace MLAPI

#endif
