#ifndef MLAPI_MULTILEVEL_H
#define MLAPI_MULTILEVEL_H

#include "ml_include.h"
#include "Teuchos_ParameterList.hpp"
#include "MLAPI_Operator.h"
#include "MLAPI_DoubleVector.h"
#include "MLAPI_InverseOperator.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_Preconditioner.h"
#include "MLAPI_Container.h"
#include "MLAPI_Workspace.h"
#include <vector>

#include "ml_agg_genP.h"

static int MLAPI_count = 0;

namespace MLAPI {

class MultiLevel : public Preconditioner {

public:

  MultiLevel(const Operator FineMatrix,
             Teuchos::ParameterList& MLList)
  {

    FineMatrix_ = FineMatrix;

    double Damping = MLList.get("aggregation: damping factor", 1.333);
    string EigenAnalysis = MLList.get("eigen-analysis: type", "Anorm");

    MaxLevels_ = MLList.get("max levels",2);
    if (MaxLevels_ <= 0) {
      cerr << "Value of `max levels' not valid (" << MaxLevels_ << ")" << endl;
      throw("invalid parameter");
    }

    H_ = new Container[MaxLevels_];

    H_[0].SetA(FineMatrix);

    double LambdaMax;
    Operator A;
    Operator C;
    Operator R;
    Operator P;
    Operator Ptent;
    Operator IminusA;
    InverseOperator S;

    int level;
    for (level = 0 ; level < MaxLevels_ - 1 ; ++level) 
    {
      if (PrintLevel()) {
        cout << endl;
        cout << "Building level " << level << "..." << endl;
        cout << endl;
      }

      A = H_[level].A();
      Ptent = BuildP(A,MLList);
      LambdaMax = MaxEigenvalue(A,EigenAnalysis,true);

      if (PrintLevel()) {
        cout << endl;
        cout << "Prolongator/Restriction smoother (level " << level 
          << ") : damping factor = " << Damping / LambdaMax << endl;
        cout << "Prolongator/Restriction smoother (level " << level
          << ") : (= " << Damping << " / " << LambdaMax << ")" << endl;
        cout << endl;
      }

#if 0
      DoubleVector Diag(A.DomainSpace());
      Diag = Diagonal(A);
      cout << MLAPI_count << endl;
      ++MLAPI_count;
      Diag = (Damping / LambdaMax) / Diag;
      Operator Dinv = Diagonal(A.DomainSpace(),A.RangeSpace(),Diag);
      Operator I = Identity(A.DomainSpace(),A.RangeSpace());
      Operator DinvA = Dinv * A;
      //Operator IminusA = I - (Damping / LambdaMax) * DinvA;
      Operator IminusA = I - DinvA;
#else
      struct ML_AGG_Matrix_Context widget = {0};
      IminusA = JacobiIterationOperator(A,Damping,LambdaMax,&widget);
#endif

      P = IminusA * Ptent;

      R = Transpose(P);
      C = RAP(R,A,P);
      // build smoothers
      S.Reshape(A,"SGS",MLList);
      // put operators in hierarchy
      H_[level].SetR(R);
      H_[level].SetP(P);
      H_[level + 1].SetA(C);
      // put smoothers in hierarchy
      H_[level].SetS(S);

      // break if coarse matrix is below specified tolerance
      if (C.DomainSpace().NumGlobalElements() < 32) {
        ++level;
        break;
      }
    }

    // set coarse solver
    S.Reshape(H_[level].A(),"Amesos",MLList);
    H_[level].SetS(S);
    MaxLevels_ = level + 1;
  }

  virtual ~MultiLevel()
  {
    delete H_;
  }

  int Solve(const DoubleVector& b_f, DoubleVector& x_f) const
  {
    SolveMultiLevel(b_f,x_f,0);
    return(0);
  }

  int SolveMultiLevel(const DoubleVector& b_f,DoubleVector& x_f, int level) const 
  {
    if (level == MaxLevels_ - 1) {
      x_f = S(level) / b_f;
      return(0);
    }

    DoubleVector r_f(P(level).RangeSpace());
    DoubleVector r_c(P(level).DomainSpace());
    DoubleVector z_c(P(level).DomainSpace());

    // apply pre-smoother
    x_f = S(level) / b_f;
    // new residual
    r_f = b_f - A(level) * x_f;
    // restrict to coarse
    r_c = R(level) * r_f;
    // solve coarse problem
    SolveMultiLevel(r_c,z_c,level + 1);
    // prolongate back and add to solution
    x_f = x_f + P(level) * z_c;
    // apply post-smoother
    S(level).ApplyInverse(b_f,x_f); 

    return(0);
  }

  const Space& DomainSpace() const {
    return(FineMatrix_.DomainSpace());
  }

  const Space& RangeSpace() const {
    return(FineMatrix_.RangeSpace());
  }

  const Operator& R(const int i) const
  {
    return(H_[i].R());
  }

  const Operator& A(const int i) const
  {
    return(H_[i].A());
  }

  const Operator& P(const int i) const
  {
    return(H_[i].P());
  }

  const InverseOperator& S(const int i) const
  {
    return(H_[i].S());
  }

private:
  Operator FineMatrix_;
  Container* H_;
  int MaxLevels_;

};

} // namespace MLAPI

#endif
