#ifndef MLAPI_MULTILEVEL_H
#define MLAPI_MULTILEVEL_H

#include "ml_include.h"
#include "Teuchos_ParameterList.hpp"
#include "MLAPI_Operator.h"
#include "MLAPI_DoubleVector.h"
#include "MLAPI_InverseOperator.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_Preconditioner.h"
#include "MLAPI_Workspace.h"
#include <vector>

#include "ml_agg_genP.h"

namespace MLAPI {

/*!
\class MultiLevel

\brief Black-box multilevel smoothed aggregation preconditioner.

\author Marzio Sala, SNL 9214

\date Last updated on 07-Jan-05

*/

class MultiLevel : public Preconditioner {

public:

  //! Constructs the hierarchy for given Operator and parameters.
  /*! Constructs the multilevel hierarchy.
   *  \param FineMatrix (In) - Fine-level matrix
   *
   *  \param MLList (In) - Teuchos list containing all the parameters
   *                that can affect the construction of the hierarchy,
   *                the definition of smoothers and coarse solver.
   *
   *  \warning: Only a limited subset of parameters supported by
   *  ML_Epetra::MultiLevelPreconditioner is recognized by this class.
   *  In particular:
   *  - only IFPACK smoothers can be defined;
   *  - only Amesos coarse solver are supported.
   *  - only increasing hierarchies.
   */
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

    A_.resize(MaxLevels_);
    R_.resize(MaxLevels_);
    P_.resize(MaxLevels_);
    S_.resize(MaxLevels_);

    // work on increasing hierarchies only.
    A_[0] = FineMatrix;

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

      // only to simplify the notation
      A = A_[level];
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
      // put operators and inverse in hierarchy
      R_[level] = R;
      P_[level] = P;
      A_[level + 1] = C;
      S_[level] = S;

      // break if coarse matrix is below specified tolerance
      if (C.DomainSpace().NumGlobalElements() < 32) {
        ++level;
        break;
      }
    }

    // set coarse solver
    S.Reshape(A_[level],"Amesos",MLList);
    S_[level] = S;
    MaxLevels_ = level + 1;
  }

  //! Destructor.
  virtual ~MultiLevel()
  { }

  //! Applies the preconditioner to b_f with starting solution x_f.
  int Solve(const DoubleVector& b_f, DoubleVector& x_f) const
  {
    SolveMultiLevel(b_f,x_f,0);
    return(0);
  }

  //! Recursively called method, core of the multi level preconditioner.
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

  //! Returns a copy of the internally stored domain space.
  const Space DomainSpace() const {
    return(FineMatrix_.DomainSpace());
  }

  //! Returns a copy of the internally stored range space.
  const Space RangeSpace() const {
    return(FineMatrix_.RangeSpace());
  }

  //! Returns a reference to the restriction operator of level \c i.
  const Operator& R(const int i) const
  {
    return(R_[i]);
  }

  //! Returns a reference to the operator of level \c i.
  const Operator& A(const int i) const
  {
    return(A_[i]);
  }

  //! Returns a reference to the prolongator operator of level \c i.
  const Operator& P(const int i) const
  {
    return(P_[i]);
  }

  //! Returns a reference to the inverse operator of level \c i.
  const InverseOperator& S(const int i) const
  {
    return(S_[i]);
  }

private:
  //! Maximum number of levels.
  int MaxLevels_;
  //! Fine-level matrix.
  Operator FineMatrix_;
  //! Contains the hierarchy of operators.
  vector<Operator> A_;
  //! Contains the hierarchy of restriction operators.
  vector<Operator> R_;
  //! Contains the hierarchy of prolongator operators.
  vector<Operator> P_;
  //! Contains the hierarchy of inverse operators.
  vector<InverseOperator> S_;

};

} // namespace MLAPI

#endif
