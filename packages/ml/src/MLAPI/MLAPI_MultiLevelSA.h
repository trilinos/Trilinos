#ifndef MLAPI_MULTILEVEL_H
#define MLAPI_MULTILEVEL_H

#include "ml_include.h"
#include "MLAPI_Operator.h"
#include "MLAPI_Operator_Utils.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_InverseOperator.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_Preconditioner.h"
#include "MLAPI_Workspace.h"
#include "MLAPI_NullSpace.h"
#include "MLAPI_Aggregation.h"
#include "MLAPI_DataBase.h"
#include "MLAPI_Eig.h"
#include <vector>

#include "ml_agg_genP.h"

namespace MLAPI {

/*!
\class MultiLevel

\brief Black-box multilevel smoothed aggregation preconditioner.

\author Marzio Sala, SNL 9214

\date Last updated on 07-Jan-05

*/

class MultiLevelSA : public Preconditioner {

public:

  //! Constructs the hierarchy for given Operator and parameters.
  /*! Constructs the multilevel hierarchy.
   *  \param FineMatrix (In) - Fine-level matrix
   *
   */
  MultiLevelSA(const Operator FineMatrix, AggregationDataBase& ADB,
               SmootherDataBase& SDB, CoarseSolverDataBase& CDB)
  {

    FineMatrix_ = FineMatrix;

    MaxLevels_           = ADB.GetMaxLevels();
    double Damping       = ADB.GetDamping();
    string EigenAnalysis = ADB.GetEigenAnalysisType();
    int MaxCoarseSize    = ADB.GetMaxCoarseSize();

    // setup the null space for the finest level
    // FIXME: now only default null space!
    // to get nullspace-stuff from GetOption
    int NullSpaceDimension = 1;
    NullSpace ThisNS(FineMatrix.GetDomainSpace(), NullSpaceDimension);
    NullSpace NextNS;     // contains the next-level null space

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
      if (GetPrintLevel()) {
        cout << endl;
        cout << "Building level " << level << "..." << endl;
        cout << endl;
      }

      // load current level into database
      ADB.SetCurrentLevel(level);

      // alias, to simplify the notation
      A = A_[level];
      BuildPtent(A, ADB, ThisNS, Ptent, NextNS);
      ThisNS = NextNS;
      
      if (EigenAnalysis == "Anorm")
        LambdaMax = MaxEigAnorm(A,true);
      else if (EigenAnalysis == "cg")
        LambdaMax = MaxEigCG(A,true);
      else if (EigenAnalysis == "power-method")
        LambdaMax = MaxEigPowerMethod(A,true);
      else
        ML_THROW("incorrect parameter (" + EigenAnalysis + ")", -1);

      if (GetPrintLevel()) {
        cout << endl;
        cout << "Prolongator/Restriction smoother (level " << level 
          << ") : damping factor = " << Damping / LambdaMax << endl;
        cout << "Prolongator/Restriction smoother (level " << level
          << ") : (= " << Damping << " / " << LambdaMax << ")" << endl;
        cout << endl;
      }

#if 0
      MultiVector Diag(A.GetDomainSpace());
      Diag = Diagonal(A);
      Diag = (Damping / LambdaMax) / Diag;
      Operator Dinv = Diagonal(A.GetDomainSpace(),A.GetRangeSpace(),Diag);
      Operator I = Identity(A.GetDomainSpace(),A.GetRangeSpace());
      Operator DinvA = Dinv * A;
      //Operator IminusA = I - (Damping / LambdaMax) * DinvA;
      Operator IminusA = I - DinvA;
#else
      struct ML_AGG_Matrix_Context widget = {0};
      IminusA = JacobiIterationOperator(A,Damping / LambdaMax,&widget);
#endif

      P = IminusA * Ptent;

      R = Transpose(P);
      C = RAP(R,A,P);
      // build smoothers
      S.Reshape(A,SDB);
      // put operators and inverse in hierarchy
      R_[level] = R;
      P_[level] = P;
      A_[level + 1] = C;
      S_[level] = S;

      // break if coarse matrix is below specified tolerance
      if (C.GetDomainSpace().GetNumGlobalElements() < MaxCoarseSize) {
        ++level;
        break;
      }
    }

    // set coarse solver
    S.Reshape(A_[level],CDB);
    S_[level] = S;
    MaxLevels_ = level + 1;

  }

  //! Destructor.
  virtual ~MultiLevelSA()
  { }

  //! Applies the preconditioner to b_f with starting solution x_f.
  int Solve(const MultiVector& b_f, MultiVector& x_f) const
  {
    SolveMultiLevelSA(b_f,x_f,0);
    return(0);
  }

  //! Applies the preconditioner to b_f with starting solution x_f.
  int Solve(const MultiVector& b_f, MultiVector& x_f,
            const int FinestLevel = 0) const
  {
    SolveMultiLevelSA(b_f,x_f,FinestLevel);
    return(0);
  }

  //! Recursively called core of the multi level preconditioner.
  int SolveMultiLevelSA(const MultiVector& b_f,MultiVector& x_f, int level) const 
  {
    if (level == MaxLevels_ - 1) {
      x_f = S(level) * b_f;
      return(0);
    }

    MultiVector r_f(P(level).GetRangeSpace());
    MultiVector r_c(P(level).GetDomainSpace());
    MultiVector z_c(P(level).GetDomainSpace());

    // apply pre-smoother
    x_f = S(level) * b_f;
    // new residual
    r_f = b_f - A(level) * x_f;
    // restrict to coarse
    r_c = R(level) * r_f;
    // solve coarse problem
    SolveMultiLevelSA(r_c,z_c,level + 1);
    // prolongate back and add to solution
    x_f = x_f + P(level) * z_c;
    // apply post-smoother
    S(level).ApplyInverse(b_f,x_f); 

    return(0);
  }

  //! Returns a copy of the internally stored domain space.
  const Space GetDomainSpace() const {
    return(FineMatrix_.GetDomainSpace());
  }

  //! Returns a copy of the internally stored range space.
  const Space GetRangeSpace() const {
    return(FineMatrix_.GetRangeSpace());
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

  //! Prints basic information about \c this preconditioner.
  std::ostream& Print(std::ostream& os, 
                      const bool verbose = true) const
  {
    if (GetMyPID() == 0) {
      os << "MLAPI::MultiLevelSA, label = `" << GetLabel() << "'" << endl;
      os << endl;
      os << "Number of levels = " << GetMaxLevels() << endl;
      os << endl;
    }
    return(os);
  }

  //! Returns the actual number of levels
  int GetMaxLevels() const
  {
    return(MaxLevels_);
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
  //! Contains the hierarchy of inverse operators.

};

} // namespace MLAPI

#endif
