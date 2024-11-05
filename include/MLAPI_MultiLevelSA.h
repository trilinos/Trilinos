#ifndef MLAPI_MULTILEVEL_H
#define MLAPI_MULTILEVEL_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

/*!
\file MLAPI_MultiLevelSA.h

\brief Standard smoothed aggregation multilevel preconditioner.

\author Marzio Sala, D-INFK/ETHZ.

\date Last updated on Mar-06.
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#include "ml_common.h"
#include "ml_agg_genP.h"
#include "MLAPI_Error.h"
#include "MLAPI_CompObject.h"
#include "MLAPI_TimeObject.h"
#include "MLAPI_Operator.h"
#include "MLAPI_Operator_Utils.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_InverseOperator.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_BaseOperator.h"
#include "MLAPI_Workspace.h"
#include "MLAPI_Aggregation.h"
#include "MLAPI_Eig.h"
#include <vector>

namespace MLAPI {

/*!
\class MultiLevelSA

\brief Black-box multilevel smoothed aggregation preconditioner.

\author Marzio Sala, SNL 9214

\date Last updated on Feb-05.

*/

class MultiLevelSA : public BaseOperator, public CompObject, public TimeObject {

public:

  // @{ \name Constructors and destructors

  //! Constructs the hierarchy for given Operator and parameters.
  MultiLevelSA(const Operator & FineMatrix, Teuchos::ParameterList& List,
               const bool ConstructNow = true) :
    IsComputed_(false)
  {
    FineMatrix_ = FineMatrix;
    List_ = List;
    if (ConstructNow) Compute();
  }

  //! Destructor.
  virtual ~MultiLevelSA()
  { }

  // @}
  // @{ \name Set and Get methods

  //! Returns a copy of the internally stored domain space.
  const Space GetOperatorDomainSpace() const
  {
    return(FineMatrix_.GetDomainSpace());
  }

  //! Returns a copy of the internally stored range space.
  const Space GetOperatorRangeSpace() const
  {
    return(FineMatrix_.GetRangeSpace());
  }

  //! Returns a copy of the internally stored domain space.
  inline const Space GetDomainSpace() const
  {
    return(FineMatrix_.GetDomainSpace());
  }

  //! Returns a copy of the internally stored range space.
  inline const Space GetRangeSpace() const
  {
    return(FineMatrix_.GetRangeSpace());
  }

  //! Returns a reference to the restriction operator of level \c i.
  inline const Operator& R(const int i) const
  {
    return(R_[i]);
  }

  //! Returns a reference to the operator of level \c i.
  inline const Operator& A(const int i) const
  {
    return(A_[i]);
  }

  //! Returns a reference to the prolongator operator of level \c i.
  inline const Operator& P(const int i) const
  {
    return(P_[i]);
  }

  //! Returns a reference to the inverse operator of level \c i.
  inline const InverseOperator& S(const int i) const
  {
    return(S_[i]);
  }

  //! Returns the actual number of levels
  inline int GetMaxLevels() const
  {
    return(MaxLevels_);
  }

  //! Returns \c true if the hierarchy has been successfully computed, \c false otherwise.
  inline bool IsComputed() const
  {
    return(IsComputed_);
  }

  // @}
  // @{ \name Mathematical methods

  //! Computes the hierarchy.
  void Compute()
  {
    ResetTimer();
    StackPush();
    IsComputed_ = false;

    // get parameter from the input list
    int         MaxLevels     = List_.get("max levels", 10);
    double      Damping       = List_.get("aggregation: damping factor", 1.3333);
    std::string      EigenAnalysis = List_.get("eigen-analysis: type", "Anorm");
    int         MaxCoarseSize = List_.get("coarse: max size", 32);
    MultiVector EmptySpace;
    MultiVector ThisNS        = List_.get("aggregation: null space", EmptySpace);
    int         NumPDEEqns    = List_.get("PDE equations", 1);
    std::string      SmootherType  = List_.get("smoother: type", "symmetric Gauss-Seidel");
    std::string      CoarseType    = List_.get("coarse: type", "Amesos-KLU");

    // build up the default null space
    if (ThisNS.GetNumVectors() == 0) {
      ThisNS.Reshape(FineMatrix_.GetDomainSpace(),NumPDEEqns);
      if (NumPDEEqns == 1)
        ThisNS = 1.0;
      else
      {
        ThisNS = 0.0;
        for (int i = 0 ; i < ThisNS.GetMyLength() ; ++i)
          for (int j = 0 ; j < NumPDEEqns ;++j)
            if (i % NumPDEEqns == j)
              ThisNS(i,j) = 1.0;
      }
    }

    MultiVector NextNS;     // contains the next-level null space

    A_.resize(MaxLevels);
    R_.resize(MaxLevels);
    P_.resize(MaxLevels);
    S_.resize(MaxLevels);

    // work on increasing hierarchies only.
    A_[0] = FineMatrix_;

    double LambdaMax;
    Operator Aop;
    Operator C;
    Operator Rop;
    Operator Pop;
    Operator Ptent;
    Operator IminusA;
    InverseOperator Sop;

    int level;

    for (level = 0 ; level < MaxLevels - 1 ; ++level) {

      // only an alias
      Aop = A_[level];

      if (level)
        List_.set("PDE equations", ThisNS.GetNumVectors());

      if (GetPrintLevel()) {
        ML_print_line("-", 80);
        std::cout << "current working level   = " << level << std::endl;
        std::cout << "number of global rows   = " << Aop.GetNumGlobalRows() << std::endl;
        std::cout << "number of global nnz    = " << Aop.GetNumGlobalNonzeros() << std::endl;
        std::cout << "threshold               = " << List_.get("aggregation: threshold", 0.0) << std::endl;
        std::cout << "number of PDE equations = " << NumPDEEqns << std::endl;
        std::cout << "null space dimension    = " << ThisNS.GetNumVectors() << std::endl;
      }

      // load current level into database
      List_.set("workspace: current level", level);

      GetPtent(Aop, List_, ThisNS, Ptent, NextNS);
      ThisNS = NextNS;

      if (Damping) {

        if (EigenAnalysis == "Anorm")
          LambdaMax = MaxEigAnorm(Aop,true);
        else if (EigenAnalysis == "cg")
          LambdaMax = MaxEigCG(Aop,true);
        else if (EigenAnalysis == "power-method")
          LambdaMax = MaxEigPowerMethod(Aop,true);
        else
          ML_THROW("incorrect parameter (" + EigenAnalysis + ")", -1);

#if 0
        MultiVector Diag = GetDiagonal(Aop);
        Diag.Reciprocal();
        Diag.Scale(Damping / LambdaMax);
        Operator Dinv = GetDiagonal(Diag);
        Operator DinvA = Dinv * Aop;
        Operator I = GetIdentity(Aop.GetDomainSpace(),Aop.GetRangeSpace());
        Operator IminusA = I - DinvA;
#else
        IminusA = GetJacobiIterationOperator(Aop,Damping / LambdaMax);
#endif
        Pop = IminusA * Ptent;
      }
      else {
        Pop = Ptent;
        LambdaMax = -1.0;
      }

      if (GetPrintLevel()) {
        std::cout << "omega                   = " << Damping << std::endl;
        if (LambdaMax != -1.0) {
          std::cout << "lambda max              = " << LambdaMax << std::endl;
          std::cout << "damping factor          = " << Damping / LambdaMax << std::endl;
        }
        std::cout << "smoother type           = " << SmootherType << std::endl;
        std::cout << "relaxation sweeps       = " << List_.get("smoother: sweeps", 1) << std::endl;
        std::cout << "smoother damping        = " << List_.get("smoother: damping factor", 0.67) << std::endl;
      }

      Rop = GetTranspose(Pop);
      C = GetRAP(Rop,Aop,Pop);
      // build smoothers
      Sop.Reshape(Aop, SmootherType, List_);

      // put operators and inverse in hierarchy
      R_[level    ] = Rop;
      P_[level    ] = Pop;
      A_[level + 1] = C;
      S_[level    ] = Sop;

      // break if coarse matrix is below specified tolerance
      if (C.GetNumGlobalRows() <= MaxCoarseSize) {
        ++level;
        break;
      }
    }

    // set coarse solver
    Sop.Reshape(A_[level], CoarseType, List_);
    S_[level] = Sop;
    MaxLevels_ = level + 1;

    // set the label
    SetLabel("SA, L = " + GetString(MaxLevels_) +
             ", smoother = " + SmootherType);

    if (GetPrintLevel()) {
      ML_print_line("-", 80);
      std::cout << "final level             = " << level << std::endl;
      std::cout << "number of global rows   = " << A_[level].GetNumGlobalRows() << std::endl;
      std::cout << "number of global nnz    = " << A_[level].GetNumGlobalNonzeros() << std::endl;
      std::cout << "coarse solver           = " << CoarseType << std::endl;
      std::cout << "time spent in constr.   = " << GetTime() << " (s)" << std::endl;
      ML_print_line("-", 80);
    }

    IsComputed_ = true;
    StackPop();

    // FIXME: update flops!
    UpdateTime();

  }

  //! Applies the preconditioner to \c b_f, returns the result in \c x_f.
  int Apply(const MultiVector& b_f, MultiVector& x_f) const
  {
    ResetTimer();
    StackPush();

    if (IsComputed() == false)
      ML_THROW("Method Compute() must be called before Apply()", -1);
    SolveMultiLevelSA(b_f,x_f,0);
    UpdateTime();

    StackPop();
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

    // reset flop counter
    S(level).SetFlops(0.0);
    A(level).SetFlops(0.0);
    R(level).SetFlops(0.0);
    P(level).SetFlops(0.0);

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
    S(level).Apply(b_f,x_f);

    UpdateFlops(2.0 * S(level).GetFlops());
    UpdateFlops(A(level).GetFlops());
    UpdateFlops(R(level).GetFlops());
    UpdateFlops(P(level).GetFlops());
    UpdateFlops(2.0 * x_f.GetGlobalLength());

    return(0);
  }

  // @}
  // @{ \name Miscellaneous methods

  //! Prints basic information about \c this preconditioner.
  std::ostream& Print(std::ostream& os,
                      const bool verbose = true) const
  {
    if (GetMyPID() == 0) {
      os << std::endl;
      os << "*** MLAPI::MultiLevelSA, label = `" << GetLabel() << "'" << std::endl;
      os << std::endl;
      os << "Number of levels = " << GetMaxLevels() << std::endl;
      os << "Flop count       = " << GetFlops() << std::endl;
      os << "Cumulative time  = " << GetTime() << std::endl;
      if (GetTime() != 0.0)
        os << "MFlops rate      = " << 1.0e-6 * GetFlops() / GetTime() << std::endl;
      else
        os << "MFlops rate      = 0.0" << std::endl;
      os << std::endl;
    }
    return(os);
  }

  // @}

private:

  //! Maximum number of levels.
  int MaxLevels_;
  //! Fine-level matrix.
  Operator FineMatrix_;
  //! Contains the hierarchy of operators.
  std::vector<Operator> A_;
  //! Contains the hierarchy of restriction operators.
  std::vector<Operator> R_;
  //! Contains the hierarchy of prolongator operators.
  std::vector<Operator> P_;
  //! Contains the hierarchy of inverse operators.
  std::vector<InverseOperator> S_;
  //! Contains a copy of the input list.
  Teuchos::ParameterList List_;
  //! \c true if the hierarchy has been successfully computed, \c false otherwise.
  bool IsComputed_;

};

} // namespace MLAPI

#endif
