#ifndef MLAPI_MULTILEVELNONSYMMETRIC_H
#define MLAPI_MULTILEVELNONSYMMETRIC_H

#include "ml_include.h"
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
\class MultiLevelNonSymmetricSA

\brief Black-box multilevel smoothed aggregation preconditioner for non-symmetric matrices.

\date Last updated on Mar-05.

*/

class MultiLevelNonSymmetricSA : public BaseOperator, public CompObject, public TimeObject {

public:

  // @{ \name Constructors and destructors

  //! Constructs the hierarchy for given Operator and parameters.
  MultiLevelNonSymmetricSA(const Operator FineMatrix, 
                           Teuchos::ParameterList& List,
                           const string RestrictionType = "classic",
                           const bool ConstructNow = true) :
    IsComputed_(false)
  {
    SetList(List);
    SetMaxLevels(List.get("max levels", 10));
    SetRestrictionType(RestrictionType);
    if (GetMaxLevels() <= 0)
      ML_THROW("Invalid number of levels, " + GetString(GetMaxLevels()), -1);

    ResizeArrays();
    A(0) = FineMatrix;

    if (ConstructNow) Compute();
  }

  //! Destructor.
  virtual ~MultiLevelNonSymmetricSA()
  { }

  // @}
  // @{ \name Set and Get methods

  //! Returns a copy of the internally stored domain space.
  const Space GetOperatorDomainSpace() const 
  {
    return(A(0).GetDomainSpace());
  }

  //! Returns a copy of the internally stored range space.
  const Space GetOperatorRangeSpace() const 
  {
    return(A(0).GetRangeSpace());
  }

  //! Returns a copy of the internally stored domain space.
  inline const Space GetDomainSpace() const 
  {
    return(A(0).GetDomainSpace());
  }

  //! Returns a copy of the internally stored range space.
  inline const Space GetRangeSpace() const 
  {
    return(A(0).GetRangeSpace());
  }

  //! Returns a reference to the restriction operator of level \c i.
  inline const Operator& R(const int i) const
  {
    return(R_[i]);
  }

  //! Returns a reference to the restriction operator of level \c i.
  inline Operator& R(const int i)
  {
    return(R_[i]);
  }

  //! Returns a reference to the operator of level \c i.
  inline const Operator& A(const int i) const
  {
    return(A_[i]);
  }

  //! Returns a reference to the operator of level \c i.
  inline Operator& A(const int i)
  {
    return(A_[i]);
  }

  //! Returns a reference to the prolongator operator of level \c i.
  inline const Operator& P(const int i) const
  {
    return(P_[i]);
  }

  //! Returns a reference to the prolongator operator of level \c i.
  inline Operator& P(const int i) 
  {
    return(P_[i]);
  }

  //! Returns a reference to the inverse operator of level \c i.
  inline const InverseOperator& S(const int i) const
  {
    return(S_[i]);
  }

  //! Returns a reference to the inverse operator of level \c i.
  inline InverseOperator& S(const int i)
  {
    return(S_[i]);
  }

  //! Returns the actual number of levels
  inline int GetMaxLevels() const
  {
    return(MaxLevels_);
  }

  //! Returns the actual number of levels
  inline void SetMaxLevels(const int MaxLevels)
  {
    MaxLevels_ = MaxLevels;
  }

  //! Returns \c true if the hierarchy has been successfully computed, \c false otherwise.
  inline bool IsComputed() const
  {
    return(IsComputed_);
  }

  //! Returns the smoothed type.
  inline string GetSmootherType()
  {
    return(List_.get("smoother: type", "symmetric Gauss-Seidel"));
  }

  //! Returns the coarse solver type.
  inline string GetCoarseType()
  {
    return(List_.get("coarse: type", "Amesos-KLU"));
  }

  //! Returns the prolongator damping factor.
  inline double GetDamping()
  {
    return(List_.get("aggregation: damping factor", 1.3333));
  }

  //! Returns the maximum coarse size.
  inline double GetMaxCoarseSize()
  {
    return(List_.get("coarse: max size", 32));
  }

  //! Reset the internally stored list
  inline void SetList(Teuchos::ParameterList& List)
  {
    List_ = List;
  }

  //! Returns the strategy to be used to build up the restriction.
  inline string GetRestrictionType() const
  {
    return(RestrictionType_);
  }

  //! Sets the strategy to be used to build up the restriction.
  inline void SetRestrictionType(const string R)
  {
    RestrictionType_ = R;
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
    MultiVector EmptySpace;
    MultiVector ThisNS     = List_.get("aggregation: null space", EmptySpace);
    int         NumPDEEqns = List_.get("PDE equations", 1);
    
    // build up the default null space
    if (ThisNS.GetNumVectors() == 0) {
      ThisNS.Reshape(A(0).GetDomainSpace(),NumPDEEqns);
      ThisNS = 0.0;
      for (int i = 0 ; i < ThisNS.GetMyLength() ; ++i)
        for (int j = 0 ; j < NumPDEEqns ;++j)
          if (i % NumPDEEqns == j)
            ThisNS(i,j) = 1.0;
    }

    MultiVector NextNS;     // contains the next-level null space

    if (GetPrintLevel()) {
      cout << endl;
      ML_print_line("-", 80);
      cout << "Non-symmetric Smoothed Aggregation" << endl;
      ML_print_line("-", 80);
    }

    int level;

    for (level = 0 ; level < GetMaxLevels() - 1 ; ++level) {

      Operator Ptent;
      double LambdaMax;

      if (level)
        List_.set("PDE equations", ThisNS.GetNumVectors());

      List_.set("workspace: current level", level);

      GetPtent(A(level), List_, ThisNS, Ptent, NextNS);
      ThisNS = NextNS;
      
      if (GetDamping() != 0.0) {

        LambdaMax = ComputeLambdaMax(A(level));

        MultiVector Diag = GetDiagonal(A(level));
        Diag.Reciprocal();
        Diag.Scale(GetDamping() / LambdaMax);
        Operator Dinv = GetDiagonal(Diag);
        Operator DinvA = Dinv * A(level);
        Operator I = GetIdentity(A(level).GetDomainSpace(),A(level).GetRangeSpace());
        Operator IminusA = I - DinvA;
        // could be replaced with the following line:
        // IminusA = GetJacobiIterationOperator(A,GetDamping() / LambdaMax);
        P(level) = IminusA * Ptent;
      }
      else {
        P(level) = Ptent;
        LambdaMax = -1.0;
      }

      if (GetPrintLevel()) 
        PrintInfo(level, NumPDEEqns, ThisNS.GetNumVectors(),
                  LambdaMax);

      BuildR(level);

      A(level + 1) = GetRAP(R(level),A(level),P(level));
      S(level).Reshape(A(level), GetSmootherType(), List_);

      // break if coarse matrix is below specified tolerance
      if (A(level + 1).GetNumGlobalRows() <= GetMaxCoarseSize()) {
        ++level;
        break;
      }
    }

    // set coarse solver
    S(level).Reshape(A(level), GetCoarseType(), List_);
    SetMaxLevels(level + 1);

    // set the label
    SetLabel("NS-SA, L = " + GetString(GetMaxLevels()) +
             ", smoother = " + GetSmootherType());

    if (GetPrintLevel()) PrintInfo(level);

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
    if (level == GetMaxLevels() - 1) {
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
      os << endl;
      os << "*** MLAPI::MultiLevelSA, label = `" << GetLabel() << "'" << endl;
      os << endl;
      os << "Number of levels = " << GetMaxLevels() << endl;
      os << "Flop count       = " << GetFlops() << endl;
      os << "Cumulative time  = " << GetTime() << endl;
      if (GetTime() != 0.0)
        os << "MFlops rate      = " << 1.0e-6 * GetFlops() / GetTime() << endl;
      else
        os << "MFlops rate      = 0.0" << endl;
      os << endl;
    }
    return(os);
  }

  // @}

private:

  double ComputeLambdaMax(Operator& A)
  {

    double LambdaMax;
    string EigenAnalysis = List_.get("eigen-analysis: type", "Anorm");

    if (EigenAnalysis == "Anorm")
      LambdaMax = MaxEigAnorm(A,true);
    else if (EigenAnalysis == "cg")
      LambdaMax = MaxEigCG(A,true);
    else if (EigenAnalysis == "power-method")
      LambdaMax = MaxEigPowerMethod(A,true);
    else
      ML_THROW("incorrect parameter (" + EigenAnalysis + ")", -1);

    return(LambdaMax);
  }

  void ResizeArrays()
  {
    A_.resize(GetMaxLevels());
    R_.resize(GetMaxLevels());
    P_.resize(GetMaxLevels());
    S_.resize(GetMaxLevels());
  }

  void PrintInfo(const int level, const int NumPDEEqns, 
                 const int NullSpaceDimension, double LambdaMax)
  {
    ML_print_line("-", 80);
    cout << "current working level   = " << level << endl;
    cout << "number of global rows   = " << A(level).GetNumGlobalRows() << endl;
    cout << "number of global nnz    = " << A(level).GetNumGlobalNonzeros() << endl;
    cout << "threshold               = " << List_.get("aggregation: threshold", 0.0) << endl;
    cout << "number of PDE equations = " << NumPDEEqns << endl;
    cout << "null space dimension    = " << NullSpaceDimension << endl;
    cout << "omega                   = " << GetDamping() << endl;
    if (LambdaMax != -1.0) {
      cout << "lambda max              = " << LambdaMax << endl;
      cout << "damping factor          = " << GetDamping() / LambdaMax << endl;
    }
    cout << "smoother type           = " << GetSmootherType() << endl;
    cout << "relaxation sweeps       = " << List_.get("smoother: sweeps", 1) << endl;
    cout << "smoother damping        = " << List_.get("smoother: damping factor", 0.67) << endl;
  }

  void PrintInfo(const int level)
  {
    ML_print_line("-", 80);
    cout << "final level             = " << level << endl;
    cout << "number of global rows   = " << A(level).GetNumGlobalRows() << endl;
    cout << "number of global nnz    = " << A(level).GetNumGlobalNonzeros() << endl;
    cout << "coarse solver           = " << GetCoarseType() << endl;
    cout << "time spent in constr.   = " << GetTime() << " (s)" << endl;
    ML_print_line("-", 80);
  }

  //! Builds the restriction for level \c level.
  void BuildR(int level) 
  {
    StackPush();

    if (GetRestrictionType() == "classic")
      R(level) = GetTranspose(P(level));
    else
      ML_THROW("Invalid choice (" + GetRestrictionType() + ")", -1);

    StackPop();
  }

  //! Maximum number of levels.
  int MaxLevels_;
  //! Contains the hierarchy of operators.
  vector<Operator> A_;
  //! Contains the hierarchy of restriction operators.
  vector<Operator> R_;
  //! Contains the hierarchy of prolongator operators.
  vector<Operator> P_;
  //! Contains the hierarchy of inverse operators.
  vector<InverseOperator> S_;
  //! Contains a copy of the input list.
  Teuchos::ParameterList List_;
  //! \c true if the hierarchy has been successfully computed, \c false otherwise.
  bool IsComputed_;
  //! Strategy to define the restriction operator
  string RestrictionType_;

};

} // namespace MLAPI

#endif
