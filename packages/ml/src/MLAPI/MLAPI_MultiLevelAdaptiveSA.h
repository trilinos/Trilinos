#ifndef MLAPI_MULTILEVELADAPTIVESA_H
#define MLAPI_MULTILEVELADAPTIVESA_H

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
\class MultiLevelAdaptiveSA

\brief Black-box multilevel adaptive smoothed aggregation preconditioner.

Still to do:
- store the structure of the aggregates for all phases.
- set threshold to zero??

\author Marzio Sala, Ray Tuminaro, Jonathan Hu, Michael Gee, Marian Brezina.

\date Last updated on Mar-05.

*/

class MultiLevelAdaptiveSA : public BaseOperator, public CompObject, public TimeObject {

public:

  // @{ \name Constructors and destructors

  //! Constructs the hierarchy for given Operator and parameters.
  MultiLevelAdaptiveSA(const Operator FineMatrix, Teuchos::ParameterList& List,
                       const bool ConstructNow = true)
  {
    FineMatrix_ = FineMatrix;
    List_ = List;
    if (ConstructNow) Compute();
  }

  // @}
  // @{ \name Hierarchy creation methods

  //! Creates an hierarchy using the provided or default null space.
  void Compute() 
  {

    ResetTimer();
    StackPush();

    // get parameter from the input list
    int         MaxLevels     = List_.get("max levels", 10);
    int         MaxCoarseSize = List_.get("coarse: max size", 32);
    MultiVector EmptySpace;
    int         NumPDEEqns    = List_.get("PDE equations", 1);
    string      SmootherType  = List_.get("smoother: type", "symmetric Gauss-Seidel");
    string      CoarseType    = List_.get("coarse: type", "Amesos-KLU");
    
    // retrive null space
    MultiVector ThisNS = GetNullSpace();

    // build up the default null space
    if (ThisNS.GetNumVectors() == 0) {
      ThisNS.Reshape(FineMatrix_.GetDomainSpace(),NumPDEEqns);
      ThisNS = 0.0;
      for (int i = 0 ; i < ThisNS.GetMyLength() ; ++i)
        for (int j = 0 ; j < NumPDEEqns ;++j)
          if (i % NumPDEEqns == j)
            ThisNS(i,j) = 1.0;

      SetNullSpace(ThisNS);
    }

    MultiVector NextNS;     // contains the next-level null space

    A_.resize(MaxLevels);
    R_.resize(MaxLevels);
    P_.resize(MaxLevels);
    S_.resize(MaxLevels);

    // work on increasing hierarchies only.
    A_[0] = FineMatrix_;

    Operator A;
    Operator C;
    Operator R;
    Operator P;
    Operator Ptent;
    InverseOperator S;

    int level;

    for (level = 0 ; level < MaxLevels - 1 ; ++level) {

      // only an alias
      A = A_[level];

      if (level)
        List_.set("PDE equations", ThisNS.GetNumVectors());

      if (GetPrintLevel()) {
      ML_print_line("-", 80);
        cout << "current working level   = " << level << endl;
        cout << "number of global rows   = " << A.GetNumGlobalRows() << endl;
        cout << "number of global nnz    = " << A.GetNumGlobalNonzeros() << endl;
        cout << "threshold               = " << List_.get("aggregation: threshold", 0.0) << endl;
        cout << "number of PDE equations = " << List_.get("PDE equations",-666) << endl;
        cout << "null space dimension    = " << ThisNS.GetNumVectors() << endl;
      }

      // load current level into database
      List_.set("workspace: current level", level);

      GetSmoothedP(A, List_, ThisNS, P, NextNS);
      ThisNS = NextNS;

      R = GetTranspose(P);
      C = GetRAP(R,A,P);
      // build smoothers
      S.Reshape(A, SmootherType, List_);

      // put operators and inverse in hierarchy
      R_[level    ] = R;
      P_[level    ] = P;
      A_[level + 1] = C;
      S_[level    ] = S;

      // break if coarse matrix is below specified tolerance
      if (C.GetNumGlobalRows() <= MaxCoarseSize) {
        ++level;
        break;
      }
    }
    List_.set("PDE equations",NumPDEEqns);

    // set coarse solver
    S.Reshape(A_[level], CoarseType, List_);
    S_[level] = S;
    MaxLevels_ = level + 1;

    // set the label
    SetLabel("SA, L = " + GetString(MaxLevels_) +
             ", smoother = " + SmootherType);

    if (GetPrintLevel()) {
      ML_print_line("-", 80);
      cout << "final level             = " << level << endl;
      cout << "number of global rows   = " << A_[level].GetNumGlobalRows() << endl;
      cout << "number of global nnz    = " << A_[level].GetNumGlobalNonzeros() << endl;
      cout << "coarse solver           = " << CoarseType << endl;
      cout << "time spent in constr.   = " << GetTime() << " (s)" << endl;
      ML_print_line("-", 80);
    }

    StackPop();
    // FIXME: update flops!
    UpdateTime();

  }

  //! Setup the adaptive multilevel hierarchy.
  /* Computes the multilevel hierarchy as specified by the user.
   *
   * \param UseDefault - (In) if \c true, the first call to Compute()
   *                     uses the default null space. Otherwise,
   *                     one null space component is computed using
   *                     SetupNullSpace().
   *
   * \param AdditionalCandidates - (In) Number of candidates, that is the
   *                     number of null space components that will be
   *                     computed using IncrementNullSpace().
   */
  void AdaptCompute(const bool UseDefault, const int AdditionalCandidates) 
  {

    if (UseDefault) 
      Compute();
    else {
      // compute the first guy, supposing that no null space
      SetupNullSpace();
      Compute();
    }

    for (int i = 0 ; i < AdditionalCandidates ; ++i) {
      IncrementNullSpace();
      Compute();
    }

  }

  //! Computes the first component of the null space.
  void SetupNullSpace() 
  {

    int    MaxLevels     = List_.get("max levels", 10);
    int    MaxCoarseSize = List_.get("coarse: max size", 1);
    string SmootherType  = List_.get("smoother: type", "symmetric Gauss-Seidel");
    string CoarseType    = List_.get("coarse: type", "Amesos-KLU");
    int    NumPDEs       = List_.get("PDE equations", 1);

    MultiVector NS(A_[0].GetDomainSpace());
    MultiVector NewNS;
    for (int v = 0 ; v < NS.GetNumVectors() ; ++v)
      NS.Random(v);

    NS = (NS + 1.0) / 2.0;

    // zero out everything except for first dof on every  node 
    for (int j=0; j < NS.GetMyLength(); ++j)
    {
      if (j % NumPDEs != 0)
        NS(j) = 0.;
    }

    vector<Operator>        A(MaxLevels);
    vector<Operator>        R(MaxLevels);
    vector<Operator>        P(MaxLevels);
    vector<InverseOperator> S(MaxLevels);

    // run pre-smoother
    MultiVector F(A_[0].GetDomainSpace());
    F = 0.0;

    S[0].Reshape(A_[0], SmootherType, List_);
    S[0].Apply(F, NS);

    double MyEnergyBefore = sqrt((A_[0] * NS) * NS);
    if (MyEnergyBefore == 0.0) {
      SetNullSpace(NewNS);
      return;
    }

    //compare last two iterates on fine level
    int SweepsBefore = List_.get("smoother: sweeps",1);
    List_.set("smoother: sweeps",1);
    S[0].Reshape(A_[0], SmootherType, List_);
    S[0].Apply(F, NS);
    double MyEnergyAfter = sqrt((A_[0] * NS) * NS);
    if (MyEnergyAfter/MyEnergyBefore < 0.1) {
      SetNullSpace(NewNS);
      return;
    }

    List_.set("smoother: sweeps",SweepsBefore);

    A[0] = A_[0];

    int level;
    for (level = 0 ; level < MaxLevels - 2 ; ++level) {

      if (level)
        List_.set("PDE equations", NS.GetNumVectors());

      if (GetPrintLevel()) {
        ML_print_line("-", 80);
        cout << "current working level   = " << level << endl;
        cout << "number of global rows   = " 
          << A[level].GetDomainSpace().GetNumGlobalElements() << endl;
        cout << "number of PDE equations = " << List_.get("PDE equations",-666) << endl;
        cout << "null space dimension    = " << NS.GetNumVectors() << endl;
      }

      GetSmoothedP(A[level], List_, NS, P[level], NewNS);
      NS = NewNS;

      R[level] = GetTranspose(P[level]);
      A[level + 1] = GetRAP(R[level],A[level],P[level]);
      Operator C = A[level + 1];
      S[level + 1].Reshape(C,SmootherType,List_);

      // break if coarse matrix is below specified size
      if (A[level + 1].GetDomainSpace().GetNumGlobalElements() <= MaxCoarseSize) {
        ++level;
        cout << "test " << level << " Max levels = " << MaxLevels << endl;
        break;
      }

      MultiVector F(C.GetDomainSpace());
      F = 0.0;
      MyEnergyBefore = sqrt((C * NS) * NS);
      S[level + 1].Apply(F, NS);
      MyEnergyAfter = sqrt((C * NS) * NS);
      cout << "Energy before smoothing = " << MyEnergyBefore << endl;
      cout << "Energy after smoothing  = " << MyEnergyAfter << endl;
      if (pow(MyEnergyAfter/MyEnergyBefore,1.0/SweepsBefore) < 0.1) {
        ++level;
        break; 
      }
    }

    List_.set("PDE equations", NumPDEs);

    if (GetPrintLevel())
      ML_print_line("-", 80);

    // interpolate candidate back to fine level
    MaxLevels = level;
    for (int level = MaxLevels ; level > 0 ; --level) {
      NS = P[level - 1] * NS;
    }

    F.Reshape(A_[0].GetDomainSpace());
    F = 0.0;
    S[0].Apply(F, NS);

    double norm = NS.NormInf();
    NS.Scale(1.0 / norm);

    SetNullSpace(NS);
  }

  // Increments the null space dimension by one.
  void IncrementNullSpace()
  {
    ResetTimer();
    StackPush();

    // get parameter from the input list
    MultiVector EmptySpace;
    int         NumPDEEqns    = List_.get("PDE equations", 1);
    string      SmootherType  = List_.get("smoother: type", "symmetric Gauss-Seidel");
    string      CoarseType    = List_.get("coarse: type", "Amesos-KLU");
    
    MultiVector InputNS = GetNullSpace();

    if (InputNS.GetNumVectors() == 0)
      ML_THROW("Empty null space not allowed", -1);

    Operator A, C, R, P;
    InverseOperator S;

    int level;

    // create the NS with the Random new candidate in it
    // NCand is dimension of previous nullspace
    int NCand = InputNS.GetNumVectors();
    MultiVector ExpandedNS = InputNS;
    ExpandedNS.Append();

    // Extract a light-weight copy of the additional component
    MultiVector AdditionalNS = Extract(ExpandedNS, NCand);
    AdditionalNS = (AdditionalNS + 1.) / 2.0;

    //  zero out in the new candidate everybody but the (Ncand+1)'th guy
    if (NCand+1 <= NumPDEEqns)
    {
      for (int i=0; i< AdditionalNS.GetMyLength(); i++)
        if ( (i+1) % (NCand+1) != 0) 
          AdditionalNS(i) = 0;
    }

   // run the current V-cycle on the candidate
   MultiVector b0(AdditionalNS.GetVectorSpace());

   int Nits = 15;
   for (int i=0; i<Nits; i++)
     SolveMultiLevelSA(b0,AdditionalNS,0);

   double NormFirst = ExpandedNS.NormInf(0);
   AdditionalNS.Scale(NormFirst / AdditionalNS.NormInf());

   // NOTE: I need this instruction, since the original pointer
   // in AdditionalNS could have been changed under the hood
   for (int i=0; i<AdditionalNS.GetMyLength(); i++)
     ExpandedNS(i,NCand) = AdditionalNS(i);

   NumPDEEqns = List_.get("PDE equations", -666);

   // ===================== //
   // cycle over all levels //
   // ===================== //
   
   for (level = 0 ; level < GetMaxLevels() - 2 ; ++level) {

      // only an alias
      A = A_[level];

      // load current level into database
      List_.set("workspace: current level", level);

      MultiVector NewNS;
   
      // creates smoothed prolongator
      GetSmoothedP(A, List_, ExpandedNS, P, NewNS);
      ExpandedNS = NewNS;

      R = GetTranspose(P);
      C = GetRAP(R,A,P);

      List_.set("PDE equations", NewNS.GetNumVectors());

      // build smoothers
      if (level != 0)
        S_[level].Reshape(A, SmootherType, List_);

      // put operators and inverse in hierarchy
      R_[level    ] = R;
      P_[level    ] = P;
      A_[level + 1] = C;
      
      S_[level+1].Reshape(C, SmootherType, List_); 
      
      AdditionalNS.Reshape(ExpandedNS.GetVectorSpace(),NCand);
      for (int i=0; i<NCand; i++)
        for (int j=0; j<AdditionalNS.GetMyLength(); j++)
          AdditionalNS(j,i) = ExpandedNS(j,i);

      Operator Pbridge;
      GetSmoothedP(C, List_, AdditionalNS, Pbridge, NewNS);

      P_[level+1] = Pbridge;
      R_[level+1] = GetTranspose(Pbridge);

      AdditionalNS = Extract(ExpandedNS, NCand);

      double MyEnergyBefore = sqrt((C * AdditionalNS) * AdditionalNS);

      // FIXME scale with something: norm of the matrix, ...;
      if (MyEnergyBefore < 1e-10) {
        ++level;
        break;
      }
        
      b0.Reshape(AdditionalNS.GetVectorSpace());
      b0 = 0.;
      
      if (level)
        Nits = 5;

      for (int i=0; i<Nits; i++)
        SolveMultiLevelSA(b0,AdditionalNS,level+1);

      // get norm of the first NS component
      double maxabs = 0.;

      for (int i=0; i<ExpandedNS.GetMyLength(); ++i)
         if (maxabs < fabs(ExpandedNS(i,0))) maxabs = fabs(ExpandedNS(i,0));

      double MyEnergyAfter = sqrt((C * AdditionalNS) * AdditionalNS);

      // scale the new candidate
      double max = AdditionalNS.NormInf();
      AdditionalNS.Scale((maxabs/max));

      for (int i=0; i<AdditionalNS.GetMyLength(); i++)
        ExpandedNS(i,NCand) = AdditionalNS(i);
      
      //TODO coarse level
      cout << "adaptedCompute: EnergyBefore=" << MyEnergyBefore << endl;
      cout << "adaptedCompute: EnergyAfter =" << MyEnergyAfter << endl;

      // FIXME: still to do:
      // - scaling of the new computed component
      // - if MyEnergyAfter is zero, take the previous guy

      if (pow(MyEnergyAfter/MyEnergyBefore,1.0/Nits) < 0.1) {
        ++level;
        break;
      }
    }
    
    --level;
    //project back to fine level
    // FIXME: put scaling on every level in here?
    for (int i=level; i>=0 ; i--) {
      AdditionalNS = P_[i] * AdditionalNS;
    }

    InputNS.Append(AdditionalNS);
   
    // reset the number of equations
    List_.set("PDE equations", NumPDEEqns);

    // set the null space in the class
    SetNullSpace(InputNS);

    StackPop();
    // FIXME: update flops!
    UpdateTime();

  }

  //! Destructor.
  virtual ~MultiLevelAdaptiveSA()
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

  inline const MultiVector& GetNullSpace() const
  {
    return(NullSpace_);
  }

  inline void SetNullSpace(MultiVector& NullSpace)
  {
    NullSpace_ = NullSpace;
  }

  // @}
  // @{ \name Mathematical methods

  //! Applies the preconditioner to \c b_f, returns the result in \c x_f.
  int Apply(const MultiVector& b_f, MultiVector& x_f) const
  {
    ResetTimer();
    StackPush();

    SolveMultiLevelSA(b_f,x_f,0);

    StackPop();
    UpdateTime();

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
    S(level).Apply(b_f,x_f); 
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

  //! Returns the smoothed prolongator operator.
  void GetSmoothedP(Operator A, Teuchos::ParameterList& List, MultiVector& NS,
                    Operator& P, MultiVector& NewNS)
  {
    double LambdaMax;
    Operator IminusA;
    string EigenAnalysis = List.get("eigen-analysis: type", "Anorm");
    double Damping       = List.get("aggregation: damping", 1.333);

    Operator Ptent;

    GetPtent(A, List, NS, Ptent, NewNS);

    if (EigenAnalysis == "Anorm")
      LambdaMax = MaxEigAnorm(A,true);
    else if (EigenAnalysis == "cg")
      LambdaMax = MaxEigCG(A,true);
    else if (EigenAnalysis == "power-method")
      LambdaMax = MaxEigPowerMethod(A,true);
    else
      ML_THROW("incorrect parameter (" + EigenAnalysis + ")", -1);

    if (GetPrintLevel()) {
      cout << "omega                   = " << Damping << endl;
      cout << "lambda max              = " << LambdaMax << endl;
      cout << "damping factor          = " << Damping / LambdaMax << endl;
    }

    if (Damping != 0.0) {
      IminusA = GetJacobiIterationOperator(A, Damping / LambdaMax);
      P = IminusA * Ptent;
    }
    else
      P = Ptent;

    return;
  }

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
  Teuchos::ParameterList List_;
  //! Contains the current null space
  MultiVector NullSpace_;

}; // class MultiLevelAdaptiveSA

} // namespace MLAPI
#endif
