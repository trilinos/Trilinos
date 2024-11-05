#ifndef MLAPI_MULTILEVELADAPTIVESA_H
#define MLAPI_MULTILEVELADAPTIVESA_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

/*!
\file MLAPI_MultiLevelAdaptiveSA.h

\brief Adaptive smoothed aggregation preconditioner.

\author Marzio Sala, Ray Tuminaro, Jonathan Hu, Michael Gee, Marian Brezina.

\date Last updated on Feb-05.
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
#include "MLAPI_MultiVector_Utils.h"
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

This class implements an adaptive smoothed aggregation preconditioner.
An example of usage is reported in file \ref ml_adaptivesa.
We note that the usage of this class is slightly different from that of
MultiLevelSA.

An instance of this class can be created as follows:
\code
int NumPDEEqns = 1;
int MaxLevels = 10;
MultiLevelAdaptiveSA Prec(FineMatrix, List, NumPDEEqns, MaxLevels);
\endcode

Important methods of this class:
- The number of PDE equations on the finest level can be queried using
  GetInputNumPDEEqns().
- GetNumPDEEqns() returns the number of PDE equations on the current level.
  This value can be set via SetNumPDEEqns().
- GetNullSpace() returns a reference to the internally stored null space;
  the null space is set using SetNullSpace().
- GetMaxLevels() returns the number of levels. If called before Compute(),
  GetMaxLevels() returns the maximum number of levels used in the
  constructor, otherwise returns the actual number of levels.
- GetSmootherType() and GetCoarseType() return the smoother and coarse type.
- The number of application of the cycle in IncrementNullSpace() is given
  by GetNumItersCoarse() and GetNumItersFine().
- Methods \c A(level), \c P(level), \c R(level) and \c S(level) return a
  reference to the internally stored operators.
- Method SetList() can be used at any time to reset the internally stored
  list.

The general usage is:
- Specify the null space using SetNullSpace(NS), where NS is a MultiVector,
  then compute the hierarchy using Compute(), or
- Compute the first component of the null space using SetupInitialNullSpace().
  This will define a single-vector null space, and store it using
  SetNullSpace(NS).
- When a non-empty null space is provided, the user can increment by one
  the dimension of the null space by calling IncrementNullSpace().
- Method AdaptCompute() performs all these operations.

\author Marzio Sala, Ray Tuminaro, Jonathan Hu, Michael Gee, Marian Brezina.

\date Last updated on Mar-05.

\todo store the structure of the aggregates for all phases.
\todo Current implementation supposes zero threshold.

*/

class MultiLevelAdaptiveSA : public BaseOperator, public CompObject, public TimeObject {

public:

  // @{ \name Constructors and destructors

  //! Constructs the hierarchy for given Operator and parameters.
  MultiLevelAdaptiveSA(const Operator & FineMatrix, Teuchos::ParameterList& List,
                       const int NumPDEEqns, const int MaxLevels = 20) :
    IsComputed_(false)
  {
    FineMatrix_ = FineMatrix;
    List_ = List;
    MaxLevels_  = MaxLevels;
    SetInputNumPDEEqns(NumPDEEqns);
    SetNumPDEEqns(NumPDEEqns);

    ResizeArrays(MaxLevels);
    A(0) = FineMatrix;
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
  inline Operator& R(const int i)
  {
    return(R_[i]);
  }

  //! Returns a reference to the restriction operator of level \c i.
  inline const Operator& R(const int i) const
  {
    return(R_[i]);
  }

  //! Returns a reference to the operator of level \c i.
  inline Operator& A(const int i)
  {
    return(A_[i]);
  }

  //! Returns a reference to the operator of level \c i.
  inline const Operator& A(const int i) const
  {
    return(A_[i]);
  }

  //! Returns a reference to the prolongator operator of level \c i.
  inline Operator& P(const int i)
  {
    return(P_[i]);
  }

  //! Returns a reference to the prolongator operator of level \c i.
  inline const Operator& P(const int i) const
  {
    return(P_[i]);
  }

  //! Returns a reference to the inverse operator of level \c i.
  inline InverseOperator& S(const int i)
  {
    return(S_[i]);
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

  //! Returns the actual number of levels
  inline void SetMaxLevels(const int MaxLevels)
  {
    MaxLevels_ = MaxLevels;
  }

  //! Gets a reference to the internally stored null space.
  inline const MultiVector GetNullSpace() const
  {
    return(NullSpace_);
  }

  //! Sets the null space multi-vector to \c NullSpace.
  inline void SetNullSpace(MultiVector& NullSpace)
  {
    NullSpace_ = NullSpace;
  }

  //! Returns \c true if the hierarchy has been successfully computed.
  inline bool IsComputed() const
  {
    return(IsComputed_);
  }

  //! Sets the internally stored list to \c List.
  inline void SetList(Teuchos::ParameterList& List)
  {
    List_ = List;
  }

  //! Returns the smoother solver type.
  inline std::string GetSmootherType()
  {
    return(List_.get("smoother: type", "symmetric Gauss-Seidel"));
  }

  //! Returns the coarse solver type.
  inline std::string GetCoarseType()
  {
    return(List_.get("coarse: type", "Amesos-KLU"));
  }

  //! Returns the number of PDE equations on the finest level.
  inline void SetInputNumPDEEqns(const int n)
  {
    NumPDEEqns_ = n;
  }

  //! Returns the number of PDE equations on the current level.
  inline int GetInputNumPDEEqns()
  {
    return(NumPDEEqns_);
  }

  //! Sets the number of PDE equations on the current level.
  inline int GetNumPDEEqns()
  {
    return(List_.get("PDE equations", 1));
  }

  inline void SetNumPDEEqns(const int NumPDEEqns)
  {
    List_.set("PDE equations", NumPDEEqns);
  }

  //! Returns the maximum allowed coarse size.
  inline int GetMaxCoarseSize()
  {
    return(List_.get("coarse: max size", 32));
  }

  //! Returns the maximum allowed reduction.
  inline double GetMaxReduction()
  {
    return(List_.get("adapt: max reduction", 0.1));
  }

  //! Returns the maximum number of applications on the coarser levels.
  inline int GetNumItersCoarse()
  {
    return(List_.get("adapt: iters coarse", 5));
  }

  //! Returns the maximum number of applications on the finest level.
  inline int GetNumItersFine()
  {
    return(List_.get("adapt: iters fine", 15));
  }

  //! Returns the multigrid preconditioner operator complexity.
  double GetComplexity()
  {
    double nnzFine = A_[0].GetNumGlobalNonzeros();
    double nnzTotal = nnzFine;
    for (int i = 1; i < GetMaxLevels(); i++) {
      nnzTotal += A_[i].GetNumGlobalNonzeros();
    }
    return nnzTotal / nnzFine;
  }


  // @}
  // @{ \name Hierarchy construction methods

  // ======================================================================
  //! Creates an hierarchy using the provided or default null space.
  // ======================================================================
  void Compute()
  {

    ResetTimer();
    StackPush();
    IsComputed_ = false;

    // get parameter from the input list
    SetNumPDEEqns(GetInputNumPDEEqns());

    // retrive null space
    MultiVector ThisNS = GetNullSpace();

    if (GetPrintLevel()) {
      std::cout << std::endl;
      ML_print_line("-", 80);
      std::cout << "Computing the hierarchy, input null space dimension = "
           << ThisNS.GetNumVectors() << std::endl;
    }

    // build up the default null space
    if (ThisNS.GetNumVectors() == 0) {
      ThisNS.Reshape(FineMatrix_.GetDomainSpace(),GetNumPDEEqns());
      ThisNS = 0.0;
      for (int i = 0 ; i < ThisNS.GetMyLength() ; ++i)
        for (int j = 0 ; j < GetNumPDEEqns() ;++j)
          if (i % GetNumPDEEqns() == j)
            ThisNS(i,j) = 1.0;

      SetNullSpace(ThisNS);
      if (GetPrintLevel()) {
        std::cout << "number of PDE equations = " << GetNumPDEEqns() << std::endl;
        std::cout << "null space dimension    = " << ThisNS.GetNumVectors() << std::endl;
      }

    }

    MultiVector NextNS;     // contains the next-level null space

    // work on increasing hierarchies only.
    A(0) = FineMatrix_;

    int level;

    for (level = 0 ; level < GetMaxLevels() - 1 ; ++level) {

      if (GetPrintLevel()) ML_print_line("-", 80);

      if (level)
        SetNumPDEEqns(ThisNS.GetNumVectors());

      // load current level into database
      List_.set("workspace: current level", level);

      GetSmoothedP(A(level), List_, ThisNS, P(level), NextNS);
      ThisNS = NextNS;

      R(level) = GetTranspose(P(level));
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
    SetLabel("SA, L = " + GetString(GetMaxLevels()) +
             ", smoother = " + GetSmootherType());

    if (GetPrintLevel()) ML_print_line("-", 80);

    IsComputed_ = true;

    StackPop();
    // FIXME: update flops!
    UpdateTime();

  }

  // ======================================================================
  //! Setup the adaptive multilevel hierarchy.
  /* Computes the multilevel hierarchy as specified by the user.
   *
   * \param UseDefaultOrSpecified - (In) if \c true, the first call to Compute()
   *                     uses either the default null space or the null space
   *                     that has been set using SetNullSpace().
   *                     If \c false, then one null space component is
   *                     computed using SetupInitialNullSpace().
   *
   * \param AdditionalCandidates - (In) Number of candidates, that is the
   *                     number of null space components that will be
   *                     computed using IncrementNullSpace(). If
   *                     \c "UseDefaultOrSpecified == false", the code computes
   *                     one additional candidate using
   *                     SetupInitialNullSpace(), and the remaining using
   *                     IncrementNullSpace().
   */
  // time is tracked within each method.
  // ======================================================================
  void AdaptCompute(const bool UseDefaultOrSpecified, int AdditionalCandidates)
  {

    StackPush();

    if (UseDefaultOrSpecified)
      Compute();
    else {
      SetupInitialNullSpace();
      Compute();
      AdditionalCandidates--;
    }

    for (int i = 0 ; i < AdditionalCandidates ; ++i) {
      IncrementNullSpace();
      Compute();
    }

    StackPop();
  }

  // ======================================================================
  //! Computes the first component of the null space.
  // ======================================================================
  void SetupInitialNullSpace()
  {
    ResetTimer();
    StackPush();

    SetNumPDEEqns(GetInputNumPDEEqns());

    if (GetPrintLevel()) {
      std::cout << std::endl;
      ML_print_line("-", 80);
      std::cout << "Computing the first null space component" << std::endl;
    }

    MultiVector NS(A_[0].GetDomainSpace());
    MultiVector NewNS;
    for (int v = 0 ; v < NS.GetNumVectors() ; ++v)
      NS.Random(v);

    NS = (NS + 1.0) / 2.0;

    // zero out everything except for first dof on every  node
    for (int j=0; j < NS.GetMyLength(); ++j)
    {
      if (j % GetNumPDEEqns() != 0)
        NS(j) = 0.;
    }

    // run pre-smoother
    MultiVector F(A(0).GetDomainSpace());
    F = 0.0;

    S(0).Reshape(A(0), GetSmootherType(), List_);
    S(0).Apply(F, NS);

    double MyEnergyBefore = sqrt((A(0) * NS) * NS);
    if (MyEnergyBefore == 0.0) {
      SetNullSpace(NewNS);
      return;
    }

    //compare last two iterates on fine level
    int SweepsBefore = List_.get("smoother: sweeps",1);
    List_.set("smoother: sweeps",1);
    S(0).Reshape(A(0), GetSmootherType(), List_);
    S(0).Apply(F, NS);
    double MyEnergyAfter = sqrt((A(0) * NS) * NS);
    if (MyEnergyAfter/MyEnergyBefore < GetMaxReduction()) {
      SetNullSpace(NewNS);
      return;
    }

    List_.set("smoother: sweeps",SweepsBefore);

    int level;
    for (level = 0 ; level < GetMaxLevels() - 2 ; ++level) {

      if (level) SetNumPDEEqns(NS.GetNumVectors());

      if (GetPrintLevel()) {
        ML_print_line("-", 80);
        std::cout << "current working level   = " << level << std::endl;
        std::cout << "number of global rows   = "
          << A(level).GetDomainSpace().GetNumGlobalElements() << std::endl;
        std::cout << "number of PDE equations = " << GetNumPDEEqns() << std::endl;
        std::cout << "null space dimension    = " << NS.GetNumVectors() << std::endl;
      }

      GetSmoothedP(A(level), List_, NS, P(level), NewNS);
      NS = NewNS;

      R(level) = GetTranspose(P(level));
      A(level + 1) = GetRAP(R(level),A(level),P(level));
      S(level + 1).Reshape(A(level + 1),GetSmootherType(),List_);

      // break if coarse matrix is below specified size
      if (A(level + 1).GetDomainSpace().GetNumGlobalElements() <= GetMaxCoarseSize()) {
        ++level;
        break;
      }

      MultiVector locF(A(level + 1).GetDomainSpace());
      locF = 0.0;
      MyEnergyBefore = sqrt((A(level + 1) * NS) * NS);
      S(level + 1).Apply(locF, NS);
      MyEnergyAfter = sqrt((A(level + 1) * NS) * NS);
      if (GetPrintLevel() == 0) {
        std::cout << "Energy before smoothing = " << MyEnergyBefore << std::endl;
        std::cout << "Energy after smoothing  = " << MyEnergyAfter << std::endl;
      }

      if (pow(MyEnergyAfter/MyEnergyBefore,1.0/SweepsBefore) < GetMaxReduction()) {
        ++level;
        break;
      }
    }

    if (GetPrintLevel())
      ML_print_line("-", 80);

    // interpolate candidate back to fine level
    int MaxLevels = level;
    for (int j = MaxLevels ; j > 0 ; --j) {
      NS = P(j - 1) * NS;
    }

    F.Reshape(A(0).GetDomainSpace());
    F = 0.0;
    S(0).Apply(F, NS);

    double norm = NS.NormInf();
    NS.Scale(1.0 / norm);

    SetNullSpace(NS);

    StackPop();
    UpdateTime();
  }

  //! Increments the null space dimension by one.
  bool IncrementNullSpace()
  {
    ResetTimer();
    StackPush();

    SetNumPDEEqns(GetInputNumPDEEqns());

    MultiVector InputNS = GetNullSpace();

    if (InputNS.GetNumVectors() == 0)
      ML_THROW("Empty null space not allowed", -1);

    if (GetPrintLevel()) {
      std::cout << std::endl;
      ML_print_line("-", 80);
      std::cout << "Incrementing the hierarchy, input null space dimension = "
           << InputNS.GetNumVectors() << std::endl;
    }

    int level;

    // =========================================================== //
    // InputNS is the currently available (and stored) null space. //
    // ExpandedNS is InputNS + AdditionalNS.                       //
    // NCand is dimension of previous nullspace.                   //
    // AdditionalNS is set to random between 0 and 1.0; however,   //
    // we might need to zero out in the new candidate everybody    //
    // but the (Ncand+1)'th guy                                    //
    // Once the new candidate is set, we run the current V-cycle   //
    // on it. NOTE: I need the final copy instruction, since the   //
    // original pointer  in AdditionalNS could have been changed   //
    // under the hood.                                             //
    // =========================================================== //

    int NCand = InputNS.GetNumVectors();
    MultiVector ExpandedNS = InputNS;
    ExpandedNS.Append();

    MultiVector AdditionalNS = Extract(ExpandedNS, NCand);
    AdditionalNS.Random();
    AdditionalNS = (AdditionalNS + 1.) / 2.0;

    if (NCand+1 <= GetNumPDEEqns())
    {
      for (int i=0; i< AdditionalNS.GetMyLength(); i++)
        if ( (i+1) % (NCand+1) != 0)
          AdditionalNS(i) = 0;
    }

    MultiVector b0(AdditionalNS.GetVectorSpace());

    for (int i=0; i< GetNumItersFine() ; i++)
      SolveMultiLevelSA(b0,AdditionalNS,0);

    double NormFirst = ExpandedNS.NormInf(0);
    AdditionalNS.Scale(NormFirst / AdditionalNS.NormInf());

    for (int i=0; i<AdditionalNS.GetMyLength(); i++)
      ExpandedNS(i,NCand) = AdditionalNS(i);

    // ===================== //
    // cycle over all levels //
    // ===================== //

    for (level = 0 ; level < GetMaxLevels() - 2 ; ++level) {

      if (GetPrintLevel()) ML_print_line("-", 80);

      List_.set("workspace: current level", level);

      // ======================================================= //
      // Create a new prolongator operator using the newly       //
      // available null space. NewNS is a temporary variable,   //
      // set to ExpandedNS for the next-level null space. We     //
      // also need to setup the smoother at this level (not on   //
      // the finest, since the finest-level matrix does not      //
      // change).                                                //
      // At this point, we stick the operators in the hierarchy. //
      // ======================================================= //

      MultiVector NewNS;

      GetSmoothedP(A(level), List_, ExpandedNS, P(level), NewNS);
      ExpandedNS = NewNS;

      R(level) = GetTranspose(P(level));
      A(level + 1) = GetRAP(R(level),A(level),P(level));

      SetNumPDEEqns(NewNS.GetNumVectors());

      if (level != 0)
        S(level).Reshape(A(level), GetSmootherType(), List_);

      S(level + 1).Reshape(A(level + 1), GetSmootherType(), List_);

      // ======================================================= //
      // Need to setup the bridge. We need to extract the NCand  //
      // components of the "old" null space, and set them in a   //
      // temporary variable, OldNS. NewNS is simply discarded.   //
      // ======================================================= //

      MultiVector OldNS = ExpandedNS;
      OldNS.Delete(NCand);

      Operator Pbridge;
      GetSmoothedP(A(level + 1), List_, OldNS, Pbridge, NewNS);

      P(level + 1) = Pbridge;
      R(level + 1) = GetTranspose(Pbridge);

      AdditionalNS = Duplicate(Extract(ExpandedNS, NCand));

      double MyEnergyBefore = sqrt((A(level + 1) * AdditionalNS) * AdditionalNS);

      // FIXME scale with something: norm of the matrix, ...;
      if (MyEnergyBefore < 1e-10) {
        ++level;
        break;
      }

      // ======================================================= //
      // run Nits_coarse cycles, using AdditionalNS as starting  //
      // solution, and a zero right hand-side.                   //
      // ======================================================= //

      b0.Reshape(AdditionalNS.GetVectorSpace());
      b0 = 0.;

      for (int i=0; i< GetNumItersCoarse() ; i++)
        SolveMultiLevelSA(b0,AdditionalNS,level+1);

      // ======================================================= //
      // Get the norm of the first null space component, then    //
      // analyze the energy after the application of the cycle.  //
      // If the energy after is zero, we have to check whether   //
      // the new guy is zero or not.                             //
      // Then, we scale the new candidate so that its largest    //
      // entry is of the same magniture of the largest entry of  //
      // the first component.                                    //
      // ======================================================= //

      double NormFirstComponent = ExpandedNS.NormInf(0);

      double MyEnergyAfter = sqrt((A(level + 1) * AdditionalNS) * AdditionalNS);

      if (MyEnergyAfter == 0.0) {
        if (AdditionalNS.NormInf() != 0.0) {
          for (int i=0; i<AdditionalNS.GetMyLength(); i++)
            ExpandedNS(i,NCand) = AdditionalNS(i);
        }
      }
      else {
        for (int i=0; i<AdditionalNS.GetMyLength(); i++) {
          ExpandedNS(i,NCand) = AdditionalNS(i);
        }
      }

      double NormExpanded = ExpandedNS.NormInf(NCand);
      ExpandedNS.Scale(NormFirstComponent / NormExpanded, NCand);

      if (GetPrintLevel() == 0) {
        std::cout << "energy before cycle =" << MyEnergyBefore << std::endl;
        std::cout << "energy after        =" << MyEnergyAfter << std::endl;
      }

      // FIXME: still to do:
      // - scaling of the new computed component

      if (pow(MyEnergyAfter/MyEnergyBefore,1.0 / GetNumItersCoarse()) < GetMaxReduction()) {
        ++level;
        break;
      }
    }

    --level;
    AdditionalNS = Extract(ExpandedNS, NCand);

    // ======================================================= //
    // project back to fine level the AdditionalNS vector.     //
    // Then, reset the number of PDE equations, and finally    //
    // set the null space of this object using SetNullSpace(). //
    // Note that at this point the hierarchy is broken, and    //
    // must be reconstructed using Compute().                  //
    // ======================================================= //

    // FIXME: put scaling on every level in here?
    for (int i=level; i>=0 ; i--) {
      AdditionalNS = P(i) * AdditionalNS;
    }

    InputNS.Append(AdditionalNS);

    SetNullSpace(InputNS);

    if (GetPrintLevel()) ML_print_line("-", 80);

    StackPop();
    UpdateTime();

    IsComputed_ = false;

    return(true);
  }

  // @}
  // @{ \name Mathematical methods

  // ======================================================================
  //! Applies the preconditioner to \c b_f, returns the result in \c x_f.
  // ======================================================================
  int Apply(const MultiVector& b_f, MultiVector& x_f) const
  {
    ResetTimer();
    StackPush();

    if (IsComputed() == false)
      ML_THROW("Method Compute() must be called", -1);

    SolveMultiLevelSA(b_f,x_f,0);

    StackPop();
    UpdateTime();

    return(0);
  }

  // ======================================================================
  //! Recursively called core of the multi level preconditioner.
  // ======================================================================
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

      for (int level = 0 ; level < GetMaxLevels() ; ++level) {
        ML_print_line("-", 80);
        std::cout << "Information for level   = " << level;
        std::cout << "number of global rows   = "
             << A(level).GetNumGlobalRows() << std::endl;
        std::cout << "number of global nnz    = "
             << A(level).GetNumGlobalNonzeros() << std::endl;
      }
      ML_print_line("-", 80);
    }
    return(os);
  }

  // @}

private:

  //! Returns the smoothed prolongator operator.
  void GetSmoothedP(const Operator & Aop, Teuchos::ParameterList& List, MultiVector& NS,
                    Operator& Pop, MultiVector& NewNS)
  {
    double LambdaMax;
    Operator IminusA;
    std::string EigenAnalysis = List.get("eigen-analysis: type", "Anorm");
    double Damping       = List.get("aggregation: damping", 1.333);

    Operator Ptent;

    GetPtent(Aop, List, NS, Ptent, NewNS);

    if (EigenAnalysis == "Anorm")
      LambdaMax = MaxEigAnorm(Aop,true);
    else if (EigenAnalysis == "cg")
      LambdaMax = MaxEigCG(Aop,true);
    else if (EigenAnalysis == "power-method")
      LambdaMax = MaxEigPowerMethod(Aop,true);
    else
      ML_THROW("incorrect parameter (" + EigenAnalysis + ")", -1);

    if (GetPrintLevel()) {
      std::cout << "omega                   = " << Damping << std::endl;
      std::cout << "lambda max              = " << LambdaMax << std::endl;
      std::cout << "damping factor          = " << Damping / LambdaMax << std::endl;
    }

    if (Damping != 0.0) {
      IminusA = GetJacobiIterationOperator(Aop, Damping / LambdaMax);
      Pop = IminusA * Ptent;
    }
    else
      Pop = Ptent;

    // fix the number of equations in Pop, so that GetRAP() will
    // get the correct number for C= RAP
    Pop.GetML_Operator()->num_PDEs = Ptent.GetML_Operator()->num_PDEs;

    return;
  }

  void ResizeArrays(const int MaxLevels)
  {
    A_.resize(MaxLevels);
    R_.resize(MaxLevels);
    P_.resize(MaxLevels);
    S_.resize(MaxLevels);
  }

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
  Teuchos::ParameterList List_;
  //! Contains the current null space
  MultiVector NullSpace_;
  //! \c true if a hierarchy has been successfully computed.
  bool IsComputed_;
  //! Number of PDE equations on the finest grid.
  int NumPDEEqns_;

}; // class MultiLevelAdaptiveSA

} // namespace MLAPI
#endif
