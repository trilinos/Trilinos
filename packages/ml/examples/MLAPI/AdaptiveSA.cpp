
//@HEADER
// ************************************************************************
// 
//               ML: A Multilevel Preconditioner Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include "ml_config.h"
#include "ml_common.h"
#ifdef HAVE_ML_MLAPI
#include "MLAPI.h"

using namespace Teuchos;
using namespace MLAPI;

void GetSmoothedP(Operator A, ParameterList& List, MultiVector& NS,
                  Operator& P, MultiVector& NewNS);


class GeneralMultiLevelSA : public BaseOperator, public CompObject, public TimeObject {

public:

  // @{ \name Constructors and destructors

  //! Constructs the hierarchy for given Operator and parameters.
  GeneralMultiLevelSA(const Operator FineMatrix, Teuchos::ParameterList& List,
                      const bool ConstructNow = true)
  {
    FineMatrix_ = FineMatrix;
    List_ = List;
    if (ConstructNow) Compute();
  }

  void Compute() {

    ResetTimer();

    // get parameter from the input list
    int         MaxLevels     = List_.get("max levels", 10);
    double      Damping       = List_.get("aggregation: damping factor", 1.3333);
    string      EigenAnalysis = List_.get("eigen-analysis: type", "Anorm");
    int         MaxCoarseSize = List_.get("coarse: max size", 32);
    MultiVector EmptySpace;
    MultiVector ThisNS        = List_.get("aggregation: null space", EmptySpace);
    int         NumPDEEqns    = List_.get("PDE equations", 1);
                cout << "blah " << NumPDEEqns << endl;
    string      SmootherType  = List_.get("smoother: type", "symmetric Gauss-Seidel");
    string      CoarseType    = List_.get("coarse: type", "Amesos-KLU");
    
    // build up the default null space
    if (ThisNS.GetNumVectors() == 0) {
      ThisNS.Reshape(FineMatrix_.GetDomainSpace(),NumPDEEqns);
      ThisNS = 0.0;
      for (int i = 0 ; i < ThisNS.GetMyLength() ; ++i)
        for (int j = 0 ; j < NumPDEEqns ;++j)
          if (i % NumPDEEqns == j)
            ThisNS(i,j) = 1.0;
    }

    MultiVector NextNS;     // contains the next-level null space

    A_.resize(MaxLevels);
    R_.resize(MaxLevels);
    P_.resize(MaxLevels);
    S_.resize(MaxLevels);

    // work on increasing hierarchies only.
    A_[0] = FineMatrix_;

    double LambdaMax;
    Operator A;
    Operator C;
    Operator R;
    Operator P;
    Operator Ptent;
    Operator IminusA;
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

      GetPtent(A, List_, ThisNS, Ptent, NextNS);
      ThisNS = NextNS;
      
      if (Damping) {

        if (EigenAnalysis == "Anorm")
          LambdaMax = MaxEigAnorm(A,true);
        else if (EigenAnalysis == "cg")
          LambdaMax = MaxEigCG(A,true);
        else if (EigenAnalysis == "power-method")
          LambdaMax = MaxEigPowerMethod(A,true);
        else
          ML_THROW("incorrect parameter (" + EigenAnalysis + ")", -1);

#if 0
        MultiVector Diag = GetDiagonal(A);
        Diag.Reciprocal();
        Diag.Scale(Damping / LambdaMax);
        Operator Dinv = GetDiagonal(Diag);
        Operator DinvA = Dinv * A;
        Operator I = GetIdentity(A.GetDomainSpace(),A.GetRangeSpace());
        Operator IminusA = I - DinvA;
#else
        IminusA = GetJacobiIterationOperator(A,Damping / LambdaMax);
#endif
        P = IminusA * Ptent;
      }
      else {
        P = Ptent;
        LambdaMax = -1.0;
      }

      if (GetPrintLevel()) {
        cout << "omega                   = " << Damping << endl;
        if (LambdaMax != -1.0) {
          cout << "lambda max              = " << LambdaMax << endl;
          cout << "damping factor          = " << Damping / LambdaMax << endl;
        }
        cout << "smoother type           = " << SmootherType << endl;
        cout << "relaxation sweeps       = " << List_.get("smoother: sweeps", 1) << endl;
        cout << "smoother damping        = " << List_.get("smoother: damping factor", 0.67) << endl;
      }

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

    // FIXME: update flops!
    UpdateTime();

    if (GetPrintLevel()) {
      ML_print_line("-", 80);
      cout << "final level             = " << level << endl;
      cout << "number of global rows   = " << A_[level].GetNumGlobalRows() << endl;
      cout << "number of global nnz    = " << A_[level].GetNumGlobalNonzeros() << endl;
      cout << "coarse solver           = " << CoarseType << endl;
      cout << "time spent in constr.   = " << GetTime() << " (s)" << endl;
      ML_print_line("-", 80);
    }

  }

  void AdaptedCompute() {

    ResetTimer();

    // get parameter from the input list
    int         MaxLevels     = List_.get("max levels", 10);
    double      Damping       = List_.get("aggregation: damping factor", 1.3333);
    string      EigenAnalysis = List_.get("eigen-analysis: type", "Anorm");
    int         MaxCoarseSize = List_.get("coarse: max size", 32);
    MultiVector EmptySpace;
    MultiVector ThisNS        = List_.get("aggregation: null space", EmptySpace);
    int         NumPDEEqns    = List_.get("PDE equations", 1);
    string      SmootherType  = List_.get("smoother: type", "symmetric Gauss-Seidel");
    string      CoarseType    = List_.get("coarse: type", "Amesos-KLU");
    
    // build up the default null space
    if (ThisNS.GetNumVectors() == 0) {
      ThisNS.Reshape(FineMatrix_.GetDomainSpace(),NumPDEEqns);
      ThisNS = 0.0;
      for (int i = 0 ; i < ThisNS.GetMyLength() ; ++i)
        for (int j = 0 ; j < NumPDEEqns ;++j)
          if (i % NumPDEEqns == j)
            ThisNS(i,j) = 1.0;
    }

    MultiVector NextNS;     // contains the next-level null space
    
    /* no need to resize, already exist
    A_.resize(MaxLevels);
    R_.resize(MaxLevels);
    P_.resize(MaxLevels);
    S_.resize(MaxLevels);*/

    // work on increasing hierarchies only.
    //A_[0] = FineMatrix_;

    double LambdaMax;
    Operator A;
    Operator C;
    Operator R;
    Operator P;
    Operator Ptent;
    Operator IminusA;
    InverseOperator S;

    int level;

    // create the NS with the Random new candidate in it
    // NCand is dimension of previous nullspace
    int NCand=ThisNS.GetNumVectors();
    MultiVector CandNS(A_[0].GetRangeSpace(),NCand+1);
    CandNS.Random();

    for (int i=0; i< CandNS.GetMyLength(); i++)
      CandNS(i,NCand) = (CandNS(i,NCand)+1.0)/2.0;
    //CandNS = (CandNS + 1.) / 2.0;

    for (int j=0; j< NCand; j++)
    for (int i=0; i< CandNS.GetMyLength(); i++)
       CandNS(i,j) = ThisNS(i,j);

    //  zero out in the new candidate everybody but the (Ncand+1)'th guy
    if (NCand+1 <= NumPDEEqns)
    {
      for (int i=0; i< CandNS.GetMyLength(); i++)
        if ( (i+1) % (NCand+1) != 0) 
          CandNS(i,NCand) = 0;
    }

   // run the current V-cycle on the candidate
   MultiVector CopyCandNS;
   CopyCandNS.Reshape(A_[0].GetRangeSpace());
   for (int i=0; i<CandNS.GetMyLength(); i++)
     CopyCandNS(i) = CandNS(i,NCand);
   MultiVector b0(CopyCandNS.GetVectorSpace());

   cout << "NumPDEEqns = " << NumPDEEqns << endl;

   /*
    MATLABStream matlab("output.m");
    CopyCandNS.SetLabel("xbefore");
    matlab << CopyCandNS;

    b0.SetLabel("rhsbefore");
    matlab << b0;
   */
   
   int Nits = 100;
   for (int i=0; i<Nits; i++)
     SolveMultiLevelSA(b0,CopyCandNS,0);

   /*
    CopyCandNS.SetLabel("xafter");
    matlab << CopyCandNS;

    b0.SetLabel("rhsafter");
    matlab << b0;

   exit(0);
   */

   for (int i=0; i<CopyCandNS.GetMyLength(); i++)
     CandNS(i,NCand) = CopyCandNS(i);

   cout << "MaxLevels_" << MaxLevels_ << endl;

   MultiVector NewNS;
   
   NumPDEEqns = List_.get("PDE equations", -666);

   for (level = 0 ; level < MaxLevels_ - 2 ; ++level) {

      // only an alias
      A = A_[level];

      // load current level into database
      List_.set("workspace: current level", level);

      GetSmoothedP(A, List_, CandNS, P, NewNS);
      CandNS = NewNS;

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
      
      CopyCandNS.Reshape(CandNS.GetVectorSpace(),NCand);
      for (int i=0; i<NCand; i++)
        for (int j=0; j<CopyCandNS.GetMyLength(); j++)
          CopyCandNS(j,i) = CandNS(j,i);

      Operator Pbridge;
      MultiVector junk;
      GetSmoothedP(C, List_, CopyCandNS, Pbridge, junk);

      P_[level+1] = Pbridge;
      R_[level+1] = GetTranspose(Pbridge);

      CopyCandNS.Reshape(CopyCandNS.GetVectorSpace());
      for (int i=0; i<CopyCandNS.GetMyLength(); i++)
        CopyCandNS(i) = CandNS(i,NCand);

      double MyEnergyBefore = sqrt((C * CopyCandNS) * CopyCandNS);

      b0.Reshape(CopyCandNS.GetVectorSpace());
      b0 = 0.;
      
      MultiVector anotherjunk = CopyCandNS;
      cout << "Nits=" << Nits << endl;

      for (i=0; i<Nits; i++)
        SolveMultiLevelSA(b0,CopyCandNS,level+1);

      anotherjunk = anotherjunk - CopyCandNS;
      cout << "Norm of DELTA=" << anotherjunk.Norm2() << endl;

      for (int i=0; i<CopyCandNS.GetMyLength(); i++)
        CandNS(i,NCand) = CopyCandNS(i);
      
      //TODO coarse level
      double MyEnergyAfter = sqrt((C * CopyCandNS) * CopyCandNS);
      cout << "adaptedCompute: EnergyBefore=" << MyEnergyBefore << endl;
      cout << "adaptedCompute: EnergyAfter =" << MyEnergyAfter << endl;


      if (pow(MyEnergyAfter/MyEnergyBefore,1.0/Nits) < 0.1)
        break;
    }
exit(0);
    //project back to fine level
    for (int i=level; i>=0 ; level--)
      CopyCandNS = P_[level] * CopyCandNS;
   
    //XXX NCand++

    List_.set("PDE equations", NumPDEEqns);

    // set coarse solver
    S.Reshape(A_[level], CoarseType, List_);
    S_[level] = S;
    MaxLevels_ = level + 1;

    // set the label
    SetLabel("SA, L = " + GetString(MaxLevels_) +
             ", smoother = " + SmootherType);

    // FIXME: update flops!
    UpdateTime();

    if (GetPrintLevel()) {
      ML_print_line("-", 80);
      cout << "final level             = " << level << endl;
      cout << "number of global rows   = " << A_[level].GetNumGlobalRows() << endl;
      cout << "number of global nnz    = " << A_[level].GetNumGlobalNonzeros() << endl;
      cout << "coarse solver           = " << CoarseType << endl;
      cout << "time spent in constr.   = " << GetTime() << " (s)" << endl;
      ML_print_line("-", 80);
    }

  }


  //! Destructor.
  virtual ~GeneralMultiLevelSA()
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

  // @}
  // @{ \name Mathematical methods

  //! Applies the preconditioner to \c b_f, returns the result in \c x_f.
  int Apply(const MultiVector& b_f, MultiVector& x_f) const
  {
    ResetTimer();
    SolveMultiLevelSA(b_f,x_f,0);
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

};



void GetSmoothedP(Operator A, ParameterList& List, MultiVector& NS,
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

MultiVector GetTentativeNullSpace(Operator FineMatrix, 
                                  Teuchos::ParameterList& List)
                                  
{
  int    MaxLevels     = List.get("max levels", 10);
  int    MaxCoarseSize = List.get("coarse: max size", 1);
  string SmootherType  = List.get("smoother: type", "symmetric Gauss-Seidel");
  string CoarseType    = List.get("coarse: type", "Amesos-KLU");
  int    NumPDEs       = List.get("PDE equations", 1);

  MultiVector NS(FineMatrix.GetDomainSpace());
  MultiVector NewNS;
  NS.Random();

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
  MultiVector F(FineMatrix.GetDomainSpace());
  F = 0.0;

  S[0].Reshape(FineMatrix, SmootherType, List);
  S[0].Apply(F, NS);

  double MyEnergyBefore = sqrt((FineMatrix * NS) * NS);
  if (MyEnergyBefore == 0.0)
    return NewNS;
  //compare last two iterates on fine level
  int SweepsBefore = List.get("smoother: sweeps",1);
  List.set("smoother: sweeps",1);
  S[0].Reshape(FineMatrix, SmootherType, List);
  S[0].Apply(F, NS);
  double MyEnergyAfter = sqrt((FineMatrix * NS) * NS);
  if (MyEnergyAfter/MyEnergyBefore < 0.1)
    return NewNS;
  List.set("smoother: sweeps",SweepsBefore);

  A[0] = FineMatrix;

  int level;
  for (level = 0 ; level < MaxLevels - 2 ; ++level) {

    if (level)
      List.set("PDE equations", NS.GetNumVectors());

    if (GetPrintLevel()) {
      ML_print_line("-", 80);
      cout << "current working level   = " << level << endl;
      cout << "number of global rows   = " 
        << A[level].GetDomainSpace().GetNumGlobalElements() << endl;
      cout << "number of PDE equations = " << List.get("PDE equations",-666) << endl;
      cout << "null space dimension    = " << NS.GetNumVectors() << endl;
    }

    GetSmoothedP(A[level], List, NS, P[level], NewNS);
    NS = NewNS;

    R[level] = GetTranspose(P[level]);
    A[level + 1] = GetRAP(R[level],A[level],P[level]);
    Operator C = A[level + 1];
    S[level + 1].Reshape(C,SmootherType,List);

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

  List.set("PDE equations", NumPDEs);

  if (GetPrintLevel())
    ML_print_line("-", 80);

  // interpolate candidate back to fine level
  MaxLevels = level;
  for (int level = MaxLevels ; level > 0 ; --level) {
    NS = P[level - 1] * NS;
  }

  F.Reshape(FineMatrix.GetDomainSpace());
  F = 0.0;
  S[0].Apply(F, NS);

  //MATLABStream matlab("output.m");
  //matlab << NS;

  return(NS);
}

// ============== //
// example driver //
// ============== //

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  // Initialize the workspace and set the output level
  Init();

  try {

    int NX = 1000;

#if 0
     //Operator NonScaledA = GetShiftedLaplacian1D(NX, 0.99999);
     Operator NonScaledA = GetShiftedLaplacian2D(NX, NX, 0.99999);
    //Operator NonScaledA = ReadMatrix(argv[1]);

    // need to get the fine space, it will be used later
    Space FineSpace = NonScaledA.GetDomainSpace();
#endif


    // define the space for fine level vectors and operators.
    Space FineSpace(2*NX);

    DistributedMatrix MatA(FineSpace, FineSpace);

    // assembel the matrix on processor 0 only
    if (GetMyPID() == 0) {
      for (int i = 0 ; i < NX ; ++i) {
        MatA.SetElement(2*i, 2*i, 2.0);
        MatA.SetElement(2*i+1, 2*i+1, 2.0);
        if (i)
        {
          MatA.SetElement(2*i, 2*(i - 1), - 1.0);
          MatA.SetElement(2*i+1, 2*(i - 1)+1, - 1.0);
        }
        if (i != NX - 1) {
          MatA.SetElement(2*i, 2*(i + 1), - 1.0);
          MatA.SetElement(2*i+1, 2*(i + 1)+1, - 1.0);
        }
      }
    }
    MatA.FillComplete();
    
    // wrap MatA as an Operator
    Operator A(FineSpace, FineSpace, &MatA, false);

#if 0
    PrintSparsity(A);
    cout << A;
#endif

#if 0
    MultiVector Scale(FineSpace);
    Scale.Random();
    Scale = Scale + 1.0001;

    Operator S = GetDiagonal(Scale);
    Operator A = GetRAP(S, NonScaledA, S);
#else
    //Operator A = NonScaledA;
#endif

#if 0
    MultiVector EI, ER, EV;
    Eig(A, ER, EI, EV);

    cout << ER << endl;
#endif

    // information about the current null space
    MultiVector NS(FineSpace);
    int i;
    for (i=0; i<NX; i++)
    {
       NS(2*i)   = 1.0;
       NS(2*i+1) = 0.0;
    }
    
    Teuchos::ParameterList List;
    List.set("PDE equations", 2);
    List.set("aggregation: null space", NS);
    List.set("max levels", 10); 
    List.set("smoother: type", "symmetric Gauss-Seidel");
    List.set("smoother: sweeps", 10);
    List.set("smoother: damping factor", 1.0);
    List.set("coarse: type", "Amesos-KLU");
    List.get("coarse: max size", 4);
 
    // ================================================== //
    // this is `phase1' ==> computation of a tentative    //
    // null space, supposing that we have no information  //
    // about any null space vector (that is, not even the //
    // constant).                                         //
    // ================================================== //
    
#if 1
    MultiVector TNS;
    TNS = GetTentativeNullSpace(A, List);
    List.set("aggregation: null space", TNS);
#endif

    // create the multilevel preconditioner
    GeneralMultiLevelSA Prec(A, List);

    // create the candNS-hierarchy
    Prec.AdaptedCompute();

    // test the solver
    MultiVector LHS(FineSpace);
    MultiVector RHS(FineSpace);

    LHS.Random();
    RHS = 0.0;

    List.set("krylov: type", "fixed point");
    Krylov(A, LHS, RHS, Prec, List);
    
    Finalize(); 
  
  }
  catch (const int e) {
    cerr << "Caught integer exception, code = " << e << endl;
  }
  catch (...) {
    cerr << "Caught exception..." << endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return(0);
  
}

#else

#include "ml_include.h"

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("The ML API requires the following configurataion options:");
  puts("\t--enable-epetra");
  puts("\t--enable-teuchos");
  puts("\t--enable-ifpack");
  puts("\t--enable-amesos");
  puts("Please check your configure line.");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(0);
}




#endif // #if defined(HAVE_ML_MLAPI)
