
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
  MultiVector NS(FineMatrix.GetDomainSpace());
  MultiVector NewNS;
  NS.Random();

  NS = (NS + 1.0) / 2.0;
 
  int    MaxLevels     = List.get("max levels", 10);
  int    MaxCoarseSize = List.get("coarse: max size", 1);
  string SmootherType  = List.get("smoother: type", "symmetric Gauss-Seidel");
  string CoarseType    = List.get("coarse: type", "Amesos-KLU");
  int    NumPDEs       = List.get("PDE equations", 1);

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
      cout << "number of PDE equations = " << NumPDEs << endl;
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

    int NX = 200;

    // Operator NonScaledA = GetShiftedLaplacian1D(NX, 0.99999);
    // Operator NonScaledA = GetShiftedLaplacian2D(NX, NX, 0.99999);
    Operator NonScaledA = ReadMatrix(argv[1]);

    // need to get the fine space, it will be used later
    Space FineSpace = NonScaledA.GetDomainSpace();

#if 1
    MultiVector Scale(FineSpace);
    Scale.Random();
    Scale = Scale + 1.0001;

    Operator S = GetDiagonal(Scale);
    Operator A = GetRAP(S, NonScaledA, S);
#else
    Operator A = NonScaledA;
#endif

#if 0
    MultiVector EI, ER, EV;
    Eig(A, ER, EI, EV);

    cout << ER << endl;
#endif

    // information about the current null space
    MultiVector NS(FineSpace);
    NS = 1.0; // constant null space
    
    Teuchos::ParameterList List;
    List.set("PDE equations", 1);
    List.set("aggregation: null space", NS);
    List.set("max levels", 10); 
    List.set("smoother: type", "symmetric Gauss-Seidel");
    List.set("smoother: sweeps", 10);
    List.set("smoother: damping factor", 1.0);
    List.set("coarse: type", "Amesos-KLU");
    List.get("coarse: max size", 32);
 
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
    MultiLevelSA Prec(A, List);

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
