
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
using namespace Trilinos_Util;
using namespace MLAPI;

MultiVector GetTentativeNullSpace(Operator FineMatrix, 
                                  Teuchos::ParameterList& List,
                                  int NullSpaceDimension)
{
  MultiVector NS(FineMatrix.GetDomainSpace(), NullSpaceDimension);
  MultiVector NewNS;
  NS.Random();
 
  int    MaxLevels     = List.get("max levels", 10);
  double Damping       = List.get("aggregation: damping", 1.333);
  string EigenAnalysis = List.get("eigen-analysis: type", "Anorm");
  int    MaxCoarseSize = List.get("coarse: max size", 32);
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

  A[0] = FineMatrix;

  double LambdaMax;
  Operator Ptent;
  Operator IminusA;

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

    // FIXME: risky, NS is in input & output...
    GetPtent(A[level], List, NS, Ptent, NewNS);
    NS = NewNS;

    if (EigenAnalysis == "Anorm")
      LambdaMax = MaxEigAnorm(A[level],true);
    else if (EigenAnalysis == "cg")
      LambdaMax = MaxEigCG(A[level],true);
    else if (EigenAnalysis == "power-method")
      LambdaMax = MaxEigPowerMethod(A[level],true);
    else
      ML_THROW("incorrect parameter (" + EigenAnalysis + ")", -1);

    if (GetPrintLevel()) {
      cout << "omega                   = " << Damping << endl;
      cout << "lambda max              = " << LambdaMax << endl;
      cout << "damping factor          = " << Damping / LambdaMax << endl;
      cout << "smoother type           = " << SmootherType << endl;
    }

    if (Damping != 0.0) {
      IminusA = GetJacobiIterationOperator(A[level],Damping / LambdaMax);
      P[level] = IminusA * Ptent;
    }
    else
      P[level] = Ptent;

    R[level] = GetTranspose(P[level]);
    A[level + 1] = GetRAP(R[level],A[level],P[level]);
    Operator C = A[level + 1];
    S[level + 1].Reshape(C,SmootherType,List);

    MultiVector F(C.GetDomainSpace());
    F = 0.0;
    double MyEnergyBefore = (C * NS) * NS;
    S[level + 1].Apply(F, NS);
    double MyEnergyAfter = (C * NS) * NS;
    cout << "Energy before smoothing = " << MyEnergyBefore << endl;
    cout << "Energy after smoothing  = " << MyEnergyAfter << endl;

    // break if coarse matrix is below specified tolerance
    if (A[level + 1].GetDomainSpace().GetNumGlobalElements() <= MaxCoarseSize) {
      ++level;
      break;
    }
  }

  MaxLevels = level;

  if (GetPrintLevel())
    ML_print_line("-", 80);

  // now run smoother on all levels but last
  for (int level = MaxLevels - 1 ; level > 0 ; --level)
    NS = P[level] * NS;

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

    SetPrintLevel(10);

    int MaxLevels = 3;
    int NumGlobalElements = (int)pow(2.0, (double)MaxLevels);

    double Pi = atan(1.0) * 4;
    double LambdaMin = 2.0 * (sin(1.0 / (NumGlobalElements + 1) * Pi / 2)) * 2 * .95;

    // define the space for fine level vectors and operators.
    Space FineSpace(NumGlobalElements);

    DistributedMatrix MatA(FineSpace, FineSpace);

    // assembel the matrix on processor 0 only
    if (GetMyPID() == 0) {
      for (int i = 0 ; i < NumGlobalElements ; ++i) {
        MatA.SetElement(i, i, 2.0); //  - LambdaMin;
        if (i)
          MatA.SetElement(i, i - 1, - 1.0);
        if (i != NumGlobalElements - 1)
          MatA.SetElement(i, i + 1, - 1.0);
      }
    }
    MatA.FillComplete();

    // wrap MatA as an Operator
    Operator A(FineSpace, FineSpace, &MatA, false);

    // information about the current null space
    int NumPDEEquations = 1;
    int NullSpaceDimension = 1;
    MultiVector NS(FineSpace, NullSpaceDimension);
    NS = 1.0; // constant null space
    
    Teuchos::ParameterList List;
    List.set("PDE equations", NumPDEEquations);
    List.set("aggregation: null space", NS);
    List.set("max levels", MaxLevels);
    List.set("smoother: type", "symmetric Gauss-Seidel");
    List.set("smoother: sweeps", 1);
    List.set("smoother: damping factor", 0.67);
    List.set("coarse: type", "Amesos-KLU");
 
    // ================================================== //
    // this is `phase1' ==> computation of a tentative    //
    // null space, supposing that we have no information  //
    // about any null space vector (that is, not even the //
    // constant).                                         //
    // ================================================== //
    
    MultiVector TNS;
    TNS = GetTentativeNullSpace(A, List, 1);
  
  
    // FIXME: to add `phase2'

    // FIXME: store the computed null space in ADB

    // create the multilevel preconditioner
    MultiLevelSA Prec(A, List);



    // test the solver
    MultiVector LHS(FineSpace);
    MultiVector RHS(FineSpace);

    LHS = 0.0;
    RHS.Random();

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
