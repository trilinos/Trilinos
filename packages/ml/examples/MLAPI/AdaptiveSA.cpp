
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
//
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

    // assemble the matrix on processor 0 only
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
    MultiVector Scale(FineSpace);
    Scale.Random();
    Scale = Scale + 1.0001;

    Operator S = GetDiagonal(Scale);
    Operator A = GetRAP(S, NonScaledA, S);
#else
    //Operator A = NonScaledA;
#endif

    Teuchos::ParameterList List;
    List.set("smoother: type", "symmetric Gauss-Seidel");
    List.set("smoother: sweeps", 10);
    List.set("smoother: damping factor", 1.0);
    List.set("coarse: type", "Amesos-KLU");
    List.set("coarse: max size", 32);
    List.set("adapt: max reduction", 0.1);
    List.get("adapt: iters fine", 15);
    List.get("adapt: iters coarse", 5);

    int NumPDEEqns = 2;
    int MaxLevels  = 10;
    MultiLevelAdaptiveSA Prec(A, List, NumPDEEqns, MaxLevels);

    // =============================================================== //
    // setup the hierarchy:                                            //
    // - `false' means that no default null space will be considered   //
    // - AdditionalCandidates = 2' means to compute two additional guy //
    // - the final null space dimension is 3.                          //
    // =============================================================== //
    
    int AdditionalCandidates = 2;
    Prec.AdaptCompute(false, AdditionalCandidates);

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

  puts("The ML API requires the following configuration options:");
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
