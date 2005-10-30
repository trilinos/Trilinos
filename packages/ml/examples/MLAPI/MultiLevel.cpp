
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

#if defined(HAVE_ML_MLAPI) && defined(HAVE_ML_GALERI)
#include "MLAPI_Space.h"
#include "MLAPI_Operator.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_Gallery.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_MultiLevelSA.h"
#include "MLAPI_Krylov.h"

using namespace Teuchos;
using namespace MLAPI;

// ============== //
// example driver //
// ============== //

int main(int argc, char *argv[])
{
  
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  try {

    // Initialize the workspace and set the output level
    Init();

    int NumGlobalElements = 10000;

    // define the space for fine level vectors and operators.
    Space FineSpace(NumGlobalElements);

    // define the linear system matrix, solution and RHS
    Operator FineMatrix = Gallery("Tridiag", FineSpace);
    MultiVector LHS(FineSpace);
    MultiVector RHS(FineSpace);

    LHS = 0.0;
    RHS.Random();

    // set parameters for aggregation and smoothers
    // NOTE: only a limited subset of the parameters accepted by
    // class ML_Epetra::MultiLevelPreconditioner is supported
    // by MLAPI::MultiLevelSA
    
    Teuchos::ParameterList MLList;
    MLList.set("max levels",3);
    MLList.set("increasing or decreasing","increasing");
    MLList.set("aggregation: type", "Uncoupled");
    MLList.set("aggregation: damping factor", 0.0);
    MLList.set("smoother: type","symmetric Gauss-Seidel");
    MLList.set("smoother: sweeps",1);
    MLList.set("smoother: damping factor",1.0);
    MLList.set("coarse: max size",3);
    MLList.set("smoother: pre or post", "both");
    MLList.set("coarse: type","Amesos-KLU");

    MultiLevelSA Prec(FineMatrix, MLList);

    // solve with GMRES (through AztecOO)
    MLList.set("krylov: solver", "gmres");
    MLList.set("krylov: max iterations", 1550);
    MLList.set("krylov: tolerance", 1e-9);
    MLList.set("krylov: output", 16);
    Krylov(FineMatrix, LHS, RHS, Prec, MLList);

    cout << Prec;
  }
  catch (const int e) {
    cerr << "Caught exception, code = " << e << endl;
  }
  catch (...) {
    cerr << "Caught exception..." << endl;
  }

  Finalize();

#ifdef HAVE_MPI
  MPI_Finalize();
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

  puts("This MLAPI example requires the following configuration options:");
  puts("\t--enable-epetra");
  puts("\t--enable-teuchos");
  puts("\t--enable-ifpack");
  puts("\t--enable-amesos");
  puts("\t--enable-galeri");
  puts("Please check your configure line.");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(0);
}

#endif // #if defined(HAVE_ML_MLAPI)
