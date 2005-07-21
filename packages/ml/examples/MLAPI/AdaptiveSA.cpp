
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
#include "MLAPI_Space.h"
#include "MLAPI_Operator.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_Gallery.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_MultiLevelAdaptiveSA.h"
#include "MLAPI_DistributedMatrix.h"
#include "MLAPI_Krylov.h"


using namespace Teuchos;
using namespace MLAPI;

// =========== //
// main driver //
// =========== //

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  try {

    // Initialize the workspace and set the output level
    Init();

    int NX = 1000;

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

    int NumPDEEqns = 2;
    int MaxLevels = 10;

    Teuchos::ParameterList List;
    List.set("additional candidates", 2);
    List.set("use default null space", true);
    List.set("krylov: type", "cg");

    MultiLevelAdaptiveSA Prec(A, List, NumPDEEqns, MaxLevels);

    // =============================================================== //
    // setup the hierarchy:                                            //
    // - `UseDefaultNullSpace' toggles the use of default candidates.  //
    // - AdditionalCandidates = 2' means to compute two additionals.   //
    // - the final null space dimension is 3.                          //
    // =============================================================== //
    
    bool UseDefaultNullSpace = true;
    int AdditionalCandidates = 1;
    Prec.AdaptCompute(UseDefaultNullSpace, AdditionalCandidates);

    MultiVector LHS(A.GetDomainSpace());
    MultiVector RHS(A.GetRangeSpace());

    LHS.Random();
    RHS = 0.0;

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
