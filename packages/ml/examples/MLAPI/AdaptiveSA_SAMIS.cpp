
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
#include "MLAPI_SAMIS.h"

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

  if (argc != 2) {
    fprintf(stderr, "Usage: `%s InputFile'\n", argv[0]);
    fprintf(stderr, "An example of input file is reported\n");
    fprintf(stderr, "in the source of this example\n");
    exit(EXIT_SUCCESS);
  }

  string InputFile = argv[1];

  // Initialize the workspace and set the output level
  Init();

  try {

    int         NumPDEEqns;
    Operator    A;

    ReadSAMISMatrix("mtx.dat", A, NumPDEEqns);

    Teuchos::ParameterList List = ReadParameterList(InputFile.c_str());
    int  MaxLevels            = List.get("max levels", 10);
    int  AdditionalCandidates = List.get("additional candidates", 2);
    int  limKer               = List.get("limit kernel", -1);

    if (AdditionalCandidates == 0 && limKer == 0)
      limKer = -1;

    // create multilevel preconditioner, do not compute hierarchy
    MultiLevelAdaptiveSA Prec(A, List, NumPDEEqns);

    if (limKer) {
      MultiVector NS;
      ReadSAMISKernel("ker.dat", NS);
      Prec.SetNullSpace(NS);
      Prec.AdaptCompute(true, AdditionalCandidates);
    }
    else {
      Prec.AdaptCompute(false, AdditionalCandidates);
    }

    MultiVector LHS(A.GetDomainSpace());
    MultiVector RHS(A.GetRangeSpace());

    LHS.Random();
    RHS = 0.0;

    List.set("krylov: type", "cg_condnum");
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
