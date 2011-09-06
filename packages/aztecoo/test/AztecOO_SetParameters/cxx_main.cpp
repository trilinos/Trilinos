//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

// AztecOO_SetParameters Test routine
#include <AztecOO.h>

#ifdef HAVE_AZTECOO_TEUCHOS
#include <Teuchos_ParameterList.hpp>
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif

int main(int argc, char* argv[]) {
  bool verbose1 = false; // user's command-line argument
  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose1 = true;

  int err;
  int returnierr = 0;

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  if (verbose1) {
  }

#ifdef HAVE_AZTECOO_TEUCHOS
  AztecOO azoo;

  Teuchos::ParameterList paramlist;
  paramlist.set("AZ_solver", AZ_cg);
  paramlist.set("max_ITER", 200);
  paramlist.set("az_Tol", 1.0e-09);
  paramlist.set("precond", AZ_Jacobi);

  paramlist.set("AZ_kspace", 2.5);//bad type

  if (verbose1==true) {
    cout << "Test: parameter 'AZ_kspace' given bad type (double), warning should"
         << "be printed to cerr"<<endl;
    err = azoo.SetParameters(paramlist, verbose1);
  }
  else {
    err = azoo.SetParameters(paramlist);
  }
  if (err != 0) {
    if (verbose1) {
      cerr << "err " << err << " returned from AztecOO::SetParameters"<<endl;
    }
    return(-1);
  }

  const int* options = azoo.GetAllAztecOptions();
  if (options[AZ_solver] != AZ_cg) {
    if (verbose1) {
      cerr << "SetParameters test failed to correctly set AZ_solver."<<endl;
    }
    return(-1);
  }

  if (options[AZ_max_iter] != 200) {
    if (verbose1) {
      cerr << "SetParameters test failed to correctly set AZ_max_iter."<<endl;
    }
    return(-1);
  }

  if (options[AZ_precond] != AZ_Jacobi) {
    if (verbose1) {
      cerr << "SetParameters test failed to correctly set AZ_precond."<<endl;
    }
    return(-1);
  }

  const double* params = azoo.GetAllAztecParams();
  if (params[AZ_tol] != 1.0e-09) {
    if (verbose1) {
      cerr << "SetParameters test failed to correctly set AZ_tol."<<endl;
      cout << "********* Test failed **********" << endl;
    }
    return(-1);
  }
#endif
// Typically we only produce output when the verbose flag is used.
// Aztecoo has a test harness script that relies on "passed" appearing
// in the output for the script, so for now we will produce the output
// below whenever the test passes.
//  if (verbose1==true) {
    cout << "********* Test passed **********" << endl;
//  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(returnierr);
}

