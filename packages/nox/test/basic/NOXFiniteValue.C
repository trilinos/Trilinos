//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

// Tests NOX's finite number test. 

#include "NOX.H"

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "NOX_StatusTest_FiniteValue.H"
#include "Teuchos_ScalarTraits.hpp"

using namespace std;

int main(int argc, char *argv[])
{

  int status = 0;

  int myPID = 0;

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myPID);
#endif

  // Check verbosity level
  bool verbose = false;
  if (argc > 1)
    if (argv[1][0]=='-' && argv[1][1]=='v')
      verbose = true;

  
  NOX::StatusTest::FiniteValue fv_test;

  double finite = 1.0e-100;
  double nan = sqrt(-1.0);
  double infinity = 1.0 / Teuchos::ScalarTraits<double>::zero();

  // Test return codes: 0 = finite, -1 = NaN, -2 = Inf
  int result = fv_test.finiteNumberTest(finite);
  
  if (result != 0)
    status = 1;  // Nonzero is failure

  if ( verbose && (myPID == 0) ) {
    if (result == 0)
      std::cout << "\nFinite value = " << finite 
	   << ", test correctly identified value as finite!" << std::endl;
    else
      std::cout << "Finite value = " << finite 
	   << ", test failed to identify value as finite!" << std::endl;
  }

  // Test return codes: 0 = finite, -1 = NaN, -2 = Inf
  result = fv_test.finiteNumberTest(nan);

  if (result != -1)
    status = 1;  // Nonzero is failure
  
  if ( verbose && (myPID == 0) ) {
    if (result == -1)
      std::cout << "NaN value = " << nan 
	   << ", test correctly identified value as nan!" << std::endl;
    else
      std::cout << "NaN value = " << nan 
	   << ", test failed to identify value as nan!" << std::endl;
  }

  // Test return codes: 0 = finite, -1 = NaN, -2 = Inf
  result = fv_test.finiteNumberTest(infinity);

  if (result != -2)
    status = 1;  // Nonzero is failure

  if ( verbose && (myPID == 0) ) {
    if (result == -2)
      std::cout << "Inf value = " << infinity
	   << ", test correctly identified value as inf!" << std::endl;
    else
      std::cout << "Inf value = " << infinity 
	   << ", test failed to identify value as inf!" << std::endl;
  }



#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (myPID == 0) {
    if (status == 0) 
      std::cout << "\nTest passed!" << std::endl;
    else
      std::cout << "\nTest Failed!" << std::endl;
  }

  // Final return value (0 = successfull, non-zero = failure)
  return status;
}
