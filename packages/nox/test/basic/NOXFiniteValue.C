//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
  double infinity = 1.0/0.0;

  // Test return codes: 0 = finite, -1 = NaN, -2 = Inf
  int result = fv_test.finiteNumberTest(finite);
  
  if (result != 0)
    status = 1;  // Nonzero is failure

  if ( verbose && (myPID == 0) ) {
    if (result == 0)
      cout << "\nFinite value = " << finite 
	   << ", test correctly identified value as finite!" << endl;
    else
      cout << "Finite value = " << finite 
	   << ", test failed to identify value as finite!" << endl;
  }

  // Test return codes: 0 = finite, -1 = NaN, -2 = Inf
  result = fv_test.finiteNumberTest(nan);

  if (result != -1)
    status = 1;  // Nonzero is failure
  
  if ( verbose && (myPID == 0) ) {
    if (result == -1)
      cout << "NaN value = " << nan 
	   << ", test correctly identified value as nan!" << endl;
    else
      cout << "NaN value = " << nan 
	   << ", test failed to identify value as nan!" << endl;
  }

  // Test return codes: 0 = finite, -1 = NaN, -2 = Inf
  result = fv_test.finiteNumberTest(infinity);

  if (result != -2)
    status = 1;  // Nonzero is failure

  if ( verbose && (myPID == 0) ) {
    if (result == -2)
      cout << "Inf value = " << infinity
	   << ", test correctly identified value as inf!" << endl;
    else
      cout << "Inf value = " << infinity 
	   << ", test failed to identify value as inf!" << endl;
  }



#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (myPID == 0) {
    if (status == 0) 
      std::cout << "\nTest passed!" << endl;
    else
      std::cout << "\nTest Failed!" << endl;
  }

  // Final return value (0 = successfull, non-zero = failure)
  return status;
}
