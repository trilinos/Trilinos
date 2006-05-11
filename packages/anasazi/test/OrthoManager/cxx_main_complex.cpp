// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2004) Sandia Corporation
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
// ***********************************************************************
// @HEADER
//
// This test is for the BasicOrthoManager using a complex operator.

#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicOrthoManager.hpp"
#include "AnasaziModalSolverUtils.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

// use Thyra with a custom, non-identity operator
// #include<thyra stuff>

using namespace Teuchos;
using namespace Anasazi;

#ifdef HAVE_COMPLEX
  typedef std::complex<double> ST;
#elif HAVE_COMPLEX_H
  typedef ::complex<double> ST;
#else
  typedef double ST;
#endif
//typedef Epetra_MultiVector MV;
//typedef Epetra_Operator OP;
//typedef MultiVecTraits<ST,MV>    MVT;
//typedef OperatorTraits<ST,MV,OP> OPT;
typedef ScalarTraits<ST>         SCT;
typedef SCT::magnitudeType       MT;

const ST ONE  = SCT::one();
const MT ZERO = SCT::magnitude(SCT::zero());

int main(int argc, char *argv[]) 
{
  int MyPID;
#ifdef HAVE_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &MyPID);
#endif
  
  int numFailed = 0;
  bool verbose = true;
  if (argc>1) {
    if (argv[1][0]=='-' && argv[1][1]=='v') {
      verbose = true;
    }
  }

  if (verbose && MyPID == 0) {
    cout << Anasazi_Version() << endl << endl;
  }
  



#ifndef HAVE_COMPLEX && HAVE_COMPLEX_H
  // no complex. quit with failure.
  if (verbose && MyPID == 0) {
    cout << "Not compiled with complex support." << endl;
    cout << "End Result: TEST FAILED" << endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }
#endif




  cout << endl;

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  cout << "End Result: TEST PASSED" << endl;
  return 0;

}	
