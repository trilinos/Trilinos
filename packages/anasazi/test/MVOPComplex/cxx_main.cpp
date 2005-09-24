//@HEADER
// ************************************************************************
// 
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
// ************************************************************************
//@HEADER
//
//  This test instantiates the Anasazi classes using a complex scalar type
//  and checks functionality.
//

#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziEpetraAdapter.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_COMPLEX
typedef std::complex<double> ScalarType;
#elif HAVE_COMPLEX_H
typedef ::complex<double> ScalarType;
#else
typedef double ScalarType
#endif

#include "MyMultiVec.hpp"
#include "MyOperator.hpp"

int main(int argc, char *argv[])
{
  int i, ierr, gerr;
  gerr = 0;

#ifdef HAVE_MPI
  // Initialize MPI and setup an Epetra communicator
  MPI_Init(&argc,&argv);
#endif

  // PID info
  int MyPID = 0;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &MyPID);
#endif
  bool verbose = 0;

  if (argc>1) {
    if (argv[1][0]=='-' && argv[1][1]=='v') {
      verbose = true;
    }
  }

  // Create a sort manager to sort for the eigenvalue with largest real part.
  Anasazi::BasicSort<ScalarType, MyMultiVec, MyOperator> mySort("LR");


#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (gerr) {
    if (verbose && MyPID==0)
      cout << "End Result: TEST FAILED" << endl;	
    return -1;
  }
  //
  // Default return value
  //
  if (verbose && MyPID==0)
    cout << "End Result: TEST PASSED" << endl;
  return 0;

}
