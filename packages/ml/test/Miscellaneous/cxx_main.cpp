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

/*
Current revision: $Revision$
Branch:           $Name$
Last modified:    $Date$
Modified by:      $Author$
*/


#include "ml_config.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "ml_include.h"

/*
   This test is meant as a catch-all for tests that don't fit into other
   categories.  I lifted code liberally from aztecoo/test/AztecOO/cxx_main.cpp.
*/

// == ========== prototypes =============================================
bool argument_is_present(const char* argument, int argc, char** argv);
int test_bug2863(Epetra_Comm& Comm, bool verbose);

// == ========== main program ===========================================

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  bool verbose = false;
  if (Comm.MyPID() == 0) {
    verbose = argument_is_present("-v", argc, argv);
  }

  if ( test_bug2863(Comm, verbose) ) {
    cout << "test_bug2863 FAILED."<<endl;
    exit(EXIT_FAILURE);
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  if (Comm.MyPID() == 0)
    cout << "TEST PASSED" << endl;

  exit(EXIT_SUCCESS);
  
} //end of main program

// == =========================================================== ==

bool argument_is_present(const char* argument,
                         int argc,
                         char** argv)
{
  if (argument == NULL || argc < 1) return(false);

  for(int i=0; i<argc; ++i) {
    if (strcmp(argument, argv[i]) == 0) {
      return(true);
    }
  }
  return(false);
}

// == =========================================================== ==

int test_bug2863(Epetra_Comm& Comm, bool verbose)
{
//This function tests the ML_random1() function in ML.
//The implementation of the Park and Miller random number
//generator was incorrect and resulted in an overflow condition.
//This is *not* a complete test of ML's RNG, but just tests a particular
//seed that was causing a failure.
//
//A more robust check is to compile ML with gcc -ftrapv and run
//the complete test suite.

#define numSeeds 3
  int seeds[numSeeds] = {1, 127773, 2147483646};

  if (verbose) printf("test_bug2863: Testing Park and Miller serial RNG.\n");

  for (int i=0; i<numSeeds; i++) {
    int j = seeds[i];
    double rand_num = ML_srandom1(&j);
    if (verbose)
      printf("test_bug2863: seed = %11d, rand_num = %e (should be in (0,1))\n",
             seeds[i],rand_num);
    if ( (rand_num > 1) || (rand_num < 0) )
      return 1;    // rand_num should be in [0,1]
  }

  return 0;
}
