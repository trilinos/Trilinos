// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
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


#include <iostream>
#include <cstdlib>
#include <string>

#ifdef HAVE_MPI
#  include <mpi.h>
#endif


int main( int argc, char* argv[] )
{

  bool success = true;
  int procRank = 0;

#ifdef HAVE_MPI
  MPI_Init (&argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &procRank );
#endif

  std::string printMsg = "";
  if (argc > 1)
    printMsg = argv[1];
  int returnVal = 0;
  if (argc > 2)
    returnVal = std::atoi(argv[2]);

  if (printMsg.length() && procRank==0)
    std::cout << "\n" << printMsg << "\n";
 
  if (procRank==0) {
    if (success)
      std::cout << "\nEnd Result: TEST PASSED" << std::endl;
    else
      std::cout << "\nEnd Result: TEST FAILED" << std::endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
 
  return returnVal;
 
}
