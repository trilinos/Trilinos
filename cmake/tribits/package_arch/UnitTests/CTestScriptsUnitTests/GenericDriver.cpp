// @HEADER
// ************************************************************************
//
//            TriBITS: Tribial Build, Integrate, and Test System
//                    Copyright 2013 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
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
