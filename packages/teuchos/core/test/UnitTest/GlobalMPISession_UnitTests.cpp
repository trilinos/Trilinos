// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// @HEADER

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Array.hpp"

#ifdef HAVE_MPI
#  include "mpi.h"
#endif

#include "Teuchos_UnitTestHarness.hpp"


//
// Unit tests for GlobalMPISession
//
// NOTE: Becuase this class is used to implement the parallel reduction
// feature of the unit test harness, we can't use that feature here and we
// have to do the global reductions across processes ourselves.
//


namespace Teuchos {


void globalReduceSuccess(bool &success, FancyOStream &out)
{
#ifdef HAVE_MPI
  int globalSumSuccessInt = -1;
  int localSuccessInt = (success ? 0 : 1);
  MPI_Allreduce(&localSuccessInt, &globalSumSuccessInt, 1,
    MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
  TEST_EQUALITY_CONST(globalSumSuccessInt, 0);
#endif
}


TEUCHOS_UNIT_TEST( GlobalMPISession, basic ) {
#ifdef HAVE_MPI
  TEST_ASSERT(GlobalMPISession::mpiIsInitialized());
  int numProcs = -1;
  ECHO(::MPI_Comm_size(MPI_COMM_WORLD, &numProcs));
  TEST_EQUALITY(GlobalMPISession::getNProc(), numProcs);
  int procRank = -1;
  ECHO(::MPI_Comm_rank(MPI_COMM_WORLD, &procRank));
  TEST_EQUALITY(GlobalMPISession::getRank(), procRank);
#else // HAVE_MPI
  TEST_ASSERT(!GlobalMPISession::mpiIsInitialized());
  TEST_EQUALITY_CONST(GlobalMPISession::getNProc(), 1);
  TEST_EQUALITY_CONST(GlobalMPISession::getRank(), 0);
#endif // HAVE_MPI
  TEST_ASSERT(!GlobalMPISession::mpiIsFinalized());
  globalReduceSuccess(success, out);
}


TEUCHOS_UNIT_TEST( GlobalMPISession, barrier ) {
  out << "*** Just make sure the basic barrier does not hang or something.\n";
  ECHO(GlobalMPISession::barrier());
  globalReduceSuccess(success, out);
}


TEUCHOS_UNIT_TEST( GlobalMPISession, sum ) {
  ECHO(const int globalSum = GlobalMPISession::sum(GlobalMPISession::getRank()+1));
  ECHO(const int n = GlobalMPISession::getNProc());
  TEST_EQUALITY(globalSum, (n*(n+1))/2);
  globalReduceSuccess(success, out);
}


TEUCHOS_UNIT_TEST( GlobalMPISession, allGather )
{
  const int numProcs = GlobalMPISession::getNProc();
  const int procRank = GlobalMPISession::getRank();
  {
    Array<int> allInts;
    ECHO(allInts.resize(numProcs-1));
    TEST_THROW(GlobalMPISession::allGather(procRank+1, allInts()), std::out_of_range);
    ECHO(allInts.resize(numProcs+1));
    TEST_THROW(GlobalMPISession::allGather(procRank+1, allInts()), std::out_of_range);
  }
  {
    Array<int> allInts_expected(numProcs);
    for (int k = 0; k < numProcs; ++k) {
      allInts_expected[k] = k+1;
    }
    Array<int> allInts(numProcs);
    ECHO(GlobalMPISession::allGather(procRank+1, allInts()));
    TEST_EQUALITY(allInts, allInts_expected);
  }
}


} // namespace Teuchos



