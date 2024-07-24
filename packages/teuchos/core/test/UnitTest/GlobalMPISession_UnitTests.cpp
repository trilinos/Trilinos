// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Array.hpp"

#ifdef HAVE_MPI
#  include "mpi.h"
#endif

#include "Teuchos_UnitTestHarness.hpp"

#ifdef HAVE_TEUCHOSCORE_KOKKOS
#  include <string>
// NOTE (mfh 18 Apr 2016) KokkosCore requires C++11, so we may include
// type_traits here.  Please do not include it unconditionally,
// because TeuchosCore does NOT require C++11.
#  include <type_traits>
#endif // HAVE_TEUCHOSCORE_KOKKOS

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


#ifdef HAVE_TEUCHOSCORE_KOKKOS
// Check whether GlobalMPISession's constructor (called in the unit
// test harness, not in this file) correctly saved a copy of the
// command-line arguments.  Command-line arguments _should_ be
// propagated to all MPI processes.
//
//
TEUCHOS_UNIT_TEST( GlobalMPISession, getArgv )
{
  auto argvCopy = GlobalMPISession::getArgv ();
  const auto numArgs = argvCopy.size ();
  static_assert (std::is_integral<std::decay<decltype (numArgs)>::type>::value,
		 "The return type of getArgv has a size() method, "
		 "but it does not return an integer.");
  if (numArgs > 0) {
    // This tests that argvCopy has an operator[](int) const.
    const auto arg0 = argvCopy[0];
    static_assert (std::is_convertible<std::decay<decltype (arg0)>::type,
		     std::string>::value,
		   "The return type of getArgv must have an operator[](int) "
		   "that returns a type convertible to std::string.");
  }
}
#endif // HAVE_TEUCHOSCORE_KOKKOS


} // namespace Teuchos



