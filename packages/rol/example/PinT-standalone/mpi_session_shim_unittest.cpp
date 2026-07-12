// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Unit tests for the standalone-ROL additions that enable real-MPI PinT:
//   (1) src/compatibility/real/mpi/ROL_GlobalMPISession.hpp   (real MPI session)
//   (2) src/compatibility/teuchos-lite/Teuchos_Time.hpp        (timer + symbol shims)
// Exercises ONLY the new code in isolation (no PinTConstraint). Run under mpirun
// (np = 1 or 2). Rank 0 prints "End Result: TEST PASSED" iff every check on every
// rank passed; the process also exits nonzero on any failure (for ctest).

#include <mpi.h>
#include <iostream>
#include <cstddef>

#include "ROL_GlobalMPISession.hpp"   // -> real/mpi under ENABLE_MPI
#include "ROL_Ptr.hpp"                // ROL::Ptr, makePtr, nullPtr
#include "Teuchos_Time.hpp"           // Teuchos::as / null, ROL::is_null, StackedTimer

namespace {

int g_fail = 0;
int g_rank = 0;

// Minimal check macro: counts failures, names the offending rank/line.
#define CHECK( cond, msg )                                                   \
  do {                                                                       \
    if( !(cond) ) {                                                          \
      ++g_fail;                                                              \
      std::cerr << "[rank " << g_rank << "] FAIL: " << (msg)                 \
                << "   (" #cond ", line " << __LINE__ << ")\n";              \
    }                                                                        \
  } while( 0 )

// ---- (1) real-MPI ROL::GlobalMPISession -------------------------------------
void test_session()
{
  int mpiRank = -1, mpiSize = -1;
  MPI_Comm_rank( MPI_COMM_WORLD, &mpiRank );
  MPI_Comm_size( MPI_COMM_WORLD, &mpiSize );

  int initFlag = 0;
  MPI_Initialized( &initFlag );
  CHECK( initFlag != 0, "MPI is initialized after the session constructor" );

  CHECK( ROL::GlobalMPISession::getRank()  == mpiRank, "getRank() matches MPI_Comm_rank" );
  CHECK( ROL::GlobalMPISession::getNProc() == mpiSize, "getNProc() matches MPI_Comm_size" );
  CHECK( ROL::GlobalMPISession::getRank() >= 0 &&
         ROL::GlobalMPISession::getRank() < ROL::GlobalMPISession::getNProc(),
         "rank is in [0, nproc)" );

  // exactly one process is rank 0 (the distinctness property we care about)
  const int amRank0 = ( ROL::GlobalMPISession::getRank() == 0 ) ? 1 : 0;
  CHECK( ROL::GlobalMPISession::sum( amRank0 ) == 1, "exactly one process is rank 0" );

  // collective reductions
  CHECK( ROL::GlobalMPISession::sum( 1 ) == mpiSize, "sum(1) == nproc (Allreduce)" );
  CHECK( ROL::GlobalMPISession::sum( mpiRank ) == mpiSize * (mpiSize - 1) / 2,
         "sum(rank) == n(n-1)/2" );

  ROL::GlobalMPISession::barrier();   // must not hang or crash
}

// ---- (2) teuchos-lite shims -------------------------------------------------
void test_shims()
{
  // Teuchos::as<Out>(in) == static_cast<Out>(in)
  CHECK( Teuchos::as<int>( 3.9 )       == 3,              "as<int>(3.9) truncates to 3" );
  CHECK( Teuchos::as<int>( -2.7 )      == -2,             "as<int>(-2.7) truncates toward zero" );
  CHECK( Teuchos::as<double>( 5 )      == 5.0,            "as<double>(5) == 5.0" );
  CHECK( Teuchos::as<std::size_t>( 7 ) == std::size_t(7), "as<size_t>(7) == 7" );

  // Teuchos::null and ROL::is_null on ROL::Ptr (= std::shared_ptr)
  ROL::Ptr<int> pnull = ROL::nullPtr;
  ROL::Ptr<int> plive = ROL::makePtr<int>( 42 );

  CHECK( pnull == Teuchos::null, "null Ptr == Teuchos::null" );
  CHECK( Teuchos::null == pnull, "Teuchos::null == null Ptr (reversed)" );
  CHECK( plive != Teuchos::null, "live Ptr != Teuchos::null" );
  CHECK( Teuchos::null != plive, "Teuchos::null != live Ptr (reversed)" );
  CHECK( ROL::is_null( pnull ),  "is_null(null) is true" );
  CHECK( !ROL::is_null( plive ), "is_null(live) is false" );
  CHECK( *plive == 42,           "live Ptr dereferences to its value" );

  // StackedTimer shim: non-null, stable static instance, start/stop are safe no-ops
  ROL::Ptr<Teuchos::StackedTimer> t1 = Teuchos::TimeMonitor::getStackedTimer();
  ROL::Ptr<Teuchos::StackedTimer> t2 = Teuchos::TimeMonitor::getStackedTimer();
  CHECK( t1 != Teuchos::null,  "getStackedTimer() returns non-null" );
  CHECK( t1.get() == t2.get(), "getStackedTimer() returns the same instance" );
  t1->start( "unit" );
  t1->stop( "unit" );          // compiles and does not crash
}

} // namespace

int main( int argc, char** argv )
{
  ROL::GlobalMPISession session( &argc, &argv );
  g_rank = ROL::GlobalMPISession::getRank();

  test_session();
  test_shims();

  // Agree on a single verdict across all ranks (also exercises sum()).
  const int globalFail = ROL::GlobalMPISession::sum( g_fail );

  if( g_rank == 0 ) {
    if( globalFail == 0 ) { std::cout << "End Result: TEST PASSED\n"; }
    else { std::cout << "End Result: TEST FAILED (" << globalFail << " checks failed)\n"; }
  }
  return ( globalFail == 0 ) ? 0 : 1;
}
