// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_GlobalMPISession_HPP
#define ROL_GlobalMPISession_HPP

// Native, Teuchos-free real-MPI ROL::GlobalMPISession. Drop-in replacement for the
// serial stub at src/compatibility/simple/mpi/ROL_GlobalMPISession.hpp. It reuses the
// SAME include guard on purpose: exactly one of compatibility/{real,simple}/mpi may be
// on the -I path (the root CMakeLists.txt selects it via ENABLE_MPI). Whichever appears
// first wins; the other is guard-suppressed.

#include <ostream>
#include <mpi.h>

namespace ROL {

struct GlobalMPISession {

  // Accepts the existing call sites:
  //   GlobalMPISession(&argc,&argv);      and   GlobalMPISession(&argc,&argv,0);
  // The 3rd arg mirrors Teuchos' optional rank-0 ostream*; accepted and ignored.
  GlobalMPISession( int* argc = nullptr,
                    char*** argv = nullptr,
                    std::ostream* /*out*/ = nullptr )
  {
    int alreadyInit = 0;
    MPI_Initialized( &alreadyInit );
    if( !alreadyInit ) { MPI_Init( argc, argv ); ownInit_ = true; }
  }

  ~GlobalMPISession()
  {
    if( ownInit_ ) {
      int finalized = 0;
      MPI_Finalized( &finalized );
      if( !finalized ) MPI_Finalize();
    }
  }

  GlobalMPISession( const GlobalMPISession& ) = delete;
  GlobalMPISession& operator=( const GlobalMPISession& ) = delete;

  static void abort() { MPI_Abort( MPI_COMM_WORLD, -1 ); }

  static void barrier() { if( initialized() ) MPI_Barrier( MPI_COMM_WORLD ); }

  static int getNProc() {
    if( !initialized() ) return 1;
    int n = 1; MPI_Comm_size( MPI_COMM_WORLD, &n ); return n;
  }

  // TRUE rank (serial stub wrongly returned the constant 1).
  static int getRank() {
    if( !initialized() ) return 0;
    int r = 0; MPI_Comm_rank( MPI_COMM_WORLD, &r ); return r;
  }

  static int sum( int localVal ) {
    if( !initialized() ) return localVal;
    int g = 0; MPI_Allreduce( &localVal, &g, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
    return g;
  }

private:
  static bool initialized() {
    int flag = 0; MPI_Initialized( &flag );
    if( !flag ) return false;
    int fin = 0; MPI_Finalized( &fin );
    return !fin;
  }
  bool ownInit_ = false;
};

} // namespace ROL

#endif // ROL_GlobalMPISession_HPP
