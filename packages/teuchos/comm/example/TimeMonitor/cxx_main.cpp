// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_Version.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

using namespace Teuchos;

// Global Timers
RCP<Time> CompTime = TimeMonitor::getNewCounter("Computational Time");
RCP<Time> FactTime = TimeMonitor::getNewCounter("Factorial Time");

// Quadratic function declaration.
double quadFunc( double x );

// Factorial function declaration.
double factFunc( int x );

int main(int argc, char* argv[])
{
  std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;

  int i;
  double x;

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  // Apply the quadratic function.
  for( i=-100; i<100; i++ ) {
    x = quadFunc( (double) i );
    (void)x; // Not used!
  }

  // Apply the factorial function.
  for( i=0; i<100; i++ ) {
    x = factFunc( i );
    (void)x; // Not used!
  }

  // Get a summary from the time monitor.
  TimeMonitor::summarize();

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}

/* Evaluate a quadratic function at point x */
double quadFunc( double x )
{
  // Construct a local time monitor, this starts the CompTime timer and will stop when leaving scope.
  Teuchos::TimeMonitor LocalTimer(*CompTime);

  // Evaluate the quadratic function.
  return ( x*x - 1.0 );
}

/* Compute the factorial of x */
double factFunc( int x )
{
  // Construct a local time monitor, this starts the FactTime timer and will stop when leaving scope.
  Teuchos::TimeMonitor LocalTimer(*FactTime);

  // Special returns for specific cases.
  if( x == 0 ) return 0.0;
  if( x == 1 ) return 1.0;

  // Evaluate the factorial function.
  return ( (double) x * factFunc(x-1) );
}
