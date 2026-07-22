// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_Version.hpp"
#include "Teuchos_as.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

using std::string;
using Teuchos::TimeMonitor;
using Teuchos::Time;
using Teuchos::RCP;
using Teuchos::ScalarTraits;
using Teuchos::as;

/* Test of Teuchos timing classes */


/* create timers for several functions */
static Time& sqrtTimer() {static RCP<Time> t = TimeMonitor::getNewTimer("square roots"); return *t;}

static Time& factTimer() {static RCP<Time> t = TimeMonitor::getNewTimer("factorials"); return *t;}

static Time& exceptTimer() {static RCP<Time> t = TimeMonitor::getNewTimer("func with std::exception"); return *t;}

static Time& localTimer() {static RCP<Time> t = TimeMonitor::getNewTimer("a function that is not called on all procs"); return *t;}

static Time& anotherTimer() {static RCP<Time> t = TimeMonitor::getNewTimer("another func"); return *t;}

static Time& yetAnotherTimer() {static RCP<Time> t = TimeMonitor::getNewTimer("yet another func"); return *t;}

static Time& yetOneMoreTimer() {static RCP<Time> t = TimeMonitor::getNewTimer("yet one more func"); return *t;}


// Prototypes of the functions we will time below.
double sqrtFunc();
double factFunc(int x);
double exceptFunc();
double localFunc();
double anotherFunc();
double yetAnotherFunc();
double yetOneMoreFunc();


int main(int argc, char* argv[])
{
  bool verbose = 0;
  int procRank = 0;
  int FailedTests = 1; // This will be set to 0, if the std::exception is caught!

#ifdef HAVE_MPI
  /* initialize MPI if we are running in parallel */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &procRank );
#endif

  // Check for verbose flag.
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  if (verbose && procRank==0)
    std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;

  try {
    // time a simple function
    for (int i=0; i<100; i++)
    {
      double x = 0.0;
      x = sqrtFunc ();
      (void) x; // forestall "variable set but not used" compiler warning
    }

    // time a reentrant function
    for (int i = 0; i < 100; ++i) {
      factFunc (100);
    }

    /* time a couple of silly functions */
    for (int i = 0; i < 100; ++i) {
      anotherFunc ();
      yetAnotherFunc ();
      yetOneMoreFunc ();
    }

    /* Time a function that will be called only on the root proc. This
     * checks that the TimeMonitor will work properly when different
     * processors have different sets of timers. */
    if (procRank == 0) {
      for (int i = 0; i < 100; ++i) {
        double x = 0.0;
        x = localFunc ();
	(void) x; // forestall "variable set but not used" compiler warning
      }
    }

    // time a function that throws an exception
    for (int i = 0; i < 100; ++i) {
      double x = 0.0;
      x = exceptFunc ();
      (void) x; // forestall "variable set but not used" compiler warning
    }
  }
  catch (std::exception& e) {
    // This _should_ only catch the one exception thrown above by the
    // first call to exceptFunc().
    if (verbose && procRank==0) {
      std::cerr << "Caught std::exception [expected]:  " << e.what() << std::endl;
    }
    // Return 0 since we caught the std::exception
    FailedTests = 0;
  }

  // Summarize timings. This must be done before finalizing MPI.
  TimeMonitor::format ().setRowsBetweenLines (3);
  if (verbose) {
    TimeMonitor::summarize ();
  }

#ifdef HAVE_MPI
  /* clean up MPI if we are running in parallel*/
  MPI_Finalize ();
#endif // HAVE_MPI

  if (FailedTests != 0) {
    std::cout << "End Result: TEST FAILED" << std::endl;
    return 1;
  }

  std::cout << "End Result: TEST PASSED" << std::endl;
  return FailedTests;
}


/* sum std::sqrt(x), x=[0, 10000). */
double sqrtFunc()
{
  /* construct a time monitor. This starts the timer. It will stop when leaving scope */
  TimeMonitor timer(sqrtTimer());

  double sum = 0.0;

  for (int i=0; i<10000; i++)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(ScalarTraits<double>::squareroot(as<double>(i)) > 1000.0, std::runtime_error,
      "throw an std::exception");
    sum += ScalarTraits<double>::squareroot(as<double>(i));
  }

  return sum;
}


/* compute log(factorial(x)) */
double factFunc(int x)
{
  /* construct a time monitor. This starts the timer. It will stop when leaving scope */
  TimeMonitor timer(factTimer());

  if (x==0) return 0;
  if (x==1) return 1;
  return std::log(as<double>(x))  + factFunc(x-1);
}


/* sum std::sqrt(x), x=[0, 10000). */
double exceptFunc()
{
  /* construct a time monitor. This starts the timer. It will stop when leaving scope */
  TimeMonitor timer(exceptTimer());

  double sum = 0.0;
  for (int i=0; i<10000; i++)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      ScalarTraits<double>::squareroot(as<double>(i)) > 60.0, std::runtime_error,
      "throw an std::exception");
    sum += ScalarTraits<double>::squareroot(as<double>(i));
  }
  return sum;
}


/* sum x, x=[0, 10000). */
double localFunc()
{
  /* construct a time monitor. This starts the timer. It will stop when leaving scope */
  TimeMonitor timer(localTimer());

  double sum = 0.0;

  for (int i=0; i<10000; i++)
  {
    sum += i;
  }

  return sum;
}


/* sum x^2, x=[0, 10000). */
double anotherFunc()
{
  /* construct a time monitor. This starts the timer. It will stop when leaving scope */
  TimeMonitor timer(anotherTimer());

  double sum = 0.0;

  for (int i=0; i<10000; i++)
  {
    sum += i*i;
  }

  return sum;
}


/* sum x^3, x=[0, 10000). */
double yetAnotherFunc()
{
  /* construct a time monitor. This starts the timer. It will stop when leaving scope */
  TimeMonitor timer(yetAnotherTimer());

  double sum = 0.0;

  for (int i=0; i<10000; i++)
  {
    sum += i*i*i;
  }

  return sum;
}


/* sum x+1, x=[0, 10000). */
double yetOneMoreFunc()
{
  /* construct a time monitor. This starts the timer. It will stop when leaving scope */
  TimeMonitor timer(yetOneMoreTimer());

  double sum = 0.0;

  for (int i=0; i<10000; i++)
  {
    sum += i+1;
  }

  return sum;
}
