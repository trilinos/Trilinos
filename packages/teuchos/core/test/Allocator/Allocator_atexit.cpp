// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_Details_Allocator.hpp>
#include <sstream>
#include <vector>

#ifdef HAVE_TEUCHOS_MPI
#  include "mpi.h"
#endif // HAVE_TEUCHOS_MPI

// This test "passes" if the atexit() hook runs and prints the right
// thing.  _exit() skips over any atexit() hooks, so it is perfect for
// making the test fail if something goes wrong before main() is done.
// However, _exit() might not exist.  See the following man page to
// learn how to test whether _exit() exists.
//
//   http://man7.org/linux/man-pages/man2/_exit.2.html

namespace { // (anonymous)
  //! Whether the calling process should print at exit.
  bool iPrint_;

  //! String prefix for printing at exit.
  std::string prefix_;

  // Demonstration of an atexit() hook that uses
  // Teuchos::Details::AllocationLogger to report maximum and current
  // memory usage at exit from main().
  void allocationLoggerHook () {
    if (iPrint_) {
      using Teuchos::Details::AllocationLogger;
      using std::cout;
      using std::endl;

      // If printing on multiple MPI processes, print to a string
      // first, to hinder interleaving of output from different
      // processes.
      std::ostringstream os;
      os << prefix_ << "Teuchos allocation tracking: "
         << "Maximum allocation (B): "
         << AllocationLogger::maxAllocInBytes ()
         << ", Current allocation (B): "
         << AllocationLogger::curAllocInBytes () << endl;
      cout << os.str ();
    }
  }
} // namespace (anonymous)

#ifdef HAVE_TEUCHOS_MPI
int main (int argc, char* argv[])
#else
int main (int, char*[])
#endif // HAVE_TEUCHOS_MPI
{
  typedef std::vector<float, Teuchos::Details::Allocator<float>  > float_vec_type;
  typedef std::vector<long, Teuchos::Details::Allocator<long>  > long_vec_type;

#ifdef HAVE_TEUCHOS_MPI
  (void) MPI_Init (&argc, &argv);
  int myRank = 0;
  (void) MPI_Comm_rank (MPI_COMM_WORLD, &myRank);
  // Only let Process 0 print.  This ensures that the atexit() hook
  // won't print multiple times, with the possibility of interleaved
  // output and a test that fails unjustifiedly (CTest determines
  // passing by searching the test output for a string).
  iPrint_ = (myRank == 0);
#else
  iPrint_ = true;
#endif // HAVE_TEUCHOS_MPI
  prefix_ = std::string ("Proc 0: ");

  // Register the atexit() "hook" (the function to be called when
  // main() exists).
  //
  // It's possible for atexit() to fail.  In that case, it returns a
  // nonzero error code.  The failure mode to expect is that too many
  // hooks have been set and the system has no more memory to store
  // them.  In that case, we simply accept that no printing will
  // happen at exit.  It doesn't make sense to retry; hooks can't be
  // _unset_, so memory will never get freed up.
  (void) atexit (allocationLoggerHook);

  const float_vec_type::size_type numEntries = 10;
  float_vec_type x_f (numEntries, 42.0);
  long_vec_type x_l (numEntries);
  std::copy (x_f.begin (), x_f.end (), x_l.begin ());

#ifdef HAVE_TEUCHOS_MPI
  (void) MPI_Finalize ();
#endif // HAVE_TEUCHOS_MPI

  return EXIT_SUCCESS;
}

