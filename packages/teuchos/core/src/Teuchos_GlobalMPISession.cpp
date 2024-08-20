// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Assert.hpp"

// The header file does not at all depend on MPI routines or types,
// so we can defer inclusion of mpi.h to here.  This also fixes Bug
// 5631:  https://software.sandia.gov/bugzilla/show_bug.cgi?id=5631
#ifdef HAVE_MPI
#  include "mpi.h"
#endif

#ifdef HAVE_TEUCHOSCORE_KOKKOS
#  include "Kokkos_Core.hpp"
#endif // HAVE_TEUCHOSCORE_KOKKOS



namespace Teuchos {


bool GlobalMPISession::haveMPIState_ = false;
bool GlobalMPISession::mpiIsFinalized_ = false;
int GlobalMPISession::rank_ = 0 ;
int GlobalMPISession::nProc_ = 1 ;

#ifdef HAVE_TEUCHOSCORE_KOKKOS

// We have to invoke the std::vector's constructor here,
// because it's a class (static) variable.
std::vector<std::string> GlobalMPISession::argvCopy_;

#endif // HAVE_TEUCHOSCORE_KOKKOS


GlobalMPISession::GlobalMPISession( int* argc, char*** argv, std::ostream *out )
{
  std::ostringstream oss;

  // Above is used to create all output before sending to *out to avoid
  // jumbled parallel output between processors

#ifdef HAVE_MPI

  int mpierr = 0;

  // Assert that MPI is not already initialized
  int mpiHasBeenStarted = 0;
  MPI_Initialized(&mpiHasBeenStarted);
  if (mpiHasBeenStarted) {
    if (out) {
      *out << "GlobalMPISession(): Error, MPI_Intialized() return true,"
           << " calling std::terminate()!\n"
           << std::flush;
    }
    std::terminate();
  }

  // Initialize MPI
  mpierr = ::MPI_Init(argc, (char ***) argv);
  if (mpierr != 0) {
    if (out) {
      *out << "GlobalMPISession(): Error, MPI_Init() returned error code="
           << mpierr << "!=0, calling std::terminate()!\n"
           << std::flush;
    }
    std::terminate();
  }

  initialize(out); // Get NProc_ and rank_

  int nameLen;
  char procName[MPI_MAX_PROCESSOR_NAME];
  mpierr = ::MPI_Get_processor_name(procName, &nameLen);
  if (mpierr != 0) {
    if (out) {
      *out << "GlobalMPISession():  Error, MPI_Get_processor_name() error code="
           << mpierr << "!=0, calling std::terminate()!\n"
           << std::flush;
    }
    std::terminate();
  }

  oss << "Teuchos::GlobalMPISession::GlobalMPISession(): started processor with name "
      << procName << " and rank " << rank_ << "!" << std::endl;

#else

  oss << "Teuchos::GlobalMPISession::GlobalMPISession(): started serial run"
      << std::endl;

#endif

#ifndef TEUCHOS_SUPPRESS_PROC_STARTUP_BANNER

  // See if we should suppress the startup banner
  bool printStartupBanner = true;
  const std::string suppress_option("--teuchos-suppress-startup-banner");
  for ( int opt_i = 0; opt_i < *argc; ++opt_i ) {
    if ( suppress_option == (*argv)[opt_i] ) {
      // We are suppressing the output!
      printStartupBanner = false;
      // Remove this option!
      // Note that (*argv)[*argc]==0 but convention so we copy it too!
      for( int i = opt_i; i < *argc; ++i )
        (*argv)[i] = (*argv)[i+1];
      --*argc;
    }
  }
  if (out && printStartupBanner) {
    *out << oss.str() << std::flush;
  }

#endif

#ifdef HAVE_TEUCHOSCORE_KOKKOS
  // mfh 15/16 Apr 2016: This is the one chance we get to save the
  // command-line arguments, so that we can (later) initialize Kokkos
  // with the correct number of threads as specified by (e.g.,) the
  // --kokkos-threads command-line argument.  We won't attempt to
  // initialize Kokkos now, because not all applications want Kokkos.
  // Some applications may also prefer to initialize Kokkos with their
  // own thread count.
  //
  // NOTE (mfh 15/16 Apr 2016): While static variables are not thread
  // safe in general, and this is not thread safe in particular, it
  // only makes sense to GlobalMPISession's constructor on a single
  // thread per MPI process anyway, because MPI_Init has the same
  // requirement.

  const int numArgs = *argc;
  argvCopy_.resize (numArgs);  
  for (int c = 0; c < numArgs; ++c) {
    argvCopy_[c] = std::string ((*argv)[c]); // deep copy
  }
#endif // HAVE_TEUCHOSCORE_KOKKOS
}

  
#ifdef HAVE_TEUCHOSCORE_KOKKOS  
std::vector<std::string> GlobalMPISession::getArgv ()
{
  return argvCopy_;
}
#endif // HAVE_TEUCHOSCORE_KOKKOS  

  
GlobalMPISession::~GlobalMPISession()
{

#ifdef HAVE_TEUCHOSCORE_KOKKOS
  try {
    if (Kokkos::is_initialized())
      Kokkos::finalize();
  }
  catch (const std::runtime_error& e) {
    std::cerr << "Kokkos::finalize failed:\n"
              << e.what() << "\n";
  }
#endif

  haveMPIState_ = false;
#ifdef HAVE_MPI
  const int mpierr = ::MPI_Finalize();
  mpiIsFinalized_ = (mpierr == 0);
  if (mpierr != 0)
    std::cerr << "Error code " << mpierr << " returned from MPI_Finalize()\n";
#else
  mpiIsFinalized_ = true;
#endif
}

void GlobalMPISession::abort() {
  justInTimeInitialize();
  #ifdef HAVE_MPI
    MPI_Abort(MPI_COMM_WORLD, MPI_ERR_UNKNOWN);
  #else
    std::abort();
  #endif
}

bool GlobalMPISession::mpiIsInitialized() {
  justInTimeInitialize();
  return haveMPIState_;
}


bool GlobalMPISession::mpiIsFinalized()
{
  return mpiIsFinalized_;
}


int GlobalMPISession::getRank()
{
  justInTimeInitialize();
  return rank_;
}


int GlobalMPISession::getNProc() {
  justInTimeInitialize();
  return nProc_;
}


void GlobalMPISession::barrier()
{
  justInTimeInitialize();
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}


int GlobalMPISession::sum(int localVal)
{
  justInTimeInitialize();
#ifdef HAVE_MPI
  int globalSum = -1;
  MPI_Allreduce(&localVal, &globalSum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  return globalSum;
#else
  return localVal;
#endif
}


void GlobalMPISession::allGather(int localVal, const ArrayView<int> &allVals)
{
  justInTimeInitialize();
  TEUCHOS_ASSERT_EQUALITY(allVals.size(), getNProc());
#ifdef HAVE_MPI
  MPI_Allgather( &localVal, 1, MPI_INT, allVals.getRawPtr(), 1, MPI_INT,
    MPI_COMM_WORLD);
#else
  allVals[0] = localVal;
#endif
}


// private


void GlobalMPISession::initialize( std::ostream *out )
{
#ifdef HAVE_MPI

  if(mpiIsFinalized_) {
    // MPI has aleady been finalized so we have a serial machine again!
    rank_ = 0;
    nProc_ = 1;
    return;
  }

  if(haveMPIState_) {
    return; // We already have what we need!
  }

  // We don't have the state of MPI so the constructor for this class must not
  // have been called.  However, if MPI has been called in another way we
  // can still get the state of MPI_COMM_WORLD here.

  int mpiHasBeenStarted = 0;
  MPI_Initialized(&mpiHasBeenStarted);

  if(!mpiHasBeenStarted)
    return;  // We have to give up and just leave NProc_ and rank_ at the default values.

  // Get the state of MPI
  // Don't throw exceptions here since this part of the code
  // is used by TEUCHOS_STANDARD_CATCH_STATEMENTS().
  // See bug #6192 <https://software.sandia.gov/bugzilla/show_bug.cgi?id=6192>.
  int mpierr = 0;
  mpierr = ::MPI_Comm_rank( MPI_COMM_WORLD, &rank_ );
  if (mpierr != 0) {
    *out << "Error code=" << mpierr << " detected in MPI_Comm_rank()"
         << std::endl;
  }

  mpierr = ::MPI_Comm_size( MPI_COMM_WORLD, &nProc_ );
  if (mpierr != 0) {
    *out << "Error code=" << mpierr << " detected in MPI_Comm_size()"
         << std::endl;
  }

  haveMPIState_ = true;
  mpiIsFinalized_ = false;

#endif // HAVE_MPI

}


void GlobalMPISession::justInTimeInitialize()
{
  if(!haveMPIState_)
    initialize(&std::cerr);
}


} // namespace Teuchos
