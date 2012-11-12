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
#include "Teuchos_Assert.hpp"

// The header file does not at all depend on MPI routines or types,
// so we can defer inclusion of mpi.h to here.  This also fixes Bug
// 5631:  https://software.sandia.gov/bugzilla/show_bug.cgi?id=5631
#ifdef HAVE_MPI
#include "mpi.h"
#endif


namespace Teuchos {

bool GlobalMPISession::haveMPIState_ = false;
bool GlobalMPISession::mpiIsFinalized_ = false;
int GlobalMPISession::rank_ = 0 ;
int GlobalMPISession::nProc_ = 1 ;

GlobalMPISession::GlobalMPISession( int* argc, char*** argv, std::ostream *out )
{
  std::ostringstream oss;

  // Above is used to create all output before sending to *out to avoid
  // jumbled parallel output between processors

#ifdef HAVE_MPI
  // initialize MPI
        int mpiHasBeenStarted = 0, mpierr = 0;
        MPI_Initialized(&mpiHasBeenStarted);
        TEUCHOS_TEST_FOR_EXCEPTION_PRINT(
    mpiHasBeenStarted, std::runtime_error
    ,"Error, you can only call this constructor once!"
    ,out
    );

  mpierr = ::MPI_Init (argc, (char ***) argv);
  TEUCHOS_TEST_FOR_EXCEPTION_PRINT(
    mpierr != 0, std::runtime_error
    ,"Error code=" << mpierr << " detected in GlobalMPISession::GlobalMPISession(argc,argv)"
    ,out
    );

  initialize(out); // Get NProc_ and rank_

  int nameLen;
        char procName[MPI_MAX_PROCESSOR_NAME];
  mpierr = ::MPI_Get_processor_name(procName,&nameLen);
  TEUCHOS_TEST_FOR_EXCEPTION_PRINT(
    mpierr != 0, std::runtime_error
    ,"Error code=" << mpierr << " detected in MPI_Get_processor_name()"
    ,out
    );

  oss << "Teuchos::GlobalMPISession::GlobalMPISession(): started processor with name "
      << procName << " and rank " << rank_ << "!" << std::endl;

#else
  oss << "Teuchos::GlobalMPISession::GlobalMPISession(): started serial run" << std::endl;
#endif
#ifndef TEUCHOS_SUPPRESS_PROC_STARTUP_BANNER
  // See if we should suppress the startup banner
  bool printStartupBanner = true;
  const std::string suppress_option("--teuchos-suppress-startup-banner");
  for( int opt_i = 0; opt_i < *argc; ++opt_i ) {
    if( suppress_option == (*argv)[opt_i] ) {
      // We are suppressing the output!
      printStartupBanner = false;
      // Remove this option!
      // Note that (*argv)[*argc]==0 but convention so we copy it too!
      for( int i = opt_i; i < *argc; ++i )
        (*argv)[i] = (*argv)[i+1];
      --*argc;
    }
  }
  if( out && printStartupBanner )
    *out << oss.str() << std::flush;
#endif
}

GlobalMPISession::~GlobalMPISession()
{
  haveMPIState_ = false;
  mpiIsFinalized_ = true;
#ifdef HAVE_MPI
  int mpierr = ::MPI_Finalize();
  TEUCHOS_TEST_FOR_EXCEPTION_PRINT(
    mpierr != 0, std::runtime_error
    ,"Error code=" << mpierr << " detected in MPI_Finalize()"
    ,&std::cerr
    );
#endif
}

bool GlobalMPISession::mpiIsInitialized() {
  if(!haveMPIState_)
    initialize(&std::cerr);
  return haveMPIState_;
}

bool GlobalMPISession::mpiIsFinalized()
{
  return mpiIsFinalized_;
}

int GlobalMPISession::getRank()
{
  if(!haveMPIState_)
    initialize(&std::cerr);
  return rank_;
}

int GlobalMPISession::getNProc() {
  if(!haveMPIState_)
    initialize(&std::cerr);
  return nProc_;
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

  if(haveMPIState_)
    return; // We already have what we need!

  // We don't have the state of MPI so the constructor for this class must not
  // have been called.  However, if MPI has been called in another way we
  // can still get the state of MPI_COMM_WORLD here.

  int mpiHasBeenStarted = 0, mpierr = 0;
  MPI_Initialized(&mpiHasBeenStarted);

  if(!mpiHasBeenStarted)
    return;  // We have to give up and just leave NProc_ and rank_ at the default values.

  // Get the state of MPI

  mpierr = ::MPI_Comm_rank( MPI_COMM_WORLD, &rank_ );
  TEUCHOS_TEST_FOR_EXCEPTION_PRINT(
    mpierr != 0, std::runtime_error
    ,"Error code=" << mpierr << " detected in MPI_Comm_rank()"
    ,out
    );

  mpierr = ::MPI_Comm_size( MPI_COMM_WORLD, &nProc_ );
  TEUCHOS_TEST_FOR_EXCEPTION_PRINT(
    mpierr != 0, std::runtime_error
    ,"Error code=" << mpierr << " detected in MPI_Comm_size()"
    ,out
    );

  haveMPIState_ = true;
  mpiIsFinalized_ = false;

#endif // HAVE_MPI

}

} // namespace Teuchos
