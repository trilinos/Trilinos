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

// Test that Teuchos::reduceAll allows aliasing input and output
// buffers, if the communicator is not an intercomm.  Motivation:
// https://github.com/trilinos/Trilinos/issues/1389

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#ifdef HAVE_TEUCHOS_MPI
#  include "mpi.h"
#  include "Teuchos_DefaultMpiComm.hpp"
#endif // HAVE_TEUCHOS_MPI
#include "Teuchos_UnitTestHarness.hpp"

namespace { // (anonymous)

#ifdef HAVE_TEUCHOS_MPI
MPI_Comm
getRawMpiCommFromTeuchosComm (const ::Teuchos::Comm<int>& comm)
{
  using ::Teuchos::Comm;
  using ::Teuchos::MpiComm;
  using ::Teuchos::RCP;

  const Comm<int>* commPtr = &comm;
  const MpiComm<int>* mpiCommPtr = dynamic_cast<const MpiComm<int>* > (commPtr);
  if (mpiCommPtr == NULL) {
    return MPI_COMM_SELF;
  }
  else {
    using ::Teuchos::OpaqueWrapper;
    using ::Teuchos::RCP;
    RCP<const OpaqueWrapper<MPI_Comm> > wrapper = mpiCommPtr->getRawMpiComm ();
    if (wrapper.is_null ()) {
      return MPI_COMM_SELF;
    }
    else {
      return *wrapper;
    }
  }
}
#endif // HAVE_TEUCHOS_MPI

TEUCHOS_UNIT_TEST( Comm, ReduceAllInPlace )
{
  using Teuchos::Comm;
  using Teuchos::outArg;
  using Teuchos::reduceAll;
  using Teuchos::REDUCE_SUM;
  using Teuchos::RCP;
  using std::endl;

  out << "Test Teuchos::reduceAll with both intracomms and intercomms" << endl;
  Teuchos::OSTab tab1 (out);

  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm ();
  const int numProcs = comm->getSize ();
  TEST_ASSERT( numProcs > 1 );
  if (numProcs > 1) {
    out << "This test requires > 1 processes in the input communicator "
      "in order to be meaningful." << endl;
    return;
  }

  typedef int reduce_type; // type of reductions' input & output, for this test
  reduce_type inputVal = 3;
  reduce_type outputVal = 0; // output of reduceAll
  const reduce_type expectedOutputVal = inputVal * numProcs;

  // Sanity check for reduceAll.
  reduceAll<int, reduce_type> (*comm, REDUCE_SUM, inputVal, outArg (outputVal));
  TEST_EQUALITY( outputVal, expectedOutputVal );

#ifdef HAVE_TEUCHOS_MPI
  MPI_Comm rawMpiComm = getRawMpiCommFromTeuchosComm (*comm);
  int isInter = 0;
  // This is a "local routine" (MPI 3.1, Section 6.6.1, p. 259).
  (void) MPI_Comm_test_inter (rawMpiComm, &isInter);
  if (isInter != 0) { // input communicator is an intercomm.
    out << "Input communicator is an intercomm; "
      "no point in continuing (not safe to exercise MPI_IN_PLACE)." << endl;
  }
  else { // input communicator is NOT an intercomm (is an intracomm)
    out << "Input communicator is NOT an intercomm; "
      "we may safely exercise MPI_IN_PLACE." << endl;

    // We hide use of MPI_IN_PLACE, by letting users alias the input
    // and output buffers.
    inputVal = 3;
    reduceAll<int, reduce_type> (*comm, REDUCE_SUM, inputVal, outArg (inputVal));
    TEST_EQUALITY( inputVal, expectedOutputVal );
  }

  // Test whether we succeeded on all processes.  This is particularly
  // important when exercising MPI_IN_PLACE, since incorrect use of
  // aliasing may cause only some processes to have incorrect results.
  int lclSuccess = success ? 1 : 0;
  int gblSuccess = 0; // output argument
  const int err = MPI_Allreduce (&lclSuccess, &gblSuccess, 1, MPI_INT, MPI_MIN, rawMpiComm);
  TEUCHOS_TEST_FOR_EXCEPTION
    (err != MPI_SUCCESS, std::logic_error, "MPI_Allreduce failed!");
  success = (gblSuccess == 1);
  TEST_ASSERT( gblSuccess == 1 );
#endif // HAVE_TEUCHOS_MPI
}

} // namespace (anonymous)
