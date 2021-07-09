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

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#ifdef HAVE_TEUCHOS_MPI
#  include "Teuchos_DefaultMpiComm.hpp"
#endif // HAVE_TEUCHOS_MPI
#include "Teuchos_UnitTestHarness.hpp"


template<class PacketType>
bool
testScatter (bool& success, std::ostream& out,
             const int root, const Teuchos::Comm<int>& comm)
{
  using Teuchos::scatter;
  using Teuchos::TypeNameTraits;
  using std::endl;
  typedef PacketType packet_type;

  // Teuchos constructs the output stream such that it only prints on
  // Process 0 anyway.
  Teuchos::OSTab tab0 (out);
  out << "Testing Teuchos::scatter<int, "
      << TypeNameTraits<packet_type>::name ()
      << "> with root = " << root << endl;
  Teuchos::OSTab tab1 (out);

  int lclSuccess = 1;
  int gblSuccess = lclSuccess;

#ifdef HAVE_TEUCHOS_MPI
  using Teuchos::MpiComm;
  int errCode = MPI_SUCCESS;

  const MpiComm<int>* mpiComm = dynamic_cast<const MpiComm<int>* > (&comm);
  TEUCHOS_TEST_FOR_EXCEPTION
    (mpiComm == NULL, std::logic_error, "Building with MPI, but default "
     "communicator is not a Teuchos::MpiComm!");
  MPI_Comm rawMpiComm = * (mpiComm->getRawMpiComm ());
#endif // HAVE_TEUCHOS_MPI

  const packet_type ZERO = Teuchos::ScalarTraits<packet_type>::zero ();
  const packet_type ONE = Teuchos::ScalarTraits<packet_type>::one ();
  const int myRank = comm.getRank ();
  const int numProcs = comm.getSize ();
  const int sendCount = 10;
  const int recvCount = 10;

  out << "Initializing receive buffer (on all processes)" << endl;

  // Set up the receive buffer on each process.
  Teuchos::Array<packet_type> recvBuf (recvCount);
  for (int i = 0; i < recvCount; ++i) {
    recvBuf[i] = ZERO;
  }

  out << "Filling send buffer (on root process only)" << endl;

  // The send buffer is only needed on the root process.  It exists
  // elsewhere for syntax's sake, but is only allocated with nonzero
  // size on the root process.
  Teuchos::Array<packet_type> sendBuf;
  if (myRank == root) {
    sendBuf.resize (sendCount * numProcs);
    try {
      // Fill the send buffer.  The root process will send (p+1)*1, (p+1)*2,
      // ..., (p+1)*sendCount to Process p.  This helps us test whether the
      // scatter worked correctly.
      for (int p = 0; p < numProcs; ++p) {
        for (int k = 0; k < sendCount; ++k) {
          const packet_type val = static_cast<packet_type> (p+1) *
            (ONE + static_cast<packet_type> (k));
          sendBuf[p*sendCount + k] = val;
        }
      }
    } catch (std::exception& e) {
      std::cerr << "Root process " << root << " threw an exception: "
                << e.what () << endl;
      lclSuccess = 0;
    }
  }

#ifdef HAVE_TEUCHOS_MPI
  // Make sure that the root process didn't throw an exception when
  // filling the send buffer, by broadcasting from it to all the other
  // processes.
  gblSuccess = lclSuccess;
  errCode = MPI_Bcast (&gblSuccess, 1, MPI_INT, root, rawMpiComm);
  TEUCHOS_TEST_FOR_EXCEPTION
    (errCode != MPI_SUCCESS, std::logic_error, "MPI_Bcast failed!");
  TEUCHOS_TEST_FOR_EXCEPTION
    (gblSuccess != 1, std::logic_error, "Filling the send buffer failed on "
     "the root (" << root << ") process!  This probably indicates a bug in "
     "the test.");
#endif // HAVE_TEUCHOS_MPI

  // Invoke the function to test.
  out << "About to invoke scatter" << endl;
  scatter<int, packet_type> (sendBuf.getRawPtr (), sendCount,
                             recvBuf.getRawPtr (), recvCount,
                             root, comm);
  out << "Done with scatter" << endl;

  // Test whether each process got the expected results.
  for (int k = 0; k < recvCount; ++k) {
    const packet_type expectedVal = static_cast<packet_type> (myRank+1) *
      (ONE + static_cast<packet_type> (k));
    TEST_EQUALITY( recvBuf[k], expectedVal );
  }
  lclSuccess = success ? 1 : 0;
  gblSuccess = lclSuccess;

  // Use an all-reduce to check whether all processes got what they
  // should have gotten.
  //
  // Don't trust that any other Teuchos communicator wrapper functions
  // work here.  Instead, if building with MPI, use raw MPI.  If not
  // building with MPI, first test that comm has only one process,
  // then assume this in the test.

#ifdef HAVE_TEUCHOS_MPI
  errCode = MPI_Allreduce (&lclSuccess, &gblSuccess, 1,
                           MPI_INT, MPI_MIN, rawMpiComm);
  TEUCHOS_TEST_FOR_EXCEPTION
    (errCode != MPI_SUCCESS, std::logic_error, "MPI_Allreduce failed!");
#else // HAVE_TEUCHOS_MPI
  TEUCHOS_TEST_FOR_EXCEPTION
    (comm.getSize () != 1, std::logic_error, "Not building with MPI, but "
     "communicator has size = " << comm.getSize () << " != 1.  We don't know "
     "how to test this case.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (comm.getRank () != root, std::logic_error, "Not building with MPI, but "
     "the one process of rank " << comm.getRank () << " is not the root "
     "process of rank " << root << ".  We don't know how to test this case.");
#endif // HAVE_TEUCHOS_MPI

  return gblSuccess == 1;
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Comm, Scatter, PacketType )
{
  using Teuchos::Comm;
  using Teuchos::RCP;
  using std::endl;
  typedef PacketType packet_type;

  Teuchos::OSTab tab0 (out);
  out << "Testing Teuchos::scatter" << endl;
  Teuchos::OSTab tab1 (out);

  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm ();
  const int numProcs = comm->getSize ();

  // Make sure that it works for all possible root processes in the
  // communicator, not just Process 0.
  for (int root = 0; root < numProcs; ++root) {
    out << "Test for root = " << root << endl;
    const bool curSuccess = testScatter<packet_type> (success, out, root, *comm);
    TEST_EQUALITY_CONST( curSuccess, true );
    success = success && curSuccess;
  }

  comm->barrier (); // make sure that everybody finished
  out << "Done with test!" << endl;
}

//
// mfh 14 Jul 2015: We only need to test int for now.  See Bug 6375.
//
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Comm, Scatter, int )
