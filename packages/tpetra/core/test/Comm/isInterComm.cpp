// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Details_isInterComm.hpp"
#ifdef HAVE_TPETRACORE_MPI
#  include "Tpetra_Details_extractMpiCommFromTeuchos.hpp"
#  include "Teuchos_DefaultMpiComm.hpp"
#endif // HAVE_TPETRACORE_MPI
#include "Teuchos_DefaultSerialComm.hpp"
#include "Tpetra_TestingUtilities.hpp"

namespace { // (anonymous)

using Teuchos::Comm;
using Teuchos::RCP;
using Teuchos::rcp;
using std::endl;
#ifdef HAVE_TPETRACORE_MPI
using Teuchos::MpiComm;
#else
using Teuchos::SerialComm;
#endif // HAVE_TPETRACORE_MPI

/// \brief Test Tpetra::Details::isInterComm.
///
/// \param out [out] Output stream; valid (writeable) only on Process
///   0 in the input communicator.
/// \param origComm [in] Communicator over which to do the test.  This
///   MUST be an intracommunicator.
void
testInterComm (bool& success,
               std::ostream& out,
               const Teuchos::Comm<int>& origCommWrapped)
{
  // lclSuccess: Local success status.  0 means a failure happened on
  //   the calling process.
  // gblSuccess [in/out] Global success status.  0 means a failure
  //   happened on some process in the input communicator.
  int lclSuccess = 1; // to be updated below
  int gblSuccess = 0; // output argument (see below)

#ifdef HAVE_TPETRACORE_MPI
  const int myRank = origCommWrapped.getRank ();
  const int numProcs = origCommWrapped.getSize ();

  out << "Test Tpetra::Details::isInterComm with MPI enabled" << endl;
  Teuchos::OSTab tab1 (out);

  MPI_Comm origComm =
    Tpetra::Details::extractMpiCommFromTeuchos (origCommWrapped);

  if (myRank == 0) {
    out << "Call MPI_Comm_test_inter to make sure that input communicator "
      "is not an intercomm" << endl;
  }
  {
    int isInter = 0;
    // This is a "local routine" (MPI 3.1, Section 6.6.1, p. 259).
    (void) MPI_Comm_test_inter (origComm, &isInter);
    lclSuccess = (isInter == 0);
    (void) MPI_Allreduce (&lclSuccess, &gblSuccess, 1,
                          MPI_INT, MPI_SUM, origComm);
    if (gblSuccess == 0) {
      if (myRank == 0) {
        out << "MPI_Comm_test_inter claims on at least one process in the "
          "original communicator, that the original communicator itself IS an "
          "intercommunicator!  This means that MPI itself is messed up." << endl;
      }
      success = false;
      return; // no sense in continuing; MPI is messed up
    }
    else {
      if (myRank == 0) {
        out << "OK" << endl;
      }
    }
  }

  if (myRank == 0) {
    out << "Test isInterComm on the input communicator "
      "(should return false)" << endl;
  }
  {
    const bool isInter = Tpetra::Details::isInterComm (origCommWrapped);
    lclSuccess = isInter ? 0 : 1;
    (void) MPI_Allreduce (&lclSuccess, &gblSuccess, 1,
                          MPI_INT, MPI_SUM, origComm);
    if (gblSuccess == 0) {
      if (myRank == 0) {
        out << "Tpetra::Details::isInterComm claims on at least one process in "
          "the original communicator, that the original communicator itself IS "
          "an intercommunicator!" << endl;
      }
      success = false;
      //return; // don't have to return immediately; let's see what happens
    }
    else {
      if (myRank == 0) {
        out << "OK" << endl;
      }
    }
  }

  if (myRank == 0) {
    out << "Duplicate the input communicator" << endl;
  }
  // The MPI standard says that it's better NOT to use MPI_COMM_WORLD
  // directly as the input to MPI_Intercomm_create.  Instead, one
  // should duplicate MPI_COMM_WORLD.
  MPI_Comm peerComm;
  (void) MPI_Comm_dup (origComm, &peerComm);

  if (myRank == 0) {
    out << "Call MPI_Comm_test_inter to make sure that the duplicated "
      "input communicator is not an intercomm" << endl;
  }
  {
    int isInter = 0;
    // This is a "local routine" (MPI 3.1, Section 6.6.1, p. 259).
    (void) MPI_Comm_test_inter (peerComm, &isInter);
    lclSuccess = (isInter == 0);
    (void) MPI_Allreduce (&lclSuccess, &gblSuccess, 1,
                          MPI_INT, MPI_SUM, origComm);
    if (gblSuccess == 0) {
      if (myRank == 0) {
        out << "MPI_Comm_test_inter claims on at least one process in the "
          "original communicator, that the result of MPI_Comm_dup on the "
          "original communicator IS an intercommunicator!" << endl;
      }
      success = false;
      return; // no sense in continuing; MPI is messed up
    }
    else {
      if (myRank == 0) {
        out << "OK" << endl;
      }
    }
  }

  if (myRank == 0) {
    out << "Test isInterComm on the duplicated input communicator "
      "(should return false)" << endl;
  }
  {
    RCP<const Comm<int> > peerCommWrapped = rcp (new MpiComm<int> (peerComm));
    const bool isInter = Tpetra::Details::isInterComm (*peerCommWrapped);
    lclSuccess = isInter ? 0 : 1;
    (void) MPI_Allreduce (&lclSuccess, &gblSuccess, 1,
                          MPI_INT, MPI_SUM, origComm);
    if (gblSuccess == 0) {
      if (myRank == 0) {
        out << "Tpetra::Details::isInterComm claims on at least one process in "
          "the original communicator, that the result of MPI_Comm_dup on the "
          "original communicator IS an intercommunicator!" << endl;
      }
      success = false;
      //return; // don't have to return immediately; let's see what happens
    }
    else {
      if (myRank == 0) {
        out << "OK" << endl;
      }
    }
  }

  if (numProcs >= 2) {
    const int color = (myRank % 2);
    const int key = 0;
    MPI_Comm splitComm;
    (void) MPI_Comm_split (peerComm, color, key, &splitComm);

    if (myRank == 0) {
      out << "Test MPI_Comm_test_inter on the result of MPI_Comm_split on the "
        "result of MPI_Comm_dup on the original communicator (should return "
        "false)" << endl;
    }
    {
      int isInter = 0;
      // This is a "local routine" (MPI 3.1, Section 6.6.1, p. 259).
      (void) MPI_Comm_test_inter (splitComm, &isInter);
      lclSuccess = (isInter == 0);
      (void) MPI_Allreduce (&lclSuccess, &gblSuccess, 1,
                            MPI_INT, MPI_SUM, origComm);
      if (gblSuccess == 0) {
        if (myRank == 0) {
          out << "MPI_Comm_test_inter claims on at least one process in the "
            "original communicator, that the result of MPI_Comm_split on the "
            "result of MPI_Comm_dup on the original communicator IS an "
            "intercommunicator!" << endl;
        }
        success = false;
        return; // no sense in continuing; MPI is messed up
      }
      else {
        if (myRank == 0) {
          out << "OK" << endl;
        }
      }
    }

    if (myRank == 0) {
      out << "Test isInterComm on the result of MPI_Comm_split on the "
        "duplicated input communicator (should return false)" << endl;
    }
    {
      RCP<const Comm<int> > splitCommWrapped = rcp (new MpiComm<int> (splitComm));
      const bool isInter = Tpetra::Details::isInterComm (*splitCommWrapped);
      lclSuccess = isInter ? 0 : 1;

      (void) MPI_Allreduce (&lclSuccess, &gblSuccess, 1,
                            MPI_INT, MPI_SUM, origComm);
      if (gblSuccess == 0) {
        if (myRank == 0) {
          out << "Tpetra::Details::isInterComm claims on at least one process in "
            "the original communicator, that the result of MPI_Comm_split on the "
            "MPI_Comm_dup of the original communicator IS an intercommunicator!"
              << endl;
        }
        success = false;
        //return; // don't have to return immediately; let's see what happens
      }
      else {
        if (myRank == 0) {
          out << "OK" << endl;
        }
      }
    }

    if (myRank == 0) {
      out << "Create an intercommunicator using MPI_Intercomm_create" << endl;
    }
    const int localLeader = 0;
    const int remoteLeader = (color == 0) ? 1 : 0;
    const int tag = 42;
    MPI_Comm interComm;
    (void) MPI_Intercomm_create (splitComm, localLeader, peerComm,
                                 remoteLeader, tag, &interComm);

    if (myRank == 0) {
      out << "Test MPI_Comm_test_inter on the result of MPI_Intercomm_create "
        "(should return true)" << endl;
    }
    {
      int isInter = 0;
      // This is a "local routine" (MPI 3.1, Section 6.6.1, p. 259).
      (void) MPI_Comm_test_inter (interComm, &isInter);
      lclSuccess = (isInter != 0);
      (void) MPI_Allreduce (&lclSuccess, &gblSuccess, 1,
                            MPI_INT, MPI_SUM, origComm);
      if (gblSuccess == 0) {
        if (myRank == 0) {
          out << "MPI_Comm_test_inter claims on at least one process in the "
            "original communicator, that the result of MPI_Intercomm_create "
            "is NOT an intercommunicator!" << endl;
        }
        success = false;
        return; // no sense in continuing; MPI is messed up
      }
      else {
        if (myRank == 0) {
          out << "OK" << endl;
        }
      }
    }

    if (myRank == 0) {
      out << "Test isInterComm on the result of MPI_Intercomm_create "
        "(should return true)" << endl;
    }
    {
      RCP<const Comm<int> > interCommWrapped;
      // mfh 22 Nov 2016: Need a default tag, else MpiComm's
      // constructor will launch a collective to compute it, which
      // causes deadlock here (all processes reach this point).
      const int defaultTag = 31;
      using Teuchos::opaqueWrapper;
      interCommWrapped = rcp (new MpiComm<int> (opaqueWrapper (interComm),
                                                defaultTag));
      const bool isInter = Tpetra::Details::isInterComm (*interCommWrapped);
      lclSuccess = isInter ? 1 : 0;
      (void) MPI_Allreduce (&lclSuccess, &gblSuccess, 1,
                            MPI_INT, MPI_SUM, origComm);
      if (gblSuccess == 0) {
        if (myRank == 0) {
          out << "Tpetra::Details::isInterComm claims on at least one process "
            "in the original communicator, that the result of "
            "MPI_Intercomm_create on the original communicator is NOT an "
            "intercommunicator!" << endl;
        }
        success = false;
        return; // no sense in continuing; isInterComm is messed up
      }
      else {
        if (myRank == 0) {
          out << "OK" << endl;
        }
      }
    }

    (void) MPI_Comm_free (&interComm);
    (void) MPI_Comm_free (&splitComm);
  }

  (void) MPI_Comm_free (&peerComm);

#else // NOT HAVE_TPETRACORE_MPI

  out << "Test Tpetra::Details::isInterComm with MPI disabled" << endl;
  Teuchos::OSTab tab1 (out);

  out << "Test Tpetra::Details::isInterComm on the input communicator "
    "(should return false)" << endl;
  const bool isInter = Tpetra::Details::isInterComm (origCommWrapped);
  lclSuccess = isInter ? 0 : 1;
  gblSuccess = lclSuccess;
  if (gblSuccess == 0) {
    out << "Tpetra::Details::isInterComm reports that the original input "
      "communicator IS an intercomm!" << endl;
    success = false;
  }

  // Not much more we can test without MPI.

#endif // HAVE_TPETRACORE_MPI

  success = (gblSuccess != 0);
}


TEUCHOS_UNIT_TEST( isInterComm, basic )
{
  out << "Testing Tpetra::Details::isInterComm" << endl;
  Teuchos::OSTab tab1 (out);

#ifdef HAVE_TPETRACORE_MPI
  RCP<const Comm<int> > comm = rcp (new MpiComm<int> (MPI_COMM_WORLD));
#else
  RCP<const Comm<int> > comm = rcp (new Teuchos::SerialComm<int> ());
#endif // HAVE_TPETRACORE_MPI

  testInterComm (success, std::cerr, *comm);
  // Just to make sure that if the test fails, we still print
  // something that the unit test framework recognizes.
  TEST_ASSERT( success );
}

} // namespace (anonymous)
