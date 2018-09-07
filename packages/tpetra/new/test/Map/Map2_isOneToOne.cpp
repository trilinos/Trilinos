/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
// @HEADER
*/

#include <Tpetra_TestingUtilities.hpp>
#include <Tpetra_Core.hpp>
#include <TpetraNew_Map.hpp>

namespace {

  using Teuchos::Array;
  using Teuchos::as;
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using std::cerr;
  using std::endl;

  using comm_type = Teuchos::Comm<int>;
  using map_type = TpetraNew::Map;
  using LO = map_type::local_ordinal_type;  
  using GO = map_type::global_ordinal_type;

  TEUCHOS_UNIT_TEST(Map_new, isOneToOne_contig_uniform_distributed)
  {
    out << "Testing Map::isOneToOne with a contiguous uniform Map" << endl;
    Teuchos::OSTab tab0 (out);


    RCP<const comm_type> comm = Tpetra::getDefaultComm ();

    const int numProcs = comm->getSize ();
    const int myRank = comm->getRank ();
    out << "Comm has " << numProcs << " process"
	<< (numProcs != 1 ? "es" : "") << endl;

    int lclSuccess = 1;
    std::ostringstream err;

    // Make a contiguous uniform Map.  Just for variety, give it a
    // number of global indices not evenly divisible by the number of
    // processes in the communicator, and an index base different than
    // the usual value of zero.
    const GO numGlobalIndices = GO (numProcs) * 10 + GO (1);
    const GO indexBase = 1;
    map_type map;
    try {
      map = map_type (numGlobalIndices, indexBase, comm,
		      Tpetra::GloballyDistributed);
      lclSuccess = 1;
    }
    catch (std::exception& e) {
      lclSuccess = 0;
      err << "On Process " << myRank << ": contiguous uniform Map "
	"constructor raised an exception with the following message: "
	<< e.what () << endl;
    }

    // If Map construction failed on _any_ process, first print all
    // the error messages on all processes, then raise an exception on
    // all processes.
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      for (int p = 0; p < numProcs; ++p) {
        if (p == myRank) {
          if (lclSuccess == 1) {
            cerr << "Process " << myRank << " succeeded" << endl;
          }
          else {
            cerr << err.str ();
          }
        }
        comm->barrier (); // Give output time to finish
        comm->barrier ();
        comm->barrier ();
      }
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Failed to construct Map!");
    }

    // Make sure that the Map is one to one on _all_ processes.
    const bool isOneToOne = map.isOneToOne ();
    lclSuccess = isOneToOne ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST(gblSuccess, 1);
  }


  TEUCHOS_UNIT_TEST(Map_new, isOneToOne_serial_contig_uniform_distributed)
  {
    out << "Testing Map::isOneToOne with a contiguous uniform Map over one process" << endl;
    Teuchos::OSTab tab0 (out);

    RCP<const comm_type> comm = rcp (new Teuchos::SerialComm<int> ());

    const int numProcs = comm->getSize ();
    const int myRank = comm->getRank ();
    out << "Comm has " << numProcs << " process"
	<< (numProcs != 1 ? "es" : "") << endl;

    int lclSuccess = 1;
    std::ostringstream err;

    // Make a contiguous uniform Map.  Just for variety, give it a
    // number of global indices not evenly divisible by the number of
    // processes in the communicator, and an index base different than
    // the usual value of zero.
    const GO numGlobalIndices = GO (numProcs) * 10 + GO (1);
    const GO indexBase = 1;
    map_type map;
    try {
      map = map_type (numGlobalIndices, indexBase, comm,
		      Tpetra::GloballyDistributed);
      lclSuccess = 1;
    }
    catch (std::exception& e) {
      lclSuccess = 0;
      err << "On Process " << myRank << ": contiguous uniform Map "
	"constructor raised an exception with the following message: "
	<< e.what () << endl;
    }

    // If Map construction failed on _any_ process, first print all
    // the error messages on all processes, then raise an exception on
    // all processes.
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      for (int p = 0; p < numProcs; ++p) {
        if (p == myRank) {
          if (lclSuccess == 1) {
            cerr << "Process " << myRank << " succeeded" << endl;
          }
          else {
            cerr << err.str ();
          }
        }
        comm->barrier (); // Give output time to finish
        comm->barrier ();
        comm->barrier ();
      }
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Failed to construct Map!");
    }

    // Make sure that the Map is one to one on _all_ processes.
    const bool isOneToOne = map.isOneToOne ();
    lclSuccess = isOneToOne ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST(gblSuccess, 1);
  }


  TEUCHOS_UNIT_TEST(Map_new, isOneToOne_contig_uniform_replicated)
  {
    out << "Testing Map::isOneToOne with a contiguous uniform replicated Map" << endl;
    Teuchos::OSTab tab0 (out);

    RCP<const comm_type> comm = Tpetra::getDefaultComm ();

    const int numProcs = comm->getSize ();
    const int myRank = comm->getRank ();
    out << "Comm has " << numProcs << " process"
	<< (numProcs != 1 ? "es" : "") << endl;

    int lclSuccess = 1;
    std::ostringstream err;

    // Make a contiguous uniform replicated Map.  Just for variety,
    // give it a number of global indices not evenly divisible by the
    // number of processes in the communicator, and an index base
    // different than the usual value of zero.
    const GO numGlobalIndices = GO (numProcs) * 10 + GO (1);
    const GO indexBase = 1;
    map_type map;
    try {
      map = map_type (numGlobalIndices, indexBase, comm,
		      Tpetra::LocallyReplicated);
      lclSuccess = 1;
    }
    catch (std::exception& e) {
      lclSuccess = 0;
      err << "On Process " << myRank << ": contiguous uniform replicated "
	"Map constructor raised an exception with the following message: "
        << e.what () << endl;
    }

    // If Map construction failed on _any_ process, first print all
    // the error messages on all processes, then raise an exception on
    // all processes.
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      for (int p = 0; p < numProcs; ++p) {
        if (p == myRank) {
          if (lclSuccess == 1) {
            cerr << "Process " << myRank << " succeeded" << endl;
          }
          else {
            cerr << err.str ();
          }
        }
        comm->barrier (); // Give output time to finish
        comm->barrier ();
        comm->barrier ();
      }
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Failed to construct Map!");
    }

    // The Map should only be one to one if numProcs == 1.  Make sure
    // that the value of isOneToOne is the same over all processes.
    const bool isOneToOne = map.isOneToOne ();
    lclSuccess = ((numProcs == 1 && isOneToOne) || (numProcs != 1 && ! isOneToOne)) ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST(gblSuccess, 1);
  }


  TEUCHOS_UNIT_TEST(Map_new, isOneToOne_serial_contig_uniform_replicated)
  {
    out << "Testing Map::isOneToOne with a contiguous uniform "
      "replicated Map over one process" << endl;
    Teuchos::OSTab tab0 (out);

    RCP<const comm_type> comm = rcp (new Teuchos::SerialComm<int> ());

    const int numProcs = comm->getSize ();
    const int myRank = comm->getRank ();
    out << "Comm has " << numProcs << " process"
	<< (numProcs != 1 ? "es" : "") << endl;

    int lclSuccess = 1;
    std::ostringstream err;

    // Make a contiguous uniform replicated Map.  Just for variety,
    // give it a number of global indices not evenly divisible by the
    // number of processes in the communicator, and an index base
    // different than the usual value of zero.
    const GO numGlobalIndices = GO (numProcs) * 10 + GO (1);
    const GO indexBase = 1;
    map_type map;
    try {
      map = map_type (numGlobalIndices, indexBase, comm,
		      Tpetra::LocallyReplicated);
      lclSuccess = 1;
    }
    catch (std::exception& e) {
      lclSuccess = 0;
      err << "On Process " << myRank << ": contiguous uniform "
	"replicated Map constructor raised an exception with the "
	"following message: " << e.what () << endl;
    }

    // If Map construction failed on _any_ process, first print all
    // the error messages on all processes, then raise an exception on
    // all processes.
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      for (int p = 0; p < numProcs; ++p) {
        if (p == myRank) {
          if (lclSuccess == 1) {
            cerr << "Process " << myRank << " succeeded" << endl;
          }
          else {
            cerr << err.str ();
          }
        }
        comm->barrier (); // Give output time to finish
        comm->barrier ();
        comm->barrier ();
      }
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Failed to construct Map!");
    }

    // The Map should only be one to one if numProcs == 1.  (This is
    // always true for a communicator over one process, but we test
    // nevertheless.)  Make sure that the value of isOneToOne is the
    // same over all processes.
    const bool isOneToOne = map.isOneToOne ();
    lclSuccess = ((numProcs == 1 && isOneToOne) || (numProcs != 1 && ! isOneToOne)) ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST(gblSuccess, 1);
  }


  TEUCHOS_UNIT_TEST(Map_new, isOneToOne_contig_nonuniform)
  {
    out << "Testing Map::isOneToOne with a contiguous nonuniform Map" << endl;
    Teuchos::OSTab tab0 (out);

    RCP<const comm_type> comm = Tpetra::getDefaultComm ();

    const int numProcs = comm->getSize ();
    const int myRank = comm->getRank ();
    out << "Comm has " << numProcs << " process" << (numProcs != 1 ? "es" : "") << endl;

    int lclSuccess = 1;
    std::ostringstream err;

    // Make a contiguous nonuniform Map.  Just for variety, give it an
    // index base different than the usual value of zero.
    const GO numGlobalIndices = Teuchos::OrdinalTraits<GO>::invalid ();
    // Just some arbitrary formula to ensure that different processes
    // have different numbers of indices.
    const LO numLocalIndices = (myRank == 0) ? 42 : (8 + (myRank % 3));
    const GO indexBase = 1;
    map_type map;
    try {
      map = map_type (numGlobalIndices, numLocalIndices, indexBase, comm);
      lclSuccess = 1;
    } catch (std::exception& e) {
      lclSuccess = 0;
      err << "Process " << myRank << ": contiguous nonuniform Map constructor "
        "raised an exception with the following message: " << e.what () << endl;
    }

    // If Map construction failed on _any_ process, first print all
    // the error messages on all processes, then raise an exception on
    // all processes.
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      for (int p = 0; p < numProcs; ++p) {
        if (p == myRank) {
          if (lclSuccess == 1) {
            cerr << "Process " << myRank << ": Map constructor did not raise "
              "an exception" << endl;
          }
          else {
            cerr << err.str ();
          }
        }
        comm->barrier (); // Give output time to finish
        comm->barrier ();
        comm->barrier ();
      }
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Failed to construct Map!");
    }

    // Make sure that the Map is one to one on _all_ processes.
    const bool isOneToOne = map.isOneToOne ();
    lclSuccess = isOneToOne ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST(gblSuccess, 1);
  }


  TEUCHOS_UNIT_TEST(Map_new, isOneToOne_serial_contig_nonuniform)
  {
    out << "Testing Map::isOneToOne with a contiguous nonuniform Map over one process" << endl;
    Teuchos::OSTab tab0 (out);

    RCP<const comm_type> comm = rcp (new Teuchos::SerialComm<int> ());

    const int numProcs = comm->getSize ();
    const int myRank = comm->getRank ();
    out << "Comm has " << numProcs << " process" << (numProcs != 1 ? "es" : "") << endl;

    int lclSuccess = 1;
    std::ostringstream err;

    // Make a contiguous nonuniform Map.  Just for variety, give it an
    // index base different than the usual value of zero.
    const GO numGlobalIndices = Teuchos::OrdinalTraits<GO>::invalid ();
    // Just some arbitrary formula to ensure that different processes
    // have different numbers of indices.
    const LO numLocalIndices = (myRank == 0) ? 42 : (8 + (myRank % 3));
    const GO indexBase = 1;
    map_type map;
    try {
      map = map_type (numGlobalIndices, numLocalIndices, indexBase, comm);
      lclSuccess = 1;
    } catch (std::exception& e) {
      lclSuccess = 0;
      err << "Process " << myRank << ": contiguous nonuniform Map constructor "
        "raised an exception with the following message: " << e.what () << endl;
    }

    // If Map construction failed on _any_ process, first print all
    // the error messages on all processes, then raise an exception on
    // all processes.
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      for (int p = 0; p < numProcs; ++p) {
        if (p == myRank) {
          if (lclSuccess == 1) {
            cerr << "Process " << myRank << ": Map constructor did not raise "
              "an exception" << endl;
          }
          else {
            cerr << err.str ();
          }
        }
        comm->barrier (); // Give output time to finish
        comm->barrier ();
        comm->barrier ();
      }
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Failed to construct Map!");
    }

    // Make sure that the Map is one to one on _all_ processes.
    const bool isOneToOne = map.isOneToOne ();
    lclSuccess = isOneToOne ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST(gblSuccess, 1);
  }


  TEUCHOS_UNIT_TEST(Map_new, isOneToOne_noncontig_oneToOne)
  {
    out << "Testing Map::isOneToOne with a noncontiguous, one-to-one Map" << endl;
    Teuchos::OSTab tab0 (out);

    RCP<const comm_type> comm = Tpetra::getDefaultComm ();

    const int numProcs = comm->getSize ();
    const int myRank = comm->getRank ();
    out << "Comm has " << numProcs << " process" << (numProcs != 1 ? "es" : "") << endl;

    int lclSuccess = 1;
    std::ostringstream err;

    // Make a noncontiguous, one-to-one Map.  The first 10 processes
    // get 10 indices, and any processes that remain get zero indices.
    // Process 0 gets indices 0, 10, 20, ..., 90.  Process 1 (if it
    // exists) gets indices 1, 11, 21, ..., 91.  Continue the pattern;
    // Process 9 (if it exists) gets indices 9, 19, 29, ..., 99.  Any
    // remaining processes get zero indices, so that the Map stays one
    // to one.
    const GO numGlobalIndices = std::min (numProcs, 10) * 10;
    const LO numLocalIndices = (myRank < 10) ? 10 : 0;
    Array<GO> myGlobalIndices (numLocalIndices);
    for (LO k = 0; k < numLocalIndices; ++k) {
      myGlobalIndices[k] = (GO (k) * GO (numLocalIndices)) + GO (myRank);
    }
    const GO indexBase = 0;

    map_type map;
    try {
      map = map_type (numGlobalIndices, myGlobalIndices (), indexBase, comm);
      lclSuccess = 1;
    } catch (std::exception& e) {
      lclSuccess = 0;
      err << "Process " << myRank << ": noncontiguous Map constructor raised "
        "an exception with the following message: " << e.what () << endl;
    }

    // If Map construction failed on _any_ process, first print all
    // the error messages on all processes, then raise an exception on
    // all processes.
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      for (int p = 0; p < numProcs; ++p) {
        if (p == myRank) {
          if (lclSuccess == 1) {
            cerr << "Process " << myRank << ": Map constructor did not raise "
              "an exception" << endl;
          }
          else {
            cerr << err.str ();
          }
        }
        comm->barrier (); // Give output time to finish
        comm->barrier ();
        comm->barrier ();
      }
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Failed to construct Map!");
    }

    // Make sure that the Map is one to one on _all_ processes.
    const bool isOneToOne = map.isOneToOne ();
    lclSuccess = isOneToOne ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST(gblSuccess, 1);
  }


  TEUCHOS_UNIT_TEST(Map_new, isOneToOne_noncontig_replicated)
  {
    out << "Testing Map::isOneToOne with a noncontiguous, replicated Map" << endl;
    Teuchos::OSTab tab0 (out);

    RCP<const comm_type> comm = Tpetra::getDefaultComm ();

    const int numProcs = comm->getSize ();
    const int myRank = comm->getRank ();
    out << "Comm has " << numProcs << " process" << (numProcs != 1 ? "es" : "") << endl;

    int lclSuccess = 1;
    std::ostringstream err;

    // Make a noncontiguous replicated Map.  All processes get indices
    // [10, 0, 5, 2].  We use an out-of-order sequence like this to
    // avoid any optimizations that could turn a contiguous increasing
    // sequence into a contiguous Map.

    const LO numLocalIndices = 4;
    const GO numGlobalIndices = numLocalIndices * numProcs;
    Array<GO> myGlobalIndices (numLocalIndices);
    myGlobalIndices[0] = 10;
    myGlobalIndices[1] = 0;
    myGlobalIndices[2] = 5;
    myGlobalIndices[3] = 2;
    const GO indexBase = 0;

    map_type map;
    try {
      map = map_type (numGlobalIndices, myGlobalIndices (), indexBase, comm);
      lclSuccess = 1;
    } catch (std::exception& e) {
      lclSuccess = 0;
      err << "Process " << myRank << ": noncontiguous Map constructor raised "
        "an exception with the following message: " << e.what () << endl;
    }

    // If Map construction failed on _any_ process, first print all
    // the error messages on all processes, then raise an exception on
    // all processes.
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      for (int p = 0; p < numProcs; ++p) {
        if (p == myRank) {
          if (lclSuccess == 1) {
            cerr << "Process " << myRank << ": Map constructor did not raise "
              "an exception" << endl;
          }
          else {
            cerr << err.str ();
          }
        }
        comm->barrier (); // Give output time to finish
        comm->barrier ();
        comm->barrier ();
      }
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Failed to construct Map!");
    }

    // The Map should only be one to one if numProcs == 1.  Make sure
    // that the value of isOneToOne is the same over all processes.
    const bool isOneToOne = map.isOneToOne ();
    lclSuccess = ((numProcs == 1 && isOneToOne) || (numProcs != 1 && ! isOneToOne)) ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST(gblSuccess, 1);
  }

  TEUCHOS_UNIT_TEST(Map_new, isOneToOne_noncontig_notOneToOne)
  {
    out << "Testing Map::isOneToOne with a noncontiguous, not one-to-one Map" << endl;
    Teuchos::OSTab tab0 (out);

    RCP<const comm_type> comm = Tpetra::getDefaultComm ();

    const int numProcs = comm->getSize ();
    const int myRank = comm->getRank ();
    out << "Comm has " << numProcs << " process" << (numProcs != 1 ? "es" : "") << endl;

    int lclSuccess = 1;
    std::ostringstream err;

    // Make a noncontiguous, not one-to-one Map.  It's only not one to
    // one if the communicator has only one process.  The first 10
    // processes get 10 elements, except for Proc 0, which gets 11.
    // Proc 0 has GIDs 0, 10, ..., 90 as well as 1, which also belongs
    // to Proc 1.
    const GO numGlobalIndices = std::min (numProcs, 10) * GO (10) + GO (1);
    const LO numLocalIndices = (myRank >= 10) ? 0 : ((myRank == 0) ? 11 : 10);
    Array<GO> myGlobalIndices (numLocalIndices);
    for (LO k = 0; k < LO (10); ++k) {
      myGlobalIndices[k] = (GO (k) * GO (numLocalIndices)) + GO (myRank);
    }
    if (myRank == 0) {
      myGlobalIndices[10] = 1; // belongs to Proc 1 as well
    }
    const GO indexBase = 0;

    map_type map;
    try {
      map = map_type (numGlobalIndices, myGlobalIndices (), indexBase, comm);
      lclSuccess = 1;
    } catch (std::exception& e) {
      lclSuccess = 0;
      err << "Process " << myRank << ": noncontiguous Map constructor raised "
        "an exception with the following message: " << e.what () << endl;
    }

    // If Map construction failed on _any_ process, first print all
    // the error messages on all processes, then raise an exception on
    // all processes.
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      for (int p = 0; p < numProcs; ++p) {
        if (p == myRank) {
          if (lclSuccess == 1) {
            cerr << "Process " << myRank << ": Map constructor did not raise "
              "an exception" << endl;
          }
          else {
            cerr << err.str ();
          }
        }
        comm->barrier (); // Give output time to finish
        comm->barrier ();
        comm->barrier ();
      }
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Failed to construct Map!");
    }

    // The Map should only be one to one if numProcs == 1.  Make sure
    // that the value of isOneToOne is the same over all processes.
    const bool isOneToOne = map.isOneToOne ();
    lclSuccess = ((numProcs == 1 && isOneToOne) || (numProcs != 1 && ! isOneToOne)) ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST(gblSuccess, 1);
  }

} // (anonymous)
