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

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_DefaultSerialComm.hpp>
#include <Tpetra_Map.hpp>

#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION
#include <Tpetra_Map_def.hpp>
#include <Tpetra_Directory_def.hpp>
#endif // HAVE_TPETRA_EXPLICIT_INSTANTIATION

using Tpetra::global_size_t;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::as;
using Teuchos::Comm;
using Teuchos::MpiComm;
using Teuchos::null;
using Teuchos::outArg;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::REDUCE_MIN;
using Teuchos::reduceAll;
using Teuchos::SerialComm;
using Teuchos::toString;
using Teuchos::tuple;
using std::cerr;
using std::endl;

// This test is only meaningful in an MPI build.
TEUCHOS_UNIT_TEST( Map, replaceCommWithSubset )
{
  typedef int local_ordinal_type;
  typedef long global_ordinal_type;
  typedef Kokkos::SerialNode node_type;
  typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;
  typedef Array<global_ordinal_type>::size_type size_type;

  RCP<const Comm<int> > origComm = rcp (new MpiComm<int> (MPI_COMM_WORLD));
  RCP<node_type> node (new node_type);
  const int numProcs = origComm->getSize ();
  const int myRank = origComm->getRank ();

  // Create a Map in which all processes have a nonzero number of elements.
  const size_type numGidsPerProc = 3;
  const size_type myNumGids = numGidsPerProc;
  Array<global_ordinal_type> myGids (myNumGids);
  for (size_type k = 0; k < myNumGids; ++k) {
    myGids[k] = as<global_ordinal_type> (myRank) *
      as<global_ordinal_type> (numGidsPerProc) +
      as<global_ordinal_type> (k);
  }
  const global_size_t globalNumElts = as<global_size_t> (numGidsPerProc) *
    as<global_size_t> (numProcs);
  const global_ordinal_type indexBase = 0;
  RCP<const map_type> origMap (new map_type (globalNumElts, myGids (),
                                             indexBase, origComm, node));

  // Create a new communicator that excludes Proc 0.
  // This will exercise recomputing the index base.
  const int color = (myRank == 0) ? 0 : 1;
  const int key = 0;
  RCP<const Comm<int> > newComm = origComm->split (color, key);
  if (myRank == 0) {
    newComm = null;
  }
  // Create the new Map distributed over the subset communicator.
  RCP<const map_type> newMap = origMap->replaceCommWithSubset (newComm);

  // Test collectively for success, so the test doesn't hang on failure.
  int localSuccess = 1;
  std::ostringstream err;
  if (myRank == 0) {
    if (! newMap.is_null ()) {
      localSuccess = 0;
      err << "removeEmptyProcesses() should have returned null, but did not."
          << endl;
    }
  } else {
    if (newMap.is_null ()) {
      localSuccess = 0;
      err << "removeEmptyProcesses() should not have returned null, but did."
          << endl;
    } else {
      RCP<const Comm<int> > newComm = newMap->getComm ();
      if (newComm->getSize () != numProcs - 1) {
        localSuccess = 0;
        err << "New communicator should have " << (numProcs - 1)
            << " processes, but has " << newComm->getSize ()
            << " processes instead." << endl;
      }

      if (newMap->getGlobalNumElements () != origMap->getGlobalNumElements () - numGidsPerProc) {
        localSuccess = 0;
        err << "New Map has " << newMap->getGlobalNumElements () << " global "
            << "elements, but should have "
            << (origMap->getGlobalNumElements () - numGidsPerProc) << "." << endl;
      }

      if (newMap->getNodeNumElements () != origMap->getNodeNumElements ()) {
        localSuccess = 0;
        err << "New Map has " << newMap->getNodeNumElements () << " local "
            << "elements, but should have " << origMap->getNodeNumElements ()
            << "." << endl;
      }

      if (newMap->getIndexBase () != as<global_ordinal_type> (numGidsPerProc)) {
        localSuccess = 0;
        err << "New Map has index base " << newMap->getIndexBase ()
            << ", but should have index base " << numGidsPerProc << "." << endl;
      }

      ArrayView<const global_ordinal_type> myNewGids =
        newMap->getNodeElementList ();
      if (myNewGids.size () != myGids.size () ||
          ! std::equal (myNewGids.begin (), myNewGids.end (), myGids.begin ())) {
        localSuccess = 0;
        err << "New Map has local GID list " << toString (myNewGids) << ", but "
            << "should have local GID list " << toString (myGids ()) << "."
            << endl;
      }
    }
  }

  int globalSuccess = 0;
  reduceAll (*origComm, REDUCE_MIN, localSuccess, outArg (globalSuccess));
  if (globalSuccess == 0) {
    if (myRank == 0) {
      cerr << "TEST FAILED" << endl
           << "Error messages from each process:" << endl << endl;
    }
    for (int p = 0; p < numProcs; ++p) {
      if (myRank == p) {
        cerr << "Process " << myRank << ": " << err.str () << endl;
      }
      origComm->barrier (); // Give time for output to finish.
      origComm->barrier ();
      origComm->barrier ();
    }
  }
  TEST_EQUALITY(globalSuccess, 1);
}


