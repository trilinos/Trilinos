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

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>

using LO = Tpetra::Map<>::local_ordinal_type;
using GO = Tpetra::Map<>::global_ordinal_type;

namespace { // anonymous
  const GO sourceMap0[] = {91, 68, 90, 56, 77, 54, 75, 51};
  const GO sourceMap1[] = {97, 74, 88, 58, 94, 71, 89, 57};
  const GO sourceMap2[] = {65, 81, 93, 66, 50, 76, 79, 52, 92, 67, 78, 53};
  const GO sourceMap3[] = {63, 83, 99, 72, 98, 73, 64, 82, 96, 69, 95, 70};
  const GO sourceMap4[] = {41, 18, 40, 6, 27, 4, 25, 1};
  const GO sourceMap5[] = {15, 31, 43, 16, 42, 17, 0, 26, 29, 2, 28, 3};
  const GO sourceMap6[] = {37, 10, 47, 24, 34, 5, 38, 8, 44, 21, 39, 7};
  const GO sourceMap7[] = {9, 30, 35, 12, 13, 33, 49, 22, 36, 11, 48, 23, 14, 32, 46, 19, 45, 20};

  const LO numSrcGids[8] = {
    8,
    8,
    12,
    12,
    8,
    12,
    12,
    18
  };

  const GO targetMap0[] = {95, 70, 94, 71, 89, 57, 92, 67, 78, 53};
  const GO targetMap1[] = {95, 70, 94, 71, 89, 57, 98, 73, 28, 3, 98, 73, 27, 4, 25, 1};
  const GO targetMap2[] = {95, 70, 92, 67, 78, 53, 64, 82, 96, 69};
  const GO targetMap3[] = {95, 70, 98, 73, 64, 82, 96, 69, 28, 3, 98, 73, 0, 26, 29, 2};
  const GO targetMap4[] = {28, 3, 27, 4, 25, 1, 42, 17, 45, 20, 42, 17, 44, 21, 39, 7};
  const GO targetMap5[] = {28, 3, 0, 26, 29, 2, 42, 17, 45, 20, 42, 17, 14, 32, 46, 19};
  const GO targetMap6[] = {45, 20, 44, 21, 39, 7, 36, 11, 48, 23};
  const GO targetMap7[] = {45, 20, 14, 32, 46, 19, 36, 11, 48, 23};

  const LO numTgtGids[8] = {
    10,
    16,
    10,
    16,
    16,
    16,
    10,
    10
  };
} // namespace (anonymous)

int main (int argc, char *argv[])
{
  using Teuchos::ArrayView;
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using std::cerr;
  using std::cout;
  using std::endl;
  typedef Tpetra::Map<>::node_type    NT;
  typedef Tpetra::Map<LO,GO,NT>       map_type;
  typedef Tpetra::Import<LO,GO,NT>    import_type;
  typedef Tpetra::Vector<LO,LO,GO,NT> vec_type;
  typedef Tpetra::global_size_t GST; // Map's constructor needs this

  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  RCP<const Teuchos::Comm<int> > commWorld = Tpetra::getDefaultComm ();
  TEUCHOS_TEST_FOR_EXCEPTION(
    commWorld->getSize () < 8, std::logic_error,
    "This test must be run with at least 8 (preferably exactly 8) "
    "MPI processes.");

  // Create a communicator containing only the first 8 processes.  We
  // do this using Comm::split, and ignore the part of the split with
  // the remaining processes (whose ranks in MPI_COMM_WORLD are >= 8).
  RCP<const Teuchos::Comm<int> > comm;
  {
    const int color = (commWorld->getRank () < 8) ? 0 : 1;
    const int key = 0; // order processes by their current rank
    comm = commWorld->split (color, key);
  }

  bool success = true; // local error state on each process
  std::ostringstream err; // collect error messages on each process
  if (commWorld->getRank () < 8) {
    const int myRank = comm->getRank ();
    const GO* rawSrcGids = NULL;
    const GO* rawTgtGids = NULL;
    if (myRank == 0) {
      rawSrcGids = sourceMap0;
      rawTgtGids = targetMap0;
    } else if (myRank == 1) {
      rawSrcGids = sourceMap1;
      rawTgtGids = targetMap1;
    } else if (myRank == 2) {
      rawSrcGids = sourceMap2;
      rawTgtGids = targetMap2;
    } else if (myRank == 3) {
      rawSrcGids = sourceMap3;
      rawTgtGids = targetMap3;
    } else if (myRank == 4) {
      rawSrcGids = sourceMap4;
      rawTgtGids = targetMap4;
    } else if (myRank == 5) {
      rawSrcGids = sourceMap5;
      rawTgtGids = targetMap5;
    } else if (myRank == 6) {
      rawSrcGids = sourceMap6;
      rawTgtGids = targetMap6;
    } else if (myRank == 7) {
      rawSrcGids = sourceMap7;
      rawTgtGids = targetMap7;
    }
    ArrayView<const GO> mySrcGids (rawSrcGids, numSrcGids[myRank]);
    ArrayView<const GO> myTgtGids (rawTgtGids, numTgtGids[myRank]);

    const GST INV = Teuchos::OrdinalTraits<GST>::invalid ();
    const GO indexBase = 0;
    RCP<const map_type> sourceMap =
      rcp (new map_type (INV, mySrcGids, indexBase, comm));
    // The source Map of an Import must be one-to-one.  NOTE (mfh 18
    // Jul 2014) Calling isOneToOne sets up the Map's Directory.  That
    // could mask a bug.  On the other hand, I have run this test
    // without calling isOneToOne.
    bool srcIsOneToOne = true;
    try {
      srcIsOneToOne = sourceMap->isOneToOne ();
    }
    catch (std::exception& e) {
      err << "Map::isOneToOne threw an exception: " << e.what () << endl;
      cerr << "Map::isOneToOne threw an exception: " << e.what () << endl;
      success = false;
    }
    catch (...) {
      err << "Map::isOneToOne threw an exception of unknown type" << endl;
      cerr << "Map::isOneToOne threw an exception of unknown type" << endl;
      success = false;
    }
    // This is really a bug in the test, not a Tpetra source bug.
    TEUCHOS_TEST_FOR_EXCEPTION( ! srcIsOneToOne, std::logic_error,
      "The source Map is not one to one.  Please report this bug to "
      "the Tpetra developers.");

    RCP<const map_type> targetMap =
      rcp (new map_type (INV, myTgtGids, indexBase, comm));
    vec_type sourceVec (sourceMap);
    sourceVec.putScalar (static_cast<LO> (1));
    vec_type targetVec (targetMap);
    RCP<import_type> importer;

    try {
      importer = rcp (new import_type (sourceMap, targetMap));
    }
    catch (std::exception& e) {
      err << "Import constructor threw an exception: " << e.what () << endl;
      cerr << "Import constructor threw an exception: " << e.what () << endl;
      success = false;
    }
    catch (...) {
      err << "Import constructor threw an exception of unknown type" << endl;
      cerr << "Import constructor threw an exception of unknown type" << endl;
      success = false;
    }

    try {
      targetVec.doImport (sourceVec, *importer, Tpetra::INSERT);
    }
    catch (std::exception& e) {
      err << "doImport threw an exception: " << e.what () << endl;
      cerr << "doImport threw an exception: " << e.what () << endl;
      success = false;
    }
    catch (...) {
      err << "doImport threw an exception of unknown type" << endl;
      cerr << "doImport threw an exception of unknown type" << endl;
      success = false;
    }
  } // if I am one of the first 8 processes in MPI_COMM_WORLD

  int lclSuccess = success ? 1 : 0;
  int gblSuccess = 0;
  reduceAll<int, int> (*commWorld, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  success = success && (gblSuccess == 1);

  const int myRank = commWorld->getRank ();
  const int numProcs = commWorld->getSize ();
  if (success) {
    if (myRank == 0) {
      cout << "Success on all processes!" << endl;
    }
  }
  else {
    if (myRank == 0) {
      cout << "At least one process was not successful!" << endl;
    }
    for (int p = 0; p < numProcs; ++p) {
      if (myRank == p) {
        std::ostringstream os;
        os << "Error messages from Process " << myRank << ":" << endl
           << err.str ();
        cout << os.str () << std::flush;
      }
      commWorld->barrier (); // Wait for output to complete
      commWorld->barrier ();
      commWorld->barrier ();
    }
  }

  if (commWorld->getRank () == 0) {
    if (success) {
      cout << "End Result: TEST PASSED" << endl;
    } else {
      cout << "End Result: TEST FAILED" << endl;
    }
  }
  return EXIT_SUCCESS;
}

