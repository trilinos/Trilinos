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

// Some Macro Magic to ensure that if CUDA and KokkosCompat is enabled
// only the .cu version of this file is actually compiled
#include <Tpetra_config.h>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Tpetra_Map.hpp>

#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION
#  include "Tpetra_Map_def.hpp"
#  include "Tpetra_Directory_def.hpp"
#  ifdef HAVE_TPETRA_FIXED_HASH_TABLE
#    include "Tpetra_Details_FixedHashTable_def.hpp"
#  else
#    include "Tpetra_HashTable_def.hpp"
#  endif // HAVE_TPETRA_FIXED_HASH_TABLE
#endif

using Tpetra::global_size_t;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::as;
using Teuchos::Comm;
using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::tuple;
using std::endl;

// This test works (and exercises the interesting case) in serial mode
// or for 1 MPI process, but it was originally written for 2 MPI
// processes.
TEUCHOS_UNIT_TEST( Map, Bug5822_StartWith3Billion )
{
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm ();
  const int numProcs = comm->getSize ();
  const int myRank = comm->getRank ();
  TEUCHOS_TEST_FOR_EXCEPTION(numProcs != 2, std::logic_error,
    "This test only makes sense to run with 2 MPI processes.");

#ifdef HAVE_TEUCHOS_LONG_LONG_INT
  typedef long long GO;
  if (sizeof (long long) <= 4) {
    out << "sizeof (long long) = " << sizeof (long long) << " <= 4.  "
      "This test only makes sense if sizeof (long long) >= 8, "
      "since the test is supposed to exercise GIDs > 2 billion.";
    return;
  }
#else // NOT HAVE_TEUCHOS_LONG_LONG_INT
  typedef long GO;
  if (sizeof (long) <= 4) {
    out << "sizeof (long) = " << sizeof (long) << " <= 4.  "
      "This test only makes sense if sizeof (long) >= 8, "
      "since the test is supposed to exercise GIDs > 2 billion.";
    return;
  }
#endif // HAVE_TEUCHOS_LONG_LONG_INT
  typedef int LO;
  typedef KokkosClassic::DefaultNode::DefaultNodeType NT;
  typedef Tpetra::Map<LO, GO, NT> map_type;

  const size_t localNumElts = 5;
  const global_size_t globalNumElts = comm->getSize () * localNumElts;
  const GO globalFirstGid = 3000000000L;
  const GO localFirstGid = globalFirstGid + 2 * localNumElts * comm->getRank ();

  // Make sure that the GIDs are not contiguous, in case Map tries to
  // detect contiguous GIDs for us.  The point is to exercise a
  // noncontiguous Map with GIDs all greater than 2 billion (so that
  // they would overflow a 32-bit integer).
  Array<GO> myGids (localNumElts);
  Array<GO> myGidsExpected (localNumElts);
  for (size_t k = 0; k < localNumElts; ++k) {
    myGids[k] = localFirstGid + as<GO> (2*k);
    // Make a copy, just to make sure that Map's constructor didn't
    // sneakily change the input ArrayView.
    myGidsExpected[k] = myGids[k];
  }
  // Proc 0: myGids = [3B, 3B+2, 3B+4, 3B+6, 3B+8].
  // Proc 1: myGids = [3B+10, 3B+12, 3B+14, 3B+16, 3B+18].

  RCP<NT> node;
  {
    Teuchos::ParameterList defaultParams;
    node = Teuchos::rcp (new NT (defaultParams));
  }
  // Tpetra::Map requires that the index base and the first GID be equal.
  const GO indexBase = globalFirstGid;
  RCP<const map_type> map (new map_type (globalNumElts, myGids (), indexBase, comm, node));

  ArrayView<const GO> myGidsFound = map->getNodeElementList ();
  TEST_COMPARE_ARRAYS( myGidsExpected (), myGidsFound () );

  Array<GO> remoteGids (5);
  Array<int> remotePids (5, -1);
  Array<LO> remoteLids (5, Teuchos::OrdinalTraits<LO>::invalid ());
  if (myRank == 0) {
    Array<int> expectedRemotePids (5);
    std::fill (expectedRemotePids.begin (), expectedRemotePids.end (), 1);
    Array<int> expectedRemoteLids (5);
    expectedRemoteLids[0] = 0;
    expectedRemoteLids[1] = 1;
    expectedRemoteLids[2] = 2;
    expectedRemoteLids[3] = 3;
    expectedRemoteLids[4] = 4;
    remoteGids[0] = globalFirstGid + 10;
    remoteGids[1] = globalFirstGid + 12;
    remoteGids[2] = globalFirstGid + 14;
    remoteGids[3] = globalFirstGid + 16;
    remoteGids[4] = globalFirstGid + 18;
    map->getRemoteIndexList (remoteGids (), remotePids (), remoteLids ());

    TEST_COMPARE_ARRAYS( remotePids (), expectedRemotePids () );
    TEST_COMPARE_ARRAYS( remoteLids (), expectedRemoteLids () );
  } else if (myRank == 1) {
    Array<int> expectedRemotePids (5);
    std::fill (expectedRemotePids.begin (), expectedRemotePids.end (), 0);
    Array<int> expectedRemoteLids (5);
    expectedRemoteLids[0] = 0;
    expectedRemoteLids[1] = 1;
    expectedRemoteLids[2] = 2;
    expectedRemoteLids[3] = 3;
    expectedRemoteLids[4] = 4;
    remoteGids[0] = globalFirstGid;
    remoteGids[1] = globalFirstGid + 2;
    remoteGids[2] = globalFirstGid + 4;
    remoteGids[3] = globalFirstGid + 6;
    remoteGids[4] = globalFirstGid + 8;
    map->getRemoteIndexList (remoteGids (), remotePids (), remoteLids ());

    TEST_COMPARE_ARRAYS( remotePids (), expectedRemotePids () );
    TEST_COMPARE_ARRAYS( remoteLids (), expectedRemoteLids () );
  }
}


