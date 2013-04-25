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
#include <Teuchos_DefaultComm.hpp>
#include <Tpetra_Map.hpp>
#include <Teuchos_TimeMonitor.hpp>

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
using std::cout;
using std::cin;

// This test works (and exercises the interesting case) in serial mode
// or for 1 MPI process, but it was originally written for 2 MPI
// processes.
TEUCHOS_UNIT_TEST( Map, Bug5822_StartWithZeroThenSkipTo3Billion )
{
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
  typedef int LO;
  typedef long long GO;
  typedef Kokkos::SerialNode NT;
  typedef Tpetra::Map<LO, GO, NT> map_type;

  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm ();
  const int myRank = comm->getRank ();

  // Proc 0 gets [0, 3B, 3B+2, 3B+4, 3B+8, 3B+10] (6 GIDs).
  // Proc 1 gets [3B+12, 3B+14, 3B+16, 3B+18, 3B+20] (5 GIDs).
  //
  // The GIDs are not contiguous in order to prevent Map from
  // detecting contiguity and bypassing the noncontiguous Map
  // case.
  const size_t localNumElts = (myRank == 0) ? 6 : 5;
  const global_size_t globalNumElts = 1 + 5*comm->getSize ();
  const GO globalFirstGid = 1L;
  const GO threeBillion = 3000000000L;

  Array<GO> myGids (localNumElts);
  // Make a copy, just to make sure that Map's constructor didn't
  // sneakily change the input ArrayView.
  Array<GO> myGidsExpected (localNumElts);
  if (myRank == 0) {
    myGids[0] = globalFirstGid;
    myGidsExpected[0] = myGids[0];
    myGids[1] = threeBillion;
    myGidsExpected[1] = myGids[1];
    for (size_t k = 2; k < localNumElts; ++k) {
      myGids[k] = myGids[k-1] + 2;
      myGidsExpected[k] = myGids[k];
    }
  } else {
    myGids[0] = threeBillion + as<GO> ((localNumElts+1) * 2);
    myGidsExpected[0] = myGids[0];
    for (size_t k = 1; k < localNumElts; ++k) {
      myGids[k] = myGids[k-1] + 2;
      myGidsExpected[k] = myGids[k];
    }
  }

  RCP<NT> node;
  {
    ParameterList junk;
    node = rcp (new NT (junk));
  }
  // if (myRank == 0) {
  //   std::cout << "type '0' and hit enter" << std::endl;
  //   int zero;
  //   cin >> zero;
  // }
  // comm->barrier();
  // Tpetra::Map requires that the index base is less than the minimum GID.
  const GO indexBase = 0L;
  out << "Constructing the map..." << endl;
  RCP<const map_type> map;
  {
    TEUCHOS_FUNC_TIME_MONITOR("Construct Map");
    map = rcp(new map_type (globalNumElts, myGids (), indexBase, comm, node));
  }

  {
    TEUCHOS_FUNC_TIME_MONITOR("Querying the map for local elements");
    ArrayView<const GO> myGidsFound = map->getNodeElementList ();
    TEST_COMPARE_ARRAYS( myGidsExpected (), myGidsFound () );
  }
  Teuchos::TimeMonitor::summarize();
#else
  out << "This test only makes sense if long long support is enabled in Teuchos.";
  return;
#endif // HAVE_TEUCHOS_LONG_LONG_INT
}




