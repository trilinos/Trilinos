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

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_Tuple.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"

using Teuchos::RCP;
using Teuchos::Array;
using Teuchos::tuple;

////
TEUCHOS_UNIT_TEST( Map, ProblematicLookup )
{
  using std::cerr;
  using std::endl;

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const int myRank = comm->getRank();
  /**********************************************************************************/
  // Map in question:
  // -----------------------------
  // SRC Map  Processor 0: Global IDs = 0 1 2 3 4 5 6
  //          Processor 1: Global IDs =                    9 10 11 12 13 14 15
  //
  // Lookup of global IDs 7 8 should return IDNotFound
  //
  if (myRank == 0) {
    cerr << "Creating Map" << endl;
  }
  comm->barrier ();
  comm->barrier ();
  comm->barrier (); // Just to make sure output finishes.

  typedef Tpetra::Map<> map_type;
  typedef map_type::local_ordinal_type LO;
  typedef map_type::global_ordinal_type GO;
  typedef int rank_type;

  RCP<const map_type> map;
  if (myRank == 0) {
    Array<GO> gids (tuple<GO> (1));
    map = Tpetra::createNonContigMap<LO, GO> (gids ().getConst () , comm);
  }
  else {
    Array<GO> gids (tuple<GO> (3));
    map = Tpetra::createNonContigMap<LO, GO> (gids ().getConst (), comm);
  }

  {
    std::ostringstream os;
    os << "Proc " << myRank << ": created Map" << endl;
    cerr << os.str ();
  }

  {
    std::ostringstream os;
    os << "Proc " << myRank << ": calling getRemoteIndexList" << endl;
    cerr << os.str ();
  }

  Array<rank_type> processRanks (1);
  Tpetra::LookupStatus lookup = map->getRemoteIndexList (tuple<GO> (2), processRanks ());

  {
    std::ostringstream os;
    os << "Proc " << myRank << ": getRemoteIndexList done" << endl;
    cerr << os.str ();
  }
  comm->barrier ();
  if (myRank == 0) {
    cerr << "getRemoteIndexList finished on all processes" << endl;
  }
  comm->barrier (); // Just to make sure output finishes.

  TEST_EQUALITY_CONST( map->isDistributed(), true )
  TEST_EQUALITY_CONST( map->isContiguous(), false )
  TEST_EQUALITY_CONST( lookup, Tpetra::IDNotPresent )
  TEST_COMPARE_ARRAYS( processRanks(), tuple<rank_type>(-1) );
}

