// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Thanks to Karen Devine for contributing this test.  We made some
// small changes to make it use the Teuchos unit test framework.
//
// Bug 6412 says that Tpetra::Map::getRemoteIndexList does not
// correctly handle repeated global indices on the same process, but
// only when the Map is noncontiguous.  Zoltan's DDirectory _does_
// correctly handle this case.

#include "Tpetra_Map.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_UnitTestHarness.hpp"

namespace { // (anonymous)

using Teuchos::arcp;
typedef Tpetra::Map<> map_type;
typedef map_type::local_ordinal_type LO;
typedef map_type::global_ordinal_type GO;

void
searchIt (bool& success, std::ostream& out, const map_type& myMap, const std::string& myName)
{
  using std::endl;
  const int me = myMap.getComm ()->getRank ();

  // Print the map elements
  out << me << " " << myName << " MINE: ";
  for (size_t i = 0; i < myMap.getLocalNumElements (); ++i) {
    out << myMap.getGlobalElement(i) << " ";
  }
  out << endl;

  // Memory for Gids for which to search
  size_t nSearch = 6;
  Teuchos::ArrayRCP<GO> searchGids = arcp (new GO[nSearch], 0, nSearch, true);
  Teuchos::ArrayRCP<int> searchRemoteRanks = arcp (new int[nSearch], 0, nSearch, true);
  Teuchos::ArrayRCP<LO> searchRemoteLids = arcp (new LO[nSearch], 0, nSearch, true);

  // Search without duplicates
  for (size_t i = 0; i < nSearch; ++i) {
    searchGids[i] = i;
  }
  myMap.getRemoteIndexList (searchGids (),
                            searchRemoteRanks (),
                            searchRemoteLids ());

  for (size_t i = 0; i < nSearch; ++i) {
    if (searchRemoteRanks[i] == -1) {
      success = false;
    }
    out << me << " " << myName
        << " NoDuplicates:  GID " << searchGids[i]
        << " RANK " << searchRemoteRanks[i]
        << " LID " << searchRemoteLids[i]
        << (searchRemoteRanks[i] == -1 ? "  BAD!" : " ")
        << endl;
  }

  // Search with duplicates in the input list of global indices.  We
  // force duplicates to show up by dividing by 2, so that 2k and 2k+1
  // both become the same global index k.
  for (size_t i = 0; i < nSearch; i++) {
    searchGids[i] = i/2;
  }
  myMap.getRemoteIndexList (searchGids (),
                            searchRemoteRanks (),
                            searchRemoteLids ());

  for (size_t i = 0; i < nSearch; ++i) {
    if (searchRemoteRanks[i] == -1) {
      success = false;
    }
    out << me << " " << myName
        << " WithDuplicates:  GID " << searchGids[i]
        << " RANK " << searchRemoteRanks[i]
        << " LID " << searchRemoteLids[i]
        << (searchRemoteRanks[i] == -1 ? "  BAD!" : " ")
        << endl;
  }
}


TEUCHOS_UNIT_TEST( Map, Bug6412 )
{
  using Teuchos::Comm;
  using Teuchos::RCP;
  typedef Tpetra::global_size_t GST;

  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm ();
  const int me = comm->getRank ();
  const int np = comm->getSize ();

  GO nGlobal = 24;   // Global number of Gids

  // Create and search Default Tpetra Map
  const map_type defaultMap (nGlobal, 0, comm);

  searchIt (success, std::cout, defaultMap, "defaultMap");

  // Create and search customized map
  // Identify locally owned GIDs:  same as default map (if nGlobal%np == 0)
  LO nLocal = nGlobal / np + (me < (nGlobal%np));
  GO myFirst = me * (nGlobal / np) + (me < (nGlobal%np) ? me : (nGlobal%np));
  Teuchos::ArrayRCP<GO> myGids = arcp (new GO[nLocal], 0, nLocal, true);
  for (LO i = 0; i < nLocal; ++i) {
    myGids[i] = myFirst + i;
  }

  // Construct customMap
  const GST dummy = Teuchos::OrdinalTraits<GST>::invalid ();
  const map_type customMap (dummy, myGids (), 0, comm);

  searchIt (success, out, customMap, "customMap");
}

} // namespace (anonymous)
