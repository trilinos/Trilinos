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

#include <Tpetra_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Tuple.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_TieBreak.hpp>

#define NUM_GLOBAL_ELEMENTS 100



namespace {

  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef Tpetra::Map<> Map;
  typedef Tpetra::Directory<> Directory;

  typedef Map::local_ordinal_type LO;
  typedef Map::global_ordinal_type GO;
  typedef Map::node_type NT;
  typedef Tpetra::global_size_t GST;

  // If shared, GID goes to rank GID%2.  If not shared (e.g.,
  // pid_and_lid.size()==1), do not change PID assignment.
  template <typename LO, typename GO>
  class ModTwoTieBreak : public Tpetra::Details::TieBreak<LO, GO> {
  public:
    std::size_t
    selectedIndex (GO GID,
                   const std::vector<std::pair<int,LO> > & pid_and_lid) const
    {
      std::size_t index = 0;
      if (pid_and_lid.size() > 1) {
        for (std::size_t i = 0; i < pid_and_lid.size (); ++i) {
          if (pid_and_lid[i].first == (GID % 2)) {
            index = i;
            break;
          }
        }
      }
      // std::cout << "TIE BREAK " << GID
      //           << " goes to " << pid_and_lid[index].first << std::endl;
      return index;
    }
  };


  TEUCHOS_UNIT_TEST(OneToOne, TieBreakAlmostAllOnOne)
  {
    using std::endl;

    // Create a noncontiguous Map with almost all of the GIDs on a
    // single processor.  All processes share GIDs 41, 42, 51, and 52.
    // In the corresponding one-to-one Map, rank 1 should have GIDs 41
    // and 51; rank 0 should have the rest.
    RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();
    const int myRank = comm->getRank ();

    out << "Creating global index lists" << endl;

    Array<GO> gidList;
    size_t num_loc_elems;
    if (myRank == 0) {
      num_loc_elems = NUM_GLOBAL_ELEMENTS;
      gidList = Array<GO> (NUM_GLOBAL_ELEMENTS);
      for (int i = 0; i < NUM_GLOBAL_ELEMENTS; ++i) {
        gidList[i] = i;
      }
    }
    else {
      num_loc_elems = 4;
      gidList = Array<GO> (num_loc_elems);
      gidList[0] = 51; // These two GIDs are in rank one's directory
      gidList[1] = 52;
      gidList[2] = 41; // These two GIDs are in rank zero's directory
      gidList[3] = 42;
    }

    out << "Building non-one-to-one Map" << endl;
    const GO indexBase = 0;
    RCP<const Map> map =
      rcp (new Map (Teuchos::OrdinalTraits<GST>::invalid (), gidList (),
                    indexBase, comm));

    out << "Calling createOneToOne" << endl;
    ModTwoTieBreak<LO,GO> tie_break;
    RCP<const Map> new_map = Tpetra::createOneToOne<LO,GO,NT> (map, tie_break);

    out << "Print the new map" << endl;
    // The "Teuchos::" stuff turns std::cout into a FancyOStream.
    new_map->describe (* (Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cout))),
                       Teuchos::VERB_EXTREME);

    out << "Checking results" << endl;

    ArrayView<const GO> my_owned = new_map->getNodeElementList ();
    const LO my_num_owned = static_cast<LO> (new_map->getNodeNumElements ());
    Array<int> nodeIDlist (num_loc_elems);
    Tpetra::LookupStatus stat =
      new_map->getRemoteIndexList (gidList (), nodeIDlist ());

    //If we pass this we didn't lose any IDs.
    TEST_EQUALITY_CONST(stat, Tpetra::AllIDsPresent);
    TEST_EQUALITY(new_map->getGlobalNumElements (), NUM_GLOBAL_ELEMENTS);

    // Make sure global indices are in the right place.  All elements
    // except 41 and 51 should be on rank 0.  Global indices 41 and 51
    // should be on rank 1.

    if (myRank == 0) {
      for (LO i = 0; i < my_num_owned; ++i) {
        TEST_INEQUALITY(my_owned[i], 41);
        TEST_INEQUALITY(my_owned[i], 51);
      }
    }
    else {
      TEST_EQUALITY(my_num_owned, 2);

      if (my_num_owned == static_cast<LO> (2)) {
        // Make a deep copy, so we can sort.  Order of the indices
        // doesn't matter; all that matters is that Proc 1 actually
        // got the indices it was supposed to get.
        Array<GO> myOwnedCopy (my_owned.begin (), my_owned.end ());
        std::sort (myOwnedCopy.begin (), myOwnedCopy.end ());
        TEST_EQUALITY(myOwnedCopy[0], 41);
        TEST_EQUALITY(myOwnedCopy[1], 51);
      }
    }

    comm->barrier ();
    out << "Done" << endl;
  }
}


