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
#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Teuchos_DefaultSerialComm.hpp"

namespace { // (anonymous)

using Teuchos::RCP;
using Teuchos::rcp;
using std::endl;

bool
falseOnSomeProcess (const Teuchos::Comm<int>& comm, const bool in)
{
  using Teuchos::outArg;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;

  const int inInt = in ? 1 : 0;
  int outInt = 0; // output argument
  reduceAll<int, int> (comm, REDUCE_MIN, inInt, outArg (outInt));
  return (outInt == 0);
}

bool
trueOnAllProcesses (const Teuchos::Comm<int>& comm, const bool in)
{
  using Teuchos::outArg;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;

  const int inInt = in ? 1 : 0;
  int outInt = 0; // output argument
  reduceAll<int, int> (comm, REDUCE_MIN, inInt, outArg (outInt));
  return (outInt == 1);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( ImportExport, IsLocallyComplete, LO, GO, NT )
{
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef Tpetra::Import<LO, GO, NT> import_type;
  typedef Tpetra::Export<LO, GO, NT> export_type;

  out << "Test {Ex,Im}port::isLocallyComplete()" << endl;
  Teuchos::OSTab tab1 (out);

  auto comm = Tpetra::getDefaultComm ();
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();
  const GO indexBase = 0;

  //const Tpetra::global_size_t INVALID =
  //  Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid ();

  out << "Test cases where both Maps live on one process" << endl;
  {
    Teuchos::OSTab tab2 (out);
    // This is equivalent to MPI_COMM_SELF.
    auto serialComm = rcp (new Teuchos::SerialComm<int> ());

    RCP<const map_type> map1 (new map_type (10, 10, indexBase, serialComm));
    RCP<const map_type> map2 (new map_type (10, 10, indexBase, serialComm));
    RCP<const map_type> map3 (new map_type (11, 11, indexBase, serialComm));

    Teuchos::Array<GO> gids4 (10);
    typedef typename Teuchos::Array<GO>::size_type size_type;
    for (size_type k = 0; k < static_cast<size_type> (10); ++k) {
      gids4[k] = static_cast<GO> ((10 - 1) - k); // reverse order
    }
    RCP<const map_type> map4 (new map_type (10, gids4 (), indexBase, serialComm));

    out << "Test identical Maps" << endl;
    export_type exp11 (map1, map1);
    TEST_ASSERT( exp11.isLocallyComplete () );
    import_type imp11 (map1, map1);
    TEST_ASSERT( imp11.isLocallyComplete () );

    out << "Test Maps which are the same, but not identical" << endl;
    export_type exp12 (map1, map2);
    TEST_ASSERT( exp12.isLocallyComplete () );
    import_type imp12 (map1, map2);
    TEST_ASSERT( imp12.isLocallyComplete () );

    out << "Test Maps [0, ..., 9], [0, ..., 10]" << endl;
    export_type exp13 (map1, map3);
    // Yes, all _source_ Map indices exist in the _target_ Map.
    TEST_ASSERT( exp13.isLocallyComplete () );
    import_type imp13 (map1, map3);
    TEST_ASSERT( ! imp13.isLocallyComplete () );

    out << "Test Maps [0, ..., 10], [0, ..., 9]" << endl;
    export_type exp31 (map3, map1);
    TEST_ASSERT( ! exp31.isLocallyComplete () );
    import_type imp31 (map3, map1);
    // Yes, all _target_ Map indices exist in the _source_ Map.
    TEST_ASSERT( imp31.isLocallyComplete () );

    out << "Test Maps [0, ..., 9], [9, 8, ..., 0]" << endl;
    export_type exp14 (map1, map4);
    TEST_ASSERT( exp14.isLocallyComplete () );
    import_type imp14 (map1, map4);
    TEST_ASSERT( imp14.isLocallyComplete () );

    out << "Test Maps [9, 8, ..., 0], [0, ..., 9]" << endl;
    export_type exp41 (map4, map1);
    TEST_ASSERT( exp41.isLocallyComplete () );
    import_type imp41 (map4, map1);
    TEST_ASSERT( imp41.isLocallyComplete () );

    out << "Test Maps [9, 8, ..., 0], [0, ..., 10]" << endl;
    export_type exp43 (map4, map3);
    // Yes, all _source_ Map indices exist in the _target_ Map.
    TEST_ASSERT( exp43.isLocallyComplete () );
    import_type imp43 (map4, map3);
    TEST_ASSERT( ! imp43.isLocallyComplete () );

    out << "Test Maps [0, ..., 10], [9, 8, ..., 0]" << endl;
    export_type exp34 (map3, map4);
    TEST_ASSERT( ! exp34.isLocallyComplete () );
    import_type imp34 (map3, map4);
    // Yes, all _target_ Map indices exist in the _source_ Map.
    TEST_ASSERT( imp34.isLocallyComplete () );

    // Make sure that the implementation looks at something other than
    // the global min and max indices.
    out << "Test Maps [0, 2, 3, 6], [0, 1, 4, 6]" << endl;
    {
      Teuchos::Array<GO> gids5 (4);
      gids5[0] = 0;
      gids5[1] = 2;
      gids5[2] = 3;
      gids5[3] = 6;
      Teuchos::Array<GO> gids6 (4);
      gids6[0] = 0;
      gids6[1] = 1;
      gids6[2] = 4;
      gids6[3] = 6;
      RCP<const map_type> map5 (new map_type (4, gids5 (), indexBase, serialComm));
      RCP<const map_type> map6 (new map_type (4, gids6 (), indexBase, serialComm));

      export_type exp56 (map5, map6);
      TEST_ASSERT( ! exp56.isLocallyComplete () );
      import_type imp56 (map5, map6);
      TEST_ASSERT( ! imp56.isLocallyComplete () );
    }
  } // test Maps with MPI_COMM_SELF / SerialComm

  if (numProcs > 1) {
    out << "Test cases where both Maps share the same communicator, "
      "with multiple processes" << endl;
    Teuchos::OSTab tab2 (out);

    out << "Create first three Maps (all contiguous)" << endl;

    RCP<const map_type> map1 (new map_type (10 * numProcs, 10, indexBase, comm));
    RCP<const map_type> map2 (new map_type (10 * numProcs, 10, indexBase, comm));
    RCP<const map_type> map3 (new map_type (11 * numProcs, 11, indexBase, comm));

    out << "Create map4 (noncontiguous, compatible with map1)" << endl;
    RCP<const map_type> map4;
    {
      const LO map1LclSize = static_cast<LO> (map1->getNodeNumElements ());
      TEST_ASSERT( trueOnAllProcesses (*comm, map1LclSize == static_cast<LO> (10)) );
      Teuchos::Array<GO> gids4 (map1LclSize);
      for (LO map1_lid = 0; map1_lid < map1LclSize; ++map1_lid) {
        //const GO map1_gid = map1->getGlobalElement (map1_lid);
        const GO map4_gid = map1->getMaxGlobalIndex () - static_cast<GO> (map1_lid);
        gids4[map1_lid] = map4_gid;
      }
      map4 = rcp (new map_type (map1LclSize * numProcs, gids4 (), indexBase, comm));
    }

    out << "Create map5 (noncontiguous, compatible with map3)" << endl;
    RCP<const map_type> map5;
    {
      const LO map3LclSize = static_cast<LO> (map3->getNodeNumElements ());
      TEST_ASSERT( trueOnAllProcesses (*comm, map3LclSize == static_cast<LO> (11)) );
      Teuchos::Array<GO> gids5 (map3LclSize);
      for (LO map3_lid = 0; map3_lid < map3LclSize; ++map3_lid) {
        //const GO map3_gid = map3->getGlobalElement (map3_lid);
        const GO map5_gid = map3->getMaxGlobalIndex () - static_cast<GO> (map3_lid);
        gids5[map3_lid] = map5_gid;
      }
      map5 = rcp (new map_type (map3LclSize * numProcs, gids5 (), indexBase, comm));
    }

    out << "Create map6 (one entry per process)" << endl;
    Teuchos::Array<GO> gids6 (1);
    gids6[0] = static_cast<GO> (myRank);
    RCP<const map_type> map6 (new map_type (numProcs, gids6 (), 0, comm));

    out << "Create map7 (one entry per process, reverse order of map6)" << endl;
    Teuchos::Array<GO> gids7 (1);
    gids7[0] = static_cast<GO> ((numProcs - 1) - myRank); // reverse order
    RCP<const map_type> map7 (new map_type (numProcs, gids7 (), 0, comm));

    out << "Create map8 (like map6, but with one process (ideally in the "
      "middle, if numProcs > 2) with no indices)" << endl;
    RCP<const map_type> map8;
    {
      Teuchos::ArrayRCP<GO> gids8;
      if (myRank == 0) { // Proc 0 always gets indexBase
        gids8 = Teuchos::arcp<GO> (1);
        gids8[0] = static_cast<GO> (0);
      }
      else if (myRank != 1) { // Proc 1 has no GIDs, always
        gids8 = Teuchos::arcp<GO> (1);
        gids8[0] = static_cast<GO> (myRank);
      }
      map8 = rcp (new map_type ((numProcs - 1), gids8 (), indexBase, comm));
    }

    {
      TEST_ASSERT( map1->isCompatible (*map2) );
      TEST_ASSERT( map1->isSameAs (*map2) );

      export_type exp11 (map1, map1);
      TEST_ASSERT( trueOnAllProcesses (*comm, exp11.isLocallyComplete ()) );
      import_type imp11 (map1, map1);
      TEST_ASSERT( trueOnAllProcesses (*comm, imp11.isLocallyComplete ()) );

      export_type exp12 (map1, map2);
      TEST_ASSERT( trueOnAllProcesses (*comm, exp12.isLocallyComplete ()) );
      import_type imp12 (map1, map2);
      TEST_ASSERT( trueOnAllProcesses (*comm, imp12.isLocallyComplete ()) );
    }

    {
      TEST_ASSERT( ! map1->isCompatible (*map3) );
      TEST_ASSERT( ! map1->isSameAs (*map3) );

      export_type exp13 (map1, map3);
      // Yes, all _source_ Map indices exist in the _target_ Map.
      TEST_ASSERT( trueOnAllProcesses (*comm, exp13.isLocallyComplete ()) );
      import_type imp13 (map1, map3);
      TEST_ASSERT( falseOnSomeProcess (*comm, imp13.isLocallyComplete ()) );

      export_type exp31 (map3, map1);
      TEST_ASSERT( falseOnSomeProcess (*comm, exp31.isLocallyComplete ()) );
      import_type imp31 (map3, map1);
      // Yes, all _target_ Map indices exist in the _source_ Map.
      TEST_ASSERT( trueOnAllProcesses (*comm, imp31.isLocallyComplete ()) );
    }

    {
      TEST_ASSERT( map1->isCompatible (*map4) );
      TEST_ASSERT( ! map1->isSameAs (*map4) );

      export_type exp14 (map1, map4);
      TEST_ASSERT( trueOnAllProcesses (*comm, exp14.getNumRemoteIDs () == 0) );
      TEST_ASSERT( trueOnAllProcesses (*comm, exp14.getNumExportIDs () == 0) );
      TEST_ASSERT( trueOnAllProcesses (*comm, exp14.isLocallyComplete ()) );

      import_type imp14 (map1, map4);
      TEST_ASSERT( trueOnAllProcesses (*comm, imp14.getNumRemoteIDs () == 0) );
      TEST_ASSERT( trueOnAllProcesses (*comm, imp14.getNumExportIDs () == 0) );
      TEST_ASSERT( trueOnAllProcesses (*comm, imp14.isLocallyComplete ()) );

      export_type exp41 (map4, map1);
      TEST_ASSERT( trueOnAllProcesses (*comm, exp41.getNumRemoteIDs () == 0) );
      TEST_ASSERT( trueOnAllProcesses (*comm, exp41.getNumExportIDs () == 0) );
      TEST_ASSERT( trueOnAllProcesses (*comm, exp41.isLocallyComplete ()) );

      import_type imp41 (map4, map1);
      TEST_ASSERT( trueOnAllProcesses (*comm, imp41.getNumRemoteIDs () == 0) );
      TEST_ASSERT( trueOnAllProcesses (*comm, imp41.getNumExportIDs () == 0) );
      TEST_ASSERT( trueOnAllProcesses (*comm, imp41.isLocallyComplete ()) );

      // isOneToOne() may cache communicated information for later
      // use, so if there's a bug in isLocallyComplete(), calling
      // isOneToOne() first may influence the above results.
      TEST_ASSERT( map4->isOneToOne () );
    }

    {
      TEST_ASSERT( ! map1->isCompatible (*map5) );
      TEST_ASSERT( ! map1->isSameAs (*map5) );

      TEST_ASSERT( map3->isCompatible (*map5) );
      TEST_ASSERT( ! map3->isSameAs (*map5) );

      export_type exp15 (map1, map5);
      // Yes, all _source_ Map indices exist in the _target_ Map.
      TEST_ASSERT( trueOnAllProcesses (*comm, exp15.isLocallyComplete ()) );
      import_type imp15 (map1, map5);
      TEST_ASSERT( falseOnSomeProcess (*comm, imp15.isLocallyComplete ()) );

      export_type exp51 (map5, map1);
      TEST_ASSERT( falseOnSomeProcess (*comm, exp51.isLocallyComplete ()) );
      import_type imp51 (map5, map1);
      // Yes, all _target_ Map indices exist in the _source_ Map.
      TEST_ASSERT( trueOnAllProcesses (*comm, imp51.isLocallyComplete ()) );
    }

    {
      export_type exp45 (map4, map5);
      // Yes, all _source_ Map indices exist in the _target_ Map.
      TEST_ASSERT( exp45.isLocallyComplete () );
      import_type imp45 (map4, map5);
      TEST_ASSERT( falseOnSomeProcess (*comm, imp45.isLocallyComplete ()) );

      export_type exp54 (map5, map4);
      // The source Map has some indices not in the target Map on any process.
      TEST_ASSERT( falseOnSomeProcess (*comm, exp54.isLocallyComplete ()) );
      import_type imp54 (map5, map4);
      // Yes, all _target_ Map indices exist in the _source_ Map.
      TEST_ASSERT( trueOnAllProcesses (*comm, imp54.isLocallyComplete ()) );
    }

    {
      export_type exp67 (map6, map7);
      TEST_ASSERT( trueOnAllProcesses (*comm, exp67.isLocallyComplete ()) );
      import_type imp67 (map6, map7);
      TEST_ASSERT( trueOnAllProcesses (*comm, imp67.isLocallyComplete ()) );

      export_type exp76 (map7, map6);
      TEST_ASSERT( trueOnAllProcesses (*comm, exp76.isLocallyComplete ()) );
      import_type imp76 (map7, map6);
      TEST_ASSERT( trueOnAllProcesses (*comm, imp76.isLocallyComplete ()) );

      export_type exp68 (map6, map8);
      // The source Map has some indices not in the target Map on any process.
      TEST_ASSERT( falseOnSomeProcess (*comm, exp68.isLocallyComplete ()) );
      import_type imp68 (map6, map8);
      // Yes, all _target_ Map indices exist in the _source_ Map.
      TEST_ASSERT( trueOnAllProcesses (*comm, imp68.isLocallyComplete ()) );

      export_type exp86 (map8, map6);
      // Yes, all _target_ Map indices exist in the _source_ Map.
      TEST_ASSERT( trueOnAllProcesses (*comm, exp86.isLocallyComplete ()) );
      import_type imp86 (map8, map6);
      // The source Map has some indices not in the target Map on any process.
      TEST_ASSERT( falseOnSomeProcess (*comm, imp86.isLocallyComplete ()) );
    }
  }
}

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( ImportExport, IsLocallyComplete, LO, GO, NT )

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_LGN( UNIT_TEST_GROUP )

} // namespace (anonymous)
