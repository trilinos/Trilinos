// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>

#include "Teuchos_as.hpp"

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_DefaultPlatform.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "Tpetra_ConfigDefs.hpp" //TODO
#include "Tpetra_DefaultPlatform.hpp" //TODO
#include "Xpetra_StridedTpetraMap.hpp"
#endif

#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_StridedEpetraMap.hpp"
#endif

#include "Xpetra_StridedMapFactory.hpp"

namespace {
#ifdef HAVE_XPETRA_EPETRA
  typedef Xpetra::StridedEpetraMap StridedEpetraMap;
#endif
  using Teuchos::Array;
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::arcp;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Teuchos::Tuple;
  using Teuchos::tuple;
  using Xpetra::global_size_t;
  using Xpetra::DefaultPlatform;
  using std::sort;
  using std::find;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::Comm;

  bool testMpi = true;
  double errorTolSlack = 1e+1;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignored and a serial comm is always used." );
    clp.setOption(
        "error-tol-slack", &errorTolSlack,
        "Slack off of machine epsilon used to check test results" );
  }

  RCP<const Comm<int> > getDefaultComm()
  {
    if (testMpi) {
      return DefaultPlatform::getDefaultPlatform().getComm();
    }
    return rcp(new Teuchos::SerialComm<int>());
  }

  //
  // UNIT TESTS
  //


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( StridedMapFactory, CreateStridedEpetraMap1, LO, GO, Node )
  {
    // create a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    //const int myImageID = comm->getRank();

#ifdef HAVE_XPETRA_EPETRA
    for (int indexBase = 0; indexBase < 2; indexBase++) {

      GO offset = 111;

      // constructor calls: (num global elements, index base)
      global_size_t numGlobalElements = 120 * numImages;
      size_t numLocalElements = 120;
      std::vector<size_t> stridedInfo;
      stridedInfo.push_back(3);
      stridedInfo.push_back(4);
      stridedInfo.push_back(5);

      Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;

      Teuchos::RCP<Xpetra::StridedMap<LO,GO,Node> > map = Xpetra::StridedMapFactory<LO,GO,Node>::Build(lib, numGlobalElements, indexBase, stridedInfo, comm, -1, offset);

      TEST_EQUALITY_CONST( map->getFixedBlockSize(), 12 );
      TEST_EQUALITY_CONST( map->isStrided(), true );
      TEST_EQUALITY_CONST( map->isBlocked(), true );
      TEST_EQUALITY_CONST( map->getMinAllGlobalIndex(), indexBase + offset );
      TEST_EQUALITY_CONST( map->getMaxAllGlobalIndex(), indexBase + offset + Teuchos::as<GO>(numGlobalElements) - 1);
      TEST_EQUALITY_CONST( map->isContiguous(), false);
      TEST_EQUALITY_CONST( map->getNodeNumElements() % 12 , 0);

      Teuchos::RCP<Xpetra::StridedEpetraMap> emap = Teuchos::rcp_dynamic_cast<Xpetra::StridedEpetraMap>(map);
      TEST_EQUALITY_CONST( emap != Teuchos::null , true);

      Teuchos::RCP<Xpetra::StridedEpetraMap> emap2 = Teuchos::null;
      for(size_t k=0; k<stridedInfo.size(); k++) {
        Teuchos::RCP<Xpetra::StridedMap<LO,GO,Node> > map2 = Xpetra::StridedMapFactory<LO,GO,Node>::Build(lib, numGlobalElements, indexBase, stridedInfo, comm, k, offset);
        TEST_EQUALITY_CONST( map2->getFixedBlockSize(), 12 );
        TEST_EQUALITY_CONST( map2->isStrided(), true );
        TEST_EQUALITY_CONST( map2->isBlocked(), true );
        TEST_EQUALITY_CONST( map2->isContiguous(), false);
        TEST_EQUALITY_CONST( map2->getNodeNumElements() % stridedInfo[k] , 0);
        TEST_EQUALITY_CONST( map2->getNodeNumElements(), numLocalElements / map2->getFixedBlockSize() * stridedInfo[k] );
        emap2 = Teuchos::rcp_dynamic_cast<Xpetra::StridedEpetraMap>(map2);
        TEST_EQUALITY_CONST( emap2 != Teuchos::null , true);
      }

    }
#endif
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( StridedMapFactory, CreateStridedTpetraMap1, LO, GO, Node )
  {
    // create a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    //const int myImageID = comm->getRank();

#ifdef HAVE_XPETRA_TPETRA
    for (int indexBase = 0; indexBase < 2; indexBase++) {
      GO offset = 111;

      // constructor calls: (num global elements, index base)
      global_size_t numGlobalElements = 120 * numImages;
      size_t numLocalElements = 120;
      std::vector<size_t> stridedInfo;
      stridedInfo.push_back(3);
      stridedInfo.push_back(4);
      stridedInfo.push_back(5);

      Xpetra::UnderlyingLib lib = Xpetra::UseTpetra;

      Teuchos::RCP<Xpetra::StridedMap<LO,GO,Node> > map = Xpetra::StridedMapFactory<LO,GO,Node>::Build(lib, numGlobalElements, indexBase, stridedInfo, comm, -1, offset);

      TEST_EQUALITY_CONST( map->getFixedBlockSize(), 12 );
      TEST_EQUALITY_CONST( map->isStrided(), true );
      TEST_EQUALITY_CONST( map->isBlocked(), true );
      TEST_EQUALITY_CONST( map->getMinAllGlobalIndex(), indexBase + offset );
      TEST_EQUALITY_CONST( map->getMaxAllGlobalIndex(), indexBase + offset + Teuchos::as<GO>(numGlobalElements) - 1);
      TEST_EQUALITY_CONST( map->isContiguous(), false);
      TEST_EQUALITY_CONST( map->getNodeNumElements() % 12 , 0);

      Teuchos::RCP<Xpetra::StridedTpetraMap<LO,GO,Node> > tmap = Teuchos::rcp_dynamic_cast<Xpetra::StridedTpetraMap<LO,GO,Node> >(map);
      TEST_EQUALITY_CONST( tmap != Teuchos::null , true);

      Teuchos::RCP<Xpetra::StridedTpetraMap<LO,GO,Node> > tmap2 = Teuchos::null;
      for(size_t k=0; k<stridedInfo.size(); k++) {
        Teuchos::RCP<Xpetra::StridedMap<LO,GO,Node> > map2 = Xpetra::StridedMapFactory<LO,GO,Node>::Build(lib, numGlobalElements, indexBase, stridedInfo, comm, k, offset);
        TEST_EQUALITY_CONST( map2->getFixedBlockSize(), 12 );
        TEST_EQUALITY_CONST( map2->isStrided(), true );
        TEST_EQUALITY_CONST( map2->isBlocked(), true );
        TEST_EQUALITY_CONST( map2->isContiguous(), false);
        TEST_EQUALITY_CONST( map2->getNodeNumElements() % stridedInfo[k] , 0);
        TEST_EQUALITY_CONST( map2->getNodeNumElements(), numLocalElements / map2->getFixedBlockSize() * stridedInfo[k] );
        tmap2 = Teuchos::rcp_dynamic_cast<Xpetra::StridedTpetraMap<LO,GO,Node> >(map);
        TEST_EQUALITY_CONST( tmap2 != Teuchos::null , true);
      }
    }
#endif
  }

  // TODO add test routines for remaining constructors of StridedMapFactory

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( StridedMapFactory, CreateStridedEpetraMap2, LO, GO, Node )
  {
    // create a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    //const int myImageID = comm->getRank();

#ifdef HAVE_XPETRA_EPETRA
    GO offset = 111;

    // constructor calls: (num global elements, index base)
    global_size_t numGlobalElements = 120 * numImages;
    size_t numLocalElements = 120;
    std::vector<size_t> stridedInfo;
    stridedInfo.push_back(3);
    stridedInfo.push_back(4);
    stridedInfo.push_back(5);

    Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;

    Teuchos::RCP<Xpetra::StridedMap<LO,GO,Node> > map = Xpetra::StridedMapFactory<LO,GO,Node>::Build(lib, numGlobalElements, 0, stridedInfo, comm, -1, offset);

    TEST_EQUALITY_CONST( map->getFixedBlockSize(), 12 );
    TEST_EQUALITY_CONST( map->isStrided(), true );
    TEST_EQUALITY_CONST( map->isBlocked(), true );
    TEST_EQUALITY_CONST( map->getMinAllGlobalIndex(), offset );
    TEST_EQUALITY_CONST( map->getMaxAllGlobalIndex(), offset + Teuchos::as<GO>(numGlobalElements) - 1);
    TEST_EQUALITY_CONST( map->isContiguous(), false);
    TEST_EQUALITY_CONST( map->getNodeNumElements() % 12 , 0);

    Teuchos::RCP<Xpetra::StridedMap<LO,GO,Node> > map2 = Xpetra::StridedMapFactory<LO,GO,Node>::Build(map, 0);
    TEST_EQUALITY_CONST( map2->getFixedBlockSize(), 12 );
    TEST_EQUALITY_CONST( map2->isStrided(), true );
    TEST_EQUALITY_CONST( map2->isBlocked(), true );
    TEST_EQUALITY_CONST( map2->getMinAllGlobalIndex(), offset );
    TEST_EQUALITY_CONST( map2->getMaxAllGlobalIndex(), offset + Teuchos::as<GO>(numGlobalElements) - 10 );
    TEST_EQUALITY_CONST( map2->isContiguous(), false);
    TEST_EQUALITY_CONST( map2->getNodeNumElements() % 3 , 0);
    TEST_EQUALITY_CONST( map2->getNodeNumElements(), numLocalElements / map2->getFixedBlockSize() * stridedInfo[0]);

    Teuchos::RCP<Xpetra::StridedMap<LO,GO,Node> > map3 = Xpetra::StridedMapFactory<LO,GO,Node>::Build(map, 1);
    TEST_EQUALITY_CONST( map3->getFixedBlockSize(), 12 );
    TEST_EQUALITY_CONST( map3->isStrided(), true );
    TEST_EQUALITY_CONST( map3->isBlocked(), true );
    TEST_EQUALITY_CONST( map3->getMinAllGlobalIndex(), offset + 3 );
    TEST_EQUALITY_CONST( map3->getMaxAllGlobalIndex(), offset + Teuchos::as<GO>(numGlobalElements) - 6 );
    TEST_EQUALITY_CONST( map3->isContiguous(), false);
    TEST_EQUALITY_CONST( map3->getNodeNumElements() % 4 , 0);
    TEST_EQUALITY_CONST( map3->getNodeNumElements(), numLocalElements / map3->getFixedBlockSize() * stridedInfo[1]);

    Teuchos::RCP<Xpetra::StridedMap<LO,GO,Node> > map4 = Xpetra::StridedMapFactory<LO,GO,Node>::Build(map, 2);
    TEST_EQUALITY_CONST( map4->getFixedBlockSize(), 12 );
    TEST_EQUALITY_CONST( map4->isStrided(), true );
    TEST_EQUALITY_CONST( map4->isBlocked(), true );
    TEST_EQUALITY_CONST( map4->getMinAllGlobalIndex(), offset + 7 );
    TEST_EQUALITY_CONST( map4->getMaxAllGlobalIndex(), offset + Teuchos::as<GO>(numGlobalElements) - 1 );
    TEST_EQUALITY_CONST( map4->isContiguous(), false);
    TEST_EQUALITY_CONST( map4->getNodeNumElements() % 5 , 0);
    TEST_EQUALITY_CONST( map4->getNodeNumElements(), numLocalElements / map4->getFixedBlockSize() * stridedInfo[2]);
#endif
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( StridedMapFactory, CreateStridedTpetraMap2, LO, GO, Node )
  {
    // create a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    //const int myImageID = comm->getRank();

#ifdef HAVE_XPETRA_TPETRA
    GO offset = 111;

    // constructor calls: (num global elements, index base)
    global_size_t numGlobalElements = 120 * numImages;
    size_t numLocalElements = 120;
    std::vector<size_t> stridedInfo;
    stridedInfo.push_back(3);
    stridedInfo.push_back(4);
    stridedInfo.push_back(5);

    Xpetra::UnderlyingLib lib = Xpetra::UseTpetra;

    Teuchos::RCP<Xpetra::StridedMap<LO,GO,Node> > map = Xpetra::StridedMapFactory<LO,GO,Node>::Build(lib, numGlobalElements, 0, stridedInfo, comm, -1, offset);

    TEST_EQUALITY_CONST( map->getFixedBlockSize(), 12 );
    TEST_EQUALITY_CONST( map->isStrided(), true );
    TEST_EQUALITY_CONST( map->isBlocked(), true );
    TEST_EQUALITY_CONST( map->getMinAllGlobalIndex(), offset );
    TEST_EQUALITY_CONST( map->getMaxAllGlobalIndex(), offset + Teuchos::as<GO>(numGlobalElements) - 1);
    TEST_EQUALITY_CONST( map->isContiguous(), false);
    TEST_EQUALITY_CONST( map->getNodeNumElements() % 12 , 0);

    Teuchos::RCP<Xpetra::StridedMap<LO,GO,Node> > map2 = Xpetra::StridedMapFactory<LO,GO,Node>::Build(map, 0);
    TEST_EQUALITY_CONST( map2->getFixedBlockSize(), 12 );
    TEST_EQUALITY_CONST( map2->isStrided(), true );
    TEST_EQUALITY_CONST( map2->isBlocked(), true );
    TEST_EQUALITY_CONST( map2->getMinAllGlobalIndex(), offset );
    TEST_EQUALITY_CONST( map2->getMaxAllGlobalIndex(), offset + Teuchos::as<GO>(numGlobalElements) - 10 );
    TEST_EQUALITY_CONST( map2->isContiguous(), false);
    TEST_EQUALITY_CONST( map2->getNodeNumElements() % 3 , 0);
    TEST_EQUALITY_CONST( map2->getNodeNumElements(), numLocalElements / map2->getFixedBlockSize() * stridedInfo[0]);

    Teuchos::RCP<Xpetra::StridedMap<LO,GO,Node> > map3 = Xpetra::StridedMapFactory<LO,GO,Node>::Build(map, 1);
    TEST_EQUALITY_CONST( map3->getFixedBlockSize(), 12 );
    TEST_EQUALITY_CONST( map3->isStrided(), true );
    TEST_EQUALITY_CONST( map3->isBlocked(), true );
    TEST_EQUALITY_CONST( map3->getMinAllGlobalIndex(), offset + 3 );
    TEST_EQUALITY_CONST( map3->getMaxAllGlobalIndex(), offset + Teuchos::as<GO>(numGlobalElements) - 6 );
    TEST_EQUALITY_CONST( map3->isContiguous(), false);
    TEST_EQUALITY_CONST( map3->getNodeNumElements() % 4 , 0);
    TEST_EQUALITY_CONST( map3->getNodeNumElements(), numLocalElements / map3->getFixedBlockSize() * stridedInfo[1]);

    Teuchos::RCP<Xpetra::StridedMap<LO,GO,Node> > map4 = Xpetra::StridedMapFactory<LO,GO,Node>::Build(map, 2);
    TEST_EQUALITY_CONST( map4->getFixedBlockSize(), 12 );
    TEST_EQUALITY_CONST( map4->isStrided(), true );
    TEST_EQUALITY_CONST( map4->isBlocked(), true );
    TEST_EQUALITY_CONST( map4->getMinAllGlobalIndex(), offset + 7 );
    TEST_EQUALITY_CONST( map4->getMaxAllGlobalIndex(), offset + Teuchos::as<GO>(numGlobalElements) - 1 );
    TEST_EQUALITY_CONST( map4->isContiguous(), false);
    TEST_EQUALITY_CONST( map4->getNodeNumElements() % 5 , 0);
    TEST_EQUALITY_CONST( map4->getNodeNumElements(), numLocalElements / map4->getFixedBlockSize() * stridedInfo[2]);
#endif
  }
  //
  // INSTANTIATIONS
  //

#define UNIT_TEST_GROUP_ORDINAL( LO, GO, Node )                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( StridedMapFactory, CreateStridedEpetraMap1, LO, GO, Node ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( StridedMapFactory, CreateStridedTpetraMap1, LO, GO, Node ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( StridedMapFactory, CreateStridedEpetraMap2, LO, GO, Node ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( StridedMapFactory, CreateStridedTpetraMap2, LO, GO, Node )

  typedef Kokkos::DefaultNode::DefaultNodeType DefaultNodeType;

  UNIT_TEST_GROUP_ORDINAL(int, int, DefaultNodeType)

}



