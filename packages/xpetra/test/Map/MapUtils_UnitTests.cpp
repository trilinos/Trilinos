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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Teuchos_Array.hpp>
#include <Teuchos_as.hpp>

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_DefaultPlatform.hpp"

#include "Xpetra_Exceptions.hpp"
#include "Xpetra_Map.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_MapUtils.hpp"
#include "Xpetra_StridedMap.hpp"
#include "Xpetra_StridedMapFactory.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraMap.hpp"
#include "Tpetra_Details_Behavior.hpp"
#endif

#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraMap.hpp"
#endif

// FINISH: add testing of operator==, operator!=, operator=, copy construct
// put these into test_same_as and test_is_compatible

namespace {

  bool testMpi = true;
  double errorTolSlack = 1e+1;

  Teuchos::RCP<const Teuchos::Comm<int> > getDefaultComm()
  {
    if (testMpi) {
      return Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
    }
    return Teuchos::rcp(new Teuchos::SerialComm<int>());
  }

  //
  // UNIT TESTS
  //
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MapUtils, ConcatenateTwoMaps_ObtainMergedMap, M, LO, GO, N )
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::RCP;

    using Map = Xpetra::Map<LO,GO,N>;
    using MapFactory = Xpetra::MapFactory<LO,GO,N>;
    using MapUtils = Xpetra::MapUtils<LO,GO,N>;

    // create a comm
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    M testMap(1, 0, comm);
    Xpetra::UnderlyingLib lib = testMap.lib();

    const GO GO_zero = Teuchos::OrdinalTraits<GO>::zero();

    const LO numLocalElementsPerMap = 3;
    const Xpetra::global_size_t numGlobalElementsPerMap =
        Teuchos::as<Xpetra::global_size_t>(comm->getSize() * numLocalElementsPerMap);

    Array<GO> gidsForMapOne;
    Array<GO> gidsForMapTwo;
    for (LO i = 0; i < numLocalElementsPerMap; ++i) {
      gidsForMapOne.push_back(Teuchos::as<GO>(comm->getRank() * numLocalElementsPerMap + i));
      gidsForMapTwo.push_back(Teuchos::as<GO>((comm->getSize() + comm->getRank() * numLocalElementsPerMap) + i));
    }

    RCP<const Map> mapOne = MapFactory::Build(lib, numGlobalElementsPerMap, gidsForMapOne(), GO_zero, comm);
    RCP<const Map> mapTwo = MapFactory::Build(lib, numGlobalElementsPerMap, gidsForMapTwo(), GO_zero, comm);

    TEST_ASSERT(Teuchos::nonnull(mapOne));
    TEST_ASSERT(Teuchos::nonnull(mapTwo));

    std::vector<RCP<const Map>> maps;
    maps.push_back(mapOne);
    maps.push_back(mapTwo);

    RCP<const Map> fullMap = MapUtils::concatenateMaps(maps);

    TEST_ASSERT(Teuchos::nonnull(fullMap));
    TEST_EQUALITY_CONST(fullMap->getNodeNumElements(), Teuchos::as<LO>(2 * numLocalElementsPerMap));
    TEST_EQUALITY_CONST(fullMap->getGlobalNumElements(), Teuchos::as<LO>(2 * numGlobalElementsPerMap));

    // Manually merge both lists of GIDs and check with GIDs in full map
    {
      Array<GO> expectedListOfGids;
      for (const auto& gid : gidsForMapOne) expectedListOfGids.push_back(gid);
      for (const auto& gid : gidsForMapTwo) expectedListOfGids.push_back(gid);

      ArrayView<const GO> gidsInFullMap = fullMap->getNodeElementList();

      TEST_COMPARE_ARRAYS(gidsInFullMap, expectedListOfGids);
    }
  }

  // TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MapUtils, AttemptToConcatenateTwoStridedMaps_ThrowException, M, LO, GO, N )
  // {
  //   using Teuchos::Array;
  //   using Teuchos::ArrayView;
  //   using Teuchos::RCP;

  //   using Map = Xpetra::Map<LO,GO,N>;
  //   using MapFactory = Xpetra::MapFactory<LO,GO,N>;
  //   using MapUtils = Xpetra::MapUtils<LO,GO,N>;
  //   using StridedMap = Xpetra::StridedMap<LO,GO,N>;
  //   using StridedMapFactory = Xpetra::StridedMapFactory<LO,GO,N>;

  //   // create a comm
  //   Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
  //   M testMap(1, 0, comm);
  //   Xpetra::UnderlyingLib lib = testMap.lib();

  //   const GO GO_zero = Teuchos::OrdinalTraits<GO>::zero();

  //   const LO numLocalElementsPerMap = 3;
  //   const Xpetra::global_size_t numGlobalElementsPerMap =
  //       Teuchos::as<Xpetra::global_size_t>(comm->getSize() * numLocalElementsPerMap);

  //   Array<GO> gidsForMapOne;
  //   Array<GO> gidsForMapTwo;
  //   for (LO i = 0; i < numLocalElementsPerMap; ++i) {
  //     gidsForMapOne.push_back(Teuchos::as<GO>(comm->getRank() * numLocalElementsPerMap + i));
  //     gidsForMapTwo.push_back(Teuchos::as<GO>((comm->getSize() + comm->getRank() * numLocalElementsPerMap) + i));
  //   }

  //   RCP<const Map> mapOne = MapFactory::Build(lib, numGlobalElementsPerMap, gidsForMapOne(), GO_zero, comm);
  //   RCP<const Map> mapTwo = MapFactory::Build(lib, numGlobalElementsPerMap, gidsForMapTwo(), GO_zero, comm);

  //   TEST_ASSERT(Teuchos::nonnull(mapOne));
  //   TEST_ASSERT(Teuchos::nonnull(mapTwo));

  //   std::vector<size_t> stridingInfo;
  //   stridingInfo.push_back(1);

  //   RCP<const StridedMap> stridedMapOne = StridedMapFactory::Build(mapOne, stridingInfo);
  //   RCP<const StridedMap> stridedMapTwo = StridedMapFactory::Build(mapTwo, stridingInfo);

  //   std::vector<RCP<const Map>> maps;
  //   maps.push_back(stridedMapOne);
  //   maps.push_back(stridedMapTwo);

  //   TEST_THROW(RCP<const Map> fullMap = MapUtils::concatenateMaps(maps), Xpetra::Exceptions::RuntimeError);
  // }

    //
  // INSTANTIATIONS
  //
#ifdef HAVE_XPETRA_TPETRA

  #define XPETRA_TPETRA_TYPES( LO, GO, N) \
    typedef typename Xpetra::TpetraMap<LO,GO,N> M##LO##GO##N;

#endif

#ifdef HAVE_XPETRA_EPETRA

  #define XPETRA_EPETRA_TYPES( LO, GO, N) \
    typedef typename Xpetra::EpetraMapT<GO,N> M##LO##GO##N;

#endif

// List of tests (which run both on Epetra and Tpetra)
#define XP_MAP_INSTANT(LO,GO,N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MapUtils, ConcatenateTwoMaps_ObtainMergedMap, M##LO##GO##N, LO, GO, N) //\
    // TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MapUtils, AttemptToConcatenateTwoStridedMaps_ThrowException, M##LO##GO##N, LO, GO, N)

#if defined(HAVE_XPETRA_TPETRA)

#include <TpetraCore_config.h>
#include <TpetraCore_ETIHelperMacros.h>

TPETRA_ETI_MANGLING_TYPEDEFS()
// no ordinal types as scalar for testing as some tests use ScalarTraits::eps...
TPETRA_INSTANTIATE_LGN ( XPETRA_TPETRA_TYPES )
TPETRA_INSTANTIATE_LGN ( XP_MAP_INSTANT )

#endif

#if defined(HAVE_XPETRA_EPETRA)

#include "Xpetra_Map.hpp" // defines EpetraNode
typedef Xpetra::EpetraNode EpetraNode;
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
XPETRA_EPETRA_TYPES(int,int,EpetraNode)
XP_MAP_INSTANT(int,int,EpetraNode)
#endif
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
typedef long long LongLong;
XPETRA_EPETRA_TYPES(int,LongLong,EpetraNode)
XP_MAP_INSTANT(int,LongLong,EpetraNode)
#endif
#endif

}
