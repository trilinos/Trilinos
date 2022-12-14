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

#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include <Xpetra_UnitTestHelpers.hpp>

#include <Xpetra_BlockedMap.hpp>
// #include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MapExtractor.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_MapFactory.hpp>
// #include <Xpetra_Matrix.hpp>

#include <Xpetra_TpetraCrsMatrix.hpp>
#include <Xpetra_TpetraMap.hpp>

#ifdef HAVE_XPETRA_EPETRA
#include <Xpetra_EpetraCrsMatrix.hpp>
#endif

namespace XpetraMapExtractorFactoryTests {

bool testMpi = true;
double errorTolSlack = 1e+1;

Teuchos::RCP<const Teuchos::Comm<int> > getDefaultComm()
{
  if (testMpi) {
    return Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
  }
  return rcp(new Teuchos::SerialComm<int>());
}

/////////////////////////////////////////////////////

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

//
// UNIT TESTS
//

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( MapExtractorFactory, ConstructFromFullAndPartialMaps, M, MA, Scalar, LO, GO, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  using Map = Xpetra::Map<LO,GO,Node>;
  using MapFactory = Xpetra::MapFactory<LO,GO,Node>;
  using MapExtractor = Xpetra::MapExtractor<Scalar,LO,GO,Node>;
  using MapExtractorFactory = Xpetra::MapExtractorFactory<Scalar,LO,GO,Node>;
  using MapUtils = Xpetra::MapUtils<LO,GO,Node>;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  M testMap(1,0,comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  // Generate partial maps
  const GO nEle1 = 63;
  const GO nEle2 = 29;
  const GO GO_zero = Teuchos::ScalarTraits<GO>::zero();

  const RCP<const Map> map1 = MapFactory::Build(lib, nEle1, GO_zero, comm);
  const RCP<const Map> map2 = MapFactory::Build(lib, nEle2, GO_zero, comm);

  TEST_ASSERT(!map1.is_null());
  TEST_ASSERT(!map2.is_null());

  TEST_EQUALITY_CONST(map1->getGlobalNumElements(), nEle1);
  TEST_EQUALITY_CONST(map2->getGlobalNumElements(), nEle2);

  // Create full map
  std::vector<RCP<const Map>> partialMaps;
  partialMaps.push_back(map1);
  partialMaps.push_back(map2);
  RCP<const Map> fullMap = MapUtils::concatenateMaps(partialMaps);

  TEST_ASSERT(!fullMap.is_null());
  TEST_EQUALITY_CONST(fullMap->getGlobalNumElements(), nEle1 + nEle2);

  RCP<const MapExtractor> mapExtractor = MapExtractorFactory::Build(fullMap, partialMaps);

  TEST_ASSERT(!mapExtractor.is_null());

  // Finally, test the MapExtractor in action
  TEST_ASSERT(mapExtractor->getMap(0)->isSameAs(*map1));
  TEST_ASSERT(mapExtractor->getMap(1)->isSameAs(*map2));
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( MapExtractorFactory, ConstructFromBlockedMap, M, MA, Scalar, LO, GO, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  using BlockedMap = Xpetra::BlockedMap<LO,GO,Node>;
  using Map = Xpetra::Map<LO,GO,Node>;
  using MapFactory = Xpetra::MapFactory<LO,GO,Node>;
  using MapExtractor = Xpetra::MapExtractor<Scalar,LO,GO,Node>;
  using MapExtractorFactory = Xpetra::MapExtractorFactory<Scalar,LO,GO,Node>;
  using MapUtils = Xpetra::MapUtils<LO,GO,Node>;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  M testMap(1,0,comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  // Generate partial maps
  const GO nEle1 = 63;
  const GO nEle2 = 29;
  const GO GO_zero = Teuchos::ScalarTraits<GO>::zero();

  const RCP<const Map> map1 = MapFactory::Build(lib, nEle1, GO_zero, comm);
  const RCP<const Map> map2 = MapFactory::Build(lib, nEle2, GO_zero, comm);

  TEST_ASSERT(!map1.is_null());
  TEST_ASSERT(!map2.is_null());

  TEST_EQUALITY_CONST(map1->getGlobalNumElements(), nEle1);
  TEST_EQUALITY_CONST(map2->getGlobalNumElements(), nEle2);

  // Create full map
  std::vector<RCP<const Map>> partialMaps;
  partialMaps.push_back(map1);
  partialMaps.push_back(map2);
  RCP<const Map> fullMap = MapUtils::concatenateMaps(partialMaps);

  TEST_ASSERT(!fullMap.is_null());
  TEST_EQUALITY_CONST(fullMap->getGlobalNumElements(), nEle1 + nEle2);

  RCP<const BlockedMap> blockedMap = rcp(new BlockedMap(fullMap, partialMaps));

  TEST_ASSERT(!blockedMap.is_null());
  TEST_EQUALITY_CONST(blockedMap->getGlobalNumElements(), nEle1 + nEle2);
  TEST_EQUALITY_CONST(blockedMap->getNumMaps(), 2);

  RCP<const MapExtractor> mapExtractor = MapExtractorFactory::Build(blockedMap);

  TEST_ASSERT(!mapExtractor.is_null());

  // Finally, test the MapExtractor in action
  TEST_ASSERT(mapExtractor->getMap(0)->isSameAs(*map1));
  TEST_ASSERT(mapExtractor->getMap(1)->isSameAs(*map2));
  TEST_ASSERT(mapExtractor->getMap()->isSameAs(*blockedMap));
}

//
// INSTANTIATIONS
//

  #define XPETRA_TPETRA_TYPES( S, LO, GO, N) \
    typedef typename Xpetra::TpetraMap<LO,GO,N> M##LO##GO##N; \
    typedef typename Xpetra::TpetraCrsMatrix<S,LO,GO,N> MA##S##LO##GO##N;

#ifdef HAVE_XPETRA_EPETRA

  #define XPETRA_EPETRA_TYPES( S, LO, GO, N) \
    typedef typename Xpetra::EpetraMapT<GO,N> M##LO##GO##N; \
    typedef typename Xpetra::EpetraCrsMatrixT<GO,N> MA##S##LO##GO##N;

#endif


#define XP_MAPEXTRACTORACTORY_INSTANT(S,LO,GO,N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( MapExtractorFactory, ConstructFromFullAndPartialMaps, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( MapExtractorFactory, ConstructFromBlockedMap, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N )

#include <TpetraCore_config.h>
#include <TpetraCore_ETIHelperMacros.h>

TPETRA_ETI_MANGLING_TYPEDEFS()
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( XPETRA_TPETRA_TYPES )
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( XP_MAPEXTRACTORACTORY_INSTANT )

#if defined(HAVE_XPETRA_EPETRA)

#include "Xpetra_Map.hpp" // defines EpetraNode
typedef Xpetra::EpetraNode EpetraNode;
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
XPETRA_EPETRA_TYPES(double,int,int,EpetraNode)
XP_MAPEXTRACTORACTORY_INSTANT(double,int,int,EpetraNode)
#endif
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
typedef long long LongLong;
XPETRA_EPETRA_TYPES(double,int,LongLong,EpetraNode)
XP_MAPEXTRACTORACTORY_INSTANT(double,int,LongLong,EpetraNode)
#endif

#endif

}
