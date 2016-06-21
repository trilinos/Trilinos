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
#include <Xpetra_UnitTestHelpers.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include <Teuchos_as.hpp>

#include <Xpetra_Map.hpp>

#include <Xpetra_MapUtils.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_BlockedMultiVector.hpp>
#include <Xpetra_Exceptions.hpp>


namespace XpetraBlockMatrixTests {

double errorTolSlack = 1e+1;

Teuchos::RCP<const Teuchos::Comm<int> > getDefaultComm()
{
  return Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
}

/////////////////////////////////////////////////////

TEUCHOS_STATIC_SETUP()
{
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
  clp.setOption(
      "error-tol-slack", &errorTolSlack,
      "Slack off of machine epsilon used to check test results" );
}

//
// Helper routines
//
template<class LocalOrdinal, class GlobalOrdinal, class Node, class MapType>
Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > CreateMap(const std::set<GlobalOrdinal>& gids, const Teuchos::Comm<int>& comm) {
  Teuchos::Array<GlobalOrdinal> mapvec;
  mapvec.reserve(gids.size());
  mapvec.assign(gids.begin(), gids.end());
  GlobalOrdinal count = Teuchos::as<GlobalOrdinal>(mapvec.size());
  GlobalOrdinal gcount;
  Teuchos::reduceAll(comm, Teuchos::REDUCE_SUM, count, Teuchos::outArg(gcount));

  Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map =
      Teuchos::rcp(new MapType(gcount,
          mapvec(),
          0,
          Teuchos::rcpFromRef(comm)));
  mapvec.clear();
  return map;
}

template<class LocalOrdinal, class GlobalOrdinal, class Node, class MapType>
Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > SplitMap(const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> & Amap, const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> & Agiven) {
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Amap.getComm();

  GlobalOrdinal count=0;
  Teuchos::Array<GlobalOrdinal> myaugids(Amap.getNodeNumElements());
  for (size_t i=0; i<Amap.getNodeNumElements(); ++i) {
    const GlobalOrdinal gid = Amap.getGlobalElement(i);
    if (Agiven.isNodeGlobalElement(gid)) continue;
    myaugids[Teuchos::as<GlobalOrdinal>(count)] = gid;
    ++count;
  }
  myaugids.resize(count);
  GlobalOrdinal gcount;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &count, &gcount);
  return Teuchos::rcp(new MapType(gcount,myaugids(),0,comm));
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class MapType>
Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > CreateMultiVector(int noBlocks, Teuchos::RCP<const Teuchos::Comm<int> > comm) {

  typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
  typedef Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVector;
  typedef Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVectorFactory;

  GlobalOrdinal nOverallDOFGidsPerProc = Teuchos::as<GlobalOrdinal>(Teuchos::ScalarTraits<GlobalOrdinal>::pow(2,noBlocks-2)) * 10;

  GlobalOrdinal procOffset = comm->getRank() * nOverallDOFGidsPerProc;

  std::set<GlobalOrdinal> myDOFGids;
  for(GlobalOrdinal i = 0; i < nOverallDOFGidsPerProc; i++)
    myDOFGids.insert(i + procOffset);

  Teuchos::RCP<Map> fullmap = CreateMap<LocalOrdinal,GlobalOrdinal,Node,MapType>(myDOFGids, *comm);

  // create Multivector
  Teuchos::RCP<MultiVector> vv = MultiVectorFactory::Build(fullmap,2,true);

  // fill multivector data (first multivector contains the GID, the second the LID as scalar)
  Teuchos::ArrayRCP< Scalar > vv1 = vv->getDataNonConst(0);
  Teuchos::ArrayRCP< Scalar > vv2 = vv->getDataNonConst(1);
  for(LocalOrdinal i = 0; i < vv->getLocalLength(); ++i) {
    vv1[i] = Teuchos::as<Scalar>(fullmap->getGlobalElement(i));
    vv2[i] = Teuchos::as<Scalar>(i);
  }

  return vv;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class MapType>
Teuchos::RCP<Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > CreateBlockedMultiVector(int noBlocks, Teuchos::RCP<const Teuchos::Comm<int> > comm) {

  typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
  typedef Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> BlockedMultiVector;
  typedef Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVectorFactory;
  typedef Xpetra::MapExtractor<Scalar,LocalOrdinal, GlobalOrdinal, Node> MapExtractor;

  GlobalOrdinal nOverallDOFGidsPerProc = Teuchos::as<GlobalOrdinal>(Teuchos::ScalarTraits<GlobalOrdinal>::pow(2,noBlocks-2)) * 10;

  GlobalOrdinal procOffset = comm->getRank() * nOverallDOFGidsPerProc;

  std::set<GlobalOrdinal> myDOFGids;
  for(GlobalOrdinal i = 0; i < nOverallDOFGidsPerProc; i++)
    myDOFGids.insert(i + procOffset);

  Teuchos::RCP<Map> fullmap = CreateMap<LocalOrdinal,GlobalOrdinal,Node,MapType>(myDOFGids, *comm);

  std::vector<Teuchos::RCP<const Map> > maps(noBlocks, Teuchos::null);
  GlobalOrdinal nPartGIDs = nOverallDOFGidsPerProc;
  Teuchos::RCP<Map> remainingpartmap = fullmap;
  for (int it=0; it<noBlocks; it++) {
    if(it == noBlocks - 1) {
      maps[0] = remainingpartmap;
      break;
    }
    // collect first half of GIDs
    nPartGIDs = nPartGIDs / 2;
    std::set<GlobalOrdinal> myHalfGIDs;
    for(GlobalOrdinal j = 0; j < nPartGIDs; j++)
      myHalfGIDs.insert(j + procOffset);

    Teuchos::RCP<Map> halfmap = CreateMap<LocalOrdinal,GlobalOrdinal,Node,MapType> (myHalfGIDs, *comm);

    Teuchos::RCP<Map> secondmap = SplitMap<LocalOrdinal,GlobalOrdinal,Node,MapType>(*remainingpartmap, *halfmap);
    remainingpartmap = halfmap;

    maps[noBlocks - 1 - it]  = secondmap;
  }

  // create map extractor (Xpetra mode)
  Teuchos::RCP<const MapExtractor> xpMapExtractor = Teuchos::rcp(new MapExtractor(fullmap, maps, false));

  // create Multivector
  Teuchos::RCP<MultiVector> vv = MultiVectorFactory::Build(fullmap,2,true);

  // fill multivector data (first multivector contains the GID, the second the LID as scalar)
  Teuchos::ArrayRCP< Scalar > vv1 = vv->getDataNonConst(0);
  Teuchos::ArrayRCP< Scalar > vv2 = vv->getDataNonConst(1);
  for(LocalOrdinal i = 0; i < vv->getLocalLength(); ++i) {
    vv1[i] = Teuchos::as<Scalar>(fullmap->getGlobalElement(i));
    vv2[i] = Teuchos::as<Scalar>(i);
  }

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv = Teuchos::rcp(new BlockedMultiVector(xpMapExtractor,vv));

  return bvv;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class MapType>
Teuchos::RCP<Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > CreateBlockedMultiVectorThyra(int noBlocks, Teuchos::RCP<const Teuchos::Comm<int> > comm) {

  typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
  typedef Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node> MapFactory;
  typedef Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> BlockedMultiVector;
  typedef Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVectorFactory;
  typedef Xpetra::MapExtractor<Scalar,LocalOrdinal, GlobalOrdinal, Node> MapExtractor;

  MapType testMap(1,0,comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  std::vector<Teuchos::RCP<const Map> > maps(noBlocks, Teuchos::null);
  maps[0] = MapFactory::Build (lib, comm->getSize() * 5, 5, 0, comm);
  for (int it=1; it<noBlocks; it++) {
    GlobalOrdinal localDofs = Teuchos::as<GlobalOrdinal>(Teuchos::ScalarTraits<GlobalOrdinal>::pow(2,it-1)*5);
    maps[it]  = MapFactory::Build (lib, comm->getSize() * localDofs, localDofs, 0, comm);
  }

  // create map extractor
  // To generate the Thyra style map extractor we do not need a full map but only the
  // information about the Map details (i.e. lib and indexBase). We can extract this
  // information from maps[0]
  Teuchos::RCP<const MapExtractor > me =
      Teuchos::rcp(new MapExtractor(maps[0], maps, true));

  // create Multivector
  Teuchos::RCP<MultiVector> vv = MultiVectorFactory::Build(me->getFullMap(),2,true);

  // fill multivector data (first multivector contains the GID, the second the LID as scalar)
  Teuchos::ArrayRCP< Scalar > vv1 = vv->getDataNonConst(0);
  Teuchos::ArrayRCP< Scalar > vv2 = vv->getDataNonConst(1);
  for(LocalOrdinal i = 0; i < vv->getLocalLength(); ++i) {
    vv1[i] = Teuchos::as<Scalar>(vv->getMap()->getGlobalElement(i));
    vv2[i] = Teuchos::as<Scalar>(i);
  }

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv = Teuchos::rcp(new BlockedMultiVector(me,vv));

  return bvv;
}

//
// UNIT TESTS
//

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedMultiVector, Constructor, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 5;

  Teuchos::RCP<MultiVector>         vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  TEST_EQUALITY(bvv->getMapExtractor()->NumMaps(),Teuchos::as<size_t>(noBlocks));
  for (size_t r = 0; r<bvv->getMapExtractor()->NumMaps(); ++r) {
    Teuchos::RCP<MultiVector> bvvi = bvv->getMultiVector(r);
    TEST_EQUALITY(bvvi->getMap()->isSameAs(*(bvv->getMapExtractor()->getMap(r))),true);
    Teuchos::ArrayRCP<const Scalar > bvvi1 = bvvi->getData(0);
    Teuchos::ArrayRCP<const Scalar > bvvi2 = bvvi->getData(1);
    for(LO l = 0; l < bvvi->getLocalLength(); ++l) {
      TEST_EQUALITY(bvvi1[l],Teuchos::as<Scalar>(bvvi->getMap()->getGlobalElement(l)));
    }
  }

  typedef typename STS::magnitudeType Magnitude;
  Teuchos::Array<Magnitude> bnorms(bvv->getNumVectors());
  Teuchos::Array<Magnitude> fnorms(vv->getNumVectors());

  TEST_NOTHROW( vv->norm1(fnorms) );
  TEST_NOTHROW( bvv->norm1(bnorms) );
  TEST_COMPARE_FLOATING_ARRAYS(fnorms,bnorms,Teuchos::ScalarTraits<Magnitude>::zero());
  TEST_NOTHROW( vv->norm2(fnorms) );
  TEST_NOTHROW( bvv->norm2(bnorms) );
  TEST_COMPARE_FLOATING_ARRAYS(fnorms,bnorms,Teuchos::ScalarTraits<Magnitude>::zero());
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedMultiVector, Norm1, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 2;

  Teuchos::RCP<MultiVector>         vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  typedef typename STS::magnitudeType Magnitude;
  Teuchos::Array<Magnitude> bnorms(bvv->getNumVectors());
  Teuchos::Array<Magnitude> fnorms(vv->getNumVectors());

  TEST_NOTHROW( vv->norm1(fnorms) );
  TEST_NOTHROW( bvv->norm1(bnorms) );
  TEST_COMPARE_FLOATING_ARRAYS(fnorms,bnorms,Teuchos::ScalarTraits<Magnitude>::zero());
  Magnitude result = Teuchos::ScalarTraits<Magnitude>::zero();
  for(GO gg = 0; gg < vv->getMap()->getGlobalNumElements(); gg++)
    result += Teuchos::as<Magnitude>(gg);
  TEST_EQUALITY( bnorms[0], result);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv2 = CreateBlockedMultiVectorThyra<Scalar, LO, GO, Node, M>(noBlocks, comm);
  Teuchos::Array<Magnitude> bnorms2(bvv2->getNumVectors());
  TEST_NOTHROW( bvv2->norm1(bnorms2) );
  TEST_EQUALITY( bnorms2[0], result);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedMultiVector, Norm2, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 2;

  Teuchos::RCP<MultiVector>         vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  typedef typename STS::magnitudeType Magnitude;
  Teuchos::Array<Magnitude> bnorms(bvv->getNumVectors());
  Teuchos::Array<Magnitude> fnorms(vv->getNumVectors());

  TEST_NOTHROW( vv->norm2(fnorms) );
  TEST_NOTHROW( bvv->norm2(bnorms) );
  TEST_COMPARE_FLOATING_ARRAYS(fnorms,bnorms,Teuchos::ScalarTraits<Magnitude>::zero());
  Magnitude result = Teuchos::ScalarTraits<Magnitude>::zero();
  for(GO gg = 0; gg < vv->getMap()->getGlobalNumElements(); gg++)
    result += Teuchos::as<Magnitude>(gg) * Teuchos::as<Magnitude>(gg);
  result = Teuchos::ScalarTraits<Magnitude>::squareroot(result);
  TEST_EQUALITY( bnorms[0], result);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv2 = CreateBlockedMultiVectorThyra<Scalar, LO, GO, Node, M>(noBlocks, comm);
  Teuchos::Array<Magnitude> bnorms2(bvv2->getNumVectors());
  TEST_NOTHROW( bvv2->norm2(bnorms2) );
  TEST_COMPARE( bnorms2[0] - result, <, 1e-10);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedMultiVector, NormInf, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 2;

  Teuchos::RCP<MultiVector>         vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  typedef typename STS::magnitudeType Magnitude;
  Teuchos::Array<Magnitude> bnorms(bvv->getNumVectors());
  Teuchos::Array<Magnitude> fnorms(vv->getNumVectors());

  TEST_NOTHROW( vv->normInf(fnorms) );
  TEST_NOTHROW( bvv->normInf(bnorms) );
  TEST_COMPARE_FLOATING_ARRAYS(fnorms,bnorms,Teuchos::ScalarTraits<Magnitude>::zero());
  Magnitude result = Teuchos::ScalarTraits<Magnitude>::zero();
  for(GO gg = 0; gg < vv->getMap()->getGlobalNumElements(); gg++)
    result = std::max(result, Teuchos::as<Magnitude>(gg));
  TEST_EQUALITY( bnorms[0], result);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv2 = CreateBlockedMultiVectorThyra<Scalar, LO, GO, Node, M>(noBlocks, comm);
  Teuchos::Array<Magnitude> bnorms2(bvv2->getNumVectors());
  TEST_NOTHROW( bvv2->normInf(bnorms2) );
  TEST_COMPARE( bnorms2[0] - result, <, 1e-10);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedMultiVector, Scale, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 2;

  Teuchos::RCP<MultiVector>         vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  typedef typename STS::magnitudeType Magnitude;
  Teuchos::Array<Magnitude> bnorms(bvv->getNumVectors());
  Teuchos::Array<Magnitude> fnorms(vv->getNumVectors());

  TEST_NOTHROW( vv->normInf(fnorms) );
  TEST_NOTHROW( bvv->normInf(bnorms) );
  TEST_COMPARE_FLOATING_ARRAYS(fnorms,bnorms,Teuchos::ScalarTraits<Magnitude>::zero());
  Magnitude myresult = Teuchos::ScalarTraits<Magnitude>::zero();
  for(GO gg = 0; gg < vv->getMap()->getGlobalNumElements(); gg++)
    myresult = std::max(myresult, Teuchos::as<Magnitude>(gg));
  TEST_EQUALITY( bnorms[0], myresult);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv2 = CreateBlockedMultiVectorThyra<Scalar, LO, GO, Node, M>(noBlocks, comm);
  Teuchos::Array<Magnitude> bnorms2(bvv2->getNumVectors());
  TEST_NOTHROW( bvv2->normInf(bnorms2) );
  TEST_COMPARE( bnorms2[0] - myresult, <, 1e-10);

  bvv->scale(Teuchos::as<Scalar>(2.0));
  vv->scale(Teuchos::as<Scalar>(2.0));
  Teuchos::Array<Magnitude> scaled_bnorms(bvv->getNumVectors());
  Teuchos::Array<Magnitude> scaled_fnorms(vv->getNumVectors());
  TEST_NOTHROW( vv->normInf(scaled_fnorms) );
  TEST_NOTHROW( bvv->normInf(scaled_bnorms) );
  TEST_COMPARE_FLOATING_ARRAYS(scaled_fnorms,scaled_bnorms,Teuchos::ScalarTraits<Magnitude>::zero());
  myresult = Teuchos::ScalarTraits<Magnitude>::zero();
  for(GO gg = 0; gg < vv->getMap()->getGlobalNumElements(); gg++)
    myresult = std::max(myresult, Teuchos::as<Magnitude>(gg));
  TEST_EQUALITY( scaled_bnorms[0], Teuchos::as<Magnitude>(2.0) * myresult);

  // create BlockedMultiVector
  bvv2 = CreateBlockedMultiVectorThyra<Scalar, LO, GO, Node, M>(noBlocks, comm);
  bvv2->scale(Teuchos::as<Scalar>(2.0));
  Teuchos::Array<Magnitude> scaled_bnorms2(bvv2->getNumVectors());
  TEST_NOTHROW( bvv2->normInf(scaled_bnorms2) );
  TEST_COMPARE( scaled_bnorms2[0] - Teuchos::as<Magnitude>(2.0) * myresult, <, 1e-10);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedMultiVector, ExtractVector, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::Map<LO, GO, Node> Map;
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  typedef Xpetra::MapExtractor<Scalar,LO,GO,Node> MapExtractor;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename STS::magnitudeType Magnitude;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 5;

  // create full vector
  Teuchos::RCP<MultiVector>         vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // extract map extractor
  Teuchos::RCP<const MapExtractor> me  = bvv->getMapExtractor();

  for(size_t r = 0; r < me->NumMaps(); ++r) {
    Teuchos::RCP<MultiVector> partB = me->ExtractVector(bvv,r);
    Teuchos::RCP<MultiVector> partV = me->ExtractVector(vv,r);
    TEST_EQUALITY(partB->getMap()->isSameAs(*(partV->getMap())),true);
    Teuchos::ArrayRCP<const Scalar > partBd = partB->getData(0);
    Teuchos::ArrayRCP<const Scalar > partVd = partV->getData(0);
    TEST_COMPARE_FLOATING_ARRAYS(partBd,partVd,Teuchos::ScalarTraits<Magnitude>::zero());
    TEST_THROW(partB = me->ExtractVector(bvv,r,true),Xpetra::Exceptions::RuntimeError);
    TEST_THROW(partV = me->ExtractVector(vv,r,true),Xpetra::Exceptions::RuntimeError);
  }

  // create a new faulty MapExtractor
  std::vector<Teuchos::RCP<const Map> > maps(noBlocks, Teuchos::null);
  for(size_t r = 0; r < me->NumMaps(); ++r) {
    maps[me->NumMaps() - 1 - r] = me->getMap(r);
  }
  Teuchos::RCP<const MapExtractor> meFaulty = Teuchos::rcp(new MapExtractor(me->getFullMap(), maps, false));

  // the faulty map extractor reverses the ordering of maps. We have five partial maps (0-4)
  // Therefore, partial map with index 2 are the same in the faulty and the original map extractor and
  // the ExtractVector call then succeeds
  for(size_t r = 0; r < me->NumMaps(); ++r) {
    TEST_NOTHROW(Teuchos::RCP<MultiVector> partV = meFaulty->ExtractVector(vv,r));
    if(r!=2)      {
      TEST_NOTHROW(Teuchos::RCP<MultiVector> partB = meFaulty->ExtractVector(bvv,r));
      Teuchos::RCP<MultiVector> partV = meFaulty->ExtractVector(vv,r);
      Teuchos::RCP<MultiVector> partB = meFaulty->ExtractVector(bvv,r);
      TEST_EQUALITY(partB->getMap()->isSameAs(*(partV->getMap())),false);
    }
    else if(r==2) {
      TEST_NOTHROW(Teuchos::RCP<MultiVector> partB = meFaulty->ExtractVector(bvv,r));
      Teuchos::RCP<MultiVector> partV = meFaulty->ExtractVector(vv,r);
      Teuchos::RCP<MultiVector> partB = meFaulty->ExtractVector(bvv,r);
      TEST_EQUALITY(partB->getMap()->isSameAs(*(partV->getMap())),true);
      Teuchos::ArrayRCP<const Scalar > partBd = partB->getData(0);
      Teuchos::ArrayRCP<const Scalar > partVd = partV->getData(0);
      TEST_COMPARE_FLOATING_ARRAYS(partBd,partVd,Teuchos::ScalarTraits<Magnitude>::zero());
    }
  }

}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedMultiVector, ExtractVectorThyra, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  typedef Xpetra::MapExtractor<Scalar,LO,GO,Node> MapExtractor;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename STS::magnitudeType Magnitude;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 5;

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv = CreateBlockedMultiVectorThyra<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create full vector
  Teuchos::RCP<MultiVector>         vv = bvv->Merge();

  // extract map extractor
  Teuchos::RCP<const MapExtractor> me  = bvv->getMapExtractor();

  for(size_t r = 0; r < me->NumMaps(); ++r) {
    Teuchos::RCP<MultiVector> partB = me->ExtractVector(bvv,r,true);
    Teuchos::RCP<MultiVector> partV = me->ExtractVector(vv,r,false);
    TEST_EQUALITY(partB->getMap()->isSameAs(*(partV->getMap())) == false || r==0,true);
    TEST_EQUALITY(partB->getMap()->getMinAllGlobalIndex(),0);
    partV = me->ExtractVector(vv,r,true);
    TEST_EQUALITY(partB->getMap()->isSameAs(*(partV->getMap())),true);
    Teuchos::ArrayRCP<const Scalar > partBd = partB->getData(0);
    Teuchos::ArrayRCP<const Scalar > partVd = partV->getData(0);
    TEST_COMPARE_FLOATING_ARRAYS(partBd,partVd,Teuchos::ScalarTraits<Magnitude>::zero());
  }

}


TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedMultiVector, InsertVector, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  typedef Xpetra::MapExtractor<Scalar,LO,GO,Node> MapExtractor;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename STS::magnitudeType Magnitude;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 5;

  // create full vector
  Teuchos::RCP<MultiVector>         vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  TEST_EQUALITY(bvv->getNumVectors(), 2);

  // extract map extractor
  Teuchos::RCP<const MapExtractor> me  = bvv->getMapExtractor();


  for(size_t r = 0; r < me->NumMaps(); ++r) {
    Teuchos::RCP<MultiVector> part = me->getVector(r,bvv->getNumVectors(),false);
    part->putScalar(STS::one());
    me->InsertVector(part,r,bvv);
  }

  for(size_t r = 0; r < me->NumMaps(); ++r) {
    Teuchos::RCP<const MultiVector> part = me->ExtractVector(bvv,r);
    Teuchos::ArrayRCP<const Scalar > partd1 = part->getData(0);
    Teuchos::ArrayRCP<const Scalar > partd2 = part->getData(1);
    for(LO l = 0; l < part->getLocalLength(); l++)
      TEST_EQUALITY(partd1[l], STS::one());
    TEST_COMPARE_FLOATING_ARRAYS(partd1,partd2,Teuchos::ScalarTraits<Magnitude>::zero());
  }

  // create malicious multivector
  Teuchos::RCP<MultiVector> part1 = me->getVector(0,23,false);
  TEST_THROW(me->InsertVector(part1,0,bvv),Xpetra::Exceptions::RuntimeError);
  Teuchos::RCP<MultiVector> part2 = me->getVector(0,2,false);
  TEST_THROW(me->InsertVector(part2,1,bvv),Xpetra::Exceptions::RuntimeError);
  TEST_THROW(Teuchos::RCP<MultiVector> part3 = me->getVector(1,2,true),Xpetra::Exceptions::RuntimeError);

}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedMultiVector, InsertVectorThyra, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  typedef Xpetra::MapExtractor<Scalar,LO,GO,Node> MapExtractor;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename STS::magnitudeType Magnitude;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 5;

  // create full vector
  Teuchos::RCP<MultiVector>         vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv = CreateBlockedMultiVectorThyra<Scalar, LO, GO, Node, M>(noBlocks, comm);

  TEST_EQUALITY(bvv->getNumVectors(), 2);

  // extract map extractor
  Teuchos::RCP<const MapExtractor> me  = bvv->getMapExtractor();
  TEST_EQUALITY(me->getThyraMode(),true);

  for(size_t r = 0; r < me->NumMaps(); ++r) {
    Teuchos::RCP<MultiVector> part = me->getVector(r,bvv->getNumVectors(),true);
    TEST_EQUALITY(part->getMap()->getMinAllGlobalIndex(),0);
    TEST_NOTHROW(part->putScalar(STS::one()));
    TEST_NOTHROW(me->InsertVector(part,r,bvv,me->getThyraMode()));
  }

  for(size_t r = 0; r < me->NumMaps(); ++r) {
    Teuchos::RCP<const MultiVector> part = me->ExtractVector(bvv,r,me->getThyraMode());
    Teuchos::ArrayRCP<const Scalar > partd1 = part->getData(0);
    Teuchos::ArrayRCP<const Scalar > partd2 = part->getData(1);
    for(LO l = 0; l < part->getLocalLength(); l++)
      TEST_EQUALITY(partd1[l], STS::one());
    TEST_COMPARE_FLOATING_ARRAYS(partd1,partd2,Teuchos::ScalarTraits<Magnitude>::zero());
  }

  // create malicious multivector
  Teuchos::RCP<MultiVector> part1 = me->getVector(0,23,true);
  TEST_THROW(me->InsertVector(part1,0,bvv),Xpetra::Exceptions::RuntimeError);
  // unfortunately, in Thyra mode there is no error thrown since the vectors in
  // block 0 and 1 have the same length (and the same GIDs)
  Teuchos::RCP<MultiVector> part2 = me->getVector(0,2,true);
  TEST_NOTHROW(me->InsertVector(part2,1,bvv,me->getThyraMode()));
  // This should throw, thought
  Teuchos::RCP<MultiVector> part3 = me->getVector(0,2,true);
  TEST_THROW(me->InsertVector(part2,2,bvv,me->getThyraMode()),Xpetra::Exceptions::RuntimeError);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedMultiVector, UpdateVector1, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename STS::magnitudeType Magnitude;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 5;

  // create full vector
  Teuchos::RCP<MultiVector>         vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv1 = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);
  Teuchos::RCP<BlockedMultiVector> bvv2 = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  TEST_EQUALITY(bvv1->getNumVectors(), 2);
  TEST_EQUALITY(bvv2->getNumVectors(), 2);

  TEST_NOTHROW(bvv1->update(-0.35*STS::one(), *bvv2, 0.7*STS::one()));
  TEST_NOTHROW(bvv1->update(-0.35*STS::one(), *bvv2, STS::one()));

  typedef typename STS::magnitudeType Magnitude;
  Teuchos::Array<Magnitude> bnorms1(bvv1->getNumVectors());
  Teuchos::Array<Magnitude> bnorms2(vv->getNumVectors());
  TEST_NOTHROW( bvv1->norm1(bnorms1) );
  TEST_EQUALITY( bnorms1[0], Teuchos::ScalarTraits<Magnitude>::zero());
  TEST_EQUALITY( bnorms1[1], Teuchos::ScalarTraits<Magnitude>::zero());
  TEST_NOTHROW( vv->norm1(bnorms1) );
  TEST_NOTHROW( bvv2->norm1(bnorms2) );
  TEST_COMPARE_FLOATING_ARRAYS(bnorms1,bnorms2,Teuchos::ScalarTraits<Magnitude>::zero());
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedMultiVector, UpdateVector1b, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  typedef Xpetra::MultiVectorFactory<Scalar, LO, GO, Node> MultiVectorFactory;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename STS::magnitudeType Magnitude;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 5;

  // create full vector
  Teuchos::RCP<MultiVector>         vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv1 = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);
  TEST_EQUALITY(bvv1->getNumVectors(), 2);

  typedef typename STS::magnitudeType Magnitude;
  Teuchos::Array<Magnitude> bnorms1(bvv1->getNumVectors());
  Teuchos::Array<Magnitude> bnorms2(vv->getNumVectors());

  TEST_NOTHROW(bvv1->update(-0.5* STS::one(), *vv, STS::one()));
  TEST_NOTHROW( bvv1->norm1(bnorms1) );
  TEST_NOTHROW( vv->norm1(bnorms2) );
  TEST_EQUALITY( bnorms1[0], 0.5 * bnorms2[0]);
  TEST_EQUALITY( bnorms1[1], 0.5 * bnorms2[1]);

  // create faulty multivector
  Teuchos::RCP<MultiVector> vvx = MultiVectorFactory::Build(bvv1->getMapExtractor()->getMap(0),2,true);
  TEST_THROW(bvv1->update(STS::one(), *vvx, STS::one()), Xpetra::Exceptions::RuntimeError);
  vvx = MultiVectorFactory::Build(bvv1->getMap(),1,true);
  TEST_THROW(bvv1->update(STS::one(), *vvx, STS::one()), Xpetra::Exceptions::RuntimeError);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedMultiVector, UpdateVector2, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename STS::magnitudeType Magnitude;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 5;

  // create full vector
  Teuchos::RCP<MultiVector>         vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv1 = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);
  Teuchos::RCP<BlockedMultiVector> bvv2 = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);
  Teuchos::RCP<BlockedMultiVector> bvv3 = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  TEST_EQUALITY(bvv1->getNumVectors(), 2);
  TEST_EQUALITY(bvv2->getNumVectors(), 2);
  TEST_EQUALITY(bvv3->getNumVectors(), 2);

  TEST_NOTHROW(bvv1->update(-0.25 * STS::one(), *bvv2, -0.25 * STS::one(), *bvv3, 0.5 * STS::one()));

  typedef typename STS::magnitudeType Magnitude;
  Teuchos::Array<Magnitude> bnorms1(bvv1->getNumVectors());
  Teuchos::Array<Magnitude> bnorms2(vv->getNumVectors());
  TEST_NOTHROW( bvv1->norm1(bnorms1) );
  TEST_EQUALITY( bnorms1[0], Teuchos::ScalarTraits<Magnitude>::zero());
  TEST_EQUALITY( bnorms1[1], Teuchos::ScalarTraits<Magnitude>::zero());
  TEST_NOTHROW( vv->norm1(bnorms1) );
  TEST_NOTHROW( bvv2->norm1(bnorms2) );
  TEST_COMPARE_FLOATING_ARRAYS(bnorms1,bnorms2,Teuchos::ScalarTraits<Magnitude>::zero());
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedMultiVector, PutScalar, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename STS::magnitudeType Magnitude;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 5;

  // create full vector
  Teuchos::RCP<MultiVector>         vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  TEST_NOTHROW(bvv->putScalar(1.0*STS::one()));

  typedef typename STS::magnitudeType Magnitude;
  Teuchos::Array<Magnitude> bnorms(bvv->getNumVectors());
  TEST_NOTHROW( bvv->norm1(bnorms) );
  TEST_EQUALITY( bnorms[0], Teuchos::as<Magnitude>(bvv->getMapExtractor()->getFullMap()->getGlobalNumElements()));
  TEST_EQUALITY( bnorms[1], Teuchos::as<Magnitude>(bvv->getMapExtractor()->getFullMap()->getGlobalNumElements()));

  TEST_NOTHROW(bvv->putScalar(3.0*STS::one()));

  for(size_t r = 0; r < bvv->getMapExtractor()->NumMaps(); ++r) {
    Teuchos::RCP<const MultiVector> part = bvv->getMapExtractor()->ExtractVector(bvv,r);
    Teuchos::ArrayRCP<const Scalar > partd1 = part->getData(0);
    Teuchos::ArrayRCP<const Scalar > partd2 = part->getData(1);
    for(LO l = 0; l < part->getLocalLength(); l++)
      TEST_EQUALITY(partd1[l], 3.0 * STS::one());
    TEST_COMPARE_FLOATING_ARRAYS(partd1,partd2,Teuchos::ScalarTraits<Magnitude>::zero());
  }
}


//
// INSTANTIATIONS
//
#ifdef HAVE_XPETRA_TPETRA

  #define XPETRA_TPETRA_TYPES( S, LO, GO, N) \
    typedef typename Xpetra::TpetraMap<LO,GO,N> M##LO##GO##N; \
    typedef typename Xpetra::TpetraMultiVector<S,LO,GO,N> MV##S##LO##GO##N;

#endif

#ifdef HAVE_XPETRA_EPETRA

  #define XPETRA_EPETRA_TYPES( S, LO, GO, N) \
    typedef typename Xpetra::EpetraMapT<GO,N> M##LO##GO##N; \
    typedef typename Xpetra::EpetraMultiVectorT<GO,N> MV##S##LO##GO##N;

#endif

#define XP_BLOCKEDMULTIVECTOR_INSTANT(S,LO,GO,N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedMultiVector, Constructor, M##LO##GO##N , MV##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedMultiVector, Norm1, M##LO##GO##N , MV##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedMultiVector, Norm2, M##LO##GO##N , MV##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedMultiVector, NormInf, M##LO##GO##N , MV##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedMultiVector, Scale, M##LO##GO##N , MV##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedMultiVector, ExtractVector, M##LO##GO##N , MV##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedMultiVector, ExtractVectorThyra, M##LO##GO##N , MV##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedMultiVector, InsertVector, M##LO##GO##N , MV##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedMultiVector, InsertVectorThyra, M##LO##GO##N , MV##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedMultiVector, UpdateVector1, M##LO##GO##N , MV##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedMultiVector, UpdateVector1b, M##LO##GO##N , MV##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedMultiVector, UpdateVector2, M##LO##GO##N , MV##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedMultiVector, PutScalar, M##LO##GO##N , MV##S##LO##GO##N, S, LO, GO, N ) \

// List of tests which run only with Tpetra
#define XP_TPETRA_BLOCKEDMULTIVECTOR_INSTANT(S,LO,GO,N)

// List of tests which run only with Epetra
#define XP_EPETRA_BLOCKEDMULTIVECTOR_INSTANT(S,LO,GO,N)

#if defined(HAVE_XPETRA_TPETRA)

#include <TpetraCore_config.h>
#include <TpetraCore_ETIHelperMacros.h>

TPETRA_ETI_MANGLING_TYPEDEFS()
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( XPETRA_TPETRA_TYPES )
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( XP_TPETRA_BLOCKEDMULTIVECTOR_INSTANT )
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( XP_BLOCKEDMULTIVECTOR_INSTANT )

#endif

#if defined(HAVE_XPETRA_EPETRA)

#include "Xpetra_Map.hpp" // defines EpetraNode
typedef Xpetra::EpetraNode EpetraNode;
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
XPETRA_EPETRA_TYPES(double,int,int,EpetraNode)
XP_EPETRA_BLOCKEDMULTIVECTOR_INSTANT(double,int,int,EpetraNode)
XP_BLOCKEDMULTIVECTOR_INSTANT(double,int,int,EpetraNode)
#endif
// EpetraExt routines are not working with 64 bit
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
typedef long long LongLong;
XPETRA_EPETRA_TYPES(double,int,LongLong,EpetraNode)
XP_EPETRA_BLOCKEDMULTIVECTOR_INSTANT(double,int,LongLong,EpetraNode)
XP_EPETRA_BLOCKEDMULTIVECTOR_INSTANT(double,int,LongLong,EpetraNode)
#endif

#endif

}
