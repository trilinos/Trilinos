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
#include <Xpetra_MapExtractor.hpp>
#include <Xpetra_BlockedMultiVector.hpp>
#include <Xpetra_ReorderedBlockedMultiVector.hpp>
#include <Xpetra_Exceptions.hpp>
#include <Xpetra_ThyraUtils.hpp>

#ifdef HAVE_XPETRA_THYRA
#include <Thyra_DefaultProductVectorSpace.hpp>
#endif

#include <set>

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
  Teuchos::Array<GlobalOrdinal> myaugids(Amap.getLocalNumElements());
  for (size_t i=0; i<Amap.getLocalNumElements(); ++i) {
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
  for(LocalOrdinal i = 0; i < Teuchos::as<LocalOrdinal>(vv->getLocalLength()); ++i) {
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
  for(LocalOrdinal i = 0; i < Teuchos::as<LocalOrdinal>(vv->getLocalLength()); ++i) {
    vv1[i] = Teuchos::as<Scalar>(fullmap->getGlobalElement(i));
    vv2[i] = Teuchos::as<Scalar>(i);
  }

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv = Teuchos::rcp(new BlockedMultiVector(xpMapExtractor,vv));

  return bvv;
}

//
// UNIT TESTS
//

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( ThyraBlockedMultiVector, Constructor, M, MA, Scalar, LO, GO, Node )
{
#ifdef HAVE_XPETRA_THYRA
  typedef Xpetra::Map<LO, GO, Node> Map;
  typedef Xpetra::BlockedMap<LO, GO, Node> BlockedMap;
  typedef Xpetra::MapFactory<LO, GO, Node> MapFactory;
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  typedef Xpetra::ThyraUtils<Scalar,LO,GO,Node> ThyraUtils;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 5;

  Teuchos::RCP<MultiVector>         vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  TEST_EQUALITY(bvv->getBlockedMap()->getNumMaps(),Teuchos::as<size_t>(noBlocks));
  for (size_t r = 0; r<bvv->getBlockedMap()->getNumMaps(); ++r) {
    Teuchos::RCP<MultiVector> bvvi = bvv->getMultiVector(r);
    TEST_EQUALITY(bvvi->getMap()->isSameAs(*(bvv->getBlockedMap()->getMap(r))),true);
    Teuchos::ArrayRCP<const Scalar > bvvi1 = bvvi->getData(0);
    Teuchos::ArrayRCP<const Scalar > bvvi2 = bvvi->getData(1);
    for(LO l = 0; l < Teuchos::as<LO>(bvvi->getLocalLength()); ++l) {
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

  M testMap(1,0,comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  Teuchos::RCP<const Map> map1 = MapFactory::Build(lib, 10, 0, comm);
  Teuchos::RCP<const Map> map2 = MapFactory::Build(lib, 15, 0, comm);
  Teuchos::RCP<const Map> map3 = MapFactory::Build(lib, 18, 0, comm);

  // create Thyra product vector space
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > vs1 = ThyraUtils::toThyra(map1);
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > vs2 = ThyraUtils::toThyra(map2);
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > vs3 = ThyraUtils::toThyra(map3);

  Teuchos::Array<Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > > vecSpacesInner(2);
  Teuchos::Array<Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > > vecSpacesOuter(2);

  vecSpacesInner[0] = vs2;
  vecSpacesInner[1] = vs3;

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > psInner =
    Thyra::productVectorSpace<Scalar>(vecSpacesInner());

  vecSpacesOuter[0] = vs1;
  vecSpacesOuter[1] = psInner;

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > ps =
    Thyra::productVectorSpace<Scalar>(vecSpacesOuter());

  Teuchos::RCP<const Thyra::ProductVectorSpaceBase<Scalar> > pps =
    Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorSpaceBase<Scalar> >(ps);
  TEST_EQUALITY(pps.is_null(),false);

  Teuchos::RCP<const Map> ppm = ThyraUtils::toXpetra(ps, comm);
  TEST_EQUALITY(ppm.is_null(),false);

  TEST_EQUALITY(ppm->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(ppm->getMaxAllGlobalIndex(), 42);

  Teuchos::RCP<const BlockedMap> ppbm = Teuchos::rcp_dynamic_cast<const BlockedMap>(ppm);
  TEST_EQUALITY(ppbm.is_null(),false);

  TEST_EQUALITY(ppbm->getMap(0,false)->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(ppbm->getMap(0,false)->getMaxAllGlobalIndex(), 9);

  TEST_EQUALITY(ppbm->getMap(1,false)->getMinAllGlobalIndex(), 10);
  TEST_EQUALITY(ppbm->getMap(1,false)->getMaxAllGlobalIndex(), 42);

  Teuchos::RCP<const Map> pp2m = ppbm->getMap(1,true);
  Teuchos::RCP<const BlockedMap> pp2bm = Teuchos::rcp_dynamic_cast<const BlockedMap>(pp2m);
  TEST_EQUALITY(pp2bm.is_null(),false);

  // there are no nested Xpetra style blocked maps!
  Teuchos::RCP<const Map> pp3m = ppbm->getMap(1,false);
  Teuchos::RCP<const BlockedMap> pp3bm = Teuchos::rcp_dynamic_cast<const BlockedMap>(pp3m);
  TEST_EQUALITY(pp3bm.is_null(),false);
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( ThyraBlockedMultiVector, ConstructorNested, M, MA, Scalar, LO, GO, Node )
{
#ifdef HAVE_XPETRA_THYRA
  typedef Xpetra::Map<LO, GO, Node> Map;
  typedef Xpetra::BlockedMap<LO, GO, Node> BlockedMap;
  typedef Xpetra::MapFactory<LO, GO, Node> MapFactory;
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::MultiVectorFactory<Scalar, LO, GO, Node> MultiVectorFactory;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;

  typedef Xpetra::ThyraUtils<Scalar,LO,GO,Node> ThyraUtils;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  M testMap(1,0,comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  Teuchos::RCP<const Map> map1 = MapFactory::Build(lib, 10, 0, comm);
  Teuchos::RCP<const Map> map2 = MapFactory::Build(lib, 15, 0, comm);
  Teuchos::RCP<const Map> map3 = MapFactory::Build(lib, 18, 0, comm);

  // create Thyra product vector space
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > vs1 = ThyraUtils::toThyra(map1);
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > vs2 = ThyraUtils::toThyra(map2);
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > vs3 = ThyraUtils::toThyra(map3);

  Teuchos::Array<Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > > vecSpacesInner(2);
  Teuchos::Array<Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > > vecSpacesOuter(2);

  vecSpacesInner[0] = vs2;
  vecSpacesInner[1] = vs3;

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > psInner =
    Thyra::productVectorSpace<Scalar>(vecSpacesInner());

  vecSpacesOuter[0] = vs1;
  vecSpacesOuter[1] = psInner;

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > ps =
    Thyra::productVectorSpace<Scalar>(vecSpacesOuter());

  Teuchos::RCP<const Thyra::ProductVectorSpaceBase<Scalar> > pps =
    Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorSpaceBase<Scalar> >(ps);
  TEST_EQUALITY(pps.is_null(),false);

  Teuchos::RCP<const Map> ppm = ThyraUtils::toXpetra(ps, comm);
  TEST_EQUALITY(ppm.is_null(),false);

  TEST_EQUALITY(ppm->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(ppm->getMaxAllGlobalIndex(), 42);

  Teuchos::RCP<const BlockedMap> ppbm = Teuchos::rcp_dynamic_cast<const BlockedMap>(ppm);
  TEST_EQUALITY(ppbm.is_null(),false);
  TEST_EQUALITY(ppbm->getThyraMode(),true);

  TEST_EQUALITY(ppbm->getMap(0,false)->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(ppbm->getMap(0,false)->getMaxAllGlobalIndex(), 9);

  TEST_EQUALITY(ppbm->getMap(1,false)->getMinAllGlobalIndex(), 10);
  TEST_EQUALITY(ppbm->getMap(1,false)->getMaxAllGlobalIndex(), 42);

  Teuchos::RCP<const Map> pp2m = ppbm->getMap(1,true);
  Teuchos::RCP<const BlockedMap> pp2bm = Teuchos::rcp_dynamic_cast<const BlockedMap>(pp2m);
  TEST_EQUALITY(pp2bm.is_null(),false);
  TEST_EQUALITY(pp2bm->getThyraMode(),true);

  // there are no nested Xpetra style blocked maps!
  Teuchos::RCP<const Map> pp3m = ppbm->getMap(1,false);
  Teuchos::RCP<const BlockedMap> pp3bm = Teuchos::rcp_dynamic_cast<const BlockedMap>(pp3m);
  TEST_EQUALITY(pp3bm.is_null(),false);
  TEST_EQUALITY(pp3bm->getThyraMode(),true);

  // create sub multivectors
  Teuchos::RCP<MultiVector> smv1 = MultiVectorFactory::Build(map1,2,true);
  {
    Teuchos::ArrayRCP< Scalar > vv1 = smv1->getDataNonConst(0);
    Teuchos::ArrayRCP< Scalar > vv2 = smv1->getDataNonConst(1);
    for(LO i = 0; i < Teuchos::as<LO>(smv1->getLocalLength()); ++i) {
      vv1[i] = Teuchos::as<Scalar>(comm->getRank() * 100 + i);
      vv2[i] = Teuchos::as<Scalar>(i);
    }
  }

  Teuchos::RCP<MultiVector> smv2 = MultiVectorFactory::Build(map2,2,true);
  {
    Teuchos::ArrayRCP< Scalar > vv1 = smv2->getDataNonConst(0);
    Teuchos::ArrayRCP< Scalar > vv2 = smv2->getDataNonConst(1);
    for(LO i = 0; i < Teuchos::as<LO>(smv2->getLocalLength()); ++i) {
      vv1[i] = Teuchos::as<Scalar>(comm->getRank() * 100 + 10 + i);
      vv2[i] = Teuchos::as<Scalar>(i);
    }
  }

  Teuchos::RCP<MultiVector> smv3 = MultiVectorFactory::Build(map3,2,true);
  {
    Teuchos::ArrayRCP< Scalar > vv1 = smv3->getDataNonConst(0);
    Teuchos::ArrayRCP< Scalar > vv2 = smv3->getDataNonConst(1);
    for(LO i = 0; i < Teuchos::as<LO>(smv3->getLocalLength()); ++i) {
      vv1[i] = Teuchos::as<Scalar>(comm->getRank() * 100 + 25 + i);
      vv2[i] = Teuchos::as<Scalar>(i);
    }
  }

  // create MultiVector in Thyra style (inner)
  Teuchos::RCP<BlockedMultiVector> bvecinner = Teuchos::rcp(new BlockedMultiVector(pp2bm, 2, true));
  TEST_EQUALITY(bvecinner.is_null(),false);
  TEST_EQUALITY(bvecinner->getBlockedMap()->getThyraMode(),true);
  TEST_EQUALITY(bvecinner->getBlockedMap()->getNumMaps(),2);
  bvecinner->setMultiVector(0,smv2,true);
  bvecinner->setMultiVector(1,smv3,true);

  // create MultiVector in Thyra style
  Teuchos::RCP<BlockedMultiVector> bvec = Teuchos::rcp(new BlockedMultiVector(ppbm, 2, true));
  TEST_EQUALITY(bvec.is_null(),false);
  TEST_EQUALITY(bvec->getBlockedMap()->getThyraMode(),true);
  TEST_EQUALITY(bvec->getBlockedMap()->getNumMaps(),2);
  bvec->setMultiVector(0,smv1,true);
  bvec->setMultiVector(1,bvecinner,true);

  // check content of MultiVector
  {
    Teuchos::ArrayRCP<const Scalar> vData  = bvec->getMultiVector(0)->getData(0);
    Teuchos::ArrayRCP<const Scalar> vData2 = bvec->getMultiVector(0)->getData(1);
    for(size_t i=0; i< bvec->getBlockedMap()->getMap(0,true)->getLocalNumElements(); i++) {
      TEST_EQUALITY(vData[i], Teuchos::as<Scalar>(comm->getRank() * 100 + i));
      TEST_EQUALITY(vData2[i], Teuchos::as<Scalar>(i));
      TEST_EQUALITY(bvec->getMultiVector(0)->getMap()->getGlobalElement(i), Teuchos::as<GO>(bvec->getMultiVector(0)->getMap()->getMinGlobalIndex() + i) );
    }
  }

#ifdef HAVE_XPETRA_DEBUG
  // sub vectors are Thyra vectors: currently the implementation throws if we try to
  // access the Xpetra version. We probably can remove the getMultiVector(...,bool)
  // routine...
  TEST_THROW(bvec->getMultiVector(1,false),Xpetra::Exceptions::RuntimeError);
#endif

  Teuchos::RCP<const MultiVector> bvecit = bvec->getMultiVector(1,true);
  TEST_EQUALITY(bvecit.is_null(),false);
  TEST_EQUALITY(bvecit->getMap()->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(bvecit->getMap()->getMaxAllGlobalIndex(), 32);
  {
    Teuchos::RCP<const BlockedMultiVector> bbvecit = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(bvecit);
    Teuchos::ArrayRCP<const Scalar> vData  = bbvecit->getMultiVector(0)->getData(0);
    Teuchos::ArrayRCP<const Scalar> vData2 = bbvecit->getMultiVector(0)->getData(1);
    for(size_t i=0; i< bbvecit->getMultiVector(0)->getMap()->getLocalNumElements(); i++) {
      TEST_EQUALITY(vData[i], Teuchos::as<Scalar>(comm->getRank() * 100 + 10 + i));
      TEST_EQUALITY(vData2[i], Teuchos::as<Scalar>(i));
      TEST_EQUALITY(bbvecit->getBlockedMap()->getMap(0,false)->getGlobalElement(i), map2->getGlobalElement(i) );
      TEST_EQUALITY(bbvecit->getBlockedMap()->getMap(0,true)->getGlobalElement(i), map2->getGlobalElement(i) );
    }
    vData  = bbvecit->getMultiVector(1)->getData(0);
    vData2 = bbvecit->getMultiVector(1)->getData(1);
    for(size_t i=0; i< bbvecit->getMultiVector(1)->getMap()->getLocalNumElements(); i++) {
      TEST_EQUALITY(vData[i], Teuchos::as<Scalar>(comm->getRank() * 100 + 25 + i));
      TEST_EQUALITY(vData2[i], Teuchos::as<Scalar>(i));
      TEST_EQUALITY(bbvecit->getBlockedMap()->getMap(1,false)->getGlobalElement(i), map3->getGlobalElement(i) + 15 );
      TEST_EQUALITY(bbvecit->getBlockedMap()->getMap(1,true)->getGlobalElement(i), map3->getGlobalElement(i) );
    }
  }
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( ThyraBlockedMultiVector, BlockedMapDeepCopy, M, MA, Scalar, LO, GO, Node )
{
#ifdef HAVE_XPETRA_THYRA
  typedef Xpetra::Map<LO, GO, Node> Map;
  typedef Xpetra::BlockedMap<LO, GO, Node> BlockedMap;
  typedef Xpetra::MapFactory<LO, GO, Node> MapFactory;

  typedef Xpetra::ThyraUtils<Scalar,LO,GO,Node> ThyraUtils;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  M testMap(1,0,comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  Teuchos::RCP<const Map> map1 = MapFactory::Build(lib, 10, 0, comm);
  Teuchos::RCP<const Map> map2 = MapFactory::Build(lib, 15, 0, comm);
  Teuchos::RCP<const Map> map3 = MapFactory::Build(lib, 18, 0, comm);

  // create Thyra product vector space
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > vs1 = ThyraUtils::toThyra(map1);
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > vs2 = ThyraUtils::toThyra(map2);
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > vs3 = ThyraUtils::toThyra(map3);

  Teuchos::Array<Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > > vecSpacesInner(2);
  Teuchos::Array<Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > > vecSpacesOuter(2);

  vecSpacesInner[0] = vs2;
  vecSpacesInner[1] = vs3;

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > psInner =
    Thyra::productVectorSpace<Scalar>(vecSpacesInner());

  vecSpacesOuter[0] = vs1;
  vecSpacesOuter[1] = psInner;

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > ps =
    Thyra::productVectorSpace<Scalar>(vecSpacesOuter());

  Teuchos::RCP<const Thyra::ProductVectorSpaceBase<Scalar> > pps =
    Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorSpaceBase<Scalar> >(ps);
  TEST_EQUALITY(pps.is_null(),false);

  Teuchos::RCP<const Map> ppm = ThyraUtils::toXpetra(ps, comm);
  TEST_EQUALITY(ppm.is_null(),false);

  Teuchos::RCP<const BlockedMap> ppbm = Teuchos::rcp_dynamic_cast<const BlockedMap>(ppm);
  TEST_EQUALITY(ppbm.is_null(),false);
  TEST_EQUALITY(ppbm->getThyraMode(),true);

  Teuchos::RCP<const BlockedMap> ppbm2 = Teuchos::rcp(new BlockedMap(*ppbm));
  TEST_EQUALITY(ppbm.is_null(),false);

  TEST_EQUALITY(ppbm->isSameAs(*ppbm2),true);

  TEST_EQUALITY(ppbm->getMap(0,false)->isSameAs(*(ppbm->getMap(0,false))),true);
  TEST_EQUALITY(ppbm->getMap(1,false)->isSameAs(*(ppbm->getMap(1,false))),true);
  TEST_EQUALITY(ppbm->getMap(0,true)->isSameAs(*(ppbm->getMap(0,true))),true);
  TEST_EQUALITY(ppbm->getMap(1,true)->isSameAs(*(ppbm->getMap(1,true))),true);

  ppbm = Teuchos::null;

  TEST_EQUALITY(ppbm2.is_null(), false);
  TEST_EQUALITY(ppbm2->getThyraMode(),true);
  TEST_EQUALITY(ppbm2->getMap(0,false)->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(ppbm2->getMap(0,false)->getMaxAllGlobalIndex(), 9);
  TEST_EQUALITY(ppbm2->getMap(1,false)->getMinAllGlobalIndex(), 10);
  TEST_EQUALITY(ppbm2->getMap(1,false)->getMaxAllGlobalIndex(), 42);

  Teuchos::RCP<const BlockedMap> subMapThyra = Teuchos::rcp_dynamic_cast<const BlockedMap>(ppbm2->getMap(1,true));
  TEST_EQUALITY(subMapThyra.is_null(), false);
  TEST_EQUALITY(subMapThyra->getNumMaps(), 2);
  TEST_EQUALITY(subMapThyra->getFullMap()->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(subMapThyra->getFullMap()->getMaxAllGlobalIndex(), 32);

  Teuchos::RCP<const BlockedMap> subMapXpetra = Teuchos::rcp_dynamic_cast<const BlockedMap>(ppbm2->getMap(1,false));
  TEST_EQUALITY(subMapXpetra.is_null(), false);
  TEST_EQUALITY(subMapXpetra->getNumMaps(), 2);
  TEST_EQUALITY(subMapXpetra->getFullMap()->getMinAllGlobalIndex(), 10);
  TEST_EQUALITY(subMapXpetra->getFullMap()->getMaxAllGlobalIndex(), 42);
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( ThyraBlockedMultiVector, BlockedVectorDeepCopy, M, MA, Scalar, LO, GO, Node )
{
#ifdef HAVE_XPETRA_THYRA
  typedef Xpetra::Map<LO, GO, Node> Map;
  typedef Xpetra::BlockedMap<LO, GO, Node> BlockedMap;
  typedef Xpetra::MapFactory<LO, GO, Node> MapFactory;
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::MultiVectorFactory<Scalar, LO, GO, Node> MultiVectorFactory;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  typedef Xpetra::ThyraUtils<Scalar,LO,GO,Node> ThyraUtils;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  M testMap(1,0,comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  Teuchos::RCP<const Map> map1 = MapFactory::Build(lib, 10, 0, comm);
  Teuchos::RCP<const Map> map2 = MapFactory::Build(lib, 15, 0, comm);
  Teuchos::RCP<const Map> map3 = MapFactory::Build(lib, 18, 0, comm);

  // create Thyra product vector space
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > vs1 = ThyraUtils::toThyra(map1);
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > vs2 = ThyraUtils::toThyra(map2);
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > vs3 = ThyraUtils::toThyra(map3);

  Teuchos::Array<Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > > vecSpacesInner(2);
  Teuchos::Array<Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > > vecSpacesOuter(2);

  vecSpacesInner[0] = vs2;
  vecSpacesInner[1] = vs3;

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > psInner =
    Thyra::productVectorSpace<Scalar>(vecSpacesInner());

  vecSpacesOuter[0] = vs1;
  vecSpacesOuter[1] = psInner;

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > ps =
    Thyra::productVectorSpace<Scalar>(vecSpacesOuter());

  Teuchos::RCP<const Thyra::ProductVectorSpaceBase<Scalar> > pps =
    Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorSpaceBase<Scalar> >(ps);
  TEST_EQUALITY(pps.is_null(),false);

  Teuchos::RCP<const Map> ppm = ThyraUtils::toXpetra(ps, comm);
  TEST_EQUALITY(ppm.is_null(),false);

  Teuchos::RCP<const BlockedMap> ppbm = Teuchos::rcp_dynamic_cast<const BlockedMap>(ppm);
  TEST_EQUALITY(ppbm.is_null(),false);
  TEST_EQUALITY(ppbm->getThyraMode(),true);

  Teuchos::RCP<const Map> pp2m = ppbm->getMap(1,true);
  Teuchos::RCP<const BlockedMap> pp2bm = Teuchos::rcp_dynamic_cast<const BlockedMap>(pp2m);
  TEST_EQUALITY(pp2bm.is_null(),false);
  TEST_EQUALITY(pp2bm->getThyraMode(),true);

  // there are no nested Xpetra style blocked maps!
  Teuchos::RCP<const Map> pp3m = ppbm->getMap(1,false);
  Teuchos::RCP<const BlockedMap> pp3bm = Teuchos::rcp_dynamic_cast<const BlockedMap>(pp3m);
  TEST_EQUALITY(pp3bm.is_null(),false);
  TEST_EQUALITY(pp3bm->getThyraMode(),true);

  // create sub multivectors
  Teuchos::RCP<MultiVector> smv1 = MultiVectorFactory::Build(map1,2,true);
  {
    Teuchos::ArrayRCP< Scalar > vv1 = smv1->getDataNonConst(0);
    Teuchos::ArrayRCP< Scalar > vv2 = smv1->getDataNonConst(1);
    for(LO i = 0; i < Teuchos::as<LO>(smv1->getLocalLength()); ++i) {
      vv1[i] = Teuchos::as<Scalar>(comm->getRank() * 100 + i);
      vv2[i] = Teuchos::as<Scalar>(i);
    }
  }

  Teuchos::RCP<MultiVector> smv2 = MultiVectorFactory::Build(map2,2,true);
  {
    Teuchos::ArrayRCP< Scalar > vv1 = smv2->getDataNonConst(0);
    Teuchos::ArrayRCP< Scalar > vv2 = smv2->getDataNonConst(1);
    for(LO i = 0; i < Teuchos::as<LO>(smv2->getLocalLength()); ++i) {
      vv1[i] = Teuchos::as<Scalar>(comm->getRank() * 100 + 10 + i);
      vv2[i] = Teuchos::as<Scalar>(i);
    }
  }

  Teuchos::RCP<MultiVector> smv3 = MultiVectorFactory::Build(map3,2,true);
  {
    Teuchos::ArrayRCP< Scalar > vv1 = smv3->getDataNonConst(0);
    Teuchos::ArrayRCP< Scalar > vv2 = smv3->getDataNonConst(1);
    for(LO i = 0; i < Teuchos::as<LO>(smv3->getLocalLength()); ++i) {
      vv1[i] = Teuchos::as<Scalar>(comm->getRank() * 100 + 25 + i);
      vv2[i] = Teuchos::as<Scalar>(i);
    }
  }

  // create MultiVector in Thyra style (inner)
  Teuchos::RCP<BlockedMultiVector> bvecinner = Teuchos::rcp(new BlockedMultiVector(pp2bm, 2, true));
  TEST_EQUALITY(bvecinner.is_null(),false);
  TEST_EQUALITY(bvecinner->getBlockedMap()->getThyraMode(),true);
  TEST_EQUALITY(bvecinner->getBlockedMap()->getNumMaps(),2);
  bvecinner->setMultiVector(0,smv2,true);
  bvecinner->setMultiVector(1,smv3,true);

  // create MultiVector in Thyra style
  Teuchos::RCP<BlockedMultiVector> bvec = Teuchos::rcp(new BlockedMultiVector(ppbm, 2, true));
  TEST_EQUALITY(bvec.is_null(),false);
  TEST_EQUALITY(bvec->getBlockedMap()->getThyraMode(),true);
  TEST_EQUALITY(bvec->getBlockedMap()->getNumMaps(),2);
  bvec->setMultiVector(0,smv1,true);
  bvec->setMultiVector(1,bvecinner,true);

  // check content of MultiVector
  {
    Teuchos::ArrayRCP<const Scalar> vData  = bvec->getMultiVector(0)->getData(0);
    Teuchos::ArrayRCP<const Scalar> vData2 = bvec->getMultiVector(0)->getData(1);
    for(size_t i=0; i< bvec->getBlockedMap()->getMap(0,true)->getLocalNumElements(); i++) {
      TEST_EQUALITY(vData[i], Teuchos::as<Scalar>(comm->getRank() * 100 + i));
      TEST_EQUALITY(vData2[i], Teuchos::as<Scalar>(i));
      TEST_EQUALITY(bvec->getMultiVector(0)->getMap()->getGlobalElement(i), Teuchos::as<GO>(bvec->getMultiVector(0)->getMap()->getMinGlobalIndex() + i) );
    }
  }

  //
  Teuchos::RCP<BlockedMultiVector> bvec2 = Teuchos::rcp(new BlockedMultiVector(pp3bm,1));
  *bvec2 = *bvec; // deep copy
  TEST_EQUALITY(bvec2.is_null(),false);
  TEST_EQUALITY(bvec2->getBlockedMap()->isSameAs(*(bvec->getBlockedMap())),true);
  TEST_EQUALITY(bvec2->getNumVectors(),bvec->getNumVectors());

  Teuchos::Array<typename STS::magnitudeType> nn(bvec->getNumVectors());
  Teuchos::Array<typename STS::magnitudeType> nn2(bvec2->getNumVectors());
  TEST_NOTHROW( bvec->norm1(nn) );

  bvec = Teuchos::null;

  TEST_NOTHROW( bvec2->norm1(nn2) );
  for (size_t t = 0; t < bvec2->getNumVectors(); t++) {
    TEST_EQUALITY(nn[t],nn2[t]);
  }
#endif
}

//
// INSTANTIATIONS
//

  #define XPETRA_TPETRA_TYPES( S, LO, GO, N) \
    typedef typename Xpetra::TpetraMap<LO,GO,N> M##LO##GO##N; \
    typedef typename Xpetra::TpetraMultiVector<S,LO,GO,N> MV##S##LO##GO##N;

#ifdef HAVE_XPETRA_EPETRA

  #define XPETRA_EPETRA_TYPES( S, LO, GO, N) \
    typedef typename Xpetra::EpetraMapT<GO,N> M##LO##GO##N; \
    typedef typename Xpetra::EpetraMultiVectorT<GO,N> MV##S##LO##GO##N;

#endif

#define XP_BLOCKEDMULTIVECTOR_INSTANT(S,LO,GO,N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( ThyraBlockedMultiVector, Constructor, M##LO##GO##N , MV##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( ThyraBlockedMultiVector, ConstructorNested, M##LO##GO##N , MV##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( ThyraBlockedMultiVector, BlockedMapDeepCopy, M##LO##GO##N , MV##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( ThyraBlockedMultiVector, BlockedVectorDeepCopy, M##LO##GO##N , MV##S##LO##GO##N, S, LO, GO, N ) \


//TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedMultiVector, ExtractVector, M##LO##GO##N , MV##S##LO##GO##N, S, LO, GO, N )
//TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedMultiVector, ExtractVectorThyra, M##LO##GO##N , MV##S##LO##GO##N, S, LO, GO, N )
//TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedMultiVector, InsertVector, M##LO##GO##N , MV##S##LO##GO##N, S, LO, GO, N )
//TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedMultiVector, InsertVectorThyra, M##LO##GO##N , MV##S##LO##GO##N, S, LO, GO, N )


// List of tests which run only with Tpetra
#define XP_TPETRA_BLOCKEDMULTIVECTOR_INSTANT(S,LO,GO,N)

// List of tests which run only with Epetra
#define XP_EPETRA_BLOCKEDMULTIVECTOR_INSTANT(S,LO,GO,N)

#include <TpetraCore_config.h>
#include <TpetraCore_ETIHelperMacros.h>

TPETRA_ETI_MANGLING_TYPEDEFS()
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( XPETRA_TPETRA_TYPES )
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( XP_TPETRA_BLOCKEDMULTIVECTOR_INSTANT )
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( XP_BLOCKEDMULTIVECTOR_INSTANT )

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
