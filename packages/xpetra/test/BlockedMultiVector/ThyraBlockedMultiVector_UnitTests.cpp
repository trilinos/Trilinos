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
#include <Xpetra_ThyraUtils.hpp>

#include <Thyra_DefaultProductVectorSpace_decl.hpp>

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

  // new stuff

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
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( ThyraBlockedMultiVector, Constructor, M##LO##GO##N , MV##S##LO##GO##N, S, LO, GO, N ) \

//TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedMultiVector, ExtractVector, M##LO##GO##N , MV##S##LO##GO##N, S, LO, GO, N )
//TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedMultiVector, ExtractVectorThyra, M##LO##GO##N , MV##S##LO##GO##N, S, LO, GO, N )
//TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedMultiVector, InsertVector, M##LO##GO##N , MV##S##LO##GO##N, S, LO, GO, N )
//TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedMultiVector, InsertVectorThyra, M##LO##GO##N , MV##S##LO##GO##N, S, LO, GO, N )


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
