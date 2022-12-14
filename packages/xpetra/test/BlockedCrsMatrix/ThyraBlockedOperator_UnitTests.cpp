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
#include <Teuchos_as.hpp>
#include <Teuchos_DefaultComm.hpp>

#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#  include "mpi.h"
#endif
#  include "Epetra_SerialComm.h"

#include <Xpetra_ConfigDefs.hpp>

#ifdef HAVE_XPETRA_EPETRAEXT
// EpetraExt
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#endif

#include "BlockedMatrixTestHelpers.hpp" // handling of Epetra block matrices (SplitMap etc...)

#include <Xpetra_DefaultPlatform.hpp>

#include <Xpetra_Map.hpp>
#include <Xpetra_MapUtils.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_BlockReorderManager.hpp>
#include <Xpetra_ReorderedBlockedCrsMatrix.hpp>
#include <Xpetra_MatrixUtils.hpp>
#include <Xpetra_IO.hpp>

#ifdef HAVE_XPETRA_THYRA
#include <Xpetra_ThyraUtils.hpp>
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultProductMultiVector.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#endif
#include <Xpetra_Exceptions.hpp>

namespace XpetraBlockMatrixTests {

using Xpetra::DefaultPlatform;


using Xpetra::viewLabel_t;

bool testMpi = true;

Teuchos::RCP<const Teuchos::Comm<int> > getDefaultComm()
{
  if (testMpi) {
    return DefaultPlatform::getDefaultPlatform().getComm();
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
}

//
// UNIT TESTS
//


TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( ThyraBlockedOperator, ThyraVectorSpace2XpetraMap_Tpetra, M, MA, Scalar, LO, GO, Node )
{
#ifdef HAVE_XPETRA_THYRA
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // TPetra version
  {
    Teuchos::RCP<const Xpetra::Map<LO,GO,Node> > map = Xpetra::MapFactory<LO,GO,Node>::Build(Xpetra::UseTpetra, 1000, 0, comm);
    TEST_EQUALITY(Teuchos::is_null(map),false);
    Teuchos::RCP<const Tpetra::Map<LO,GO,Node > > tMap = Xpetra::toTpetra(map);
    TEST_EQUALITY(Teuchos::is_null(tMap),false);

    // transform to Thyra...
    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > thyraVectorSpace = Thyra::createVectorSpace<Scalar, LO, GO, Node>(tMap);
    TEST_EQUALITY(Teuchos::is_null(thyraVectorSpace),false);

    // transform back to Xpetra...
    Teuchos::RCP<Xpetra::Map<LO,GO,Node> > xMap = Xpetra::ThyraUtils<Scalar,LO,GO,Node>::toXpetra(thyraVectorSpace,comm);

    TEST_EQUALITY(Teuchos::is_null(xMap),false);

    TEST_EQUALITY(xMap->isCompatible(*map),true);
    TEST_EQUALITY(xMap->isSameAs(*map),true);
  }
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( ThyraBlockedOperator, ThyraVectorSpace2XpetraMap_Epetra, M, MA, Scalar, LO, GO, Node )
{
#ifdef HAVE_XPETRA_THYRA
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // Epetra version
#ifdef HAVE_XPETRA_EPETRA
  {
    Teuchos::RCP<const Xpetra::Map<LO,GO,Node> > map = Xpetra::MapFactory<LO,GO,Node>::Build(Xpetra::UseEpetra, 1000, 0, comm);
    TEST_EQUALITY(Teuchos::is_null(map),false);
    const Epetra_Map ret = Xpetra::toEpetra(map);
    Teuchos::RCP<const Epetra_Map> eMap = Teuchos::rcp(&ret,false);
    TEST_EQUALITY(Teuchos::is_null(eMap),false);

    // transform to Thyra...
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > thyraVectorSpace = Thyra::create_VectorSpace(eMap);
    TEST_EQUALITY(Teuchos::is_null(thyraVectorSpace),false);

    // transform back to Xpetra...
    Teuchos::RCP<Xpetra::Map<LO,GO,Node> > xMap = Xpetra::ThyraUtils<Scalar,LO,GO,Node>::toXpetra(thyraVectorSpace,comm);

    TEST_EQUALITY(Teuchos::is_null(xMap),false);

    TEST_EQUALITY(xMap->isCompatible(*map),true);
    TEST_EQUALITY(xMap->isSameAs(*map),true);
  }
#endif
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( ThyraBlockedOperator, ThyraShrinkMaps, M, MA, Scalar, LO, GO, Node )
{
#ifdef HAVE_XPETRA_THYRA
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  M testMap(1,0,comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  typedef Xpetra::Map<LO, GO, Node> MapClass;
  typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;
  typedef Xpetra::MapUtils<LO, GO, Node> MapUtilsClass;

  // generate non-overlapping map
  Teuchos::Array<GO> myGIDs;
  for(int i = 0; i < 10; i++) {
    myGIDs.push_back(comm->getRank() * 100 + i * 3);
  }
  const Teuchos::RCP<const MapClass> map = MapFactoryClass::Build (lib,Teuchos::OrdinalTraits<GO>::invalid(),myGIDs(),0,comm);

  const Teuchos::RCP<const MapClass> thMap = MapUtilsClass::shrinkMapGIDs(*map,*map);

  TEST_EQUALITY(thMap->getGlobalNumElements() , Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 10));
  TEST_EQUALITY(thMap->getLocalNumElements() , 10);
  TEST_EQUALITY(thMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(thMap->getMaxLocalIndex(), 9);
  TEST_EQUALITY(thMap->getMinGlobalIndex(), comm->getRank() * 10);
  TEST_EQUALITY(thMap->getMaxGlobalIndex(), comm->getRank() * 10 + 9);
  TEST_EQUALITY(thMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(thMap->getMaxAllGlobalIndex(), comm->getSize() * 10 - 1);
  //TEST_EQUALITY(thMap->isContiguous(), true);

  // generate overlapping map
  Teuchos::Array<GO> myovlGIDs;
  if(comm->getRank() > 0) {
    myovlGIDs.push_back((comm->getRank()-1) * 100 + 21);
    myovlGIDs.push_back((comm->getRank()-1) * 100 + 24);
    myovlGIDs.push_back((comm->getRank()-1) * 100 + 27);
  }
  for(int i = 0; i < 10; i++) {
    myovlGIDs.push_back(comm->getRank() * 100 + i * 3);
  }

  const Teuchos::RCP<const MapClass> map2 = MapFactoryClass::Build (lib,Teuchos::OrdinalTraits<GO>::invalid(),myovlGIDs(),0,comm);

  const Teuchos::RCP<const MapClass> thMap2 = MapUtilsClass::shrinkMapGIDs(*map2,*map);

  TEST_EQUALITY(thMap2->getMinGlobalIndex() , std::max(0,comm->getRank() * 10 - 3));
  TEST_EQUALITY(thMap2->getMaxGlobalIndex() , comm->getRank() * 10 + 9);
  TEST_EQUALITY(thMap2->getLocalNumElements() , (comm->getRank() > 0) ? 13 : 10 );
  TEST_EQUALITY(thMap2->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(thMap2->getMaxAllGlobalIndex(), comm->getSize() * 10 - 1);
  TEST_EQUALITY(thMap2->getGlobalNumElements() , Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 13 - 3));
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( ThyraBlockedOperator, ThyraOperator2XpetraCrsMat, M, MA, Scalar, LO, GO, Node )
{
#ifdef HAVE_XPETRA_THYRA
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  M testMap(1,0,comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;

  // generate the matrix
  LO nEle = 63;
  const Teuchos::RCP<const Xpetra::Map<LO,GO,Node> > map = MapFactoryClass::Build(lib, nEle, 0, comm);

  Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > matrix =
               Xpetra::CrsMatrixFactory<Scalar,LO,GO,Node>::Build(map, 10);

  LO NumMyElements = map->getLocalNumElements();
  Teuchos::ArrayView<const GO> MyGlobalElements = map->getLocalElementList();

  for (LO i = 0; i < NumMyElements; ++i) {
      matrix->insertGlobalValues(MyGlobalElements[i],
                              Teuchos::tuple<GO>(MyGlobalElements[i]),
                              Teuchos::tuple<Scalar>(1.0) );
  }

  matrix->fillComplete();

  // create a Thyra operator from Xpetra::CrsMatrix
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > thyraOp =
      Xpetra::ThyraUtils<Scalar,LO,GO,Node>::toThyra(matrix);

  // transform Thyra operator 2 Xpetra::CrsMatrix
  Teuchos::RCP<const Xpetra::Matrix<Scalar, LO, GO, Node> > xMat =
      Xpetra::ThyraUtils<Scalar,LO,GO,Node>::toXpetra(thyraOp);
  Teuchos::RCP<const Xpetra::CrsMatrix<Scalar, LO, GO, Node> > xCrsMat =
    Teuchos::rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LO, GO, Node> >(xMat)->getCrsMatrix();
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xCrsMat));

  TEST_EQUALITY(xCrsMat->getFrobeniusNorm()   ,matrix->getFrobeniusNorm());
  TEST_EQUALITY(xCrsMat->getGlobalNumRows()   ,matrix->getGlobalNumRows());
  TEST_EQUALITY(xCrsMat->getGlobalNumCols()   ,matrix->getGlobalNumCols());
  TEST_EQUALITY(xCrsMat->getLocalNumRows()     ,matrix->getLocalNumRows()  );
  TEST_EQUALITY(xCrsMat->getGlobalNumEntries(),matrix->getGlobalNumEntries());
  TEST_EQUALITY(xCrsMat->getLocalNumEntries()  ,matrix->getLocalNumEntries());

#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( ThyraBlockedOperator, ThyraBlockedOperator2XpetraBlockedCrsMat, M, MA, Scalar, LO, GO, Node )
{
#ifdef HAVE_XPETRA_THYRA
  typedef Xpetra::Map<LO, GO, Node> MapClass;
  typedef Xpetra::StridedMap<LO, GO, Node> StridedMapClass;
  typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;
  typedef Xpetra::StridedMapFactory<LO, GO, Node> StridedMapFactoryClass;
  typedef Xpetra::MapExtractor<Scalar,LO,GO,Node> MapExtractorClass;
  typedef Xpetra::MapExtractorFactory<Scalar,LO,GO,Node> MapExtractorFactoryClass;
  typedef Xpetra::Matrix<Scalar, LO, GO, Node> MatrixClass;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  M testMap(1,0,comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  const Teuchos::RCP<const MapClass> pointmap = MapFactoryClass::Build(lib, 12, 0, comm);

  // generate local maps for loading matrices
  Teuchos::Array<GO> velgidvec; // global strided maps
  Teuchos::Array<GO> pregidvec;
  Teuchos::Array<GO> fullgidvec; // full global map
  for (LO i=0; i<Teuchos::as<LO>(pointmap->getLocalNumElements()); i++)
  {
    // loop over all local ids in pointmap

    // get corresponding global id
    GO gid = pointmap->getGlobalElement(i);

    // store global strided gids
    velgidvec.push_back(3*gid);
    velgidvec.push_back(3*gid+1);
    pregidvec.push_back(3*gid+2);

    // gid for full map
    fullgidvec.push_back(3*gid);
    fullgidvec.push_back(3*gid+1);
    fullgidvec.push_back(3*gid+2);
  }

  std::vector<size_t> stridingInfo;
  stridingInfo.push_back(2);
  stridingInfo.push_back(1);

  const Teuchos::RCP<const StridedMapClass> velmap = StridedMapFactoryClass::Build (lib,
      Teuchos::OrdinalTraits<GO>::invalid(),
      velgidvec(),
      0,
      stridingInfo,
      comm, 0);
  const Teuchos::RCP<const StridedMapClass> premap = StridedMapFactoryClass::Build (lib,
      Teuchos::OrdinalTraits<GO>::invalid(),
      pregidvec(),
      0,
      stridingInfo,
      comm,1);
  const Teuchos::RCP<const StridedMapClass> fullmap = StridedMapFactoryClass::Build (lib,
      Teuchos::OrdinalTraits<GO>::invalid(),
      fullgidvec(),
      0,
      stridingInfo,
      comm,-1);

  std::string tname = typeid(Scalar).name();
  if(tname.find("complex")!=std::string::npos) {
    std::cout << "Skip test for scalar=" << tname << std::endl;
    return;
  }

  Teuchos::RCP<MatrixClass> A = Xpetra::IO<Scalar,LO,GO,Node>::Read("A.mat", fullmap->getMap());

  std::vector<Teuchos::RCP<const MapClass> > xmaps;
  xmaps.push_back(velmap);
  xmaps.push_back(premap);

  Teuchos::RCP<const MapExtractorClass> map_extractor = MapExtractorFactoryClass::Build(fullmap,xmaps);

  Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> > bOp =
      Xpetra::MatrixUtils<Scalar, LO, GO, Node>::SplitMatrix(*A,map_extractor,map_extractor);

  bOp->fillComplete();
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(bOp));

  // create Thyra operator
  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > thOp =
      Xpetra::ThyraUtils<Scalar,LO,GO,Node>::toThyra(bOp);
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thOp));

  Teuchos::RCP<Thyra::BlockedLinearOpBase<Scalar> > thbOp =
      Teuchos::rcp_dynamic_cast<Thyra::BlockedLinearOpBase<Scalar> >(thOp);
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thbOp));

  Teuchos::RCP<const Thyra::ProductVectorSpaceBase<Scalar> > productRange  = thbOp->productRange();
  Teuchos::RCP<const Thyra::ProductVectorSpaceBase<Scalar> > productDomain = thbOp->productDomain();

  TEST_EQUALITY(productRange->numBlocks()  ,2);
  TEST_EQUALITY(productDomain->numBlocks() ,2);
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productRange->dim()) ,fullmap->getGlobalNumElements());
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productDomain->dim()) ,fullmap->getGlobalNumElements());
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productRange->getBlock(0)->dim())  ,velmap->getGlobalNumElements());
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productDomain->getBlock(0)->dim()) ,velmap->getGlobalNumElements());
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productRange->getBlock(1)->dim())  ,premap->getGlobalNumElements());
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productDomain->getBlock(1)->dim()) ,premap->getGlobalNumElements());
#endif // end HAVE_XPETRA_THYRA
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( ThyraBlockedOperator, XpetraBlockedCrsMatConstructor, M, MA, Scalar, LO, GO, Node )
{
#ifdef HAVE_XPETRA_THYRA

  typedef Xpetra::Map<LO, GO, Node> MapClass;
  typedef Xpetra::StridedMap<LO, GO, Node> StridedMapClass;
  typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;
  typedef Xpetra::StridedMapFactory<LO, GO, Node> StridedMapFactoryClass;
  typedef Xpetra::MapExtractor<Scalar,LO,GO,Node> MapExtractorClass;
  typedef Xpetra::MapExtractorFactory<Scalar,LO,GO,Node> MapExtractorFactoryClass;
  typedef Xpetra::Matrix<Scalar, LO, GO, Node> MatrixClass;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  M testMap(1,0,comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  const Teuchos::RCP<const MapClass> pointmap = MapFactoryClass::Build(lib, 12, 0, comm);

  // generate local maps for loading matrices
  Teuchos::Array<GO> velgidvec; // global strided maps
  Teuchos::Array<GO> pregidvec;
  Teuchos::Array<GO> fullgidvec; // full global map
  for (LO i=0; i<Teuchos::as<LO>(pointmap->getLocalNumElements()); i++)
  {
    // loop over all local ids in pointmap

    // get corresponding global id
    GO gid = pointmap->getGlobalElement(i);

    // store global strided gids
    velgidvec.push_back(3*gid);
    velgidvec.push_back(3*gid+1);
    pregidvec.push_back(3*gid+2);

    // gid for full map
    fullgidvec.push_back(3*gid);
    fullgidvec.push_back(3*gid+1);
    fullgidvec.push_back(3*gid+2);
  }

  std::vector<size_t> stridingInfo;
  stridingInfo.push_back(2);
  stridingInfo.push_back(1);

  const Teuchos::RCP<const StridedMapClass> velmap = StridedMapFactoryClass::Build (lib,
      Teuchos::OrdinalTraits<GO>::invalid(),
      velgidvec(),
      0,
      stridingInfo,
      comm, 0);
  const Teuchos::RCP<const StridedMapClass> premap = StridedMapFactoryClass::Build (lib,
      Teuchos::OrdinalTraits<GO>::invalid(),
      pregidvec(),
      0,
      stridingInfo,
      comm,1);
  const Teuchos::RCP<const StridedMapClass> fullmap = StridedMapFactoryClass::Build (lib,
      Teuchos::OrdinalTraits<GO>::invalid(),
      fullgidvec(),
      0,
      stridingInfo,
      comm,-1);

  std::string tname = typeid(Scalar).name();
  if(tname.find("complex")!=std::string::npos) {
    std::cout << "Skip test for scalar=" << tname << std::endl;
    return;
  }

  Teuchos::RCP<MatrixClass> A = Xpetra::IO<Scalar,LO,GO,Node>::Read("A.mat", fullmap->getMap());

  std::vector<Teuchos::RCP<const MapClass> > xmaps;
  xmaps.push_back(velmap);
  xmaps.push_back(premap);

  Teuchos::RCP<const MapExtractorClass> map_extractor = MapExtractorFactoryClass::Build(fullmap,xmaps);

  Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> > bOp =
      Xpetra::MatrixUtils<Scalar, LO, GO, Node>::SplitMatrix(*A,map_extractor,map_extractor);

  bOp->fillComplete();
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(bOp));
  TEST_EQUALITY(bOp->getRangeMapExtractor()->getThyraMode(),false);
  TEST_EQUALITY(bOp->getDomainMapExtractor()->getThyraMode(),false);

  // create Thyra operator
  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > thOp =
      Xpetra::ThyraUtils<Scalar,LO,GO,Node>::toThyra(bOp);
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thOp));

  Teuchos::RCP<Thyra::BlockedLinearOpBase<Scalar> > thbOp =
      Teuchos::rcp_dynamic_cast<Thyra::BlockedLinearOpBase<Scalar> >(thOp);
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thbOp));

  Teuchos::RCP<const Thyra::ProductVectorSpaceBase<Scalar> > productRange  = thbOp->productRange();
  Teuchos::RCP<const Thyra::ProductVectorSpaceBase<Scalar> > productDomain = thbOp->productDomain();

  TEST_EQUALITY(productRange->numBlocks()  ,2);
  TEST_EQUALITY(productDomain->numBlocks() ,2);
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productRange->dim())  ,fullmap->getGlobalNumElements());
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productDomain->dim()) ,fullmap->getGlobalNumElements());
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productRange->getBlock(0)->dim())  ,velmap->getGlobalNumElements());
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productDomain->getBlock(0)->dim()) ,velmap->getGlobalNumElements());
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productRange->getBlock(1)->dim())  ,premap->getGlobalNumElements());
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productDomain->getBlock(1)->dim()) ,premap->getGlobalNumElements());

  //construct a Xpetra::BlockedCrsMatrix object
  Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> > bOp2 =
      Teuchos::rcp(new Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node>(thbOp,comm));

  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(bOp2));
  TEST_EQUALITY(bOp2->getGlobalNumRows(),bOp->getGlobalNumRows());
  TEST_EQUALITY(bOp2->getGlobalNumCols(),bOp->getGlobalNumCols());
  TEST_EQUALITY(bOp2->getLocalNumRows(),bOp->getLocalNumRows());
  TEST_EQUALITY(bOp2->getGlobalNumEntries(),bOp->getGlobalNumEntries());
  TEST_EQUALITY(bOp2->getLocalNumEntries(),bOp->getLocalNumEntries());
  TEST_EQUALITY(bOp2->isFillComplete(),bOp->isFillComplete());
  TEST_EQUALITY(bOp2->getRangeMapExtractor()->getThyraMode(),false);
  TEST_EQUALITY(bOp2->getDomainMapExtractor()->getThyraMode(),false);
  TEST_EQUALITY(bOp2->getDomainMap()->isCompatible(*bOp->getDomainMap()),true);
  TEST_EQUALITY(bOp2->getRangeMap()->isCompatible(*bOp->getRangeMap()),true);
  TEST_EQUALITY(bOp2->getDomainMap(0)->isCompatible(*bOp->getDomainMap(0)),true);
  TEST_EQUALITY(bOp2->getRangeMap(0)->isCompatible(*bOp->getRangeMap(0)),true);
  TEST_EQUALITY(bOp2->getDomainMap(1)->isCompatible(*bOp->getDomainMap(1)),true);
  TEST_EQUALITY(bOp2->getRangeMap(1)->isCompatible(*bOp->getRangeMap(1)),true);
  TEST_EQUALITY(bOp2->getDomainMap()->isSameAs(*bOp->getDomainMap()),true);
  TEST_EQUALITY(bOp2->getRangeMap()->isSameAs(*bOp->getRangeMap()),true);
  TEST_EQUALITY(bOp2->getDomainMap(0)->isSameAs(*bOp->getDomainMap(0)),true);
  TEST_EQUALITY(bOp2->getRangeMap(0)->isSameAs(*bOp->getRangeMap(0)),true);
  TEST_EQUALITY(bOp2->getDomainMap(1)->isSameAs(*bOp->getDomainMap(1)),true);
  TEST_EQUALITY(bOp2->getRangeMap(1)->isSameAs(*bOp->getRangeMap(1)),true);

#endif // end HAVE_XPETRA_THYRA
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( ThyraBlockedOperator, SplitMatrixForThyra, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::Map<LO, GO, Node> MapClass;
  typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;
  typedef Xpetra::Vector<Scalar, LO, GO, Node> VectorClass;
  typedef Xpetra::VectorFactory<Scalar, LO, GO, Node> VectorFactoryClass;
  typedef Xpetra::MapExtractor<Scalar,LO,GO,Node> MapExtractorClass;
  typedef Xpetra::MapExtractorFactory<Scalar,LO,GO,Node> MapExtractorFactoryClass;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  M testMap(1,0,comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  // generate problem
  GO nEle = 63;
  const Teuchos::RCP<const MapClass> map = MapFactoryClass::Build(lib, nEle, 0, comm);

  LO NumMyElements = map->getLocalNumElements();
  GO NumGlobalElements = map->getGlobalNumElements();
  Teuchos::ArrayView<const GO> MyGlobalElements = map->getLocalElementList();

  Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > A =
      Xpetra::CrsMatrixFactory<Scalar,LO,GO,Node>::Build(map, 3);
  TEUCHOS_TEST_FOR_EXCEPTION(A->isFillComplete() == true || A->isFillActive() == false, std::runtime_error, "");

  for (LO i = 0; i < NumMyElements; i++) {
     if (MyGlobalElements[i] == 0) {
       A->insertGlobalValues(MyGlobalElements[i],
                             Teuchos::tuple<GO>(MyGlobalElements[i], MyGlobalElements[i] +1),
                             Teuchos::tuple<Scalar> (Teuchos::as<Scalar>(i)*STS::one(), -1.0));
     }
     else if (MyGlobalElements[i] == NumGlobalElements - 1) {
       A->insertGlobalValues(MyGlobalElements[i],
                             Teuchos::tuple<GO>(MyGlobalElements[i] -1, MyGlobalElements[i]),
                             Teuchos::tuple<Scalar> (-1.0, Teuchos::as<Scalar>(i)*STS::one()));
     }
     else {
       A->insertGlobalValues(MyGlobalElements[i],
                             Teuchos::tuple<GO>(MyGlobalElements[i] -1, MyGlobalElements[i], MyGlobalElements[i] +1),
                             Teuchos::tuple<Scalar> (-1.0, Teuchos::as<Scalar>(i)*STS::one(), -1.0));
     }
  }

  A->fillComplete();
  TEUCHOS_TEST_FOR_EXCEPTION(A->isFillComplete() == false || A->isFillActive() == true, std::runtime_error, "");

  Teuchos::RCP<Xpetra::Matrix<Scalar, LO, GO, Node> > mat =
      Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LO, GO, Node>(A));

  Teuchos::Array<GO> gids1;
  Teuchos::Array<GO> gids2;
  for(LO i=0; i<NumMyElements; i++) {
    if(i % 3 < 2)
      gids1.push_back(map->getGlobalElement(i));
    else
      gids2.push_back(map->getGlobalElement(i));
  }

  const Teuchos::RCP<const MapClass> map1 = MapFactoryClass::Build (lib,
      Teuchos::OrdinalTraits<GO>::invalid(),
      gids1.view(0,gids1.size()),
      0,
      comm);
  const Teuchos::RCP<const MapClass> map2 = MapFactoryClass::Build (lib,
      Teuchos::OrdinalTraits<GO>::invalid(),
      gids2.view(0,gids2.size()),
      0,
      comm);

  std::vector<Teuchos::RCP<const MapClass> > xmaps;
  xmaps.push_back(map1);
  xmaps.push_back(map2);

  Teuchos::RCP<const MapExtractorClass> map_extractor = MapExtractorFactoryClass::Build(map,xmaps);

  Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> > bOp =
      Xpetra::MatrixUtils<Scalar, LO, GO, Node>::SplitMatrix(*mat,map_extractor,map_extractor,Teuchos::null,true);

  // build gloabl vector with one entries
  Teuchos::RCP<VectorClass> ones_A   = VectorFactoryClass::Build(map, true);
  Teuchos::RCP<VectorClass> exp  = VectorFactoryClass::Build(map, true);
  Teuchos::RCP<VectorClass> ones_bOp = VectorFactoryClass::Build(bOp->getRangeMap(), true);
  Teuchos::RCP<VectorClass> res  = VectorFactoryClass::Build(bOp->getRangeMap(), true);
  ones_A->putScalar(STS::one());
  ones_bOp->putScalar(STS::one());

  A->apply(*ones_A, *exp);
  bOp->apply(*ones_bOp, *res);

  TEST_EQUALITY(res->norm1(),exp->norm1());
  TEST_EQUALITY(res->normInf(),exp->normInf());
  TEST_EQUALITY(bOp->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(bOp->getDomainMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(bOp->getRangeMap(0)->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(bOp->getRangeMap(1)->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(bOp->getDomainMap(0)->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(bOp->getDomainMap(1)->getMinAllGlobalIndex(), 0);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(ThyraBlockedOperator, NestedBlockOperator2ThyraBlockedCrsMat, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::Map< LO, GO, Node> MapClass;
  typedef Xpetra::Vector<Scalar, LO, GO, Node> VectorClass;
  typedef Xpetra::VectorFactory<Scalar, LO, GO, Node> VectorFactoryClass;
  typedef Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> BlockedCrsMatrixClass;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 8;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrixThyra<Scalar,LO,GO,Node,M>(noBlocks, *comm);


  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 6 [3 2] ]");

  Teuchos::RCP<const BlockedCrsMatrixClass> brop =
      Teuchos::rcp_dynamic_cast<const BlockedCrsMatrixClass>(buildReorderedBlockedCrsMatrix(brm, bop));

  Teuchos::RCP<const MapClass> fullmap = brop->getRowMap();

  TEST_EQUALITY(brop->Rows(),2);
  TEST_EQUALITY(brop->Cols(),2);
  TEST_EQUALITY(brop->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop->getDomainMapExtractor()->getThyraMode(), true);

  Teuchos::RCP<VectorClass> ones_A   = VectorFactoryClass::Build(fullmap, true);
  ones_A->putScalar(STS::one());

  Teuchos::RCP<VectorClass> exp   = VectorFactoryClass::Build(fullmap, true);
  exp->putScalar(STS::zero());


#ifdef HAVE_XPETRA_THYRA
  // create Thyra operator
  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > thOp =
      Xpetra::ThyraUtils<Scalar,LO,GO,Node>::toThyra(Teuchos::rcp_const_cast<BlockedCrsMatrixClass>(brop));
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thOp));

  Teuchos::RCP<Thyra::BlockedLinearOpBase<Scalar> > thbOp =
      Teuchos::rcp_dynamic_cast<Thyra::BlockedLinearOpBase<Scalar> >(thOp);
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thbOp));

  Teuchos::RCP<const Thyra::ProductVectorSpaceBase<Scalar> > productRange  = thbOp->productRange();
  Teuchos::RCP<const Thyra::ProductVectorSpaceBase<Scalar> > productDomain = thbOp->productDomain();

  TEST_EQUALITY(productRange->numBlocks()  ,2);
  TEST_EQUALITY(productDomain->numBlocks() ,2);
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productRange->dim())  ,fullmap->getGlobalNumElements());
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productDomain->dim()) ,fullmap->getGlobalNumElements());
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productRange->getBlock(0)->dim())  ,Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 160));
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productDomain->getBlock(0)->dim()) ,Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 160));
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productRange->getBlock(1)->dim())  ,Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 30));
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productDomain->getBlock(1)->dim()) ,Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 30));

  Teuchos::RCP<const Thyra::BlockedLinearOpBase<Scalar> > thbOp11 =
      Teuchos::rcp_dynamic_cast<const Thyra::BlockedLinearOpBase<Scalar> >(thbOp->getBlock(1,1));
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thbOp11));

  Teuchos::RCP<const Thyra::ProductVectorSpaceBase<Scalar> > productRange11  = thbOp11->productRange();
  Teuchos::RCP<const Thyra::ProductVectorSpaceBase<Scalar> > productDomain11 = thbOp11->productDomain();
  TEST_EQUALITY(productRange11->numBlocks()  ,2);
  TEST_EQUALITY(productDomain11->numBlocks() ,2);
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productRange11->dim())  ,Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 30));
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productDomain11->dim()) ,Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 30));
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productRange11->getBlock(0)->dim())  ,Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 20));
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productDomain11->getBlock(0)->dim()) ,Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 20));
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productRange11->getBlock(1)->dim())  ,Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 10));
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productDomain11->getBlock(1)->dim()) ,Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 10));

#if 0
  // create Thyra product vector space
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > vs1 =
      Xpetra::ThyraUtils<Scalar,LO,GO,Node>::toThyra(brop->getRangeMapExtractor()->getMap(0,true));
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > vs2 =
      Xpetra::ThyraUtils<Scalar,LO,GO,Node>::toThyra(brop->getRangeMapExtractor()->getMap(1,true));
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > vs3 =
      Xpetra::ThyraUtils<Scalar,LO,GO,Node>::toThyra(brop->getRangeMapExtractor()->getMap(2,true));

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

  Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > X = Thyra::createMembers(ps,1);

  Thyra::assign<Scalar>(X.ptr(), Teuchos::as<Scalar>(0.0));

  Xpetra::ThyraUtils<Scalar,LO,GO,Node>::updateThyra(ones_A, brop->getDomainMapExtractor(), Teuchos::rcp_dynamic_cast<Thyra::ProductMultiVectorBase<Scalar> >(X));

  Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > B = Thyra::createMembers(ps,1);
  Thyra::assign<Scalar>(B.ptr(), Teuchos::as<Scalar>(0.0));

  thbOp->apply( Thyra::NOTRANS, *X, B.ptr(), STS::one(), STS::zero());

  Teuchos::RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> > xres =
        Xpetra::ThyraUtils<Scalar,LO,GO,Node>::toXpetra(B, comm);

  Teuchos::ArrayRCP<const Scalar> xdata = xres->getData(0);
  bool bCheck = true;
  for(int i=0; i<xres->getLocalLength(); i++) {
    if (i>=0  && i< 160) { if(xdata[i] != (Scalar) 7.0){ std::cout << i << " " << xdata[i] << std::endl; bCheck = false; }}
    if (i>=160 && i< 180) { if(xdata[i] != (Scalar) 4.0){std::cout << i << " " << xdata[i] << std::endl; bCheck = false; }}
    if (i>=180 && i< 190) { if(xdata[i] != (Scalar) 3.0){std::cout << i << " " << xdata[i] << std::endl; bCheck = false; }}
  }
  TEST_EQUALITY(bCheck, true);

  brop->apply(*ones_A, *exp);
  exp->update(-STS::one(),*xres,STS::one());
  TEST_EQUALITY(exp->norm1(), STS::zero());
#endif // not working (no nested Thyra operators!)
#endif // HAVE_XPETRA_THYRA
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(ThyraBlockedOperator, ReorderBlockOperator2ThyraBlockedCrsMat, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::Map< LO, GO, Node> MapClass;
  typedef Xpetra::Vector<Scalar, LO, GO, Node> VectorClass;
  typedef Xpetra::VectorFactory<Scalar, LO, GO, Node> VectorFactoryClass;
  typedef Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> BlockedCrsMatrixClass;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 8;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrixThyra<Scalar,LO,GO,Node,M>(noBlocks, *comm);


  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 6 3 2]");

  Teuchos::RCP<const BlockedCrsMatrixClass> brop =
      Teuchos::rcp_dynamic_cast<const BlockedCrsMatrixClass>(buildReorderedBlockedCrsMatrix(brm, bop));

  Teuchos::RCP<const MapClass> fullmap = brop->getRowMap();

  TEST_EQUALITY(brop->Rows(),3);
  TEST_EQUALITY(brop->Cols(),3);
  TEST_EQUALITY(brop->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop->getDomainMapExtractor()->getThyraMode(), true);

  Teuchos::RCP<VectorClass> ones_A   = VectorFactoryClass::Build(fullmap, true);
  ones_A->putScalar(STS::one());

  Teuchos::RCP<VectorClass> exp   = VectorFactoryClass::Build(fullmap, true);
  exp->putScalar(STS::zero());


#ifdef HAVE_XPETRA_THYRA
  // create Thyra operator
  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > thOp =
      Xpetra::ThyraUtils<Scalar,LO,GO,Node>::toThyra(Teuchos::rcp_const_cast<BlockedCrsMatrixClass>(brop));
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thOp));

  Teuchos::RCP<Thyra::BlockedLinearOpBase<Scalar> > thbOp =
      Teuchos::rcp_dynamic_cast<Thyra::BlockedLinearOpBase<Scalar> >(thOp);
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thbOp));

  Teuchos::RCP<const Thyra::ProductVectorSpaceBase<Scalar> > productRange  = thbOp->productRange();
  Teuchos::RCP<const Thyra::ProductVectorSpaceBase<Scalar> > productDomain = thbOp->productDomain();

  TEST_EQUALITY(productRange->numBlocks()  ,3);
  TEST_EQUALITY(productDomain->numBlocks() ,3);
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productRange->dim())  ,fullmap->getGlobalNumElements());
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productDomain->dim()) ,fullmap->getGlobalNumElements());
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productRange->getBlock(0)->dim())  ,Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 160));
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productDomain->getBlock(0)->dim()) ,Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 160));
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productRange->getBlock(1)->dim())  ,Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 20));
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productDomain->getBlock(1)->dim()) ,Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 20));
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productRange->getBlock(2)->dim())  ,Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 10));
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productDomain->getBlock(2)->dim()) ,Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 10));

  // create Thyra product vector space
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > vs1 =
      Xpetra::ThyraUtils<Scalar,LO,GO,Node>::toThyra(brop->getRangeMapExtractor()->getMap(0,true));
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > vs2 =
      Xpetra::ThyraUtils<Scalar,LO,GO,Node>::toThyra(brop->getRangeMapExtractor()->getMap(1,true));
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > vs3 =
      Xpetra::ThyraUtils<Scalar,LO,GO,Node>::toThyra(brop->getRangeMapExtractor()->getMap(2,true));

  Teuchos::Array<Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > > vecSpaces(3);
  vecSpaces[0] = vs1;
  vecSpaces[1] = vs2;
  vecSpaces[2] = vs3;

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > ps =
    Thyra::productVectorSpace<Scalar>(vecSpaces());

  Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > X = Thyra::createMembers(ps,1);

  Thyra::assign<Scalar>(X.ptr(), Teuchos::as<Scalar>(0.0));

  Xpetra::ThyraUtils<Scalar,LO,GO,Node>::updateThyra(ones_A, brop->getDomainMapExtractor(), Teuchos::rcp_dynamic_cast<Thyra::ProductMultiVectorBase<Scalar> >(X));

  Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > B = Thyra::createMembers(ps,1);
  Thyra::assign<Scalar>(B.ptr(), Teuchos::as<Scalar>(0.0));

  thbOp->apply( Thyra::NOTRANS, *X, B.ptr(), STS::one(), STS::zero());

  Teuchos::RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> > xres =
        Xpetra::ThyraUtils<Scalar,LO,GO,Node>::toXpetra(B, comm);

  Teuchos::RCP<Xpetra::BlockedMultiVector<Scalar,LO,GO,Node> > xresb =
      Teuchos::rcp_dynamic_cast<Xpetra::BlockedMultiVector<Scalar,LO,GO,Node>>(xres);
  TEST_EQUALITY(xresb.is_null(), false);

  Teuchos::RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> > xresmerged = xresb->Merge();

  Teuchos::ArrayRCP<const Scalar> xdata = xresmerged->getData(0);
  bool bCheck = true;
  for(int i=0; i<Teuchos::as<int>(xresmerged->getLocalLength()); i++) {
    if (i>=0  && i< 160) { if(xdata[i] != (Scalar) 7.0){ std::cout << i << " " << xdata[i] << std::endl; bCheck = false; }}
    if (i>=160 && i< 180) { if(xdata[i] != (Scalar) 4.0){std::cout << i << " " << xdata[i] << std::endl; bCheck = false; }}
    if (i>=180 && i< 190) { if(xdata[i] != (Scalar) 3.0){std::cout << i << " " << xdata[i] << std::endl; bCheck = false; }}
  }
  TEST_EQUALITY(bCheck, true);

  brop->apply(*ones_A, *exp);
  exp->update(-STS::one(),*xresmerged,STS::one());
  TEST_EQUALITY(exp->norm1(), STS::zero());
#endif // HAVE_XPETRA_THYRA
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( ThyraBlockedOperator, ReadWriteMatrixMatrixMarket, M, MA, Scalar, LO, GO, Node )
{
#ifdef HAVE_XPETRA_THYRA
  typedef Xpetra::Map<LO, GO, Node> MapClass;
  typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;
  typedef Xpetra::Vector<Scalar, LO, GO, Node> VectorClass;
  typedef Xpetra::VectorFactory<Scalar, LO, GO, Node> VectorFactoryClass;
  typedef Xpetra::MapExtractor<Scalar,LO,GO,Node> MapExtractorClass;
  typedef Xpetra::MapExtractorFactory<Scalar,LO,GO,Node> MapExtractorFactoryClass;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  M testMap(1,0,comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  // generate problem
  GO nEle = 63;
  const Teuchos::RCP<const MapClass> map = MapFactoryClass::Build(lib, nEle, 0, comm);

  LO NumMyElements = map->getLocalNumElements();
  GO NumGlobalElements = map->getGlobalNumElements();
  Teuchos::ArrayView<const GO> MyGlobalElements = map->getLocalElementList();

  Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > A =
      Xpetra::CrsMatrixFactory<Scalar,LO,GO,Node>::Build(map, 3);
  TEUCHOS_TEST_FOR_EXCEPTION(A->isFillComplete() == true || A->isFillActive() == false, std::runtime_error, "");

  for (LO i = 0; i < NumMyElements; i++) {
     if (MyGlobalElements[i] == 0) {
       A->insertGlobalValues(MyGlobalElements[i],
                             Teuchos::tuple<GO>(MyGlobalElements[i], MyGlobalElements[i] +1),
                             Teuchos::tuple<Scalar> (Teuchos::as<Scalar>(i)*STS::one(), -1.0));
     }
     else if (MyGlobalElements[i] == NumGlobalElements - 1) {
       A->insertGlobalValues(MyGlobalElements[i],
                             Teuchos::tuple<GO>(MyGlobalElements[i] -1, MyGlobalElements[i]),
                             Teuchos::tuple<Scalar> (-1.0, Teuchos::as<Scalar>(i)*STS::one()));
     }
     else {
       A->insertGlobalValues(MyGlobalElements[i],
                             Teuchos::tuple<GO>(MyGlobalElements[i] -1, MyGlobalElements[i], MyGlobalElements[i] +1),
                             Teuchos::tuple<Scalar> (-1.0, Teuchos::as<Scalar>(i)*STS::one(), -1.0));
     }
  }

  A->fillComplete();
  TEUCHOS_TEST_FOR_EXCEPTION(A->isFillComplete() == false || A->isFillActive() == true, std::runtime_error, "");

  Teuchos::RCP<Xpetra::Matrix<Scalar, LO, GO, Node> > mat =
      Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LO, GO, Node>(A));

  Teuchos::Array<GO> gids1;
  Teuchos::Array<GO> gids2;
  for(LO i=0; i<NumMyElements; i++) {
    if(i % 3 < 2)
      gids1.push_back(map->getGlobalElement(i));
    else
      gids2.push_back(map->getGlobalElement(i));
  }

  const Teuchos::RCP<const MapClass> map1 = MapFactoryClass::Build (lib,
      Teuchos::OrdinalTraits<GO>::invalid(),
      gids1.view(0,gids1.size()),
      0,
      comm);
  const Teuchos::RCP<const MapClass> map2 = MapFactoryClass::Build (lib,
      Teuchos::OrdinalTraits<GO>::invalid(),
      gids2.view(0,gids2.size()),
      0,
      comm);

  std::vector<Teuchos::RCP<const MapClass> > xmaps;
  xmaps.push_back(map1);
  xmaps.push_back(map2);

  Teuchos::RCP<const MapExtractorClass> map_extractor = MapExtractorFactoryClass::Build(map,xmaps);

  Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> > bMat =
      Xpetra::MatrixUtils<Scalar, LO, GO, Node>::SplitMatrix(*mat,map_extractor,map_extractor,Teuchos::null,true);


  // Write matrices out, read fine A back in, and check that the read was ok
  // by using a matvec with a random vector.
  // JJH: 22-Feb-2016 Append scalar type to file name. The theory is that for dashboard
  //      tests with multiple Scalar instantiations of this test, a test with Scalar type
  //      A could try to read in the results of the test with Scalar type B, simply because
  //      the test with type B overwrote A's output matrix file.  A better solution would be
  //      to write to a file stream, but this would involve writing new interfaces to Epetra's
  //      file I/O capabilities.
  std::string tname = "MATRIX";
  tname = tname + typeid(Scalar).name();
  tname = tname + typeid(LO).name();
  tname = tname + typeid(GO).name();
#ifdef HAVE_MUELU_KOKKOSCORE
  std::string nn = Kokkos::Compat::KokkosDeviceWrapperNode<typename Node::execution_space>::name();
  nn.erase(std::remove(nn.begin(), nn.end(), '/'), nn.end());
  tname = tname + nn;
#endif
  tname = "_" + tname;

  const bool alwaysWriteMaps = true;
  Xpetra::IO<Scalar, LO, GO, Node>::WriteBlockedCrsMatrix(tname, *bMat, alwaysWriteMaps);
  Teuchos::RCP<const Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> > bMat2 = Xpetra::IO<Scalar, LO, GO, Node>::ReadBlockedCrsMatrix(tname, lib, comm);

  TEST_EQUALITY(bMat->getMatrix(0,0)->getGlobalNumEntries(),bMat2->getMatrix(0,0)->getGlobalNumEntries());
  TEST_EQUALITY(bMat->getMatrix(0,1)->getGlobalNumEntries(),bMat2->getMatrix(0,1)->getGlobalNumEntries());
  TEST_EQUALITY(bMat->getMatrix(1,0)->getGlobalNumEntries(),bMat2->getMatrix(1,0)->getGlobalNumEntries());
  TEST_EQUALITY(bMat->getMatrix(1,1)->getGlobalNumEntries(),bMat2->getMatrix(1,1)->getGlobalNumEntries());

  TEST_EQUALITY(bMat->getMatrix(0,0)->getLocalNumEntries(),bMat2->getMatrix(0,0)->getLocalNumEntries());
  TEST_EQUALITY(bMat->getMatrix(0,1)->getLocalNumEntries(),bMat2->getMatrix(0,1)->getLocalNumEntries());
  TEST_EQUALITY(bMat->getMatrix(1,0)->getLocalNumEntries(),bMat2->getMatrix(1,0)->getLocalNumEntries());
  TEST_EQUALITY(bMat->getMatrix(1,1)->getLocalNumEntries(),bMat2->getMatrix(1,1)->getLocalNumEntries());

  TEST_EQUALITY(bMat->getMatrix(0,0)->getFrobeniusNorm(),bMat2->getMatrix(0,0)->getFrobeniusNorm());
  TEST_EQUALITY(bMat->getMatrix(0,1)->getFrobeniusNorm(),bMat2->getMatrix(0,1)->getFrobeniusNorm());
  TEST_EQUALITY(bMat->getMatrix(1,0)->getFrobeniusNorm(),bMat2->getMatrix(1,0)->getFrobeniusNorm());
  TEST_EQUALITY(bMat->getMatrix(1,1)->getFrobeniusNorm(),bMat2->getMatrix(1,1)->getFrobeniusNorm());

  TEST_EQUALITY(bMat->getRangeMapExtractor()->getMap(0)->isSameAs(*(bMat2->getRangeMapExtractor()->getMap(0))),true);
  TEST_EQUALITY(bMat->getDomainMapExtractor()->getMap(0)->isSameAs(*(bMat2->getDomainMapExtractor()->getMap(0))),true);

  TEST_EQUALITY(bMat->getRangeMapExtractor()->getFullMap()->isSameAs(*(bMat2->getRangeMapExtractor()->getFullMap())),true);
  TEST_EQUALITY(bMat->getDomainMapExtractor()->getFullMap()->isSameAs(*(bMat2->getDomainMapExtractor()->getFullMap())),true);

  // these tests are false with Tpetra? TODO check me: why only in Tpetra?
  // bMat2 is always in Xpetra mode so far. This is, since the Read routine and Write routine for the MapExtractor do not really
  // consider the Thyra mode so far.
  //TEST_EQUALITY(bMat->getRangeMapExtractor()->getMap(1)->isSameAs(*(bMat2->getRangeMapExtractor()->getMap(1))),true);
  //TEST_EQUALITY(bMat->getDomainMapExtractor()->getMap(1)->isSameAs(*(bMat2->getDomainMapExtractor()->getMap(1))),true);

  TEST_EQUALITY(bMat->getMatrix(0,0)->getRowMap()->isSameAs(*(bMat2->getMatrix(0,0)->getRowMap())),true);
  TEST_EQUALITY(bMat->getMatrix(0,1)->getRowMap()->isSameAs(*(bMat2->getMatrix(0,1)->getRowMap())),true);
  TEST_EQUALITY(bMat->getMatrix(1,0)->getRowMap()->isSameAs(*(bMat2->getMatrix(1,0)->getRowMap())),true);
  TEST_EQUALITY(bMat->getMatrix(1,1)->getRowMap()->isSameAs(*(bMat2->getMatrix(1,1)->getRowMap())),true);

  TEST_EQUALITY(bMat->getMatrix(0,0)->getColMap()->isSameAs(*(bMat2->getMatrix(0,0)->getColMap())),true);
  TEST_EQUALITY(bMat->getMatrix(0,1)->getColMap()->isSameAs(*(bMat2->getMatrix(0,1)->getColMap())),true);
  // the following test fails with Teptra. Why?
  //TEST_EQUALITY(bMat->getMatrix(1,0)->getColMap()->isSameAs(*(bMat2->getMatrix(1,0)->getColMap())),true);
  TEST_EQUALITY(bMat->getMatrix(1,1)->getColMap()->isSameAs(*(bMat2->getMatrix(1,1)->getColMap())),true);


  // build gloabl vector with one entries
  Teuchos::RCP<VectorClass> ones_A   = VectorFactoryClass::Build(bMat->getRangeMap(), true);
  Teuchos::RCP<VectorClass> exp  = VectorFactoryClass::Build(bMat->getRangeMap(), true);
  Teuchos::RCP<VectorClass> ones_bOp = VectorFactoryClass::Build(bMat2->getRangeMap(), true);
  Teuchos::RCP<VectorClass> res  = VectorFactoryClass::Build(bMat2->getRangeMap(), true);
  ones_A->putScalar(STS::one());
  ones_bOp->putScalar(STS::one());

  bMat->apply(*ones_A, *exp);
  bMat2->apply(*ones_bOp, *res);

  TEST_EQUALITY(res->norm2(),exp->norm2());
  TEST_EQUALITY(res->normInf(),exp->normInf());
  TEST_EQUALITY(bMat->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(bMat->getDomainMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(bMat->getRangeMap(0)->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(bMat->getRangeMap(1)->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(bMat->getDomainMap(0)->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(bMat->getDomainMap(1)->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(bMat2->getRangeMapExtractor()->getThyraMode(), false);  // thyra mode is not correctly transferred!!
  TEST_EQUALITY(bMat2->getDomainMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(bMat2->getRangeMap(0)->getMinAllGlobalIndex(), 0);
  TEST_INEQUALITY(bMat2->getRangeMap(1)->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(bMat2->getDomainMap(0)->getMinAllGlobalIndex(), 0);
  TEST_INEQUALITY(bMat2->getDomainMap(1)->getMinAllGlobalIndex(), 0);
#endif
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

// List of tests which run both with Epetra and Tpetra
#define XP_MATRIX_INSTANT(S,LO,GO,N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( ThyraBlockedOperator, ThyraOperator2XpetraCrsMat, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( ThyraBlockedOperator, ThyraShrinkMaps,            M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( ThyraBlockedOperator, SplitMatrixForThyra, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( ThyraBlockedOperator, ReadWriteMatrixMatrixMarket, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( ThyraBlockedOperator, XpetraBlockedCrsMatConstructor, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( ThyraBlockedOperator, ThyraBlockedOperator2XpetraBlockedCrsMat, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( ThyraBlockedOperator, NestedBlockOperator2ThyraBlockedCrsMat, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( ThyraBlockedOperator, ReorderBlockOperator2ThyraBlockedCrsMat, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N )

// List of tests which run only with Tpetra
#define XP_TPETRA_MATRIX_INSTANT(S,LO,GO,N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( ThyraBlockedOperator, ThyraVectorSpace2XpetraMap_Tpetra, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N )

// List of tests which run only with Epetra
#define XP_EPETRA_MATRIX_INSTANT(S,LO,GO,N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( ThyraBlockedOperator, ThyraVectorSpace2XpetraMap_Epetra, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N )


#include <TpetraCore_config.h>
#include <TpetraCore_ETIHelperMacros.h>

TPETRA_ETI_MANGLING_TYPEDEFS()
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( XPETRA_TPETRA_TYPES )
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( XP_MATRIX_INSTANT )
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( XP_TPETRA_MATRIX_INSTANT )


#if defined(HAVE_XPETRA_EPETRA)

#include "Xpetra_Map.hpp" // defines EpetraNode
typedef Xpetra::EpetraNode EpetraNode;
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
XPETRA_EPETRA_TYPES(double,int,int,EpetraNode)
XP_MATRIX_INSTANT(double,int,int,EpetraNode)
XP_EPETRA_MATRIX_INSTANT(double,int,int,EpetraNode)
#endif
// EpetraExt routines are not working with 64 bit
/*#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
typedef long long LongLong;
XPETRA_EPETRA_TYPES(double,int,LongLong,EpetraNode)
XP_EPETRA_MATRIX_INSTANT(double,int,LongLong,EpetraNode)
#endif*/

#endif

}
