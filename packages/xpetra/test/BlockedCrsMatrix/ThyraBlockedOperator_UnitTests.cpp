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
/*
 * BlockedCrsMatrix_UnitTests.cpp
 *
 *  Created on: Aug 22, 2011
 *      Author: wiesner
 */
#include <Teuchos_UnitTestHarness.hpp>
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
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>

#ifdef HAVE_XPETRA_THYRA
#include <Xpetra_ThyraUtils.hpp>
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


TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( ThyraBlockedOperator, ThyraVectorSpace2XpetraMap, Scalar, LO, GO, Node )
{
#ifdef HAVE_XPETRA_THYRA
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // TPetra version
#ifdef HAVE_XPETRA_TPETRA
#ifdef HAVE_XPETRA_TPETRA_INST_INT_INT // test only if LO=GO=int active
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
#endif

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

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( ThyraBlockedOperator, ThyraOperator2XpetraCrsMat, Scalar, LO, GO, Node )
{
#ifdef HAVE_XPETRA_THYRA
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  std::vector<Xpetra::UnderlyingLib> libs;

  // TPetra version
#ifdef HAVE_XPETRA_TPETRA
#ifdef HAVE_XPETRA_TPETRA_INST_INT_INT // only if LO=GO=int enabled
  libs.push_back(Xpetra::UseTpetra);
#endif
#endif

  // Epetra version
#ifdef HAVE_XPETRA_EPETRA
  libs.push_back(Xpetra::UseEpetra);
#endif

  typedef Xpetra::Map<LO, GO, Node> MapClass;
  typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;

  for(size_t ll = 0; ll<libs.size(); ++ll) {
    Xpetra::UnderlyingLib lib = libs[ll];

    // generate the matrix
    LO nEle = 63;
    const Teuchos::RCP<const MapClass> map = MapFactoryClass::Build(lib, nEle, 0, comm);

    Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > matrix =
                 Xpetra::CrsMatrixFactory<Scalar,LO,GO,Node>::Build(map, 10);

    LO NumMyElements = map->getNodeNumElements();
    Teuchos::ArrayView<const GO> MyGlobalElements = map->getNodeElementList();

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
    Teuchos::RCP<const Xpetra::CrsMatrix<Scalar, LO, GO, Node> > xCrsMat =
        Xpetra::ThyraUtils<Scalar,LO,GO,Node>::toXpetra(thyraOp);
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xCrsMat));

    TEST_EQUALITY(xCrsMat->getFrobeniusNorm()   ,matrix->getFrobeniusNorm());
    TEST_EQUALITY(xCrsMat->getGlobalNumRows()   ,matrix->getGlobalNumRows());
    TEST_EQUALITY(xCrsMat->getGlobalNumCols()   ,matrix->getGlobalNumCols());
    TEST_EQUALITY(xCrsMat->getNodeNumRows()     ,matrix->getNodeNumRows()  );
    TEST_EQUALITY(xCrsMat->getGlobalNumEntries(),matrix->getGlobalNumEntries());
    TEST_EQUALITY(xCrsMat->getNodeNumEntries()  ,matrix->getNodeNumEntries());
  } //end for
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( ThyraBlockedOperator, ThyraBlockedOperator2XpetraBlockedCrsMat, Scalar, LO, GO, Node )
{
#ifdef HAVE_XPETRA_THYRA
#ifdef HAVE_XPETRA_EPETRAEXT

  Teuchos::RCP<Epetra_Comm> Comm;
#ifdef HAVE_MPI
  Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  Comm = Teuchos::rcp(new Epetra_SerialComm);
#endif

  // 1) load all matrices
  Epetra_Map pointmap(36,0,*Comm);

  // generate local maps for loading matrices
  std::vector<int> velgidvec; // global strided maps
  std::vector<int> pregidvec;
  std::vector<int> fullgidvec; // full global map
  for (int i=0; i<pointmap.NumMyElements(); i++)
  {
    // loop over all local ids in pointmap

    // get corresponding global id
    int gid = pointmap.GID(i);

    // store global strided gids
    velgidvec.push_back(3*gid);
    velgidvec.push_back(3*gid+1);
    pregidvec.push_back(3*gid+2);

    // gid for full map
    fullgidvec.push_back(3*gid);
    fullgidvec.push_back(3*gid+1);
    fullgidvec.push_back(3*gid+2);
  }

  // generate strided maps
  Teuchos::RCP<const Epetra_Map> velmap = Teuchos::rcp(new const Epetra_Map(-1, velgidvec.size(), &velgidvec[0], 0, *Comm));
  Teuchos::RCP<const Epetra_Map> premap = Teuchos::rcp(new const Epetra_Map(-1, pregidvec.size(), &pregidvec[0], 0, *Comm));

  // generate full map
  const Teuchos::RCP<const Epetra_Map> fullmap = Teuchos::rcp(new const Epetra_Map(-1, fullgidvec.size(), &fullgidvec[0], 0, *Comm));

  // read in matrices
  Epetra_CrsMatrix* ptrA = 0;
  EpetraExt::MatrixMarketFileToCrsMatrix("A.mat",*fullmap,*fullmap,*fullmap,ptrA);

  Teuchos::RCP<Epetra_CrsMatrix> fullA = Teuchos::rcp(ptrA);
  Teuchos::RCP<Epetra_Vector>    x     = Teuchos::rcp(new Epetra_Vector(*fullmap));
  x->PutScalar(1.0);

  // split fullA into A11,..., A22
  Teuchos::RCP<Epetra_CrsMatrix> A11;
  Teuchos::RCP<Epetra_CrsMatrix> A12;
  Teuchos::RCP<Epetra_CrsMatrix> A21;
  Teuchos::RCP<Epetra_CrsMatrix> A22;

  TEST_EQUALITY(SplitMatrix2x2(fullA,*velmap,*premap,A11,A12,A21,A22),true);

  // build Xpetra objects from Epetra_CrsMatrix objects
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xfuA = Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<GO,Node>(fullA));
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xA11 = Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<GO,Node>(A11));
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xA12 = Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<GO,Node>(A12));
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xA21 = Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<GO,Node>(A21));
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xA22 = Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<GO,Node>(A22));

  // build map extractor
  Teuchos::RCP<Xpetra::EpetraMapT<GO,Node> > xfullmap = Teuchos::rcp(new Xpetra::EpetraMapT<GO,Node>(fullmap));
  Teuchos::RCP<Xpetra::EpetraMapT<GO,Node> > xvelmap  = Teuchos::rcp(new Xpetra::EpetraMapT<GO,Node>(velmap ));
  Teuchos::RCP<Xpetra::EpetraMapT<GO,Node> > xpremap  = Teuchos::rcp(new Xpetra::EpetraMapT<GO,Node>(premap ));

  std::vector<Teuchos::RCP<const Xpetra::Map<LO,GO,Node> > > xmaps;
  xmaps.push_back(xvelmap);
  xmaps.push_back(xpremap);

  Teuchos::RCP<const Xpetra::MapExtractor<Scalar,LO,GO,Node> > map_extractor = Xpetra::MapExtractorFactory<Scalar,LO,GO,Node>::Build(xfullmap,xmaps);

  // build blocked operator
  Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> > bOp = Teuchos::rcp(new Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node>(map_extractor,map_extractor,10));
  bOp->setMatrix(0,0,xA11);
  bOp->setMatrix(0,1,xA12);
  bOp->setMatrix(1,0,xA21);
  bOp->setMatrix(1,1,xA22);

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
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productRange->dim()) ,xfullmap->getGlobalNumElements());
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productDomain->dim()) ,xfullmap->getGlobalNumElements());
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productRange->getBlock(0)->dim())  ,xvelmap->getGlobalNumElements());
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productDomain->getBlock(0)->dim()) ,xvelmap->getGlobalNumElements());
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productRange->getBlock(1)->dim())  ,xpremap->getGlobalNumElements());
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productDomain->getBlock(1)->dim()) ,xpremap->getGlobalNumElements());

#endif // end HAVE_XPETRA_EPETRAEXT
#endif // end HAVE_XPETRA_THYRA
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( ThyraBlockedOperator, XpetraBlockedCrsMatConstructor, Scalar, LO, GO, Node )
{
#ifdef HAVE_XPETRA_THYRA
#ifdef HAVE_XPETRA_EPETRAEXT
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  const Teuchos::RCP<const Epetra_Comm> epComm = Xpetra::toEpetra(comm);

  // 1) load all matrices
  Epetra_Map pointmap(36,0,*epComm);

  // generate local maps for loading matrices
  std::vector<int> velgidvec; // global strided maps
  std::vector<int> pregidvec;
  std::vector<int> fullgidvec; // full global map
  for (int i=0; i<pointmap.NumMyElements(); i++)
  {
    // loop over all local ids in pointmap

    // get corresponding global id
    int gid = pointmap.GID(i);

    // store global strided gids
    velgidvec.push_back(3*gid);
    velgidvec.push_back(3*gid+1);
    pregidvec.push_back(3*gid+2);

    // gid for full map
    fullgidvec.push_back(3*gid);
    fullgidvec.push_back(3*gid+1);
    fullgidvec.push_back(3*gid+2);
  }

  // generate strided maps
  Teuchos::RCP<const Epetra_Map> velmap = Teuchos::rcp(new const Epetra_Map(-1, velgidvec.size(), &velgidvec[0], 0, *epComm));
  Teuchos::RCP<const Epetra_Map> premap = Teuchos::rcp(new const Epetra_Map(-1, pregidvec.size(), &pregidvec[0], 0, *epComm));

  // generate full map
  const Teuchos::RCP<const Epetra_Map> fullmap = Teuchos::rcp(new const Epetra_Map(-1, fullgidvec.size(), &fullgidvec[0], 0, *epComm));

  // read in matrices
  Epetra_CrsMatrix* ptrA = 0;
  EpetraExt::MatrixMarketFileToCrsMatrix("A.mat",*fullmap,*fullmap,*fullmap,ptrA);

  Teuchos::RCP<Epetra_CrsMatrix> fullA = Teuchos::rcp(ptrA);
  Teuchos::RCP<Epetra_Vector>    x     = Teuchos::rcp(new Epetra_Vector(*fullmap));
  x->PutScalar(1.0);

  // split fullA into A11,..., A22
  Teuchos::RCP<Epetra_CrsMatrix> A11;
  Teuchos::RCP<Epetra_CrsMatrix> A12;
  Teuchos::RCP<Epetra_CrsMatrix> A21;
  Teuchos::RCP<Epetra_CrsMatrix> A22;

  TEST_EQUALITY(SplitMatrix2x2(fullA,*velmap,*premap,A11,A12,A21,A22),true);

  // build Xpetra objects from Epetra_CrsMatrix objects
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xfuA = Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<GO,Node>(fullA));
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xA11 = Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<GO,Node>(A11));
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xA12 = Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<GO,Node>(A12));
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xA21 = Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<GO,Node>(A21));
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xA22 = Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<GO,Node>(A22));

  // build map extractor
  Teuchos::RCP<Xpetra::EpetraMapT<GO,Node> > xfullmap = Teuchos::rcp(new Xpetra::EpetraMapT<GO,Node>(fullmap));
  Teuchos::RCP<Xpetra::EpetraMapT<GO,Node> > xvelmap  = Teuchos::rcp(new Xpetra::EpetraMapT<GO,Node>(velmap ));
  Teuchos::RCP<Xpetra::EpetraMapT<GO,Node> > xpremap  = Teuchos::rcp(new Xpetra::EpetraMapT<GO,Node>(premap ));

  std::vector<Teuchos::RCP<const Xpetra::Map<LO,GO,Node> > > xmaps;
  xmaps.push_back(xvelmap);
  xmaps.push_back(xpremap);

  Teuchos::RCP<const Xpetra::MapExtractor<Scalar,LO,GO,Node> > map_extractor = Xpetra::MapExtractorFactory<Scalar,LO,GO,Node>::Build(xfullmap,xmaps);

  // build blocked operator
  Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> > bOp = Teuchos::rcp(new Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node>(map_extractor,map_extractor,10));
  bOp->setMatrix(0,0,xA11);
  bOp->setMatrix(0,1,xA12);
  bOp->setMatrix(1,0,xA21);
  bOp->setMatrix(1,1,xA22);

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
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productRange->dim())  ,xfullmap->getGlobalNumElements());
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productDomain->dim()) ,xfullmap->getGlobalNumElements());
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productRange->getBlock(0)->dim())  ,xvelmap->getGlobalNumElements());
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productDomain->getBlock(0)->dim()) ,xvelmap->getGlobalNumElements());
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productRange->getBlock(1)->dim())  ,xpremap->getGlobalNumElements());
  TEST_EQUALITY(Teuchos::as<Xpetra::global_size_t>(productDomain->getBlock(1)->dim()) ,xpremap->getGlobalNumElements());

  //construct a Xpetra::BlockedCrsMatrix object
  Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> > bOp2 =
      Teuchos::rcp(new Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node>(thbOp,comm));

  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(bOp2));
  TEST_EQUALITY(bOp2->getGlobalNumRows(),bOp->getGlobalNumRows());
  TEST_EQUALITY(bOp2->getGlobalNumCols(),bOp->getGlobalNumCols());
  TEST_EQUALITY(bOp2->getNodeNumRows(),bOp->getNodeNumRows());
  TEST_EQUALITY(bOp2->getGlobalNumEntries(),bOp->getGlobalNumEntries());
  TEST_EQUALITY(bOp2->getNodeNumEntries(),bOp->getNodeNumEntries());
  TEST_EQUALITY(bOp2->isFillComplete(),bOp->isFillComplete());
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

#endif // end HAVE_XPETRA_EPETRAEXT
#endif // end HAVE_XPETRA_THYRA
}

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP_ORDINAL( SC, LO, GO, Node )                     \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( ThyraBlockedOperator, ThyraVectorSpace2XpetraMap, SC, LO, GO, Node ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( ThyraBlockedOperator, ThyraOperator2XpetraCrsMat, SC, LO, GO, Node ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( ThyraBlockedOperator, ThyraBlockedOperator2XpetraBlockedCrsMat, SC, LO, GO, Node ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( ThyraBlockedOperator, XpetraBlockedCrsMatConstructor, SC, LO, GO, Node)

// TODO: fix me
#ifdef HAVE_XPETRA_SERIAL
typedef Kokkos::Compat::KokkosSerialWrapperNode DefaultNodeType;
// all tests work only for LO=GO=int (and/or Epetra)
// So, it is sufficient to define only a TEST_GROUP for LO=GO=int
UNIT_TEST_GROUP_ORDINAL(double, int, int, DefaultNodeType)
#endif


}
