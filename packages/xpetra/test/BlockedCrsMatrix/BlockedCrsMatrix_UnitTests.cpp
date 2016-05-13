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
#include <Xpetra_UnitTestHelpers.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>

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

// Epetra routines to split matrix and maps
#include "BlockedMatrixTestHelpers.hpp"

#include <Xpetra_DefaultPlatform.hpp>
#include <Teuchos_as.hpp>

#include <Xpetra_Map.hpp>

#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_Exceptions.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_MatrixUtils.hpp>
#include <Xpetra_IO.hpp>
#include <Xpetra_BlockReorderManager.hpp>
#include <Xpetra_ReorderedBlockedCrsMatrix.hpp>

namespace XpetraBlockMatrixTests {

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

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedCrsMatrix, SplitMatrix, M, MA, Scalar, LO, GO, Node )
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

  LO NumMyElements = map->getNodeNumElements();
  GO NumGlobalElements = map->getGlobalNumElements();
  Teuchos::ArrayView<const GO> MyGlobalElements = map->getNodeElementList();

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
      Xpetra::MatrixUtils<Scalar, LO, GO, Node>::SplitMatrix(*mat,map_extractor,map_extractor);

  // build gloabl vector with one entries
  Teuchos::RCP<VectorClass> ones = VectorFactoryClass::Build(map, true);
  Teuchos::RCP<VectorClass> exp  = VectorFactoryClass::Build(map, true);
  Teuchos::RCP<VectorClass> res  = VectorFactoryClass::Build(map, true);
  Teuchos::RCP<VectorClass> rnd  = VectorFactoryClass::Build(map, true);
  ones->putScalar(STS::one());
  rnd->randomize();

  A->apply(*ones, *exp);
  bOp->apply(*ones, *res);
  res->update(-STS::one(),*exp,STS::one());
  TEUCHOS_TEST_COMPARE(res->norm2(), <, 1e-16, out, success);
  TEUCHOS_TEST_COMPARE(res->normInf(), <, 1e-16, out, success);

  A->apply(*rnd, *exp);
  bOp->apply(*rnd, *res);
  res->update(-STS::one(),*exp,STS::one());
  TEUCHOS_TEST_COMPARE(res->norm2(), <, 5e-14, out, success);
  TEUCHOS_TEST_COMPARE(res->normInf(), <, 5e-14, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedCrsMatrix, CreateBlockedDiagonalOp, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> BlockedCrsMatrixClass;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 4;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar,LO,GO,Node,M>(noBlocks, *comm);

  GO goNumRows = Teuchos::as<GO>(Teuchos::ScalarTraits<GO>::pow(2,noBlocks-2)) * 10 * comm->getSize();

  TEST_EQUALITY(bop->Rows(),4);
  TEST_EQUALITY(bop->Cols(),4);
  TEST_EQUALITY(bop->getRangeMap()->getGlobalNumElements(),goNumRows);
  TEST_EQUALITY(bop->getDomainMap()->getGlobalNumElements(),goNumRows);
  TEST_EQUALITY(bop->getMatrix(0,0)->getRowMap()->getMinGlobalIndex(),comm->getRank() * 40);
  TEST_EQUALITY(bop->getMatrix(0,0)->getRowMap()->getMaxGlobalIndex(),comm->getRank() * 40 + 4);
  TEST_EQUALITY(bop->getMatrix(1,1)->getRowMap()->getMinGlobalIndex(),comm->getRank() * 40 + 5);
  TEST_EQUALITY(bop->getMatrix(1,1)->getRowMap()->getMaxGlobalIndex(),comm->getRank() * 40 + 9);
  TEST_EQUALITY(bop->getMatrix(2,2)->getRowMap()->getMinGlobalIndex(),comm->getRank() * 40 + 10);
  TEST_EQUALITY(bop->getMatrix(2,2)->getRowMap()->getMaxGlobalIndex(),comm->getRank() * 40 + 19);
  TEST_EQUALITY(bop->getMatrix(3,3)->getRowMap()->getMinGlobalIndex(),comm->getRank() * 40 + 20);
  TEST_EQUALITY(bop->getMatrix(3,3)->getRowMap()->getMaxGlobalIndex(),comm->getRank() * 40 + 39);
  TEST_EQUALITY(bop->getMatrix(0,0)->getRowMap()->getGlobalNumElements(),comm->getSize() * 5);
  TEST_EQUALITY(bop->getMatrix(1,1)->getRowMap()->getGlobalNumElements(),comm->getSize() * 5);
  TEST_EQUALITY(bop->getMatrix(2,2)->getRowMap()->getGlobalNumElements(),comm->getSize() * 10);
  TEST_EQUALITY(bop->getMatrix(3,3)->getRowMap()->getGlobalNumElements(),comm->getSize() * 20);
  TEST_EQUALITY(bop->getMatrix(0,0)->getRowMap()->getNodeNumElements(),5);
  TEST_EQUALITY(bop->getMatrix(1,1)->getRowMap()->getNodeNumElements(),5);
  TEST_EQUALITY(bop->getMatrix(2,2)->getRowMap()->getNodeNumElements(),10);
  TEST_EQUALITY(bop->getMatrix(3,3)->getRowMap()->getNodeNumElements(),20);
  TEST_EQUALITY(bop->getMatrix(0,0)->getColMap()->getMinGlobalIndex(),comm->getRank() * 40);
  TEST_EQUALITY(bop->getMatrix(0,0)->getColMap()->getMaxGlobalIndex(),comm->getRank() * 40 + 4);
  TEST_EQUALITY(bop->getMatrix(1,1)->getColMap()->getMinGlobalIndex(),comm->getRank() * 40 + 5);
  TEST_EQUALITY(bop->getMatrix(1,1)->getColMap()->getMaxGlobalIndex(),comm->getRank() * 40 + 9);
  TEST_EQUALITY(bop->getMatrix(2,2)->getColMap()->getMinGlobalIndex(),comm->getRank() * 40 + 10);
  TEST_EQUALITY(bop->getMatrix(2,2)->getColMap()->getMaxGlobalIndex(),comm->getRank() * 40 + 19);
  TEST_EQUALITY(bop->getMatrix(3,3)->getColMap()->getMinGlobalIndex(),comm->getRank() * 40 + 20);
  TEST_EQUALITY(bop->getMatrix(3,3)->getColMap()->getMaxGlobalIndex(),comm->getRank() * 40 + 39);
  TEST_EQUALITY(bop->getMatrix(0,0)->getColMap()->getGlobalNumElements(),comm->getSize() * 5);
  TEST_EQUALITY(bop->getMatrix(1,1)->getColMap()->getGlobalNumElements(),comm->getSize() * 5);
  TEST_EQUALITY(bop->getMatrix(2,2)->getColMap()->getGlobalNumElements(),comm->getSize() * 10);
  TEST_EQUALITY(bop->getMatrix(3,3)->getColMap()->getGlobalNumElements(),comm->getSize() * 20);
  TEST_EQUALITY(bop->getMatrix(0,0)->getColMap()->getNodeNumElements(),5);
  TEST_EQUALITY(bop->getMatrix(1,1)->getColMap()->getNodeNumElements(),5);
  TEST_EQUALITY(bop->getMatrix(2,2)->getColMap()->getNodeNumElements(),10);
  TEST_EQUALITY(bop->getMatrix(3,3)->getColMap()->getNodeNumElements(),20);

  TEST_EQUALITY(bop->getMatrix(0,1)->getRowMap()->getMinGlobalIndex(),comm->getRank() * 40);
  TEST_EQUALITY(bop->getMatrix(0,1)->getRowMap()->getMaxGlobalIndex(),comm->getRank() * 40 + 4);
  TEST_EQUALITY(bop->getMatrix(0,2)->getRowMap()->getMinGlobalIndex(),comm->getRank() * 40);
  TEST_EQUALITY(bop->getMatrix(0,2)->getRowMap()->getMaxGlobalIndex(),comm->getRank() * 40 + 4);
  TEST_EQUALITY(bop->getMatrix(0,3)->getRowMap()->getMinGlobalIndex(),comm->getRank() * 40);
  TEST_EQUALITY(bop->getMatrix(0,3)->getRowMap()->getMaxGlobalIndex(),comm->getRank() * 40 + 4);
  TEST_EQUALITY(bop->getMatrix(1,0)->getColMap()->getMinGlobalIndex(),0); // TODO
  TEST_EQUALITY(bop->getMatrix(1,0)->getDomainMap()->getMaxGlobalIndex(),comm->getRank() * 40 + 4);
  //TEST_EQUALITY(bop->getMatrix(1,0)->getColMap()->getMaxGlobalIndex(),comm->getRank() * 40);
  TEST_EQUALITY(bop->getMatrix(2,0)->getColMap()->getMinGlobalIndex(),0); // TODO
  TEST_EQUALITY(bop->getMatrix(2,0)->getDomainMap()->getMaxGlobalIndex(),comm->getRank() * 40 + 4);
  //TEST_EQUALITY(bop->getMatrix(2,0)->getColMap()->getMaxGlobalIndex(),comm->getRank() * 40);
  TEST_EQUALITY(bop->getMatrix(3,0)->getColMap()->getMinGlobalIndex(),0); // TODO
  TEST_EQUALITY(bop->getMatrix(3,0)->getDomainMap()->getMaxGlobalIndex(),comm->getRank() * 40 + 4);
  //TEST_EQUALITY(bop->getMatrix(3,0)->getColMap()->getMaxGlobalIndex(),comm->getRank() * 40);

  TEST_EQUALITY(bop->getMatrix(2,1)->getRowMap()->getMinGlobalIndex(),comm->getRank() * 40 + 10);
  TEST_EQUALITY(bop->getMatrix(2,1)->getRowMap()->getMaxGlobalIndex(),comm->getRank() * 40 + 19);
  TEST_EQUALITY(bop->getMatrix(2,2)->getRowMap()->getMinGlobalIndex(),comm->getRank() * 40 + 10);
  TEST_EQUALITY(bop->getMatrix(2,2)->getRowMap()->getMaxGlobalIndex(),comm->getRank() * 40 + 19);
  TEST_EQUALITY(bop->getMatrix(2,3)->getRowMap()->getMinGlobalIndex(),comm->getRank() * 40 + 10);
  TEST_EQUALITY(bop->getMatrix(2,3)->getRowMap()->getMaxGlobalIndex(),comm->getRank() * 40 + 19);

  TEST_EQUALITY(bop->getMatrix(0,0)->isFillComplete(),true);
  TEST_EQUALITY(bop->getMatrix(2,2)->isFillComplete(),true);
  TEST_EQUALITY(bop->getMatrix(2,3)->isFillComplete(),true);
  TEST_EQUALITY(bop->getMatrix(1,0)->isFillComplete(),true);
  TEST_EQUALITY(bop->getMatrix(3,1)->isFillComplete(),true);

  TEST_EQUALITY(bop->getRangeMapExtractor()->getThyraMode(),false);
  TEST_EQUALITY(bop->getDomainMapExtractor()->getThyraMode(),false);

  TEST_THROW(bop->getRangeMap(0,true), Xpetra::Exceptions::RuntimeError);

  TEST_EQUALITY(bop->getRangeMapExtractor()->getMap(2)->getMinAllGlobalIndex(),10);
  TEST_EQUALITY(bop->getRangeMapExtractor()->getMap(2)->getMaxAllGlobalIndex(),(comm->getSize() - 1) * 40 + 19);

  TEST_EQUALITY(bop->getRangeMapExtractor()->getMap(3,false)->getMinGlobalIndex(),comm->getRank() * 40 + 20);
  TEST_EQUALITY(bop->getRangeMapExtractor()->getMap(3,false)->getMaxGlobalIndex(),comm->getRank() * 40 + 39);
  TEST_EQUALITY(bop->getRangeMapExtractor()->getMap(3,false)->getMinAllGlobalIndex(),20);
  TEST_EQUALITY(bop->getRangeMapExtractor()->getMap(3,false)->getMaxAllGlobalIndex(),comm->getSize() * 40 - 1);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedCrsMatrix, CreateBlockedDiagonalOpThyra, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> BlockedCrsMatrixClass;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 4;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrixThyra<Scalar,LO,GO,Node,M>(noBlocks, *comm);

  GO goNumRows = Teuchos::as<GO>(Teuchos::ScalarTraits<GO>::pow(2,noBlocks-2)) * 10 * comm->getSize();

  TEST_EQUALITY(bop->Rows(),4);
  TEST_EQUALITY(bop->Cols(),4);
  TEST_EQUALITY(bop->getRangeMap()->getGlobalNumElements(),goNumRows);
  TEST_EQUALITY(bop->getDomainMap()->getGlobalNumElements(),goNumRows);
  // Thyra GIDs
  TEST_EQUALITY(bop->getMatrix(0,0)->getRowMap()->getMinGlobalIndex(),comm->getRank() * 5);
  TEST_EQUALITY(bop->getMatrix(0,0)->getRowMap()->getMaxGlobalIndex(),comm->getRank() * 5 + 4);
  TEST_EQUALITY(bop->getMatrix(1,1)->getRowMap()->getMinGlobalIndex(),comm->getRank() * 5);
  TEST_EQUALITY(bop->getMatrix(1,1)->getRowMap()->getMaxGlobalIndex(),comm->getRank() * 5 + 4);
  TEST_EQUALITY(bop->getMatrix(2,2)->getRowMap()->getMinGlobalIndex(),comm->getRank() * 10);
  TEST_EQUALITY(bop->getMatrix(2,2)->getRowMap()->getMaxGlobalIndex(),comm->getRank() * 10 + 9);
  TEST_EQUALITY(bop->getMatrix(3,3)->getRowMap()->getMinGlobalIndex(),comm->getRank() * 20);
  TEST_EQUALITY(bop->getMatrix(3,3)->getRowMap()->getMaxGlobalIndex(),comm->getRank() * 20 + 19);
  TEST_EQUALITY(bop->getMatrix(0,0)->getRowMap()->getGlobalNumElements(),comm->getSize() * 5);
  TEST_EQUALITY(bop->getMatrix(1,1)->getRowMap()->getGlobalNumElements(),comm->getSize() * 5);
  TEST_EQUALITY(bop->getMatrix(2,2)->getRowMap()->getGlobalNumElements(),comm->getSize() * 10);
  TEST_EQUALITY(bop->getMatrix(3,3)->getRowMap()->getGlobalNumElements(),comm->getSize() * 20);
  TEST_EQUALITY(bop->getMatrix(0,0)->getRowMap()->getNodeNumElements(),5);
  TEST_EQUALITY(bop->getMatrix(1,1)->getRowMap()->getNodeNumElements(),5);
  TEST_EQUALITY(bop->getMatrix(2,2)->getRowMap()->getNodeNumElements(),10);
  TEST_EQUALITY(bop->getMatrix(3,3)->getRowMap()->getNodeNumElements(),20);
  TEST_EQUALITY(bop->getMatrix(0,0)->getColMap()->getMinGlobalIndex(),comm->getRank() * 5);
  TEST_EQUALITY(bop->getMatrix(0,0)->getColMap()->getMaxGlobalIndex(),comm->getRank() * 5 + 4);
  TEST_EQUALITY(bop->getMatrix(1,1)->getColMap()->getMinGlobalIndex(),comm->getRank() * 5);
  TEST_EQUALITY(bop->getMatrix(1,1)->getColMap()->getMaxGlobalIndex(),comm->getRank() * 5 + 4);
  TEST_EQUALITY(bop->getMatrix(2,2)->getColMap()->getMinGlobalIndex(),comm->getRank() * 10);
  TEST_EQUALITY(bop->getMatrix(2,2)->getColMap()->getMaxGlobalIndex(),comm->getRank() * 10 + 9);
  TEST_EQUALITY(bop->getMatrix(3,3)->getColMap()->getMinGlobalIndex(),comm->getRank() * 20);
  TEST_EQUALITY(bop->getMatrix(3,3)->getColMap()->getMaxGlobalIndex(),comm->getRank() * 20 + 19);
  TEST_EQUALITY(bop->getMatrix(0,0)->getColMap()->getGlobalNumElements(),comm->getSize() * 5);
  TEST_EQUALITY(bop->getMatrix(1,1)->getColMap()->getGlobalNumElements(),comm->getSize() * 5);
  TEST_EQUALITY(bop->getMatrix(2,2)->getColMap()->getGlobalNumElements(),comm->getSize() * 10);
  TEST_EQUALITY(bop->getMatrix(3,3)->getColMap()->getGlobalNumElements(),comm->getSize() * 20);
  TEST_EQUALITY(bop->getMatrix(0,0)->getColMap()->getNodeNumElements(),5);
  TEST_EQUALITY(bop->getMatrix(1,1)->getColMap()->getNodeNumElements(),5);
  TEST_EQUALITY(bop->getMatrix(2,2)->getColMap()->getNodeNumElements(),10);
  TEST_EQUALITY(bop->getMatrix(3,3)->getColMap()->getNodeNumElements(),20);

  TEST_EQUALITY(bop->getMatrix(0,1)->getRowMap()->getMinGlobalIndex(),comm->getRank() * 5);
  TEST_EQUALITY(bop->getMatrix(0,1)->getRowMap()->getMaxGlobalIndex(),comm->getRank() * 5 + 4);
  TEST_EQUALITY(bop->getMatrix(0,2)->getRowMap()->getMinGlobalIndex(),comm->getRank() * 5);
  TEST_EQUALITY(bop->getMatrix(0,2)->getRowMap()->getMaxGlobalIndex(),comm->getRank() * 5 + 4);
  TEST_EQUALITY(bop->getMatrix(0,3)->getRowMap()->getMinGlobalIndex(),comm->getRank() * 5);
  TEST_EQUALITY(bop->getMatrix(0,3)->getRowMap()->getMaxGlobalIndex(),comm->getRank() * 5 + 4);
  TEST_EQUALITY(bop->getMatrix(1,0)->getColMap()->getMinGlobalIndex(),0); // TODO
  TEST_EQUALITY(bop->getMatrix(1,0)->getDomainMap()->getMaxGlobalIndex(),comm->getRank() * 5 + 4);
  //TEST_EQUALITY(bop->getMatrix(1,0)->getColMap()->getMaxGlobalIndex(),comm->getRank() * 5);
  TEST_EQUALITY(bop->getMatrix(2,0)->getColMap()->getMinGlobalIndex(),0); // TODO
  TEST_EQUALITY(bop->getMatrix(2,0)->getDomainMap()->getMaxGlobalIndex(),comm->getRank() * 5 + 4);
  //TEST_EQUALITY(bop->getMatrix(2,0)->getColMap()->getMaxGlobalIndex(),comm->getRank() * 5);
  TEST_EQUALITY(bop->getMatrix(3,0)->getColMap()->getMinGlobalIndex(),0); // TODO
  TEST_EQUALITY(bop->getMatrix(3,0)->getDomainMap()->getMaxGlobalIndex(),comm->getRank() * 5 + 4);
  //TEST_EQUALITY(bop->getMatrix(3,0)->getColMap()->getMaxGlobalIndex(),comm->getRank() * 5);

  TEST_EQUALITY(bop->getMatrix(2,1)->getRowMap()->getMinGlobalIndex(),comm->getRank() * 10);
  TEST_EQUALITY(bop->getMatrix(2,1)->getRowMap()->getMaxGlobalIndex(),comm->getRank() * 10 + 9);
  TEST_EQUALITY(bop->getMatrix(2,2)->getRowMap()->getMinGlobalIndex(),comm->getRank() * 10);
  TEST_EQUALITY(bop->getMatrix(2,2)->getRowMap()->getMaxGlobalIndex(),comm->getRank() * 10 + 9);
  TEST_EQUALITY(bop->getMatrix(2,3)->getRowMap()->getMinGlobalIndex(),comm->getRank() * 10);
  TEST_EQUALITY(bop->getMatrix(2,3)->getRowMap()->getMaxGlobalIndex(),comm->getRank() * 10 + 9);

  TEST_EQUALITY(bop->getMatrix(0,0)->isFillComplete(),true);
  TEST_EQUALITY(bop->getMatrix(2,2)->isFillComplete(),true);
  TEST_EQUALITY(bop->getMatrix(2,3)->isFillComplete(),true);
  TEST_EQUALITY(bop->getMatrix(1,0)->isFillComplete(),true);
  TEST_EQUALITY(bop->getMatrix(3,1)->isFillComplete(),true);

  // check Xpetra replacement maps
  TEST_EQUALITY(bop->getRangeMap(0,false)->getMinGlobalIndex(),comm->getRank() * 5 + 0);
  TEST_EQUALITY(bop->getRangeMap(0,false)->getMaxGlobalIndex(),comm->getRank() * 5 + 4);
  TEST_EQUALITY(bop->getRangeMap(1,false)->getMinGlobalIndex(),comm->getSize() * 5 + comm->getRank() * 5);
  TEST_EQUALITY(bop->getRangeMap(1,false)->getMaxGlobalIndex(),comm->getSize() * 5 + comm->getRank() * 5 + 4);
  TEST_EQUALITY(bop->getRangeMap(2,false)->getMinGlobalIndex(),comm->getSize() * 10 + comm->getRank() * 10);
  TEST_EQUALITY(bop->getRangeMap(2,false)->getMaxGlobalIndex(),comm->getSize() * 10 + comm->getRank() * 10 + 9);
  TEST_EQUALITY(bop->getRangeMap(3,false)->getMinGlobalIndex(),comm->getSize() * 20 + comm->getRank() * 20);
  TEST_EQUALITY(bop->getRangeMap(3,false)->getMaxGlobalIndex(),comm->getSize() * 20 + comm->getRank() * 20 + 19);

  // check Thyra maps
  TEST_EQUALITY(bop->getRangeMap(0)->getMinGlobalIndex(),comm->getRank() * 5 + 0);
  TEST_EQUALITY(bop->getRangeMap(0)->getMaxGlobalIndex(),comm->getRank() * 5 + 4);
  TEST_EQUALITY(bop->getRangeMap(1)->getMinGlobalIndex(),comm->getRank() * 5 + 0);
  TEST_EQUALITY(bop->getRangeMap(1)->getMaxGlobalIndex(),comm->getRank() * 5 + 4);
  TEST_EQUALITY(bop->getRangeMap(2)->getMinGlobalIndex(),comm->getRank() * 10 + 0);
  TEST_EQUALITY(bop->getRangeMap(2)->getMaxGlobalIndex(),comm->getRank() * 10 + 9);
  TEST_EQUALITY(bop->getRangeMap(3)->getMinGlobalIndex(),comm->getRank() * 20 + 0);
  TEST_EQUALITY(bop->getRangeMap(3)->getMaxGlobalIndex(),comm->getRank() * 20 + 19);

  TEST_EQUALITY(bop->getRangeMap(0)->getMinAllGlobalIndex(),0);
  TEST_EQUALITY(bop->getRangeMap(0)->getMaxAllGlobalIndex(),comm->getSize() * 5 - 1);
  TEST_EQUALITY(bop->getRangeMap(1)->getMinAllGlobalIndex(),0);
  TEST_EQUALITY(bop->getRangeMap(1)->getMaxAllGlobalIndex(),comm->getSize() * 5 - 1);
  TEST_EQUALITY(bop->getRangeMap(2)->getMinAllGlobalIndex(),0);
  TEST_EQUALITY(bop->getRangeMap(2)->getMaxAllGlobalIndex(),comm->getSize() * 10 - 1);
  TEST_EQUALITY(bop->getRangeMap(3)->getMinAllGlobalIndex(),0);
  TEST_EQUALITY(bop->getRangeMap(3)->getMaxAllGlobalIndex(),comm->getSize() * 20 - 1);

  TEST_EQUALITY(bop->getRangeMapExtractor()->getThyraMode(),true);
  TEST_EQUALITY(bop->getDomainMapExtractor()->getThyraMode(),true);

  // check Xpetra replacement submaps
  TEST_EQUALITY(bop->getRangeMapExtractor()->getMap(3,false)->getMinGlobalIndex(),comm->getSize() * 20 + comm->getRank() * 20);
  TEST_EQUALITY(bop->getRangeMapExtractor()->getMap(3,false)->getMaxGlobalIndex(),comm->getSize() * 20 + comm->getRank() * 20 + 19);
  TEST_EQUALITY(bop->getRangeMapExtractor()->getMap(3,false)->getMinAllGlobalIndex(),comm->getSize() * 20);
  TEST_EQUALITY(bop->getRangeMapExtractor()->getMap(3,false)->getMaxAllGlobalIndex(),comm->getSize() * 40 - 1);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedCrsMatrix, ReorderBlockOperator, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> BlockedCrsMatrixClass;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 8;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar,LO,GO,Node,M>(noBlocks, *comm);


  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ [ 0 [ [1 2] 3] ] 4 [ 5 6 7] ]");

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> > brop =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> >(buildReorderedBlockedCrsMatrix(brm, bop));

  GO goNumRows = Teuchos::as<GO>(Teuchos::ScalarTraits<GO>::pow(2,noBlocks-2)) * 10 * comm->getSize();

  TEST_EQUALITY(brop->Rows(),3);
  TEST_EQUALITY(brop->Cols(),3);
  TEST_EQUALITY(brop->getRangeMap()->getGlobalNumElements(),goNumRows);
  TEST_EQUALITY(brop->getDomainMap()->getGlobalNumElements(),goNumRows);

  // block 00
  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> > brop00 =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> >(brop->getMatrix(0,0));

  GO goNumRows00 = Teuchos::as<GO>(Teuchos::ScalarTraits<GO>::pow(2,2)) * 10 * comm->getSize();

  TEST_EQUALITY(brop00->Rows(),2);
  TEST_EQUALITY(brop00->Cols(),2);
  TEST_EQUALITY(brop00->getRangeMap()->getGlobalNumElements(),goNumRows00);
  TEST_EQUALITY(brop00->getDomainMap()->getGlobalNumElements(),goNumRows00);

  // block 11
  Teuchos::RCP<const Xpetra::Matrix<Scalar,LO,GO,Node> > brop11 = brop->getMatrix(1,1);

  GO goNumRows11 = Teuchos::as<GO>(40 * comm->getSize());
  TEST_EQUALITY(brop11->getRangeMap()->getGlobalNumElements(),goNumRows11);
  TEST_EQUALITY(brop11->getDomainMap()->getGlobalNumElements(),goNumRows11);
  TEST_EQUALITY(brop11->getRangeMap()->getMinGlobalIndex(),comm->getRank() * 640 + 40);
  TEST_EQUALITY(brop11->getRangeMap()->getMaxGlobalIndex(),comm->getRank() * 640 + 79);

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> > brop11test =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> >(brop11);
  TEST_INEQUALITY(brop11test,Teuchos::null);

  TEST_EQUALITY(brop11test->Rows(), 1);
  TEST_EQUALITY(brop11test->Cols(), 1);
  TEST_EQUALITY(brop11test->getRangeMap()->getGlobalNumElements(),goNumRows11);
  TEST_EQUALITY(brop11test->getDomainMap()->getGlobalNumElements(),goNumRows11);
  TEST_EQUALITY(brop11test->getRangeMap()->getMinGlobalIndex(),comm->getRank() * 640 + 40);
  TEST_EQUALITY(brop11test->getRangeMap()->getMaxGlobalIndex(),comm->getRank() * 640 + 79);

  // block 22
  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> > brop22 =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> >(brop->getMatrix(2,2));

  GO goNumRows22 = Teuchos::as<GO>(560 * comm->getSize());

  TEST_EQUALITY(brop22->Rows(),3);
  TEST_EQUALITY(brop22->Cols(),3);
  TEST_EQUALITY(brop22->getRangeMap()->getGlobalNumElements(),goNumRows22);
  TEST_EQUALITY(brop22->getDomainMap()->getGlobalNumElements(),goNumRows22);
  TEST_EQUALITY(brop22->getRangeMap()->getMinGlobalIndex(),comm->getRank() * 640 + 80);
  TEST_EQUALITY(brop22->getRangeMap()->getMaxGlobalIndex(),comm->getRank() * 640 + 639);
  TEST_EQUALITY(brop22->getMatrix(0,0)->getRangeMap()->getMinGlobalIndex(),comm->getRank() * 640 + 80);
  TEST_EQUALITY(brop22->getMatrix(0,0)->getRangeMap()->getMaxGlobalIndex(),comm->getRank() * 640 + 159);
  TEST_EQUALITY(brop22->getMatrix(1,1)->getRangeMap()->getMinGlobalIndex(),comm->getRank() * 640 + 160);
  TEST_EQUALITY(brop22->getMatrix(1,1)->getRangeMap()->getMaxGlobalIndex(),comm->getRank() * 640 + 319);
  TEST_EQUALITY(brop22->getMatrix(2,2)->getRangeMap()->getMinGlobalIndex(),comm->getRank() * 640 + 320);
  TEST_EQUALITY(brop22->getMatrix(2,2)->getRangeMap()->getMaxGlobalIndex(),comm->getRank() * 640 + 639);

  // block 00_11
  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> > brop00_11 =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> >(brop00->getMatrix(1,1));

  GO goNumRows00_11 = Teuchos::as<GO>(35 * comm->getSize());

  TEST_EQUALITY(brop00_11->Rows(),2);
  TEST_EQUALITY(brop00_11->Cols(),2);
  TEST_EQUALITY(brop00_11->getRangeMap()->getGlobalNumElements(),goNumRows00_11);
  TEST_EQUALITY(brop00_11->getDomainMap()->getGlobalNumElements(),goNumRows00_11);
  TEST_EQUALITY(brop00_11->getRangeMap()->getMinGlobalIndex(),comm->getRank() * 640 + 5);
  TEST_EQUALITY(brop00_11->getRangeMap()->getMaxGlobalIndex(),comm->getRank() * 640 + 39);
  TEST_EQUALITY(brop00_11->getMatrix(0,0)->getRangeMap()->getGlobalNumElements(),15 * comm->getSize());
  TEST_EQUALITY(brop00_11->getMatrix(0,0)->getDomainMap()->getGlobalNumElements(),15 * comm->getSize());
  TEST_EQUALITY(brop00_11->getMatrix(0,0)->getRangeMap()->getMinGlobalIndex(),comm->getRank() * 640 + 5);
  TEST_EQUALITY(brop00_11->getMatrix(0,0)->getRangeMap()->getMaxGlobalIndex(),comm->getRank() * 640 + 19);
  TEST_EQUALITY(brop00_11->getMatrix(1,1)->getRangeMap()->getGlobalNumElements(),20 * comm->getSize());
  TEST_EQUALITY(brop00_11->getMatrix(1,1)->getDomainMap()->getGlobalNumElements(),20 * comm->getSize());
  TEST_EQUALITY(brop00_11->getMatrix(1,1)->getRangeMap()->getMinGlobalIndex(),comm->getRank() * 640 + 20);
  TEST_EQUALITY(brop00_11->getMatrix(1,1)->getRangeMap()->getMaxGlobalIndex(),comm->getRank() * 640 + 39);

  // block 01
  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> > brop01 =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> >(brop->getMatrix(0,1));

  TEST_EQUALITY(brop01->Rows(),2);
  TEST_EQUALITY(brop01->Cols(),1);
  TEST_EQUALITY(brop01->getRangeMap()->getGlobalNumElements(),comm->getSize() * 40);
  TEST_EQUALITY(brop01->getDomainMap()->getGlobalNumElements(),comm->getSize() * 40);
  TEST_EQUALITY(brop01->getRangeMap()->getMinGlobalIndex(),comm->getRank() * 640 + 0);
  TEST_EQUALITY(brop01->getRangeMap()->getMaxGlobalIndex(),comm->getRank() * 640 + 39);
  TEST_EQUALITY(brop01->getDomainMap()->getMinGlobalIndex(),comm->getRank() * 640 + 40);
  TEST_EQUALITY(brop01->getDomainMap()->getMaxGlobalIndex(),comm->getRank() * 640 + 79);

}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedCrsMatrix, ReorderBlockOperator2, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> BlockedCrsMatrixClass;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 8;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar,LO,GO,Node,M>(noBlocks, *comm);


  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 4 3 1 7 ]");

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> > brop =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> >(buildReorderedBlockedCrsMatrix(brm, bop));

  TEST_EQUALITY(brop->Rows(),4);
  TEST_EQUALITY(brop->Cols(),4);
  TEST_EQUALITY(brop->getRangeMap()->getGlobalNumElements(),comm->getSize() * 385);
  TEST_EQUALITY(brop->getDomainMap()->getGlobalNumElements(),comm->getSize() * 385);

  // block 00
  Teuchos::RCP<const Xpetra::Matrix<Scalar,LO,GO,Node> > brop00 = brop->getMatrix(0,0);

  TEST_EQUALITY(brop00->getRangeMap()->getGlobalNumElements(),comm->getSize() * 40);
  TEST_EQUALITY(brop00->getDomainMap()->getGlobalNumElements(),comm->getSize() * 40);
  TEST_EQUALITY(brop00->getRangeMap()->getMinGlobalIndex(),comm->getRank() * 640 + 40);
  TEST_EQUALITY(brop00->getRangeMap()->getMaxGlobalIndex(),comm->getRank() * 640 + 79);

  // block 11
  Teuchos::RCP<const Xpetra::Matrix<Scalar,LO,GO,Node> > brop11 = brop->getMatrix(1,1);

  TEST_EQUALITY(brop11->getRangeMap()->getGlobalNumElements(),comm->getSize() * 20);
  TEST_EQUALITY(brop11->getDomainMap()->getGlobalNumElements(),comm->getSize() * 20);
  TEST_EQUALITY(brop11->getRangeMap()->getMinGlobalIndex(),comm->getRank() * 640 + 20);
  TEST_EQUALITY(brop11->getRangeMap()->getMaxGlobalIndex(),comm->getRank() * 640 + 39);

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> > brop11test =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> >(brop11);

  TEST_EQUALITY(brop11test->Rows(), 1);
  TEST_EQUALITY(brop11test->Cols(), 1);
  TEST_EQUALITY(brop11test->getRangeMap()->getGlobalNumElements(),comm->getSize() * 20);
  TEST_EQUALITY(brop11test->getDomainMap()->getGlobalNumElements(),comm->getSize() * 20);
  TEST_EQUALITY(brop11test->getRangeMap()->getMinGlobalIndex(),comm->getRank() * 640 + 20);
  TEST_EQUALITY(brop11test->getRangeMap()->getMaxGlobalIndex(),comm->getRank() * 640 + 39);

  // block 22
  Teuchos::RCP<const Xpetra::Matrix<Scalar,LO,GO,Node> > brop22 = brop->getMatrix(2,2);

  TEST_EQUALITY(brop22->getRangeMap()->getGlobalNumElements(),comm->getSize() * 5);
  TEST_EQUALITY(brop22->getDomainMap()->getGlobalNumElements(),comm->getSize() * 5);
  TEST_EQUALITY(brop22->getRangeMap()->getMinGlobalIndex(),comm->getRank() * 640 + 5);
  TEST_EQUALITY(brop22->getRangeMap()->getMaxGlobalIndex(),comm->getRank() * 640 + 9);

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> > brop22test =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> >(brop22);

  TEST_EQUALITY(brop22test->Rows(), 1);
  TEST_EQUALITY(brop22test->Cols(), 1);
  TEST_EQUALITY(brop22test->getRangeMap()->getGlobalNumElements(),comm->getSize() * 5);
  TEST_EQUALITY(brop22test->getDomainMap()->getGlobalNumElements(),comm->getSize() * 5);
  TEST_EQUALITY(brop22test->getRangeMap()->getMinGlobalIndex(),comm->getRank() * 640 + 5);
  TEST_EQUALITY(brop22test->getRangeMap()->getMaxGlobalIndex(),comm->getRank() * 640 + 9);

  // block 33
  Teuchos::RCP<const Xpetra::Matrix<Scalar,LO,GO,Node> > brop33 = brop->getMatrix(3,3);

  TEST_EQUALITY(brop33->getRangeMap()->getGlobalNumElements(),comm->getSize() * 320);
  TEST_EQUALITY(brop33->getDomainMap()->getGlobalNumElements(),comm->getSize() * 320);
  TEST_EQUALITY(brop33->getRangeMap()->getMinGlobalIndex(),comm->getRank() * 640 + 320);
  TEST_EQUALITY(brop33->getRangeMap()->getMaxGlobalIndex(),comm->getRank() * 640 + 639);

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> > brop33test =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> >(brop33);

  TEST_EQUALITY(brop33test->Rows(), 1);
  TEST_EQUALITY(brop33test->Cols(), 1);
  TEST_EQUALITY(brop33test->getRangeMap()->getGlobalNumElements(),comm->getSize() * 320);
  TEST_EQUALITY(brop33test->getDomainMap()->getGlobalNumElements(),comm->getSize() * 320);
  TEST_EQUALITY(brop33test->getRangeMap()->getMinGlobalIndex(),comm->getRank() * 640 + 320);
  TEST_EQUALITY(brop33test->getRangeMap()->getMaxGlobalIndex(),comm->getRank() * 640 + 639);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedCrsMatrix, ReorderBlockOperatorThyra, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> BlockedCrsMatrixClass;
  typedef Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> ReorderedBlockedCrsMatrix;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 8;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrixThyra<Scalar,LO,GO,Node,M>(noBlocks, *comm);


  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ [ 0 [ [1 2] 3] ] 4 [ 5 6 7] ]");

  Teuchos::RCP<const ReorderedBlockedCrsMatrix> brop =
      Teuchos::rcp_dynamic_cast<const ReorderedBlockedCrsMatrix>(buildReorderedBlockedCrsMatrix(brm, bop));

  GO goNumRows = Teuchos::as<GO>(Teuchos::ScalarTraits<GO>::pow(2,noBlocks-2)) * 10 * comm->getSize();

  TEST_EQUALITY(brop->Rows(),3);
  TEST_EQUALITY(brop->Cols(),3);
  TEST_EQUALITY(brop->getRangeMap()->getGlobalNumElements(),goNumRows);
  TEST_EQUALITY(brop->getDomainMap()->getGlobalNumElements(),goNumRows);

  TEST_EQUALITY(brop->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop->getDomainMapExtractor()->getThyraMode(), true);

  // block 00
  Teuchos::RCP<const ReorderedBlockedCrsMatrix> brop00 =
      Teuchos::rcp_dynamic_cast<const ReorderedBlockedCrsMatrix>(brop->getMatrix(0,0));

  GO goNumRows00 = Teuchos::as<GO>(Teuchos::ScalarTraits<GO>::pow(2,2)) * 10 * comm->getSize();

  TEST_EQUALITY(brop00->Rows(),2);
  TEST_EQUALITY(brop00->Cols(),2);
  TEST_EQUALITY(brop00->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop00->getDomainMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop00->getRangeMap()->getGlobalNumElements(),goNumRows00);
  TEST_EQUALITY(brop00->getDomainMap()->getGlobalNumElements(),goNumRows00);
  TEST_EQUALITY(brop00->getRangeMap()->getMinGlobalIndex(),comm->getRank() * 5 + 0);
  TEST_EQUALITY(brop00->getRangeMap()->getMaxGlobalIndex(),comm->getSize() * 20 + comm->getRank() * 20 + 19);
  TEST_EQUALITY(brop00->getRangeMap()->getMinAllGlobalIndex(),0);
  TEST_EQUALITY(brop00->getRangeMap()->getMaxAllGlobalIndex(),comm->getSize() * 40 - 1);
  // Thyra maps (these might have duplicate GID entries!)
  TEST_EQUALITY(brop00->getRangeMap(0,true)->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(brop00->getRangeMap(0,true)->getMaxAllGlobalIndex(), comm->getSize()*5 - 1);
  TEST_EQUALITY(brop00->getRangeMap(1,true)->getGlobalNumElements(), comm->getSize() * 35);
  TEST_EQUALITY(brop00->getRangeMap(1,true)->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(brop00->getRangeMap(1,true)->getMaxAllGlobalIndex(), comm->getSize() * 20 - 1);
  // Xpetra maps
  TEST_EQUALITY(brop00->getRangeMap(0,false)->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(brop00->getRangeMap(0,false)->getMaxAllGlobalIndex(), comm->getSize()*5 - 1);
  TEST_EQUALITY(brop00->getRangeMap(1,false)->getGlobalNumElements(), comm->getSize() * 35);
  TEST_EQUALITY(brop00->getRangeMap(1,false)->getMinAllGlobalIndex(), comm->getSize()*5);
  TEST_EQUALITY(brop00->getRangeMap(1,false)->getMaxAllGlobalIndex(), comm->getSize()*5 + comm->getSize() * 35 - 1);


  // block 11
  Teuchos::RCP<const Xpetra::Matrix<Scalar,LO,GO,Node> > brop11 = brop->getMatrix(1,1);

  // Thyra GIDs for the matrix
  GO goNumRows11 = Teuchos::as<GO>(40 * comm->getSize());
  TEST_EQUALITY(brop11->getRangeMap()->getGlobalNumElements(),goNumRows11);
  TEST_EQUALITY(brop11->getDomainMap()->getGlobalNumElements(),goNumRows11);
  TEST_EQUALITY(brop11->getRangeMap()->getMinGlobalIndex(),comm->getSize() * 40 + comm->getRank() * 40 + 0);
  TEST_EQUALITY(brop11->getRangeMap()->getMaxGlobalIndex(),comm->getSize() * 40 + comm->getRank() * 40 + 39);
  TEST_EQUALITY(brop11->getRangeMap()->getMinAllGlobalIndex(),comm->getSize() * 40 + 0);
  TEST_EQUALITY(brop11->getRangeMap()->getMaxAllGlobalIndex(),comm->getSize() * 40 + comm->getSize() * 40 - 1);
  // Xpetra GIDs
  TEST_EQUALITY(brop->getRangeMap(1,false)->getMinAllGlobalIndex(),comm->getSize() * 40);
  TEST_EQUALITY(brop->getRangeMap(1,false)->getMaxAllGlobalIndex(),2 * comm->getSize() * 40 - 1);

  Teuchos::RCP<const ReorderedBlockedCrsMatrix> brop11test =
      Teuchos::rcp_dynamic_cast<const ReorderedBlockedCrsMatrix>(brop11);

  TEST_EQUALITY(brop11test->Rows(), 1);
  TEST_EQUALITY(brop11test->Cols(), 1);
  TEST_EQUALITY(brop11test->getRangeMap()->getMinGlobalIndex(),comm->getSize() * 40 + comm->getRank() * 40);
  TEST_EQUALITY(brop11test->getRangeMap()->getMaxGlobalIndex(),comm->getSize() * 40 + comm->getRank() * 40 + 39);
  TEST_EQUALITY(brop11test->getMatrix(0,0)->getRangeMap()->getMinGlobalIndex(),comm->getRank() * 40 + 0);
  TEST_EQUALITY(brop11test->getMatrix(0,0)->getRangeMap()->getMaxGlobalIndex(),comm->getRank() * 40 + 39);
  TEST_EQUALITY(brop11test->getMatrix(0,0)->getRangeMap()->getMinAllGlobalIndex(),0);
  TEST_EQUALITY(brop11test->getMatrix(0,0)->getRangeMap()->getMaxAllGlobalIndex(),comm->getSize() * 40 - 1);

  // block 22
  Teuchos::RCP<const ReorderedBlockedCrsMatrix> brop22 =
      Teuchos::rcp_dynamic_cast<const ReorderedBlockedCrsMatrix>(brop->getMatrix(2,2));

  GO goNumRows22 = Teuchos::as<GO>(560 * comm->getSize());

  TEST_EQUALITY(brop22->Rows(),3);
  TEST_EQUALITY(brop22->Cols(),3);
  TEST_EQUALITY(brop22->getRangeMap()->getGlobalNumElements(),goNumRows22);
  TEST_EQUALITY(brop22->getDomainMap()->getGlobalNumElements(),goNumRows22);
  // Xpetra replacement GIDs
  TEST_EQUALITY(brop22->getRangeMap()->getMinAllGlobalIndex(),comm->getSize() * 80);
  TEST_EQUALITY(brop22->getRangeMap()->getMaxAllGlobalIndex(),comm->getSize() * 80 + comm->getSize() * 560 - 1);
  TEST_EQUALITY(brop22->getRangeMap()->getMinGlobalIndex(),comm->getRank() * 80 + comm->getSize() * 80);
  TEST_EQUALITY(brop22->getRangeMap()->getMaxGlobalIndex(),comm->getSize() * 80 + comm->getSize() * 240 + comm->getRank() * 320 + 319);
  // Xpetra GIDs
  TEST_EQUALITY(brop22->getRangeMap(0,false)->getMinGlobalIndex(),comm->getSize() * 80 + comm->getRank() * 80);
  TEST_EQUALITY(brop22->getRangeMap(0,false)->getMaxGlobalIndex(),comm->getSize() * 80 + comm->getRank() * 80 + 79);
  TEST_EQUALITY(brop22->getRangeMap(1,false)->getMinGlobalIndex(),comm->getSize() * 160 + comm->getRank() * 160);
  TEST_EQUALITY(brop22->getRangeMap(1,false)->getMaxGlobalIndex(),comm->getSize() * 160 + comm->getRank() * 160 + 159);
  TEST_EQUALITY(brop22->getRangeMap(2,false)->getMinGlobalIndex(),comm->getSize() * 320 + comm->getRank() * 320);
  TEST_EQUALITY(brop22->getRangeMap(2,false)->getMaxGlobalIndex(),comm->getSize() * 320 + comm->getRank() * 320 + 319);

  // block 00_11
  /*Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> > brop00_11 =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> >(brop00->getMatrix(1,1));

  GO goNumRows00_11 = Teuchos::as<GO>(35 * comm->getSize());

  TEST_EQUALITY(brop00_11->Rows(),2);
  TEST_EQUALITY(brop00_11->Cols(),2);
  TEST_EQUALITY(brop00_11->getRangeMap()->getGlobalNumElements(),goNumRows00_11);
  TEST_EQUALITY(brop00_11->getDomainMap()->getGlobalNumElements(),goNumRows00_11);
  TEST_EQUALITY(brop00_11->getRangeMap()->getMinGlobalIndex(),comm->getRank() * 640 + 5);
  TEST_EQUALITY(brop00_11->getRangeMap()->getMaxGlobalIndex(),comm->getRank() * 640 + 39);
  TEST_EQUALITY(brop00_11->getMatrix(0,0)->getRangeMap()->getGlobalNumElements(),15 * comm->getSize());
  TEST_EQUALITY(brop00_11->getMatrix(0,0)->getDomainMap()->getGlobalNumElements(),15 * comm->getSize());
  TEST_EQUALITY(brop00_11->getMatrix(0,0)->getRangeMap()->getMinGlobalIndex(),comm->getRank() * 640 + 5);
  TEST_EQUALITY(brop00_11->getMatrix(0,0)->getRangeMap()->getMaxGlobalIndex(),comm->getRank() * 640 + 19);
  TEST_EQUALITY(brop00_11->getMatrix(1,1)->getRangeMap()->getGlobalNumElements(),20 * comm->getSize());
  TEST_EQUALITY(brop00_11->getMatrix(1,1)->getDomainMap()->getGlobalNumElements(),20 * comm->getSize());
  TEST_EQUALITY(brop00_11->getMatrix(1,1)->getRangeMap()->getMinGlobalIndex(),comm->getRank() * 640 + 20);
  TEST_EQUALITY(brop00_11->getMatrix(1,1)->getRangeMap()->getMaxGlobalIndex(),comm->getRank() * 640 + 39);
*/
  // block 01
  Teuchos::RCP<const ReorderedBlockedCrsMatrix> brop01 =
      Teuchos::rcp_dynamic_cast<const ReorderedBlockedCrsMatrix>(brop->getMatrix(0,1));

  // Xpetra like maps
  TEST_EQUALITY(brop01->Rows(),2);
  TEST_EQUALITY(brop01->Cols(),1);
  TEST_EQUALITY(brop01->getRangeMap()->getGlobalNumElements(),comm->getSize() * 40);
  TEST_EQUALITY(brop01->getDomainMap()->getGlobalNumElements(),comm->getSize() * 40);
  TEST_EQUALITY(brop01->getRangeMap()->getMinGlobalIndex(),comm->getRank() * 5 + 0);
  TEST_EQUALITY(brop01->getRangeMap()->getMaxGlobalIndex(),comm->getSize() * 5 + comm->getSize() * 15 + comm->getRank() * 20 + 19);
  TEST_EQUALITY(brop01->getDomainMap()->getMinGlobalIndex(),comm->getSize() * 40 + comm->getRank() * 40);
  TEST_EQUALITY(brop01->getDomainMap()->getMaxGlobalIndex(),comm->getSize() * 40 + comm->getRank() * 40 + 39);

}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedCrsMatrix, ReorderBlockOperator2Thyra, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> BlockedCrsMatrixClass;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 8;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrixThyra<Scalar,LO,GO,Node,M>(noBlocks, *comm);

  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 4 3 1 7 ]");

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> > brop =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> >(buildReorderedBlockedCrsMatrix(brm, bop));

  TEST_EQUALITY(brop->Rows(),4);
  TEST_EQUALITY(brop->Cols(),4);
  TEST_EQUALITY(brop->getRangeMap()->getGlobalNumElements(),comm->getSize() * 385);
  TEST_EQUALITY(brop->getDomainMap()->getGlobalNumElements(),comm->getSize() * 385);

  // block 00
  Teuchos::RCP<const Xpetra::Matrix<Scalar,LO,GO,Node> > brop00 = brop->getMatrix(0,0);

  TEST_EQUALITY(brop00->getRangeMap()->getGlobalNumElements(),comm->getSize() * 40);
  TEST_EQUALITY(brop00->getDomainMap()->getGlobalNumElements(),comm->getSize() * 40);
  // Thyra GIDs
  TEST_EQUALITY(brop00->getRangeMap()->getMinGlobalIndex(),comm->getSize() * 40 + comm->getRank() * 40);
  TEST_EQUALITY(brop00->getRangeMap()->getMaxGlobalIndex(),comm->getSize() * 40 + comm->getRank() * 40 + 39);
  TEST_EQUALITY(brop00->getRangeMap()->getMinAllGlobalIndex(),comm->getSize() * 40);
  TEST_EQUALITY(brop00->getRangeMap()->getMaxAllGlobalIndex(),comm->getSize() * 40 + comm->getSize() * 40 - 1);
  // Xpetra GIDs
  TEST_EQUALITY(brop->getDomainMap(0,false)->getMinGlobalIndex(), comm->getSize() * 40 + comm->getRank() * 40);
  TEST_EQUALITY(brop->getDomainMap(0,false)->getMaxGlobalIndex(), comm->getSize() * 40 + comm->getRank() * 40 + 39);

  // block 11
  Teuchos::RCP<const Xpetra::Matrix<Scalar,LO,GO,Node> > brop11 = brop->getMatrix(1,1);

  TEST_EQUALITY(brop11->getRangeMap()->getGlobalNumElements(),comm->getSize() * 20);
  TEST_EQUALITY(brop11->getDomainMap()->getGlobalNumElements(),comm->getSize() * 20);
  // Thyra GIDs (+ Xpetra shift)
  TEST_EQUALITY(brop11->getRangeMap()->getMinGlobalIndex(),comm->getSize() * 20 + comm->getRank() * 20);
  TEST_EQUALITY(brop11->getRangeMap()->getMaxGlobalIndex(),comm->getSize() * 20 + comm->getRank() * 20 + 19);
  TEST_EQUALITY(brop11->getRangeMap()->getMinAllGlobalIndex(),comm->getSize() * 20 + 0);
  TEST_EQUALITY(brop11->getRangeMap()->getMaxAllGlobalIndex(),comm->getSize() * 20 + comm->getSize() * 20 - 1);

  // Xpetra GIDs
  TEST_EQUALITY(brop->getDomainMap(1,false)->getMinGlobalIndex(), comm->getSize() * 20 + comm->getRank() * 20);
  TEST_EQUALITY(brop->getDomainMap(1,false)->getMaxGlobalIndex(), comm->getSize() * 20 + comm->getRank() * 20 + 19);


  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> > brop11test =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> >(brop11);

  TEST_EQUALITY(brop11test->Rows(), 1);
  TEST_EQUALITY(brop11test->Cols(), 1);
  // Thyra GIDs
  TEST_EQUALITY(brop11test->getMatrix(0,0)->getRangeMap()->getMinGlobalIndex(),comm->getRank() * 20);
  TEST_EQUALITY(brop11test->getMatrix(0,0)->getRangeMap()->getMaxGlobalIndex(),comm->getRank() * 20 + 19);
  TEST_EQUALITY(brop11test->getMatrix(0,0)->getRangeMap()->getMinAllGlobalIndex(),0);
  TEST_EQUALITY(brop11test->getMatrix(0,0)->getRangeMap()->getMaxAllGlobalIndex(),comm->getSize() * 20 - 1);

  // block 22
  Teuchos::RCP<const Xpetra::Matrix<Scalar,LO,GO,Node> > brop22 = brop->getMatrix(2,2);

  TEST_EQUALITY(brop22->getRangeMap()->getGlobalNumElements(),comm->getSize() * 5);
  TEST_EQUALITY(brop22->getDomainMap()->getGlobalNumElements(),comm->getSize() * 5);
  TEST_EQUALITY(brop22->getRangeMap()->getMinGlobalIndex(),comm->getSize() * 5 + comm->getRank() * 5);
  TEST_EQUALITY(brop22->getRangeMap()->getMaxGlobalIndex(),comm->getSize() * 5 + comm->getRank() * 5 + 4);
  TEST_EQUALITY(brop22->getRangeMap()->getMinAllGlobalIndex(),comm->getSize() * 5 + 0);
  TEST_EQUALITY(brop22->getRangeMap()->getMaxAllGlobalIndex(),comm->getSize() * 5 + comm->getSize() * 5 - 1);
  // Xpetra GIDs
  TEST_EQUALITY(brop->getDomainMap(2,false)->getMinGlobalIndex(), comm->getSize() * 5 + comm->getRank() * 5);
  TEST_EQUALITY(brop->getDomainMap(2,false)->getMaxGlobalIndex(), comm->getSize() * 5 + comm->getRank() * 5 + 4);

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> > brop22test =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> >(brop22);
  TEST_EQUALITY(brop22test->Rows(), 1);
  TEST_EQUALITY(brop22test->Cols(), 1);

  // Thyra GIDs
  TEST_EQUALITY(brop22test->getMatrix(0,0)->getRangeMap()->getMinGlobalIndex(),comm->getRank() * 5);
  TEST_EQUALITY(brop22test->getMatrix(0,0)->getRangeMap()->getMaxGlobalIndex(),comm->getRank() * 5 + 4);
  TEST_EQUALITY(brop22test->getMatrix(0,0)->getRangeMap()->getMinAllGlobalIndex(),0);
  TEST_EQUALITY(brop22test->getMatrix(0,0)->getRangeMap()->getMaxAllGlobalIndex(),comm->getSize() * 5 - 1);


  // block 33
  Teuchos::RCP<const Xpetra::Matrix<Scalar,LO,GO,Node> > brop33 = brop->getMatrix(3,3);

  TEST_EQUALITY(brop33->getRangeMap()->getGlobalNumElements(),comm->getSize() * 320);
  TEST_EQUALITY(brop33->getDomainMap()->getGlobalNumElements(),comm->getSize() * 320);
  TEST_EQUALITY(brop33->getRangeMap()->getMinGlobalIndex(),comm->getSize() * 320 + comm->getRank() * 320);
  TEST_EQUALITY(brop33->getRangeMap()->getMaxGlobalIndex(),comm->getSize() * 320 + comm->getRank() * 320 + 319);
  TEST_EQUALITY(brop33->getRangeMap()->getMinAllGlobalIndex(),comm->getSize() * 320 + 0);
  TEST_EQUALITY(brop33->getRangeMap()->getMaxAllGlobalIndex(),comm->getSize() * 320 + comm->getSize() * 320 - 1);
  // Xpetra GIDs
  TEST_EQUALITY(brop->getDomainMap(3,false)->getMinGlobalIndex(), comm->getSize() * 320 + comm->getRank() * 320);
  TEST_EQUALITY(brop->getDomainMap(3,false)->getMaxGlobalIndex(), comm->getSize() * 320 + comm->getRank() * 320 + 319);

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> > brop33test =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> >(brop33);
  TEST_EQUALITY(brop33test->Rows(), 1);
  TEST_EQUALITY(brop33test->Cols(), 1);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedCrsMatrix, ReorderBlockOperatorApply, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::Vector<Scalar, LO, GO, Node> VectorClass;
  typedef Xpetra::VectorFactory<Scalar, LO, GO, Node> VectorFactoryClass;
  typedef Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> BlockedCrsMatrixClass;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 8;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar,LO,GO,Node,M>(noBlocks, *comm);


  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 0 [ 1 [ 2 3 4 ] 5 ] [6 7] ]");

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> > brop =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> >(buildReorderedBlockedCrsMatrix(brm, bop));

  TEST_EQUALITY(brop->Rows(),3);
  TEST_EQUALITY(brop->Cols(),3);
  TEST_EQUALITY(brop->getRangeMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(brop->getDomainMapExtractor()->getThyraMode(), false);

  // build gloabl vector with one entries
  Teuchos::RCP<VectorClass> ones = VectorFactoryClass::Build(bop->getRangeMap(), true);
  Teuchos::RCP<VectorClass> exp  = VectorFactoryClass::Build(bop->getRangeMap(), true);
  Teuchos::RCP<VectorClass> res  = VectorFactoryClass::Build(bop->getRangeMap(), true);
  Teuchos::RCP<VectorClass> rnd  = VectorFactoryClass::Build(bop->getRangeMap(), true);
  ones->putScalar(STS::one());
  rnd->randomize();

  bop->apply(*ones, *exp);
  brop->apply(*ones, *res);
  res->update(-STS::one(),*exp,STS::one());
  TEUCHOS_TEST_COMPARE(res->norm2(), <, 1e-16, out, success);
  TEUCHOS_TEST_COMPARE(res->normInf(), <, 1e-16, out, success);

  bop->apply(*rnd, *exp);
  brop->apply(*rnd, *res);
  res->update(-STS::one(),*exp,STS::one());
  TEUCHOS_TEST_COMPARE(res->norm2(), <, 5e-14, out, success);
  TEUCHOS_TEST_COMPARE(res->normInf(), <, 5e-14, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedCrsMatrix, ReorderBlockOperatorApply2, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::Vector<Scalar, LO, GO, Node> VectorClass;
  typedef Xpetra::VectorFactory<Scalar, LO, GO, Node> VectorFactoryClass;
  typedef Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> BlockedCrsMatrixClass;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 8;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar,LO,GO,Node,M>(noBlocks, *comm);


  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 6 3 2 ]");

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> > brop =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> >(buildReorderedBlockedCrsMatrix(brm, bop));

  TEST_EQUALITY(brop->Rows(),3);
  TEST_EQUALITY(brop->Cols(),3);
  TEST_EQUALITY(brop->getRangeMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(brop->getDomainMapExtractor()->getThyraMode(), false);

  // build gloabl vector with one entries
  Teuchos::RCP<VectorClass> ones = VectorFactoryClass::Build(brop->getRangeMap(), true);
  Teuchos::RCP<VectorClass> res  = VectorFactoryClass::Build(brop->getRangeMap(), true);
  ones->putScalar(STS::one());

  brop->apply(*ones, *res);

  TEST_EQUALITY(res->norm1(), STS::magnitude(Teuchos::as<Scalar>(comm->getSize()) * Teuchos::as<Scalar>(1230)));
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedCrsMatrix, ReorderBlockOperatorApplyThyra, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::Vector<Scalar, LO, GO, Node> VectorClass;
  typedef Xpetra::VectorFactory<Scalar, LO, GO, Node> VectorFactoryClass;
  typedef Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> BlockedCrsMatrixClass;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 8;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrixThyra<Scalar,LO,GO,Node,M>(noBlocks, *comm);

  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 0 [ 1 [ 2 3 4 ] 5 ] [6 7] ]");

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> > brop =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> >(buildReorderedBlockedCrsMatrix(brm, bop));

  TEST_EQUALITY(brop->Rows(),3);
  TEST_EQUALITY(brop->Cols(),3);
  TEST_EQUALITY(brop->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop->getDomainMapExtractor()->getThyraMode(), true);

  // build gloabl vector with one entries
  Teuchos::RCP<VectorClass> ones = VectorFactoryClass::Build(bop->getRangeMap(), true);
  Teuchos::RCP<VectorClass> exp  = VectorFactoryClass::Build(bop->getRangeMap(), true);
  Teuchos::RCP<VectorClass> res  = VectorFactoryClass::Build(bop->getRangeMap(), true);
  Teuchos::RCP<VectorClass> rnd  = VectorFactoryClass::Build(bop->getRangeMap(), true);
  ones->putScalar(STS::one());
  rnd->randomize();

  bop->apply(*ones, *exp);
  brop->apply(*ones, *res);
  res->update(-STS::one(),*exp,STS::one());
  TEUCHOS_TEST_COMPARE(res->norm2(), <, 1e-16, out, success);
  TEUCHOS_TEST_COMPARE(res->normInf(), <, 1e-16, out, success);

  bop->apply(*rnd, *exp);
  brop->apply(*rnd, *res);
  res->update(-STS::one(),*exp,STS::one());
  TEUCHOS_TEST_COMPARE(res->norm2(), <, 5e-14, out, success);
  TEUCHOS_TEST_COMPARE(res->normInf(), <, 5e-14, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedCrsMatrix, ReorderBlockOperatorApply2Thyra, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::Vector<Scalar, LO, GO, Node> VectorClass;
  typedef Xpetra::VectorFactory<Scalar, LO, GO, Node> VectorFactoryClass;
  typedef Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> BlockedCrsMatrixClass;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 8;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrixThyra<Scalar,LO,GO,Node,M>(noBlocks, *comm);


  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 6 3 2 ]");

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> > brop =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> >(buildReorderedBlockedCrsMatrix(brm, bop));

  TEST_EQUALITY(brop->Rows(),3);
  TEST_EQUALITY(brop->Cols(),3);
  TEST_EQUALITY(brop->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop->getDomainMapExtractor()->getThyraMode(), true);

  // build gloabl vector with one entries
  Teuchos::RCP<VectorClass> ones = VectorFactoryClass::Build(brop->getRangeMap(), true);
  Teuchos::RCP<VectorClass> res  = VectorFactoryClass::Build(brop->getRangeMap(), true);
  ones->putScalar(STS::one());

  brop->apply(*ones, *res);

  TEST_EQUALITY(res->norm1(), STS::magnitude(Teuchos::as<Scalar>(comm->getSize()) * Teuchos::as<Scalar>(1230)));
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedCrsMatrix, ReadWriteBlockedMatrix, M, MA, Scalar, LO, GO, Node )
{
  // TODO: it seems that the Tpetra matrix reader is only working for standard maps??

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

  LO NumMyElements = map->getNodeNumElements();
  GO NumGlobalElements = map->getGlobalNumElements();
  Teuchos::ArrayView<const GO> MyGlobalElements = map->getNodeElementList();

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
      Xpetra::MatrixUtils<Scalar, LO, GO, Node>::SplitMatrix(*mat,map_extractor,map_extractor);

  // Write matrices out, read fine A back in, and check that the read was ok
  // by using a matvec with a random vector.
  // JJH: 22-Feb-2016 Append scalar type to file name. The theory is that for dashboard
  //      tests with multiple Scalar instantiations of this test, a test with Scalar type
  //      A could try to read in the results of the test with Scalar type B, simply because
  //      the test with type B overwrote A's output matrix file.  A better solution would be
  //      to write to a file stream, but this would involve writing new interfaces to Epetra's
  //      file I/O capabilities.
  std::string tname = "BLOCKEDMATRIX";
  tname = tname + typeid(Scalar).name();
  tname = tname + typeid(LO).name();
  tname = tname + typeid(GO).name();
#ifdef HAVE_MUELU_KOKKOSCORE
  std::string nn = Kokkos::Compat::KokkosDeviceWrapperNode<typename Node::execution_space>::name();
  nn.erase(std::remove(nn.begin(), nn.end(), '/'), nn.end());
  tname = tname + nn;
#endif
  tname = "_" + tname;

  Xpetra::IO<Scalar, LO, GO, Node>::WriteBlockedCrsMatrix(tname, *bMat);
  Teuchos::RCP<const Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> > bMat2 = Xpetra::IO<Scalar, LO, GO, Node>::ReadBlockedCrsMatrix(tname, lib, comm);

  TEST_EQUALITY(bMat->getMatrix(0,0)->getGlobalNumEntries(),bMat2->getMatrix(0,0)->getGlobalNumEntries());
  TEST_EQUALITY(bMat->getMatrix(0,1)->getGlobalNumEntries(),bMat2->getMatrix(0,1)->getGlobalNumEntries());
  TEST_EQUALITY(bMat->getMatrix(1,0)->getGlobalNumEntries(),bMat2->getMatrix(1,0)->getGlobalNumEntries());
  TEST_EQUALITY(bMat->getMatrix(1,1)->getGlobalNumEntries(),bMat2->getMatrix(1,1)->getGlobalNumEntries());

  TEST_EQUALITY(bMat->getMatrix(0,0)->getNodeNumEntries(),bMat2->getMatrix(0,0)->getNodeNumEntries());
  TEST_EQUALITY(bMat->getMatrix(0,1)->getNodeNumEntries(),bMat2->getMatrix(0,1)->getNodeNumEntries());
  TEST_EQUALITY(bMat->getMatrix(1,0)->getNodeNumEntries(),bMat2->getMatrix(1,0)->getNodeNumEntries());
  TEST_EQUALITY(bMat->getMatrix(1,1)->getNodeNumEntries(),bMat2->getMatrix(1,1)->getNodeNumEntries());

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
  TEST_EQUALITY(bMat->getRangeMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(bMat->getDomainMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(bMat->getRangeMap(0)->getMinAllGlobalIndex(), 0);
  TEST_INEQUALITY(bMat->getRangeMap(1)->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(bMat->getDomainMap(0)->getMinAllGlobalIndex(), 0);
  TEST_INEQUALITY(bMat->getDomainMap(1)->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(bMat2->getRangeMapExtractor()->getThyraMode(), false);  // thyra mode is not correctly transferred!!
  TEST_EQUALITY(bMat2->getDomainMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(bMat2->getRangeMap(0)->getMinAllGlobalIndex(), 0);
  TEST_INEQUALITY(bMat2->getRangeMap(1)->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(bMat2->getDomainMap(0)->getMinAllGlobalIndex(), 0);
  TEST_INEQUALITY(bMat2->getDomainMap(1)->getMinAllGlobalIndex(), 0);
}

/// simple test routine for the apply function of BlockedCrsMatrix
TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedCrsMatrix, Apply, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::Map<LO, GO, Node> MapClass;
  typedef Xpetra::StridedMap<LO, GO, Node> StridedMapClass;
  typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;
  typedef Xpetra::StridedMapFactory<LO, GO, Node> StridedMapFactoryClass;
  typedef Xpetra::Vector<Scalar, LO, GO, Node> VectorClass;
  typedef Xpetra::VectorFactory<Scalar, LO, GO, Node> VectorFactoryClass;
  typedef Xpetra::MapExtractor<Scalar,LO,GO,Node> MapExtractorClass;
  typedef Xpetra::MapExtractorFactory<Scalar,LO,GO,Node> MapExtractorFactoryClass;
  typedef Xpetra::Matrix<Scalar, LO, GO, Node> MatrixClass;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  M testMap(1,0,comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  const Teuchos::RCP<const MapClass> pointmap = MapFactoryClass::Build(lib, 12, 0, comm);

  // generate local maps for loading matrices
  Teuchos::Array<GO> velgidvec; // global strided maps
  Teuchos::Array<GO> pregidvec;
  Teuchos::Array<GO> fullgidvec; // full global map
  for (LO i=0; i<Teuchos::as<LO>(pointmap->getNodeNumElements()); i++)
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

  // build gloabl vector with one entries
  Teuchos::RCP<VectorClass> ones = VectorFactoryClass::Build(fullmap, true);
  Teuchos::RCP<VectorClass> exp  = VectorFactoryClass::Build(fullmap, true);
  Teuchos::RCP<VectorClass> res  = VectorFactoryClass::Build(fullmap, true);
  Teuchos::RCP<VectorClass> rnd  = VectorFactoryClass::Build(fullmap, true);
  ones->putScalar(STS::one());
  rnd->randomize();

  A->apply(*ones, *exp);
  bOp->apply(*ones, *res);
  res->update(-STS::one(),*exp,STS::one());
  TEUCHOS_TEST_COMPARE(res->norm2(), <, 1e-16, out, success);
  TEUCHOS_TEST_COMPARE(res->normInf(), <, 1e-16, out, success);

  A->apply(*rnd, *exp);
  bOp->apply(*rnd, *res);
  res->update(-STS::one(),*exp,STS::one());
  TEUCHOS_TEST_COMPARE(res->norm2(), <, 5e-14, out, success);
  TEUCHOS_TEST_COMPARE(res->normInf(), <, 5e-14, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedCrsMatrix, getLocalDiagCopy, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> BlockedCrsMatrixClass;
  typedef Xpetra::Vector<Scalar,LO,GO,Node> VectorClass;
  typedef Xpetra::VectorFactory<Scalar,LO,GO,Node> VectorFactoryClass;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 5;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar,LO,GO,Node,M>(noBlocks, *comm);
  TEST_EQUALITY(bop->Rows(), noBlocks);
  TEST_EQUALITY(bop->Cols(), noBlocks);


  Teuchos::RCP<VectorClass> vorig = VectorFactoryClass::Build(bop->getRangeMap(),true);

  bop->getLocalDiagCopy(*vorig);

  Teuchos::ArrayRCP< const Scalar > vdataorig = vorig->getData(0);
  bool bCheck = true;
  for(int i=0; i<5; i++)  if(vdataorig[i] != Teuchos::as<Scalar>(1.0)) bCheck = false;
  for(int i=5; i<10; i++) if(vdataorig[i] != Teuchos::as<Scalar>(2.0)) bCheck = false;
  for(int i=10; i<20; i++) if(vdataorig[i] != Teuchos::as<Scalar>(3.0)) bCheck = false;
  for(int i=20; i<40; i++) if(vdataorig[i] != Teuchos::as<Scalar>(4.0)) bCheck = false;
  for(int i=40; i<80; i++) if(vdataorig[i] != Teuchos::as<Scalar>(5.0)) bCheck = false;
  TEST_EQUALITY(bCheck, true);

  // reordered blocked operator (Xpetra style)
  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 4 [3 2] 1 0]");

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> > brop =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> >(buildReorderedBlockedCrsMatrix(brm, bop));

  GO goNumRows = Teuchos::as<GO>(Teuchos::ScalarTraits<GO>::pow(2,noBlocks-2)) * 10 * comm->getSize();

  TEST_EQUALITY(brop->Rows(), 4);
  TEST_EQUALITY(brop->Cols(), 4);
  TEST_EQUALITY(brop->getRangeMap()->getGlobalNumElements(),goNumRows);
  TEST_EQUALITY(brop->getDomainMap()->getGlobalNumElements(),goNumRows);

  Teuchos::RCP<VectorClass> v = VectorFactoryClass::Build(brop->getRangeMap(),true);

  brop->getLocalDiagCopy(*v);

  Teuchos::ArrayRCP< const Scalar > vdata = v->getData(0);
  bCheck = true;
  for(int i=0; i<40; i++)  if(vdata[i] != Teuchos::as<Scalar>(5.0)) bCheck = false;
  for(int i=40; i<60; i++) if(vdata[i] != Teuchos::as<Scalar>(4.0)) bCheck = false;
  for(int i=60; i<70; i++) if(vdata[i] != Teuchos::as<Scalar>(3.0)) bCheck = false;
  for(int i=70; i<75; i++) if(vdata[i] != Teuchos::as<Scalar>(2.0)) bCheck = false;
  for(int i=75; i<80; i++) if(vdata[i] != Teuchos::as<Scalar>(1.0)) bCheck = false;
  TEST_EQUALITY(bCheck, true);

  // Thyra style (reordered) operator
  Teuchos::RCP<const BlockedCrsMatrixClass> btop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrixThyra<Scalar,LO,GO,Node,M>(noBlocks, *comm);

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> > brtop =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> >(buildReorderedBlockedCrsMatrix(brm, btop));

  TEST_EQUALITY(brtop->Rows(), 4);
  TEST_EQUALITY(brtop->Cols(), 4);
  TEST_EQUALITY(brtop->getRangeMap()->getGlobalNumElements(),goNumRows);
  TEST_EQUALITY(brtop->getDomainMap()->getGlobalNumElements(),goNumRows);

  Teuchos::RCP<VectorClass> v2 = VectorFactoryClass::Build(brtop->getRangeMap(),true);

  brtop->getLocalDiagCopy(*v2);

  Teuchos::ArrayRCP< const Scalar > vdata2 = v2->getData(0);
  bCheck = true;
  for(int i=0; i<40; i++)  if(vdata2[i] != Teuchos::as<Scalar>(5.0)) bCheck = false;
  for(int i=40; i<60; i++) if(vdata2[i] != Teuchos::as<Scalar>(4.0)) bCheck = false;
  for(int i=60; i<70; i++) if(vdata2[i] != Teuchos::as<Scalar>(3.0)) bCheck = false;
  for(int i=70; i<75; i++) if(vdata2[i] != Teuchos::as<Scalar>(2.0)) bCheck = false;
  for(int i=75; i<80; i++) if(vdata2[i] != Teuchos::as<Scalar>(1.0)) bCheck = false;
  TEST_EQUALITY(bCheck, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedCrsMatrix, deepCopy, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> BlockedCrsMatrixClass;
  typedef Xpetra::Vector<Scalar,LO,GO,Node> VectorClass;
  typedef Xpetra::VectorFactory<Scalar,LO,GO,Node> VectorFactoryClass;
  typedef Xpetra::Matrix<Scalar,LO,GO,Node> MatrixClass;
  typedef Xpetra::MatrixFactory<Scalar,LO,GO,Node> MatrixFactoryClass;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 5;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar,LO,GO,Node,M>(noBlocks, *comm);
  TEST_EQUALITY(bop->Rows(), noBlocks);
  TEST_EQUALITY(bop->Cols(), noBlocks);

  Teuchos::RCP<const MatrixClass> A2 = MatrixFactoryClass::BuildCopy(bop);
  Teuchos::RCP<const BlockedCrsMatrixClass> bop2 = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrixClass>(A2);
  TEST_EQUALITY(bop2->Rows(), noBlocks);
  TEST_EQUALITY(bop2->Cols(), noBlocks);

  TEST_EQUALITY(bop2->getRangeMapExtractor()->NumMaps(), bop->getRangeMapExtractor()->NumMaps());
  TEST_EQUALITY(bop2->getDomainMapExtractor()->NumMaps(), bop->getDomainMapExtractor()->NumMaps());
  TEST_EQUALITY(bop2->getGlobalMaxNumRowEntries(), bop->getGlobalMaxNumRowEntries());
  TEST_EQUALITY(bop2->getGlobalNumEntries(), bop->getGlobalNumEntries());
  TEST_EQUALITY(bop2->getGlobalNumRows(), bop->getGlobalNumRows());
  TEST_EQUALITY(bop2->getGlobalNumCols(), bop->getGlobalNumCols());

  Teuchos::RCP<VectorClass> v1 = VectorFactoryClass::Build(bop->getRangeMap(),true);
  Teuchos::RCP<VectorClass> v2 = VectorFactoryClass::Build(bop2->getRangeMap(),true);
  bop->getLocalDiagCopy(*v1);
  bop2->getLocalDiagCopy(*v2);

  v1->update(-Teuchos::ScalarTraits<Scalar>::one(),*v2,Teuchos::ScalarTraits<Scalar>::one());
  TEUCHOS_TEST_COMPARE(v1->norm2(), <, 1e-16, out, success);
  TEUCHOS_TEST_COMPARE(v1->normInf(), <, 1e-16, out, success);

  v1 = Teuchos::null;
  v2 = Teuchos::null;

  bop = Teuchos::null;

  TEST_EQUALITY(bop2->getRangeMapExtractor()->NumMaps(), noBlocks);
  TEST_EQUALITY(bop2->getDomainMapExtractor()->NumMaps(), noBlocks);
  TEST_EQUALITY(bop2->getRangeMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(bop2->getDomainMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(bop2->getGlobalMaxNumRowEntries(), 1);
  TEST_EQUALITY(bop2->getGlobalNumRows(), comm->getSize() * 80);
  TEST_EQUALITY(bop2->getGlobalNumCols(), comm->getSize() * 80);
  TEST_EQUALITY(bop2->getMatrix(0,0)!=Teuchos::null, true);
  TEST_EQUALITY(bop2->getMatrix(0,0)->getRangeMap()->getGlobalNumElements(), comm->getSize() * 5);

  // Thyra blocked operator
  Teuchos::RCP<const BlockedCrsMatrixClass> bop3 = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrixThyra<Scalar,LO,GO,Node,M>(noBlocks, *comm);
  TEST_EQUALITY(bop3->Rows(), noBlocks);
  TEST_EQUALITY(bop3->Cols(), noBlocks);

  Teuchos::RCP<const MatrixClass> A4 = MatrixFactoryClass::BuildCopy(bop3);
  Teuchos::RCP<const BlockedCrsMatrixClass> bop4 = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrixClass>(A4);
  TEST_EQUALITY(bop4->Rows(), noBlocks);
  TEST_EQUALITY(bop4->Cols(), noBlocks);

  bop3 = Teuchos::null;

  TEST_EQUALITY(bop4->getRangeMapExtractor()->NumMaps(), noBlocks);
  TEST_EQUALITY(bop4->getDomainMapExtractor()->NumMaps(), noBlocks);
  TEST_EQUALITY(bop4->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(bop4->getDomainMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(bop4->getGlobalMaxNumRowEntries(), 1);
  TEST_EQUALITY(bop4->getGlobalNumRows(), comm->getSize() * 80);
  TEST_EQUALITY(bop4->getGlobalNumCols(), comm->getSize() * 80);
  TEST_EQUALITY(bop4->getMatrix(0,0)!=Teuchos::null, true);
  TEST_EQUALITY(bop4->getMatrix(0,0)->getRangeMap()->getGlobalNumElements(), comm->getSize() * 5);

  // Nested Xpetra blocked operator
  bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar,LO,GO,Node,M>(noBlocks, *comm);
  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 1 [ [ 2 4 0 ] 3] ]");

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> > brop =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> >(buildReorderedBlockedCrsMatrix(brm, bop));

  TEST_EQUALITY(brop->Rows(), 2);
  TEST_EQUALITY(brop->Cols(), 2);

  Teuchos::RCP<const MatrixClass> A = MatrixFactoryClass::BuildCopy(brop);
  Teuchos::RCP<const BlockedCrsMatrixClass> brop2 = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrixClass>(A);
  brop = Teuchos::null;
  TEST_EQUALITY(brop2->Rows(), 2);
  TEST_EQUALITY(brop2->Cols(), 2);
  Teuchos::RCP<const BlockedCrsMatrixClass> brop200 = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrixClass>(brop2->getMatrix(0,0));
  Teuchos::RCP<const BlockedCrsMatrixClass> brop211 = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrixClass>(brop2->getMatrix(1,1));
  Teuchos::RCP<const BlockedCrsMatrixClass> brop21100 = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrixClass>(brop211->getMatrix(0,0));
  TEST_EQUALITY(brop200->Rows(), 1);
  TEST_EQUALITY(brop200->Cols(), 1);
  TEST_EQUALITY(brop211->Rows(), 2);
  TEST_EQUALITY(brop211->Cols(), 2);
  TEST_EQUALITY(brop21100->Rows(), 3);
  TEST_EQUALITY(brop21100->Cols(), 3);
  TEST_EQUALITY(brop21100->getMatrix(0,0)->getRangeMap()->getGlobalNumElements(),comm->getSize() * 10);
  TEST_EQUALITY(brop21100->getMatrix(1,1)->getRangeMap()->getGlobalNumElements(),comm->getSize() * 40);
  TEST_EQUALITY(brop2->getRangeMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(brop2->getDomainMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(brop200->getRangeMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(brop200->getDomainMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(brop211->getRangeMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(brop211->getDomainMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(brop21100->getRangeMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(brop21100->getDomainMapExtractor()->getThyraMode(), false);

  // Nested Thyra blocked operator
  bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrixThyra<Scalar,LO,GO,Node,M>(noBlocks, *comm);

  brop = Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> >(buildReorderedBlockedCrsMatrix(brm, bop));

  TEST_EQUALITY(brop->Rows(), 2);
  TEST_EQUALITY(brop->Cols(), 2);

  A = MatrixFactoryClass::BuildCopy(brop);
  brop2 = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrixClass>(A);
  brop = Teuchos::null;
  TEST_EQUALITY(brop2->Rows(), 2);
  TEST_EQUALITY(brop2->Cols(), 2);
  brop200 = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrixClass>(brop2->getMatrix(0,0));
  brop211 = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrixClass>(brop2->getMatrix(1,1));
  brop21100 = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrixClass>(brop211->getMatrix(0,0));
  TEST_EQUALITY(brop200->Rows(), 1);
  TEST_EQUALITY(brop200->Cols(), 1);
  TEST_EQUALITY(brop211->Rows(), 2);
  TEST_EQUALITY(brop211->Cols(), 2);
  TEST_EQUALITY(brop21100->Rows(), 3);
  TEST_EQUALITY(brop21100->Cols(), 3);
  TEST_EQUALITY(brop21100->getMatrix(0,0)->getRangeMap()->getGlobalNumElements(),comm->getSize() * 10);
  TEST_EQUALITY(brop21100->getMatrix(1,1)->getRangeMap()->getGlobalNumElements(),comm->getSize() * 40);
  TEST_EQUALITY(brop2->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop2->getDomainMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop200->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop200->getDomainMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop211->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop211->getDomainMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop21100->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop21100->getDomainMapExtractor()->getThyraMode(), true);
}


/// simple test routine for the apply function of BlockedCrsMatrix
TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedCrsMatrix, MatrixMatrixMult, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::Map<LO, GO, Node> MapClass;
  typedef Xpetra::StridedMap<LO, GO, Node> StridedMapClass;
  typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;
  typedef Xpetra::StridedMapFactory<LO, GO, Node> StridedMapFactoryClass;
  typedef Xpetra::Vector<Scalar, LO, GO, Node> VectorClass;
  typedef Xpetra::VectorFactory<Scalar, LO, GO, Node> VectorFactoryClass;
  typedef Xpetra::MapExtractor<Scalar,LO,GO,Node> MapExtractorClass;
  typedef Xpetra::MapExtractorFactory<Scalar,LO,GO,Node> MapExtractorFactoryClass;
  typedef Xpetra::Matrix<Scalar, LO, GO, Node> MatrixClass;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  M testMap(1,0,comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  const Teuchos::RCP<const MapClass> pointmap = MapFactoryClass::Build(lib, 12, 0, comm);

  // generate local maps for loading matrices
  Teuchos::Array<GO> velgidvec; // global strided maps
  Teuchos::Array<GO> pregidvec;
  Teuchos::Array<GO> fullgidvec; // full global map
  for (LO i=0; i<Teuchos::as<LO>(pointmap->getNodeNumElements()); i++)
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

  Teuchos::RCP<MatrixClass> A2 = Xpetra::IO<Scalar,LO,GO,Node>::Read("A.mat", fullmap->getMap());

  std::vector<Teuchos::RCP<const MapClass> > xmaps2;
  xmaps2.push_back(velmap);
  xmaps2.push_back(premap);

  Teuchos::RCP<const MapExtractorClass> map_extractor2 = MapExtractorFactoryClass::Build(fullmap,xmaps2);

  Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> > bOp2 =
      Xpetra::MatrixUtils<Scalar, LO, GO, Node>::SplitMatrix(*A2,map_extractor2,map_extractor2);

  // matrix-matrix multiplication of standard matrices
  Teuchos::RCP<Xpetra::Matrix<Scalar,LO,GO,Node> > fuAfuA_2 = Xpetra::MatrixMatrix<Scalar,LO,GO,Node>::Multiply(*A,false,*A2,false,out);
  fuAfuA_2->describe(out);

  // matrix-matrix multiplication of blocked operators
  Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> > bOpbOp_2 = Xpetra::MatrixMatrix<Scalar,LO,GO,Node>::TwoMatrixMultiplyBlock(*bOp,false,*bOp2,false,out);

  // build gloabl vector with one entries
  Teuchos::RCP<VectorClass> ones = VectorFactoryClass::Build(fullmap, true);
  Teuchos::RCP<VectorClass> exp  = VectorFactoryClass::Build(fullmap, true);
  Teuchos::RCP<VectorClass> res  = VectorFactoryClass::Build(fullmap, true);
  Teuchos::RCP<VectorClass> rnd  = VectorFactoryClass::Build(fullmap, true);
  ones->putScalar(STS::one());
  rnd->randomize();

  fuAfuA_2->apply(*ones, *exp);
  bOpbOp_2->apply(*ones, *res);
  res->update(-STS::one(),*exp,STS::one());
  TEUCHOS_TEST_COMPARE(res->norm2(), <, 1e-16, out, success);
  TEUCHOS_TEST_COMPARE(res->normInf(), <, 1e-16, out, success);

  A->apply(*rnd, *exp);
  bOp->apply(*rnd, *res);
  res->update(-STS::one(),*exp,STS::one());
  TEUCHOS_TEST_COMPARE(res->norm2(), <, 5e-14, out, success);
  TEUCHOS_TEST_COMPARE(res->normInf(), <, 5e-14, out, success);

  TEUCHOS_TEST_EQUALITY(fuAfuA_2->getGlobalNumEntries(),312,out,success);
  TEUCHOS_TEST_EQUALITY(bOpbOp_2->getGlobalNumEntries(),312,out,success);

  Teuchos::RCP<const MapClass> rgMap0 = bOpbOp_2->getRangeMap(0);
  Teuchos::RCP<const StridedMapClass> strRgMap0 = Teuchos::rcp_dynamic_cast<const StridedMapClass>(rgMap0);
  TEUCHOS_TEST_EQUALITY(strRgMap0==Teuchos::null, false, out, success );
  std::vector<size_t> strInfoData = strRgMap0->getStridingData();
  TEUCHOS_TEST_EQUALITY(strInfoData[0], 2, out, success );
  TEUCHOS_TEST_EQUALITY(strInfoData[1], 1, out, success );
  TEUCHOS_TEST_EQUALITY(strRgMap0->getFixedBlockSize(), 3, out, success );
  TEUCHOS_TEST_EQUALITY(strRgMap0->getStridedBlockId(), 0, out, success );

  Teuchos::RCP<const MapClass> rgMap = bOpbOp_2->getRangeMap();
  Teuchos::RCP<const StridedMapClass> strRgMap = Teuchos::rcp_dynamic_cast<const StridedMapClass>(rgMap);
  TEUCHOS_TEST_EQUALITY(strRgMap==Teuchos::null, false, out, success );
  strInfoData = strRgMap->getStridingData();
  TEUCHOS_TEST_EQUALITY(strInfoData[0], 2, out, success );
  TEUCHOS_TEST_EQUALITY(strInfoData[1], 1, out, success );
  TEUCHOS_TEST_EQUALITY(strRgMap->getFixedBlockSize(), 3, out, success );
  TEUCHOS_TEST_EQUALITY(strRgMap->getStridedBlockId(), -1, out, success );

  Teuchos::RCP<const MapClass> doMap0 = bOpbOp_2->getDomainMap(0);
  Teuchos::RCP<const StridedMapClass> strDoMap0 = Teuchos::rcp_dynamic_cast<const StridedMapClass>(doMap0);
  TEUCHOS_TEST_EQUALITY(strDoMap0==Teuchos::null, false, out, success );
  strInfoData = strDoMap0->getStridingData();
  TEUCHOS_TEST_EQUALITY(strInfoData[0], 2, out, success );
  TEUCHOS_TEST_EQUALITY(strInfoData[1], 1, out, success );
  TEUCHOS_TEST_EQUALITY(strDoMap0->getFixedBlockSize(), 3, out, success );
  TEUCHOS_TEST_EQUALITY(strDoMap0->getStridedBlockId(), 0, out, success );

  Teuchos::RCP<const MapClass> doMap = bOpbOp_2->getDomainMap();
  Teuchos::RCP<const StridedMapClass> strDoMap = Teuchos::rcp_dynamic_cast<const StridedMapClass>(doMap);
  TEUCHOS_TEST_EQUALITY(strDoMap==Teuchos::null, false, out, success );
  strInfoData = strDoMap->getStridingData();
  TEUCHOS_TEST_EQUALITY(strInfoData[0], 2, out, success );
  TEUCHOS_TEST_EQUALITY(strInfoData[1], 1, out, success );
  TEUCHOS_TEST_EQUALITY(strDoMap->getFixedBlockSize(), 3, out, success );
  TEUCHOS_TEST_EQUALITY(strDoMap->getStridedBlockId(), -1, out, success );

}

/// simple test for matrix-matrix multiplication for a 2x2 blocked matrix with a 2x1 blocked matrix
/*TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockedCrsMatrix, EpetraMatrixMatrixMult2x1, Scalar, LO, GO, Node)
{
#ifdef HAVE_XPETRA_EPETRAEXT
  RCP<const Comm<int> > comm = getDefaultComm();

  // build maps
  RCP<Epetra_Map> rowmap1     = Teuchos::rcp(new Epetra_Map(24,0,*Xpetra::toEpetra(comm)));
  RCP<Epetra_Map> rowmap2     = Teuchos::rcp(new Epetra_Map(12,24,*Xpetra::toEpetra(comm)));
  RCP<Epetra_Map> dommap1     = Teuchos::rcp(new Epetra_Map(8,0,*Xpetra::toEpetra(comm)));
  RCP<Epetra_Map> dommap2     = Teuchos::rcp(new Epetra_Map(4,8,*Xpetra::toEpetra(comm)));

  std::vector<RCP<const Epetra_Map> > rowmaps;
  rowmaps.push_back(rowmap1); rowmaps.push_back(rowmap2);
  std::vector<RCP<const Epetra_Map> > dommaps;
  dommaps.push_back(dommap1); dommaps.push_back(dommap2);

  RCP<Epetra_Map> fullrowmap = MergeMaps(rowmaps);
  RCP<Epetra_Map> fulldommap = MergeMaps(dommaps);

  // read in matrices in matrix market format
  Epetra_CrsMatrix* ptrA   = 0;
  Epetra_CrsMatrix* ptrP   = 0;
  EpetraExt::MatrixMarketFileToCrsMatrix("A.mat",*fullrowmap,*fullrowmap,*fullrowmap,ptrA);
  EpetraExt::MatrixMarketFileToCrsMatrix("P.mat",*fullrowmap,*fullrowmap,*fulldommap,ptrP);
  Teuchos::RCP<Epetra_CrsMatrix> epA = Teuchos::rcp(ptrA);
  Teuchos::RCP<Epetra_CrsMatrix> epP = Teuchos::rcp(ptrP);

  // split fullA into A11,..., A22
  Teuchos::RCP<Epetra_CrsMatrix> epA11;
  Teuchos::RCP<Epetra_CrsMatrix> epA12;
  Teuchos::RCP<Epetra_CrsMatrix> epA21;
  Teuchos::RCP<Epetra_CrsMatrix> epA22;

  SplitMatrix2x2(epA,*rowmap1,*rowmap2,epA11,epA12,epA21,epA22);

  Teuchos::RCP<Epetra_CrsMatrix> epP1;
  Teuchos::RCP<Epetra_CrsMatrix> epP2;

  SplitMatrix2x1(epP,*rowmap1,*rowmap2,*fulldommap,epP1,epP2);

  ////////////////// transform Epetra stuff to Xpetra

  // build Xpetra objects from Epetra_CrsMatrix objects
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO> > xA11 = Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<GO,Node>(epA11));
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO> > xA12 = Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<GO,Node>(epA12));
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO> > xA21 = Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<GO,Node>(epA21));
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO> > xA22 = Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<GO,Node>(epA22));

  // build Xpetra objects from Epetra_CrsMatrix objects
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO> > xP1 = Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<GO,Node>(epP1));
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO> > xP2 = Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<GO,Node>(epP2));

  Teuchos::RCP<Xpetra::EpetraMapT<GO, Node> > xfullrowmap = Teuchos::rcp(new Xpetra::EpetraMapT<GO, Node> (fullrowmap));
  Teuchos::RCP<Xpetra::EpetraMapT<GO, Node> > xfulldommap = Teuchos::rcp(new Xpetra::EpetraMapT<GO, Node> (fulldommap));
  Teuchos::RCP<Xpetra::EpetraMapT<GO, Node> > xrowmap1  = Teuchos::rcp(new Xpetra::EpetraMapT<GO, Node> (rowmap1 ));
  Teuchos::RCP<Xpetra::EpetraMapT<GO, Node> > xrowmap2  = Teuchos::rcp(new Xpetra::EpetraMapT<GO, Node> (rowmap2 ));
  Teuchos::RCP<Xpetra::EpetraMapT<GO, Node> > xdommap1  = Teuchos::rcp(new Xpetra::EpetraMapT<GO, Node> (dommap1 ));
  Teuchos::RCP<Xpetra::EpetraMapT<GO, Node> > xdommap2  = Teuchos::rcp(new Xpetra::EpetraMapT<GO, Node> (dommap2 ));

  // build map extractor objects
  std::vector<Teuchos::RCP<const Xpetra::Map<LO,GO> > > xrowmaps;
  xrowmaps.push_back(xrowmap1);
  xrowmaps.push_back(xrowmap2);

  Teuchos::RCP<const Xpetra::MapExtractor<Scalar,LO,GO> > map_extractor = Xpetra::MapExtractorFactory<Scalar,LO,GO>::Build(xfullrowmap,xrowmaps);

  std::vector<Teuchos::RCP<const Xpetra::Map<LO,GO> > > xdommaps;
  xdommaps.push_back(xfulldommap);

  Teuchos::RCP<const Xpetra::MapExtractor<Scalar,LO,GO> > map_domextractor = Xpetra::MapExtractorFactory<Scalar,LO,GO>::Build(xfulldommap,xdommaps);

  // build blocked operators

  // build 2x2 blocked operator
  Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar,LO,GO> > bA = Teuchos::rcp(new Xpetra::BlockedCrsMatrix<Scalar,LO,GO>(map_extractor,map_extractor,10));
  bA->setMatrix(0,0,xA11);
  bA->setMatrix(0,1,xA12);
  bA->setMatrix(1,0,xA21);
  bA->setMatrix(1,1,xA22);
  bA->fillComplete();

  // build 2x1 blocked operator
  Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar,LO,GO> > bP = Teuchos::rcp(new Xpetra::BlockedCrsMatrix<Scalar,LO,GO>(map_extractor,map_domextractor,10));
  bP->setMatrix(0,0,xP1);
  bP->setMatrix(1,0,xP2);
  bP->fillComplete();

  RCP<Xpetra::BlockedCrsMatrix<Scalar,LO,GO> > bAbP = MueLu::Utils<Scalar,LO,GO>::TwoMatrixMultiplyBlock(bA,false,bP,false);

  TEUCHOS_TEST_EQUALITY(bAbP->Rows(), 2, out, success );
  TEUCHOS_TEST_EQUALITY(bAbP->Cols(), 1, out, success );

  RCP<Xpetra::CrsMatrix<Scalar,LO,GO> > bAbPmerged = bAbP->Merge();

  // tests
  Teuchos::RCP<Xpetra::Vector<Scalar,LO,GO> > onevector =  Xpetra::VectorFactory<Scalar,LO,GO>::Build(bAbP->getDomainMap() ,true);
  Teuchos::RCP<Xpetra::Vector<Scalar,LO,GO> > resvector =  Xpetra::VectorFactory<Scalar,LO,GO>::Build(bAbP->getRangeMap() ,true);
  Teuchos::RCP<Xpetra::Vector<Scalar,LO,GO> > resvector2=  Xpetra::VectorFactory<Scalar,LO,GO>::Build(bAbPmerged->getRangeMap() ,true);
  onevector->putScalar(1.0);
  bAbP->apply(*onevector,*resvector);
  bAbPmerged->apply(*onevector,*resvector2);

  resvector2->update(1.0,*resvector,-1.0);
  TEUCHOS_TEST_COMPARE(resvector2->norm2(), <, 1e-16, out, success);
#endif
}*/

//
// INSTANTIATIONS
//
#ifdef HAVE_XPETRA_TPETRA

  #define XPETRA_TPETRA_TYPES( S, LO, GO, N) \
    typedef typename Xpetra::TpetraMap<LO,GO,N> M##LO##GO##N; \
    typedef typename Xpetra::TpetraCrsMatrix<S,LO,GO,N> MA##S##LO##GO##N;

#endif

#ifdef HAVE_XPETRA_EPETRA

  #define XPETRA_EPETRA_TYPES( S, LO, GO, N) \
    typedef typename Xpetra::EpetraMapT<GO,N> M##LO##GO##N; \
    typedef typename Xpetra::EpetraCrsMatrixT<GO,N> MA##S##LO##GO##N;

#endif


#define XP_MATRIX_INSTANT(S,LO,GO,N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedCrsMatrix, SplitMatrix, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedCrsMatrix, CreateBlockedDiagonalOp, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedCrsMatrix, CreateBlockedDiagonalOpThyra, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedCrsMatrix, ReorderBlockOperator, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedCrsMatrix, ReorderBlockOperator2, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedCrsMatrix, ReorderBlockOperatorThyra, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedCrsMatrix, ReorderBlockOperator2Thyra, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedCrsMatrix, ReorderBlockOperatorApply, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedCrsMatrix, ReorderBlockOperatorApply2, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedCrsMatrix, ReorderBlockOperatorApplyThyra, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedCrsMatrix, ReorderBlockOperatorApply2Thyra, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedCrsMatrix, Apply, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedCrsMatrix, getLocalDiagCopy, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedCrsMatrix, deepCopy, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedCrsMatrix, MatrixMatrixMult, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N )

// List of tests which run only with Tpetra
#define XP_TPETRA_MATRIX_INSTANT(S,LO,GO,N)

// List of tests which run only with Epetra
#define XP_EPETRA_MATRIX_INSTANT(S,LO,GO,N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedCrsMatrix, ReadWriteBlockedMatrix, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N )

// TODO reactivate these tests after moving MM multiplication code to xpetra...
//TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockedCrsMatrix, EpetraMatrixMatrixMult2x1, SC, LO, GO, Node )

#if defined(HAVE_XPETRA_TPETRA)

#include <TpetraCore_config.h>
#include <TpetraCore_ETIHelperMacros.h>

TPETRA_ETI_MANGLING_TYPEDEFS()
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( XPETRA_TPETRA_TYPES )
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( XP_TPETRA_MATRIX_INSTANT )
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( XP_MATRIX_INSTANT )

#endif

#if defined(HAVE_XPETRA_EPETRA)

#include "Xpetra_Map.hpp" // defines EpetraNode
typedef Xpetra::EpetraNode EpetraNode;
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
XPETRA_EPETRA_TYPES(double,int,int,EpetraNode)
XP_EPETRA_MATRIX_INSTANT(double,int,int,EpetraNode)
XP_MATRIX_INSTANT(double,int,int,EpetraNode)
#endif
// EpetraExt routines are not working with 64 bit
/*#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
typedef long long LongLong;
XPETRA_EPETRA_TYPES(double,int,LongLong,EpetraNode)
XP_EPETRA_MATRIX_INSTANT(double,int,LongLong,EpetraNode)
#endif*/

#endif

}
