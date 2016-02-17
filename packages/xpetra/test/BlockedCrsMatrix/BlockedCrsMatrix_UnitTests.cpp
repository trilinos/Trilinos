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
#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_Exceptions.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_MatrixUtils.hpp>
#include <Xpetra_IO.hpp>

//#include <MueLu_Utilities.hpp> //TODO: Xpetra tests should not use MueLu

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
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedCrsMatrix, Apply, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedCrsMatrix, MatrixMatrixMult, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N )

// List of tests which run only with Tpetra
#define XP_TPETRA_MATRIX_INSTANT(S,LO,GO,N)

// List of tests which run only with Epetra
#define XP_EPETRA_MATRIX_INSTANT(S,LO,GO,N)

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
