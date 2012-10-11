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
/*
 * MatrixMatrix_UnitTests.cpp
 *
 *  Created on: Oct 9, 2012
 *      Author: wiesner
 */

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <Xpetra_ConfigDefs.hpp>

#ifdef HAVE_XPETRA_EPETRAEXT
// EpetraExt
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#endif

#include <Xpetra_DefaultPlatform.hpp>
#include <Teuchos_as.hpp>

#include <Xpetra_Map.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#ifdef HAVE_XPETRA_TPETRA
#include <Xpetra_TpetraCrsMatrix.hpp>
#include <MatrixMarket_Tpetra.hpp>
#endif
#ifdef HAVE_XPETRA_EPETRA
#include <Xpetra_EpetraCrsMatrix.hpp>
#endif
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Exceptions.hpp>

#include <XpetraExt_MatrixMatrix.hpp>

//#include <MueLu_Utilities.hpp> //TODO: Xpetra tests should not use MueLu

namespace {

  using Teuchos::Array;
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::arcp;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Teuchos::Tuple;
  using Teuchos::tuple;
  using std::sort;
  using std::find;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::Comm;

  using Xpetra::DefaultPlatform;
  using Xpetra::Matrix;
  using Xpetra::CrsMatrix;
#ifdef HAVE_XPETRA_TPETRA
  using Xpetra::TpetraCrsMatrix; //TMP
#endif
  using Xpetra::Map;

  using Xpetra::viewLabel_t;

  bool testMpi = true;
  double errorTolSlack = 1e+1;



  RCP<const Comm<int> > getDefaultComm()
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
    clp.setOption(
                  "error-tol-slack", &errorTolSlack,
                  "Slack off of machine epsilon used to check test results" );
  }

  //
  // UNIT TESTS
  //

  /// unit test for matrix-matrix multiplication (both for Epetra and Tpetra)
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MatrixMatrix, Multiply, Scalar, LO, GO, Node )
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    typedef Matrix<Scalar, LO, GO, Node> Matrix;
    typedef CrsMatrix<Scalar, LO, GO, Node> CrsMatrix;
    typedef Xpetra::Map<LO, GO, Node> MapClass;
    typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;

    //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
    //yAB->describe(*fos, Teuchos::VERB_EXTREME);

#ifdef HAVE_XPETRA_EPETRA
    { // Epetra test
      // get a comm and node
      RCP<const Comm<int> > comm = getDefaultComm();

      Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;

      // generate problem
      LO nEle = 6;
      const RCP<const MapClass> map = MapFactoryClass::Build(lib, nEle, 0, comm);
      const RCP<const Xpetra::EpetraMap> XepMap = Teuchos::rcp_dynamic_cast<const Xpetra::EpetraMap>(map);

      /////////////////////////////////////// transform Xpetra::Map objects to Epetra
      // this is needed for AztecOO
      const Teuchos::RCP<const Epetra_Map> epMap = Teuchos::rcpFromRef(XepMap->getEpetra_Map());
      /////////////////////////////////////// import problem matrix and RHS from files (-> Epetra)

      // read in problem
      Epetra_CrsMatrix * ptrA  = 0;
      Epetra_CrsMatrix * ptrB  = 0;
      Epetra_CrsMatrix * ptrAB = 0;
      Epetra_CrsMatrix * ptrAtB = 0;
      Epetra_CrsMatrix * ptrABt = 0;
      Epetra_CrsMatrix * ptrAtBt = 0;
      EpetraExt::MatrixMarketFileToCrsMatrix("A.mat",*epMap,*epMap,*epMap,ptrA);
      EpetraExt::MatrixMarketFileToCrsMatrix("B.mat",*epMap,*epMap,*epMap,ptrB);
      EpetraExt::MatrixMarketFileToCrsMatrix("AB.mat",*epMap,*epMap,*epMap,ptrAB);
      EpetraExt::MatrixMarketFileToCrsMatrix("AtB.mat",*epMap,*epMap,*epMap,ptrAtB);
      EpetraExt::MatrixMarketFileToCrsMatrix("ABt.mat",*epMap,*epMap,*epMap,ptrABt);
      EpetraExt::MatrixMarketFileToCrsMatrix("AtBt.mat",*epMap,*epMap,*epMap,ptrAtBt);
      if(ptrA == NULL || ptrB == NULL || ptrAB == NULL ||
         ptrAtB == NULL || ptrABt == NULL || ptrAtBt == NULL)
        std::cout << "matrix not found" << std::endl;
      RCP<Epetra_CrsMatrix> epA = Teuchos::rcp(ptrA);
      RCP<Epetra_CrsMatrix> epB = Teuchos::rcp(ptrB);
      RCP<Epetra_CrsMatrix> epAB = Teuchos::rcp(ptrAB);
      RCP<Epetra_CrsMatrix> epAtB = Teuchos::rcp(ptrAtB);
      RCP<Epetra_CrsMatrix> epABt = Teuchos::rcp(ptrABt);
      RCP<Epetra_CrsMatrix> epAtBt = Teuchos::rcp(ptrAtBt);

      /////////////////////////////////////// transform Epetra objects to Xpetra (needed for MueLu)

      // build Xpetra objects from Epetra_CrsMatrix objects
      Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xAmat = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(epA));
      Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xBmat = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(epB));
      Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xABmat = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(epAB));
      Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xAtBmat = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(epAtB));
      Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xABtmat = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(epABt));
      Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xAtBtmat = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(epAtBt));

      Teuchos::RCP<Xpetra::Matrix<Scalar,LO,GO,Node> > xA = Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar,LO,GO,Node>(xAmat));
      Teuchos::RCP<Xpetra::Matrix<Scalar,LO,GO,Node> > xB = Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar,LO,GO,Node>(xBmat));
      Teuchos::RCP<Xpetra::Matrix<Scalar,LO,GO,Node> > xAB= Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar,LO,GO,Node>(xABmat));
      Teuchos::RCP<Xpetra::Matrix<Scalar,LO,GO,Node> > xAtB= Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar,LO,GO,Node>(xAtBmat));
      Teuchos::RCP<Xpetra::Matrix<Scalar,LO,GO,Node> > xABt= Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar,LO,GO,Node>(xABtmat));
      Teuchos::RCP<Xpetra::Matrix<Scalar,LO,GO,Node> > xAtBt= Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar,LO,GO,Node>(xAtBtmat));

      Teuchos::RCP<Xpetra::Matrix<Scalar,LO,GO,Node> > yAB= Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar,LO,GO,Node>(map, 6));

      //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
      //yAB->describe(*fos,Teuchos::VERB_EXTREME);


      Xpetra::MatrixMatrix::Multiply (
        *xA,
        false,
        *xB,
        false,
        *yAB);
      Xpetra::MatrixMatrix::Add (*xAB,false,1.0,*yAB,-1.0);
      TEUCHOS_TEST_EQUALITY(yAB->getFrobeniusNorm(), 0, out, success );

      Xpetra::MatrixMatrix::Multiply (
        *xA,
        true,
        *xB,
        false,
        *yAB);
      Xpetra::MatrixMatrix::Add (*xAtB,false,1.0,*yAB,-1.0);
      TEUCHOS_TEST_EQUALITY(yAB->getFrobeniusNorm(), 0, out, success );

      Xpetra::MatrixMatrix::Multiply (
        *xA,
        false,
        *xB,
        true,
        *yAB);
      Xpetra::MatrixMatrix::Add (*xABt,false,1.0,*yAB,-1.0);
      TEUCHOS_TEST_EQUALITY(yAB->getFrobeniusNorm(), 0, out, success );

      Xpetra::MatrixMatrix::Multiply (
        *xA,
        true,
        *xB,
        true,
        *yAB);
      Xpetra::MatrixMatrix::Add (*xAtBt,false,1.0,*yAB,-1.0);
      TEUCHOS_TEST_EQUALITY(yAB->getFrobeniusNorm(), 0, out, success );
    } // end Epetra test
#endif


#ifdef HAVE_XPETRA_TPETRA
    { // Tpetra test

      // get a comm and node
      RCP<const Comm<int> > comm = getDefaultComm();

      Xpetra::UnderlyingLib lib = Xpetra::UseTpetra;

      // define node
      Teuchos::ParameterList defaultParams;
      Teuchos::RCP<Kokkos::SerialNode> pNode = Teuchos::rcp (new Kokkos::SerialNode (defaultParams));

      // define map
      LO nEle = 6;
      const RCP<const MapClass> map = MapFactoryClass::Build(lib, nEle, 0, comm);

      // read in matrices
      typedef Tpetra::MatrixMarket::Reader<Tpetra::CrsMatrix<Scalar, LO, GO, Node> > reader_type;

      Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO,Node> > tpA = reader_type::readSparseFile("A.mat",comm,pNode );
      Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO,Node> > tpB = reader_type::readSparseFile("B.mat",comm,pNode );
      Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO,Node> > tpAB = reader_type::readSparseFile("AB.mat",comm,pNode );
      Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO,Node> > tpAtB = reader_type::readSparseFile("AtB.mat",comm,pNode );
      Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO,Node> > tpABt = reader_type::readSparseFile("ABt.mat",comm,pNode );
      Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO,Node> > tpAtBt = reader_type::readSparseFile("AtBt.mat",comm,pNode );

      // transform to Xpetra
      Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xAmat = Teuchos::rcp(new Xpetra::TpetraCrsMatrix<Scalar,LO,GO,Node>(tpA));
      Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xBmat = Teuchos::rcp(new Xpetra::TpetraCrsMatrix<Scalar,LO,GO,Node>(tpB));
      Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xABmat = Teuchos::rcp(new Xpetra::TpetraCrsMatrix<Scalar,LO,GO,Node>(tpAB));
      Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xAtBmat = Teuchos::rcp(new Xpetra::TpetraCrsMatrix<Scalar,LO,GO,Node>(tpAtB));
      Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xABtmat = Teuchos::rcp(new Xpetra::TpetraCrsMatrix<Scalar,LO,GO,Node>(tpABt));
      Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node> > xAtBtmat = Teuchos::rcp(new Xpetra::TpetraCrsMatrix<Scalar,LO,GO,Node>(tpAtBt));

      Teuchos::RCP<Xpetra::Matrix<Scalar,LO,GO,Node> > xA = Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar,LO,GO,Node>(xAmat));
      Teuchos::RCP<Xpetra::Matrix<Scalar,LO,GO,Node> > xB = Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar,LO,GO,Node>(xBmat));
      Teuchos::RCP<Xpetra::Matrix<Scalar,LO,GO,Node> > xAB= Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar,LO,GO,Node>(xABmat));
      Teuchos::RCP<Xpetra::Matrix<Scalar,LO,GO,Node> > xAtB= Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar,LO,GO,Node>(xAtBmat));
      Teuchos::RCP<Xpetra::Matrix<Scalar,LO,GO,Node> > xABt= Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar,LO,GO,Node>(xABtmat));
      Teuchos::RCP<Xpetra::Matrix<Scalar,LO,GO,Node> > xAtBt= Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar,LO,GO,Node>(xAtBtmat));

      Teuchos::RCP<Xpetra::Matrix<Scalar,LO,GO,Node> > yAB= Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar,LO,GO,Node>(map, 10));
      Teuchos::RCP<Xpetra::Matrix<Scalar,LO,GO,Node> > yAtB= Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar,LO,GO,Node>(map, 10));
      Teuchos::RCP<Xpetra::Matrix<Scalar,LO,GO,Node> > yABt= Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar,LO,GO,Node>(map, 10));
      Teuchos::RCP<Xpetra::Matrix<Scalar,LO,GO,Node> > yAtBt= Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar,LO,GO,Node>(map, 10));

      Xpetra::MatrixMatrix::Multiply (
        *xA,
        false,
        *xB,
        false,
        *yAB);
      TEUCHOS_TEST_EQUALITY(xAB->getFrobeniusNorm(), yAB->getFrobeniusNorm(), out, success );
      TEUCHOS_TEST_EQUALITY(xAB->getGlobalNumDiags(), yAB->getGlobalNumDiags(), out, success );
      TEUCHOS_TEST_EQUALITY(xAB->getNodeNumEntries(), yAB->getNodeNumEntries(), out, success );

      Xpetra::MatrixMatrix::Multiply (
          *xA,
          true,
          *xB,
          false,
          *yAtB);
        TEUCHOS_TEST_EQUALITY(xAtB->getFrobeniusNorm(), yAtB->getFrobeniusNorm(), out, success );
        TEUCHOS_TEST_EQUALITY(xAtB->getGlobalNumDiags(), yAtB->getGlobalNumDiags(), out, success );
        TEUCHOS_TEST_EQUALITY(xAtB->getNodeNumEntries(), yAtB->getNodeNumEntries(), out, success );

      Xpetra::MatrixMatrix::Multiply (
          *xA,
          false,
          *xB,
          true,
          *yABt);
        TEUCHOS_TEST_EQUALITY(xABt->getFrobeniusNorm(), yABt->getFrobeniusNorm(), out, success );
        TEUCHOS_TEST_EQUALITY(xABt->getGlobalNumDiags(), yABt->getGlobalNumDiags(), out, success );
        TEUCHOS_TEST_EQUALITY(xABt->getNodeNumEntries(), yABt->getNodeNumEntries(), out, success );

      Xpetra::MatrixMatrix::Multiply (
          *xA,
          true,
          *xB,
          true,
          *yAtBt);
        TEUCHOS_TEST_EQUALITY(xAtBt->getFrobeniusNorm(), yAtBt->getFrobeniusNorm(), out, success );
        TEUCHOS_TEST_EQUALITY(xAtBt->getGlobalNumDiags(), yAtBt->getGlobalNumDiags(), out, success );
        TEUCHOS_TEST_EQUALITY(xAtBt->getNodeNumEntries(), yAtBt->getNodeNumEntries(), out, success );
    }// end Tpetra test
#endif
  }

  //
  // INSTANTIATIONS
  //
#define UNIT_TEST_GROUP_ORDINAL( SC, LO, GO, Node )                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MatrixMatrix, Multiply, SC, LO, GO, Node )

  typedef Kokkos::DefaultNode::DefaultNodeType DefaultNodeType;

  UNIT_TEST_GROUP_ORDINAL(double, int, int, DefaultNodeType)

} // end namespace

