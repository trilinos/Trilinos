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
#include <Xpetra_TpetraCrsMatrix.hpp>
#include <Xpetra_TpetraBlockCrsMatrix.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_BlockCrsMatrix.hpp>
#include <Tpetra_BlockCrsMatrix_Helpers.hpp>
#ifdef HAVE_XPETRA_EPETRA
#include <Xpetra_EpetraCrsMatrix.hpp>
#endif
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Exceptions.hpp>

//#include <XpetraExt_MatrixMatrix.hpp>
#include <Xpetra_MatrixMatrix.hpp>

//#include <MueLu_Utilities.hpp> //TODO: Xpetra tests should not use MueLu

#ifdef XPETRA_TEST_USE_LONGLONG_GO
#define MatrixMarketFileToCrsMatrix MatrixMarketFileToCrsMatrix64
#endif

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
  TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( MatrixMatrix, Multiply_Epetra, M, MA, Scalar, LO, GO, Node )
  {
#if defined(HAVE_XPETRA_EPETRA) && defined(HAVE_XPETRA_EPETRAEXT)
    typedef Xpetra::Map<LO, GO, Node> MapClass;
    typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;
    typedef Xpetra::CrsMatrix<Scalar,LO,GO,Node> CrsMatrixClass;
    typedef Xpetra::Matrix<Scalar,LO,GO,Node> MatrixClass;
    typedef Xpetra::CrsMatrixWrap<Scalar,LO,GO,Node> CrsMatrixWrapClass;

    //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
    //yAB->describe(*fos, Teuchos::VERB_EXTREME);

    { // Epetra test
      // get a comm and node
      RCP<const Comm<int> > comm = getDefaultComm();

      Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;

      // generate problem
      LO nEle = 6;
      const RCP<const MapClass> map = MapFactoryClass::Build(lib, nEle, 0, comm);
      // TODO get rid of this...
      //#ifndef XPETRA_TEST_USE_LONGLONG_GO
      const RCP<const Xpetra::EpetraMapT<GO, Node> > XepMap = Teuchos::rcp_dynamic_cast<const Xpetra::EpetraMapT<GO, Node> >(map);
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
      TEUCHOS_TEST_FOR_EXCEPTION(ptrA == NULL || ptrB == NULL || ptrAB == NULL || ptrAtB == NULL || ptrABt == NULL || ptrAtBt == NULL,std::logic_error,"Could not open one or more matrix files");
      RCP<Epetra_CrsMatrix> epA = Teuchos::rcp(ptrA);
      RCP<Epetra_CrsMatrix> epB = Teuchos::rcp(ptrB);
      RCP<Epetra_CrsMatrix> epAB = Teuchos::rcp(ptrAB);
      RCP<Epetra_CrsMatrix> epAtB = Teuchos::rcp(ptrAtB);
      RCP<Epetra_CrsMatrix> epABt = Teuchos::rcp(ptrABt);
      RCP<Epetra_CrsMatrix> epAtBt = Teuchos::rcp(ptrAtBt);

      /////////////////////////////////////// transform Epetra objects to Xpetra (needed for MueLu)

      // build Xpetra objects from Epetra_CrsMatrix objects
      Teuchos::RCP<CrsMatrixClass> xAmat = Teuchos::rcp(new MA(epA));
      Teuchos::RCP<CrsMatrixClass> xBmat = Teuchos::rcp(new MA(epB));
      Teuchos::RCP<CrsMatrixClass> xABmat = Teuchos::rcp(new MA(epAB));
      Teuchos::RCP<CrsMatrixClass> xAtBmat = Teuchos::rcp(new MA(epAtB));
      Teuchos::RCP<CrsMatrixClass> xABtmat = Teuchos::rcp(new MA(epABt));
      Teuchos::RCP<CrsMatrixClass> xAtBtmat = Teuchos::rcp(new MA(epAtBt));

      Teuchos::RCP<MatrixClass> xA = Teuchos::rcp(new CrsMatrixWrapClass(xAmat));
      Teuchos::RCP<MatrixClass> xB = Teuchos::rcp(new CrsMatrixWrapClass(xBmat));
      Teuchos::RCP<MatrixClass> xAB= Teuchos::rcp(new CrsMatrixWrapClass(xABmat));
      Teuchos::RCP<MatrixClass> xAtB= Teuchos::rcp(new CrsMatrixWrapClass(xAtBmat));
      Teuchos::RCP<MatrixClass> xABt= Teuchos::rcp(new CrsMatrixWrapClass(xABtmat));
      Teuchos::RCP<MatrixClass> xAtBt= Teuchos::rcp(new CrsMatrixWrapClass(xAtBtmat));

      Teuchos::RCP<MatrixClass> yAB= Teuchos::rcp(new CrsMatrixWrapClass(map, 6));

      //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
      //yAB->describe(*fos,Teuchos::VERB_EXTREME);


      //Xpetra::MatrixMatrix::Multiply<Scalar, LO, GO, Node> (
      Xpetra::MatrixMatrix<Scalar, LO, GO, Node>::Multiply (
        *xA,
        false,
        *xB,
        false,
        *yAB);
      //Xpetra::MatrixMatrix::Add<Scalar, LO, GO, Node>(*xAB,false,1.0,*yAB,-1.0);
      TEUCHOS_TEST_EQUALITY(yAB->getFrobeniusNorm(), xAB->getFrobeniusNorm(), out, success );

      Xpetra::MatrixMatrix<Scalar, LO, GO, Node>::Multiply (
        *xA,
        true,
        *xB,
        false,
        *yAB);
      //Xpetra::MatrixMatrix::Add<Scalar, LO, GO, Node>(*xAtB,false,1.0,*yAB,-1.0);
      TEUCHOS_TEST_EQUALITY(yAB->getFrobeniusNorm(), xAtB->getFrobeniusNorm(), out, success );

      Xpetra::MatrixMatrix<Scalar, LO, GO, Node>::Multiply (
        *xA,
        false,
        *xB,
        true,
        *yAB);
      //Xpetra::MatrixMatrix::Add<Scalar, LO, GO, Node> (*xABt,false,1.0,*yAB,-1.0);
      TEUCHOS_TEST_EQUALITY(yAB->getFrobeniusNorm(), xABt->getFrobeniusNorm(), out, success );

      Xpetra::MatrixMatrix<Scalar, LO, GO, Node>::Multiply (
        *xA,
        true,
        *xB,
        true,
        *yAB);
      //Xpetra::MatrixMatrix::Add<Scalar, LO, GO, Node> (*xAtBt,false,1.0,*yAB,-1.0);
      TEUCHOS_TEST_EQUALITY(yAB->getFrobeniusNorm(), xAtBt->getFrobeniusNorm(), out, success );
    } // end Epetra test
#endif
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( MatrixMatrix, Multiply_Epetra64, M, MA, Scalar, LO, GO, Node )
  {
#if defined(HAVE_XPETRA_EPETRA) && defined(HAVE_XPETRA_EPETRAEXT)
    typedef Xpetra::Map<LO, GO, Node> MapClass;
    typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;
    typedef Xpetra::CrsMatrix<Scalar,LO,GO,Node> CrsMatrixClass;
    typedef Xpetra::Matrix<Scalar,LO,GO,Node> MatrixClass;
    typedef Xpetra::CrsMatrixWrap<Scalar,LO,GO,Node> CrsMatrixWrapClass;

    //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
    //yAB->describe(*fos, Teuchos::VERB_EXTREME);

    { // Epetra test
      // get a comm and node
      RCP<const Comm<int> > comm = getDefaultComm();

      Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;

      // generate problem
      LO nEle = 6;
      const RCP<const MapClass> map = MapFactoryClass::Build(lib, nEle, 0, comm);
      const RCP<const Xpetra::EpetraMapT<GO, Node> > XepMap = Teuchos::rcp_dynamic_cast<const Xpetra::EpetraMapT<GO, Node> >(map);
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
      EpetraExt::MatrixMarketFileToCrsMatrix64("A.mat",*epMap,*epMap,*epMap,ptrA);
      EpetraExt::MatrixMarketFileToCrsMatrix64("B.mat",*epMap,*epMap,*epMap,ptrB);
      EpetraExt::MatrixMarketFileToCrsMatrix64("AB.mat",*epMap,*epMap,*epMap,ptrAB);
      EpetraExt::MatrixMarketFileToCrsMatrix64("AtB.mat",*epMap,*epMap,*epMap,ptrAtB);
      EpetraExt::MatrixMarketFileToCrsMatrix64("ABt.mat",*epMap,*epMap,*epMap,ptrABt);
      EpetraExt::MatrixMarketFileToCrsMatrix64("AtBt.mat",*epMap,*epMap,*epMap,ptrAtBt);
      TEUCHOS_TEST_FOR_EXCEPTION(ptrA == NULL || ptrB == NULL || ptrAB == NULL || ptrAtB == NULL || ptrABt == NULL || ptrAtBt == NULL,std::logic_error,"Could not open one or more matrix files");
      RCP<Epetra_CrsMatrix> epA = Teuchos::rcp(ptrA);
      RCP<Epetra_CrsMatrix> epB = Teuchos::rcp(ptrB);
      RCP<Epetra_CrsMatrix> epAB = Teuchos::rcp(ptrAB);
      RCP<Epetra_CrsMatrix> epAtB = Teuchos::rcp(ptrAtB);
      RCP<Epetra_CrsMatrix> epABt = Teuchos::rcp(ptrABt);
      RCP<Epetra_CrsMatrix> epAtBt = Teuchos::rcp(ptrAtBt);

      /////////////////////////////////////// transform Epetra objects to Xpetra (needed for MueLu)

      // build Xpetra objects from Epetra_CrsMatrix objects
      Teuchos::RCP<CrsMatrixClass> xAmat = Teuchos::rcp(new MA(epA));
      Teuchos::RCP<CrsMatrixClass> xBmat = Teuchos::rcp(new MA(epB));
      Teuchos::RCP<CrsMatrixClass> xABmat = Teuchos::rcp(new MA(epAB));
      Teuchos::RCP<CrsMatrixClass> xAtBmat = Teuchos::rcp(new MA(epAtB));
      Teuchos::RCP<CrsMatrixClass> xABtmat = Teuchos::rcp(new MA(epABt));
      Teuchos::RCP<CrsMatrixClass> xAtBtmat = Teuchos::rcp(new MA(epAtBt));

      Teuchos::RCP<MatrixClass> xA = Teuchos::rcp(new CrsMatrixWrapClass(xAmat));
      Teuchos::RCP<MatrixClass> xB = Teuchos::rcp(new CrsMatrixWrapClass(xBmat));
      Teuchos::RCP<MatrixClass> xAB= Teuchos::rcp(new CrsMatrixWrapClass(xABmat));
      Teuchos::RCP<MatrixClass> xAtB= Teuchos::rcp(new CrsMatrixWrapClass(xAtBmat));
      Teuchos::RCP<MatrixClass> xABt= Teuchos::rcp(new CrsMatrixWrapClass(xABtmat));
      Teuchos::RCP<MatrixClass> xAtBt= Teuchos::rcp(new CrsMatrixWrapClass(xAtBtmat));

      Teuchos::RCP<MatrixClass> yAB= Teuchos::rcp(new CrsMatrixWrapClass(map, 6));

      //Xpetra::MatrixMatrix::Multiply<Scalar, LO, GO, Node> (
      Xpetra::MatrixMatrix<Scalar, LO, GO, Node>::Multiply (
        *xA,
        false,
        *xB,
        false,
        *yAB);
      //Xpetra::MatrixMatrix::Add<Scalar, LO, GO, Node>(*xAB,false,1.0,*yAB,-1.0);
      TEUCHOS_TEST_EQUALITY(yAB->getFrobeniusNorm(), xAB->getFrobeniusNorm(), out, success );

      Xpetra::MatrixMatrix<Scalar, LO, GO, Node>::Multiply (
        *xA,
        true,
        *xB,
        false,
        *yAB);
      //Xpetra::MatrixMatrix::Add<Scalar, LO, GO, Node>(*xAtB,false,1.0,*yAB,-1.0);
      TEUCHOS_TEST_EQUALITY(yAB->getFrobeniusNorm(), xAtB->getFrobeniusNorm(), out, success );

      Xpetra::MatrixMatrix<Scalar, LO, GO, Node>::Multiply (
        *xA,
        false,
        *xB,
        true,
        *yAB);
      //Xpetra::MatrixMatrix::Add<Scalar, LO, GO, Node> (*xABt,false,1.0,*yAB,-1.0);
      TEUCHOS_TEST_EQUALITY(yAB->getFrobeniusNorm(), xABt->getFrobeniusNorm(), out, success );

      Xpetra::MatrixMatrix<Scalar, LO, GO, Node>::Multiply (
        *xA,
        true,
        *xB,
        true,
        *yAB);
      //Xpetra::MatrixMatrix::Add<Scalar, LO, GO, Node> (*xAtBt,false,1.0,*yAB,-1.0);
      TEUCHOS_TEST_EQUALITY(yAB->getFrobeniusNorm(), xAtBt->getFrobeniusNorm(), out, success );
    } // end Epetra test
#endif
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( MatrixMatrix, Multiply_Tpetra, M, MA, Scalar, LO, GO, Node )
  {
#if defined(HAVE_TPETRA_INST_COMPLEX_FLOAT) || defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE)
    // The matrix reader does not work with complex scalars. Skip all tests then.
    return;
#endif
    typedef Xpetra::Map<LO, GO, Node> MapClass;
    typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;
    typedef Xpetra::CrsMatrix<Scalar,LO,GO,Node> CrsMatrixClass;
    typedef Xpetra::Matrix<Scalar,LO,GO,Node> MatrixClass;
    typedef Xpetra::CrsMatrixWrap<Scalar,LO,GO,Node> CrsMatrixWrapClass;

    //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
    //yAB->describe(*fos, Teuchos::VERB_EXTREME);

    { // Tpetra test

      // get a comm and node
      RCP<const Comm<int> > comm = getDefaultComm();
      Xpetra::UnderlyingLib lib = Xpetra::UseTpetra;

      // define map
      LO nEle = 6;
      const RCP<const MapClass> map = MapFactoryClass::Build(lib, nEle, 0, comm);

      // read in matrices
      typedef Tpetra::MatrixMarket::Reader<Tpetra::CrsMatrix<Scalar, LO, GO, Node> > reader_type;

      Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO,Node> > tpA = reader_type::readSparseFile("A.mat",comm);
      Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO,Node> > tpB = reader_type::readSparseFile("B.mat",comm);
      Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO,Node> > tpAB = reader_type::readSparseFile("AB.mat",comm);
      Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO,Node> > tpAtB = reader_type::readSparseFile("AtB.mat",comm);
      Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO,Node> > tpABt = reader_type::readSparseFile("ABt.mat",comm);
      Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO,Node> > tpAtBt = reader_type::readSparseFile("AtBt.mat",comm);

      // transform to Xpetra
      Teuchos::RCP<CrsMatrixClass> xAmat = Teuchos::rcp(new MA(tpA));
      Teuchos::RCP<CrsMatrixClass> xBmat = Teuchos::rcp(new MA(tpB));
      Teuchos::RCP<CrsMatrixClass> xABmat = Teuchos::rcp(new MA(tpAB));
      Teuchos::RCP<CrsMatrixClass> xAtBmat = Teuchos::rcp(new MA(tpAtB));
      Teuchos::RCP<CrsMatrixClass> xABtmat = Teuchos::rcp(new MA(tpABt));
      Teuchos::RCP<CrsMatrixClass> xAtBtmat = Teuchos::rcp(new MA(tpAtBt));

      Teuchos::RCP<MatrixClass> xA = Teuchos::rcp(new CrsMatrixWrapClass(xAmat));
      Teuchos::RCP<MatrixClass> xB = Teuchos::rcp(new CrsMatrixWrapClass(xBmat));
      Teuchos::RCP<MatrixClass> xAB= Teuchos::rcp(new CrsMatrixWrapClass(xABmat));
      Teuchos::RCP<MatrixClass> xAtB= Teuchos::rcp(new CrsMatrixWrapClass(xAtBmat));
      Teuchos::RCP<MatrixClass> xABt= Teuchos::rcp(new CrsMatrixWrapClass(xABtmat));
      Teuchos::RCP<MatrixClass> xAtBt= Teuchos::rcp(new CrsMatrixWrapClass(xAtBtmat));

      Teuchos::RCP<MatrixClass> yAB= Teuchos::rcp(new CrsMatrixWrapClass(map, 10));
      Teuchos::RCP<MatrixClass> yAtB= Teuchos::rcp(new CrsMatrixWrapClass(map, 10));
      Teuchos::RCP<MatrixClass> yABt= Teuchos::rcp(new CrsMatrixWrapClass(map, 10));
      Teuchos::RCP<MatrixClass> yAtBt= Teuchos::rcp(new CrsMatrixWrapClass(map, 10));

      Xpetra::MatrixMatrix<Scalar, LO, GO, Node>::Multiply (
        *xA,
        false,
        *xB,
        false,
        *yAB);
      TEUCHOS_TEST_EQUALITY(xAB->getFrobeniusNorm(), yAB->getFrobeniusNorm(), out, success );
      TEUCHOS_TEST_EQUALITY(xAB->getLocalNumEntries(), yAB->getLocalNumEntries(), out, success );

      Xpetra::MatrixMatrix<Scalar, LO, GO, Node>::Multiply (
          *xA,
          true,
          *xB,
          false,
          *yAtB);
        TEUCHOS_TEST_EQUALITY(xAtB->getFrobeniusNorm(), yAtB->getFrobeniusNorm(), out, success );
        TEUCHOS_TEST_EQUALITY(xAtB->getLocalNumEntries(), yAtB->getLocalNumEntries(), out, success );

      Xpetra::MatrixMatrix<Scalar, LO, GO, Node>::Multiply (
          *xA,
          false,
          *xB,
          true,
          *yABt);
        TEUCHOS_TEST_EQUALITY(xABt->getFrobeniusNorm(), yABt->getFrobeniusNorm(), out, success );
        TEUCHOS_TEST_EQUALITY(xABt->getLocalNumEntries(), yABt->getLocalNumEntries(), out, success );

      Xpetra::MatrixMatrix<Scalar, LO, GO, Node>::Multiply (
          *xA,
          true,
          *xB,
          true,
          *yAtBt);
        TEUCHOS_TEST_EQUALITY(xAtBt->getFrobeniusNorm(), yAtBt->getFrobeniusNorm(), out, success );
        TEUCHOS_TEST_EQUALITY(xAtBt->getLocalNumEntries(), yAtBt->getLocalNumEntries(), out, success );
    }// end Tpetra test
  }


  TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( MatrixMatrix, BlockCrs, M, MB, Scalar, LO, GO, Node )
  {
#if defined(HAVE_TPETRA_INST_COMPLEX_FLOAT) || defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE)
    // The matrix reader does not work with complex scalars. Skip all tests then.
    return;
#endif
    typedef Tpetra::BlockCrsMatrix<Scalar, LO, GO, Node> BCM;
    typedef Tpetra::CrsGraph<LO, GO, Node> graph_type;
    typedef Tpetra::Map<LO, GO, Node> map_type;

    typedef Xpetra::CrsMatrix<Scalar,LO,GO,Node> CrsMatrixClass;
    typedef Xpetra::TpetraBlockCrsMatrix<Scalar,LO,GO,Node> BlockCrsMatrixClass;
    typedef Xpetra::Matrix<Scalar,LO,GO,Node> MatrixClass;
    typedef Xpetra::CrsMatrixWrap<Scalar,LO,GO,Node> CrsMatrixWrapClass;
    using helpers = Xpetra::Helpers<Scalar,LO,GO,Node>;

    RCP<const Comm<int> > comm = getDefaultComm ();
    const GO INVALID = Teuchos::OrdinalTraits<GO>::invalid ();

    const size_t numLocalMeshPoints = 12;
    const GO indexBase = 1;
    // mfh 16 May 2014: Tpetra::CrsGraph still needs the row Map as an
    // RCP.  Later interface changes will let us pass in the Map by
    // const reference and assume view semantics.
    RCP<const map_type> meshRowMapPtr =
      rcp (new map_type (INVALID, numLocalMeshPoints, indexBase, comm));

    // Test w/ an empty graph
    const LO blockSize = 4;
    graph_type graph (meshRowMapPtr, 0);
    graph.fillComplete ();

    // Make the matrix (Tpetra)
    RCP<BCM> blockMat = rcp(new BCM (graph, blockSize));
    RCP<CrsMatrixClass> bmc = rcp(new BlockCrsMatrixClass(blockMat));
    RCP<CrsMatrixWrapClass> wrap = rcp(new CrsMatrixWrapClass(bmc));
    RCP<MatrixClass>  mat = wrap;

    // Now for the checks
    TEUCHOS_TEST_EQUALITY(helpers::isTpetraBlockCrs(mat), true, out, success);
    TEUCHOS_TEST_EQUALITY(helpers::isTpetraCrs(mat), false, out, success);

  }


  //
  // INSTANTIATIONS
  //

  #define XPETRA_TPETRA_TYPES( S, LO, GO, N) \
    typedef typename Xpetra::TpetraMap<LO,GO,N> M##LO##GO##N; \
    typedef typename Xpetra::TpetraCrsMatrix<S,LO,GO,N> MA##S##LO##GO##N; \
    typedef typename Xpetra::TpetraBlockCrsMatrix<S,LO,GO,N> MB##S##LO##GO##N; \

#ifdef HAVE_XPETRA_EPETRA

  #define XPETRA_EPETRA_TYPES( S, LO, GO, N) \
    typedef typename Xpetra::EpetraMapT<GO,N> M##LO##GO##N; \
    typedef typename Xpetra::EpetraCrsMatrixT<GO,N> MA##S##LO##GO##N;

#endif

  // List of tests which run only with Tpetra
#define XP_TPETRA_MATRIX_INSTANT(S,LO,GO,N) \
      TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( MatrixMatrix, Multiply_Tpetra, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( MatrixMatrix, BlockCrs, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N ) 

  // List of tests which run only with Epetra
#define XP_EPETRA_MATRIX_INSTANT(S,LO,GO,N) \
      TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( MatrixMatrix, Multiply_Epetra, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N ) \

  // List of tests which run only with Epetra64
#define XP_EPETRA64_MATRIX_INSTANT(S,LO,GO,N) \
      TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( MatrixMatrix, Multiply_Epetra64, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N )

#include <TpetraCore_config.h>
#include <TpetraCore_ETIHelperMacros.h>

TPETRA_ETI_MANGLING_TYPEDEFS()
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( XPETRA_TPETRA_TYPES )
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( XP_TPETRA_MATRIX_INSTANT )


#if defined(HAVE_XPETRA_EPETRA)

#include "Xpetra_Map.hpp" // defines EpetraNode
typedef Xpetra::EpetraNode EpetraNode;
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
XPETRA_EPETRA_TYPES(double,int,int,EpetraNode)
XP_EPETRA_MATRIX_INSTANT(double,int,int,EpetraNode)
#endif
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
typedef long long LongLong;
XPETRA_EPETRA_TYPES(double,int,LongLong,EpetraNode)
XP_EPETRA64_MATRIX_INSTANT(double,int,LongLong,EpetraNode)
#endif

#endif

} // end namespace

