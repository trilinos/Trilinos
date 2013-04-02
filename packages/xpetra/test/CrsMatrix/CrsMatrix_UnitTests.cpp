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
#ifdef HAVE_XPETRA_TPETRA
#include <Xpetra_TpetraCrsMatrix.hpp>
#endif
#ifdef HAVE_XPETRA_EPETRA
#include <Xpetra_EpetraCrsMatrix.hpp>
#endif
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_Exceptions.hpp>

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


  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, Apply, Scalar, LO, GO, Node )
  {
#ifdef HAVE_XPETRA_EPETRA

    typedef Xpetra::Map<LO, GO, Node> MapClass;
    typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;

    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();

    Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;

    // generate problem
    LO nEle = 63;
    const RCP<const MapClass> map = MapFactoryClass::Build(lib, nEle, 0, comm);

    RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > matrix =
        Xpetra::CrsMatrixFactory<Scalar,LO,GO,Node>::Build(map, 10);

    LO NumMyElements = map->getNodeNumElements();
    Teuchos::ArrayView<const GO> MyGlobalElements = map->getNodeElementList();

    for (LO i = 0; i < NumMyElements; ++i) {
        matrix->insertGlobalValues(MyGlobalElements[i],
                                Teuchos::tuple<GO>(MyGlobalElements[i]),
                                Teuchos::tuple<Scalar>(1.0) );
    }

    matrix->fillComplete();

    RCP<Xpetra::Vector<Scalar, LO, GO, Node> > vec =
        Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(map);

    vec->putScalar(1.0);

    RCP<Xpetra::Vector<Scalar, LO, GO, Node> > vec_sol =
        Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(matrix->getRangeMap());

    vec_sol->putScalar(0.0);

    matrix->apply(*vec, *vec_sol, Teuchos::NO_TRANS, 1.0, 0.0);

    vec_sol->putScalar(2.0);

    matrix->apply(*vec, *vec_sol, Teuchos::NO_TRANS, 1.0, -0.5);

    TEUCHOS_TEST_COMPARE(vec_sol->norm2(), <, 1e-16, out, success);
#endif
  }

#ifdef HAVE_XPETRA_EXPERIMENTAL
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, TpetraDeepCopy, Scalar, LO, GO, Node )
  {
#ifdef HAVE_XPETRA_TPETRA

    typedef Xpetra::Map<LO, GO, Node> MapClass;
    typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;

    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();

    Xpetra::UnderlyingLib lib = Xpetra::UseTpetra;

    // generate problem
    LO nEle = 63;
    const RCP<const MapClass> map = MapFactoryClass::Build(lib, nEle, 0, comm);

    RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > A =
        Xpetra::CrsMatrixFactory<Scalar,LO,GO,Node>::Build(map, 10);

    LO NumMyElements = map->getNodeNumElements();
    Teuchos::ArrayView<const GO> MyGlobalElements = map->getNodeElementList();

    for (LO i = 0; i < NumMyElements; ++i) {
        A->insertGlobalValues(MyGlobalElements[i],
                                Teuchos::tuple<GO>(MyGlobalElements[i]),
                                Teuchos::tuple<Scalar>(1.0) );
    }

    A->fillComplete();

    RCP<Xpetra::Vector<Scalar, LO, GO, Node> > v = Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(map);
    v->setSeed(8675309);
    v->randomize(true);

    RCP<Xpetra::Vector<Scalar, LO, GO, Node> > r = Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(A->getRangeMap());
    RCP<Xpetra::Vector<Scalar, LO, GO, Node> > rcopy = Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(A->getRangeMap());

    A->apply(*v, *r, Teuchos::NO_TRANS, 1.0, 0.0);

    RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > Acopy(new Xpetra::TpetraCrsMatrix<Scalar,LO,GO,Node>(*(Teuchos::rcp_static_cast<Xpetra::TpetraCrsMatrix<Scalar,LO,GO,Node> >(A))));
    A = Teuchos::null;

    Acopy->apply(*v, *rcopy, Teuchos::NO_TRANS, 1.0, 0.0);

    Teuchos::ArrayRCP<Scalar> rdata = r->getDataNonConst(0), rdatacopy = rcopy->getDataNonConst(0);
    Scalar s = Teuchos::ScalarTraits<Scalar>::zero();
    for (LO i = 0; i < NumMyElements; i++)
        s += Teuchos::ScalarTraits<Scalar>::magnitude(rdata[i] - rdatacopy[i]);
    TEUCHOS_TEST_COMPARE(s, <, 1e-16, out, success);
#endif
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, EpetraDeepCopy, Scalar, LO, GO, Node )
  {
#ifdef HAVE_XPETRA_EPETRA

    typedef Xpetra::Map<LO, GO, Node> MapClass;
    typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;

    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();

    Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;

    // generate problem
    LO nEle = 63;
    const RCP<const MapClass> map = MapFactoryClass::Build(lib, nEle, 0, comm);

    RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > A =
        Xpetra::CrsMatrixFactory<Scalar,LO,GO,Node>::Build(map, 10);

    LO NumMyElements = map->getNodeNumElements();
    Teuchos::ArrayView<const GO> MyGlobalElements = map->getNodeElementList();

    for (LO i = 0; i < NumMyElements; ++i) {
        A->insertGlobalValues(MyGlobalElements[i],
                                Teuchos::tuple<GO>(MyGlobalElements[i]),
                                Teuchos::tuple<Scalar>(1.0) );
    }

    A->fillComplete();

    RCP<Xpetra::Vector<Scalar, LO, GO, Node> > v = Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(map);
    v->setSeed(8675309);
    v->randomize(true);

    RCP<Xpetra::Vector<Scalar, LO, GO, Node> > r = Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(A->getRangeMap());
    RCP<Xpetra::Vector<Scalar, LO, GO, Node> > rcopy = Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(A->getRangeMap());

    A->apply(*v, *r, Teuchos::NO_TRANS, 1.0, 0.0);

    RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > Acopy(new Xpetra::EpetraCrsMatrix(*(Teuchos::rcp_static_cast<Xpetra::EpetraCrsMatrix>(A))));
    A = Teuchos::null;

    Acopy->apply(*v, *rcopy, Teuchos::NO_TRANS, 1.0, 0.0);

    Teuchos::ArrayRCP<Scalar> rdata = r->getDataNonConst(0), rdatacopy = rcopy->getDataNonConst(0);
    Scalar s = Teuchos::ScalarTraits<Scalar>::zero();
    for (LO i = 0; i < NumMyElements; i++)
        s += Teuchos::ScalarTraits<Scalar>::magnitude(rdata[i] - rdatacopy[i]);
    TEUCHOS_TEST_COMPARE(s, <, 1e-16, out, success);
#endif
  }
#endif // ifdef HAVE_XPETRA_EXPERIMENTAL




  //
  // INSTANTIATIONS
  //

#define UNIT_TEST_GROUP_ORDINAL( SC, LO, GO, Node )                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, Apply, SC, LO, GO, Node )
#define UNIT_TEST_GROUP_ORDINAL1( SC, LO, GO, Node )                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, TpetraDeepCopy, SC, LO, GO, Node )
#define UNIT_TEST_GROUP_ORDINAL2( SC, LO, GO, Node )                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, EpetraDeepCopy, SC, LO, GO, Node )

  typedef Kokkos::DefaultNode::DefaultNodeType DefaultNodeType;

  UNIT_TEST_GROUP_ORDINAL(double, int, int, DefaultNodeType)
#ifdef HAVE_XPETRA_EXPERIMENTAL
  UNIT_TEST_GROUP_ORDINAL1(double, int, int, DefaultNodeType)
  UNIT_TEST_GROUP_ORDINAL2(double, int, int, DefaultNodeType)
#endif

}

