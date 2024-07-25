// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Xpetra_UnitTestHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <KokkosSparse_CrsMatrix.hpp>

#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include <Xpetra_Parameters.hpp>

#include <Xpetra_Map.hpp>
#include <Xpetra_MapExtractor.hpp>
#include <Xpetra_MultiVectorFactory.hpp>  // taw: include MultiVectorFactory before VectorFactory for BlockedMultiVector
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_Exceptions.hpp>
#include "Teuchos_ScalarTraits.hpp"

namespace {

using Teuchos::RCP;

bool testMpi         = true;
double errorTolSlack = 1e+1;

RCP<const Teuchos::Comm<int> > getDefaultComm() {
  if (testMpi) {
    return Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
  }
  return rcp(new Teuchos::SerialComm<int>());
}

/////////////////////////////////////////////////////

TEUCHOS_STATIC_SETUP() {
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
  clp.setOption(
      "test-mpi", "test-serial", &testMpi,
      "Test MPI (if available) or force test of serial.  In a serial build,"
      " this option is ignored and a serial comm is always used.");
  clp.setOption(
      "error-tol-slack", &errorTolSlack,
      "Slack off of machine epsilon used to check test results");
}

//
// UNIT TESTS
//

TEUCHOS_UNIT_TEST_TEMPLATE_5_DECL(CrsMatrix, Constructor_Vector, M, Scalar, LO, GO, Node) {
  using MapClass        = Xpetra::Map<LO, GO, Node>;
  using MapFactoryClass = Xpetra::MapFactory<LO, GO, Node>;
  using STS             = Teuchos::ScalarTraits<Scalar>;
  typedef typename STS::magnitudeType MT;

  // get a comm and node
  RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  M testMap(1, 0, comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  // generate problem
  const LO nEle                 = 63;
  const RCP<const MapClass> map = MapFactoryClass::Build(lib, nEle, 0, comm);

  const RCP<Xpetra::Vector<Scalar, LO, GO, Node> > vec = Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(map);
  vec->randomize();
  RCP<Xpetra::Matrix<Scalar, LO, GO, Node> > mat = Xpetra::MatrixFactory<Scalar, LO, GO, Node>::Build(vec.getConst());

  const MT tol = 1e-12;

  TEST_ASSERT(!mat.is_null());
  TEST_EQUALITY(nEle, mat->getGlobalNumEntries());
  TEST_FLOATING_EQUALITY(vec->norm2(), mat->getFrobeniusNorm(), tol);

  const RCP<Xpetra::Vector<Scalar, LO, GO, Node> > diagonal = Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(map);

  mat->getLocalDiagCopy(*diagonal);
  TEST_FLOATING_EQUALITY(vec->norm2(), diagonal->norm2(), tol);
}

TEUCHOS_UNIT_TEST_TEMPLATE_5_DECL(CrsMatrix, Apply, M, Scalar, LO, GO, Node) {
  typedef Xpetra::Map<LO, GO, Node> MapClass;
  typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;

  // get a comm and node
  RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  M testMap(1, 0, comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  // generate problem
  LO nEle                       = 63;
  const RCP<const MapClass> map = MapFactoryClass::Build(lib, nEle, 0, comm);

  RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > matrix =
      Xpetra::CrsMatrixFactory<Scalar, LO, GO, Node>::Build(map, 10);

  LO NumMyElements                              = map->getLocalNumElements();
  Teuchos::ArrayView<const GO> MyGlobalElements = map->getLocalElementList();

  for (LO i = 0; i < NumMyElements; ++i) {
    matrix->insertGlobalValues(MyGlobalElements[i],
                               Teuchos::tuple<GO>(MyGlobalElements[i]),
                               Teuchos::tuple<Scalar>(1.0));
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
}

TEUCHOS_UNIT_TEST_TEMPLATE_5_DECL(CrsMatrix, ReplaceGlobalAndLocalValues, M, Scalar, LO, GO, Node) {
  typedef Xpetra::Map<LO, GO, Node> MapClass;
  typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;

  // get a comm and node
  RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  M testMap(1, 0, comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  // generate problem
  LO nEle                       = 63;
  const RCP<const MapClass> map = MapFactoryClass::Build(lib, nEle, 0, comm);

  LO NumMyElements                              = map->getLocalNumElements();
  GO NumGlobalElements                          = map->getGlobalNumElements();
  Teuchos::ArrayView<const GO> MyGlobalElements = map->getLocalElementList();

  RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > A =
      Xpetra::CrsMatrixFactory<Scalar, LO, GO, Node>::Build(map, 3);
  TEUCHOS_TEST_FOR_EXCEPTION(A->isFillComplete() == true || A->isFillActive() == false, std::runtime_error, "");

  for (size_t i = 0; i < static_cast<size_t>(NumMyElements); i++) {
    if (MyGlobalElements[i] == 0) {
      A->insertGlobalValues(MyGlobalElements[i],
                            Teuchos::tuple<GO>(MyGlobalElements[i], MyGlobalElements[i] + 1),
                            Teuchos::tuple<Scalar>(2.0, -1.0));
    } else if (MyGlobalElements[i] == NumGlobalElements - 1) {
      A->insertGlobalValues(MyGlobalElements[i],
                            Teuchos::tuple<GO>(MyGlobalElements[i] - 1, MyGlobalElements[i]),
                            Teuchos::tuple<Scalar>(-1.0, 2.0));
    } else {
      A->insertGlobalValues(MyGlobalElements[i],
                            Teuchos::tuple<GO>(MyGlobalElements[i] - 1, MyGlobalElements[i], MyGlobalElements[i] + 1),
                            Teuchos::tuple<Scalar>(-1.0, 2.0, -1.0));
    }
  }

  // test only in serial
  if (comm->getSize() == 1) {
    Teuchos::ArrayView<const GO> indices;
    Teuchos::ArrayView<const Scalar> values;

    A->getGlobalRowView(0, indices, values);

    const size_t nnz = indices.size();
    TEUCHOS_TEST_FOR_EXCEPTION(nnz != 2, std::runtime_error, "Wrong number of nonzeros in row with row gid 0.");
    Teuchos::Array<Scalar> newValues(nnz, 0.0);

    for (size_t j = 0; j < nnz; j++) {
      newValues[j] = 42;
    }

    A->replaceGlobalValues(0, indices, newValues);
  }

  A->fillComplete();
  TEUCHOS_TEST_FOR_EXCEPTION(A->isFillComplete() == false || A->isFillActive() == true, std::runtime_error, "");

  // edit matrix

  const LO nodeNumRows = A->getLocalNumRows();
  A->resumeFill();
  TEUCHOS_TEST_FOR_EXCEPTION(A->isFillComplete() == true || A->isFillActive() == false, std::runtime_error, "");

  for (LO i = 0; i < nodeNumRows; i++) {
    Teuchos::ArrayView<const LO> indices;
    Teuchos::ArrayView<const Scalar> values;

    A->getLocalRowView(i, indices, values);

    const size_t nnz = indices.size();
    TEUCHOS_TEST_FOR_EXCEPTION(nnz != 2 && nnz != 3, std::runtime_error, "Wrong number of nonzeros in row with row gid 0.");
    Teuchos::Array<Scalar> newValues(nnz, 0.0);

    for (size_t j = 0; j < nnz; j++) {
      if (indices[j] == i) /* diagonal term */
        newValues[j] = 4;
      else
        newValues[j] = -2;
    }

    A->replaceLocalValues(i, indices, newValues);
  }

  A->fillComplete();
  TEUCHOS_TEST_FOR_EXCEPTION(A->isFillComplete() == false || A->isFillActive() == true, std::runtime_error, "");
}

TEUCHOS_UNIT_TEST_TEMPLATE_5_DECL(CrsMatrix, leftScale, M, Scalar, LO, GO, Node) {
  typedef Xpetra::Map<LO, GO, Node> MapClass;
  typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;
  typedef Xpetra::Vector<Scalar, LO, GO, Node> VectorClass;
  typedef Xpetra::VectorFactory<Scalar, LO, GO, Node> VectorFactoryClass;
  const Scalar ZERO = Teuchos::ScalarTraits<Scalar>::zero();

  // get a comm and node
  RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  M testMap(1, 0, comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  // generate problem
  LO nEle                       = 63;
  const RCP<const MapClass> map = MapFactoryClass::Build(lib, nEle, 0, comm);

  TEST_EQUALITY(Teuchos::as<LO>(map->getGlobalNumElements()), nEle);

  LO NumMyElements                              = map->getLocalNumElements();
  GO NumGlobalElements                          = map->getGlobalNumElements();
  Teuchos::ArrayView<const GO> MyGlobalElements = map->getLocalElementList();

  RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > A =
      Xpetra::CrsMatrixFactory<Scalar, LO, GO, Node>::Build(map, 3);
  TEUCHOS_TEST_FOR_EXCEPTION(A->isFillComplete() == true || A->isFillActive() == false, std::runtime_error, "");

  for (size_t i = 0; i < static_cast<size_t>(NumMyElements); i++) {
    if (MyGlobalElements[i] == 0) {
      A->insertGlobalValues(MyGlobalElements[i],
                            Teuchos::tuple<GO>(MyGlobalElements[i], MyGlobalElements[i] + 1),
                            Teuchos::tuple<Scalar>(2.0, -1.0));
    } else if (MyGlobalElements[i] == NumGlobalElements - 1) {
      A->insertGlobalValues(MyGlobalElements[i],
                            Teuchos::tuple<GO>(MyGlobalElements[i] - 1, MyGlobalElements[i]),
                            Teuchos::tuple<Scalar>(-1.0, 2.0));
    } else {
      A->insertGlobalValues(MyGlobalElements[i],
                            Teuchos::tuple<GO>(MyGlobalElements[i] - 1, MyGlobalElements[i], MyGlobalElements[i] + 1),
                            Teuchos::tuple<Scalar>(-1.0, 2.0, -1.0));
    }
  }
  A->fillComplete();
  TEST_EQUALITY(A->isFillComplete(), true);

  RCP<VectorClass> s = VectorFactoryClass::Build(map, true);
  {
    Teuchos::ArrayRCP<Scalar> sd = s->getDataNonConst(0);
    for (LO i = 0; i < NumMyElements; ++i) {
      sd[i] = Teuchos::as<Scalar>(map->getGlobalElement(i));
    }
  }

  A->leftScale(*s);

#ifdef HAVE_XPETRA_TPETRA
  Kokkos::fence();
#endif

  for (size_t i = 0; i < static_cast<size_t>(NumMyElements); i++) {
    if (MyGlobalElements[i] == 0) {
      Teuchos::ArrayView<const LO> indices;
      Teuchos::ArrayView<const Scalar> values;
      A->getLocalRowView(map->getLocalElement(MyGlobalElements[i]), indices, values);
      TEST_EQUALITY(indices.size(), 2);
      TEST_EQUALITY(values[0], ZERO);
      TEST_EQUALITY(values[1], ZERO);
    } else if (MyGlobalElements[i] == NumGlobalElements - 1) {
      if (map->isNodeGlobalElement(MyGlobalElements[i])) {
        Teuchos::ArrayView<const LO> indices;
        Teuchos::ArrayView<const Scalar> values;
        A->getLocalRowView(map->getLocalElement(MyGlobalElements[i]), indices, values);
        TEST_EQUALITY(indices.size(), 2);
        TEST_EQUALITY(values[0], Teuchos::as<Scalar>(-1.0) * Teuchos::as<Scalar>(MyGlobalElements[i]));
        TEST_EQUALITY(values[1], Teuchos::as<Scalar>(2.0) * Teuchos::as<Scalar>(MyGlobalElements[i]));
      }
    } else {
      // skip check on first row of each processor due to different local ordering (communication overlap)
      if (map->isNodeGlobalElement(MyGlobalElements[i]) && map->getLocalElement(MyGlobalElements[i]) != 0) {
        Teuchos::ArrayView<const LO> indices;
        Teuchos::ArrayView<const Scalar> values;
        A->getLocalRowView(map->getLocalElement(MyGlobalElements[i]), indices, values);
        TEST_EQUALITY(indices.size(), 3);
        TEST_EQUALITY(values[0], Teuchos::as<Scalar>(-1.0) * Teuchos::as<Scalar>(MyGlobalElements[i]));
        TEST_EQUALITY(values[1], Teuchos::as<Scalar>(2.0) * Teuchos::as<Scalar>(MyGlobalElements[i]));
        TEST_EQUALITY(values[2], Teuchos::as<Scalar>(-1.0) * Teuchos::as<Scalar>(MyGlobalElements[i]));
      }
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_5_DECL(CrsMatrix, rightScale, M, Scalar, LO, GO, Node) {
  typedef Xpetra::Map<LO, GO, Node> MapClass;
  typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;
  typedef Xpetra::Vector<Scalar, LO, GO, Node> VectorClass;
  typedef Xpetra::VectorFactory<Scalar, LO, GO, Node> VectorFactoryClass;
  const Scalar ZERO = Teuchos::ScalarTraits<Scalar>::zero();
  const Scalar ONE  = Teuchos::ScalarTraits<Scalar>::one();

  // get a comm and node
  RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  M testMap(1, 0, comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  // generate problem
  LO nEle                       = 63;
  const RCP<const MapClass> map = MapFactoryClass::Build(lib, nEle, 0, comm);

  TEST_EQUALITY(Teuchos::as<LO>(map->getGlobalNumElements()), nEle);

  LO NumMyElements                              = map->getLocalNumElements();
  GO NumGlobalElements                          = map->getGlobalNumElements();
  Teuchos::ArrayView<const GO> MyGlobalElements = map->getLocalElementList();

  RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > A =
      Xpetra::CrsMatrixFactory<Scalar, LO, GO, Node>::Build(map, 3);
  TEUCHOS_TEST_FOR_EXCEPTION(A->isFillComplete() == true || A->isFillActive() == false, std::runtime_error, "");

  for (size_t i = 0; i < static_cast<size_t>(NumMyElements); i++) {
    if (MyGlobalElements[i] == 0) {
      A->insertGlobalValues(MyGlobalElements[i],
                            Teuchos::tuple<GO>(MyGlobalElements[i], MyGlobalElements[i] + 1),
                            Teuchos::tuple<Scalar>(2.0, -1.0));
    } else if (MyGlobalElements[i] == NumGlobalElements - 1) {
      A->insertGlobalValues(MyGlobalElements[i],
                            Teuchos::tuple<GO>(MyGlobalElements[i] - 1, MyGlobalElements[i]),
                            Teuchos::tuple<Scalar>(-1.0, 2.0));
    } else {
      A->insertGlobalValues(MyGlobalElements[i],
                            Teuchos::tuple<GO>(MyGlobalElements[i] - 1, MyGlobalElements[i], MyGlobalElements[i] + 1),
                            Teuchos::tuple<Scalar>(-1.0, 2.0, -1.0));
    }
  }
  A->fillComplete();
  TEST_EQUALITY(A->isFillComplete(), true);

  RCP<VectorClass> s = VectorFactoryClass::Build(map, true);
  {
    Teuchos::ArrayRCP<Scalar> sd = s->getDataNonConst(0);
    for (LO i = 0; i < NumMyElements; ++i) {
      sd[i] = Teuchos::as<Scalar>(map->getGlobalElement(i));
    }
  }

  A->rightScale(*s);

#ifdef HAVE_XPETRA_TPETRA
  Kokkos::fence();
#endif

  for (size_t i = 0; i < static_cast<size_t>(NumMyElements); i++) {
    if (MyGlobalElements[i] == 0) {
      Teuchos::ArrayView<const LO> indices;
      Teuchos::ArrayView<const Scalar> values;
      A->getLocalRowView(map->getLocalElement(MyGlobalElements[i]), indices, values);
      TEST_EQUALITY(indices.size(), 2);
      TEST_EQUALITY(values[0], ZERO);
      TEST_EQUALITY(values[1], -ONE);
    } else if (MyGlobalElements[i] == NumGlobalElements - 1) {
      if (map->isNodeGlobalElement(MyGlobalElements[i])) {
        Teuchos::ArrayView<const LO> indices;
        Teuchos::ArrayView<const Scalar> values;
        A->getLocalRowView(map->getLocalElement(MyGlobalElements[i]), indices, values);
        TEST_EQUALITY(indices.size(), 2);
        TEST_EQUALITY(values[0], Teuchos::as<Scalar>(-1.0) * Teuchos::as<Scalar>(MyGlobalElements[i] - 1));
        TEST_EQUALITY(values[1], Teuchos::as<Scalar>(2.0) * Teuchos::as<Scalar>(MyGlobalElements[i]));
      }
    } else {
      // skip check on first row of each processor due to different local ordering (communication overlap)
      if (map->isNodeGlobalElement(MyGlobalElements[i]) && map->getLocalElement(MyGlobalElements[i]) != 0) {
        Teuchos::ArrayView<const LO> indices;
        Teuchos::ArrayView<const Scalar> values;
        A->getLocalRowView(map->getLocalElement(MyGlobalElements[i]), indices, values);
        TEST_EQUALITY(indices.size(), 3);
        TEST_EQUALITY(values[0], Teuchos::as<Scalar>(-1.0) * Teuchos::as<Scalar>(MyGlobalElements[i] - 1));
        TEST_EQUALITY(values[1], Teuchos::as<Scalar>(2.0) * Teuchos::as<Scalar>(MyGlobalElements[i]));
        TEST_EQUALITY(values[2], Teuchos::as<Scalar>(-1.0) * Teuchos::as<Scalar>(MyGlobalElements[i] + 1));
      }
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_5_DECL(CrsMatrix, replaceDiagonal, M, Scalar, LO, GO, Node) {
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename STS::magnitudeType MT;
  const Scalar SC_ONE = STS::one();
  const GO GO_ONE     = Teuchos::OrdinalTraits<GO>::one();
  // get a comm and node
  RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  M testMap(1, 0, comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  // generate problem
  GO nEle = Teuchos::as<GO>(2 * comm->getSize());
  const RCP<const Xpetra::Map<LO, GO, Node> > map =
      Xpetra::MapFactory<LO, GO, Node>::Build(lib, nEle, 2, 0, comm);

  RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > A =
      Xpetra::CrsMatrixFactory<Scalar, LO, GO, Node>::Build(map, 3);
  const Scalar rankAsScalar = static_cast<Scalar>(static_cast<MT>(comm->getRank()));

  Teuchos::Array<Scalar> vals = {{SC_ONE, rankAsScalar + SC_ONE, SC_ONE}};
  for (size_t lclRowIdx = 0; lclRowIdx < 2; ++lclRowIdx) {
    const GO gblRowIdx      = Teuchos::as<GO>(2 * comm->getRank() + lclRowIdx);
    Teuchos::Array<GO> cols = {{gblRowIdx - GO_ONE, gblRowIdx, gblRowIdx + GO_ONE}};

    if ((comm->getRank() == 0) && (lclRowIdx == 0)) {  // First row of the matrix
      A->insertGlobalValues(gblRowIdx, cols(1, 2), vals(1, 2));
    } else if ((comm->getRank() == comm->getSize() - 1) && (lclRowIdx == 1)) {  // Last row of the matrix
      A->insertGlobalValues(gblRowIdx, cols(0, 2), vals(0, 2));
    } else {
      A->insertGlobalValues(gblRowIdx, cols(), vals());
    }
  }

  A->fillComplete();
  TEST_ASSERT(A->isFillComplete());

  comm->barrier();
  {
    /* Replace the diagonal of the matrix by the ID of the owning MPI rank
     *
     * 1. Create map
     * 2. Create vector with new diagonal values
     * 3. Replace the diagonal
     * 4. Test for
     *    - successful replacement of diagonal values
     *    - unchanged off-diagonal values (not implemented yet)
     */

    // Create vector with new diagonal values
    RCP<Xpetra::Vector<Scalar, LO, GO, Node> > newDiag =
        Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(A->getRowMap(), true);
    newDiag->putScalar(rankAsScalar);

    // Replace the diagonal
    A->replaceDiag(*newDiag);

    // Tests
    {
      RCP<Xpetra::Vector<Scalar, LO, GO, Node> > diagCopy =
          Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(A->getRowMap(), true);
      A->getLocalDiagCopy(*diagCopy);

      Teuchos::ArrayRCP<const Scalar> diagCopyData = diagCopy->getData(0);

      for (size_t i = 0; i < static_cast<size_t>(diagCopyData.size()); ++i)
        TEST_EQUALITY_CONST(diagCopyData[i], rankAsScalar);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, Constructor_Epetra, Scalar, LO, GO, Node) {
#ifdef HAVE_XPETRA_EPETRA

  // get a comm and node
  RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  {
    typedef Xpetra::EpetraMapT<GO, Node> mm;
    TEST_NOTHROW(mm(10, 0, comm));
    typedef Xpetra::EpetraCrsMatrixT<GO, Node> mx;
    TEST_NOTHROW(mx(Teuchos::rcp(new mm(10, 0, comm)), 0));
  }

#if defined(HAVE_XPETRA_TPETRA) && defined(HAVE_TPETRA_INST_PTHREAD)
  {
    typedef Xpetra::EpetraMapT<GO, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> mm;
    TEST_THROW(mm(10, 0, comm), Xpetra::Exceptions::RuntimeError);
    typedef Xpetra::EpetraCrsMatrixT<GO, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> mx;
    TEST_THROW(mx(Teuchos::null, 0), Xpetra::Exceptions::RuntimeError);
  }
#endif
#if defined(HAVE_XPETRA_TPETRA) && defined(HAVE_TPETRA_INST_CUDA)
  {
    typedef Xpetra::EpetraMapT<GO, Tpetra::KokkosCompat::KokkosCudaWrapperNode> mm;
    TEST_THROW(mm(10, 0, comm), Xpetra::Exceptions::RuntimeError);
    typedef Xpetra::EpetraCrsMatrixT<GO, Tpetra::KokkosCompat::KokkosCudaWrapperNode> mx;
    TEST_THROW(mx(Teuchos::null, 0), Xpetra::Exceptions::RuntimeError);
  }
#endif
#if defined(HAVE_XPETRA_TPETRA) && defined(HAVE_TPETRA_INST_HIP)
  {
    typedef Xpetra::EpetraMapT<GO, Tpetra::KokkosCompat::KokkosHIPWrapperNode> mm;
    TEST_THROW(mm(10, 0, comm), Xpetra::Exceptions::RuntimeError);
    typedef Xpetra::EpetraCrsMatrixT<GO, Tpetra::KokkosCompat::KokkosHIPWrapperNode> mx;
    TEST_THROW(mx(Teuchos::null, 0), Xpetra::Exceptions::RuntimeError);
  }
#endif

#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, Epetra_ReplaceLocalValues, Scalar, LO, GO, Node) {
#ifdef HAVE_XPETRA_EPETRA

  typedef Xpetra::Map<LO, GO, Node> MapClass;
  typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;

  // get a comm and node
  RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;

  // generate problem
  LO nEle                       = 63;
  const RCP<const MapClass> map = MapFactoryClass::Build(lib, nEle, 0, comm);

  RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > matrix =
      Xpetra::CrsMatrixFactory<Scalar, LO, GO, Node>::Build(map, 10);

  LO NumMyElements                              = map->getLocalNumElements();
  Teuchos::ArrayView<const GO> MyGlobalElements = map->getLocalElementList();

  for (LO i = 0; i < NumMyElements; ++i) {
    matrix->insertGlobalValues(MyGlobalElements[i],
                               Teuchos::tuple<GO>(MyGlobalElements[i]),
                               Teuchos::tuple<Scalar>(1.0));
  }

  matrix->fillComplete();
  matrix->resumeFill();

  Teuchos::Array<LO> indout(1, 0);
  Teuchos::Array<Scalar> valout(1, 5.0);
  matrix->replaceLocalValues(0, indout.view(0, indout.size()), valout.view(0, valout.size()));
  matrix->fillComplete();

  RCP<Xpetra::Vector<Scalar, LO, GO, Node> > vec =
      Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(map);

  vec->putScalar(1.0);

  RCP<Xpetra::Vector<Scalar, LO, GO, Node> > vec_sol =
      Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(matrix->getRangeMap());

  vec_sol->putScalar(0.0);

  matrix->apply(*vec, *vec_sol, Teuchos::NO_TRANS, 1.0, 0.0);

  RCP<Xpetra::Vector<Scalar, LO, GO, Node> > vectest =
      Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(map);
  vectest->putScalar(1.0);
  Teuchos::ArrayRCP<Scalar> vectestData = vectest->getDataNonConst(0);
  vectestData[0]                        = 5.0;

  vec_sol->update(-1.0, *vectest, 1.0);

  TEUCHOS_TEST_COMPARE(vec_sol->norm2(), <, 1e-16, out, success);
#endif
}

// just a copy of the Epetra_ReplaceLocalValues test for Tpetra
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, Tpetra_ReplaceLocalValues, Scalar, LO, GO, Node) {
#ifdef HAVE_XPETRA_TPETRA
  using std::endl;
  using Teuchos::outArg;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  typedef Xpetra::Map<LO, GO, Node> map_type;
  typedef Xpetra::MapFactory<LO, GO, Node> map_factory_type;
  typedef Xpetra::CrsMatrixFactory<Scalar, LO, GO, Node> crs_matrix_factory_type;
  typedef Xpetra::CrsMatrix<Scalar, LO, GO, Node> crs_matrix_type;
  typedef Xpetra::VectorFactory<Scalar, LO, GO, Node> vec_factory_type;
  typedef Xpetra::Vector<Scalar, LO, GO, Node> vec_type;
  typedef typename Teuchos::Array<LO>::size_type size_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename STS::magnitudeType MT;
  typedef Teuchos::ScalarTraits<MT> STM;
  const Scalar ONE  = STS::one();
  const Scalar FIVE = ONE + ONE + ONE + ONE + ONE;

  out << "Tpetra replaceLocalValues test" << endl;

  // get a comm and node
  RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  Xpetra::UnderlyingLib lib = Xpetra::UseTpetra;
  // Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

  out << "Create Map and matrix" << endl;

  // generate problem
  LO nEle                 = 63;
  RCP<const map_type> map = map_factory_type::Build(lib, nEle, 0, comm);

  RCP<crs_matrix_type> matrix                   = crs_matrix_factory_type::Build(map, 10);
  const LO NumMyElements                        = map->getLocalNumElements();
  Teuchos::ArrayView<const GO> MyGlobalElements = map->getLocalElementList();

  // Make the matrix the identity matrix.
  out << "Fill matrix by calling insertGlobalValues" << endl;
  for (LO i = 0; i < NumMyElements; ++i) {
    matrix->insertGlobalValues(MyGlobalElements[i],
                               Teuchos::tuple<GO>(MyGlobalElements[i]),
                               Teuchos::tuple<Scalar>(1.0));
  }

  out << "Call fillComplete and resumeFill on matrix" << endl;
  matrix->fillComplete();
  matrix->resumeFill();

  // Change the 0,0 local entry, on each process, to be 5.0.
  out << "Modify entries of the matrix using replaceLocalValues, "
         "and test the result before calling fillComplete"
      << endl;
  Teuchos::Array<LO> indout(1, 0);
  Teuchos::Array<Scalar> valout(1, 5.0);

  // Every process should have a local row index 0.
  TEST_ASSERT(map->isNodeLocalElement(0));
  TEST_ASSERT(!matrix->getColMap().is_null());

  if (map->isNodeLocalElement(0) && !matrix->getColMap().is_null()) {
    bool validLocalColumnIndices = true;
    for (size_type k = 0; k < indout.size(); ++k) {
      if (!matrix->getColMap()->isNodeLocalElement(indout[k])) {
        validLocalColumnIndices = false;
        break;
      }
    }
    // Every process should have a local column index 0.
    TEST_ASSERT(validLocalColumnIndices);
    if (validLocalColumnIndices) {
      // Make sure that we are changing the first diagonal entry on
      // this process.  We determine whether a matrix is diagonal
      // using global indices.
      TEST_ASSERT(matrix->getColMap()->getGlobalElement(indout[0]) ==
                  map->getGlobalElement(0));
      // Replace the local (0,0) entry with valout[0].  We know from
      // the above test that the local (0,0) entry is the first
      // diagonal entry on the calling process.
      matrix->replaceLocalValues(0, indout.view(0, indout.size()),
                                 valout.view(0, valout.size()));
    }

    // Make sure that replaceLocalValues worked, by getting the
    // values in the local row 0.
    const size_t numEnt = matrix->getNumEntriesInLocalRow(0);
    TEST_EQUALITY_CONST(numEnt, static_cast<size_t>(1));

    if (numEnt == static_cast<size_t>(1)) {
      Teuchos::Array<LO> ind(numEnt);
      Teuchos::Array<Scalar> val(numEnt);
      size_t numEntOut = 0;
      matrix->getLocalRowCopy(0, ind(), val(), numEntOut);
      TEST_EQUALITY(numEnt, numEntOut);

      if (numEntOut == static_cast<size_t>(1)) {
        TEST_EQUALITY(ind[0], 0);
        TEST_EQUALITY(val[0], FIVE);
      }
    }
  }

  // Make sure that all processes got this far.
  int lclSuccess = success ? 1 : 0;
  int gblSuccess = 0;
  reduceAll<int, int>(*comm, REDUCE_MIN, lclSuccess, outArg(gblSuccess));
  success = success && (gblSuccess == 1);
  TEST_EQUALITY_CONST(gblSuccess, 1);

  out << "Call fillComplete on matrix for the second time" << endl;
  matrix->fillComplete();

  out << "Test the result of replaceLocalValues after fillComplete" << endl;
  if (map->isNodeLocalElement(0)) {
    // Make sure that replaceLocalValues worked, by getting the
    // values in the local row 0.
    const size_t numEnt = matrix->getNumEntriesInLocalRow(0);
    TEST_EQUALITY_CONST(numEnt, static_cast<size_t>(1));

    if (numEnt == static_cast<size_t>(1)) {
      Teuchos::Array<LO> ind(numEnt);
      Teuchos::Array<Scalar> val(numEnt);
      size_t numEntOut = 0;
      matrix->getLocalRowCopy(0, ind(), val(), numEntOut);
      TEST_EQUALITY(numEnt, numEntOut);

      if (numEntOut == static_cast<size_t>(1)) {
        TEST_EQUALITY(ind[0], 0);
        TEST_EQUALITY(val[0], FIVE);
      }
    }
  }

  RCP<vec_type> vec = vec_factory_type::Build(map);
  vec->putScalar(1.0);
  out << "Test that vec->putScalar(1.0) filled vec with ones" << endl;
  {
    const MT N             = static_cast<MT>(vec->getGlobalLength());
    const MT expectedNorm2 = STM::squareroot(N);
    const MT actualNorm2   = vec->norm2();
    TEST_EQUALITY(actualNorm2, expectedNorm2);
  }

  RCP<const map_type> rangeMap = matrix->getRangeMap();
  TEST_ASSERT(!rangeMap.is_null());
  RCP<vec_type> vec_sol = vec_factory_type::Build(rangeMap);
  vec_sol->putScalar(0.0);
  out << "Test that vec_sol->putScalar(0.0) filled vec with zeros" << endl;
  {
    const MT expectedNorm2 = STM::zero();
    const MT actualNorm2   = vec_sol->norm2();
    TEST_EQUALITY(actualNorm2, expectedNorm2);
  }

  // Compute vec_sol := matrix*vec.  The result _should_ be a vector
  // of ones everywhere, except for the entry at local index zero
  // (on every process), which should be 5.0.
  matrix->apply(*vec, *vec_sol, Teuchos::NO_TRANS, 1.0, 0.0);
  if (rangeMap->getLocalNumElements() > 0) {
    // Test this both for a const view and for a nonconst view.
    // This may also be a test for {T,X}petra::MultiVector::getData
    // and {T,X}petra::MultiVector::getDataNonConst.

    // Create the const view.
    Teuchos::ArrayRCP<const Scalar> outData = vec_sol->getData(0);
    TEST_ASSERT(outData.size() == static_cast<size_type>(rangeMap->getLocalNumElements()));
    if (outData.size() == static_cast<size_type>(rangeMap->getLocalNumElements()) &&
        outData.size() > static_cast<size_type>(0)) {
      TEST_EQUALITY(outData[0], FIVE);
    }
    if (rangeMap->getLocalNumElements() > static_cast<size_t>(1)) {
      bool allOnes = true;
      for (size_t k = 1; k < size_t(rangeMap->getLocalNumElements()); ++k) {
        if (outData[k] != ONE) {
          allOnes = false;
        }
      }
      TEST_ASSERT(allOnes);
    }

    // Invalidate the const view, before creating a nonconst view.
    outData = Teuchos::null;
    // Create the nonconst view.
    Teuchos::ArrayRCP<Scalar> outDataNonConst = vec_sol->getDataNonConst(0);
    TEST_ASSERT(outDataNonConst.size() == static_cast<size_type>(rangeMap->getLocalNumElements()));
    if (outDataNonConst.size() == static_cast<size_type>(rangeMap->getLocalNumElements()) &&
        outDataNonConst.size() > static_cast<size_type>(0)) {
      TEST_EQUALITY(outDataNonConst[0], FIVE);
    }
    if (rangeMap->getLocalNumElements() > static_cast<size_t>(1)) {
      bool allOnes = true;
      for (size_type k = 1; k < static_cast<size_type>(rangeMap->getLocalNumElements()); ++k) {
        if (outDataNonConst[k] != ONE) {
          allOnes = false;
        }
      }
      TEST_ASSERT(allOnes);
    }
  }

  lclSuccess = success ? 1 : 0;
  gblSuccess = 0;
  reduceAll<int, int>(*comm, REDUCE_MIN, lclSuccess, outArg(gblSuccess));
  success = success && (gblSuccess == 1);
  TEST_EQUALITY_CONST(gblSuccess, 1);

  if (gblSuccess == 1) {
    out << "Vector result is correct" << endl;
  }

  RCP<vec_type> vectest = vec_factory_type::Build(map);
  vectest->putScalar(1.0);
  {
    Teuchos::ArrayRCP<Scalar> vectestData = vectest->getDataNonConst(0);
    vectestData[0]                        = 5.0;
  }

  vec_sol->update(-1.0, *vectest, 1.0);

  // FIXME (mfh 23 Feb 2020) Tolerance should depend on Scalar here.
  TEUCHOS_TEST_COMPARE(vec_sol->norm2(), <, 1e-16, out, success);

  // Make sure that all processes got this far.
  lclSuccess = success ? 1 : 0;
  gblSuccess = 0;
  reduceAll<int, int>(*comm, REDUCE_MIN, lclSuccess, outArg(gblSuccess));
  success = success && (gblSuccess == 1);
  TEST_EQUALITY_CONST(gblSuccess, 1);
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, TpetraDeepCopy, Scalar, LO, GO, Node) {
#ifdef HAVE_XPETRA_TPETRA
  typedef Xpetra::Map<LO, GO, Node> MapClass;
  typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename STS::magnitudeType magnitude_type;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;

  // get a comm and node
  RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  {
    using std::cerr;
    using std::endl;

    std::ostringstream os;
    const int myRank = comm->getRank();
    os << "Process " << myRank << ": (CrsMatrix, TpetraDeepCopy) test" << endl;
    cerr << os.str();
  }

  Xpetra::UnderlyingLib lib = Xpetra::UseTpetra;

  // Create a Map, which will be the row, domain, and range Map of the matrix A.
  LO nEle                       = 63;
  const RCP<const MapClass> map = MapFactoryClass::Build(lib, nEle, 0, comm);

  {
    using std::cerr;
    using std::endl;

    std::ostringstream os;
    const int myRank = comm->getRank();
    os << "Process " << myRank << ": Creating matrix" << endl;
    cerr << os.str();
  }

  RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > A =
      Xpetra::CrsMatrixFactory<Scalar, LO, GO, Node>::Build(map, 10);

  LO NumMyElements                              = map->getLocalNumElements();
  Teuchos::ArrayView<const GO> MyGlobalElements = map->getLocalElementList();

  {
    using std::cerr;
    using std::endl;

    std::ostringstream os;
    const int myRank = comm->getRank();
    os << "Process " << myRank << ": Filling matrix" << endl;
    cerr << os.str();
  }

  // Make A the identity matrix.
  for (LO i = 0; i < NumMyElements; ++i) {
    A->insertGlobalValues(MyGlobalElements[i],
                          Teuchos::tuple<GO>(MyGlobalElements[i]),
                          Teuchos::tuple<Scalar>(1.0));
  }

  {
    using std::cerr;
    using std::endl;

    std::ostringstream os;
    const int myRank = comm->getRank();
    os << "Process " << myRank << ": Calling fillComplete on matrix A" << endl;
    cerr << os.str();
  }

  A->fillComplete();

  {
    using std::cerr;
    using std::endl;

    std::ostringstream os;
    const int myRank = comm->getRank();
    os << "Process " << myRank << ": Building an input Vector" << endl;
    cerr << os.str();
  }

  RCP<Xpetra::Vector<Scalar, LO, GO, Node> > v = Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(A->getRangeMap());
  v->setSeed(8675309);
  v->randomize(true);

  {
    using std::cerr;
    using std::endl;

    std::ostringstream os;
    const int myRank = comm->getRank();
    os << "Process " << myRank << ": Testing that Xpetra::Vector::operator= does a deep copy" << endl;
    cerr << os.str();
  }

  // Remember the norm of v, to make sure that neither apply() call changes it.
  const magnitude_type v_norm = v->norm2();

  // Keep a copy of v, to test that neither apply() call changes it.
  RCP<Xpetra::Vector<Scalar, LO, GO, Node> > vcopy =
      Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(map);
  // Xpetra's operator= does a deep copy, like Epetra, but unlike
  // Tpetra (as of early 2014).
  *vcopy = *v;

  // Make sure that vcopy and v have the same norm.  It's OK for the
  // norms to be slightly different, due to nondeterminism in
  // parallel collectives.
  const magnitude_type vcopy_norm = vcopy->norm2();

  {
    using std::cerr;
    using std::endl;

    std::ostringstream os;
    const int myRank = comm->getRank();
    os << "Process " << myRank << ": v->norm2() = " << v_norm
       << ", vcopy->norm2() = " << vcopy_norm << endl;
    cerr << os.str();
  }

  const magnitude_type norm_tol =
      static_cast<magnitude_type>(map->getGlobalNumElements()) * STM::eps();
  TEUCHOS_TEST_COMPARE(STM::magnitude(v_norm - vcopy_norm), <, norm_tol, out, success);

  // Make sure that if you change vcopy, v doesn't change.
  // That is, vcopy must be a true deep copy of v.
  {
    Teuchos::ArrayRCP<Scalar> vcopy_data = vcopy->getDataNonConst(0);
    if (NumMyElements != 0) {
      vcopy_data[0] += static_cast<magnitude_type>(10000.0);
    }
    // Destroy the view, so that the changes get written back to the Vector.
    vcopy_data = Teuchos::null;

    // Adding 10000 to an entry had better change the 2-norm by at least sqrt(10000) = 100.
    const magnitude_type norm_tol2 = static_cast<magnitude_type>(100.0);
    TEUCHOS_TEST_COMPARE(STM::magnitude(vcopy_norm - vcopy->norm2()), >, norm_tol2, out, success);

    // Restore the original vcopy, by doing a deep copy again.
    // Xpetra's operator= does a deep copy, like Epetra, but unlike
    // Tpetra (as of early 2014).
    *vcopy = *v;

    // Make sure the original copy got restored.
    TEUCHOS_TEST_COMPARE(STM::magnitude(vcopy_norm - vcopy->norm2()), <, norm_tol, out, success);
  }

  // r and rcopy are distinct Vectors with the same Map, namely the
  // range Map of A.  All the Vectors v, r, and rcopy are distinct.
  RCP<Xpetra::Vector<Scalar, LO, GO, Node> > r =
      Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(A->getRangeMap());
  RCP<Xpetra::Vector<Scalar, LO, GO, Node> > rcopy =
      Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(A->getRangeMap());

  {
    using std::cerr;
    using std::endl;

    std::ostringstream os;
    const int myRank = comm->getRank();
    os << "Process " << myRank << ": Applying matrix A" << endl;
    cerr << os.str();
  }

  // r := A * v.
  A->apply(*v, *r, Teuchos::NO_TRANS, STS::one(), STS::zero());

  // Since A is the identity matrix, after the above line finishes,
  // r and v should be exactly equal.  This should be true even in
  // finite-precision arithmetic.  Test this here.
  {
    RCP<Xpetra::Vector<Scalar, LO, GO, Node> > diff =
        Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(A->getRangeMap());

    // diff := 1*v + (-1)*r.
    diff->update(STS::one(), *v, -STS::one(), *r, STS::zero());
    Teuchos::Array<magnitude_type> norms(1);
    diff->norm2(norms());
    // The norm of v - r must be _exactly_ zero.
    TEST_EQUALITY(norms[0], STM::zero());
  }

  // Make sure that the above apply() call didn't change v, by
  // testing against vcopy.
  {
    RCP<Xpetra::Vector<Scalar, LO, GO, Node> > diff =
        Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(A->getRangeMap());

    // diff := 1*v + (-1)*vcopy.
    diff->update(STS::one(), *v, -STS::one(), *vcopy, STS::zero());
    Teuchos::Array<magnitude_type> norms(1);
    diff->norm2(norms());
    // The norm of v - vcopy must be _exactly_ zero.
    TEST_EQUALITY(norms[0], STM::zero());
  }
  // TODO Make sure norm of v didn't change

  {
    using std::cerr;
    using std::endl;

    std::ostringstream os;
    const int myRank = comm->getRank();
    os << "Process " << myRank << ": Deep-copying matrix" << endl;
    cerr << os.str();
  }

  using Teuchos::rcp_static_cast;
  // NOTE (mfh 24 Apr 2014): This invokes the
  // Xpetra::TpetraCrsMatrix copy constructor, on line 329 of
  // Xpetra_TpetraCrsMatrix.hpp (as of 24 Apr 2014).  That in turn
  // calls Tpetra::CrsMatrix::clone().
  RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > Acopy =
      rcp(new Xpetra::TpetraCrsMatrix<Scalar, LO, GO, Node>(*(rcp_static_cast<Xpetra::TpetraCrsMatrix<Scalar, LO, GO, Node> >(A))));

  // Make sure that A and Acopy have the same gross properties.  For
  // example, they must be both fill complete and locally indexed,
  // and their four Maps must match.

  const bool bothFillComplete = A->isFillComplete() && Acopy->isFillComplete();
  TEST_EQUALITY_CONST(bothFillComplete, true);

  const bool bothLocallyIndexed = A->isLocallyIndexed() && Acopy->isLocallyIndexed();
  TEST_EQUALITY_CONST(bothLocallyIndexed, true);

  const bool bothNotGloballyIndexed = !A->isGloballyIndexed() && !Acopy->isGloballyIndexed();
  TEST_EQUALITY_CONST(bothNotGloballyIndexed, true);

  const bool rowMapsMatch = A->getRowMap()->isSameAs(*(Acopy->getRowMap()));
  TEST_EQUALITY_CONST(rowMapsMatch, true);

  const bool colMapsMatch = A->getColMap()->isSameAs(*(Acopy->getColMap()));
  TEST_EQUALITY_CONST(colMapsMatch, true);

  const bool domainMapsMatch = A->getDomainMap()->isSameAs(*(Acopy->getDomainMap()));
  TEST_EQUALITY_CONST(domainMapsMatch, true);

  const bool rangeMapsMatch = A->getRangeMap()->isSameAs(*(Acopy->getRangeMap()));
  TEST_EQUALITY_CONST(rangeMapsMatch, true);

  TEST_EQUALITY(A->getGlobalNumRows(), Acopy->getGlobalNumRows());
  TEST_EQUALITY(A->getGlobalNumCols(), Acopy->getGlobalNumCols());
  TEST_EQUALITY(A->getGlobalNumEntries(), Acopy->getGlobalNumEntries());
  TEST_EQUALITY(A->getGlobalMaxNumRowEntries(), Acopy->getGlobalMaxNumRowEntries());

  // FIXME (mfh 24 Apr 2014) Need to test separately on each MPI
  // process and do an all-reduce to check if all got it right.
  TEST_EQUALITY(A->getLocalNumRows(), Acopy->getLocalNumRows());
  TEST_EQUALITY(A->getLocalNumCols(), Acopy->getLocalNumCols());
  TEST_EQUALITY(A->getLocalNumEntries(), Acopy->getLocalNumEntries());
  TEST_EQUALITY(A->getLocalMaxNumRowEntries(), Acopy->getLocalMaxNumRowEntries());

  // Acopy and A should be identically the same.  We can verify this
  // in two ways.  First, we can directly compare the rows of both
  // matrices on each process.  Second, we can repeat the apply()
  // operation with Acopy and verify that it produces the same
  // result.  We will take both approaches here.

  // This test only makes sense if the row Maps of A and Acopy
  // match.  Otherwise, some row indices that are valid for one
  // matrix might not be valid for the other.
  if (rowMapsMatch) {
    typedef typename Teuchos::ArrayView<const GO>::size_type size_type;
    // mfh 24 Apr 2014: Apparently, Xpetra::CrsMatrix implements
    // neither getGlobalRowCopy nor getNumEntriesInGlobalRow.
    //
    // Teuchos::Array<GO> A_ginds (A->getLocalMaxNumRowEntries ());
    Teuchos::Array<LO> A_linds(A->getLocalMaxNumRowEntries());
    Teuchos::Array<Scalar> A_vals(A->getLocalMaxNumRowEntries());

    // mfh 24 Apr 2014: Apparently, Xpetra::CrsMatrix implements
    // neither getGlobalRowCopy nor getNumEntriesInGlobalRow.
    //
    // Teuchos::Array<GO> Acopy_ginds (Acopy->getLocalMaxNumRowEntries ());
    Teuchos::Array<LO> Acopy_linds(Acopy->getLocalMaxNumRowEntries());
    Teuchos::Array<Scalar> Acopy_vals(Acopy->getLocalMaxNumRowEntries());

    for (size_type k = 0; k < static_cast<size_type>(NumMyElements); ++k) {
      const LO lrow = static_cast<LO>(k);
      // const GO grow = MyGlobalElements[k];
      size_t A_numEnt     = 0;
      size_t Acopy_numEnt = 0;

      // mfh 24 Apr 2014: Apparently, Xpetra::CrsMatrix implements
      // neither getGlobalRowCopy nor getNumEntriesInGlobalRow.
      //
      // A->getGlobalRowCopy (grow, A_ginds (), A_vals (), A_numEnt);
      // Acopy->getGlobalRowCopy (grow, Acopy_ginds (), Acopy_vals (), Acopy_numEnt);
      // TEST_COMPARE_ARRAYS( A_ginds (0, A_numEnt), Acopy_ginds (0, Acopy_numEnt) );
      // TEST_COMPARE_ARRAYS( A_vals (0, A_numEnt), Acopy_vals (0, Acopy_numEnt) );

      TEST_EQUALITY(A->getNumEntriesInLocalRow(lrow), Acopy->getNumEntriesInLocalRow(lrow));

      // mfh 24 Apr 2014: Apparently, Xpetra::CrsMatrix implements
      // neither getGlobalRowCopy nor getNumEntriesInGlobalRow.
      //
      // TEST_EQUALITY(A->getNumEntriesInGlobalRow (grow), Acopy->getNumEntriesInGlobalRow (grow));

      A->getLocalRowCopy(lrow, A_linds(), A_vals(), A_numEnt);
      Acopy->getLocalRowCopy(lrow, Acopy_linds(), Acopy_vals(), Acopy_numEnt);

      TEST_COMPARE_ARRAYS(A_linds(0, A_numEnt), Acopy_linds(0, Acopy_numEnt));
      TEST_COMPARE_ARRAYS(A_vals(0, A_numEnt), Acopy_vals(0, Acopy_numEnt));
    }
  }

  // Make sure that A doesn't exist anymore.
  A = Teuchos::null;

  {
    using std::cerr;
    using std::endl;

    std::ostringstream os;
    const int myRank = comm->getRank();
    os << "Process " << myRank << ": Applying matrix Acopy" << endl;
    cerr << os.str();
  }

  Acopy->apply(*v, *rcopy, Teuchos::NO_TRANS, STS::one(), STS::zero());

  {
    // Since A is the identity matrix and Acopy is a deep copy of A,
    // after this line finishes, rcopy and v should be exactly
    // equal.  This should be true even in finite-precision
    // arithmetic.  Test this here.
    RCP<Xpetra::Vector<Scalar, LO, GO, Node> > diff =
        Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(Acopy->getRangeMap());

    // diff := 1*v + (-1)*rcopy.
    diff->update(STS::one(), *v, -STS::one(), *rcopy, STS::zero());
    Teuchos::Array<magnitude_type> norms(1);
    diff->norm2(norms());
    // The norm of v - r must be _exactly_ zero.
    TEST_EQUALITY(norms[0], STM::zero());
  }

  // Repeat the above test in a different way.
  Teuchos::ArrayRCP<const Scalar> rdata     = r->getData(0);
  Teuchos::ArrayRCP<const Scalar> rdatacopy = rcopy->getData(0);
  magnitude_type s                          = STM::zero();
  for (LO i = 0; i < NumMyElements; ++i) {
    s += Teuchos::ScalarTraits<Scalar>::magnitude(rdata[i] - rdatacopy[i]);
  }
  TEUCHOS_TEST_COMPARE(s, <, 1e-16, out, success);
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, EpetraDeepCopy, Scalar, LO, GO, Node) {
#ifdef HAVE_XPETRA_EPETRA

  typedef Xpetra::Map<LO, GO, Node> MapClass;
  typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;

  // get a comm and node
  RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;

  // generate problem
  LO nEle                       = 63;
  const RCP<const MapClass> map = MapFactoryClass::Build(lib, nEle, 0, comm);

  RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > A =
      Xpetra::CrsMatrixFactory<Scalar, LO, GO, Node>::Build(map, 10);

  LO NumMyElements                              = map->getLocalNumElements();
  Teuchos::ArrayView<const GO> MyGlobalElements = map->getLocalElementList();

  for (LO i = 0; i < NumMyElements; ++i) {
    A->insertGlobalValues(MyGlobalElements[i],
                          Teuchos::tuple<GO>(MyGlobalElements[i]),
                          Teuchos::tuple<Scalar>(1.0));
  }

  A->fillComplete();

  RCP<Xpetra::Vector<Scalar, LO, GO, Node> > v = Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(map);
  v->setSeed(8675309);
  v->randomize(true);

  RCP<Xpetra::Vector<Scalar, LO, GO, Node> > r     = Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(A->getRangeMap());
  RCP<Xpetra::Vector<Scalar, LO, GO, Node> > rcopy = Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(A->getRangeMap());

  A->apply(*v, *r, Teuchos::NO_TRANS, 1.0, 0.0);

  RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > Acopy(new Xpetra::EpetraCrsMatrixT<GO, Node>(*(Teuchos::rcp_static_cast<Xpetra::EpetraCrsMatrixT<GO, Node> >(A))));
  A = Teuchos::null;

  Acopy->apply(*v, *rcopy, Teuchos::NO_TRANS, 1.0, 0.0);

  Teuchos::ArrayRCP<Scalar> rdata = r->getDataNonConst(0), rdatacopy = rcopy->getDataNonConst(0);
  Scalar s = Teuchos::ScalarTraits<Scalar>::zero();
  for (LO i = 0; i < NumMyElements; i++) {
    s += Teuchos::ScalarTraits<Scalar>::magnitude(rdata[i] - rdatacopy[i]);
  }
  TEUCHOS_TEST_COMPARE(s, <, 1e-16, out, success);
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_5_DECL(CrsMatrix, GetLocalMatrix, M, Scalar, LO, GO, Node) {
  typedef Xpetra::Map<LO, GO, Node> MapClass;
  typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;
  typedef typename Xpetra::CrsMatrix<Scalar, LO, GO, Node>::local_matrix_type local_matrix_type;
  // typedef typename local_matrix_type::size_type size_type;
  typedef typename local_matrix_type::value_type value_type;
  typedef typename local_matrix_type::ordinal_type ordinal_type;
  using STS        = Teuchos::ScalarTraits<Scalar>;
  const Scalar ONE = STS::one();
  using mag_type   = typename STS::magnitudeType;

  // get a comm and node
  RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  M testMap(1, 0, comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  // generate problem
  LO nEle                       = 63;
  const RCP<const MapClass> map = MapFactoryClass::Build(lib, nEle, 0, comm);

  RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > A =
      Xpetra::CrsMatrixFactory<Scalar, LO, GO, Node>::Build(map, 10);

  LO NumMyElements                              = map->getLocalNumElements();
  Teuchos::ArrayView<const GO> MyGlobalElements = map->getLocalElementList();

  for (LO i = 0; i < NumMyElements; ++i) {
    A->insertGlobalValues(MyGlobalElements[i],
                          Teuchos::tuple<GO>(MyGlobalElements[i]),
                          Teuchos::tuple<Scalar>(1.0));
  }

  A->fillComplete();

  {
    // access data after fill complete!
    auto view2 = A->getLocalMatrixHost();
    TEST_EQUALITY(Teuchos::as<size_t>(view2.numRows()), A->getLocalNumRows());
    TEST_EQUALITY(Teuchos::as<size_t>(view2.numCols()), A->getLocalNumCols());
    TEST_EQUALITY(Teuchos::as<size_t>(view2.nnz()), A->getLocalNumEntries());

    // check that the local_matrix_type taken the second time is the same
    auto view3 = A->getLocalMatrixHost();
    // The row pointer is only identical for Tpetra. For Epetra, the
    // rowptr has a different type than what the local matrix wants,
    // so we are copying to a new Kokkos view..
    if (map->lib() == Xpetra::UseTpetra)
      TEST_EQUALITY(view2.graph.row_map.data(), view3.graph.row_map.data());
    TEST_EQUALITY(view2.graph.entries.data(), view3.graph.entries.data());
    TEST_EQUALITY(view2.values.data(), view3.values.data());

    for (LO r = 0; r < view2.numRows(); ++r) {
      // extract data from current row r
      auto rowview = view2.row(r);

      for (LO c = 0; c < rowview.length; ++c) {
        Scalar vv = rowview.value(c);
        LO cc     = rowview.colidx(c);
        TEST_EQUALITY(rowview.length, 1);
        TEST_EQUALITY(cc, r);
        TEST_EQUALITY(vv, ONE);
      }
    }

    Teuchos::ArrayView<const LO> indices;
    Teuchos::ArrayView<const Scalar> values;
    A->getLocalRowView(0, indices, values);
    TEST_EQUALITY(indices.size(), 1);
    TEST_EQUALITY(values[0], ONE);

    /////////////////////////////////////////

    // check whether later changes are updated in view!
    ordinal_type nColIdx = 0;
    value_type value     = 42.0;
    view2.replaceValues(0, &nColIdx, 1, &value);

    A->getLocalRowView(0, indices, values);
    TEST_EQUALITY(indices.size(), 1);

    // NOTE (mfh 23 Feb 2020) You can't always convert double to
    // Scalar directly; e.g., with Scalar=complex<float>.
    const Scalar FORTY_TWO = Scalar(mag_type(42.0));
    TEST_EQUALITY(values[0], FORTY_TWO);  // changes in the view also changes matrix values
  }

  A->setAllToScalar(-123.4);

  {
    auto view2 = A->getLocalMatrixHost();
    TEST_EQUALITY(Teuchos::as<size_t>(view2.numRows()), A->getLocalNumRows());
    TEST_EQUALITY(Teuchos::as<size_t>(view2.numCols()), A->getLocalNumCols());
    TEST_EQUALITY(Teuchos::as<size_t>(view2.nnz()), A->getLocalNumEntries());

    for (LO r = 0; r < view2.numRows(); ++r) {
      // extract data from current row r
      auto rowview = view2.row(r);

      for (LO c = 0; c < rowview.length; ++c) {
        Scalar vv = rowview.value(c);
        LO cc     = rowview.colidx(c);
        TEST_EQUALITY(rowview.length, 1);
        TEST_EQUALITY(cc, r);
        // NOTE (mfh 23 Feb 2020) You can't always convert double to
        // Scalar directly; e.g., with Scalar=complex<float>.
        const Scalar expected_vv = Scalar(mag_type(-123.4));
        TEST_EQUALITY(vv, expected_vv);
      }
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_5_DECL(CrsMatrix, ConstructMatrixKokkos, M, Scalar, LO, GO, Node) {
#ifdef HAVE_XPETRA_TPETRA  // Note: get Kokkos interface for Epetra is only available if Tpetra is also enabled!
  std::cout << "Run ConstructMatrixKokkos test" << std::endl;
  // Kokkos::initialize();
  typedef Xpetra::Map<LO, GO, Node> MapClass;
  typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;
  typedef typename Xpetra::CrsMatrix<Scalar, LO, GO, Node> CrsMatrixClass;
  typedef typename Xpetra::CrsMatrixFactory<Scalar, LO, GO, Node> CrsMatrixFactoryClass;
  typedef typename CrsMatrixClass::local_matrix_type local_matrix_type;

  // get a comm and node
  RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
  int rank                            = comm->getRank();
  int numProcs                        = comm->getSize();

  M testMap(1, 0, comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  // choose memory space of local_matrix_type (could be e.g. on the device)
  typedef typename local_matrix_type::memory_space MemorySpace;
  typedef typename local_matrix_type::size_type size_type;
  typedef typename local_matrix_type::value_type value_type;
  typedef typename local_matrix_type::ordinal_type ordinal_type;

  // typedefs for the data of the local CRS matrix
  typedef Kokkos::View<size_type *, MemorySpace> ptr_type;
  typedef Kokkos::View<ordinal_type *, MemorySpace> ind_type;
  typedef Kokkos::View<value_type *, MemorySpace> val_type;

  // build local Kokkos::CrsMatrix
  LO numRows, numCols, nnz;
  if (rank == 0 && numProcs == 1) { /* special case: only one proc */
    numRows = 3;
    numCols = 3;
    nnz     = 7;
  } else if (rank == 0 || rank == numProcs - 1) { /* first or last processor */
    numRows = 3;
    numCols = 4;
    nnz     = 8;
  } else { /* all other processors processor */
    numRows = 3;
    numCols = 5;
    nnz     = 9;
  }

  // Create the output Views.
  ptr_type ptr = ptr_type("ptr", numRows + 1);
  ind_type ind = ind_type("ind", nnz);
  val_type val = val_type("val", nnz);

  // build local Kokkos::CrsMatrix
  if (rank == 0 && numProcs == 1) { /* special case: only one processor */
    const size_type ptrRaw[]    = {0, 2, 5, 7};
    const ordinal_type indRaw[] = {0, 1,
                                   0, 1, 2,
                                   1, 2};
    const value_type valRaw[]   = {2.0, -1.0,
                                   -1.0, 2.0, -1.0,
                                   -1.0, 2.0};
    // Wrap the above three arrays in unmanaged Views, so we can use deep_copy.
    typename ptr_type::HostMirror::const_type ptrIn(ptrRaw, numRows + 1);
    typename ind_type::HostMirror::const_type indIn(indRaw, nnz);
    typename val_type::HostMirror::const_type valIn(valRaw, nnz);

    Kokkos::deep_copy(ptr, ptrIn);
    Kokkos::deep_copy(ind, indIn);
    Kokkos::deep_copy(val, valIn);
  } else if (rank == 0 && numProcs > 1) { /* first processor */
    const size_type ptrRaw[]    = {0, 2, 5, 8};
    const ordinal_type indRaw[] = {0, 1,
                                   0, 1, 2,
                                   1, 2, 3};
    const value_type valRaw[]   = {2.0, -1.0,
                                   -1.0, 2.0, -1.0,
                                   -1.0, 2.0, -1.0};
    // Wrap the above three arrays in unmanaged Views, so we can use deep_copy.
    typename ptr_type::HostMirror::const_type ptrIn(ptrRaw, numRows + 1);
    typename ind_type::HostMirror::const_type indIn(indRaw, nnz);
    typename val_type::HostMirror::const_type valIn(valRaw, nnz);

    Kokkos::deep_copy(ptr, ptrIn);
    Kokkos::deep_copy(ind, indIn);
    Kokkos::deep_copy(val, valIn);
  } else if (rank == numProcs - 1) { /* last processor */
    const size_type ptrRaw[]    = {0, 3, 6, 8};
    const ordinal_type indRaw[] = {3, 0, 1,
                                   0, 1, 2,
                                   1, 2};
    const value_type valRaw[]   = {-1.0, 2.0, -1.0,
                                   -1.0, 2.0, -1.0,
                                   -1.0, 2.0};
    // Wrap the above three arrays in unmanaged Views, so we can use deep_copy.
    typename ptr_type::HostMirror::const_type ptrIn(ptrRaw, numRows + 1);
    typename ind_type::HostMirror::const_type indIn(indRaw, nnz);
    typename val_type::HostMirror::const_type valIn(valRaw, nnz);

    Kokkos::deep_copy(ptr, ptrIn);
    Kokkos::deep_copy(ind, indIn);
    Kokkos::deep_copy(val, valIn);
  } else { /* all other processors processor */
    const size_type ptrRaw[]    = {0, 3, 6, 9};
    const ordinal_type indRaw[] = {3, 0, 1,
                                   0, 1, 2,
                                   1, 2, 4};
    const value_type valRaw[]   = {-1.0, 2.0, -1.0,
                                   -1.0, 2.0, -1.0,
                                   -1.0, 2.0, -1.0};
    // Wrap the above three arrays in unmanaged Views, so we can use deep_copy.
    typename ptr_type::HostMirror::const_type ptrIn(ptrRaw, numRows + 1);
    typename ind_type::HostMirror::const_type indIn(indRaw, nnz);
    typename val_type::HostMirror::const_type valIn(valRaw, nnz);

    Kokkos::deep_copy(ptr, ptrIn);
    Kokkos::deep_copy(ind, indIn);
    Kokkos::deep_copy(val, valIn);
  }

  // reconstruct row and column map
  // std::vector<GO> rowMapGids;  // vector for collecting row map GIDs
  // std::vector<GO> colMapGids;  // vector for collecting column map GIDs
  Teuchos::Array<GO> rowMapGids;  // vector for collecting row map GIDs
  Teuchos::Array<GO> colMapGids;  // vector for collecting column map GID
  if (rank == 0 && numProcs == 1) {
    rowMapGids.push_back(0);
    colMapGids.push_back(0);
    rowMapGids.push_back(1);
    colMapGids.push_back(1);
    rowMapGids.push_back(2);
    colMapGids.push_back(2);
  } else if (rank == 0 && numProcs > 1) { /* first processor */
    rowMapGids.push_back(0);
    colMapGids.push_back(0);
    rowMapGids.push_back(1);
    colMapGids.push_back(1);
    rowMapGids.push_back(2);
    colMapGids.push_back(2);
    colMapGids.push_back(3);
  } else if (rank == numProcs - 1) { /* last processor */
    rowMapGids.push_back(rank * 3 + 0);
    colMapGids.push_back(rank * 3 + 0);
    rowMapGids.push_back(rank * 3 + 1);
    colMapGids.push_back(rank * 3 + 1);
    rowMapGids.push_back(rank * 3 + 2);
    colMapGids.push_back(rank * 3 + 2);
    colMapGids.push_back(rank * 3 - 1);
  } else { /* all other processors */

    rowMapGids.push_back(rank * 3 + 0);
    colMapGids.push_back(rank * 3 + 0);
    rowMapGids.push_back(rank * 3 + 1);
    colMapGids.push_back(rank * 3 + 1);
    rowMapGids.push_back(rank * 3 + 2);
    colMapGids.push_back(rank * 3 + 2);
    colMapGids.push_back(rank * 3 - 1);
    colMapGids.push_back(rank * 3 + 3);
  }

  // create Xpetra::Map objects for row map and column map
  const GO INVALID = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();
  // Teuchos::ArrayView<GO> rowMapGidsView(&rowMapGids[0], rowMapGids.size());
  // Teuchos::ArrayView<GO> colMapGidsView(&colMapGids[0], colMapGids.size());
  Teuchos::RCP<const MapClass> rowMap = MapFactoryClass::Build(lib, INVALID, rowMapGids.view(0, rowMapGids.size()), 0, comm);
  Teuchos::RCP<const MapClass> colMap = MapFactoryClass::Build(lib, INVALID, colMapGids.view(0, colMapGids.size()), 0, comm);

  Teuchos::RCP<CrsMatrixClass> mat;
  {
    // create local CrsMatrix
    local_matrix_type lclMatrix = local_matrix_type("A", numRows, numCols, nnz, val, ptr, ind);

    mat = CrsMatrixFactoryClass::Build(rowMap, colMap, lclMatrix);
    if (mat == Teuchos::null) std::cout << "mat is Teuchos::null..." << std::endl;

    // Blank out views to reduce ref count
    ptr = ptr_type();
    ind = ind_type();
    val = val_type();
  }

  TEST_EQUALITY(mat->isFillComplete(), true);
  TEST_EQUALITY(mat->getGlobalMaxNumRowEntries(), 3);
  TEST_EQUALITY(mat->getLocalMaxNumRowEntries(), 3);
  TEST_EQUALITY(mat->getLocalNumRows(), 3);
  TEST_EQUALITY(static_cast<size_t>(mat->getGlobalNumCols()),
                static_cast<size_t>(3 * numProcs));
  TEST_EQUALITY(static_cast<size_t>(mat->getGlobalNumRows()),
                static_cast<size_t>(3 * numProcs));
  TEST_EQUALITY(static_cast<size_t>(mat->getGlobalNumEntries()),
                static_cast<size_t>(9 * numProcs - 2));

  size_t numLocalRows = mat->getLocalNumRows();
  for (size_t row = 0; row < numLocalRows; row++) {
    GO grid = mat->getRowMap()->getGlobalElement(row);
    // extract row information from input matrix
    Teuchos::ArrayView<const LO> indices;
    Teuchos::ArrayView<const Scalar> vals;
    mat->getLocalRowView(row, indices, vals);
    for (size_t col = 0; col < size_t(indices.size()); col++) {
      using STS        = Teuchos::ScalarTraits<Scalar>;
      const Scalar ONE = STS::one();
      const Scalar TWO = ONE + ONE;

      if (grid == colMap->getGlobalElement(indices[col])) {
        TEST_EQUALITY(vals[col], TWO);
      } else {
        TEST_EQUALITY(vals[col], -ONE);
      }
    }
  }
  // Kokkos::finalize();
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, TpetraGraphAndValuesConstructor, Scalar, LO, GO, Node) {
#ifdef HAVE_XPETRA_TPETRA
  typedef Xpetra::Map<LO, GO, Node> MapClass;
  typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename STS::magnitudeType magnitude_type;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;
  using Matrix = Xpetra::Matrix<Scalar, LO, GO, Node>;
  using IST    = typename Xpetra::CrsMatrix<Scalar, LO, GO, Node>::impl_scalar_type;

  // get a comm and node
  RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  {
    using std::cerr;
    using std::endl;

    std::ostringstream os;
    const int myRank = comm->getRank();
    os << "Process " << myRank << ": (CrsMatrix, TpetraDeepCopy) test" << endl;
    cerr << os.str();
  }

  Xpetra::UnderlyingLib lib = Xpetra::UseTpetra;

  // Create a Map, which will be the row, domain, and range Map of the matrix A.
  LO nEle                       = 63;
  const RCP<const MapClass> map = MapFactoryClass::Build(lib, nEle, 0, comm);

  {
    using std::cerr;
    using std::endl;

    std::ostringstream os;
    const int myRank = comm->getRank();
    os << "Process " << myRank << ": Creating matrix" << endl;
    cerr << os.str();
  }

  RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > A =
      Xpetra::CrsMatrixFactory<Scalar, LO, GO, Node>::Build(map, 10);

  LO NumMyElements                              = map->getLocalNumElements();
  Teuchos::ArrayView<const GO> MyGlobalElements = map->getLocalElementList();

  {
    using std::cerr;
    using std::endl;

    std::ostringstream os;
    const int myRank = comm->getRank();
    os << "Process " << myRank << ": Filling matrix" << endl;
    cerr << os.str();
  }

  // Make A the identity matrix.
  for (LO i = 0; i < NumMyElements; ++i) {
    A->insertGlobalValues(MyGlobalElements[i],
                          Teuchos::tuple<GO>(MyGlobalElements[i]),
                          Teuchos::tuple<Scalar>(1.0));
  }

  {
    using std::cerr;
    using std::endl;

    std::ostringstream os;
    const int myRank = comm->getRank();
    os << "Process " << myRank << ": Calling fillComplete on matrix A" << endl;
    cerr << os.str();
  }

  A->fillComplete();

  // Create new values array containing all ones.
  typename Xpetra::CrsMatrix<Scalar, LO, GO, Node>::local_matrix_type::values_type new_values("new_values", A->getLocalNumEntries());
  Kokkos::deep_copy(new_values, Teuchos::ScalarTraits<IST>::one());

  // Use the MatrixFactory to build a copy
  RCP<Matrix> Acopy = Xpetra::MatrixFactory<Scalar, LO, GO, Node>::Build(A->getCrsGraph(), new_values);
  Acopy->fillComplete(A->getDomainMap(), A->getRangeMap());

  {
    using std::cerr;
    using std::endl;

    std::ostringstream os;
    const int myRank = comm->getRank();
    os << "Process " << myRank << ": Building an input Vector" << endl;
    cerr << os.str();
  }

  RCP<Xpetra::Vector<Scalar, LO, GO, Node> > v = Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(A->getRangeMap());
  v->setSeed(8675309);
  v->randomize(true);

  // r and rcopy are distinct Vectors with the same Map, namely the
  // range Map of A.  All the Vectors v, r, and rcopy are distinct.
  RCP<Xpetra::Vector<Scalar, LO, GO, Node> > r =
      Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(A->getRangeMap());
  RCP<Xpetra::Vector<Scalar, LO, GO, Node> > rcopy =
      Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(A->getRangeMap());

  {
    using std::cerr;
    using std::endl;

    std::ostringstream os;
    const int myRank = comm->getRank();
    os << "Process " << myRank << ": Applying matrix A" << endl;
    cerr << os.str();
  }

  // r := A * v.
  A->apply(*v, *r, Teuchos::NO_TRANS, STS::one(), STS::zero());
  Acopy->apply(*v, *rcopy, Teuchos::NO_TRANS, STS::one(), STS::zero());

  // Since A is the identity matrix, after the above line finishes,
  // r and rcopy should be exactly equal.  This should be true even in
  // finite-precision arithmetic.  Test this here.
  {
    RCP<Xpetra::Vector<Scalar, LO, GO, Node> > diff =
        Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(A->getRangeMap());

    // diff := 1*rcopy + (-1)*r.
    diff->update(STS::one(), *rcopy, -STS::one(), *r, STS::zero());
    Teuchos::Array<magnitude_type> norms(1);
    diff->norm2(norms());
    // The norm of v - r must be _exactly_ zero.
    TEST_EQUALITY(norms[0], STM::zero());
  }

#endif
}

//
// INSTANTIATIONS
//

#ifdef HAVE_XPETRA_TPETRA

#define XPETRA_TPETRA_TYPES(SC, LO, GO, Node) \
  typedef typename Xpetra::TpetraMap<LO, GO, Node> M##LO##GO##Node;

#endif

#ifdef HAVE_XPETRA_EPETRA

#define XPETRA_EPETRA_TYPES(SC, LO, GO, Node) \
  typedef typename Xpetra::EpetraMapT<GO, Node> M##LO##GO##Node;

#endif

// for common tests (Epetra and Tpetra...)
#define UNIT_TEST_GROUP_ORDINAL(SC, LO, GO, Node)                                                                 \
  TEUCHOS_UNIT_TEST_TEMPLATE_5_INSTANT(CrsMatrix, Constructor_Vector, M##LO##GO##Node, SC, LO, GO, Node)          \
  TEUCHOS_UNIT_TEST_TEMPLATE_5_INSTANT(CrsMatrix, Apply, M##LO##GO##Node, SC, LO, GO, Node)                       \
  TEUCHOS_UNIT_TEST_TEMPLATE_5_INSTANT(CrsMatrix, ReplaceGlobalAndLocalValues, M##LO##GO##Node, SC, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_5_INSTANT(CrsMatrix, leftScale, M##LO##GO##Node, SC, LO, GO, Node)                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_5_INSTANT(CrsMatrix, rightScale, M##LO##GO##Node, SC, LO, GO, Node)                  \
  TEUCHOS_UNIT_TEST_TEMPLATE_5_INSTANT(CrsMatrix, replaceDiagonal, M##LO##GO##Node, SC, LO, GO, Node)
// for Tpetra tests only
#define UNIT_TEST_GROUP_ORDINAL_TPETRAONLY(SC, LO, GO, Node)                                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, TpetraDeepCopy, SC, LO, GO, Node)            \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, Tpetra_ReplaceLocalValues, SC, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, TpetraGraphAndValuesConstructor, SC, LO, GO, Node)
// for Epetra tests only
#define UNIT_TEST_GROUP_ORDINAL_EPETRAONLY(SC, LO, GO, Node)                                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, Constructor_Epetra, SC, LO, GO, Node)        \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, Epetra_ReplaceLocalValues, SC, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, EpetraDeepCopy, SC, LO, GO, Node)
// for Kokkos-specific tests
#define UNIT_TEST_GROUP_ORDINAL_KOKKOS(SC, LO, GO, Node)                                             \
  TEUCHOS_UNIT_TEST_TEMPLATE_5_INSTANT(CrsMatrix, GetLocalMatrix, M##LO##GO##Node, SC, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_5_INSTANT(CrsMatrix, ConstructMatrixKokkos, M##LO##GO##Node, SC, LO, GO, Node)

#if defined(HAVE_XPETRA_TPETRA)

#include <TpetraCore_config.h>
#include <TpetraCore_ETIHelperMacros.h>

TPETRA_ETI_MANGLING_TYPEDEFS()
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(XPETRA_TPETRA_TYPES)
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(UNIT_TEST_GROUP_ORDINAL)
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(UNIT_TEST_GROUP_ORDINAL_TPETRAONLY)
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(UNIT_TEST_GROUP_ORDINAL_KOKKOS)

#endif

#if defined(HAVE_XPETRA_EPETRA)
#include "Xpetra_Map.hpp"  // defines EpetraNode
typedef Xpetra::EpetraNode EpetraNode;
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
XPETRA_EPETRA_TYPES(double, int, int, EpetraNode)
UNIT_TEST_GROUP_ORDINAL_EPETRAONLY(double, int, int, EpetraNode)
UNIT_TEST_GROUP_ORDINAL_KOKKOS(double, int, int, EpetraNode)
#endif
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
typedef long long LongLong;
XPETRA_EPETRA_TYPES(double, int, LongLong, EpetraNode)
UNIT_TEST_GROUP_ORDINAL_EPETRAONLY(double, int, LongLong, EpetraNode)
UNIT_TEST_GROUP_ORDINAL_KOKKOS(double, int, LongLong, EpetraNode)
#endif
#endif

}  // namespace
