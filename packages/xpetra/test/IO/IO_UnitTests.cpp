// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_as.hpp>
#include <Xpetra_UnitTestHelpers.hpp>
#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_DefaultPlatform.hpp"
#include <Xpetra_IO.hpp>

namespace {

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(IO, MMMissingRows, M, MA, Scalar, LO, GO, Node) {
  using Teuchos::as;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
  TEUCHOS_ASSERT_EQUALITY(comm->getSize(), 1);

  if (Teuchos::ScalarTraits<Scalar>::isComplex)
    return;

  M testMap(1, 0, comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  auto A = Xpetra::IO<Scalar, LO, GO, Node>::Read("test.mtx", lib, comm, false);
  TEUCHOS_ASSERT_EQUALITY(A->getGlobalNumRows(), 5);
  TEUCHOS_ASSERT_EQUALITY(A->getGlobalNumCols(), 5);
  TEUCHOS_ASSERT_EQUALITY(A->getGlobalNumEntries(), 3);

  auto colmap = A->getColMap();
  auto crsA   = Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrixWrap<Scalar, LO, GO, Node> >(A, true)->getCrsMatrix();
  Teuchos::ArrayView<const LO> indices;
  Teuchos::ArrayView<const Scalar> values;
  crsA->getLocalRowView(0, indices, values);
  TEST_EQUALITY(indices.size(), 2);
  TEST_EQUALITY(colmap->getGlobalElement(indices[0]), 0);
  TEST_EQUALITY(colmap->getGlobalElement(indices[1]), 3);
  TEST_EQUALITY(values[0], as<Scalar>(2.));
  TEST_EQUALITY(values[1], as<Scalar>(3.));

  crsA->getLocalRowView(1, indices, values);
  TEST_EQUALITY(indices.size(), 1);
  TEST_EQUALITY(colmap->getGlobalElement(indices[0]), 4);
  TEST_EQUALITY(values[0], as<Scalar>(4.));
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(IO, BinaryMissingRows, M, MA, Scalar, LO, GO, Node) {
  using Teuchos::as;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
  TEUCHOS_ASSERT_EQUALITY(comm->getSize(), 1);

  M testMap(1, 0, comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  auto A = Xpetra::IO<Scalar, LO, GO, Node>::Read("test.mtx.bin", lib, comm, true);
  TEUCHOS_ASSERT_EQUALITY(A->getGlobalNumRows(), 5);
  TEUCHOS_ASSERT_EQUALITY(A->getGlobalNumCols(), 5);
  TEUCHOS_ASSERT_EQUALITY(A->getGlobalNumEntries(), 3);

  auto colmap = A->getColMap();
  auto crsA   = Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrixWrap<Scalar, LO, GO, Node> >(A, true)->getCrsMatrix();
  Teuchos::ArrayView<const LO> indices;
  Teuchos::ArrayView<const Scalar> values;
  crsA->getLocalRowView(0, indices, values);
  TEST_EQUALITY(indices.size(), 2);
  TEST_EQUALITY(colmap->getGlobalElement(indices[0]), 0);
  TEST_EQUALITY(colmap->getGlobalElement(indices[1]), 3);
  TEST_EQUALITY(values[0], as<Scalar>(2.));
  TEST_EQUALITY(values[1], as<Scalar>(3.));

  crsA->getLocalRowView(1, indices, values);
  TEST_EQUALITY(indices.size(), 1);
  TEST_EQUALITY(colmap->getGlobalElement(indices[0]), 4);
  TEST_EQUALITY(values[0], as<Scalar>(4.));
}

//
// INSTANTIATIONS
//
#ifdef HAVE_XPETRA_TPETRA

#define XPETRA_TPETRA_TYPES(S, LO, GO, N)                     \
  typedef typename Xpetra::TpetraMap<LO, GO, N> M##LO##GO##N; \
  typedef typename Xpetra::TpetraCrsMatrix<S, LO, GO, N> MA##S##LO##GO##N;

#endif

#ifdef HAVE_XPETRA_EPETRA

#define XPETRA_EPETRA_TYPES(S, LO, GO, N)                  \
  typedef typename Xpetra::EpetraMapT<GO, N> M##LO##GO##N; \
  typedef typename Xpetra::EpetraCrsMatrixT<GO, N> MA##S##LO##GO##N;

#endif

// list of all tests which run both with Epetra and Tpetra
#define XP_IO_INSTANT(S, LO, GO, N)                                                                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(IO, MMMissingRows, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N) \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(IO, BinaryMissingRows, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)

#if defined(HAVE_XPETRA_TPETRA)

#include <TpetraCore_config.h>
#include <TpetraCore_ETIHelperMacros.h>

TPETRA_ETI_MANGLING_TYPEDEFS()
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(XPETRA_TPETRA_TYPES)
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(XP_IO_INSTANT)

#endif

#if defined(HAVE_XPETRA_EPETRA)

#include "Xpetra_Map.hpp"  // defines EpetraNode
typedef Xpetra::EpetraNode EpetraNode;
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
XPETRA_EPETRA_TYPES(double, int, int, EpetraNode)
XP_IO_INSTANT(double, int, int, EpetraNode)
#endif
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
typedef long long LongLong;
XPETRA_EPETRA_TYPES(double, int, LongLong, EpetraNode)
XP_IO_INSTANT(double, int, LongLong, EpetraNode)
#endif

#endif

}  // namespace
