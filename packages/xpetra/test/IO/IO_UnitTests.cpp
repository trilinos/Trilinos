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
#include <Tpetra_BinaryIO.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <cstdio>
#include <sstream>

namespace {

std::string makeBinaryFilename(const std::string& stem,
                               const Teuchos::Comm<int>& comm) {
  std::ostringstream os;
  os << stem << "_p" << comm.getSize() << ".bin";
  return os.str();
}

void cleanupFile(const std::string& filename,
                 const Teuchos::RCP<const Teuchos::Comm<int> >& comm) {
  comm->barrier();
  if (comm->getRank() == 0) {
    std::remove(filename.c_str());
  }
  comm->barrier();
}

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
  Xpetra::UnderlyingLib lib  = testMap.lib();
  const std::string filename = makeBinaryFilename("xpetra_io_binary_missing_rows", *comm);

  using tpetra_map_type    = Tpetra::Map<LO, GO, Node>;
  using tpetra_matrix_type = Tpetra::CrsMatrix<Scalar, LO, GO, Node>;
  using binary_io_type     = Tpetra::BinaryIO<Scalar, LO, GO, Node>;

  auto tpetraRowMap = Teuchos::rcp(new tpetra_map_type(5, static_cast<GO>(0), comm));
  auto tpetraAWrite = Teuchos::rcp(new tpetra_matrix_type(tpetraRowMap, 2));
  Teuchos::Array<GO> cols;
  Teuchos::Array<Scalar> vals;
  cols.push_back(static_cast<GO>(0));
  cols.push_back(static_cast<GO>(3));
  vals.push_back(as<Scalar>(2.));
  vals.push_back(as<Scalar>(3.));
  tpetraAWrite->insertGlobalValues(static_cast<GO>(0), cols(), vals());
  cols.resize(1);
  vals.resize(1);
  cols[0] = static_cast<GO>(4);
  vals[0] = as<Scalar>(4.);
  tpetraAWrite->insertGlobalValues(static_cast<GO>(1), cols(), vals());
  tpetraAWrite->fillComplete(tpetraRowMap, tpetraRowMap);
  binary_io_type::writeSparseFile(filename, *tpetraAWrite);

  auto A = Xpetra::IO<Scalar, LO, GO, Node>::Read(filename, lib, comm, true);
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

  cleanupFile(filename, comm);
}

//
// INSTANTIATIONS
//

#define XPETRA_TPETRA_TYPES(S, LO, GO, N)                     \
  typedef typename Xpetra::TpetraMap<LO, GO, N> M##LO##GO##N; \
  typedef typename Xpetra::TpetraCrsMatrix<S, LO, GO, N> MA##S##LO##GO##N;

// list of all tests which run both with Epetra and Tpetra
#define XP_IO_INSTANT(S, LO, GO, N)                                                                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(IO, MMMissingRows, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N) \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(IO, BinaryMissingRows, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)

#include <TpetraCore_config.h>
#include <TpetraCore_ETIHelperMacros.h>

TPETRA_ETI_MANGLING_TYPEDEFS()
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(XPETRA_TPETRA_TYPES)
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(XP_IO_INSTANT)

}  // namespace
