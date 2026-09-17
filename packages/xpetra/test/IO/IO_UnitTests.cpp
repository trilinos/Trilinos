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
#include <fstream>
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

void writeLegacyBinaryMissingRowsFile(const std::string& filename,
                                      const Teuchos::RCP<const Teuchos::Comm<int> >& comm) {
  comm->barrier();
  if (comm->getRank() == 0) {
    std::ofstream out(filename.c_str(), std::ios::binary | std::ios::trunc);
    TEUCHOS_TEST_FOR_EXCEPTION(!out.good(), std::runtime_error, "Failed to open " << filename << " for writing.");

    const int m   = 5;
    const int n   = 5;
    const int nnz = 3;
    out.write(reinterpret_cast<const char*>(&m), sizeof(m));
    out.write(reinterpret_cast<const char*>(&n), sizeof(n));
    out.write(reinterpret_cast<const char*>(&nnz), sizeof(nnz));

    for (int row = 0; row < m; ++row) {
      const int rownnz = (row == 0 ? 2 : (row == 1 ? 1 : 0));
      out.write(reinterpret_cast<const char*>(&row), sizeof(row));
      out.write(reinterpret_cast<const char*>(&rownnz), sizeof(rownnz));
      if (row == 0) {
        const int columns[2] = {0, 3};
        const double values[2] = {2.0, 3.0};
        out.write(reinterpret_cast<const char*>(columns), sizeof(columns));
        out.write(reinterpret_cast<const char*>(values), sizeof(values));
      } else if (row == 1) {
        const int column = 4;
        const double value = 4.0;
        out.write(reinterpret_cast<const char*>(&column), sizeof(column));
        out.write(reinterpret_cast<const char*>(&value), sizeof(value));
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION(!out.good(), std::runtime_error, "Failed to write " << filename << ".");
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

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(IO, BinaryLegacyConversion, M, MA, Scalar, LO, GO, Node) {
  using Teuchos::as;

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
  TEUCHOS_ASSERT_EQUALITY(comm->getSize(), 1);

  M testMap(1, 0, comm);
  Xpetra::UnderlyingLib lib = testMap.lib();
  const std::string legacyFilename = makeBinaryFilename("xpetra_io_legacy_missing_rows", *comm);
  const std::string convertedFilename = makeBinaryFilename("xpetra_io_legacy_missing_rows_converted", *comm);

  writeLegacyBinaryMissingRowsFile(legacyFilename, comm);
  Xpetra::IO<Scalar, LO, GO, Node>::ConvertLegacyBinaryToBinary(legacyFilename, convertedFilename, lib, comm);

  auto A = Xpetra::IO<Scalar, LO, GO, Node>::Read(convertedFilename, lib, comm, true);
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

  cleanupFile(legacyFilename, comm);
  cleanupFile(convertedFilename, comm);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(IO, BinaryLegacyFallback, M, MA, Scalar, LO, GO, Node) {
  using Teuchos::as;

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
  TEUCHOS_ASSERT_EQUALITY(comm->getSize(), 1);

  M testMap(1, 0, comm);
  Xpetra::UnderlyingLib lib = testMap.lib();
  const std::string legacyFilename = makeBinaryFilename("xpetra_io_legacy_fallback", *comm);

  writeLegacyBinaryMissingRowsFile(legacyFilename, comm);

  auto A = Xpetra::IO<Scalar, LO, GO, Node>::Read(legacyFilename, lib, comm, true);
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

  cleanupFile(legacyFilename, comm);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(IO, BinaryNonLegacyFallbackRejects, M, MA, Scalar, LO, GO, Node) {
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
  TEUCHOS_ASSERT_EQUALITY(comm->getSize(), 1);

  M testMap(1, 0, comm);
  Xpetra::UnderlyingLib lib = testMap.lib();
  using io_type = Xpetra::IO<Scalar, LO, GO, Node>;

  TEST_THROW(io_type::Read("test.mtx", lib, comm, true), std::exception);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(IO, BinaryCustomColMap, M, MA, Scalar, LO, GO, Node) {
  using Teuchos::as;

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
  TEUCHOS_ASSERT_EQUALITY(comm->getSize(), 1);

  const std::string filename = makeBinaryFilename("xpetra_io_binary_custom_colmap", *comm);

  using tpetra_map_type    = Tpetra::Map<LO, GO, Node>;
  using tpetra_matrix_type = Tpetra::CrsMatrix<Scalar, LO, GO, Node>;
  using binary_io_type     = Tpetra::BinaryIO<Scalar, LO, GO, Node>;

  auto tpetraRowMap = Teuchos::rcp(new tpetra_map_type(5, static_cast<GO>(0), comm));
  Teuchos::Array<GO> colGids;
  colGids.push_back(static_cast<GO>(0));
  colGids.push_back(static_cast<GO>(3));
  colGids.push_back(static_cast<GO>(4));
  auto tpetraColMap = Teuchos::rcp(new tpetra_map_type(5, colGids(), static_cast<GO>(0), comm));
  auto tpetraAWrite = Teuchos::rcp(new tpetra_matrix_type(tpetraRowMap, tpetraColMap, 2));
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

  Teuchos::RCP<const Xpetra::Map<LO, GO, Node> > rowMap = Teuchos::rcp(new Xpetra::TpetraMap<LO, GO, Node>(tpetraRowMap));
  Teuchos::RCP<const Xpetra::Map<LO, GO, Node> > colMap = Teuchos::rcp(new Xpetra::TpetraMap<LO, GO, Node>(tpetraColMap));
  auto A                                                = Xpetra::IO<Scalar, LO, GO, Node>::Read(filename, rowMap, colMap, rowMap, rowMap, true, true);

  auto xpetraColMap = A->getColMap();
  TEST_ASSERT(colMap->isSameAs(*xpetraColMap));

  auto crsA = Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrixWrap<Scalar, LO, GO, Node> >(A, true)->getCrsMatrix();
  Teuchos::ArrayView<const LO> indices;
  Teuchos::ArrayView<const Scalar> values;
  crsA->getLocalRowView(0, indices, values);
  TEST_EQUALITY(indices.size(), 2);
  TEST_EQUALITY(xpetraColMap->getGlobalElement(indices[0]), 0);
  TEST_EQUALITY(xpetraColMap->getGlobalElement(indices[1]), 3);
  TEST_EQUALITY(values[0], as<Scalar>(2.));
  TEST_EQUALITY(values[1], as<Scalar>(3.));

  crsA->getLocalRowView(1, indices, values);
  TEST_EQUALITY(indices.size(), 1);
  TEST_EQUALITY(xpetraColMap->getGlobalElement(indices[0]), 4);
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
#define XP_IO_INSTANT(S, LO, GO, N)                                                                             \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(IO, MMMissingRows, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)         \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(IO, BinaryMissingRows, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)     \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(IO, BinaryLegacyConversion, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N) \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(IO, BinaryLegacyFallback, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)   \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(IO, BinaryNonLegacyFallbackRejects, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N) \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(IO, BinaryCustomColMap, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)

#include <TpetraCore_config.h>
#include <TpetraCore_ETIHelperMacros.h>

TPETRA_ETI_MANGLING_TYPEDEFS()
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(XPETRA_TPETRA_TYPES)
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(XP_IO_INSTANT)

}  // namespace
