// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_as.hpp>
#include <Xpetra_UnitTestHelpers.hpp>
#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_DefaultPlatform.hpp"
#include <Xpetra_IO.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>

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

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IO, LegacyBinaryFallbackDistributed, Scalar, LO, GO, Node) {
  using Teuchos::as;

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
  Xpetra::UnderlyingLib lib = Xpetra::UseTpetra;
  const std::string filename = makeBinaryFilename("xpetra_io_legacy_fallback_distributed", *comm);

  writeLegacyBinaryMissingRowsFile(filename, comm);

  auto A = Xpetra::IO<Scalar, LO, GO, Node>::Read(filename, lib, comm, true);
  TEUCHOS_ASSERT_EQUALITY(A->getGlobalNumRows(), 5);
  TEUCHOS_ASSERT_EQUALITY(A->getGlobalNumCols(), 5);
  TEUCHOS_ASSERT_EQUALITY(A->getGlobalNumEntries(), 3);

  auto colmap = A->getColMap();
  auto crsA   = Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrixWrap<Scalar, LO, GO, Node> >(A, true)->getCrsMatrix();
  Teuchos::ArrayView<const LO> indices;
  Teuchos::ArrayView<const Scalar> values;

  const auto rowMap = A->getRowMap();
  for (size_t lclRow = 0; lclRow < rowMap->getLocalNumElements(); ++lclRow) {
    const GO gblRow = rowMap->getGlobalElement(static_cast<LO>(lclRow));
    crsA->getLocalRowView(static_cast<LO>(lclRow), indices, values);
    if (gblRow == static_cast<GO>(0)) {
      TEST_EQUALITY(indices.size(), 2);
      TEST_EQUALITY(colmap->getGlobalElement(indices[0]), 0);
      TEST_EQUALITY(colmap->getGlobalElement(indices[1]), 3);
      TEST_EQUALITY(values[0], as<Scalar>(2.));
      TEST_EQUALITY(values[1], as<Scalar>(3.));
    } else if (gblRow == static_cast<GO>(1)) {
      TEST_EQUALITY(indices.size(), 1);
      TEST_EQUALITY(colmap->getGlobalElement(indices[0]), 4);
      TEST_EQUALITY(values[0], as<Scalar>(4.));
    } else {
      TEST_EQUALITY(indices.size(), 0);
    }
  }

  cleanupFile(filename, comm);
}

#define XP_IO_LEGACY_INSTANT(S, LO, GO, N) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IO, LegacyBinaryFallbackDistributed, S, LO, GO, N)

#include <TpetraCore_config.h>
#include <TpetraCore_ETIHelperMacros.h>

TPETRA_ETI_MANGLING_TYPEDEFS()
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(XP_IO_LEGACY_INSTANT)

}  // namespace
