// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <TpetraCore_ETIHelperMacros.h>
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_TestingUtilities.hpp>
#include <Tpetra_Filter.hpp>
#include <Teuchos_ScalarTraits.hpp>

namespace {  // anonymous

using Teuchos::ArrayView;
using Teuchos::Comm;
using Teuchos::OrdinalTraits;
using Teuchos::RCP;
using Teuchos::rcp;
using Tpetra::TestingUtilities::getDefaultComm;

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, Filter, SC, LO, GO, NT) {
  using MapType       = Tpetra::Map<LO, GO, NT>;
  using CrsMatrixType = Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using magnitudeType = typename Teuchos::ScalarTraits<SC>::magnitudeType;

  const auto ONE = Teuchos::ScalarTraits<magnitudeType>::one();

#if KOKKOS_VERSION >= 40799
  using ATS      = KokkosKernels::ArithTraits<SC>;
  using impl_SC  = typename ATS::val_type;
  using impl_ATS = KokkosKernels::ArithTraits<impl_SC>;
#else
  using ATS      = Kokkos::ArithTraits<SC>;
  using impl_SC  = typename ATS::val_type;
  using impl_ATS = Kokkos::ArithTraits<impl_SC>;
#endif

  RCP<const Comm<int> > comm = getDefaultComm();

  const int numProcs = comm->getSize();

  const size_t lclNumRows                = 4;
  const Tpetra::global_size_t gblNumRows = lclNumRows * numProcs;
  const GO indexBase                     = 0;
  auto rowMap                            = rcp(new MapType(gblNumRows, lclNumRows, indexBase, comm));

  CrsMatrixType matrix(rowMap, lclNumRows);

  if (rowMap->getLocalNumElements() != 0) {
    LO k = 1;
    for (LO lclRow = rowMap->getMinLocalIndex();
         lclRow <= rowMap->getMaxLocalIndex(); ++lclRow) {
      const GO gblRow = rowMap->getGlobalElement(lclRow);
      matrix.insertGlobalValues(gblRow,
                                Teuchos::tuple<GO>(gblRow),
                                Teuchos::tuple<SC>(k * ONE));
      ++k;
    }
  }
  matrix.fillComplete();

  TEST_EQUALITY(matrix.getLocalNumEntries(), lclNumRows);
  TEST_EQUALITY(matrix.getGlobalNumEntries(), numProcs * lclNumRows);

  auto filteredMatrix = Tpetra::applyFilter_vals(
      matrix, KOKKOS_LAMBDA(typename CrsMatrixType::impl_scalar_type val) {
        return impl_ATS::magnitude(val) <= 2.0;
      });

  TEST_EQUALITY(filteredMatrix->getLocalNumEntries(), lclNumRows - 2);
  TEST_EQUALITY(filteredMatrix->getGlobalNumEntries(), numProcs * (lclNumRows - 2));
}

#define UNIT_TEST_GROUP_SC_LO_GO_NO(SC, LO, GO, NT) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, Filter, SC, LO, GO, NT)

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(UNIT_TEST_GROUP_SC_LO_GO_NO)

}  // namespace
