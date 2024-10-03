// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// This bug is described at
// https://github.com/trilinos/Trilinos/issues/8447
//
#include <Tpetra_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <TpetraCore_ETIHelperMacros.h>
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_TestingUtilities.hpp>

#include <iostream>
#include <vector>

namespace { // anonymous

using Teuchos::ArrayView;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Comm;
using Teuchos::OrdinalTraits;
using Tpetra::TestingUtilities::getDefaultComm;

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, Bug8447, SC, LO, GO, NT)
{

  typedef Tpetra::global_size_t GST; // Map's constructor needs this

  typedef Tpetra::Map<LO,GO,NT> MapType;
  typedef Tpetra::Import<LO,GO,NT> ImportType;
  typedef Tpetra::CrsMatrix<SC,LO,GO,NT> CrsMatrixType;
  typedef typename CrsMatrixType::impl_scalar_type implScalarType;
  typedef typename CrsMatrixType::local_matrix_device_type lclMatrixType;

  RCP<const Comm<int> > comm = getDefaultComm();
  TEUCHOS_TEST_FOR_EXCEPTION(
    comm->getSize () != 2, std::logic_error,
    "This test must be run with exactly 2 MPI processes.");

  const int myPID = comm->getRank();

  std::vector<GO> globalIDs, globalIDs2;
  std::vector<LO> rowptr;
  std::vector<LO> indices;
  std::vector<implScalarType> values;
  const SC ONE = Teuchos::ScalarTraits<implScalarType>::one ();
  if (myPID == 0) {
    globalIDs.push_back(0);
    globalIDs.push_back(1);

    globalIDs2.push_back(0);
    globalIDs2.push_back(1);
    globalIDs2.push_back(2);
    globalIDs2.push_back(3);

    rowptr.push_back(0);
    rowptr.push_back(2);
    rowptr.push_back(4);

    indices.push_back(0);
    indices.push_back(1);
    indices.push_back(0);
    indices.push_back(1);

    values.push_back(ONE);
    values.push_back(ONE);
    values.push_back(ONE);
    values.push_back(ONE);
  }
  else {
    globalIDs.push_back(2);
    globalIDs.push_back(3);

    rowptr.push_back(0);
    rowptr.push_back(2);
    rowptr.push_back(4);

    indices.push_back(0);
    indices.push_back(1);
    indices.push_back(0);
    indices.push_back(1);

    values.push_back(ONE);
    values.push_back(ONE);
    values.push_back(ONE);
    values.push_back(ONE);
  }

  const GST INVALID = OrdinalTraits<GST>::invalid ();
  const GO indexBase = 0;
  RCP<const MapType> rowMap =
    rcp (new MapType (INVALID, ArrayView<GO> (globalIDs), indexBase, comm));
  RCP<const MapType> rowMap2 =
    rcp (new MapType (INVALID, ArrayView<GO> (globalIDs2), indexBase, comm));
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! rowMap->isOneToOne (), std::logic_error,
    "In this test, the row Map is supposed to be one to one.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! rowMap2->isOneToOne (), std::logic_error,
    "In this test, the row Map is supposed to be one to one.");


  lclMatrixType lclMatrix = lclMatrixType("lclA", rowptr.size()-1, rowptr.size()-1, values.size(), values.data(), rowptr.data(), indices.data());

  RCP<CrsMatrixType> matrix = rcp (new CrsMatrixType (lclMatrix, rowMap, rowMap));

  RCP<const ImportType> importer = rcp(new ImportType (rowMap, rowMap2));
  RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList());

  // This line used to fail.
  RCP<CrsMatrixType> matrix2 = Tpetra::importAndFillCompleteCrsMatrix<CrsMatrixType>(matrix, *importer, *importer, rowMap2, rowMap2, params);

}

#define UNIT_TEST_GROUP_SC_LO_GO_NO( SC, LO, GO, NT )                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, Bug8447, SC, LO, GO, NT)

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(UNIT_TEST_GROUP_SC_LO_GO_NO)

} // anonymous
