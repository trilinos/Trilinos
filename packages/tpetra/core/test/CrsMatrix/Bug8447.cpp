/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

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
  typedef typename CrsMatrixType::local_matrix_type lclMatrixType;

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
