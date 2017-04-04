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
// https://software.sandia.gov/bugzilla/show_bug.cgi?id=6171
//
// The relevant description of the problem this test exposed is:
//
// Vectors x and y start out containing 1's in all their rows. After a
//
//     CrsMatrix->apply(x, y, Teuchos::TRANS, alpha, beta)
//
// with alpha = 0 and beta = 1, some of the entries in y become zero when
// none of these entries should have changed for the stated values of alpha and
// beta. I've attached an example to demonstrate the problem.
//
// Some Macro Magic to ensure that if CUDA and KokkosCompat is enabled
// only the .cu version of this file is actually compiled
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

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, Bug6171, SC, LO, GO, NT)
{

  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType MagnitudeType;
  typedef Tpetra::global_size_t GST; // Map's constructor needs this

  typedef Tpetra::Map<LO,GO,NT> MapType;
  typedef Tpetra::CrsMatrix<SC,LO,GO,NT> CrsMatrixType;
  typedef Tpetra::Vector<SC,LO,GO,NT> VectorType;

  RCP<const Comm<int> > comm = getDefaultComm();
  TEUCHOS_TEST_FOR_EXCEPTION(
    comm->getSize () != 2, std::logic_error,
    "This test must be run with exactly 2 MPI processes.");

  const int myPID = comm->getRank();

  std::vector<GO> globalIDs;
  std::vector< std::vector<GO> > indices;
  std::vector< std::vector<SC> > values;
  if (myPID == 0) {
    globalIDs.push_back(0);
    globalIDs.push_back(1);
    globalIDs.push_back(2);
    values.resize(3);
    values[0].resize(2, 1);
    values[1].resize(2, 1);
    values[2].resize(2, 1);
    indices.resize(3);
    indices[0].push_back(0); indices[0].push_back(4);
    indices[1].push_back(0); indices[1].push_back(1);
    indices[2].push_back(1); indices[2].push_back(3);
  }
  else {
    globalIDs.push_back(3);
    globalIDs.push_back(4);
    values.resize(2);
    values[0].resize(2, 1);
    values[1].resize(2, 1);
    indices.resize(2);
    indices[0].push_back(2); indices[0].push_back(3);
    indices[1].push_back(2); indices[1].push_back(4);
  }

  const GST INVALID = OrdinalTraits<GST>::invalid ();
  const GO indexBase = 0;
  RCP<const MapType> rowMap =
    rcp (new MapType (INVALID, ArrayView<GO> (globalIDs), indexBase, comm));
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! rowMap->isOneToOne (), std::logic_error,
    "In this test, the row Map is supposed to be one to one.");

  RCP<CrsMatrixType> matrix = rcp (new CrsMatrixType (rowMap, 0));
  for (size_t i = 0; i < static_cast<size_t> (globalIDs.size ()); ++i) {
    matrix->insertGlobalValues (globalIDs[i],
                                ArrayView<const GO> (indices[i]),
                                ArrayView<const SC> (values[i]));
  }
  matrix->fillComplete ();

  VectorType x(rowMap), y(rowMap);
  x.putScalar (1.0);
  y.putScalar (1.0);
  const MagnitudeType normBefore = y.norm2 ();
  if (myPID == 0) {
    std::cout << "norm of y before matrix->apply = " << normBefore << std::endl;
  }
  SC alpha (0.0), beta (1.0);
  matrix->apply (x, y, Teuchos::TRANS, alpha, beta);
  const MagnitudeType normAfter = y.norm2 ();
  if (myPID == 0) {
    std::cout << "norm of y after matrix->apply = " << normAfter << std::endl;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
    normBefore != normAfter, std::logic_error,
    "normBefore = " << normBefore << " != normAfter = " << normAfter << ".");

}

#define UNIT_TEST_GROUP_SC_LO_GO_NO( SC, LO, GO, NT )                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, Bug6171, SC, LO, GO, NT)

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(UNIT_TEST_GROUP_SC_LO_GO_NO)

} // anonymous
