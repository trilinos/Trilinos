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

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"

// Test construction of two matrices using the same local matrix but 
// different domain and range maps; ensure they work the same way wrt apply

namespace {

template <typename map_t>
Teuchos::RCP<const map_t> getCyclicMap(
  size_t nIndices, 
  int mapNumProc, 
  const Teuchos::RCP<const Teuchos::Comm<int> > &comm)
{
  // Return a map that is cyclic (like dealing rows to processors)
  Teuchos::Array<typename map_t::global_ordinal_type> indices(nIndices);
  size_t cnt = 0;
  int me = comm->getRank();
  int np = comm->getSize();
  if (mapNumProc > np) mapNumProc = np; // corner case: bad input
  if (mapNumProc <= 0) mapNumProc = 1;  // corner case: np is too small

  for (size_t i = 0; i < nIndices; i++) 
    if (me == int(i % np)) indices[cnt++] = i;

  Tpetra::global_size_t dummy =
          Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();

  return rcp(new map_t(dummy, indices(0,cnt), 0, comm));
}

template <typename vector_t>
void initVector(vector_t &y)
{
  // initialize vector entries to their global element number
  auto data = y.getLocalViewHost(Tpetra::Access::OverwriteAll);
  for (size_t i = 0; i < y.getMap()->getNodeNumElements(); i++)
    data(i, 0) = y.getMap()->getGlobalElement(i);
}

///////////////////////////////////////////////////////////////////////////////

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, Bug9391, SC, LO, GO, NT)
{
  // Generate a matrix (ABlock) with default range and domain maps
  // Construct identical matrix (ACyclic) with cyclic range and domain maps
  // Do SpMV with both; compare norms of results

  using map_t = Tpetra::Map<>;
  using matrix_t = Tpetra::CrsMatrix<SC>;
  using vector_t = Tpetra::Vector<SC>;

  auto comm = Tpetra::getDefaultComm();
  int np = comm->getSize();

  // Build ABlock matrix -- a simple matrix with block rowMap and
  // domainMap == rangeMap == rowMap
  int nnzPerRow = 5;
  int nGlobalRowsAndCols = 100;
  Teuchos::RCP<map_t> blockMap = rcp(new map_t(nGlobalRowsAndCols, 0, comm));

  GO row;
  Teuchos::Array<GO> cols(nnzPerRow);
  Teuchos::Array<SC> vals(nnzPerRow);

  matrix_t ABlock(blockMap, nnzPerRow);
  for (size_t i = 0; i < blockMap->getNodeNumElements(); i++) {
    row = blockMap->getGlobalElement(i);
    int cnt = 0;
    for (int j = -2; j <= 2; j++) {
      if ((row+j >= 0) && (row+j < nGlobalRowsAndCols)) {
        cols[cnt] = row+j;
        vals[cnt] = row;
        cnt++;
      }
    }
    ABlock.insertGlobalValues(row, cols(0, cnt), vals(0, cnt));
  }
  ABlock.fillComplete();

  // Make ACyclic:  same matrix with different Domain and Range maps
  // Use cyclic domain and range maps
  // To make domain and range maps differ for square matrices,
  // keep some processors empty in the cyclic maps

  Teuchos::RCP<const map_t> domainMapCyclic = 
               getCyclicMap<map_t>(nGlobalRowsAndCols, np-1, comm);
  Teuchos::RCP<const map_t> rangeMapCyclic =
               getCyclicMap<map_t>(nGlobalRowsAndCols, np-2, comm);

  auto lclMatrix = ABlock.getLocalMatrixHost();
  matrix_t ACyclic(ABlock.getLocalMatrixHost(),
                   ABlock.getRowMap(), ABlock.getColMap(),
                   domainMapCyclic, rangeMapCyclic);

  // Create vectors and run SpMV for each matrix; compare norms
  SC alpha = SC(2);
  SC beta = SC(10);

  vector_t xBlock(blockMap);
  vector_t yBlock(blockMap);
  xBlock.putScalar(SC(1));
  initVector(yBlock);
  ABlock.apply(xBlock, yBlock, Teuchos::NO_TRANS, alpha, beta);

  vector_t xCyclic(domainMapCyclic);
  vector_t yCyclic(rangeMapCyclic);
  xCyclic.putScalar(SC(1));
  initVector(yCyclic);
  ACyclic.apply(xCyclic, yCyclic, Teuchos::NO_TRANS, alpha, beta);

  using magnitude_t = typename Teuchos::ScalarTraits<SC>::magnitudeType;
//  const magnitude_t tol = magnitude_t(10) 
//                        * Teuchos::ScalarTraits<magnitude_t>::eps();
  const magnitude_t tol = 0.000005;

  TEUCHOS_TEST_FLOATING_EQUALITY(yBlock.norm1(), yCyclic.norm1(), tol,
                                 out, success);
  TEUCHOS_TEST_FLOATING_EQUALITY(yBlock.norm2(), yCyclic.norm2(), tol,
                                 out, success);
  TEUCHOS_TEST_FLOATING_EQUALITY(yBlock.normInf(), yCyclic.normInf(), tol,
                                 out, success);
}

#define UNIT_TEST_GROUP( SC, LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, Bug9391, SC, LO, GO, NT)

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_SLGN( UNIT_TEST_GROUP )

} // anonymous
