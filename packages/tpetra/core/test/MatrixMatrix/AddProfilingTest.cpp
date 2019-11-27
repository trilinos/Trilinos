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
#include <Tpetra_TestingUtilities.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include "TpetraExt_MatrixMatrix.hpp"
#include "Tpetra_MatrixIO.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrixMultiplyOp.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"
#include "impl/Kokkos_Timer.hpp"
#include "TpetraExt_MatrixMatrix.hpp"

#include <cmath>
#include <unordered_set>

namespace Tpetra
{
namespace AddProfiling
{

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Comm;

static const int numRows = 1000;

//Produce a random matrix with given nnz per global row
template<typename SC, typename LO, typename GO, typename NT>
RCP<Tpetra::CrsMatrix<SC, LO, GO, NT>> getTestMatrix(RCP<Tpetra::Map<LO, GO, NT>>& domainRangeMap,
    RCP<Tpetra::Map<LO, GO, NT>>& colMap, int seed, RCP<const Comm<int>>& comm)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  const LO maxNnzPerRow = 100;
  //use the range map as row map (non-overlapping distribution)
  auto mat = rcp(new Tpetra::CrsMatrix<SC, LO, GO, NT>(domainRangeMap, colMap, maxNnzPerRow));
  //get consistent results between trials
  srand(comm->getRank() * 7 + 42 + seed);
  ArrayView<const GO> ownedRows = domainRangeMap->getNodeElementList();
  ArrayView<const GO> ownedCols = colMap->getNodeElementList();
  for(GO row : ownedRows)
  {
    Teuchos::Array<SC> vals(maxNnzPerRow);
    Teuchos::Array<GO> inds(maxNnzPerRow);
    for(int j = 0; j < maxNnzPerRow; j++)
    {
      //get random between 0.0 and 10.0
      vals[j] = 10.0 * rand() / RAND_MAX;
      inds[j] = ownedCols[rand() % ownedCols.size()];
    }
    //No guarantee that the indices are unique, so use sumInto
    mat->sumIntoGlobalValues(row, inds(), vals());
  }
  mat->fillComplete(domainRangeMap, domainRangeMap);
  return mat;
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_AddProfiling, sorted, SC, LO, GO, NT)
{
  typedef Tpetra::CrsMatrix<SC, LO, GO, NT> crs_matrix_type;
  typedef Tpetra::Map<LO, GO, NT> map_type;
  RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
  if(comm->getRank() == 0)
    std::cout << "Running sorted add test on " << comm->getSize() << " MPI ranks.\n";
  RCP<map_type> rowMap = rcp(new map_type(numRows, 0, comm));
  RCP<map_type> colMap = rcp(new map_type(numRows, 0, comm));
  RCP<crs_matrix_type> A = getTestMatrix<SC, LO, GO, NT>(rowMap, colMap, 1, comm);
  RCP<crs_matrix_type> B = getTestMatrix<SC, LO, GO, NT>(rowMap, colMap, 2, comm);
  Kokkos::Impl::Timer addTimer;
  auto one = Teuchos::ScalarTraits<SC>::one();
  RCP<crs_matrix_type> C = MatrixMatrix::add(one, false, *A, one, false, *B);
  double tkernel = addTimer.seconds();
  std::cout << "sorted (kernel): addition took " << tkernel << "s.\n";
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_AddProfiling, different_col_maps, SC, LO, GO, NT)
{
  typedef Tpetra::CrsMatrix<SC, LO, GO, NT> crs_matrix_type;
  typedef Tpetra::Map<LO, GO, NT> map_type;
  RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
  if(comm->getRank() == 0)
    std::cout << "Running \"different col maps\" add test on " << comm->getSize() << " MPI ranks.\n";
  RCP<map_type> rowMap = rcp(new map_type(numRows, 0, comm));
  //Number of random global columns to include in column map
  GO colMapSize = numRows / comm->getSize();
  auto getRandColMap = [&] (int seed) -> Teuchos::Array<GO>
  {
    srand(seed);
    std::unordered_set<GO> cols;
    for(GO i = 0; i < colMapSize; i++)
    {
      GO col;
      do
      {
        col = rand() % numRows;
      }
      while(cols.find(col) != cols.end());
      cols.insert(col);
    }
    Teuchos::Array<GO> colList(colMapSize);
    {
      int i = 0;
      for(auto c : cols)
        colList[i++] = c;
    }
    return colList;
  };
  Teuchos::Array<GO> colMapARaw = getRandColMap(comm->getRank() + 12342);
  Teuchos::Array<GO> colMapBRaw = getRandColMap(comm->getRank() + 87345);
  RCP<map_type> colMapA = rcp(new map_type(numRows, colMapARaw(), 0, comm));
  RCP<map_type> colMapB = rcp(new map_type(numRows, colMapBRaw(), 0, comm));
  RCP<crs_matrix_type> A = getTestMatrix<SC, LO, GO, NT>(rowMap, colMapA, 1, comm);
  RCP<crs_matrix_type> B = getTestMatrix<SC, LO, GO, NT>(rowMap, colMapB, 2, comm);
  Kokkos::Impl::Timer addTimer;
  auto one = Teuchos::ScalarTraits<SC>::one();
  RCP<crs_matrix_type> C = MatrixMatrix::add(one, false, *A, one, false, *B);
  double tkernel = addTimer.seconds();
  std::cout << "mismatched col maps (kernel): addition took " << tkernel << "s.\n";
}

#define UNIT_TEST_GROUP_SC_LO_GO_NO( SC, LO, GO, NT )			\
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_AddProfiling, sorted, SC, LO, GO, NT) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_AddProfiling, different_col_maps, SC, LO, GO, NT)

  TPETRA_ETI_MANGLING_TYPEDEFS()
  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( UNIT_TEST_GROUP_SC_LO_GO_NO )

} //namespace AddProfiling
} //namespace Tpetra

