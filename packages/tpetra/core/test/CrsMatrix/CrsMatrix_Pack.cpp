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

#include <Teuchos_UnitTestHarness.hpp>
#include <Tpetra_ConfigDefs.hpp>
#include <TpetraCore_ETIHelperMacros.h>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_TestingUtilities.hpp>
#include <Tpetra_Details_packCrsMatrix.hpp>
#include <Tpetra_Distributor.hpp>

namespace { // anonymous

using std::endl;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Comm;
using Teuchos::Array;
using Teuchos::ArrayRCP;
using Teuchos::ArrayView;
using Teuchos::OrdinalTraits;
using Tpetra::Map;
using Tpetra::CrsMatrix;
using Tpetra::DefaultPlatform;

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, Pack1, SC, LO, GO, NT)
{
  using Tpetra::Details::packCrsMatrix;
  using Tpetra::Distributor;
  typedef Map<LO, GO, NT> MapType;
  typedef CrsMatrix<SC, LO, GO, NT, false> MatrixType;
  typedef typename MatrixType::local_matrix_type LocalMatrixType;
  RCP<const Comm<int> > comm = DefaultPlatform::getDefaultPlatform().getComm();

  const int world_size = comm->getSize();
  const int world_rank = comm->getRank();

  if (world_size != 2) {
    out << "This test must be run with exactly 2 MPI processes, but you ran "
        << "it with " << world_size << " process"
        << (world_size != 1 ? "es" : "") << "." << endl;
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::runtime_error, "This test must be run with exactly "
      "2 MPI processes, but you ran it with " << world_size << " process"
      << (world_size != 1 ? "es" : "") << ".");
  }
  out << "Proc " << world_rank << ": Running test" << endl;
  comm->barrier();

  // Create a 4x4 matrix
  // Rows 0,1 on P0
  // Rows 2,3 on P1
  Array<GO> row_gids, col_gids;
  GO start = static_cast<GO>(2*world_rank);
  row_gids.push_back(start);
  row_gids.push_back(start+1);

  // All columns on each
  for (GO i=0; i<4; i++) col_gids.push_back(i);

  out << "Proc " << world_rank << ": Creating row and column Maps" << endl;
  comm->barrier ();

  const GO INVALID = OrdinalTraits<GO>::invalid();
  RCP<const MapType> row_map = rcp(new MapType(INVALID, row_gids(), 0, comm));
  RCP<const MapType> col_map = rcp(new MapType(INVALID, col_gids(), 0, comm));

  out << "Proc " << world_rank << ": Creating matrix" << endl;
  comm->barrier();

  size_t count = 4;
  MatrixType A(row_map, col_map, count);

  out << "Proc " << world_rank << ": Filling matrix" << endl;
  comm->barrier();

  Array<LO> columns(4);
  Array<SC> entries(4, 7);
  for (int j=0; j<4; ++j) columns[j] = j;
  for (int i=0; i<2; ++i) A.insertLocalValues(i, columns(), entries());

  out << "Proc " << world_rank << ": Calling fillComplete" << endl;
  comm->barrier();

  A.fillComplete(col_map, row_map);
  comm->barrier();

  out << "Proc " << world_rank << ": Calling packCrsMatrix" << endl;
  comm->barrier();

  // Prepare arguments for pack
  LO ids[] = {0, 1};
  ArrayView<const LO> exportLIDs(ids, 2);

  Array<char> exports;

  size_t numpacks[] = {0, 0};
  ArrayView<size_t> numPacketsPerLID(numpacks, 2);

  size_t constantNumPackets;
  Distributor distor(comm);

  LocalMatrixType A_lcl = A.getLocalMatrix();

  packCrsMatrix<SC,LO,GO,NT,false>(exportLIDs, exports, numPacketsPerLID,
                  constantNumPackets, distor, *col_map, A_lcl);

  out << "Proc " << world_rank << ": Done with test" << endl;
}

#define UNIT_TEST_GROUP_SC_LO_GO_NO( SC, LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, Pack1, SC, LO, GO, NT)

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(UNIT_TEST_GROUP_SC_LO_GO_NO)

} // namespace (anonymous)
