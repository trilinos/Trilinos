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

#include <Tpetra_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <TpetraCore_ETIHelperMacros.h>
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_TestingUtilities.hpp>

#include <iostream>
#include <vector>

namespace { // anonymous

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
create_matrix()
{
  // Create matrix on ran 0
  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::Array;
  using GST = Tpetra::global_size_t;
  using map_type = Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
  using matrix_type = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  auto comm = Tpetra::getDefaultComm();
  const int rank = comm->getRank();
  const GST invalid = Teuchos::OrdinalTraits<GST>::invalid();


  // Map with everything on rank zero
  const GlobalOrdinal num_rows = 25;
  const LocalOrdinal num_ent_per_row = 5;
  RCP<const map_type> map = rcp(
    new map_type(invalid, (rank == 0 ? num_rows : 0), 0, comm)
  );

  // Matrix using mapOne -- all entries on rank 0
  RCP<matrix_type> matrix = rcp(new matrix_type(map, num_ent_per_row));

  if (rank == 0)
  {
    Array<GlobalOrdinal> gids(num_ent_per_row);
    Array<Scalar> values(num_ent_per_row);
    for (int i = 0; i < num_ent_per_row; i++) values[i] = 1.;
    for (GlobalOrdinal i = 0; i < num_rows; i++) {
      // five-point stencil
      int nnz = 0;

      gids[nnz++] = i;

      if (((i-1) >= 0) && (((i-1)%num_ent_per_row) < (i%num_ent_per_row)))
        gids[nnz++] = i-1; // west nbor

      if (((i+1) < num_rows) && (((i+1) % num_ent_per_row) > (i%num_ent_per_row)))
        gids[nnz++] = i+1; // east nbor

      if ((i-num_ent_per_row) >= 0)
        gids[nnz++] = i-num_ent_per_row; // south nbor

      if ((i+num_ent_per_row) < num_rows)
       gids[nnz++] = i+num_ent_per_row; // north nbor

      matrix->insertGlobalValues(i, gids(0, nnz), values(0, nnz));
    }
  }
  matrix->fillComplete();
  return matrix;
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, Bug10008, SC, LO, GO, NT)
{
  // matrix whose entries live on rank 0
  auto A = create_matrix<SC, LO, GO, NT>();

  // Map across all processors
  auto proc_zero_map = A->getRowMap();
  auto comm = proc_zero_map->getComm();
  using map_type = typename std::decay<decltype(*proc_zero_map)>::type;
  auto num_rows = proc_zero_map->getGlobalNumElements();
  Teuchos::RCP<const map_type> dist_map = Teuchos::rcp(new map_type(num_rows, 0, comm));

  // Exporter from one to many
  using export_type = Tpetra::Export<LO, GO, NT>;
  export_type exporter(proc_zero_map, dist_map);

  // Import to spread matrix across all processors
  using matrix_type = typename std::decay<decltype(*A)>::type;
  Teuchos::RCP<matrix_type> B = rcp(new matrix_type(dist_map, 0));
  B->doExport(*A, exporter, Tpetra::INSERT);
  B->fillComplete();
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, Bug10008_1, SC, LO, GO, NT)
{
  // Demonstrate that importAndFillComplete fails for one-to-many 
  // transfers (entire matrix initially on one rank, redistributed to all ranks)

  // matrix whose entries live on rank 0
  auto A = create_matrix<SC, LO, GO, NT>();

  // Map across all processors
  auto proc_zero_map = A->getRowMap();
  auto comm = proc_zero_map->getComm();
  using map_type = typename std::decay<decltype(*proc_zero_map)>::type;
  auto num_rows = proc_zero_map->getGlobalNumElements();
  Teuchos::RCP<const map_type> dist_map = Teuchos::rcp(new map_type(num_rows, 0, comm));

  // Importer from one to many 
  using import_type = Tpetra::Import<LO, GO, NT>;
  import_type importer(proc_zero_map, dist_map);

  // Import to spread matrix across all processors
  using matrix_type = typename std::decay<decltype(*A)>::type;
  Teuchos::RCP<matrix_type> B = rcp(new matrix_type(dist_map, 0));
  B->doImport(*A, importer, Tpetra::INSERT);
  B->fillComplete();
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, Bug10008_2, SC, LO, GO, NT)
{
  // matrix whose entries live on rank 0
  auto A = create_matrix<SC, LO, GO, NT>();

  // Map across all processors
  auto proc_zero_map = A->getRowMap();
  auto comm = proc_zero_map->getComm();
  using map_type = typename std::decay<decltype(*proc_zero_map)>::type;
  auto num_rows = proc_zero_map->getGlobalNumElements();
  Teuchos::RCP<const map_type> dist_map = Teuchos::rcp(new map_type(num_rows, 0, comm));

  // Exporter from one to many
  using export_type = Tpetra::Export<LO, GO, NT>;
  export_type exporter(proc_zero_map, dist_map);

  // Import to spread matrix across all processors
  using matrix_type = typename std::decay<decltype(*A)>::type;
  Teuchos::RCP<matrix_type> B = rcp(new matrix_type(dist_map, 0));
  // #10008 reports that this call to importAndFillComplete throws and hangs
  A->exportAndFillComplete(B, exporter, dist_map, dist_map);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, Bug10008_3, SC, LO, GO, NT)
{
  // matrix whose entries live on rank 0
  auto A = create_matrix<SC, LO, GO, NT>();

  // Map across all processors
  auto proc_zero_map = A->getRowMap();
  auto comm = proc_zero_map->getComm();
  using map_type = typename std::decay<decltype(*proc_zero_map)>::type;
  auto num_rows = proc_zero_map->getGlobalNumElements();
  Teuchos::RCP<const map_type> dist_map = Teuchos::rcp(new map_type(num_rows, 0, comm));

  // Importer from one to many
  using import_type = Tpetra::Import<LO, GO, NT>;
  import_type importer(proc_zero_map, dist_map);

  // Import to spread matrix across all processors
  using matrix_type = typename std::decay<decltype(*A)>::type;
  Teuchos::RCP<matrix_type> B = rcp(new matrix_type(dist_map, 0));
  // #10008 reports that this call to importAndFillComplete throws and hangs
  A->importAndFillComplete(B, importer, dist_map, dist_map);
}

#define UNIT_TEST_GROUP_SC_LO_GO_NO( SC, LO, GO, NT )                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, Bug10008, SC, LO, GO, NT) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, Bug10008_1, SC, LO, GO, NT) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, Bug10008_2, SC, LO, GO, NT) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, Bug10008_3, SC, LO, GO, NT)

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(UNIT_TEST_GROUP_SC_LO_GO_NO)

} // anonymous
