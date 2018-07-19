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
#include "TpetraCore_ETIHelperMacros.h"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Distributor.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Teuchos_CommHelpers.hpp"

namespace { // anonymous

// Test for Issue #1752
//
// Description of Bug
// ==================
// In the "slow" path of the Kokkos version of the 4 argument
// Distributor::doPosts, there was a bug in which the variable numPacketsTo_p
// was not properly incremented from its initial value of 0.  Not incrementing
// numPacketsTo_p fooled Distributor in to thinking there were not any packets
// to send to participated processes.  When Distributor::doWaits was later
// called, it would hang, waiting on the sends that were never initiated because
// of numPacketsTo_p being equal to 0.
//
// See https://github.com/trilinos/Trilinos/issues/1752 for more details
//
// This test creates a matrix whose rows are uniquely owned, but not evenly
// (uniformly) distributed among processors.  This matrix is then imported, in
// both forward and reverse modes, in to matrices that are uniformly
// distributed.

////////////////////////////////////////////////////////////////////////////////
// Fill the matrix
template<class CrsMatrixType>
void fill_and_complete(CrsMatrixType& matrix)
{
  using Teuchos::tuple;
  typedef typename CrsMatrixType::scalar_type ST;
  typedef typename CrsMatrixType::local_ordinal_type LO;
  typedef typename CrsMatrixType::global_ordinal_type GO;
  typedef Tpetra::global_size_t GST;

  const ST neg_one = static_cast<ST>(-1.0);
  const ST two = static_cast<ST>(2.0);

  // Fill the sparse matrix, one row at a time.
  auto map = matrix.getRowMap();
  auto my_num_rows = map->getNodeNumElements();
  auto num_rows = map->getGlobalNumElements();
  for (LO i=0; i<static_cast<LO>(my_num_rows); i++) {
    auto gbl_row = map->getGlobalElement(i);
    // A(0, 0:1) = [2, -1]
    if (gbl_row == 0) {
      matrix.insertGlobalValues(
          gbl_row,
          tuple<GO>(gbl_row, gbl_row+1),
          tuple<ST> (two, neg_one));
    }
    // A(N-1, N-2:N-1) = [-1, 2]
    else if (static_cast<GST>(gbl_row) == num_rows - 1) {
      matrix.insertGlobalValues(
          gbl_row,
          tuple<GO>(gbl_row - 1, gbl_row),
          tuple<ST>(neg_one, two));
    }
    // A(i, i-1:i+1) = [-1, 2, -1]
    else {
      matrix.insertGlobalValues(
          gbl_row,
          tuple<GO>(gbl_row - 1, gbl_row, gbl_row + 1),
          tuple<ST>(neg_one, two, neg_one));
    }
  }
  matrix.fillComplete();
}

////////////////////////////////////////////////////////////////////////////////
// Check that matrix has the expected values
template<class CrsMatrixType>
std::pair<int, std::string> check_matrix(CrsMatrixType& matrix)
{
  using Teuchos::ArrayView;
  typedef typename CrsMatrixType::scalar_type ST;
  typedef typename CrsMatrixType::local_ordinal_type LO;
  typedef typename CrsMatrixType::global_ordinal_type GO;
  typedef typename ArrayView<LO>::size_type size_type;
  typedef Tpetra::global_size_t GST;

  int ierr = 0;
  std::ostringstream os;

  const ST neg_one = static_cast<ST>(-1.0);
  const ST two = static_cast<ST>(2.0);

  // Fill the sparse matrix, one row at a time.
  auto map = matrix.getRowMap();
  auto my_num_rows = map->getNodeNumElements();
  auto num_rows = matrix.getGlobalNumRows();

  for (LO i=0; i<static_cast<LO>(my_num_rows); i++) {
    auto gbl_row = map->getGlobalElement(i);
    ArrayView<const LO> cols;
    ArrayView<const ST> vals;
    matrix.getLocalRowView(i, cols, vals);

    std::map<GO,ST> expected;
    if (gbl_row == 0) {
      // A(0, 0:1) = [2, -1]
      expected[gbl_row] = two;
      expected[gbl_row+1] = neg_one;
    }
    else if (static_cast<GST>(gbl_row) == num_rows - 1) {
      // A(N-1, N-2:N-1) = [-1, 2]
      expected[gbl_row-1] = neg_one;
      expected[gbl_row] = two;
    }
    else {
      // A(i, i-1:i+1) = [-1, 2, -1]
      expected[gbl_row-1] = neg_one;
      expected[gbl_row] = two;
      expected[gbl_row+1] = neg_one;
    }

    if (static_cast<size_type>(expected.size()) != cols.size()) {
      ierr++;
      os << " Error: expected row " << gbl_row
         << " to have " << expected.size()
         << " non-zero columns, but got "
         << cols.size() << std::endl;
      continue;
    }

    for (typename ArrayView<const ST>::size_type j=0; j<cols.size(); j++) {
      auto gbl_col = matrix.getColMap()->getGlobalElement(cols[j]);
      if (vals[j] != expected[gbl_col]) {
        ierr++;
        os << " Error: expected entry [" << gbl_row << "," << gbl_col << "]"
           << " to be " << expected[gbl_col]
           << ", but got " << vals[j] << std::endl;
        continue;
      }
    }
  }

  return std::make_pair(ierr, os.str());
}

////////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST(Distributor, ReverseDistributeToNonuniformMap)
{
  using std::endl;
  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::reduceAll;
  using Teuchos::outArg;

  typedef Tpetra::CrsMatrix<> matrix_type;
  // typedef typename matrix_type::scalar_type ST; // unused
  typedef typename matrix_type::local_ordinal_type LO;
  typedef typename matrix_type::global_ordinal_type GO;
  typedef typename matrix_type::node_type NT;
  typedef Tpetra::global_size_t GST;
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef Tpetra::Export<LO, GO> export_type;
  typedef Tpetra::Import<LO, GO> import_type;

  int gblSuccess = 0; // output argument

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  auto world_size = comm->getSize();
  auto my_rank = comm->getRank();

  TEUCHOS_TEST_FOR_EXCEPTION(
      world_size <= 1,
      std::invalid_argument,
      "Test is only meaningful with number of processors > 1");

  size_t my_num_rows = 4;
  GST num_rows = static_cast<GST>(my_num_rows * world_size);

  // Create a unique map
  const GO idx_b = 0;
  RCP<const map_type> unique_map;
  {
    Teuchos::Array<GO> unique_idx;
    const GST invalid = Teuchos::OrdinalTraits<GST>::invalid();
    for (size_t i=0; i<my_num_rows; i++)
      unique_idx.push_back(static_cast<GO>(my_rank + world_size * i));
    unique_map = rcp(new map_type(invalid, unique_idx(), idx_b, comm));
  }

  // Create a default map
  RCP<const map_type> default_map = rcp(new map_type(num_rows, idx_b, comm));

  // Create matrix using the unique map
  out << std::right << std::setfill('*') <<  std::setw(80)
      << "Proc " << my_rank << ": Creating matrix with unique map" << endl;
  matrix_type unique_mtx(unique_map, 3, Tpetra::StaticProfile);
  fill_and_complete(unique_mtx);

  // Sanity check the unique matrix
  {
    auto check = check_matrix(unique_mtx);
    int lclSuccess = check.first != 0 ? 0 : 1;
    reduceAll<int, int> (*comm, Teuchos::REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      if (my_rank == 0)
        out << "Initial check_matrix resulted in the following errors: " << endl;
      Tpetra::Details::gathervPrint(out, check.second, *comm);
      if (my_rank == 0)
        out << endl << "Abandoning test; no point in continuing." << endl;
      return;
    }
  }

  // Import from the unique to the default map
  import_type importer(unique_map, default_map);

  // Matrix built with default map
  matrix_type default_mtx_fwd(default_map, 3, Tpetra::DynamicProfile);

  // Do a forward import (an import operation using an Import plan) from the
  // unique matrix to the default matrix.  i.e., communicate entries
  // in the unique matrix in to their corresponding locations in the
  // default matrix.
  out << std::right << std::setfill('*') <<  std::setw(80)
      << "Proc " << my_rank << ": Doing forward mode import" << endl;
  default_mtx_fwd.doImport(unique_mtx, importer, Tpetra::INSERT);
  default_mtx_fwd.fillComplete(unique_map, unique_map);

  // Sanity check the matrix
  {
    auto check = check_matrix(default_mtx_fwd);
    int lclSuccess = check.first != 0 ? 0 : 1;
    reduceAll<int, int> (*comm, Teuchos::REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      if (my_rank == 0)
        out << "Forward mode import resulted in the following errors: " << endl;
      Tpetra::Details::gathervPrint(out, check.second, *comm);
      if (my_rank == 0)
        out << endl << "Abandoning test; no point in continuing." << endl;
      return;
    }
  }

  // Matrix built with default map
  matrix_type default_mtx_rev(default_map, 3, Tpetra::DynamicProfile);

  // Do a reverse mode import (an import operation using an Export plan) from the
  // unique matrix to the default matrix.  i.e., communicate entries
  // in the unique matrix in to their corresponding locations in the
  // default matrix.
  out << std::right << std::setfill('*') <<  std::setw(80)
      << "Proc " << my_rank << ": Doing reverse mode import" << endl;
  export_type exporter(default_map, unique_map);
  default_mtx_rev.doImport(unique_mtx, exporter, Tpetra::INSERT);
  default_mtx_rev.fillComplete(unique_map, unique_map);

  // Sanity check the matrix
  {
    auto check = check_matrix(default_mtx_rev);
    int lclSuccess = check.first != 0 ? 0 : 1;
    reduceAll<int, int> (*comm, Teuchos::REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      if (my_rank == 0)
        out << "Reverse mode import resulted in the following errors: " << endl;
      Tpetra::Details::gathervPrint(out, check.second, *comm);
      if (my_rank == 0)
        out << endl << "Abandoning test; no point in continuing." << endl;
      return;
    }
  }

  TEST_EQUALITY_CONST(gblSuccess, 1);
}

} // namespace (anonymous)
