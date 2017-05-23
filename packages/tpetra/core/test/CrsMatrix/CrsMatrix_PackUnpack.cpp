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
#include <set>
#include <random>
#include <Teuchos_UnitTestHarness.hpp>
#include <Tpetra_ConfigDefs.hpp>
#include <TpetraCore_ETIHelperMacros.h>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_TestingUtilities.hpp>
#include <Tpetra_Details_packCrsMatrix.hpp>
#include <Tpetra_Details_unpackCrsMatrix.hpp>
#include <Tpetra_Distributor.hpp>
#include "Tpetra_Details_gathervPrint.hpp"

namespace { // anonymous

template<class T>
bool
essentially_equal(T a, T b) {
  typedef Kokkos::ArithTraits<T> KAT;

  const auto eps = KAT::eps ();
  return KAT::abs(a - b) <= ( (KAT::abs(a) > KAT::abs(b) ? KAT::abs(b) : KAT::abs(a)) * eps);
}

template<class MatrixType>
void
generate_test_matrix(Teuchos::RCP<MatrixType>& A,
                     Teuchos::RCP<const Teuchos::Comm<int> >& comm)
{
  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::Array;
  using Teuchos::OrdinalTraits;
  typedef typename MatrixType::scalar_type SC;
  typedef typename MatrixType::local_ordinal_type LO;
  typedef typename MatrixType::global_ordinal_type GO;
  typedef typename MatrixType::node_type NT;
  typedef Tpetra::Map<LO, GO, NT> MapType;

  const int world_rank = comm->getRank();

  const int num_row_per_proc = 4;
  Array<GO> row_gids(num_row_per_proc);

  int start = num_row_per_proc * world_rank;
  for (int i=0; i<num_row_per_proc; i++)
    row_gids[i] = static_cast<GO>(start+i);

  // Create random, unique column GIDs.
  const int num_nz_cols = 100;
  std::random_device rand_dev;
  std::mt19937 generator(rand_dev());
  std::uniform_int_distribution<GO>  distr(1, 2000);
  std::set<GO> col_gids_set;
  col_gids_set.insert(0);
  while (col_gids_set.size() < num_nz_cols) {
    col_gids_set.insert(distr(generator));
  }
  Array<GO> col_gids(col_gids_set.begin(), col_gids_set.end());

  comm->barrier();

  const GO INVALID = OrdinalTraits<GO>::invalid();
  RCP<const MapType> row_map = rcp(new MapType(INVALID, row_gids(), 0, comm));
  RCP<const MapType> col_map = rcp(new MapType(INVALID, col_gids(), 0, comm));

  comm->barrier();

  A = rcp(new MatrixType(row_map, col_map, num_nz_cols));

  comm->barrier();

  Array<LO> columns(num_nz_cols);
  Array<SC> entries(num_nz_cols);
  for (LO j=0; j<static_cast<LO>(num_nz_cols); ++j) columns[j] = j;
  for (int i=0; i<num_row_per_proc; ++i) {
    //LO lcl_row = static_cast<LO>(start + i); // unused
    for (int j=0; j<entries.size(); ++j)
      entries[j] = static_cast<SC>(j+1);
    A->insertLocalValues(i, columns(), entries());
  }

  comm->barrier();

  A->fillComplete(col_map, row_map);
  comm->barrier();

}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, PackUnpack, SC, LO, GO, NT)
{
  using std::endl;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Comm;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Tpetra::CrsMatrix;
  using Tpetra::DefaultPlatform;

  using Tpetra::Details::packCrsMatrix;
  using Tpetra::Details::unpackCrsMatrixAndCombine;
  using Tpetra::Distributor;
  typedef CrsMatrix<SC, LO, GO, NT, false> MatrixType;
  typedef typename MatrixType::local_matrix_type LocalMatrixType;

  RCP<const Comm<int> > comm = DefaultPlatform::getDefaultPlatform().getComm();

  const int world_size = comm->getSize();
  const int world_rank = comm->getRank();

  if (world_size%2 != 0) {
    out << "This test must be run with multiples of 2 MPI processes, but was "
        << "run with " << world_size << " process"
        << (world_size != 1 ? "es" : "") << "." << endl;
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::runtime_error, "This test must be run with multiples of "
      "2 MPI processes, but was run with " << world_size << " process"
      << (world_size != 1 ? "es" : "") << ".");
  }
  out << "Proc " << world_rank << ": Running test" << endl;
  comm->barrier();

  out << "Proc " << world_rank << ": Creating matrix" << endl;
  RCP<MatrixType> A;
  generate_test_matrix(A, comm);
  auto col_map = A->getColMap();
  auto row_map = A->getRowMap();

  out << "Proc " << world_rank << ": Calling packCrsMatrix" << endl;
  comm->barrier();

  // Prepare arguments for pack
  size_t num_loc_rows = A->getNodeNumRows();
  std::vector<LO> ids(num_loc_rows);
  for (size_t i=0; i<num_loc_rows; i++) ids[i] = static_cast<LO>(i);
  ArrayView<const LO> exportLIDs(ids);
  Array<char> exports;
  std::vector<size_t> numpacks(num_loc_rows, 0);
  ArrayView<size_t> numPacketsPerLID(numpacks);

  size_t constantNumPackets;
  Distributor distor(comm);

  LocalMatrixType A_lcl = A->getLocalMatrix();

  std::unique_ptr<std::string> errStr;
  const bool lcl_pack_OK =
    packCrsMatrix (A_lcl, col_map->getLocalMap(), errStr,
        exports, numPacketsPerLID, constantNumPackets, exportLIDs,
         world_rank, distor);
  TEST_ASSERT( lcl_pack_OK );

  // Now make sure that the pack is correct by building a matrix of zeros and
  // replacing its values with the exported values from above (using the
  // Tpetra::REPLACE combine mode).  Then, the values in the create matrix
  // should be the same as the above matrix.
  out << "Proc " << world_rank << ": Building second matrix" << endl;
  comm->barrier();
  auto graph = A->getCrsGraph();
  RCP<MatrixType> B = rcp(new MatrixType(graph));
  B->setAllToScalar(static_cast<SC>(0.));
  B->fillComplete();
  LocalMatrixType B_lcl = B->getLocalMatrix();

#ifdef KOKKOS_HAVE_SERIAL
  typedef typename NT::device_type device_type;
  typedef typename device_type::execution_space ES;
  const bool atomic_updates = ! std::is_same<ES, Kokkos::Serial>::value;
#else
  const bool atomic_updates = true;
#endif // KOKKOS_HAVE_SERIAL

  out << "Proc " << world_rank << ": Calling unpackCrsMatrixAndCombine with "
      << "CombineMode=Tpetra::REPLACE" << endl;
  comm->barrier();
  bool lcl_unpack_OK;
  lcl_unpack_OK = unpackCrsMatrixAndCombine(
      B_lcl, col_map->getLocalMap(), errStr,
      exportLIDs, exports, numPacketsPerLID,
      constantNumPackets, world_rank, distor, Tpetra::REPLACE, atomic_updates);
  TEST_ASSERT( lcl_unpack_OK );

  // Loop through rows and compare matrix values
  out << "Proc " << world_rank << ": comparing matrices after first "
      << "unpackCrsMatrixAndCombine" << endl;
  for (LO loc_row=0; loc_row<static_cast<LO>(num_loc_rows); loc_row++) {
    ArrayView<const LO> A_indices;
    ArrayView<const SC> A_values;
    A->getLocalRowView(loc_row, A_indices, A_values);

    ArrayView<const LO> B_indices;
    ArrayView<const SC> B_values;
    B->getLocalRowView(loc_row, B_indices, B_values);

    TEST_ASSERT(A_indices.size() == B_indices.size())

    {
      int errors = 0;
      std::ostringstream os;
      for (int i=0; i<A_indices.size(); i++) {
        if (!essentially_equal<SC>(A_values[i], B_values[i])) {
          os << "ERROR: Proc " << world_rank << ", row " << loc_row
              << ", A[" << i << "]=" << A_values[i] << ", but "
              <<   "B[" << i << "]=" << B_values[i] << "!\n";
          errors += 1;
        }
      }
      if (errors) {
        Tpetra::Details::gathervPrint(out, os.str(), *(col_map->getComm()));
        TEST_ASSERT(false);
      }
    }
  }

  // Now let's do the same thing, but use the Tpetra::ADD combine mode.
  // The resultant matrix should then be twice the original.
  out << "Proc " << world_rank << ": Calling unpackCrsMatrixAndCombinex with "
      << "CombineMode=Tpetra::ADD" << endl;
  comm->barrier();
  lcl_unpack_OK = unpackCrsMatrixAndCombine(
      B_lcl, col_map->getLocalMap(), errStr,
      exportLIDs, exports, numPacketsPerLID,
      constantNumPackets, world_rank, distor, Tpetra::ADD, atomic_updates);
  TEST_ASSERT( lcl_unpack_OK );

  // Loop through rows and compare matrix values
  out << "Proc " << world_rank << ": comparing matrices after first "
      << "unpackCrsMatrixAndCombine" << endl;
  for (LO loc_row=0; loc_row<static_cast<LO>(num_loc_rows); loc_row++) {
    ArrayView<const LO> A_indices;
    ArrayView<const SC> A_values;
    A->getLocalRowView(loc_row, A_indices, A_values);

    ArrayView<const LO> B_indices;
    ArrayView<const SC> B_values;
    B->getLocalRowView(loc_row, B_indices, B_values);

    TEST_ASSERT(A_indices.size() == B_indices.size())

    {
      int errors = 0;
      std::ostringstream os;
      for (int i=0; i<A_indices.size(); i++) {
        if (!essentially_equal<SC>(2.0 * A_values[i], B_values[i])) {
          os << "ERROR: Proc " << world_rank << ", row " << loc_row
              << ", 2*A[" << i << "]=" << 2.0 * A_values[i] << ", but "
              <<     "B[" << i << "]=" <<   B_values[i] << "!\n";
          errors += 1;
        }
      }
      if (errors) {
        Tpetra::Details::gathervPrint(out, os.str(), *(col_map->getComm()));
        TEST_ASSERT(false);
      }
    }
  }

  int lclSuccess = success ? 1 : 0;
  int gblSuccess = 0; // output argument

  using Teuchos::outArg;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  reduceAll<int, int> (* (col_map->getComm ()), REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY_CONST( gblSuccess, 1 );

  if (gblSuccess != 1) {
    if (world_rank == 0) {
      out << "Error in packCrsMatrix!" << endl;
    }
    using ::Tpetra::Details::gathervPrint;
    const std::string errStr2 =
      errStr.get () == NULL ? std::string ("") : *errStr;
    gathervPrint (out, errStr2, * (col_map->getComm ()));
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, PackError, SC, LO, GO, NT)
{
  using Tpetra::Details::gathervPrint;
  using Tpetra::Details::packCrsMatrix;
  using Tpetra::Details::unpackCrsMatrixAndCombine;
  using Tpetra::CrsMatrix;
  using Tpetra::DefaultPlatform;
  using Tpetra::Distributor;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::Comm;
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using std::endl;

  typedef CrsMatrix<SC, LO, GO, NT, false> MatrixType;
  typedef typename MatrixType::local_matrix_type LocalMatrixType;

  RCP<const Comm<int> > comm = DefaultPlatform::getDefaultPlatform().getComm();

  const int world_size = comm->getSize();
  const int world_rank = comm->getRank();

  if (world_size%2 != 0) {
    out << "This test must be run with multiples of 2 MPI processes, but was "
        << "run with " << world_size << " process"
        << (world_size != 1 ? "es" : "") << "." << endl;
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::runtime_error, "This test must be run with multiples of "
      "2 MPI processes, but was run with " << world_size << " process"
      << (world_size != 1 ? "es" : "") << ".");
  }
  out << "Proc " << world_rank << ": Running test" << endl;
  comm->barrier();

  out << "Proc " << world_rank << ": Creating matrix" << endl;
  RCP<MatrixType> A;
  generate_test_matrix(A, comm);
  auto col_map = A->getColMap();
  auto row_map = A->getRowMap();

  out << "Proc " << world_rank << ": Calling packCrsMatrix" << endl;
  comm->barrier();

  // Prepare arguments for pack.  This test is similar to the PackUnpack test,
  // but incorrect exportLIDs are sent in to induce a packing error.
  int lclSuccess = success ? 1 : 0;
  int gblSuccess = 0; // output argument
  std::ostringstream errStrm; // for error string local to each process

  bool lcl_pack_OK = false;
  {
    size_t num_loc_rows = A->getNodeNumRows();
    std::vector<LO> ids(num_loc_rows);
    // ids[i] should equal i, but we set it to i+2
    for (size_t i=0; i<num_loc_rows; i++) ids[i] = static_cast<LO>(i+2);
    ArrayView<const LO> exportLIDs(ids);
    Array<char> exports;
    std::vector<size_t> numpacks(num_loc_rows, 0);
    ArrayView<size_t> numPacketsPerLID(numpacks);
    size_t constantNumPackets;
    Distributor distor(comm);
    LocalMatrixType A_lcl = A->getLocalMatrix();
    std::unique_ptr<std::string> errStr;
    lcl_pack_OK =
      packCrsMatrix (A_lcl, col_map->getLocalMap(), errStr,
                     exports, numPacketsPerLID, constantNumPackets, exportLIDs,
                     world_rank, distor);
    if (lcl_pack_OK) {
      // Local pack should not be OK!  We requested bad local IDs be exported!
      errStrm << "Proc " << world_rank
              << ": packCrsMatrix returned OK, but bad local IDs were requested!"
              << endl;
      lclSuccess = 0;
    }
  }

  reduceAll<int, int> (* (col_map->getComm ()), REDUCE_MIN,
                       lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY( gblSuccess, 1 );
  if (gblSuccess != 1) {
    out << "packCrsMatrix failed to notice bad export IDs on some process!" << endl;
    gathervPrint (out, errStrm.str (), * (col_map->getComm ()));
  }

  {
    // Let's try this again, but send in the wrong number of exportLIDs
    size_t num_loc_rows = A->getNodeNumRows();
    // Note the -1!
    out << "Proc " << world_rank << ": Allocating ids... ";
    std::vector<LO> ids(num_loc_rows-1);
    out << "done" << endl;
    for (size_t i=0; i<num_loc_rows-1; i++) ids[i] = static_cast<LO>(i);
    out << "Proc " << world_rank << ": Creating ids ArrayView... ";
    ArrayView<const LO> exportLIDs(ids);
    out << "done" << endl;
    Array<char> exports;
    std::vector<size_t> numpacks(num_loc_rows, 0);
    ArrayView<size_t> numPacketsPerLID(numpacks);
    size_t constantNumPackets;
    Distributor distor(comm);
    LocalMatrixType A_lcl = A->getLocalMatrix();
    std::unique_ptr<std::string> errStr;
    out << "Proc " << world_rank << ": Calling packCrsMatrix" << endl;
    lcl_pack_OK =
      packCrsMatrix (A_lcl, col_map->getLocalMap(), errStr,
                     exports, numPacketsPerLID, constantNumPackets, exportLIDs,
                     world_rank, distor);
    if (lcl_pack_OK) {
      // Local pack should not be OK!  We requested too few local IDs be exported!
      errStrm << "Proc " << world_rank
              << ": packCrsMatrix returned OK, but too few IDs "
              << "were requested to be exported!"
              << endl;
      lclSuccess = 0;
    }
    reduceAll<int, int> (* (col_map->getComm ()), REDUCE_MIN,
                         lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
  }
}

#define UNIT_TEST_GROUP_SC_LO_GO_NO( SC, LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, PackUnpack, SC, LO, GO, NT) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, PackError, SC, LO, GO, NT)

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(UNIT_TEST_GROUP_SC_LO_GO_NO)

} // namespace (anonymous)
