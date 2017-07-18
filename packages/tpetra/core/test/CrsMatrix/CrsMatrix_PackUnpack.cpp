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
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Distributor.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Tpetra_Details_packCrsMatrix.hpp"
#include "Tpetra_Details_unpackCrsMatrix.hpp"
#include "Teuchos_CommHelpers.hpp"
#include <random>
#include <set>

namespace { // anonymous

template<class T>
bool
essentially_equal(T a, T b) {
  typedef Kokkos::ArithTraits<T> KAT;

  const auto eps = KAT::eps ();
  return KAT::abs(a - b) <= ( (KAT::abs(a) > KAT::abs(b) ? KAT::abs(b) : KAT::abs(a)) * eps);
}

template<class MatrixType>
Teuchos::RCP<MatrixType>
generate_test_matrix (const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
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

  const LO num_row_per_proc = 4;
  Array<GO> row_gids(num_row_per_proc);

  GO start = static_cast<GO>(num_row_per_proc * world_rank);
  for (LO i=0; i<num_row_per_proc; ++i) {
    row_gids[i] = static_cast<GO>(start+i);
  }

  // Create random, unique column GIDs.
  const LO num_nz_cols = 100;
  std::random_device rand_dev;
  std::mt19937 generator(rand_dev());
  std::uniform_int_distribution<GO>  distr(1, 2000);
  std::set<GO> col_gids_set;
  typedef typename std::set<GO>::size_type SGO;
  SGO num_gids_in_set = static_cast<SGO>(num_nz_cols);
  col_gids_set.insert(0);
  int num_passes = 0;
  while (col_gids_set.size() < num_gids_in_set && num_passes <= 2000) {
    col_gids_set.insert(distr(generator));
    num_passes += 1;
  }
  if (col_gids_set.size() != num_gids_in_set) {
    throw std::runtime_error("Random column IDs not generate");
  }

  Array<GO> col_gids(col_gids_set.begin(), col_gids_set.end());

  const GO INVALID = OrdinalTraits<GO>::invalid();
  RCP<const MapType> row_map = rcp(new MapType(INVALID, row_gids(), 0, comm));
  RCP<const MapType> col_map = rcp(new MapType(INVALID, col_gids(), 0, comm));

  auto A = rcp (new MatrixType (row_map, col_map, num_nz_cols));

  Array<LO> columns(num_nz_cols);
  Array<SC> entries(num_nz_cols);
  for (LO j=0; j<num_nz_cols; ++j) columns[j] = j;
  for (LO i=0; i<num_row_per_proc; ++i) {
    //LO lcl_row = static_cast<LO>(start + i); // unused
    for (LO j=0; j<num_nz_cols; ++j) {
      entries[j] = static_cast<SC>(j+1);
    }
    A->insertLocalValues(i, columns(), entries());
  }
  A->fillComplete(col_map, row_map);

  return A;
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, PackUnpack, SC, LO, GO, NT)
{
  using std::endl;
  using Tpetra::CrsMatrix;
  using Tpetra::Details::gathervPrint;
  using Tpetra::Details::packCrsMatrix;
  using Tpetra::Details::unpackCrsMatrixAndCombine;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::Comm;
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::REDUCE_MIN;
  using Teuchos::REDUCE_SUM;
  using Teuchos::reduceAll;
  typedef CrsMatrix<SC, LO, GO, NT, false> MatrixType;
  typedef typename MatrixType::local_matrix_type LocalMatrixType;
  typedef typename NT::device_type device_type;
  typedef typename device_type::execution_space execution_space;

  int lclSuccess = 1; // to be revised below
  int gblSuccess = 0; // output argument

  RCP<const Comm<int> > comm =
    Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
  const int world_size = comm->getSize();
  const int world_rank = comm->getRank();

  // MPI_Comm_size returns the same value on all processes in the
  // communicator, so we don't need to do a test before exiting early.
  TEST_ASSERT( world_size % 2 == 0 );
  if (world_size % 2 != 0) {
    out << "This test must be run with multiples of 2 MPI processes, but was "
        << "run with " << world_size << " process"
        << (world_size != 1 ? "es" : "") << "." << endl;
    return;
  }

  out << "Creating matrix" << endl;

  auto A = generate_test_matrix<MatrixType> (comm);
  auto col_map = A->getColMap();
  auto row_map = A->getRowMap();

  out << "Preparing arguments for packCrsMatrix" << endl;

  LO num_loc_rows = static_cast<LO>(A->getNodeNumRows());
  Teuchos::Array<LO> exportLIDs (num_loc_rows); // input argument
  for (LO i=0; i < num_loc_rows; ++i) {
    exportLIDs[i] = static_cast<LO>(i); // pack all the rows
  }
  Array<char> exports; // output argument; to be realloc'd
  Array<size_t> numPacketsPerLID (num_loc_rows, 0); // output argument
  size_t constantNumPackets; // output argument
  Tpetra::Distributor distor (comm); // argument required, but not used
  LocalMatrixType A_lcl = A->getLocalMatrix (); // local matrix; input argument
  std::unique_ptr<std::string> errStr; // error string output argument

  out << "Calling packCrsMatrix" << endl;

  const bool lcl_pack_OK =
    packCrsMatrix (A_lcl, col_map->getLocalMap(), errStr,
                   exports, numPacketsPerLID (), constantNumPackets,
                   exportLIDs (), world_rank, distor);
  TEST_ASSERT( lcl_pack_OK );
  lclSuccess = success ? 1 : 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY_CONST( gblSuccess, 1 );
  if (gblSuccess != 1) {
    if (world_rank == 0) {
      out << "packCrsMatrix reported an error!" << endl;
    }
    const std::string errStr2 =
      errStr.get () == NULL ? std::string ("") : *errStr;
    gathervPrint (out, errStr2, *comm);
    out << endl << "Abandoning test; no point in continuing." << endl;
    return;
  }

  // Now make sure that the pack is correct by building a matrix of zeros and
  // replacing its values with the exported values from above (using the
  // Tpetra::REPLACE combine mode).  Then, the values in the create matrix
  // should be the same as the above matrix.
  out << "Building second matrix" << endl;
  auto graph = A->getCrsGraph();
  RCP<MatrixType> B = rcp(new MatrixType(graph));
  B->setAllToScalar(static_cast<SC>(0.));
  B->fillComplete();
  LocalMatrixType B_lcl = B->getLocalMatrix();

#ifdef KOKKOS_HAVE_SERIAL
  typedef typename device_type::execution_space ES;
  const bool atomic_updates = ! std::is_same<ES, Kokkos::Serial>::value;
#else
  const bool atomic_updates = true;
#endif // KOKKOS_HAVE_SERIAL

  out << "Calling unpackCrsMatrixAndCombine with "
      << "CombineMode=Tpetra::REPLACE" << endl;
  bool lcl_unpack_OK =
    unpackCrsMatrixAndCombine (B_lcl, col_map->getLocalMap (), errStr,
                               exportLIDs (), exports, numPacketsPerLID (),
                               constantNumPackets, world_rank, distor, Tpetra::REPLACE,
                               atomic_updates);
  TEST_ASSERT( lcl_unpack_OK );
  lclSuccess = success ? 1 : 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY_CONST( gblSuccess, 1 );
  if (gblSuccess != 1) {
    if (world_rank == 0) {
      out << "unpackCrsMatrixAndCombine reported an error!" << endl;
    }
    const std::string errStr2 =
      errStr.get () == NULL ? std::string ("") : *errStr;
    gathervPrint (out, errStr2, *comm);
    return; // no point in continuing
  }

  // The test below uses the host Tpetra::CrsMatrix interface to
  // compare matrix values.  Thus, we need to do a fence before
  // comparing matrix values, in order to ensure that changes made on
  // device are visible on host.
  execution_space::fence ();

  out << "Comparing matrices after unpackCrsMatrixAndCombine "
    "with CombineMode=REPLACE" << endl;
  {
    std::ostringstream errStrm;
    int lclNumErrors = 0;
    for (LO lclRow=0; lclRow<num_loc_rows; ++lclRow) {
      ArrayView<const LO> A_indices;
      ArrayView<const SC> A_values;
      A->getLocalRowView(lclRow, A_indices, A_values);

      ArrayView<const LO> B_indices;
      ArrayView<const SC> B_values;
      B->getLocalRowView(lclRow, B_indices, B_values);

      TEST_EQUALITY( A_indices.size (), B_indices.size () );

      int curNumErrors = 0;
      LO num_indices = static_cast<LO>(A_indices.size());
      for (LO i=0; i<num_indices; i++) {
        if (! essentially_equal<SC>(A_values[i], B_values[i])) {
          errStrm << "ERROR: Proc " << world_rank << ", row " << lclRow
                  << ", A[" << i << "]=" << A_values[i] << ", but "
                  <<   "B[" << i << "]=" << B_values[i] << "!\n";
          ++curNumErrors;
        }
      }
      lclNumErrors += curNumErrors;
    }
    TEST_ASSERT( lclNumErrors == 0 );

    int gblNumErrors = 0;
    reduceAll<int, int> (*comm, REDUCE_SUM, lclNumErrors, outArg (gblNumErrors));
    TEST_EQUALITY_CONST( gblNumErrors, 0 );
    if (gblNumErrors != 0) {
      if (world_rank == 0) {
        out << "unpackCrsMatrixAndCombine comparison found " << gblNumErrors
            << " error" << (gblNumErrors != 1 ? "s" : "") << "!" << endl;
      }
      gathervPrint (out, errStrm.str (), *comm);
      return; // no point in continuing
    }

    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      if (world_rank == 0) {
        out << "unpackCrsMatrixAndCombine comparison claims zero errors, "
          "but success is false on at least one process!" << endl;
      }
      gathervPrint (out, errStrm.str (), *comm);
      return; // no point in continuing
    }
  }

  // Now let's do the same thing, but use the Tpetra::ADD combine mode.
  // The resultant matrix should then be twice the original.
  out << "Calling unpackCrsMatrixAndCombine with "
      << "CombineMode=Tpetra::ADD" << endl;
  lcl_unpack_OK =
    unpackCrsMatrixAndCombine (B_lcl, col_map->getLocalMap (), errStr,
                               exportLIDs (), exports, numPacketsPerLID (),
                               constantNumPackets, world_rank, distor, Tpetra::ADD,
                               atomic_updates);
  TEST_ASSERT( lcl_unpack_OK );
  lclSuccess = success ? 1 : 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY_CONST( gblSuccess, 1 );
  if (gblSuccess != 1) {
    if (world_rank == 0) {
      out << "unpackCrsMatrixAndCombine reported an error!" << endl;
    }
    const std::string errStr2 =
      errStr.get () == NULL ? std::string ("") : *errStr;
    gathervPrint (out, errStr2, *comm);
    out << endl << "Abandoning test; no point in continuing." << endl;
    return;
  }

  // Loop through rows and compare matrix values
  out << "Comparing matrices after unpackCrsMatrixAndCombine "
      "with CombineMode=ADD" << endl;
  {
    std::ostringstream errStrm;
    int lclNumErrors = 0;
    for (LO loc_row=0; loc_row<num_loc_rows; ++loc_row) {
      ArrayView<const LO> A_indices;
      ArrayView<const SC> A_values;
      A->getLocalRowView(loc_row, A_indices, A_values);

      ArrayView<const LO> B_indices;
      ArrayView<const SC> B_values;
      B->getLocalRowView(loc_row, B_indices, B_values);

      TEST_EQUALITY( A_indices.size (), B_indices.size () );

      int curNumErrors = 0;
      LO num_indices = static_cast<LO>(A_indices.size());
      for (LO i=0; i<num_indices; i++) {
        if (! essentially_equal<SC>(2.0 * A_values[i], B_values[i])) {
          errStrm << "ERROR: Proc " << world_rank << ", row " << loc_row
                  << ", 2*A[" << i << "]=" << 2.0 * A_values[i] << ", but "
                  <<     "B[" << i << "]=" <<   B_values[i] << "!\n";
          ++curNumErrors;
        }
      }
      lclNumErrors += curNumErrors;
    }
    TEST_ASSERT( lclNumErrors == 0 );

    int gblNumErrors = 0;
    reduceAll<int, int> (*comm, REDUCE_SUM, lclNumErrors, outArg (gblNumErrors));
    TEST_EQUALITY_CONST( gblNumErrors, 0 );
    if (gblNumErrors != 0) {
      if (world_rank == 0) {
        out << "unpackCrsMatrixAndCombine comparison found " << gblNumErrors
            << " error" << (gblNumErrors != 1 ? "s" : "") << "!" << endl;
      }
      gathervPrint (out, errStrm.str (), *comm);
      return; // no point in continuing
    }

    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      if (world_rank == 0) {
        out << "unpackCrsMatrixAndCombine comparison claims zero errors, "
          "but success is false on at least one process!" << endl;
      }
      gathervPrint (out, errStrm.str (), *comm);
      return; // no point in continuing
    }
  }

  lclSuccess = success ? 1 : 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY_CONST( gblSuccess, 1 );

  if (gblSuccess != 1) {
    if (world_rank == 0) {
      out << "Error in packCrsMatrix!" << endl;
    }
    using ::Tpetra::Details::gathervPrint;
    const std::string errStr2 =
      errStr.get () == NULL ? std::string ("") : *errStr;
    gathervPrint (out, errStr2, *comm);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, PackError, SC, LO, GO, NT)
{
  using Tpetra::Details::gathervPrint;
  using Tpetra::Details::packCrsMatrix;
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

  out << "Creating matrix" << endl;

  auto A = generate_test_matrix<MatrixType> (comm);
  auto col_map = A->getColMap();
  auto row_map = A->getRowMap();

  out << "Calling packCrsMatrix" << endl;

  // Prepare arguments for pack.  This test is similar to the PackUnpack test,
  // but incorrect exportLIDs are sent in to induce a packing error.
  int lclSuccess = success ? 1 : 0;
  int gblSuccess = 0; // output argument
  std::ostringstream errStrm; // for error string local to each process

  bool lcl_pack_OK = false;
  {
    LO num_loc_rows = static_cast<LO>(A->getNodeNumRows());
    Array<LO> exportLIDs(num_loc_rows);
    // exportLIDs[i] should equal i, but we set it to i+2
    for (LO i=0; i<num_loc_rows; i++) {
      exportLIDs[i] = i + 2;
    }

    Array<char> exports;
    Array<size_t> numPacketsPerLID(num_loc_rows, 0);
    size_t constantNumPackets;
    Distributor distor(comm);
    LocalMatrixType A_lcl = A->getLocalMatrix();
    std::unique_ptr<std::string> errStr;
    lcl_pack_OK =
      packCrsMatrix (A_lcl, col_map->getLocalMap(), errStr,
                     exports, numPacketsPerLID (), constantNumPackets,
                     exportLIDs (), world_rank, distor);
    if (lcl_pack_OK) {
      // Local pack should not be OK!  We requested bad local IDs be exported!
      errStrm << "Proc " << world_rank
              << ": packCrsMatrix returned OK, but bad local IDs were requested!"
              << endl;
      lclSuccess = 0;
    }
  }

  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY( gblSuccess, 1 );
  if (gblSuccess != 1) {
    out << "packCrsMatrix failed to notice bad export IDs on some process!" << endl;
    gathervPrint (out, errStrm.str (), *comm);
  }

  {
    // Let's try this again, but send in the wrong number of exportLIDs
    LO num_loc_rows = static_cast<LO>(A->getNodeNumRows());
    // Note the -1!
    out << "Allocating ids... ";
    Array<LO> exportLIDs(num_loc_rows-1);
    for (LO i=0; i < num_loc_rows-1; ++i) {
      exportLIDs[i] = i;
    }
    out << "done" << endl;

    Array<char> exports;
    Array<size_t> numPacketsPerLID(num_loc_rows, 0);
    size_t constantNumPackets;
    Distributor distor(comm);
    LocalMatrixType A_lcl = A->getLocalMatrix();
    std::unique_ptr<std::string> errStr;
    out << "Calling packCrsMatrix" << endl;
    lcl_pack_OK =
      packCrsMatrix (A_lcl, col_map->getLocalMap(), errStr,
                     exports, numPacketsPerLID (), constantNumPackets,
                     exportLIDs (), world_rank, distor);
    if (lcl_pack_OK) {
      // Local pack should not be OK!  We requested too few local IDs be exported!
      errStrm << "Proc " << world_rank
              << ": packCrsMatrix returned OK, but too few IDs "
              << "were requested to be exported!"
              << endl;
      lclSuccess = 0;
    }
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
  }
}

#define UNIT_TEST_GROUP_SC_LO_GO_NO( SC, LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, PackUnpack, SC, LO, GO, NT) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, PackError, SC, LO, GO, NT)

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(UNIT_TEST_GROUP_SC_LO_GO_NO)

} // namespace (anonymous)
