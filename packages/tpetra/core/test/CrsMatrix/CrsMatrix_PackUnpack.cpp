// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "TpetraCore_ETIHelperMacros.h"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Tpetra_Details_packCrsMatrix.hpp"
#include "Tpetra_Details_unpackCrsMatrixAndCombine.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Kokkos_ArithTraits.hpp"
#include <random>
#include <set>

namespace { // anonymous

#define NUM_ROW_PER_PROC 4
#define NUM_NZ_COLS 100

using Tpetra::TestingUtilities::getDefaultComm;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::Comm;
using Teuchos::outArg;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::REDUCE_MIN;
using Teuchos::REDUCE_SUM;
using Teuchos::reduceAll;
using Tpetra::Details::gathervPrint;
using Tpetra::Details::packCrsMatrix;
using Tpetra::Details::unpackCrsMatrixAndCombine;
using std::endl;

template<class T>
bool
essentially_equal(T a, T b) {
  typedef Kokkos::ArithTraits<T> KAT;
  const auto eps = KAT::eps();
  return KAT::abs(a - b) <= ( (KAT::abs(a) > KAT::abs(b) ? KAT::abs(b) : KAT::abs(a)) * eps);
}

template<class CrsMatrixType>
Teuchos::RCP<CrsMatrixType>
generate_test_matrix (const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
{
  using crs_matrix_type = CrsMatrixType;
  using SC = typename crs_matrix_type::scalar_type;
  using LO = typename crs_matrix_type::local_ordinal_type;
  using GO = typename crs_matrix_type::global_ordinal_type;
  using NT = typename crs_matrix_type::node_type;
  using MapType = Tpetra::Map<LO, GO, NT>;

  const int world_rank = comm->getRank();

  const LO num_row_per_proc = NUM_ROW_PER_PROC; // 4;
  Array<GO> row_gids(num_row_per_proc);

  GO start = static_cast<GO>(num_row_per_proc * world_rank);
  for (LO i=0; i<num_row_per_proc; ++i) {
    row_gids[i] = static_cast<GO>(start+i);
  }

  // Create random, unique column GIDs.
  const LO num_nz_cols = NUM_NZ_COLS; //100;
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
  TEUCHOS_TEST_FOR_EXCEPTION(col_gids_set.size() != num_gids_in_set,
      std::runtime_error, "Random column IDs not generated");

  Array<GO> col_gids(col_gids_set.begin(), col_gids_set.end());

  const GO INVALID = Teuchos::OrdinalTraits<GO>::invalid();
  RCP<const MapType> row_map = rcp(new MapType(INVALID, row_gids(), 0, comm));
  RCP<const MapType> col_map = rcp(new MapType(INVALID, col_gids(), 0, comm));

  auto A = rcp (new crs_matrix_type (row_map, col_map, num_nz_cols));

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

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, PackThenUnpackAndCombine, SC, LO, GO, NT)
{
  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using device_type = typename NT::device_type;
  using execution_space = typename device_type::execution_space;

  int lclSuccess = 1; // to be revised below
  int gblSuccess = 0; // output argument

  RCP<const Comm<int> > comm = getDefaultComm();
  const int world_rank = comm->getRank();

  out << "Creating matrix" << endl;

  auto A = generate_test_matrix<crs_matrix_type> (comm);
  auto col_map = A->getColMap();
  auto row_map = A->getRowMap();

  out << "Preparing arguments for packCrsMatrix" << endl;

  LO num_loc_rows = static_cast<LO>(A->getLocalNumRows());
  Array<LO> exportLIDs (num_loc_rows); // input argument
  for (LO i=0; i < num_loc_rows; ++i) {
    exportLIDs[i] = static_cast<LO>(i); // pack all the rows
  }
  Array<char> exports; // output argument; to be realloc'd
  Array<size_t> numPacketsPerLID (num_loc_rows, 0); // output argument
  size_t constantNumPackets; // output argument

  out << "Calling packCrsMatrix" << endl;

  {
    int local_op_ok;
    std::ostringstream msg;
    try {
      packCrsMatrix (*A, exports, numPacketsPerLID (),
                     exportLIDs ().getConst (),
                     constantNumPackets);
      local_op_ok = 1;
    }
    catch (std::exception& e) {
      local_op_ok = 0;
      msg << e.what();
    }
    TEST_ASSERT(local_op_ok == 1);
    lclSuccess = success ? 1 : 0;
    reduceAll (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      if (world_rank == 0) {
        out << "packCrsMatrix reported an error!" << endl;
      }
      gathervPrint (out, msg.str(), *comm);
      out << endl << "Abandoning test; no point in continuing." << endl;
      return;
    }
  }

  // Now make sure that the pack is correct by building a matrix of zeros and
  // replacing its values with the exported values from above (using the
  // Tpetra::REPLACE combine mode).  Then, the values in the create matrix
  // should be the same as the above matrix.
  out << "Building second matrix" << endl;
  auto graph = A->getCrsGraph();
  RCP<crs_matrix_type> B (new crs_matrix_type (graph));
  B->setAllToScalar(SC{-1.});
  B->fillComplete();

  out << "Calling unpackCrsMatrixAndCombine with "
      << "CombineMode=Tpetra::REPLACE" << endl;

  {
    int local_op_ok;
    std::ostringstream msg;
    try {
      unpackCrsMatrixAndCombine (*B, exports, numPacketsPerLID (),
                                 exportLIDs ().getConst (),
                                 constantNumPackets,
                                 Tpetra::REPLACE);
      local_op_ok = 1;
    }
    catch (std::exception& e) {
      local_op_ok = 0;
      msg << e.what();
    }

    TEST_ASSERT(local_op_ok == 1);
    lclSuccess = success ? 1 : 0;
    reduceAll (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      if (world_rank == 0) {
        out << "unpackCrsMatrixAndCombine reported an error!" << endl;
      }
      gathervPrint(out, msg.str(), *comm);
      return; // no point in continuing
    }
  }

  // The test below uses the host Tpetra::CrsMatrix interface to
  // compare matrix values.  Thus, we need to do a fence before
  // comparing matrix values, in order to ensure that changes made on
  // device are visible on host.
  execution_space ().fence ();

  out << "Comparing matrices after unpackCrsMatrixAndCombine "
    "with CombineMode=REPLACE" << endl;
  {
    std::ostringstream errStrm;
    int lclNumErrors = 0;
    for (LO lclRow=0; lclRow<num_loc_rows; ++lclRow) {
      typename crs_matrix_type::local_inds_host_view_type A_indices;
      typename crs_matrix_type::values_host_view_type A_values;
      A->getLocalRowView(lclRow, A_indices, A_values);

      typename crs_matrix_type::local_inds_host_view_type B_indices;
      typename crs_matrix_type::values_host_view_type B_values;
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
        else
        {
          errStrm << "INFO: Proc " << world_rank << ", row " << lclRow
                  << ", A[" << i << "]=" << A_values[i] << ", and "
                  <<   "B[" << i << "]=" << B_values[i] << "!\n";
        }
      }
      lclNumErrors += curNumErrors;
    }
    TEST_ASSERT( lclNumErrors == 0 );

    int gblNumErrors = 0;
    reduceAll (*comm, REDUCE_SUM, lclNumErrors, outArg (gblNumErrors));
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
    reduceAll (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
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
  {
    int local_op_ok;
    std::ostringstream msg;
    try {
      unpackCrsMatrixAndCombine (*B, exports, numPacketsPerLID (),
                                 exportLIDs ().getConst (),
                                 constantNumPackets,
                                 Tpetra::ADD);
      local_op_ok = 1;
    }
    catch (std::exception& e) {
      local_op_ok = 0;
      msg << e.what();
    }
    TEST_ASSERT(local_op_ok == 1);
    lclSuccess = success ? 1 : 0;
    reduceAll (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      if (world_rank == 0) {
        out << "unpackCrsMatrixAndCombine reported an error!" << endl;
      }
      gathervPrint (out, msg.str(), *comm);
      out << endl << "Abandoning test; no point in continuing." << endl;
      return;
    }
  }

  // Loop through rows and compare matrix values
  out << "Comparing matrices after unpackCrsMatrixAndCombine "
      "with CombineMode=ADD" << endl;
  {
    std::ostringstream errStrm;
    int lclNumErrors = 0;
    for (LO loc_row=0; loc_row<num_loc_rows; ++loc_row) {
      typename crs_matrix_type::local_inds_host_view_type A_indices;
      typename crs_matrix_type::values_host_view_type A_values;
      A->getLocalRowView(loc_row, A_indices, A_values);

      typename crs_matrix_type::local_inds_host_view_type B_indices;
      typename crs_matrix_type::values_host_view_type B_values;
      B->getLocalRowView(loc_row, B_indices, B_values);
//      std::cout << "A_values: " << A_values << "\n";
//      std::cout << "B_values: " << B_values << "\n";
//      std::cout << std::flush;

      TEST_EQUALITY( A_indices.size (), B_indices.size () );

      int curNumErrors = 0;
      LO num_indices = static_cast<LO>(A_indices.size());
      for (LO i=0; i<num_indices; i++) {
        if (! essentially_equal<SC>(SC(2.0) * A_values[i], B_values[i])) {
          errStrm << "ERROR: Proc " << world_rank << ", row " << loc_row
                  << ", 2*A[" << i << "]=" << SC(2.0) * A_values[i] << ", but "
                  <<     "B[" << i << "]=" <<   B_values[i] << "!\n";
          ++curNumErrors;
        }
      }
      lclNumErrors += curNumErrors;
    }
    TEST_ASSERT( lclNumErrors == 0 );

    int gblNumErrors = 0;
    reduceAll (*comm, REDUCE_SUM, lclNumErrors, outArg (gblNumErrors));
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
    reduceAll (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
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

}

// PackWithError sends intentionally bad inputs to pack/unpack to make sure
// that CrsMatrix will detect the bad inputs and return the correct
// error diagnostics.
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, PackWithError, SC, LO, GO, NT)
{
  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;

  RCP<const Comm<int> > comm = getDefaultComm();
  const int world_rank = comm->getRank();

  out << "Creating matrix" << endl;

  auto A = generate_test_matrix<crs_matrix_type> (comm);
  auto col_map = A->getColMap();
  auto row_map = A->getRowMap();

  out << "Calling packCrsMatrix" << endl;

  // Prepare arguments for pack.  This test is similar to the
  // PackThenUnpackAndCombine test,
  // but incorrect exportLIDs are sent in to induce a packing error.
  int lclSuccess = success ? 1 : 0;
  int gblSuccess = 0; // output argument
  std::ostringstream errStrm; // for error string local to each process

  {
    LO num_loc_rows = static_cast<LO>(A->getLocalNumRows());
    Array<LO> exportLIDs(num_loc_rows);
    // exportLIDs[i] should equal i, but we set it to i+2
    for (LO i=0; i<num_loc_rows; i++) {
      exportLIDs[i] = i + 2;
    }

    Array<char> exports;
    Array<size_t> numPacketsPerLID(num_loc_rows, 0);
    size_t constantNumPackets;
    {
      int local_op_ok;
      std::ostringstream msg;
      try {
        packCrsMatrix (*A, exports, numPacketsPerLID (),
                       exportLIDs ().getConst (),
                       constantNumPackets);
        local_op_ok = 1;
      }
      catch (std::exception& e) {
        local_op_ok = 0;
        msg << e.what();
      }
      if (local_op_ok == 1) {
        // Local pack should not be OK!  We requested bad local IDs be exported!
        errStrm << "Proc " << world_rank
                << ": packCrsMatrix returned OK, but bad local IDs were requested!"
                << endl;
        lclSuccess = 0;
      }
    }
  }

  reduceAll (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY( gblSuccess, 1 );
  if (gblSuccess != 1) {
    out << "packCrsMatrix failed to notice bad export IDs on some process!" << endl;
    gathervPrint (out, errStrm.str (), *comm);
  }

  {
    // Let's try this again, but send in the wrong number of exportLIDs
    LO num_loc_rows = static_cast<LO>(A->getLocalNumRows());
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
    out << "Calling packCrsMatrix" << endl;
    {
      int local_op_ok;
      std::ostringstream msg;
      try {
        packCrsMatrix (*A, exports, numPacketsPerLID (),
                       exportLIDs ().getConst (),
                       constantNumPackets);
        local_op_ok = 1;
      }
      catch (std::exception& e) {
        local_op_ok = 0;
        msg << e.what();
      }
      if (local_op_ok == 1) {
        // Local pack should not be OK!  We requested too few local IDs be exported!
        errStrm << "Proc " << world_rank
                << ": packCrsMatrix returned OK, but too few IDs "
                << "were requested to be exported!"
                << endl;
        lclSuccess = 0;
      }
    }
    reduceAll (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
  }
}

// PackPartial exercises the use case where only *some* of the rows of the
// matrix are packed.  This test exists because an errant check in pack/unpack
// checked that the number of rows to pack was equal to the number of rows in
// the matrix.  However, this is not a requirement, and that particular check
// caused existing code to fail.  See Issues #1374 and #1408
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, PackPartial, SC, LO, GO, NT)
{
  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using device_type = typename NT::device_type;
  using execution_space = typename device_type::execution_space;

  int lclSuccess = 1; // to be revised below
  int gblSuccess = 0; // output argument

  RCP<const Comm<int> > comm = getDefaultComm();
  const int world_rank = comm->getRank();

  out << "Creating matrix" << endl;

  auto A = generate_test_matrix<crs_matrix_type> (comm);
  auto col_map = A->getColMap();
  auto row_map = A->getRowMap();

  out << "Preparing arguments for packCrsMatrix" << endl;

  LO num_loc_rows = static_cast<LO>(A->getLocalNumRows());
  Array<LO> exportLIDs (num_loc_rows); // input argument
  for (LO i=0; i < num_loc_rows; ++i) {
    if (i%2 != 0) continue;  // pack only even rows!
    exportLIDs[i] = static_cast<LO>(i);
  }
  Array<char> exports; // output argument; to be realloc'd
  Array<size_t> numPacketsPerLID (num_loc_rows, 0); // output argument
  size_t constantNumPackets; // output argument

  out << "Calling packCrsMatrix" << endl;

  {
    int local_op_ok;
    std::ostringstream msg;
    try {
      packCrsMatrix (*A, exports, numPacketsPerLID (),
                     exportLIDs ().getConst (),
                     constantNumPackets);
      local_op_ok = 1;
    }
    catch (std::exception& e) {
      local_op_ok = 0;
      msg << e.what();
    }
    TEST_ASSERT(local_op_ok == 1);
    lclSuccess = success ? 1 : 0;
    reduceAll (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      if (world_rank == 0) {
        out << "packCrsMatrix reported an error!" << endl;
      }
      gathervPrint (out, msg.str(), *comm);
      out << endl << "Abandoning test; no point in continuing." << endl;
      return;
    }
  }

  // Now make sure that the pack is correct by building a matrix of zeros and
  // replacing its values with the exported values from above (using the
  // Tpetra::REPLACE combine mode).  Then, the values in the create matrix
  // should be the same as the above matrix in the rows that were packed.
  out << "Building second matrix" << endl;
  auto graph = A->getCrsGraph ();
  RCP<crs_matrix_type> B (new crs_matrix_type (graph));
  B->setAllToScalar (SC {});
  B->fillComplete ();

  out << "Calling unpackCrsMatrixAndCombine with "
      << "CombineMode=Tpetra::REPLACE" << endl;

  {
    int local_op_ok;
    std::ostringstream msg;
    try {
      unpackCrsMatrixAndCombine (*B, exports, numPacketsPerLID (),
                                 exportLIDs ().getConst (),
                                 constantNumPackets,
                                 Tpetra::REPLACE);
      local_op_ok = 1;
    }
    catch (std::exception& e) {
      local_op_ok = 0;
      msg << e.what();
    }

    TEST_ASSERT(local_op_ok == 1);
    lclSuccess = success ? 1 : 0;
    reduceAll (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      if (world_rank == 0) {
        out << "unpackCrsMatrixAndCombine reported an error!" << endl;
      }
      gathervPrint(out, msg.str(), *comm);
      return; // no point in continuing
    }
  }

  // The test below uses the host Tpetra::CrsMatrix interface to
  // compare matrix values.  Thus, we need to do a fence before
  // comparing matrix values, in order to ensure that changes made on
  // device are visible on host.
  execution_space ().fence ();

  out << "Comparing matrices after unpackCrsMatrixAndCombine "
    "with CombineMode=REPLACE" << endl;
  {
    std::ostringstream errStrm;
    int lclNumErrors = 0;
    for (LO lclRow=0; lclRow<num_loc_rows; ++lclRow) {
      typename crs_matrix_type::local_inds_host_view_type A_indices;
      typename crs_matrix_type::values_host_view_type A_values;
      A->getLocalRowView(lclRow, A_indices, A_values);

      typename crs_matrix_type::local_inds_host_view_type B_indices;
      typename crs_matrix_type::values_host_view_type B_values;
      B->getLocalRowView(lclRow, B_indices, B_values);

      TEST_EQUALITY( A_indices.size (), B_indices.size () );

      int curNumErrors = 0;
      LO num_indices = static_cast<LO>(A_indices.size());
      for (LO i=0; i<num_indices; i++) {
        SC a_val;
        if (lclRow%2 == 0) {
          // Only even rows were packed!
          a_val = A_values[i];
        } else {
          // Odd rows should still be zero!
          a_val = 0.;
        }
        if (! essentially_equal<SC>(a_val, B_values[i])) {
          errStrm << "ERROR: Proc " << world_rank << ", row " << lclRow
                  << ", A[" << i << "]=" << a_val << ", but "
                  <<   "B[" << i << "]=" << B_values[i] << "!\n";
          ++curNumErrors;
        }
      }
      lclNumErrors += curNumErrors;
    }
    TEST_ASSERT( lclNumErrors == 0 );

    int gblNumErrors = 0;
    reduceAll (*comm, REDUCE_SUM, lclNumErrors, outArg (gblNumErrors));
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
    reduceAll (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
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
}

#define UNIT_TEST_GROUP( SC, LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, PackThenUnpackAndCombine, SC, LO, GO, NT) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, PackWithError, SC, LO, GO, NT) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, PackPartial, SC, LO, GO, NT)

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(UNIT_TEST_GROUP)

} // namespace (anonymous)
