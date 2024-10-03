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
#include "Tpetra_Distributor.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Tpetra_Details_packCrsMatrix.hpp"
#include "Tpetra_Details_unpackCrsMatrixAndCombine.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Kokkos_ArithTraits.hpp"
#include <random>
#include <set>

namespace { // anonymous


using Tpetra::TestingUtilities::getDefaultComm;
using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::Comm;
using Teuchos::outArg;
using Tpetra::Details::gathervPrint;
using Tpetra::Details::packCrsMatrix;
using Tpetra::Details::unpackCrsMatrixAndCombine;

// Create and return a simple example CrsMatrix, with row distribution
// over the given Map.
//
// CrsMatrixType: The type of the Tpetra::CrsMatrix specialization to use.
template<class SC, class LO, class GO, class NT>
RCP<const Tpetra::CrsMatrix<SC, LO,GO,NT>>
generate_crs_matrix(const RCP<const Tpetra::Map<LO,GO,NT>>& map)
{
  using Teuchos::tuple;
  using matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;

  const SC two     = static_cast<SC>( 2.0);
  const SC neg_one = static_cast<SC>(-1.0);

  // Create a Tpetra::Matrix using the Map, with dynamic allocation.
  auto A = rcp(new matrix_type(map, 3));

  // const size_t numMyElements = map->getLocalNumElements ();
  // The list of global elements owned by this MPI process.
  const Tpetra::global_size_t num_gbl_inds = map->getGlobalNumElements();
  auto my_gbl_elems = map->getLocalElementList();
  for (auto it : my_gbl_elems) {
    const auto gbl_row = map->getGlobalElement(it);
    if (gbl_row == 0)
      A->insertGlobalValues(gbl_row, tuple(gbl_row, gbl_row+1), tuple(two, neg_one));
    else if (static_cast<Tpetra::global_size_t>(gbl_row) == num_gbl_inds - 1)
      A->insertGlobalValues(gbl_row, tuple(gbl_row-1, gbl_row), tuple(neg_one, two));
    else
      A->insertGlobalValues(gbl_row, tuple(gbl_row - 1, gbl_row, gbl_row + 1), tuple(neg_one, two, neg_one));
  }
  // Finish up the matrix.
  A->fillComplete();
  return A;
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, ImportToStaticMatrix, SC, LO, GO, NT)
{
  // Set up Tpetra typedefs.
  using map_type = Tpetra::Map<LO, GO, NT>;
  using matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using device_type = typename NT::device_type;
  using execution_space = typename device_type::execution_space;
  const char prefix[] = "ImportToStaticMatrix: ";

  Tpetra::Details::Behavior::disable_verbose_behavior();

  auto comm = getDefaultComm();
  const auto my_rank = comm->getRank();

  auto loc_success = 1; // to be revised below
  auto gbl_success = 0; // output argument

  // The global number of rows in the matrix A to create.  We scale
  // this relative to the number of (MPI) processes, so that no matter
  // how many MPI processes you run, every process will have 10 rows.
  const Tpetra::global_size_t num_gbl_inds = 10 * comm->getSize();
  const GO idx_base = 0;

  // Construct a Map that is global, but puts all the equations on MPI Proc 0.
  RCP<const map_type> proc_zero_map;
  {
    const size_t num_loc_inds = (my_rank == 0) ? num_gbl_inds : 0;
    proc_zero_map = rcp(new map_type(num_gbl_inds, num_loc_inds, idx_base, comm));
  }

  // Construct a Map that puts approximately the same number of
  // equations on each processor.
  auto global_map = rcp(new map_type(num_gbl_inds, idx_base, comm, Tpetra::GloballyDistributed));

  // Create a sparse matrix using procZeroMap.
  auto A = generate_crs_matrix<SC,LO,GO,NT>(proc_zero_map);
  comm->barrier();

  // We've created a sparse matrix that lives entirely on Process 0.
  // Now we want to distribute it over all the processes.
  //
  // Redistribute the matrix.  Since both the source and target Maps
  // are one-to-one, we could use either an Import or an Export.  If
  // only the source Map were one-to-one, we would have to use an
  // Import; if only the target Map were one-to-one, we would have to
  // use an Export.  We do not allow redistribution using Import or
  // Export if neither source nor target Map is one-to-one.
  RCP<matrix_type> B;
  {
    using export_type = Tpetra::Export<LO,GO,NT>;
    export_type exporter(proc_zero_map, global_map);
    comm->barrier();

    Tpetra::Details::Behavior::enable_verbose_behavior();
    // Make a new sparse matrix whose row map is the global Map.
    out << prefix << "Creating empty matrix from global map.\n";
    B = rcp(new matrix_type(global_map, 0));
    out << prefix << "Empty matrix from global map created.\n";

    // Redistribute the data, NOT in place, from matrix A (which lives
    // entirely on Proc 0) to matrix B (which is distributed evenly over
    // the processes).
    out << prefix << "Performing export operation\n";
    B->doExport(*A, exporter, Tpetra::INSERT);
    out << prefix << "Export operation done.\n";
  }
  B->fillComplete();

  auto loc_num_errs = 0;

  // The test below uses the host Tpetra::CrsMatrix interface to
  // compare matrix values.  Thus, we need to do a fence before
  // comparing matrix values, in order to ensure that changes made on
  // device are visible on host.
  execution_space().fence();

  { // test output
    std::ostringstream errStrm;
    TEST_ASSERT(loc_num_errs == 0);

    auto gbl_num_errs = 0;
    Teuchos::reduceAll<int, int>(*comm,Teuchos::REDUCE_SUM,loc_num_errs,outArg(gbl_num_errs));
    TEST_EQUALITY_CONST(gbl_num_errs, 0);
    if (gbl_num_errs != 0) {
      if (my_rank == 0) {
        out << "export in to static matrix failed with " << gbl_num_errs
            << " error" << (gbl_num_errs != 1 ? "s" : "") << "!\n";
      }
      gathervPrint(out, errStrm.str (), *comm);
      return; // no point in continuing
    }

    loc_success = success ? 1 : 0;
    Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_MIN, loc_success, outArg(gbl_success));
    TEST_EQUALITY_CONST(gbl_success, 1 );
    if (gbl_success != 1) {
      if (my_rank == 0) {
        out << "export in to static matrix comparison claims zero errors, "
          "but success is false on at least one process!\n";
      }
      gathervPrint (out, errStrm.str (), *comm);
      return; // no point in continuing
    }
  }
}


TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, ImportToStaticMatrixLocal, SC, LO, GO, NT)
{
  // Set up Tpetra typedefs.
  using map_type = Tpetra::Map<LO, GO, NT>;
  using matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using device_type = typename NT::device_type;
  using execution_space = typename device_type::execution_space;
  const char prefix[] = "ImportToStaticMatrixLocal: ";

  Tpetra::Details::Behavior::disable_verbose_behavior();

  auto comm = getDefaultComm();
  const auto my_rank = comm->getRank();

  auto loc_success = 1; // to be revised below
  auto gbl_success = 0; // output argument

  // The global number of rows in the matrix A to create.  We scale
  // this relative to the number of (MPI) processes, so that no matter
  // how many MPI processes you run, every process will have 10 rows.
  const Tpetra::global_size_t num_gbl_inds = 10 * comm->getSize();
  const GO idx_base = 0;

  // Construct a Map that is global, but puts all the equations on MPI Proc 0.
  RCP<const map_type> proc_zero_map;
  {
    const size_t num_loc_inds = (my_rank == 0) ? num_gbl_inds : 0;
    proc_zero_map = rcp(new map_type(num_gbl_inds, num_loc_inds, idx_base, comm));
  }

  // Create a sparse matrix using procZeroMap.
  auto A = generate_crs_matrix<SC,LO,GO,NT>(proc_zero_map);
  comm->barrier();

  // We've created a sparse matrix that lives entirely on Process 0.
  // Now we want to distribute it over all the processes.
  //
  // Redistribute the matrix.  Since both the source and target Maps
  // are one-to-one, we could use either an Import or an Export.  If
  // only the source Map were one-to-one, we would have to use an
  // Import; if only the target Map were one-to-one, we would have to
  // use an Export.  We do not allow redistribution using Import or
  // Export if neither source nor target Map is one-to-one.
  {
    using export_type = Tpetra::Export<LO,GO,NT>;
    export_type exporter(proc_zero_map, proc_zero_map);
    comm->barrier();

    Tpetra::Details::Behavior::enable_verbose_behavior();
    // Make a new sparse matrix whose row map is the global Map.
    out << prefix << "Creating empty matrix from proc 0 map.\n";
    auto B = matrix_type(proc_zero_map, 0);
    out << prefix << "Empty matrix from proc 0 map created.\n";

    // Redistribute the data from matrix A (which lives
    // entirely on Proc 0) to matrix B which also lives on processor 0
    out << prefix << "Performing export operation\n";
    B.doExport(*A, exporter, Tpetra::INSERT);
    out << prefix << "Export operation done.\n";
    B.fillComplete();
  }

  {
    auto global_map = rcp(new map_type(num_gbl_inds, idx_base, comm, Tpetra::GloballyDistributed));
    using export_type = Tpetra::Export<LO,GO,NT>;
    export_type exporter1(proc_zero_map, global_map);
    comm->barrier();
    out << prefix << "Creating empty matrix from global map.\n";
    auto B = matrix_type(global_map, 0);
    B.doExport(*A, exporter1, Tpetra::INSERT);
    B.fillComplete();

    using export_type = Tpetra::Export<LO,GO,NT>;
    export_type exporter2(global_map, global_map);
    comm->barrier();
    out << prefix << "Creating another empty matrix from global map.\n";
    auto C = matrix_type(global_map, 0);
    out << prefix << "Exporting between global map matrices\n";
    C.doExport(B, exporter2, Tpetra::INSERT);
    C.fillComplete();
  }

  auto loc_num_errs = 0;

  // The test below uses the host Tpetra::CrsMatrix interface to
  // compare matrix values.  Thus, we need to do a fence before
  // comparing matrix values, in order to ensure that changes made on
  // device are visible on host.
  execution_space().fence();

  { // test output
    std::ostringstream errStrm;
    TEST_ASSERT(loc_num_errs == 0);

    auto gbl_num_errs = 0;
    Teuchos::reduceAll<int, int>(*comm,Teuchos::REDUCE_SUM,loc_num_errs,outArg(gbl_num_errs));
    TEST_EQUALITY_CONST(gbl_num_errs, 0);
    if (gbl_num_errs != 0) {
      if (my_rank == 0) {
        out << "export in to static matrix failed with " << gbl_num_errs
            << " error" << (gbl_num_errs != 1 ? "s" : "") << "!\n";
      }
      gathervPrint(out, errStrm.str (), *comm);
      return; // no point in continuing
    }

    loc_success = success ? 1 : 0;
    Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_MIN, loc_success, outArg(gbl_success));
    TEST_EQUALITY_CONST(gbl_success, 1 );
    if (gbl_success != 1) {
      if (my_rank == 0) {
        out << "export in to static matrix comparison claims zero errors, "
          "but success is false on at least one process!\n";
      }
      gathervPrint (out, errStrm.str (), *comm);
      return; // no point in continuing
    }
  }
}

#define UNIT_TEST_GROUP( SC, LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, ImportToStaticMatrix, SC, LO, GO, NT) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, ImportToStaticMatrixLocal, SC, LO, GO, NT)

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_SLGN(UNIT_TEST_GROUP)

} // namespace (anonymous)

int
main (int argc, char* argv[])
{
  // Initialize MPI (if enabled) before initializing Kokkos.  This
  // lets MPI control things like pinning processes to sockets.
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  const int errCode =
    Teuchos::UnitTestRepository::runUnitTestsFromMain (argc, argv);
  return errCode;
}
