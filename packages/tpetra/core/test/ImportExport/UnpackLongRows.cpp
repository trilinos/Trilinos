// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <numeric>

#include "Kokkos_Core.hpp"

#include "Tpetra_TestingUtilities.hpp"
#include "TpetraCore_ETIHelperMacros.h"

#include "Tpetra_Map.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Distributor.hpp"

#include "Teuchos_CommHelpers.hpp"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_PerformanceMonitorBase.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_YamlParameterListCoreHelpers.hpp>

#define INFO(X) std::cout << "==> info [" << comm->getRank() << "/" << comm->getSize() << "]: " << X << std::flush

namespace { // anonymous

// Create a crs matrix with a Poisson-like structure and extra dense rows at its
// end. Add column entries to each row so that the matrix remains symmetric.
//
// The final matrix will have the following form
//
//    4 -1  0 -1  0    ...      0 1 ... 1
//   -1  4 -1  0 -1  0   ...    0 1 ... 1
//    0 -1  4 -1  0 -1 0  ...   0 1 ... 1
//   -1  0 -1  4 -1  0 -1 0 ... 0 1 ... 1
//    .        .                .    .
//    .          .              .    .
//    .            .            .    .
//    0  ... -1  0 -1  4 -1  0 -1 1 ... 1
//    0   ...   -1  0 -1  4 -1  0 1 ... 1
//    0     ...    -1  0 -1  4 -1 1 ... 1
//    0       ...     -1  0 -1  4 1 ... 1
//    1  1          ...                 1
//    .              .                  .
//    .              .                  .
//    1  1          ...                 1
//
// The final matrix is distributed evenly over processors in the communicator,
// but constructed in an arbitrary way such that each rank contributes to an
// overlap region on the rank's boundary and to the dense region.
//
// This test is constructed as such to isolate slow downs in export operations
// experienced by Aria when constructing a linear system that has a Poisson
// structure with additional dense rows.
template<class graph_type>
void
generate_graphs(
  Teuchos::RCP<const Teuchos::Comm<int>> const &comm,
  Teuchos::RCP<graph_type>& owned,
  Teuchos::RCP<graph_type>& shared,
  int const rows_per_rank,
  int const overlap,
  int const dense_rows
)
{

  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::Array;
  using map_type = typename graph_type::map_type;
  using GO = typename graph_type::global_ordinal_type;

  auto const rank = comm->getRank();
  auto const procs = comm->getSize();

  using gsize_t = Tpetra::global_size_t;
  const gsize_t global_rows = rows_per_rank * procs + dense_rows;

  INFO("GENERATING GRAPHS\n");
  // one-to-one map for entries on my rank
  RCP<map_type> owned_map;
  {
    Array<GO> indices((rank != procs - 1) ? rows_per_rank : rows_per_rank + dense_rows);
    GO start = rank * rows_per_rank;
    std::iota(indices.begin(), indices.end(), start);
    owned_map = rcp(new map_type(global_rows, indices(), 0, comm));
  }

  // overlapping map for shared entries
  RCP<map_type> shared_map;
  {
    Array<GO> indices;
    if (rank == 0)
    {
      indices.resize(overlap + dense_rows);
      // overlap in to rank 1
      std::iota(indices.begin(), indices.begin()+overlap, GO(rows_per_rank));
      // dense rows at end of matrix
      std::iota(indices.begin()+overlap, indices.end(), GO(procs * rows_per_rank));
    }
    else if (rank == procs - 1)
    {
      indices.resize(overlap);
      // overlap in to rank (procs - 1)
      std::iota(indices.begin(), indices.end(), GO(rank * rows_per_rank - overlap));
    }
    else
    {
      indices.resize(2 * overlap + dense_rows);
      // overlap in to rank (rank - 1)
      std::iota(indices.begin(), indices.begin()+overlap, GO(rank * rows_per_rank - overlap));
      // overlap in to rank (rank + 1)
      std::iota(indices.begin()+overlap, indices.begin()+2*overlap, GO((rank+1) * rows_per_rank));
      // dense rows at end of matrix
      std::iota(indices.begin()+2*overlap, indices.end(), GO(procs * rows_per_rank));
    }
    auto invalid = Teuchos::OrdinalTraits<gsize_t>::invalid();
    shared_map = rcp(new map_type(invalid, indices(), 0, comm));
  }

  owned = rcp(new graph_type(owned_map, rows_per_rank + dense_rows));
  shared = rcp(new graph_type(shared_map, rows_per_rank + dense_rows));

  {
    using Teuchos::tuple;
    auto rows_to_fill = GO(rows_per_rank + overlap);
    if (rank > 0 && rank < procs - 1) rows_to_fill += overlap;
    auto start = (rank == 0) ? GO(0) : GO(rank * rows_per_rank - overlap);
    for (GO row=start; row<start+rows_to_fill; row++)
    {
      Array<GO> columns;

      if (row == 0)
      {
        // [4, -1, 0, -1]
        auto my_cols = tuple(row, row + 1, row + 3);
        columns.assign(my_cols.begin(), my_cols.end());
      }
      else if (row == 1 || row == 2)
      {
        // 1: [-1, 4, -1, 0, -1]
        // 2: [0, -1, 4, -1, 0, -1]
        auto my_cols = tuple(row - 1, row, row + 1, row + 3);
        columns.assign(my_cols.begin(), my_cols.end());
      }
      else if (
        gsize_t(row) == global_rows - 3 - dense_rows ||
        gsize_t(row) == global_rows - 2 - dense_rows
      )
      {
        // -3: [-1, 0, -1, 4, -1, 0]
        // -2:    [-1, 0, -1, 4, -1]
        auto my_cols = tuple(row - 3, row - 1, row, row + 1);
        columns.assign(my_cols.begin(), my_cols.end());
      }
      else if (gsize_t(row) == global_rows - 1 - dense_rows)
      {
        // [-1, 0, -1, 4]
        auto my_cols = tuple(row - 3, row - 1, row);
        columns.assign(my_cols.begin(), my_cols.end());
      }
      else
      {
        // [-1, 0, -1, 4, -1, 0, -1]
        auto my_cols = tuple(row - 3, row - 1, row, row + 1, row + 3);
        columns.assign(my_cols.begin(), my_cols.end());
      }

      // Fill in columns at end of row associated with the extra dense rows to
      // assure symmetry of the final matrix
      for (int i=0; i<dense_rows; i++)
      {
        columns.push_back(global_rows - dense_rows + i);
      }

      auto start_this_rank = GO(rank * rows_per_rank);
      auto end_this_rank = GO(start_this_rank + rows_per_rank);
      if (rank == procs - 1) end_this_rank += dense_rows;

      if (owned_map->isNodeGlobalElement(row))
      {
        owned->insertGlobalIndices(row, columns());
      }
      else if (shared_map->isNodeGlobalElement(row))
      {
        shared->insertGlobalIndices(row, columns());
      }
      else
      {
        TEUCHOS_TEST_FOR_EXCEPTION(
          true,
          std::logic_error,
          "Row " << row << " is not owned by anyone!"
        );
      }
    }
  }

  {
    for (int i=0; i<dense_rows; i++)
    {
      // fill in the dense rows with my contribution
      auto row = global_rows - dense_rows + i;
      auto n = rows_per_rank;
      if (rank == procs - 1) n += dense_rows;
      Array<GO> columns(n);
      std::iota(columns.begin(), columns.end(), GO(rank * rows_per_rank));
      if (rank == procs - 1)
      {
        TEUCHOS_TEST_FOR_EXCEPTION(
          !owned_map->isNodeGlobalElement(row),
          std::logic_error,
          "==> Error [" << rank << "/" << procs << "]: " <<
          "the global row " << row << " is not in this owned map\n" <<
          owned_map->getLocalElementList()
        );
        owned->insertGlobalIndices(row, columns());
      }
      else
      {
        TEUCHOS_TEST_FOR_EXCEPTION(
          !shared_map->isNodeGlobalElement(row),
          std::logic_error,
          "==> Error [" << rank << "/" << procs << "]: " <<
          "the global row " << row << " is not in this shared map\n" <<
          shared_map->getLocalElementList()
        );
        shared->insertGlobalIndices(row, columns());
      }
    }
  }
  shared->fillComplete();
  comm->barrier();

  // We've created a sparse matrix containing owned indices and another
  // containing shared. The owned matrix has one-to-one map, while the shared
  // matrix is overlapping. Here, export the entries from the shared matrix in
  // to the owned matrix.
  //
  // Since only the target Map is one-to-one, we have to use an Export.
  {
    using export_type = typename graph_type::export_type;
    export_type exporter(shared_map, owned_map);
    comm->barrier();

    owned->doExport(*shared, exporter, Tpetra::ADD);
  }
  owned->fillComplete();

  INFO("GENERATING GRAPHS DONE\n");
}


template<class matrix_type>
Teuchos::RCP<const matrix_type>
generate_matrix(
  Teuchos::RCP<const Teuchos::Comm<int>> const &comm,
  Teuchos::RCP<const typename matrix_type::crs_graph_type> const &g_owned,
  Teuchos::RCP<const typename matrix_type::crs_graph_type> const &g_shared,
  int const rows_per_rank,
  int const overlap,
  int const dense_rows
)
{

  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::Time;
  using Teuchos::Array;
  using Teuchos::TimeMonitor;
  using real = typename matrix_type::scalar_type;
  using GO = typename matrix_type::global_ordinal_type;
  using LO = typename matrix_type::local_ordinal_type;

  auto const rank = comm->getRank();
  auto const procs = comm->getSize();

  INFO("GENERATING MATRIX\n");
  using gsize_t = Tpetra::global_size_t;
  const gsize_t global_rows = rows_per_rank * procs + dense_rows;

  const real mone = static_cast<real>(-1.0);
  const real one = static_cast<real>(1.0);
  const real four = static_cast<real>(4.0);

  // one-to-one map for entries on my rank
  auto owned_map = g_owned->getRowMap();
  auto mtx_owned = rcp(new matrix_type(g_owned));

  auto shared_map = g_shared->getRowMap();
  auto mtx_shared = rcp(new matrix_type(g_shared));
  INFO("MATRIX CREATED\n");

  {
    INFO("FILLING MATRIX\n");
    RCP<Time> stage = TimeMonitor::getNewCounter("UnpackLongRows::fill");
    TimeMonitor tm(*stage);

    using Teuchos::tuple;
    auto rows_to_fill = GO(rows_per_rank + overlap);
    if (rank > 0 && rank < procs - 1) rows_to_fill += overlap;
    auto start = (rank == 0) ? GO(0) : GO(rank * rows_per_rank - overlap);
    for (GO row=start; row<start+rows_to_fill; row++)
    {
      Array<GO> columns;
      Array<real> values;

      if (row == 0)
      {
        // [4, -1, 0, -1]
        auto my_cols = tuple(row, row + 1, row + 3);
        auto my_vals = tuple(four, mone, mone);
        columns.assign(my_cols.begin(), my_cols.end());
        values.assign(my_vals.begin(), my_vals.end());
      }
      else if (row == 1 || row == 2)
      {
        // 1: [-1, 4, -1, 0, -1]
        // 2: [0, -1, 4, -1, 0, -1]
        auto my_cols = tuple(row - 1, row, row + 1, row + 3);
        auto my_vals = tuple(mone, four, mone, mone);
        columns.assign(my_cols.begin(), my_cols.end());
        values.assign(my_vals.begin(), my_vals.end());
      }
      else if (
        gsize_t(row) == global_rows - 3 - dense_rows ||
        gsize_t(row) == global_rows - 2 - dense_rows
      )
      {
        // -3: [-1, 0, -1, 4, -1, 0]
        // -2:    [-1, 0, -1, 4, -1]
        auto my_cols = tuple(row - 3, row - 1, row, row + 1);
        auto my_vals = tuple(mone, mone, four, mone);
        columns.assign(my_cols.begin(), my_cols.end());
        values.assign(my_vals.begin(), my_vals.end());
      }
      else if (gsize_t(row) == global_rows - 1 - dense_rows)
      {
        // [-1, 0, -1, 4]
        auto my_cols = tuple(row - 3, row - 1, row);
        auto my_vals = tuple(mone, mone, four);
        columns.assign(my_cols.begin(), my_cols.end());
        values.assign(my_vals.begin(), my_vals.end());
      }
      else
      {
        // [-1, 0, -1, 4, -1, 0, -1]
        auto my_cols = tuple(row - 3, row - 1, row, row + 1, row + 3);
        auto my_vals = tuple(mone, mone, four, mone, mone);
        columns.assign(my_cols.begin(), my_cols.end());
        values.assign(my_vals.begin(), my_vals.end());
      }

      // Fill in columns at end of row associated with the extra dense rows to
      // assure symmetry of the final matrix
      for (int i=0; i<dense_rows; i++)
      {
        real value = one;
        if (gsize_t(row) >= global_rows - dense_rows) value /= real(procs);
        columns.push_back(global_rows - dense_rows + i);
        values.push_back(value);
      }

      auto start_this_rank = GO(rank * rows_per_rank);
      auto end_this_rank = GO(start_this_rank + rows_per_rank);
      if (rank == procs - 1) end_this_rank += dense_rows;

      if (
        (rank != procs - 1 && row > end_this_rank - overlap) ||
        (rank != 0 && row < start_this_rank + overlap)
      )
      {
        // Scale values in overlap region
        std::transform(
          values.begin(),
          values.end(),
          values.begin(),
          [](real & x){ return .5 * x; }
        );
      }

      Array<LO> local_columns(columns.size());
      if (owned_map->isNodeGlobalElement(row))
      {
        auto local_row = owned_map->getLocalElement(row);
        auto colmap = mtx_owned->getColMap();
        std::transform(
          columns.begin(),
          columns.end(),
          local_columns.begin(),
          [&](GO & col){ return colmap->getLocalElement(col); }
        );
        mtx_owned->sumIntoLocalValues(local_row, local_columns(), values());
      }
      else if (shared_map->isNodeGlobalElement(row))
      {
        auto local_row = shared_map->getLocalElement(row);
        auto colmap = mtx_shared->getColMap();
        std::transform(
          columns.begin(),
          columns.end(),
          local_columns.begin(),
          [&](GO & col){ return colmap->getLocalElement(col); }
        );
        mtx_shared->sumIntoLocalValues(local_row, local_columns(), values());
      }
      else
      {
        TEUCHOS_TEST_FOR_EXCEPTION(
          true,
          std::logic_error,
          "Row " << row << " is not owned by anyone!"
        );
      }
    }
    INFO("MATRIX FILLED\n");
  }

  {
    INFO("FILLING DENSE ROWS\n");
    RCP<Time> stage = TimeMonitor::getNewCounter("UnpackLongRows::fill_dense_rows");
    TimeMonitor tm(*stage);

    for (int i=0; i<dense_rows; i++)
    {
      // fill in the dense rows with my contribution
      auto row = global_rows - dense_rows + i;
      auto n = rows_per_rank;
      if (rank == procs - 1) n += dense_rows;
      Array<GO> columns(n);
      Array<LO> local_columns(n);
      std::iota(columns.begin(), columns.end(), GO(rank * rows_per_rank));
      Array<real> values(n, one);
      if (rank == procs - 1)
      {
        TEUCHOS_TEST_FOR_EXCEPTION(
          !owned_map->isNodeGlobalElement(row),
          std::logic_error,
          "==> Error [" << rank << "/" << procs << "]: " <<
          "the global row " << row << " is not in this owned map\n" <<
          owned_map->getLocalElementList()
        );
        auto local_row = owned_map->getLocalElement(row);
        auto colmap = mtx_owned->getColMap();
        std::transform(
          columns.begin(),
          columns.end(),
          local_columns.begin(),
          [&](GO & col){ return colmap->getLocalElement(col); }
        );
        mtx_owned->sumIntoLocalValues(local_row, local_columns(), values());
      }
      else
      {
        TEUCHOS_TEST_FOR_EXCEPTION(
          !shared_map->isNodeGlobalElement(row),
          std::logic_error,
          "==> Error [" << rank << "/" << procs << "]: " <<
          "the global row " << row << " is not in this shared map\n" <<
          shared_map->getLocalElementList()
        );
        auto local_row = shared_map->getLocalElement(row);
        auto colmap = mtx_shared->getColMap();
        std::transform(
          columns.begin(),
          columns.end(),
          local_columns.begin(),
          [&](GO & col){ return colmap->getLocalElement(col); }
        );
        mtx_shared->sumIntoLocalValues(local_row, local_columns(), values());
      }
    }
    INFO("DENSE ROWS FILLED\n");
  }
  comm->barrier();

  // We've created a sparse matrix containing owned values and another
  // containing shared. The owned matrix has one-to-one map, while the shared
  // matrix is overlapping. Here, export the entries from the shared matrix in
  // to the owned matrix.
  //
  // Since only the target Map is one-to-one, we have to use an Export.
  {
    INFO("EXPORTING\n");
    RCP<Time> stage = TimeMonitor::getNewCounter("UnpackLongRows::export");
    TimeMonitor tm(*stage);
    using export_type = typename matrix_type::export_type;
    export_type exporter(shared_map, owned_map);
    comm->barrier();

    mtx_owned->doExport(*mtx_shared, exporter, Tpetra::ADD);
    INFO("EXPORT DONE\n");
  }
  {
    INFO("CALLING FILLCOMPLETE\n");
    RCP<Time> stage = TimeMonitor::getNewCounter("UnpackLongRows::fill_complete");
    TimeMonitor tm(*stage);
    mtx_owned->fillComplete();
    INFO("FILLCOMPLETE DONE\n");
  }

  INFO("GENERATING MATRIX DONE\n");
  return mtx_owned;
}

Teuchos::ParameterList
get_timer_stats(const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
{
  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::ParameterList;
  using Teuchos::TimeMonitor;
  typedef std::vector<std::string>::size_type size_type;
  typedef Teuchos::stat_map_type stat_map_type;

  // bool alwaysWriteLocal = false;
  // bool writeGlobalStats = true;
  // bool writeZeroTimers = false;

  stat_map_type stat_data;
  std::vector<std::string> stat_names;
  TimeMonitor::computeGlobalTimerStatistics(stat_data, stat_names, comm.ptr(), Teuchos::Union);

  ParameterList p1;
  for (auto it=stat_data.begin(); it!=stat_data.end(); ++it) {
    ParameterList px;
    ParameterList counts;
    ParameterList times;
    const std::vector<std::pair<double, double> >& cur_data = it->second;
    for (size_type ix = 0; ix < cur_data.size(); ++ix) {
      times.set(stat_names[ix], cur_data[ix].first);
      counts.set(stat_names[ix], static_cast<int>(cur_data[ix].second));
    }
    px.set("Total times", times);
    px.set("Call counts", counts);
    p1.set(it->first, px);
  }

  return p1;
}

template<class local_ordinal, class global_ordinal, class node, class buffer_device>
void
time_single_row_unpack()
{
  std::vector<size_t> test_row_lengths{10, 100, 1000, 10000};
  auto sizeof_int = sizeof(int);
  for (auto test_row_length : test_row_lengths)
  {
    {
      // Parallel kokkos unpack
      Kokkos::View<int*, buffer_device> a("unpacked", test_row_length);
      Kokkos::View<char*, buffer_device> c("packed", test_row_length * sizeof_int);
      std::ostringstream os;
      os << "Standalone test: parallel: " << test_row_length;
      Teuchos::RCP<Teuchos::Time> st = Teuchos::TimeMonitor::getNewCounter(os.str());
      Teuchos::TimeMonitor tm(*st);
      Kokkos::parallel_for(
        test_row_length,
        KOKKOS_LAMBDA(const size_t i)
        {
          auto start_a = i;
          auto start_c = i * sizeof_int;
          memcpy(c.data() + start_c, a.data() + start_a, sizeof_int);
        }
      );
    }

    {
      // memcpy
      Kokkos::View<int*, buffer_device> a("unpacked", test_row_length);
      Kokkos::View<char*, buffer_device> c("packed", test_row_length * sizeof_int);
      std::ostringstream os;
      os << "Standalone test: one row: " << test_row_length;
      Teuchos::RCP<Teuchos::Time> st = Teuchos::TimeMonitor::getNewCounter(os.str());
      Teuchos::TimeMonitor tm(*st);
      Kokkos::parallel_for(
        1,
        KOKKOS_LAMBDA(const size_t i)
        {
          memcpy(c.data(), a.data(), test_row_length * sizeof_int);
        }
      );
    }
  }
}

} // namespace (anonymous)

int
main(int argc, char* argv[])
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  using local_ordinal = Tpetra::Map<>::local_ordinal_type;
  using global_ordinal = Tpetra::Map<>::global_ordinal_type;
  using node_type = Tpetra::Map<>::node_type;
#if defined (HAVE_TPETRA_INST_DOUBLE)
  using real = double;
#elif defined (HAVE_TPETRA_INST_FLOAT)
  using real = float;
#else
#  error "Tpetra: Must enable at least one Scalar type in {double, float} in order to build this test."
#endif

  Tpetra::Details::Behavior::disable_timing();

  RCP<Teuchos::oblackholestream> blackhole = rcp (new Teuchos::oblackholestream);
  Teuchos::GlobalMPISession mpi_session (&argc, &argv, blackhole.getRawPtr());
  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();
  const int rank = comm->getRank();
  const int procs = comm->getSize();

  RCP<Teuchos::FancyOStream> proc_zero_out = (rank == 0) ?
    Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout)) :
    Teuchos::getFancyOStream(blackhole);

  //
  // Default values of command-line options.
  //
  int rows_per_rank = 5000;
  int overlap = 5;
  int dense_rows = 5;
  Teuchos::CommandLineProcessor parser(false, true);
  {
    parser.setOption("rows-per-rank", &rows_per_rank, "Rows per rank in owned matrix");
    parser.setOption("overlap", &overlap, "Number of overlapping rows in overlapping matrix");
    parser.setOption("dense-rows", &dense_rows, "Number of dense rows in owned matrix");

    auto p = parser.parse(argc, argv);
    if (p == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
    {
      return 0;
    }
    else if (p != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
    {
      *proc_zero_out << "End Result: TEST FAILED" << std::endl;
      return EXIT_FAILURE;
    }
  }

  using graph_type = Tpetra::CrsGraph<local_ordinal, global_ordinal, node_type>;
  RCP<graph_type> g_owned;
  RCP<graph_type> g_shared;
  generate_graphs(comm, g_owned, g_shared, rows_per_rank, overlap, dense_rows);

  Tpetra::Details::Behavior::enable_timing();

  using matrix_type = Tpetra::CrsMatrix<real, local_ordinal, global_ordinal, node_type>;
  RCP<matrix_type> m_owned = rcp(new matrix_type(g_owned));
  RCP<matrix_type> m_shared = rcp(new matrix_type(g_shared));

  {
    RCP<Teuchos::Time> g_stage = Teuchos::TimeMonitor::getNewCounter("UnpackLongRows::run");
    Teuchos::TimeMonitor g_tm(*g_stage);
    auto matrix = generate_matrix<matrix_type>(
      comm, g_owned, g_shared, rows_per_rank, overlap, dense_rows
    );

    using dist_object_type = Tpetra::DistObject<char, local_ordinal, global_ordinal, node_type>;
    using bdt = typename dist_object_type::buffer_device_type;
    time_single_row_unpack<local_ordinal, global_ordinal, node_type, bdt>();
  }

  Teuchos::ParameterList p0("UnpackLongRows");
  p0.set("rows_per_rank", rows_per_rank);
  p0.set("overlap", overlap);
  p0.set("dense_rows", dense_rows);
  p0.set("processors", procs);

  auto pl = get_timer_stats(comm);
  p0.set("Timing", pl);

  std::ostringstream f;
  f << "UnpackLongRows_"
    << procs << "_"
    << rows_per_rank << "_"
    << overlap << "_"
    << dense_rows << ".xml";
  Teuchos::writeParameterListToXmlFile(p0, f.str());
  Teuchos::TimeMonitor::clearCounters();

  // Always successful - we just want to time the thing
  bool success = true;
  *proc_zero_out << "End Result: TEST " << (success ? "PASSED" : "FAILED") << std::endl;
  return success ? EXIT_SUCCESS : EXIT_FAILURE;

}

#undef INFO
