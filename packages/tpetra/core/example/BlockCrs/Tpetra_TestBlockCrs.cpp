// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Core.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_RCP.hpp"

#include "Tpetra_TestBlockCrsMeshDatabase.hpp"

#include <memory>
#include <type_traits>

namespace { // (anonymous)

  int cubeRootRoundedToInt (const double input) {
    constexpr double power = 1.0 / 3.0;
    return static_cast<int> (std::pow (input, power));
  }

  int cubeRootLowerBound (const int P) { // P >= 1
    double input = static_cast<double> (P);
    int estimate = cubeRootRoundedToInt (input);
    int cube = estimate * estimate * estimate;

    if (cube > P && input >= 1.0 && estimate >= 1.0) {
      // If rounding resulted in too big of a result, reduce input by
      // one until result is a lower bound.  do-while avoids redundant
      // condition check.
      do {
        input = input - 1.0;
        estimate = cubeRootRoundedToInt (input);
        cube = estimate * estimate * estimate;
      }
      while (cube > P && input >= 1.0 && estimate >= 1.0);
    }
    else if (cube < P) {
      // See if we can increment the result to get it closer.
      do {
        input = input + 1.0;
        estimate = cubeRootRoundedToInt (input);
        cube = estimate * estimate * estimate;
      }
      while (cube < P);

      if (cube > P && estimate > 1) {
        estimate = estimate - 1;
      }
    }

    return estimate; // estimate*estimate*estimate <= P
  }

#if 0
  void
  runTest (const Teuchos::RCP<const Teuchos::Comm<int> >& commTest,
           const LO num_global_elements_i,
           const LO num_global_elements_j,
           const LO num_global_elements_k,
           const process_rank_type num_procs_i,
           const process_rank_type num_procs_j,
           const process_rank_type num_procs_k,
           const LO blocksize,
           const LO nrhs,
           const LO repeat,
           const bool debug,
           const bool verbose)
  {

  }
#endif // 0

} // namespace (anonymous)

using namespace BlockCrsTest;

int main (int argc, char *argv[])
{
  using Teuchos::RCP;
  using Teuchos::TimeMonitor;
  using std::endl;
  typedef int process_rank_type;

  auto out = Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cout));

  // Initialize MPI and Kokkos
  Tpetra::initialize (&argc, &argv);
  {
    auto commWorld = Tpetra::getDefaultComm ();

    // Don't try to read environment variables until after MPI_Init has
    // been called.  I'm not sure whether all the processes can see them
    // before.
    const bool debug = Tpetra::Details::Behavior::debug ();
    const bool verbose = Tpetra::Details::Behavior::verbose ();
    std::unique_ptr<std::string> debugPrefix;
    if (debug) {
      std::ostringstream os;
      os << "Proc " << commWorld->getRank () << ": ";
      debugPrefix = std::unique_ptr<std::string> (new std::string (os.str ()));
    }

    // By default, distribute the problem over a defaultProcGridDim^3 cube.
    const int defaultProcGridDim = cubeRootLowerBound (commWorld->getSize ());
    if (debug) {
      std::ostringstream os;
      os << *debugPrefix << "defaultProcGridDim: " << defaultProcGridDim << endl;
      std::cerr << os.str ();
    }

    // Command-line input
    process_rank_type num_procs_i = defaultProcGridDim;
    process_rank_type num_procs_j = defaultProcGridDim;
    process_rank_type num_procs_k = defaultProcGridDim;

    constexpr LO defaultNumElementsMultiplierAlongEachDim = 10;

    LO num_global_elements_i = defaultNumElementsMultiplierAlongEachDim * defaultProcGridDim;
    LO num_global_elements_j = defaultNumElementsMultiplierAlongEachDim * defaultProcGridDim;
    LO num_global_elements_k = defaultNumElementsMultiplierAlongEachDim * defaultProcGridDim;

    LO blocksize = 5, nrhs = 1, repeat = 100;
    bool dump_sparse_matrix = false;

    Teuchos::CommandLineProcessor clp (false);
    clp.setDocString ("Tpetra::BlockCrsMatrix performance test using 3-D 7-point stencil.\n");
    clp.setOption ("num-elements-i", &num_global_elements_i, "Number of cells in the I dimension.");
    clp.setOption ("num-elements-j", &num_global_elements_j, "Number of cells in the J dimension.");
    clp.setOption ("num-elements-k", &num_global_elements_k, "Number of cells in the K dimension.");
    clp.setOption ("num-procs-i", &num_procs_i,
                   "Process grid of (npi,npj,npk); npi*npj*npk <= number of MPI ranks.");
    clp.setOption ("num-procs-j", &num_procs_j,
                   "Process grid of (npi,npj,npk); npi*npj*npk <= number of MPI ranks.");
    clp.setOption ("num-procs-k", &num_procs_k,
                   "Process grid of (npi,npj,npk); npi*npj*npk <= number of MPI ranks.");
    clp.setOption ("blocksize", &blocksize,
                   "Block size. The # of DOFs coupled in a multiphysics flow problem.");
    clp.setOption ("nrhs", &nrhs,
                   "Number of right-hand sides for which to solve.");
    clp.setOption ("repeat", &repeat,
                   "Number of iterations of matvec operations to measure performance.");
    clp.setOption ("dump-sparse-matrix", "dont-dump-sparse-matrix", &dump_sparse_matrix,
                   "If true, dump the test sparse matrix to a MatrixMarket file "
                   "in the current directory.  This is a debugging option and may "
                   "take a lot of disk space and time.");

    {
      bool returnEarly = false;
      bool returnSuccess = true;

      auto r_parse = clp.parse(argc, argv, commWorld->getRank() == 0 ? &std::cout : NULL);
      if (r_parse == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        returnEarly = true;
        returnSuccess = true; // printing help doesn't mean the test failed
      }
      else if (r_parse != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
        returnEarly = true;
        returnSuccess = false;
        if (commWorld->getRank () == 0) {
          std::cout << "Failed to parse command-line arguments." << endl;
        }
      }
      else if (num_procs_i*num_procs_j*num_procs_k > commWorld->getSize ()) {
        returnEarly = true;
        returnSuccess = false;
        if (commWorld->getRank () == 0) {
          std::cout << "Invalid process grid: "
                    << num_procs_i << " x " << num_procs_j << " x " << num_procs_k << " > "
                    << " number of MPI processes "<< commWorld->getSize () << endl;
        }
      }

      if (! returnSuccess && commWorld->getRank () == 0) {
        std::cout << "End Result: TEST FAILED" << endl;
      }
      if (returnEarly) {
        Tpetra::finalize();
        return returnSuccess ? EXIT_SUCCESS : EXIT_FAILURE;
      }
    }

    // Create a communicator just for the test, that excludes processes
    // outside the process grid.
    const bool myProcessIsIncludedInTheTest =
      commWorld->getRank () < num_procs_i*num_procs_j*num_procs_k;
    auto commTest = commWorld->split (myProcessIsIncludedInTheTest ? 0 : 1,
                                      commWorld->getRank ());
    if (! myProcessIsIncludedInTheTest) {
      Tpetra::finalize ();
      return EXIT_SUCCESS;
    }

    ////////////////////////////////////////////////////////////////
    // At this point, the only processes left are those in commTest.
    ////////////////////////////////////////////////////////////////
    const bool myProcessPrintsDuringTheTest = commTest->getRank () == 0;

    // When multiple GPUs are used, we would like to see a different mpi 
    // uses a different device. An easy check is device configuration for
    // 4 mpi nodes as the typical number of GPUs per node is 4.
    if (commTest->getRank () < 4) {
      exec_space{}.print_configuration (std::cout, false);
    }

    if (debug) {
      std::ostringstream os;
      os << *debugPrefix << "Make MeshDatabase" << endl;
      std::cerr << os.str ();
    }
    MeshDatabase mesh (commTest,
                       num_global_elements_i,
                       num_global_elements_j,
                       num_global_elements_k,
                       num_procs_i,
                       num_procs_j,
                       num_procs_k);

    const auto sb = mesh.getStructuredBlock();
    const auto part = mesh.getStructuredBlockPart();

    const auto num_owned_elements = mesh.getNumOwnedElements();
    const auto num_owned_and_remote_elements = mesh.getNumElements();

    if (verbose) {
      mesh.print(std::cout);
    }

    if (debug) {
      std::ostringstream os;
      os << *debugPrefix << "Start global timer" << endl;
      std::cerr << os.str ();
    }
    {
      TimeMonitor timerGlobal(*TimeMonitor::getNewTimer("X) Global"));

      // Build Tpetra Maps
      // -----------------

      const Tpetra::global_size_t TpetraComputeGlobalNumElements
        = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();

      // - all internal views are allocated on device; mirror as mesh database is constructed on host
      const auto mesh_gids_host = mesh.getElementGlobalIDs();
      const auto mesh_gids =
        Kokkos::create_mirror_view_and_copy (mem_space(), mesh_gids_host);

      // for convenience, separate the access to owned and remote gids
      const auto owned_gids =
        Kokkos::subview (mesh_gids,
                         local_ordinal_range_type (0, num_owned_elements));
      const auto remote_gids =
        Kokkos::subview (mesh_gids,
                         local_ordinal_range_type (num_owned_elements,
                                                   num_owned_and_remote_elements));

      if (debug) {
        std::ostringstream os;
        os << *debugPrefix << "Make row and column Maps" << endl;
        std::cerr << os.str ();
      }
      RCP<const map_type> row_map (new map_type (TpetraComputeGlobalNumElements, owned_gids, 0, commTest));
      RCP<const map_type> col_map (new map_type (TpetraComputeGlobalNumElements, mesh_gids, 0, commTest));

      if (verbose) {
        row_map->describe(*out);
        col_map->describe(*out);
      }

      if (debug) {
        std::ostringstream os;
        os << *debugPrefix << "Make local graph" << endl;
        std::cerr << os.str ();
      }
      // Graph Construction
      // ------------------
      // local graph is constructed on device space
      typedef tpetra_crs_graph_type::local_graph_device_type local_graph_device_type;
      typedef local_graph_device_type::row_map_type::non_const_type rowptr_view_type;
      typedef typename local_graph_device_type::entries_type colidx_view_type;

      rowptr_view_type rowptr;
      colidx_view_type colidx;
      {
        TimeMonitor timerLocalGraphConstruction(*TimeMonitor::getNewTimer("0) LocalGraphConstruction"));

        rowptr = rowptr_view_type("rowptr", num_owned_elements + 1);

        local_ordinal_range_type owned_range_i, owned_range_j, owned_range_k;
        local_ordinal_range_type remote_range_i, remote_range_j, remote_range_k;

        part.getOwnedRange(owned_range_i, owned_range_j, owned_range_k);
        part.getRemoteRange(sb, remote_range_i, remote_range_j, remote_range_k);

        // count # of nonzeros per row
        {
          LocalGraphConstruction local_graph_construction
            (sb, owned_gids,
             remote_range_i, remote_range_j, remote_range_k,
             rowptr);
          local_graph_construction.run();
        }

        // the last entry of rowptr is the total number of nonzeros in the local graph
        // mirror to host to use the information in constructing colidx
        auto nnz = Kokkos::subview(rowptr, num_owned_elements);
        const auto nnz_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), nnz);

        // allocate colidx
        colidx = colidx_view_type("colidx", nnz_host());

        // fill
        {
          LocalGraphFill local_graph_fill(sb, owned_gids, remote_gids,
                                          owned_range_i, owned_range_j, owned_range_k,
                                          remote_range_i, remote_range_j, remote_range_k,
                                          rowptr, colidx);
          local_graph_fill.run();
        }

      } // end local graph timer

      if (debug) {
        std::ostringstream os;
        os << *debugPrefix << "Make global graph" << endl;
        std::cerr << os.str ();
      }
      // Call fillComplete on the bcrs_graph to finalize it
      RCP<tpetra_crs_graph_type> bcrs_graph;
      {
        TimeMonitor timerGlobalGraphConstruction(*TimeMonitor::getNewTimer("1) GlobalGraphConstruction"));
        rowptr_view_type rowptr_tpetra = 
          rowptr_view_type(Kokkos::ViewAllocateWithoutInitializing("rowptr_tpetra"), rowptr.extent(0));
        colidx_view_type colidx_tpetra =
          colidx_view_type(Kokkos::ViewAllocateWithoutInitializing("colidx_tpetra"), colidx.extent(0));
        Kokkos::deep_copy(rowptr_tpetra, rowptr);
        Kokkos::deep_copy(colidx_tpetra, colidx);
        bcrs_graph = rcp(new tpetra_crs_graph_type(row_map, col_map, 
                                                   local_graph_device_type(colidx_tpetra, rowptr_tpetra),
                                                   Teuchos::null));
      } // end global graph timer

      if (verbose) {
        bcrs_graph->describe(*out, Teuchos::VERB_EXTREME);
      }

      if (debug) {
        std::ostringstream os;
        os << *debugPrefix << "Make BlockCrsMatrix" << endl;
        std::cerr << os.str ();
      }
      // Create BlockCrsMatrix
      RCP<tpetra_blockcrs_matrix_type> A_bcrs (new tpetra_blockcrs_matrix_type (*bcrs_graph, blocksize));

      if (debug) {
        std::ostringstream os;
        os << *debugPrefix << "Fill BlockCrsMatrix" << endl;
        std::cerr << os.str ();
      }
      typedef tpetra_blockcrs_matrix_type::little_block_type block_type;
      Kokkos::View<block_type*, device_type> blocks;
      {
        TimeMonitor timerLocalBlockCrsFill(*TimeMonitor::getNewTimer("2) LocalBlockCrsFill"));

        // Tpetra BlockCrsMatrix only has high level access functions
        // To fill this on device, we need an access to the meta data of blocks
        const auto rowptr_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), rowptr);
        const auto colidx_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), colidx);

        blocks = Kokkos::View<block_type*, device_type>("blocks", rowptr_host(num_owned_elements));

        const auto blocks_host = Kokkos::create_mirror_view(blocks);
        // This MUST run on host, since it invokes a host-only method,
        // getLocalBlock.  This means we must NOT use KOKKOS_LAMBDA,
        // since that would build the lambda for both host AND device.

        /// without UVM, the getLocalBlockDeviceNonConst cannot be called within the parallel for 
        /// even though it is host execution space as the method can involve kernel launch 
        /// for memory transfers.
        // Kokkos::parallel_for
        //   (Kokkos::RangePolicy<host_space, LO> (0, num_owned_elements),
        //    [&] (const LO row) {
        for (LO row=0;row<LO(num_owned_elements);++row) {
          const auto beg = rowptr_host(row);
          const auto end = rowptr_host(row+1);
          typedef typename std::remove_const<decltype (beg) >::type offset_type;
          for (offset_type loc = beg; loc < end; ++loc) {
            blocks_host(loc) = A_bcrs->getLocalBlockDeviceNonConst(row, colidx_host(loc));
          }
        }
        //   });

        Kokkos::deep_copy(blocks, blocks_host);

        Kokkos::parallel_for
          (Kokkos::RangePolicy<exec_space, LO> (0, num_owned_elements),
           KOKKOS_LAMBDA (const LO row) {
            const auto beg = rowptr(row);
            const auto end = rowptr(row+1);
            typedef typename std::remove_const<decltype (beg) >::type offset_type;
            for (offset_type loc = beg; loc < end; ++loc) {
              const GO gid_row = mesh_gids(row);
              const GO gid_col = mesh_gids(colidx(loc));

              LO i0, j0, k0, i1, j1, k1;
              sb.idx_to_ijk(gid_row, i0, j0, k0);
              sb.idx_to_ijk(gid_col, i1, j1, k1);

              const LO diff_i = i0 - i1;
              const LO diff_j = j0 - j1;
              const LO diff_k = j0 - k1;

              auto block = blocks(loc);
              for (LO l0=0;l0<blocksize;++l0)
                for (LO l1=0;l1<blocksize;++l1)
                  block(l0, l1) = get_block_crs_entry<value_type>(i0, j0, k0,
                                                                  diff_i, diff_j, diff_k,
                                                                  l0, l1);
            }
          });
      }

      {
        TimeMonitor timerBlockCrsFillComplete(*TimeMonitor::getNewTimer("3) BlockCrsMatrix FillComplete - currently do nothing"));
        // this function does not exist
        //A_bcrs->fillComplete();
      }

      if (verbose) {
        A_bcrs->describe(*out, Teuchos::VERB_EXTREME);
      }

      if (debug) {
        std::ostringstream os;
        os << *debugPrefix << "Make and fill MultiVector" << endl;
        std::cerr << os.str ();
      }
      // Create MultiVector
      RCP<tpetra_multivector_type> X = Teuchos::rcp(new tpetra_multivector_type(A_bcrs->getDomainMap(), nrhs));
      {
        TimeMonitor timerMultiVectorFill(*TimeMonitor::getNewTimer("4) MultiVectorFill"));

        auto value = X->getLocalViewDevice(Tpetra::Access::OverwriteAll);
        auto map = X->getMap()->getLocalMap();
        Kokkos::parallel_for
          (value.extent(0), KOKKOS_LAMBDA(const LO i) {
            const GO gid = map.getGlobalElement(i);
            for (LO j=0,jend=value.extent(1);j<jend;++j)
              value(i, j) = get_multi_vector_entry<value_type>(gid, j);
          });
      }

      if (verbose) {
        X->describe(*out, Teuchos::VERB_EXTREME);
      }

      if (debug) {
        std::ostringstream os;
        os << *debugPrefix << "Make B_bcrs MultiVector" << endl;
        std::cerr << os.str ();
      }
      RCP<tpetra_multivector_type> B_bcrs = Teuchos::rcp(new tpetra_multivector_type(A_bcrs->getRangeMap(),  nrhs));

      if (debug) {
        std::ostringstream os;
        os << *debugPrefix << "Benchmark BlockCrsMatrix sparse matrix-vector multiply" << endl;
        std::cerr << os.str ();
      }
      // matrix vector multiplication
      {
        for (LO iter=0;iter<repeat;++iter) {
          TimeMonitor timerBlockCrsApply(*TimeMonitor::getNewTimer("5) BlockCrs Apply"));
          A_bcrs->apply(*X, *B_bcrs);
        }
      }

      if (debug) {
        std::ostringstream os;
        os << *debugPrefix << "Convert BlockCrsMatrix to CrsMatrix" << endl;
        std::cerr << os.str ();
      }
      // direct conversion: block crs -> point crs
      RCP<tpetra_crs_matrix_type> A_crs;
      {
        TimeMonitor timerConvertBlockCrsToPointCrs(*TimeMonitor::getNewTimer("6) Conversion from BlockCrs to PointCrs"));

        // construct row map and column map for a point crs matrix
        // point-wise row map can be obtained from A_bcrs->getDomainMap().
        // A constructor exist for crs matrix with a local matrix and a row map.
        // see, Tpetra_CrsMatrix_decl.hpp, line 504
        //     CrsMatrix (const local_matrix_device_type& lclMatrix,
        //                const Teuchos::RCP<const map_type>& rowMap,
        //                const Teuchos::RCP<const map_type>& colMap = Teuchos::null,
        //                const Teuchos::RCP<const map_type>& domainMap = Teuchos::null,
        //                const Teuchos::RCP<const map_type>& rangeMap = Teuchos::null,
        //                const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);
        // However, this does not work with the following use case.
        //   A_crs = rcp(new tpetra_crs_matrix_type(local_matrix, A_bcrs->getDomainMap()));

        // here, we create pointwise row and column maps manually.
        decltype(mesh_gids) crs_gids("crs_gids", mesh_gids.extent(0)*blocksize);
        Kokkos::parallel_for
          (num_owned_and_remote_elements,
           KOKKOS_LAMBDA(const LO idx) {
            for (LO l=0;l<blocksize;++l)
              crs_gids(idx*blocksize+l) = mesh_gids(idx)*blocksize+l;
          });
        const auto owned_crs_gids = Kokkos::subview(crs_gids, local_ordinal_range_type(0, num_owned_elements*blocksize));

        RCP<const map_type> row_crs_map (new map_type (TpetraComputeGlobalNumElements,
                                                       owned_crs_gids, 0, commTest));
        RCP<const map_type> col_crs_map (new map_type (TpetraComputeGlobalNumElements,
                                                       crs_gids, 0, commTest));

        rowptr_view_type crs_rowptr = rowptr_view_type("crs_rowptr", num_owned_elements*blocksize+1);
        colidx_view_type crs_colidx = colidx_view_type("crs_colidx", colidx.extent(0)*blocksize*blocksize);
        typename tpetra_crs_matrix_type::local_matrix_device_type::values_type
          crs_values("crs_values", colidx.extent(0)*blocksize*blocksize);

        Kokkos::parallel_for
          (Kokkos::RangePolicy<exec_space, LO> (0, num_owned_elements),
           KOKKOS_LAMBDA (const LO &idx) {
            const GO nnz_per_block_row = rowptr(idx+1)-rowptr(idx); // FIXME could be LO if no duplicates
            const GO nnz_per_point_row = nnz_per_block_row*blocksize;
            const GO crs_rowptr_begin  = idx*blocksize;
            const GO crs_colidx_begin  = rowptr(idx)*blocksize*blocksize;

            for (LO i=0;i<(blocksize+1);++i) {
              crs_rowptr(crs_rowptr_begin+i) = crs_colidx_begin + i*nnz_per_point_row;
            }

            GO loc = crs_colidx_begin;
            // loop over the rows in a block
            for (LO l0=0;l0<blocksize;++l0) {
              // loop over the block row
              typedef typename std::decay<decltype (rowptr(idx)) >::type offset_type;

              for (offset_type jj = rowptr(idx); jj < rowptr(idx+1); ++jj) {
                const auto block = blocks(jj);
                // loop over the cols in a block
                const GO offset = colidx(jj)*blocksize;
                for (LO l1=0;l1<blocksize;++l1) {
                  crs_colidx(loc) = offset+l1;
                  crs_values(loc) = block(l0,l1);
                  ++loc;
                }
              }
            }
          });

        typename tpetra_crs_matrix_type::local_matrix_device_type
          local_matrix("local_crs_matrix",
                       num_owned_and_remote_elements*blocksize,
                       crs_values,
                       local_graph_device_type(crs_colidx, crs_rowptr));

        A_crs = rcp(new tpetra_crs_matrix_type(row_crs_map,
                                               col_crs_map,
                                               local_matrix,
                                               Teuchos::null));

      } // end conversion timer

      if (verbose) {
        A_crs->describe(*out, Teuchos::VERB_EXTREME);
      }

      if (debug) {
        std::ostringstream os;
        os << *debugPrefix << "Make B_crs MultiVector" << endl;
        std::cerr << os.str ();
      }
      RCP<tpetra_multivector_type> B_crs = Teuchos::rcp(new tpetra_multivector_type(A_bcrs->getRangeMap(),  nrhs));

      if (debug) {
        std::ostringstream os;
        os << *debugPrefix << "Benchmark CrsMatrix sparse matrix-vector multiply" << endl;
        std::cerr << os.str ();
      }
      // perform on point crs matrix
      {
        for (LO iter=0;iter<repeat;++iter) {
          TimeMonitor timerPointCrsApply(*TimeMonitor::getNewTimer("7) PointCrs Apply"));
          A_crs->apply(*X, *B_crs);
        }
      }

      if (debug) {
        std::ostringstream os;
        os << *debugPrefix << "Check results" << endl;
        std::cerr << os.str ();
      }
      // Check B_bcrs against B_crs
      {
        B_crs->update(-1.0, *B_bcrs, 1.0);

        Kokkos::View<typename tpetra_multivector_type::dot_type*, host_space> norm2 ("norm2", nrhs);
        Kokkos::View<typename tpetra_multivector_type::dot_type*, host_space> diff2 ("diff2", nrhs);
        B_bcrs->dot(*B_bcrs, norm2);
        B_crs->dot(*B_crs, diff2);

        if (myProcessPrintsDuringTheTest) {
          for (LO i = 0; i < nrhs; ++i) {
            std::cout << "Column = " << i << "  Error norm = "
                      << std::sqrt(diff2(i)/norm2(i)) << endl;
          }
        }
      }

      if (dump_sparse_matrix) {
        if (debug) {
          std::ostringstream os;
          os << *debugPrefix << "Write CrsMatrix to Matrix Market file" << endl;
          std::cerr << os.str ();
        }
        // no function to export block crs
        TimeMonitor timerMatrixMarket(*TimeMonitor::getNewTimer("8) Export MatrixMarket "));
        std::ofstream ofs("BlockCrsTestMatrix.out", std::ofstream::out);
        Tpetra::MatrixMarket::Writer<tpetra_crs_matrix_type>::writeSparse(ofs, A_crs);
        ofs.close();
      }
    } // end global timer

    // Print out timing results.
    TimeMonitor::report(commTest.ptr(), std::cout);

    // This tells the Trilinos test framework that the test passed.
    if (commTest->getRank() == 0) {
      std::cout << "End Result: TEST PASSED" << endl;
    }
  }
  Tpetra::finalize ();
  return EXIT_SUCCESS;
}
