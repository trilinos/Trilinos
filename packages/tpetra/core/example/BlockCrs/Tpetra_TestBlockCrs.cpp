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

#include <MatrixMarket_Tpetra.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_FancyOStream.hpp>

#include "Tpetra_TestBlockCrsMeshDatabase.hpp"

//#define PRINT_VERBOSE 1

using namespace BlockCrsTest;

int main (int argc, char *argv[]) 
{
  using Teuchos::RCP;
  using Teuchos::TimeMonitor;

  const auto INVALID_GLOBAL_ORDINAL = Teuchos::OrdinalTraits<global_ordinal_type>::invalid();
  auto out = Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cout));
  
  // MPI boilerplate
  Tpetra::initialize(&argc, &argv);
  RCP<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform ().getComm();

  // Command-line input
  LO num_global_elements_i = 2, num_global_elements_j = 2, num_global_elements_k = 2;
  LO num_procs_i = 1, num_procs_j = 1, num_procs_k = 1;
  LO blocksize = 5, nrhs = 1, repeat = 100;

  Teuchos::CommandLineProcessor clp(false);
  clp.setDocString("BlockCrs performance test using 3D 7 point stencil.\n");
  clp.setOption("num-elements-i", &num_global_elements_i, "Number of cells in the I dimension.");
  clp.setOption("num-elements-j", &num_global_elements_j, "Number of cells in the J dimension.");
  clp.setOption("num-elements-k", &num_global_elements_k, "Number of cells in the K dimension.");
  clp.setOption("num-procs-i", &num_procs_i, "Processor grid of (npi,npj,npk); npi*npj*npk should be equal to the number of MPI ranks.");
  clp.setOption("num-procs-j", &num_procs_j, "Processor grid of (npi,npj,npk); npi*npj*npk should be equal to the number of MPI ranks.");
  clp.setOption("num-procs-k", &num_procs_k, "Processor grid of (npi,npj,npk); npi*npj*npk should be equal to the number of MPI ranks.");
  clp.setOption("blocksize", &blocksize, "Block size. The # of DOFs coupled in a multiphysics flow problem.");
  clp.setOption("nrhs", &nrhs, "Number of right hand sides to solve for.");
  clp.setOption("repeat", &repeat, "Number of iterations of matvec operations to measure performance.");

  {
    auto r_parse = clp.parse(argc, argv, comm->getRank() == 0 ? &std::cout : NULL);
    if (r_parse == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) return 0;
    if (r_parse != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL  ) return -1;

    if (num_procs_i*num_procs_j*num_procs_k != comm->getSize()) {
      Tpetra::finalize();
      if (comm->getRank() == 0) {
        std::cout << "Invalid processor grid ["<<num_procs_i<<"x"<<num_procs_j<<"x"<<num_procs_k<<"]"
                  << " does not match to the mpi rank size ("<<comm->getSize()<<")" << std::endl;
        std::cout << "End Result: TEST PASSED" << std::endl;
      }
      return -1;
    }
  }
  
  if (comm->getRank() == 0) {
    if (std::is_same<exec_space,Kokkos::Serial>::value)
      std::cout << "Kokkos::Serial " << std::endl;
    else
      exec_space::print_configuration(std::cout, false);
  }

  MeshDatabase mesh(comm,
                    num_global_elements_i,
                    num_global_elements_j,
                    num_global_elements_k,
                    num_procs_i,
                    num_procs_j,
                    num_procs_k);

  const auto sb = mesh.getStructuredBlock();
  const auto part = mesh.getStructuredBlockPart();

  const auto num_owned_elements = mesh.getNumOwnedElements();
  const auto num_remote_elements = mesh.getNumRemoteElements();
  const auto num_owned_and_remote_elements = mesh.getNumElements();

#if PRINT_VERBOSE
  mesh.print(std::cout);
#endif
  
  {
    TimeMonitor timerGlobal(*TimeMonitor::getNewTimer("X) Global"));
    
    // Build Tpetra Maps
    // -----------------
    
    const Tpetra::global_size_t TpetraComputeGlobalNumElements 
      = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();  
    
    // - all internal views are allocated on device; mirror as mesh database is constructed on host  
    const auto mesh_gids_host = mesh.getElementGlobalIDs();
    const auto mesh_gids = Kokkos::create_mirror_view(typename exec_space::memory_space(), mesh.getElementGlobalIDs());
    Kokkos::deep_copy(mesh_gids, mesh_gids_host);

    // for convenience, separate the access to owned and remote gids
    const auto owned_gids  = Kokkos::subview(mesh_gids, local_ordinal_range_type(                 0, num_owned_elements));
    const auto remote_gids = Kokkos::subview(mesh_gids, local_ordinal_range_type(num_owned_elements, num_owned_and_remote_elements));
    
    RCP<const map_type> row_map = rcp(new map_type(TpetraComputeGlobalNumElements, owned_gids, 0, comm));
    RCP<const map_type> col_map = rcp(new map_type(TpetraComputeGlobalNumElements, mesh_gids, 0, comm));
    
#if PRINT_VERBOSE
    row_map->describe(*out);
    col_map->describe(*out);
#endif
    
    // Graph Construction
    // ------------------
    // local graph is constructed on device space
    typedef tpetra_crs_graph_type::local_graph_type local_graph_type;
    typedef local_graph_type::row_map_type::non_const_type rowptr_view_type;
    typedef typename local_graph_type::entries_type colidx_view_type;

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
      const auto nnz_host = Kokkos::create_mirror_view(nnz);
      Kokkos::deep_copy(nnz_host, nnz);
      
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
    
    // Call fillComplete on the bcrs_graph to finalize it
    RCP<tpetra_crs_graph_type> bcrs_graph;
    { 
      TimeMonitor timerGlobalGraphConstruction(*TimeMonitor::getNewTimer("1) GlobalGraphConstruction"));
      bcrs_graph = rcp(new tpetra_crs_graph_type(row_map, col_map, local_graph_type(colidx, rowptr),
                                                Teuchos::null));
    } // end global graph timer
    
#if PRINT_VERBOSE
    bcrs_graph->describe(*out, Teuchos::VERB_EXTREME);
#endif

    // Create BlockCrsMatrix
    RCP<tpetra_blockcrs_matrix_type> A_bcrs = Teuchos::rcp(new tpetra_blockcrs_matrix_type(*bcrs_graph, blocksize));

    typedef tpetra_blockcrs_matrix_type::little_block_type block_type;
    Kokkos::View<block_type*,exec_space> blocks;
    {
      TimeMonitor timerLocalBlockCrsFill(*TimeMonitor::getNewTimer("2) LocalBlockCrsFill"));
      
      // Tpetra BlockCrsMatrix only has high level access functions
      // To fill this on device, we need an access to the meta data of blocks
      const auto rowptr_host = Kokkos::create_mirror_view(rowptr);      
      const auto colidx_host = Kokkos::create_mirror_view(colidx);
      
      Kokkos::deep_copy(rowptr_host, rowptr);
      Kokkos::deep_copy(colidx_host, colidx);

      blocks = Kokkos::View<block_type*,exec_space>("blocks", rowptr_host(num_owned_elements));

      const auto blocks_host = Kokkos::create_mirror_view(blocks);
      Kokkos::parallel_for
        (Kokkos::RangePolicy<host_space>(0,num_owned_elements), KOKKOS_LAMBDA(const LO row) {
          for (LO loc=rowptr_host(row);loc<rowptr_host(row+1);++loc) 
            blocks_host(loc) = A_bcrs->getLocalBlock(row, colidx(loc));
        });    
      Kokkos::deep_copy(blocks, blocks_host);

      Kokkos::parallel_for
        (num_owned_elements, KOKKOS_LAMBDA(const LO row) {
          for (LO loc=rowptr(row);loc<rowptr(row+1);++loc) {
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
#if PRINT_VERBOSE
    A_bcrs->describe(*out, Teuchos::VERB_EXTREME);
#endif

    // Create MultiVector
    RCP<tpetra_multivector_type> X = Teuchos::rcp(new tpetra_multivector_type(A_bcrs->getDomainMap(), nrhs));
    {
      TimeMonitor timerMultiVectorFill(*TimeMonitor::getNewTimer("4) MultiVectorFill"));
      
      auto value = X->getLocalView<typename exec_space::memory_space>();
      auto map = X->getMap()->getLocalMap();
      Kokkos::parallel_for
        (value.extent(0), KOKKOS_LAMBDA(const LO i) {
          const GO gid = map.getGlobalElement(i);
          for (LO j=0,jend=value.extent(1);j<jend;++j) 
            value(i, j) = get_multi_vector_entry<value_type>(gid, j);
        });
    }
    
#if PRINT_VERBOSE
    X->describe(*out, Teuchos::VERB_EXTREME);
#endif

    RCP<tpetra_multivector_type> B_bcrs = Teuchos::rcp(new tpetra_multivector_type(A_bcrs->getRangeMap(),  nrhs));

    // matrix vector multiplication
    {
      for (LO iter=0;iter<repeat;++iter) {
        TimeMonitor timerBlockCrsApply(*TimeMonitor::getNewTimer("5) BlockCrs Apply"));        
        A_bcrs->apply(*X, *B_bcrs);
      }
    }
      
    // direct conversion: block crs -> point crs 
    RCP<tpetra_crs_matrix_type> A_crs;
    {
      TimeMonitor timerConvertBlockCrsToPointCrs(*TimeMonitor::getNewTimer("6) Conversion from BlockCrs to PointCrs"));

      // construct row map and column map for a point crs matrix
      // point-wise row map can be obtained from A_bcrs->getDomainMap().
      // A constructor exist for crs matrix with a local matrix and a row map.
      // see, Tpetra_CrsMatrix_decl.hpp, line 504
      //     CrsMatrix (const local_matrix_type& lclMatrix,
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
      
      RCP<const map_type> row_crs_map = rcp(new map_type(TpetraComputeGlobalNumElements, owned_crs_gids, 0, comm));
      RCP<const map_type> col_crs_map = rcp(new map_type(TpetraComputeGlobalNumElements, crs_gids, 0, comm));

      rowptr_view_type crs_rowptr = rowptr_view_type("crs_rowptr", num_owned_elements*blocksize+1);      
      colidx_view_type crs_colidx = colidx_view_type("crs_colidx", colidx.extent(0)*blocksize*blocksize);
      typename tpetra_crs_matrix_type::local_matrix_type::values_type 
        crs_values("crs_values", colidx.extent(0)*blocksize*blocksize);        

      Kokkos::parallel_for
        (num_owned_elements,
         KOKKOS_LAMBDA(const LO &idx) {
          const GO nnz_per_block_row = rowptr(idx+1)-rowptr(idx);
          const GO nnz_per_point_row = nnz_per_block_row*blocksize;
          const GO crs_rowptr_begin  = idx*blocksize;
          const GO crs_colidx_begin  = rowptr(idx)*blocksize*blocksize;

          for (LO i=0;i<(blocksize+1);++i) 
            crs_rowptr(crs_rowptr_begin+i) = crs_colidx_begin + i*nnz_per_point_row;
          
          GO loc = crs_colidx_begin;
          // loop over the rows in a block
          for (LO l0=0;l0<blocksize;++l0) {
            // loop over the block row
            for (GO jj=rowptr(idx);jj<rowptr(idx+1);++jj) {
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

      typename tpetra_crs_matrix_type::local_matrix_type 
        local_matrix("local_crs_matrix",
                     num_owned_and_remote_elements*blocksize,
                     crs_values,
                     local_graph_type(crs_colidx, crs_rowptr));
      
      A_crs = rcp(new tpetra_crs_matrix_type(row_crs_map,
                                             col_crs_map,
                                             local_matrix, 
                                             Teuchos::null));
      
    } // end conversion timer

#if PRINT_VERBOSE
    A_crs->describe(*out, Teuchos::VERB_EXTREME);
#endif

    RCP<tpetra_multivector_type> B_crs = Teuchos::rcp(new tpetra_multivector_type(A_bcrs->getRangeMap(),  nrhs));

    // perform on point crs matrix
    {
      for (LO iter=0;iter<repeat;++iter) {
        TimeMonitor timerPointCrsApply(*TimeMonitor::getNewTimer("7) PointCrs Apply"));
        A_crs->apply(*X, *B_crs);
      }
    }

    // veryfy B_bcrs and B_crs are identical
    {
      B_crs->update(-1.0, *B_bcrs, 1.0);

      Kokkos::View<typename tpetra_multivector_type::dot_type*,host_space> 
        norm2("norm2", nrhs), diff2("diff2", nrhs);;
      B_bcrs->dot(*B_bcrs, norm2);
      B_crs->dot(*B_crs, diff2);
      if (comm->getRank() == 0) 
        for (LO i=0;i<nrhs;++i) 
          std::cout << "Column = " << i << "  Error norm = " << std::sqrt(diff2(i)/norm2(i)) << "\n";
    }
    // Save crs_matrix as a MatrixMarket file.
    {
      // no function to export block crs 
      TimeMonitor timerMatrixMarket(*TimeMonitor::getNewTimer("8) Export MatrixMarket "));
      std::ofstream ofs("BlockCrsTestMatrix.out", std::ofstream::out);
      Tpetra::MatrixMarket::Writer<tpetra_crs_matrix_type>::writeSparse(ofs, A_crs);
      ofs.close();
    }      

  } // end global timer

  // Print out timing results.
  TimeMonitor::report(comm.ptr(), std::cout, "");

  Tpetra::finalize();
  
  // This tells the Trilinos test framework that the test passed.
  if (comm->getRank() == 0) {
    std::cout << "End Result: TEST PASSED" << std::endl;
  }
  
  return EXIT_SUCCESS;
}  // END main()
