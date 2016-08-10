#ifndef __TACHO_EXAMPLE_CHOL_SUPERNODES_BY_BLOCKS_HPP__
#define __TACHO_EXAMPLE_CHOL_SUPERNODES_BY_BLOCKS_HPP__

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

//#define TACHO_EXECUTE_TASKS_SERIAL

#include "Tacho_Util.hpp"

#include "Tacho_CrsMatrixBase.hpp"
#include "Tacho_CrsMatrixView.hpp"
#include "Tacho_CrsRowView.hpp"

#include "Tacho_DenseMatrixBase.hpp"
#include "Tacho_DenseMatrixView.hpp"
#include "Tacho_CrsMatrixViewExt.hpp"

#include "Tacho_CrsMatrixTools.hpp"
#include "Tacho_DenseMatrixTools.hpp"

#include "Tacho_MatrixMarket.hpp"
#include "Tacho_CrsData.hpp"

#include "Tacho_GraphTools.hpp"

#include "Tacho_GraphTools_Scotch.hpp"
#include "Tacho_GraphTools_CAMD.hpp"

#include "Tacho_TaskView.hpp"
#include "Tacho_TaskFactory.hpp"

#include "Tacho_Gemm.hpp"
#include "Tacho_Herk.hpp"
#include "Tacho_Trsm.hpp"
#include "Tacho_Chol.hpp"
#include "Tacho_TriSolve.hpp"

#ifdef HAVE_SHYLUTACHO_VTUNE
#include "ittnotify.h"
#endif

#ifdef HAVE_SHYLUTACHO_MKL
#include "mkl_service.h"
#endif

namespace Tacho {

  template<typename DeviceSpaceType>
  int exampleCholSuperNodesByBlocks(const std::string file_input,
                                    const int treecut,
                                    const int prunecut,
                                    const int max_concurrency,
                                    const int memory_pool_grain_size,
                                    const int mkl_nthreads,
                                    const int nrhs,
                                    const int mb,
                                    const int nb,
                                    const bool verbose_blocks,
                                    const bool verbose) {
    typedef typename
      Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;

    const bool detail = false;
    std::cout << "DeviceSpace::  "; DeviceSpaceType::print_configuration(std::cout, detail);
    std::cout << "HostSpace::    ";   HostSpaceType::print_configuration(std::cout, detail);
    std::cout << std::endl;

    // for simple test, let's use host space only here, for device it needs mirroring.

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,HostSpaceType> CrsMatrixBaseHostType;
    typedef CrsMatrixViewExt<CrsMatrixBaseHostType> CrsMatrixViewHostType;

    typedef GraphTools<ordinal_type,size_type,HostSpaceType> GraphToolsHostType;

    typedef GraphTools_Scotch<ordinal_type,size_type,HostSpaceType> GraphToolsHostType_Scotch;
    typedef GraphTools_CAMD<ordinal_type,size_type,HostSpaceType> GraphToolsHostType_CAMD;

    typedef Kokkos::TaskPolicy<DeviceSpaceType> PolicyType;

    typedef TaskView<CrsMatrixViewHostType> CrsTaskViewHostType;
    typedef CrsMatrixBase<CrsTaskViewHostType,ordinal_type,size_type,HostSpaceType> CrsHierBaseHostType;
    typedef CrsMatrixView<CrsHierBaseHostType> CrsHierViewHostType;
    typedef TaskView<CrsHierViewHostType> CrsTaskHierViewHostType;

    typedef DenseMatrixBase<value_type,ordinal_type,size_type,HostSpaceType> DenseMatrixBaseHostType;
    typedef DenseMatrixView<DenseMatrixBaseHostType> DenseMatrixViewHostType;

    typedef TaskView<DenseMatrixViewHostType> DenseTaskViewHostType;

    typedef DenseMatrixBase<DenseTaskViewHostType,ordinal_type,size_type,HostSpaceType> DenseHierBaseHostType;
    typedef DenseMatrixView<DenseHierBaseHostType> DenseHierViewHostType;
    typedef TaskView<DenseHierViewHostType> DenseTaskHierViewHostType;

    typedef Kokkos::pair<size_type,size_type> range_type;
    typedef Kokkos::Future<int,HostSpaceType> future_type;

#ifdef HAVE_SHYLUTACHO_VTUNE
    __itt_pause();
#endif

    int r_val = 0;

    Kokkos::Impl::Timer timer;

    /// Phase 0 : Read input matrix 
    /// ------------------------------------------------------------------------------------

    ///
    /// Read from matrix market
    ///
    ///     input  - file
    ///     output - AA
    ///
    CrsMatrixBaseHostType AA("AA");
    timer.reset();
    {
      std::ifstream in;
      in.open(file_input);
      if (!in.good()) {
        std::cout << "Failed in open the file: " << file_input << std::endl;
        return -1;
      }

      const auto extension = file_input.substr(file_input.find_last_of(".") + 1);
      if        (extension == "mtx") {
        std::cout << "CholSuperNodesByBlocks:: Input matrix is MatrixMarket format" << std::endl;
        MatrixMarket::read(AA, in);
      } else if (extension == "crs") {
        std::cout << "CholSuperNodesByBlocks:: Input matrix is CRS data format" << std::endl;
        CrsData::read(AA, in);
        CrsMatrixTools::sortColumnsPerRow(AA);
      }
    }
    const double t_read = timer.seconds();

    if (verbose)
      AA.showMe(std::cout) << std::endl;


    /// Phase 1 : Reorder matrix
    /// ------------------------------------------------------------------------------------

    ///
    /// Create a graph structure for Scotch and CAMD (rptr, cidx)
    ///
    ///     rptr and cidx are need to be set up for Scotch and CAMD
    ///
    typename GraphToolsHostType::size_type_array rptr("Graph::RowPtrArray", AA.NumRows() + 1);
    typename GraphToolsHostType::ordinal_type_array cidx("Graph::ColIndexArray", AA.NumNonZeros());

    ///
    /// Run Scotch
    ///
    ///     input  - rptr, cidx, AA
    ///     output - S (perm, iperm, nblks, range, tree), AA_scotch (permuted)
    ///
    timer.reset();
    GraphToolsHostType::getGraph(rptr, cidx, AA);
    double t_graph = timer.seconds();

    GraphToolsHostType_Scotch S;
    S.setGraph(AA.NumRows(), rptr, cidx);
    S.setSeed(0);
    S.setTreeLevel();
    S.setStrategy( SCOTCH_STRATSPEED
                   | SCOTCH_STRATLEVELMAX
                   | SCOTCH_STRATLEVELMIN
                   | SCOTCH_STRATLEAFSIMPLE
                   | SCOTCH_STRATSEPASIMPLE
                   );

    timer.reset();
    S.computeOrdering(treecut);
    S.pruneTree(prunecut);
    const double t_scotch = timer.seconds();

    if (verbose)
      S.showMe(std::cout) << std::endl;

    CrsMatrixBaseHostType AA_scotch("AA_scotch");
    AA_scotch.createConfTo(AA);
    
    CrsMatrixTools::copy(AA_scotch,
                         S.PermVector(),
                         S.InvPermVector(),
                         AA);

    if (verbose)
      AA_scotch.showMe(std::cout) << std::endl;

    ///
    /// Run CAMD
    ///
    ///     input  - rptr, cidx, AA, S
    ///     output - C (perm, iperm), AA_camd (permuted)
    ///
    timer.reset();
    GraphToolsHostType::getGraph(rptr, cidx, AA_scotch);
    t_graph += timer.seconds();

    GraphToolsHostType_CAMD C;
    C.setGraph(AA_scotch.NumRows(),
               rptr, cidx,
               S.NumBlocks(),
               S.RangeVector());

    timer.reset();
    C.computeOrdering();
    const double t_camd = timer.seconds();

    if (verbose)
      C.showMe(std::cout) << std::endl;

    CrsMatrixBaseHostType AA_camd("AA_camd");
    AA_camd.createConfTo(AA_scotch);

    CrsMatrixTools::copy(AA_camd,
                         C.PermVector(),
                         C.InvPermVector(),
                         AA_scotch);

    if (verbose)
      AA_camd.showMe(std::cout) << std::endl;

    ///
    /// Assign reordered matrix
    ///
    ///     input  - either AA_scotch or AA_camd
    ///     output - AA_reordered
    ///
    CrsMatrixBaseHostType AA_reordered = AA_camd;

    ///
    /// Clean tempoerary matrices
    ///
    ///     input  - AA_scotch, AA_reordered, C, rptr, cidx
    ///     output - none
    ///
    AA_scotch = CrsMatrixBaseHostType();
    AA_camd   = CrsMatrixBaseHostType();

    C = GraphToolsHostType_CAMD();
    rptr = typename GraphToolsHostType::size_type_array();
    cidx = typename GraphToolsHostType::ordinal_type_array();

    /// Phase 2 : Construct task policy (task policy is always constructed before any future
    ///           is declared.
    /// ------------------------------------------------------------------------------------

    ///
    /// Create task policy
    ///
    ///     input  - max_task_size
    ///     output - policy
    ///
    const size_type max_task_size = (3*sizeof(CrsTaskViewHostType)+sizeof(PolicyType)+128);

    timer.reset();

    PolicyType policy( typename PolicyType::memory_space(),
                       max_task_size*max_concurrency,
                       memory_pool_grain_size );

    const double t_policy = timer.seconds();

    /// Phase 3 : Perform symbolic factorization 
    /// ------------------------------------------------------------------------------------

    ///
    /// Symbolic factorization
    ///
    ///     input  - AA_reordered
    ///     output - HA_factor
    ///
    CrsHierBaseHostType HA_factor("HA_factor");
    Kokkos::View<typename CrsMatrixViewHostType::row_view_type*,HostSpaceType> rows;

    timer.reset();
    {
      CrsHierBaseHostType HA_reordered("HA_reordered");
      
      // find fills in block matrix
      typename CrsHierBaseHostType::size_type_array ap;
      typename CrsHierBaseHostType::ordinal_type_array aj;
      {
        // to apply symbolic factorization on the matrix of blocks, we need a graph of 
        // the block matrix (full matrix).
        const bool full = true;
        CrsMatrixTools::createHierMatrix(HA_reordered,
                                         AA_reordered,
                                         S.NumBlocks(),
                                         S.RangeVector(),
                                         S.TreeVector(),
                                         full);
        
        {
          const size_type nblocks = HA_reordered.NumNonZeros();
          Kokkos::View<size_type*,HostSpaceType> offs("offs", nblocks + 1);
          offs(0) = 0;
          for (size_type k=0;k<nblocks;++k) {
            const auto &block = HA_reordered.Value(k);
            offs(k+1) = offs(k) + block.NumRows();
          }
          rows = Kokkos::View<typename CrsMatrixViewHostType::row_view_type*,HostSpaceType>("RowViewsInBlocks", offs(nblocks));
          Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, nblocks),
                               [&](const size_type k) {
                                 auto &block = HA_reordered.Value(k);
                                 block.setRowViewArray(Kokkos::subview(rows, range_type(offs(k), offs(k+1))));
                               } );
        }
        CrsMatrixTools::filterEmptyBlocks(HA_reordered);
        CrsMatrixTools::createSymbolicStructure(ap, aj, HA_reordered);
      }

      // fill block information 
      typename CrsHierBaseHostType::value_type_array ax("ax", aj.dimension(0));
      {
        const auto range = S.RangeVector();        
        const ordinal_type m = S.NumBlocks();
        Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, m),
                             [&](const ordinal_type i) {
                               const auto beg = ap(i);
                               const auto end = ap(i+1);
                               for (auto idx=beg;idx<end;++idx) {
                                 const auto j = aj(idx);
                                 ax(idx).setView(AA_reordered, range(i), (range(i+1) - range(i)),
                                                 /**/          range(j), (range(j+1) - range(j)));
                               }
                             } );
      }

      // construct hierachical matrix
      {
        const ordinal_type m = S.NumBlocks();
        const auto nnz = aj.dimension(0);
        const auto ap_begin = Kokkos::subview(ap, Kokkos::pair<ordinal_type,ordinal_type>(0,m));
        const auto ap_end   = Kokkos::subview(ap, Kokkos::pair<ordinal_type,ordinal_type>(1,m+1));
        HA_factor = CrsHierBaseHostType("HA_factor", m, m, nnz, ap_begin, ap_end, aj, ax);
        HA_factor.setNumNonZeros();
      }
      
      // copy row view structure to block symbolic factors
      {
        const size_type m = HA_reordered.NumRows();
        Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, m),
                             [&](const ordinal_type i) {
                               const auto cols_a = HA_reordered.ColsInRow(i);
                               const auto vals_a = HA_reordered.ValuesInRow(i);

                               const auto cols_f = HA_factor.ColsInRow(i);
                               const auto vals_f = HA_factor.ValuesInRow(i);
                               
                               const size_type nnz_in_a = HA_reordered.NumNonZerosInRow(i);
                               const size_type nnz_in_f = HA_factor.NumNonZerosInRow(i);

                               for (size_type idx_a=0,idx_f=0;idx_a<nnz_in_a && idx_f<nnz_in_f;) {
                                 const auto j_a = cols_a(idx_a);
                                 const auto j_f = cols_f(idx_f);
                                 
                                 if (j_a == j_f)
                                   vals_f(idx_f) = vals_a(idx_a);

                                 idx_a += (j_a <= j_f);
                                 idx_f += (j_a >= j_f);
                               }
                             } );
      }
    }
    const double t_symbolic = timer.seconds();


    /// Phase 4 : Perform numeric factorization 
    /// ------------------------------------------------------------------------------------

#ifdef HAVE_SHYLUTACHO_MKL
    mkl_set_num_threads(mkl_nthreads);
#endif

    ///
    /// Allocate blocks 
    ///
    ///    input  - HA_factor
    ///    output - mats, blks
    ///
    Kokkos::View<value_type*,HostSpaceType> mats;
    Kokkos::View<DenseTaskViewHostType*,HostSpaceType> blks;
    
    timer.reset();
    {
      const size_type nblocks = HA_factor.NumNonZeros();
      { 
        Kokkos::View<size_type*,HostSpaceType> offs("offs", nblocks + 1);
        offs(0) = 0;
        for (size_type k=0;k<nblocks;++k) {
          const auto &block = HA_factor.Value(k);
          offs(k+1) = offs(k) + block.NumRows()*block.NumCols();
        }
        
        mats = Kokkos::View<value_type*,HostSpaceType>("MatsInBlocks", offs(nblocks));
        Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, nblocks),
                             [&](const ordinal_type k) {
                               auto &block = HA_factor.Value(k);
                               block.Flat().setExternalMatrix(block.NumRows(),
                                                              block.NumCols(),
                                                              -1, -1,
                                                              Kokkos::subview(mats, range_type(offs(k), offs(k+1))));
                               block.copyToFlat();
                             } );
      }
      if (mb) {
        Kokkos::View<size_type*,HostSpaceType> offs("offs", nblocks + 1);
        offs(0) = 0;
        for (size_type k=0;k<nblocks;++k) {
          const auto &block = HA_factor.Value(k);
          ordinal_type hm, hn;
          DenseMatrixTools::getDimensionOfHierMatrix(hm, hn, block.Flat(), mb, mb);
          offs(k+1) = offs(k) + hm*hn;
        }
        
        blks = Kokkos::View<DenseTaskViewHostType*,HostSpaceType>("DenseBlocksInCrsBlocks", offs(nblocks));
        Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, nblocks),
                             [&](const ordinal_type k) {
                               auto &block = HA_factor.Value(k);
                               ordinal_type hm, hn;
                               DenseMatrixTools::getDimensionOfHierMatrix(hm, hn, block.Flat(), mb, mb);
                               block.Hier().setExternalMatrix(hm, hn,
                                                              -1, -1,
                                                              Kokkos::subview(blks, range_type(offs(k), offs(k+1))));
                               DenseMatrixTools::getHierMatrix(block.Hier(),
                                                               block.Flat(),
                                                               mb, mb);
                             } );
      }

    }
    const double t_blocks = timer.seconds();

    if (verbose || verbose_blocks)
      for (auto k=0;k<HA_factor.NumNonZeros();++k) 
        HA_factor.Value(k).showMe(std::cout) << std::endl;
      
    {
      const size_type nblocks = HA_factor.NumNonZeros();
      size_type nnz_blocks = 0, size_blocks = 0, max_blk_size = 0, max_blk_nrows = 0, max_blk_ncols = 0;
      for (size_type k=0;k<nblocks;++k) {
        const auto &block = HA_factor.Value(k);
        nnz_blocks  += block.NumNonZeros();
        
        const auto current_blk_size = block.NumRows()*block.NumCols();
        size_blocks += current_blk_size;
        
        if (max_blk_size < current_blk_size) {
          max_blk_nrows = block.NumRows();
          max_blk_ncols = block.NumCols();
          max_blk_size  = current_blk_size;
        }
      }
      std::cout << "CholSuperNodesByBlocks:: "
                << "size of blocks = " << size_blocks
                << ", max block = " << max_blk_nrows << " x " << max_blk_ncols
                << std::endl;
    }

    ///
    /// Perform numeric factorization
    ///
    ///    input  - HA_factor
    ///    output - HA_factor
    ///
    timer.reset();    
    {
      CrsTaskHierViewHostType TA_factor(HA_factor);
#ifdef HAVE_SHYLUTACHO_VTUNE
      __itt_resume();
#endif
      {
        future_type future;
        if (mb) {
          // call nested block version
          std::cout << "CholSuperNodesByBlocks:: use DenseByBlocks with external LAPACK and BLAS" << std::endl;
          future = policy.host_spawn(Chol<Uplo::Upper,AlgoChol::ByBlocks,Variant::Three>
                                     ::createTaskFunctor(policy,
                                                         TA_factor));
        } else {
          // call plain block version
          std::cout << "CholSuperNodesByBlocks:: use external LAPACK and BLAS" << std::endl;
          future = policy.host_spawn(Chol<Uplo::Upper,AlgoChol::ByBlocks,Variant::Two>
                                     ::createTaskFunctor(policy,
                                                         TA_factor));
        }
        TACHO_TEST_FOR_EXCEPTION(future.is_null(), std::runtime_error,
                                 ">> host_spawn returns a null future");
        Kokkos::wait(policy);
        TACHO_TEST_FOR_ABORT(future.get(), "Fail to perform CholeskySuperNodesByBlocks");
      }
#ifdef HAVE_SHYLUTACHO_VTUNE
      __itt_pause();
#endif
    }
    const double t_chol = timer.seconds();    

    /// Phase 4 : Solve problem
    /// ------------------------------------------------------------------------------------
    
    ///
    /// Solution check
    ///
    ///    input  - AA_reordered, BB, XX
    ///    output - RR
    ///
    double t_solve = 0;
    double error = 0, norm = 1;
    if (nrhs) {
      const auto m = AA_reordered.NumRows();
      DenseMatrixBaseHostType BB("BB", m, nrhs), XX("XX"), RR("RR");
      XX.createConfTo(BB);
      RR.createConfTo(BB);
      
      srand(time(NULL));
      for (ordinal_type rhs=0;rhs<nrhs;++rhs) {
        for (ordinal_type i=0;i<m;++i)
          XX.Value(i, rhs) = ((value_type)rand()/(RAND_MAX));
        
        // matvec
        HostSpaceType::execution_space::fence();
        Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, m),
                             [&](const ordinal_type i) {
                               const auto nnz  = AA_reordered.NumNonZerosInRow(i);
                               const auto cols = AA_reordered.ColsInRow(i);
                               const auto vals = AA_reordered.ValuesInRow(i);
                               
                               value_type tmp = 0;
                               for (ordinal_type j=0;j<nnz;++j)
                                 tmp += vals(j)*XX.Value(cols(j), rhs);
                               BB.Value(i, rhs) = tmp;
                             } );
        HostSpaceType::execution_space::fence();
      }
      if (verbose) {
        XX.showMe(std::cout) << std::endl;
        BB.showMe(std::cout) << std::endl;
      }
      DenseMatrixTools::copy(RR, XX); // keep solution on RR
      DenseMatrixTools::copy(XX, BB); // copy BB into XX

      DenseHierBaseHostType HX("HX");

      DenseMatrixTools::createHierMatrix(HX, XX,
                                         S.NumBlocks(),
                                         S.RangeVector(),
                                         nb);

      CrsTaskHierViewHostType TA_factor(HA_factor);
      DenseTaskHierViewHostType TX(HX);

      timer.reset();
      {
        auto future_forward_solve
          = policy.host_spawn(TriSolve<Uplo::Upper,Trans::ConjTranspose,
                              AlgoTriSolve::ByBlocks,Variant::Two>
                              ::createTaskFunctor(policy,
                                                  Diag::NonUnit, TA_factor, TX));
        TACHO_TEST_FOR_EXCEPTION(future_forward_solve.is_null(), std::runtime_error,
                                 ">> host_spawn returns a null future");


        auto future_backward_solve
          = policy.host_spawn(TriSolve<Uplo::Upper,Trans::NoTranspose,
                              AlgoTriSolve::ByBlocks,Variant::Two>
                              ::createTaskFunctor(policy,
                                                  Diag::NonUnit, TA_factor, TX),
                              future_forward_solve);
        TACHO_TEST_FOR_EXCEPTION(future_backward_solve.is_null(), std::runtime_error,
                                 ">> host_spawn returns a null future");

        Kokkos::wait(policy);

        TACHO_TEST_FOR_EXCEPTION(future_forward_solve.get(),  std::logic_error, "Fail to perform TriSolveSuperNodesByBlocks (forward)");
        TACHO_TEST_FOR_EXCEPTION(future_backward_solve.get(), std::logic_error, "Fail to perform TriSolveSuperNodesByBlocks (backward)");
      }
      t_solve = timer.seconds();

      if (verbose) {
        XX.showMe(std::cout) << std::endl;
        BB.showMe(std::cout) << std::endl;
      }
      
      for (ordinal_type rhs=0;rhs<nrhs;++rhs) {
        for (ordinal_type i=0;i<m;++i) {
          {
            const auto val = Util::abs(XX.Value(i, rhs) - RR.Value(i, rhs));
            error += val*val;
          }
          {
            const auto val = Util::abs(RR.Value(i, rhs));
            norm  += val*val;
          }
        }
      }
      error = std::sqrt(error);
      norm  = std::sqrt(norm);

      std::cout << std::scientific;
      std::cout << "CholSuperNodesByBlocks:: error = " << error
                << ", norm = " << norm
                << ", rel error = " << (error/norm)
                << std::endl;
      std::cout.unsetf(std::ios::scientific);
    }

    ///
    /// Print out
    ///
    {
      const auto prec = std::cout.precision();
      std::cout.precision(4);

      std::cout << std::scientific;
      std::cout << "CholSuperNodesByBlocks:: Given    matrix = " << AA.NumRows() << " x " << AA.NumCols() << ", nnz = " << AA.NumNonZeros() << std::endl;
      std::cout << "CholSuperNodesByBlocks:: Hier     matrix = " << HA_factor.NumRows() << " x " << HA_factor.NumCols() << ", nnz = " << HA_factor.NumNonZeros() << std::endl;

      std::cout << "CholSuperNodesByBlocks:: "
                << "read = " << t_read << " [sec], "
                << "graph generation = " << (t_graph/2.0) << " [sec] "
                << "scotch reordering = " << t_scotch << " [sec] "
                << "camd reordering = " << t_camd << " [sec] "
                << std::endl
                << "CholSuperNodesByBlocks:: "
                << "symbolic factorization = " << t_symbolic << " [sec] "
                << std::endl
                << "CholSuperNodesByBlocks:: "
                << "policy creation = " << t_policy << " [sec] "
                << "block specification = " << t_blocks << " [sec] "
                << std::endl
                << "CholSuperNodesByBlocks:: "
                << "Chol = " << t_chol << " [sec] ";
      if (nrhs)
        std::cout << "Solve = " << t_solve << " [sec] ";

      std::cout << std::endl;
      
      std::cout.unsetf(std::ios::scientific);
      std::cout.precision(prec);
    }

    return r_val;
  }
}

#endif
