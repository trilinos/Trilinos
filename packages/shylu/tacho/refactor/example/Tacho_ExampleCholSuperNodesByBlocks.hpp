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

namespace Tacho {

  template<typename DeviceSpaceType>
  int exampleCholSuperNodesByBlocks(const std::string file_input,
                                    const int treecut,
                                    const int prunecut,
                                    const int max_concurrency,
                                    const int max_task_dependence,
                                    const int team_size,
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

    // for simple test, let's use host space only here, for device it needs mirroring.

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,HostSpaceType> CrsMatrixBaseHostType;
    typedef CrsMatrixViewExt<CrsMatrixBaseHostType> CrsMatrixViewHostType;

    typedef GraphTools<ordinal_type,size_type,HostSpaceType> GraphToolsHostType;

    typedef GraphTools_Scotch<ordinal_type,size_type,HostSpaceType> GraphToolsHostType_Scotch;
    typedef GraphTools_CAMD<ordinal_type,size_type,HostSpaceType> GraphToolsHostType_CAMD;

    typedef Kokkos::Experimental::TaskPolicy<DeviceSpaceType> PolicyType;

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
    typedef Kokkos::Experimental::Future<int,HostSpaceType> future_type;

#ifdef HAVE_SHYLUTACHO_VTUNE
    __itt_pause();
#endif

    int r_val = 0;

    Kokkos::Impl::Timer timer;

    ///
    /// Read from matrix market
    ///
    ///     input  - file
    ///     output - AA_host
    ///
    CrsMatrixBaseHostType AA_host("AA_host");
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
        MatrixMarket::read(AA_host, in);
      } else if (extension == "crs") {
        std::cout << "CholSuperNodesByBlocks:: Input matrix is CRS data format" << std::endl;
        CrsData::read(AA_host, in);
        CrsMatrixTools::sortColumnsPerRow(AA_host);
      }
    }
    double t_read = timer.seconds();

    if (verbose)
      AA_host.showMe(std::cout) << std::endl;

    ///
    /// Create a graph structure for Scotch and CAMD (rptr, cidx)
    ///
    ///     rptr and cidx are need to be set up for Scotch and CAMD
    ///
    typename GraphToolsHostType::size_type_array rptr("Graph::RowPtrArray", AA_host.NumRows() + 1);
    typename GraphToolsHostType::ordinal_type_array cidx("Graph::ColIndexArray", AA_host.NumNonZeros());

    ///
    /// Run Scotch
    ///
    ///     input  - rptr, cidx, A_host
    ///     output - S (perm, iperm, nblks, range, tree), AA_scotch_host (permuted)
    ///
    timer.reset();
    GraphToolsHostType::getGraph(rptr, cidx, AA_host);
    double t_graph = timer.seconds();

    GraphToolsHostType_Scotch S;
    S.setGraph(AA_host.NumRows(), rptr, cidx);
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
    double t_scotch = timer.seconds();

    if (verbose)
      S.showMe(std::cout) << std::endl;

    CrsMatrixBaseHostType AA_scotch_host("AA_scotch_host");
    AA_scotch_host.createConfTo(AA_host);

    CrsMatrixTools::copy(AA_scotch_host,
                         S.PermVector(),
                         S.InvPermVector(),
                         AA_host);

    if (verbose)
      AA_scotch_host.showMe(std::cout) << std::endl;

    ///
    /// Run CAMD
    ///
    ///     input  - rptr, cidx, A_host
    ///     output - S (perm, iperm, nblks, range, tree), AA_scotch_host (permuted)
    ///
    timer.reset();
    GraphToolsHostType::getGraph(rptr, cidx, AA_scotch_host);
    t_graph += timer.seconds();

    GraphToolsHostType_CAMD C;
    C.setGraph(AA_scotch_host.NumRows(),
               rptr, cidx,
               S.NumBlocks(),
               S.RangeVector());

    timer.reset();
    C.computeOrdering();
    double t_camd = timer.seconds();

    if (verbose)
      C.showMe(std::cout) << std::endl;

    CrsMatrixBaseHostType AA_camd_host("AA_camd_host");
    AA_camd_host.createConfTo(AA_scotch_host);

    CrsMatrixTools::copy(AA_camd_host,
                         C.PermVector(),
                         C.InvPermVector(),
                         AA_scotch_host);

    if (verbose)
      AA_camd_host.showMe(std::cout) << std::endl;

    CrsMatrixBaseHostType AA_reordered_host = AA_camd_host;

    ///
    /// Clean tempoerary matrices
    ///
    ///     input  - AA_scotch_host, AA_reordered_host, C, rptr, cidx
    ///     output - none
    ///
    AA_scotch_host = CrsMatrixBaseHostType();
    AA_camd_host   = CrsMatrixBaseHostType();

    C = GraphToolsHostType_CAMD();
    rptr = typename GraphToolsHostType::size_type_array();
    cidx = typename GraphToolsHostType::ordinal_type_array();

    ///
    /// Symbolic factorization
    ///
    ///     input  -
    ///     output - S (perm, iperm, nblks, range, tree), AA_scotch_host (permuted)
    ///
    CrsHierBaseHostType HA_factor_host("HA_factor_host");
    Kokkos::View<typename CrsMatrixViewHostType::row_view_type*,HostSpaceType> rows;

    timer.reset();
    {
      typename CrsHierBaseHostType::size_type_array ap;
      typename CrsHierBaseHostType::ordinal_type_array aj;
      {
        CrsHierBaseHostType HA_reordered_host("HA_reordered_host");
        
        // create block matrices wrapping reordred matrix
        S.pruneTree(prunecut);
        
        const bool full = true;
        CrsMatrixTools::createHierMatrix(HA_reordered_host,
                                         AA_reordered_host,
                                         S.NumBlocks(),
                                         S.RangeVector(),
                                         S.TreeVector(),
                                         full);
        
        {
          const size_type nblocks = HA_reordered_host.NumNonZeros();
          Kokkos::View<size_type*,HostSpaceType> apblk("apblk", nblocks + 1);
          apblk(0) = 0;
          for (size_type k=0;k<nblocks;++k) {
            const auto &block = HA_reordered_host.Value(k);
            apblk(k+1) = apblk(k) + block.NumRows();
          }
          rows = Kokkos::View<typename CrsMatrixViewHostType::row_view_type*,HostSpaceType>("RowViewsInBlocks", apblk(nblocks));
          Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, nblocks),
                               [&](const ordinal_type k) {
                                 auto &block = HA_reordered_host.Value(k);
                                 block.setRowViewArray(Kokkos::subview(rows, range_type(apblk(k), apblk(k+1))));
                               } );
        }
        CrsMatrixTools::filterEmptyBlocks(HA_reordered_host);
        CrsMatrixTools::createSymbolicStructure(ap, aj, HA_reordered_host);
      }
      
      typename CrsHierBaseHostType::value_type_array ax("ax", aj.dimension(0));
      {
        auto flat = AA_reordered_host;
        auto range = S.RangeVector();
        
        const auto m = S.NumBlocks();
        for (auto i=0;i<m;++i) {
          const auto row_begin = ap(i);
          const auto row_end = ap(i+1);
          for (auto idx=row_begin;idx<row_end;++idx) {
            const auto j = aj(idx);
            ax(idx).setView(flat, range(i), (range(i+1) - range(i)),
                            /**/  range(j), (range(j+1) - range(j)));
          }
        }
      }
      {
        const auto m = S.NumBlocks();
        const auto nnz = aj.dimension(0);
        HA_factor_host = CrsHierBaseHostType("HA_factor_host",
                                             m, m, nnz,
                                             Kokkos::subview(ap, Kokkos::pair<ordinal_type,ordinal_type>(0,m)),
                                             Kokkos::subview(ap, Kokkos::pair<ordinal_type,ordinal_type>(1,m+1)),
                                             aj,
                                             ax);
        HA_factor_host.setNumNonZeros();
      }
      
      // copy row view structure to block symbolic factors
      {
        const size_type nblocks = HA_factor_host.NumNonZeros();
        Kokkos::View<size_type*,HostSpaceType> apblk("apblk", nblocks + 1);
        apblk(0) = 0;
        for (size_type k=0;k<nblocks;++k) {
          const auto &block = HA_factor_host.Value(k);
          apblk(k+1) = apblk(k) + block.NumRows();
        }
        rows = Kokkos::View<typename CrsMatrixViewHostType::row_view_type*,HostSpaceType>("RowViewsInBlocks", apblk(nblocks));
        Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, nblocks),
                             [&](const ordinal_type k) {
                               auto &block = HA_factor_host.Value(k);
                               block.setRowViewArray(Kokkos::subview(rows, range_type(apblk(k), apblk(k+1))));
                             } );
      }
    }
    double t_symbolic = timer.seconds();

    ///
    /// Create task policy
    ///
    ///     input  - max_task_size
    ///     output - policy
    ///
    const size_type max_task_size = (3*sizeof(CrsTaskViewHostType)+sizeof(PolicyType)+128);

    timer.reset();
    PolicyType policy(max_concurrency,
                      max_task_size,
                      max_task_dependence,
                      team_size);
    double t_policy = timer.seconds();

    ///
    /// Task parallel execution
    ///
    ///    input  - AA_factor_host, rowviews
    ///    output - HA_factor_host, AA_factor_host, B_factor_host
    ///
    double t_blocks = 0, t_chol = 0;

    Kokkos::View<value_type*,HostSpaceType> mats;
    Kokkos::View<DenseTaskViewHostType*,HostSpaceType> blks;
    {
      const size_type nblocks = HA_factor_host.NumNonZeros();
      timer.reset();
      { 
        Kokkos::View<size_type*,HostSpaceType> apblk("apblk", nblocks + 1);
        apblk(0) = 0;
        for (size_type k=0;k<nblocks;++k) {
          const auto &block = HA_factor_host.Value(k);
          apblk(k+1) = apblk(k) + block.NumRows()*block.NumCols();
        }

        mats = Kokkos::View<value_type*,HostSpaceType>("MatsInBlocks", apblk(nblocks));
        Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, nblocks),
                             [&](const ordinal_type k) {
                               auto &block = HA_factor_host.Value(k);
                               block.Flat().setExternalMatrix(block.NumRows(),
                                                              block.NumCols(),
                                                              -1, -1,
                                                              Kokkos::subview(mats, range_type(apblk(k), apblk(k+1))));
                               block.copyToFlat();
                             } );
      }
      if (mb) {
        Kokkos::View<size_type*,HostSpaceType> apblk("apblk", nblocks + 1);
        apblk(0) = 0;
        for (size_type k=0;k<nblocks;++k) {
          const auto &block = HA_factor_host.Value(k);
          ordinal_type hm, hn;
          DenseMatrixTools::getDimensionOfHierMatrix(hm, hn, block.Flat(), mb, mb);
          apblk(k+1) = apblk(k) + hm*hn;
        }

        blks = Kokkos::View<DenseTaskViewHostType*,HostSpaceType>("DenseBlocksInCrsBlocks", apblk(nblocks));

        Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, nblocks),
                             [&](const ordinal_type k) {
                               auto &block = HA_factor_host.Value(k);
                               ordinal_type hm, hn;
                               DenseMatrixTools::getDimensionOfHierMatrix(hm, hn, block.Flat(), mb, mb);
                               block.Hier().setExternalMatrix(hm, hn,
                                                              -1, -1,
                                                              Kokkos::subview(blks, range_type(apblk(k), apblk(k+1))));
                               DenseMatrixTools::getHierMatrix(block.Hier(),
                                                               block.Flat(),
                                                               mb, mb);
                             } );
      }
      t_blocks = timer.seconds();

      if (verbose || verbose_blocks)
        for (auto k=0;k<HA_factor_host.NumNonZeros();++k) 
          HA_factor_host.Value(k).showMe(std::cout) << std::endl;

      {
        const size_type nblocks = HA_factor_host.NumNonZeros();
        size_type nnz_blocks = 0, size_blocks = 0, max_blk_size = 0, max_blk_nnz = 0, max_blk_nrows = 0, max_blk_ncols = 0;
        for (size_type k=0;k<nblocks;++k) {
          const auto &block = HA_factor_host.Value(k);
          nnz_blocks  += block.NumNonZeros();
          const auto current_blk_size = block.NumRows()*block.NumCols();
          size_blocks += current_blk_size;

          if (max_blk_size < current_blk_size) {
            max_blk_nnz   = block.NumNonZeros();
            max_blk_nrows = block.NumRows();
            max_blk_ncols = block.NumCols();
            max_blk_size = current_blk_size;
          }
        }
        std::cout << "CholSuperNodesByBlocks:: nnz in blocks = " << nnz_blocks
                  << ", size of blocks = " << size_blocks
                  << ", ratio = " << (double(nnz_blocks)/size_blocks)
                  << std::endl;
        std::cout << "CholSuperNodesByBlocks:: max block info, nnz = " << max_blk_nnz
                  << ", dim = " << max_blk_nrows << " x " << max_blk_ncols
                  << ", ratio = " << (double(max_blk_nnz)/max_blk_size)
                  << std::endl;
      }

      CrsTaskHierViewHostType TA_factor_host(HA_factor_host);
#ifdef HAVE_SHYLUTACHO_VTUNE
      __itt_resume();
#endif
      timer.reset();
      {
        future_type future;
        if (mb) {
          // call nested block version
          std::cout << "CholSuperNodesByBlocks:: use DenseByBlocks with external LAPACK and BLAS" << std::endl;
          future = policy.proc_create_team(Chol<Uplo::Upper,AlgoChol::ByBlocks,Variant::Three>
                                           ::createTaskFunctor(policy,
                                                               TA_factor_host),
                                           0);
          policy.spawn(future);
        } else {
          // call plain block version
          std::cout << "CholSuperNodesByBlocks:: use external LAPACK and BLAS" << std::endl;
          future = policy.proc_create_team(Chol<Uplo::Upper,AlgoChol::ByBlocks,Variant::Two>
                                           ::createTaskFunctor(policy,
                                                               TA_factor_host),
                                           0);
          policy.spawn(future);
        }
        Kokkos::Experimental::wait(policy);
        TACHO_TEST_FOR_ABORT(future.get(), "Fail to perform CholeskySuperNodesByBlocks");
      }
      t_chol = timer.seconds();
#ifdef HAVE_SHYLUTACHO_VTUNE
      __itt_pause();
#endif
    }

    ///
    /// Solution check
    ///
    ///    input  - AA_reordered_host
    ///    output -
    ///
    double t_solve = 0;
    double error = 0, norm = 1;
    if (nrhs) {
      DenseMatrixBaseHostType 
        BB_host("BB_host", AA_reordered_host.NumRows(), nrhs), 
        XX_host("XX_host"), RR_host("RR_host");
      XX_host.createConfTo(BB_host);
      RR_host.createConfTo(BB_host);
      
      srand(time(NULL));
      const ordinal_type m = BB_host.NumRows();
      for (ordinal_type rhs=0;rhs<nrhs;++rhs) {
        for (ordinal_type i=0;i<m;++i)
          XX_host.Value(i, rhs) = ((value_type)rand()/(RAND_MAX));

        // matvec
        HostSpaceType::execution_space::fence();
        Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, m),
                             [&](const ordinal_type i) {
                               const auto nnz  = AA_reordered_host.NumNonZerosInRow(i);
                               const auto cols = AA_reordered_host.ColsInRow(i);
                               const auto vals = AA_reordered_host.ValuesInRow(i);

                               value_type tmp = 0;
                               for (ordinal_type j=0;j<nnz;++j)
                                 tmp += vals(j)*XX_host.Value(cols(j), rhs);
                               BB_host.Value(i, rhs) = tmp;
                             } );
        HostSpaceType::execution_space::fence();
      }
      if (verbose) {
        XX_host.showMe(std::cout) << std::endl;
        BB_host.showMe(std::cout) << std::endl;
      }
      DenseMatrixTools::copy(RR_host, XX_host); // keep solution on RR
      DenseMatrixTools::copy(XX_host, BB_host); // copy BB into XX

      DenseHierBaseHostType HX_host("HX_host");

      DenseMatrixTools::createHierMatrix(HX_host, XX_host,
                                         S.NumBlocks(),
                                         S.RangeVector(),
                                         nb);

      CrsTaskHierViewHostType TA_factor_host(HA_factor_host);
      DenseTaskHierViewHostType TX_host(HX_host);

      timer.reset();
      {
        auto future_forward_solve
          = policy.proc_create_team(TriSolve<Uplo::Upper,Trans::ConjTranspose,
                                    AlgoTriSolve::ByBlocks,Variant::Two>
                                    ::createTaskFunctor(policy,
                                                        Diag::NonUnit, TA_factor_host, TX_host),
                                    0);
        policy.spawn(future_forward_solve);

        auto future_backward_solve
          = policy.proc_create_team(TriSolve<Uplo::Upper,Trans::NoTranspose,
                                    AlgoTriSolve::ByBlocks,Variant::Two>
                                    ::createTaskFunctor(policy,
                                                        Diag::NonUnit, TA_factor_host, TX_host),
                                    1);

        policy.add_dependence(future_backward_solve, future_forward_solve);
        policy.spawn(future_backward_solve);

        Kokkos::Experimental::wait(policy);

        TACHO_TEST_FOR_ABORT(future_forward_solve.get(),  "Fail to perform TriSolveSuperNodesByBlocks (forward)");
        TACHO_TEST_FOR_ABORT(future_backward_solve.get(), "Fail to perform TriSolveSuperNodesByBlocks (backward)");
      }
      t_solve = timer.seconds();

      if (verbose) {
        XX_host.showMe(std::cout) << std::endl;
        BB_host.showMe(std::cout) << std::endl;
      }
      
      for (ordinal_type rhs=0;rhs<nrhs;++rhs) {
        for (ordinal_type i=0;i<m;++i) {
          error += Util::abs(XX_host.Value(i, rhs) - RR_host.Value(i, rhs));
          norm  += Util::abs(RR_host.Value(i, rhs));
        }
      }

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
      std::cout << "CholSuperNodesByBlocks:: Given    matrix = " << AA_host.NumRows() << " x " << AA_host.NumCols() << ", nnz = " << AA_host.NumNonZeros() << std::endl;
      std::cout << "CholSuperNodesByBlocks:: Hier     matrix = " << HA_factor_host.NumRows() << " x " << HA_factor_host.NumCols() << ", nnz = " << HA_factor_host.NumNonZeros() << std::endl;

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
        //<< "hier creation = " << t_hier << " [sec] "
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
