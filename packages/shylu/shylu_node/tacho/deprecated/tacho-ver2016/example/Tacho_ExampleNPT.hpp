#ifndef __TACHO_EXAMPLE_NPT_HPP__
#define __TACHO_EXAMPLE_NPT_HPP__

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "TachoExp_Util.hpp"
//#include "TachoExp_CrsMatrixBase.hpp"
//#include "TachoExp_MatrixMarket.hpp"
//#include "Tacho_CrsMatrixView.hpp"
//#include "Tacho_CrsRowView.hpp"

#ifdef HAVE_SHYLU_NODETACHO_VTUNE
#include "ittnotify.h"
#endif

#ifdef HAVE_SHYLU_NODETACHO_MKL
#include "mkl_service.h"
#endif

namespace Tacho {
  namespace Experimental {
    template<typename DeviceSpaceType>
    int exampleNPT(const std::string file_input,
                   const int max_concurrency,
                   const int memory_pool_grain_size,
                   const int mkl_nthreads,
                   const int nrhs,
                   const int mb,
                   const int nb,
                   const bool verbose_blocks,
                   const bool verbose) {
      typedef Kokkos::DefaultHostExecutionSpace HostSpaceType;

      const bool detail = false;
      std::cout << "DeviceSpace::  "; DeviceSpaceType::print_configuration(std::cout, detail);
      std::cout << "HostSpace::    ";   HostSpaceType::print_configuration(std::cout, detail);
      std::cout << std::endl;

      //typedef CrsMatrixBase<ValueType,HostSpaceType> CrsMatrixBaseHostType;

#ifdef HAVE_SHYLU_NODETACHO_VTUNE
      __itt_pause();
#endif

      int r_val = 0;

      Kokkos::Impl::Timer timer;

      std::cout << " partiiton id = " << Partition::Top::tag 
                << "  " << (Partition::Top() == Partition::Top())
                << std::endl;
      
      // /// Phase 0 : Read input matrix 
      // /// ------------------------------------------------------------------------------------

      // ///
      // /// Read from matrix market
      // ///
      // ///     input  - file
      // ///     output - AA
      // ///
      // CrsMatrixBaseHostType AA("AA");
      // timer.reset();
      // {
      //   const auto extension = file_input.substr(file_input.find_last_of(".") + 1);
      //   if        (extension == "mtx") {
      //     std::cout << "NPT:: Input matrix is MatrixMarket format" << std::endl;
      //     AA = MatrixMarket<ValueType>::read(file_input);
      //   }
      // }
      // const double t_read = timer.seconds();
      
      // if (verbose) {
      //   std::ofstream out("tacho_npt_input.mtx");
      //   MatrixMarket<ValueType>::write(out, AA);
      //   // AA.showMe(std::cout) << std::endl;
      // }
      
      // /// Phase 1 : Reorder matrix
      // /// ------------------------------------------------------------------------------------

      // ///
      // /// Create a graph structure for Scotch and CAMD (rptr, cidx)
      // ///
      // ///     rptr and cidx are need to be set up for Scotch and CAMD
      // ///
      // typename GraphToolsHostType::size_type_array rptr("Graph::RowPtrArray", AA.NumRows() + 1);
      // typename GraphToolsHostType::ordinal_type_array cidx("Graph::ColIndexArray", AA.NumNonZeros());
      
      // ///
      // /// Run Scotch
      // ///
      // ///     input  - rptr, cidx, AA
      // ///     output - S (perm, iperm, nblks, range, tree), AA_scotch (permuted)
      // ///
      // timer.reset();
      // GraphToolsHostType::getGraph(rptr, cidx, AA);
      // double t_graph = timer.seconds();

      // GraphToolsHostType_Scotch S;
      // S.setGraph(AA.NumRows(), rptr, cidx);
      // S.setSeed(0);
      // S.setTreeLevel();
      // S.setStrategy( SCOTCH_STRATSPEED
      //                | SCOTCH_STRATLEVELMAX
      //                | SCOTCH_STRATLEVELMIN
      //                | SCOTCH_STRATLEAFSIMPLE
      //                | SCOTCH_STRATSEPASIMPLE
      //                );

      // timer.reset();
      // S.computeOrdering(treecut);
      // S.pruneTree(prunecut);
      // const double t_scotch = timer.seconds();

      // if (verbose) {
      //   std::ofstream out("tacho_scotch.txt");
      //   S.showMe(out) << std::endl;
      // }

      // CrsMatrixBaseHostType AA_scotch("AA_scotch");
      // AA_scotch.createConfTo(AA);
    
      // CrsMatrixTools::copy(AA_scotch,
      //                      S.PermVector(),
      //                      S.InvPermVector(),
      //                      AA);

      // if (verbose) {
      //   std::ofstream out("tacho_s.mtx");
      //   MatrixMarket::write(out, AA_scotch);
      //   //AA_scotch.showMe(std::cout) << std::endl;
      // }

      return r_val;
    }
  }
}

#endif














//       ///
//       /// Assign reordered matrix
//       ///
//       ///     input  - either AA_scotch or AA_camd
//       ///     output - AA_reordered
//       ///
//       CrsMatrixBaseHostType AA_reordered = AA_camd;
//       {
//         std::ofstream out("tacho_r.mtx");
//         MatrixMarket::write(out, AA_reordered);
//       }

//       ///
//       /// Clean tempoerary matrices
//       ///
//       ///     input  - AA_scotch, AA_reordered, C, rptr, cidx
//       ///     output - none
//       ///
//       AA_scotch = CrsMatrixBaseHostType();
//       AA_camd   = CrsMatrixBaseHostType();

//       C = GraphToolsHostType_CAMD();
//       rptr = typename GraphToolsHostType::size_type_array();
//       cidx = typename GraphToolsHostType::ordinal_type_array();

//       /// Phase 2 : Construct task policy (task policy is always constructed before any future
//       ///           is declared.
//       /// ------------------------------------------------------------------------------------

//       ///
//       /// Create task policy
//       ///
//       ///     input  - max_task_size
//       ///     output - policy
//       ///
//       const SizeType max_task_size = (3*sizeof(CrsTaskViewHostType)+sizeof(PolicyType)+128);

//       timer.reset();

//       PolicyType policy( typename PolicyType::memory_space(),
//                          max_task_size*max_concurrency,
//                          memory_pool_grain_size );

//       const double t_policy = timer.seconds();

//       /// Phase 3 : Perform symbolic factorization 
//       /// ------------------------------------------------------------------------------------

//       ///
//       /// Symbolic factorization
//       ///
//       ///     input  - AA_reordered
//       ///     output - HA_factor
//       ///
//       CrsHierBaseHostType HA_factor("HA_factor");
//       Kokkos::View<typename CrsMatrixViewHostType::row_view_type*,HostSpaceType> rows;

//       timer.reset();
//       {
//         CrsHierBaseHostType HA_reordered("HA_reordered");
      
//         // find fills in block matrix
//         typename CrsHierBaseHostType::size_type_array ap;
//         typename CrsHierBaseHostType::ordinal_type_array aj;
//         {
//           // to apply symbolic factorization on the matrix of blocks, we need a graph of 
//           // the block matrix (full matrix).
//           const bool full = true;
//           CrsMatrixTools::createHierMatrix(HA_reordered,
//                                            AA_reordered,
//                                            S.NumBlocks(),
//                                            S.RangeVector(),
//                                            S.TreeVector(),
//                                            full);
        
//           {
//             const SizeType nblocks = HA_reordered.NumNonZeros();
//             Kokkos::View<SizeType*,HostSpaceType> offs("offs", nblocks + 1);
//             offs(0) = 0;
//             for (SizeType k=0;k<nblocks;++k) {
//               const auto &block = HA_reordered.Value(k);
//               offs(k+1) = offs(k) + block.NumRows();
//             }
//             rows = Kokkos::View<typename CrsMatrixViewHostType::row_view_type*,HostSpaceType>("RowViewsInBlocks", offs(nblocks));
//             Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, nblocks),
//                                  [&](const SizeType k) {
//                                    auto &block = HA_reordered.Value(k);
//                                    block.setRowViewArray(Kokkos::subview(rows, range_type(offs(k), offs(k+1))));
//                                  } );
//           }
//           CrsMatrixTools::filterEmptyBlocks(HA_reordered);
//           CrsMatrixTools::createSymbolicStructure(ap, aj, HA_reordered);
//         }

//         // fill block information 
//         typename CrsHierBaseHostType::value_type_array ax("ax", aj.dimension(0));
//         {
//           const auto range = S.RangeVector();        
//           const OrdinalType m = S.NumBlocks();
//           Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, m),
//                                [&](const OrdinalType i) {
//                                  const auto beg = ap(i);
//                                  const auto end = ap(i+1);
//                                  for (auto idx=beg;idx<end;++idx) {
//                                    const auto j = aj(idx);
//                                    ax(idx).setView(AA_reordered, range(i), (range(i+1) - range(i)),
//                                                    /**/          range(j), (range(j+1) - range(j)));
//                                  }
//                                } );
//         }

//         // construct hierachical matrix
//         {
//           const OrdinalType m = S.NumBlocks();
//           const auto nnz = aj.dimension(0);
//           const auto ap_begin = Kokkos::subview(ap, Kokkos::pair<OrdinalType,OrdinalType>(0,m));
//           const auto ap_end   = Kokkos::subview(ap, Kokkos::pair<OrdinalType,OrdinalType>(1,m+1));
//           HA_factor = CrsHierBaseHostType("HA_factor", m, m, nnz, ap_begin, ap_end, aj, ax);
//           HA_factor.setNumNonZeros();
//         }
      
//         // copy row view structure to block symbolic factors
//         {
//           const SizeType m = HA_reordered.NumRows();
//           Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, m),
//                                [&](const OrdinalType i) {
//                                  const auto cols_a = HA_reordered.ColsInRow(i);
//                                  const auto vals_a = HA_reordered.ValuesInRow(i);

//                                  const auto cols_f = HA_factor.ColsInRow(i);
//                                  const auto vals_f = HA_factor.ValuesInRow(i);
                               
//                                  const SizeType nnz_in_a = HA_reordered.NumNonZerosInRow(i);
//                                  const SizeType nnz_in_f = HA_factor.NumNonZerosInRow(i);

//                                  for (SizeType idx_a=0,idx_f=0;idx_a<nnz_in_a && idx_f<nnz_in_f;) {
//                                    const auto j_a = cols_a(idx_a);
//                                    const auto j_f = cols_f(idx_f);
                                 
//                                    if (j_a == j_f)
//                                      vals_f(idx_f) = vals_a(idx_a);

//                                    idx_a += (j_a <= j_f);
//                                    idx_f += (j_a >= j_f);
//                                  }
//                                } );
//         }
//         // if one wants to see block structure for given matrix
//         // if (verbose || verbose_blocks) {
//         //   std::ofstream out("tacho_r.blk");              
//         //   for (auto k=0;k<HA_reordered.NumNonZeros();++k) 
//         //     HA_reordered.Value(k).showMe(out) << std::endl;
//         // }
//       }
//       const double t_symbolic = timer.seconds();


//       /// Phase 4 : Perform numeric factorization 
//       /// ------------------------------------------------------------------------------------

// #ifdef HAVE_SHYLU_NODETACHO_MKL
//       mkl_set_num_threads(mkl_nthreads);
// #endif

//       ///
//       /// Allocate blocks 
//       ///
//       ///    input  - HA_factor
//       ///    output - mats, blks
//       ///
//       Kokkos::View<ValueType*,HostSpaceType> mats;
//       Kokkos::View<DenseTaskViewHostType*,HostSpaceType> blks;
    
//       timer.reset();
//       {
//         const SizeType nblocks = HA_factor.NumNonZeros();
//         { 
//           Kokkos::View<SizeType*,HostSpaceType> offs("offs", nblocks + 1);
//           offs(0) = 0;
//           for (SizeType k=0;k<nblocks;++k) {
//             const auto &block = HA_factor.Value(k);
//             offs(k+1) = offs(k) + block.NumRows()*block.NumCols();
//           }
        
//           mats = Kokkos::View<ValueType*,HostSpaceType>("MatsInBlocks", offs(nblocks));
//           Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, nblocks),
//                                [&](const OrdinalType k) {
//                                  auto &block = HA_factor.Value(k);
//                                  block.Flat().setExternalMatrix(block.NumRows(),
//                                                                 block.NumCols(),
//                                                                 -1, -1,
//                                                                 Kokkos::subview(mats, range_type(offs(k), offs(k+1))));
//                                  block.copyToFlat();
//                                } );
//         }
//         if (mb) {
//           Kokkos::View<SizeType*,HostSpaceType> offs("offs", nblocks + 1);
//           offs(0) = 0;
//           for (SizeType k=0;k<nblocks;++k) {
//             const auto &block = HA_factor.Value(k);
//             OrdinalType hm, hn;
//             DenseMatrixTools::getDimensionOfHierMatrix(hm, hn, block.Flat(), mb, mb);
//             offs(k+1) = offs(k) + hm*hn;
//           }
        
//           blks = Kokkos::View<DenseTaskViewHostType*,HostSpaceType>("DenseBlocksInCrsBlocks", offs(nblocks));
//           Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, nblocks),
//                                [&](const OrdinalType k) {
//                                  auto &block = HA_factor.Value(k);
//                                  OrdinalType hm, hn;
//                                  DenseMatrixTools::getDimensionOfHierMatrix(hm, hn, block.Flat(), mb, mb);
//                                  block.Hier().setExternalMatrix(hm, hn,
//                                                                 -1, -1,
//                                                                 Kokkos::subview(blks, range_type(offs(k), offs(k+1))));
//                                  Impl::DenseMatrixTools::Serial::getHierMatrix(block.Hier(),
//                                                                                block.Flat(),
//                                                                                mb, mb);
//                                } );
//         }

//       }
//       const double t_blocks = timer.seconds();

//       if (verbose || verbose_blocks) {
//         std::ofstream out("tacho_f.blk");              
//         for (auto k=0;k<HA_factor.NumNonZeros();++k) 
//           HA_factor.Value(k).showMe(out) << std::endl;
//       }

//       {
//         const SizeType nblocks = HA_factor.NumNonZeros();
//         SizeType nnz_blocks = 0, size_blocks = 0, max_blk_size = 0, max_blk_nrows = 0, max_blk_ncols = 0;
//         for (SizeType k=0;k<nblocks;++k) {
//           const auto &block = HA_factor.Value(k);
//           nnz_blocks  += block.NumNonZeros();
        
//           const auto current_blk_size = block.NumRows()*block.NumCols();
//           size_blocks += current_blk_size;
        
//           if (max_blk_size < current_blk_size) {
//             max_blk_nrows = block.NumRows();
//             max_blk_ncols = block.NumCols();
//             max_blk_size  = current_blk_size;
//           }
//         }
//         std::cout << "CholSuperNodesByBlocks:: "
//                   << "# of blocks = " << nblocks
//                   << ", size of blocks = " << size_blocks
//                   << ", max block = " << max_blk_nrows << " x " << max_blk_ncols
//                   << std::endl;
//       }

//       ///
//       /// Perform numeric factorization
//       ///
//       ///    input  - HA_factor
//       ///    output - HA_factor
//       ///

//       Stat statChol;
//       {
//         CrsTaskHierViewHostType TA_factor(HA_factor);
//         if (mb) {
//           std::cout << "CholSuperNodesByBlocks:: stat using DenseByBlocks with external LAPACK and BLAS" << std::endl;
//           statChol = Chol<Uplo::Upper,AlgoChol::ByBlocks,Variant::Three>::stat(TA_factor);
//         } else {
//           std::cout << "CholSuperNodesByBlocks:: stat using external LAPACK and BLAS" << std::endl;
//           statChol = Chol<Uplo::Upper,AlgoChol::ByBlocks,Variant::Two>::stat(TA_factor);
//         }
//       }

//       timer.reset();    
//       {
//         CrsTaskHierViewHostType TA_factor(HA_factor);
// #ifdef HAVE_SHYLU_NODETACHO_VTUNE
//         __itt_resume();
// #endif
//         {
//           future_type future;
//           if (mb) {
//             // call nested block version
//             std::cout << "CholSuperNodesByBlocks:: use DenseByBlocks with external LAPACK and BLAS" << std::endl;
//             future = policy.host_spawn(Chol<Uplo::Upper,AlgoChol::ByBlocks,Variant::Three>
//                                        ::createTaskFunctor(policy,
//                                                            TA_factor));
//           } else {
//             // call plain block version
//             std::cout << "CholSuperNodesByBlocks:: use external LAPACK and BLAS" << std::endl;
//             future = policy.host_spawn(Chol<Uplo::Upper,AlgoChol::ByBlocks,Variant::Two>
//                                        ::createTaskFunctor(policy,
//                                                            TA_factor));
//           }
//           TACHO_TEST_FOR_EXCEPTION(future.is_null(), std::runtime_error,
//                                    ">> host_spawn returns a null future");
//           Kokkos::wait(policy);
//           TACHO_TEST_FOR_ABORT(future.get(), "Fail to perform CholeskySuperNodesByBlocks");
//         }
// #ifdef HAVE_SHYLU_NODETACHO_VTUNE
//         __itt_pause();
// #endif
//       }
//       const double t_chol = timer.seconds();    

//       if (verbose) {
//         std::ofstream out("tacho_f.mtx");
//         out.precision(16);
//         out << "%%MatrixMarket matrix coordinate "
//             << "real "
//             << "general "
//             << std::endl;
//         SizeType nnz = 0, cnt = 0;
//         const SizeType nnz_blk = HA_factor.NumNonZeros();
//         const ValueType eps = 1.0e-11;
//         for (OrdinalType k=0;k<nnz_blk;++k) {
//           const auto blk = HA_factor.Value(k);
//           const auto flat = blk.Flat();
//           const OrdinalType 
//             offm = blk.OffsetRows()+1,
//             offn = blk.OffsetCols()+1,
//             m = blk.NumRows(),
//             n = blk.NumCols();
        
//           if (offm == offn) {
//             // print upper part only
//             for (OrdinalType j=0;j<m;++j)
//               for (OrdinalType i=0;i<(j+1);++i,++cnt)
//                 if (std::abs(flat.Value(i,j)) > eps) {
//                   out << (offm+i) << " " << (offn+j) << " " << flat.Value(i,j) << std::endl;
//                   ++nnz;
//                 }
//           } else {
//             // print entire matrix
//             for (OrdinalType j=0;j<n;++j)
//               for (OrdinalType i=0;i<m;++i,++cnt)
//                 if (std::abs(flat.Value(i,j)) > eps) {
//                   out << (offm+i) << " " << (offn+j) << " " << flat.Value(i,j) << std::endl;          
//                   ++nnz;
//                 }
//           }
//         }
//         out << "## " << AA_reordered.NumRows() << " " << AA_reordered.NumCols() << " " << nnz << " ## cnt " << cnt << std::endl; 
//       }

//       /// Phase 4 : Solve problem
//       /// ------------------------------------------------------------------------------------
    
//       ///
//       /// Solution check
//       ///
//       ///    input  - AA_reordered, BB, XX
//       ///    output - RR
//       ///
//       double t_solve = 0;
//       double error = 0, norm = 0;
//       if (nrhs) {
//         const auto m = AA_reordered.NumRows();
//         DenseMatrixBaseHostType BB("BB", m, nrhs), XX("XX"), RR("RR");
//         XX.createConfTo(BB);
//         RR.createConfTo(BB);
      
//         srand(time(NULL));
//         for (OrdinalType rhs=0;rhs<nrhs;++rhs) {
//           for (OrdinalType i=0;i<m;++i)
//             XX.Value(i, rhs) = ((ValueType)rand()/(RAND_MAX));
        
//           // matvec
//           HostSpaceType::execution_space::fence();
//           Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, m),
//                                [&](const OrdinalType i) {
//                                  const auto nnz  = AA_reordered.NumNonZerosInRow(i);
//                                  const auto cols = AA_reordered.ColsInRow(i);
//                                  const auto vals = AA_reordered.ValuesInRow(i);
                               
//                                  ValueType tmp = 0;
//                                  for (OrdinalType j=0;j<nnz;++j)
//                                    tmp += vals(j)*XX.Value(cols(j), rhs);
//                                  BB.Value(i, rhs) = tmp;
//                                } );
//           HostSpaceType::execution_space::fence();
//         }
//         if (verbose) {
//           std::ofstream out("tacho_x_and_b_before.txt");
//           XX.showMe(out) << std::endl;
//           BB.showMe(out) << std::endl;
//         }
//         DenseMatrixTools::copy(RR, XX); // keep solution on RR
//         DenseMatrixTools::copy(XX, BB); // copy BB into XX

//         DenseHierBaseHostType HX("HX");

//         DenseMatrixTools::createHierMatrix(HX, XX,
//                                            S.NumBlocks(),
//                                            S.RangeVector(),
//                                            nb);

//         CrsTaskHierViewHostType TA_factor(HA_factor);
//         DenseTaskHierViewHostType TX(HX);

//         timer.reset();
//         {
//           auto future_forward_solve
//             = policy.host_spawn(TriSolve<Uplo::Upper,Trans::ConjTranspose,
//                                 AlgoTriSolve::ByBlocks,Variant::Two>
//                                 ::createTaskFunctor(policy,
//                                                     Diag::NonUnit, TA_factor, TX));
//           TACHO_TEST_FOR_EXCEPTION(future_forward_solve.is_null(), std::runtime_error,
//                                    ">> host_spawn returns a null future");


//           auto future_backward_solve
//             = policy.host_spawn(TriSolve<Uplo::Upper,Trans::NoTranspose,
//                                 AlgoTriSolve::ByBlocks,Variant::Two>
//                                 ::createTaskFunctor(policy,
//                                                     Diag::NonUnit, TA_factor, TX),
//                                 future_forward_solve);
//           TACHO_TEST_FOR_EXCEPTION(future_backward_solve.is_null(), std::runtime_error,
//                                    ">> host_spawn returns a null future");

//           Kokkos::wait(policy);

//           TACHO_TEST_FOR_EXCEPTION(future_forward_solve.get(),  std::logic_error, "Fail to perform TriSolveSuperNodesByBlocks (forward)");
//           TACHO_TEST_FOR_EXCEPTION(future_backward_solve.get(), std::logic_error, "Fail to perform TriSolveSuperNodesByBlocks (backward)");
//         }
//         t_solve = timer.seconds();

//         if (verbose) {
//           std::ofstream out("tacho_x_and_b_after.txt");
//           XX.showMe(out) << std::endl;
//           BB.showMe(out) << std::endl;
//         }
      
//         for (OrdinalType rhs=0;rhs<nrhs;++rhs) {
//           for (OrdinalType i=0;i<m;++i) {
//             {
//               const auto val = Util::abs(XX.Value(i, rhs) - RR.Value(i, rhs));
//               error += val*val;
//             }
//             {
//               const auto val = Util::abs(RR.Value(i, rhs));
//               norm  += val*val;
//             }
//           }
//         }
//         error = std::sqrt(error);
//         norm  = std::sqrt(norm);

//         std::cout << std::scientific;
//         std::cout << "CholSuperNodesByBlocks:: error = " << error
//                   << ", norm = " << norm
//                   << ", rel error = " << (error/norm)
//                   << std::endl;
//         std::cout.unsetf(std::ios::scientific);
//       }

//       ///
//       /// Print out
//       ///
//       {
//         const auto prec = std::cout.precision();
//         std::cout.precision(4);

//         std::cout << std::scientific;
//         std::cout << "CholSuperNodesByBlocks:: Given    matrix = " << AA.NumRows() << " x " << AA.NumCols() << ", nnz = " << AA.NumNonZeros() << std::endl;
//         std::cout << "CholSuperNodesByBlocks:: Hier     matrix = " << HA_factor.NumRows() << " x " << HA_factor.NumCols() << ", nnz = " << HA_factor.NumNonZeros() << std::endl;

//         std::cout << "CholSuperNodesByBlocks:: "
//                   << "read = " << t_read << " [sec], "
//                   << "graph generation = " << (t_graph/2.0) << " [sec] "
//                   << "scotch reordering = " << t_scotch << " [sec] "
//                   << "camd reordering = " << t_camd << " [sec] "
//                   << std::endl
//                   << "CholSuperNodesByBlocks:: "
//                   << "symbolic factorization = " << t_symbolic << " [sec] "
//                   << std::endl
//                   << "CholSuperNodesByBlocks:: "
//                   << "policy creation = " << t_policy << " [sec] "
//                   << "block specification = " << t_blocks << " [sec] "
//                   << std::endl
//                   << "CholSuperNodesByBlocks:: "
//                   << "Chol = " 
//                   << statChol.flop/1.0e9 << " [GFLOP] " 
//                   << t_chol << " [sec] " 
//                   << (statChol.flop/t_chol)/1.0e9 << " [GFLOPs] ";

//         if (nrhs)
//           std::cout << "Solve = " << t_solve << " [sec] ";

//         std::cout << std::endl;
      
//         std::cout.unsetf(std::ios::scientific);
//         std::cout.precision(prec);
//       }
