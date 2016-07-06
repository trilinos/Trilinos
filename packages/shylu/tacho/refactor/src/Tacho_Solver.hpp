#ifndef __TACHO_SOLVER_HPP__
#define __TACHO_SOLVER_HPP__

/// \file Tacho_Solver.hpp
/// \brief Front solver interface
/// \author Kyungjoo Kim (kyukim@sandia.gov)


#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

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

#include "Tacho_SymbolicFactorization.hpp"

#include "Tacho_TaskView.hpp"
#include "Tacho_TaskFactory.hpp"

#include "Tacho_Gemm.hpp"
#include "Tacho_Herk.hpp"
#include "Tacho_Trsm.hpp"
#include "Tacho_Chol.hpp"
#include "Tacho_TriSolve.hpp"

namespace Tacho {

  template<typename ValueType,
           typename OrdinalType,
           typename SizeType,
           typename DeviceSpaceType>
  class Solver {
  public:
    // basic typedefs
    typedef ValueType   value_type;
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;

    // Host space objects
    typedef typename Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType;

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,HostSpaceType> CrsMatrixBaseHostType;
    typedef CrsMatrixViewExt<CrsMatrixBaseHostType> CrsMatrixViewHostType;

    typedef GraphTools<ordinal_type,ordinal_type,HostSpaceType> GraphToolsHostType;

    typedef GraphTools_Scotch<ordinal_type,ordinal_type,HostSpaceType> GraphToolsHostType_Scotch;
    typedef GraphTools_CAMD<ordinal_type,ordinal_type,HostSpaceType> GraphToolsHostType_CAMD;

    typedef SymbolicFactorization<CrsMatrixBaseHostType> SymbolicFactorizationType;

    // Policy on device
    typedef Kokkos::Experimental::TaskPolicy<DeviceSpaceType> PolicyType;

    // Hierarchical block matrices
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

    struct PhaseMask { 
      static constexpr int Reorder = 1;
      static constexpr int Analyze = 2;
      static constexpr int Factorize = 4;
      static constexpr int Solve = 8;
      static constexpr int MaxPhase = 16;
    };

    struct ProblemSetMask { 
      static constexpr int Matrix = 1;
      static constexpr int LeftHandSide = 2;
      static constexpr int RightHandSide = 4;
      static constexpr int WorkspaceAllocated = 8;
      static constexpr int SetProblem = (Matrix | LeftHandSide | RightHandSide);
    };

  private:
    std::string _label;

    // Kokkos policy
    size_type _max_concurrency;
    PolicyType _policy;
    
    // // algorithm
    // int _algo;

    // input matrices
    CrsMatrixBaseHostType _AA;
    DenseMatrixBaseHostType _BB, _XX, _YY;

    // extract permutation and blocking information
    CrsMatrixBaseHostType _AA_reorder;    
    ordinal_type _nblks;
    Kokkos::View<ordinal_type*,HostSpaceType> _perm, _peri, _range, _tree;
    
    // symbolic factorization pattern
    CrsMatrixBaseHostType _FF;

    // block hierachical matrix
    CrsHierBaseHostType   _HF;
    ordinal_type _mb;
    Kokkos::View<typename CrsMatrixViewHostType::row_view_type*,HostSpaceType> _rows;
    Kokkos::View<value_type*,HostSpaceType> _mats;
    Kokkos::View<DenseTaskViewHostType*,HostSpaceType> _blks;

    // block hierarchical matrix for right hand side
    DenseHierBaseHostType _HY;
    ordinal_type _nb;

    // problem set: 3 bits are used to indicate problem set
    ordinal_type _prob_set;

    // verbose output
    double _t[PhaseMask::MaxPhase];
    //bool _verbose;

  public:
    Solver(const std::string label = "Tacho::Solver")
      : _label(label)
    {}

    Solver(const Solver &b) = default;

    CrsMatrixBaseHostType   getA() const { return _AA; }
    DenseMatrixBaseHostType getB() const { return _BB; }
    DenseMatrixBaseHostType getX() const { return _XX; }

    CrsMatrixBaseHostType getFlatU() const { return _FF; }
    CrsHierBaseHostType   getHierU() const { return _HF; }

    void setPolicy(const size_type max_concurrency = 1000000) {
      const size_type max_task_size = (3*sizeof(CrsTaskViewHostType)+sizeof(PolicyType)+128);
      const size_type max_task_dependence = 3;
      const size_type team_size = 1;

      _policy = PolicyType(max_concurrency,
                           max_task_size,
                           max_task_dependence,
                           team_size);
    }

    void setProblem(const CrsMatrixBaseHostType AA,
                    const DenseMatrixBaseHostType BB,
                    const DenseMatrixBaseHostType XX) {
      _AA = AA;
      _XX = XX;
      _BB = BB;
      _prob_set = ProblemSetMask::SetProblem;
    }

    void setMatrix(const CrsMatrixBaseHostType AA) {
      _AA = AA;
      _prob_set |= ProblemSetMask::Matrix;
    }
    void setLeftHandSide(const DenseMatrixBaseHostType XX) {
      _XX = XX;
      _prob_set |=  ProblemSetMask::LeftHandSide;
      _prob_set &= ~ProblemSetMask::WorkspaceAllocated;
    }
    void setRightHandSide(const DenseMatrixBaseHostType BB) {
      _BB = BB;
      _prob_set |= ProblemSetMask::RightHandSide;
    }

    // void setAlgorithm(const int algo) {
    //   switch (algo) {
    //   case AlgoChol::SuperNodes:
    //   case AlgoChol::SuperNodesByBlocks: 
    //     break;
    //   default: {
    //     TACHO_TEST_FOR_ABORT(true, "Fail to perform CholeskySuperNodesByBlocks");
    //   }
    //   }
    //   _algo = algo;
    // }

    void setBlocksize(const size_type mb, 
                      const size_type nb) {
      _mb = mb;
      _nb = nb;
    }

    int reorder(int prunecut = 0) {
      {
        const auto prob_set = _prob_set & ProblemSetMask::Matrix;
        TACHO_TEST_FOR_EXCEPTION(prob_set != ProblemSetMask::Matrix, std::runtime_error, 
                                 "Matrix (A) is not properly set");
      }

      // graph conversion
      Kokkos::View<ordinal_type*,HostSpaceType> rptr("Tacho::Solver::Graph::rptr", _AA.NumRows() + 1);
      Kokkos::View<ordinal_type*,HostSpaceType> cidx("Tacho::Solver::Graph::cidx", _AA.NumNonZeros());
      
      GraphToolsHostType::getGraph(rptr, cidx, _AA);

      GraphToolsHostType_Scotch S;
      S.setGraph(_AA.NumRows(), rptr, cidx);
      S.setSeed();
      S.setTreeLevel();
      S.setStrategy( /**/ SCOTCH_STRATSPEED
                     //|    SCOTCH_STRATLEVELMAX
                     //|    SCOTCH_STRATLEVELMIN
                     |    SCOTCH_STRATLEAFSIMPLE
                     |    SCOTCH_STRATSEPASIMPLE
                     );

      const int treecut = 0; // let's not use treecut anymore
      S.computeOrdering(treecut);

      S.pruneTree(prunecut);

      CrsMatrixBaseHostType AA_scotch("AA_scotch");
      AA_scotch.createConfTo(_AA);

      CrsMatrixTools::copy(AA_scotch,
                           S.PermVector(),
                           S.InvPermVector(),
                           _AA);

      GraphToolsHostType::getGraph(rptr, cidx, AA_scotch);

#define USE_CAMD
#ifdef USE_CAMD
      GraphToolsHostType_CAMD C;
      C.setGraph(AA_scotch.NumRows(),
                 rptr, cidx,
                 S.NumBlocks(),
                 S.RangeVector());
      C.computeOrdering();

      CrsMatrixBaseHostType AA_camd("AA_camd");
      AA_camd.createConfTo(AA_scotch);

      CrsMatrixTools::copy(AA_camd,
                           C.PermVector(),
                           C.InvPermVector(),
                           AA_scotch);
#endif
      {
        const auto m = _AA.NumRows();
        _nblks = S.NumBlocks();
        _perm  = Kokkos::View<ordinal_type*,HostSpaceType>("Tacho::Solver::perm",  m);
        _peri  = Kokkos::View<ordinal_type*,HostSpaceType>("Tacho::Solver::peri",  m);
        _range = Kokkos::View<ordinal_type*,HostSpaceType>("Tacho::Solver::range", m);
        _tree  = Kokkos::View<ordinal_type*,HostSpaceType>("Tacho::Solver::tree",  m);
      
        const auto s_perm = S.PermVector();
        const auto s_peri = S.InvPermVector();
#ifdef USE_CAMD
        const auto c_perm = C.PermVector();
        const auto c_peri = C.InvPermVector();
#endif
        const auto s_range = S.RangeVector();
        const auto s_tree = S.TreeVector();

        for (auto i=0;i<m;++i) {
#ifdef USE_CAMD
          _perm(i)  = c_perm(s_perm(i));
          _peri(i)  = s_peri(c_peri(i));
#else
          _perm(i)  = s_perm(i);
          _peri(i)  = s_peri(i);
#endif
          _range(i) = s_range(i);
          _tree(i)  = s_tree(i);
        }

#ifdef USE_CAMD
        _AA_reorder = AA_camd;
#else
        _AA_reorder = AA_scotch;    
#endif
#undef USE_CAMD
      }

      return 0;
    }

    int analyze() {
      {
        const auto prob_set = _prob_set & ProblemSetMask::Matrix;
        TACHO_TEST_FOR_EXCEPTION(prob_set != ProblemSetMask::Matrix, std::runtime_error, 
                                 "Matrix (A) is not properly set");
      }

      const int rows_per_team = 4096;
      const int fill_level = -1; // direct factorization
      SymbolicFactorizationType::createNonZeroPattern(_FF,
                                                      fill_level,
                                                      Uplo::Upper,
                                                      _AA_reorder,
                                                      rows_per_team);

      CrsMatrixTools::createHierMatrix(_HF,
                                       _FF,
                                       _nblks,
                                       _range,
                                       _tree);

      // ** sparse blocks construction
      {
        const auto nblocks = _HF.NumNonZeros();

        // scan ap
        Kokkos::View<size_type*,HostSpaceType> ap("ap", nblocks + 1);
        ap(0) = 0;
        for (auto k=0;k<nblocks;++k) {
          const auto &block = _HF.Value(k);
          ap(k+1) = ap(k) + block.NumRows();
        }
        
        _rows = Kokkos::View<typename CrsMatrixViewHostType::row_view_type*,HostSpaceType>("Tacho::Solver::rows", ap(nblocks));
        Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, nblocks),
                             [&](const ordinal_type k) {
                               auto &block = _HF.Value(k);
                               block.setRowViewArray(Kokkos::subview(_rows, range_type(ap(k), ap(k+1))));
                             } );
        // this does not work
        // CrsMatrixTools::filterEmptyBlocks(_HF);
      }

      // ** check hier nnz and flat nnz
      {
        size_type nnz = 0;
        const auto nblocks = _HF.NumNonZeros();        
        for (auto k=0;k<nblocks;++k) {
          const auto &block = _HF.Value(k);          
          nnz += block.NumNonZeros();
        }
        TACHO_TEST_FOR_EXCEPTION(nnz != _FF.NumNonZeros(), std::runtime_error, 
                                 "Matrix of blocks does not cover the flat matrix (Scotch tree is not binary)");
      }

      // ** dense block construction
      {
        const auto nblocks = _HF.NumNonZeros();

        // scan ap
        Kokkos::View<size_type*,HostSpaceType> ap("ap", nblocks + 1);
        ap(0) = 0;
        for (auto k=0;k<nblocks;++k) {
          const auto &block = _HF.Value(k);
          const size_type dense_matrix_size = block.NumRows()*block.NumCols();
          ap(k+1) = ap(k) + dense_matrix_size;
        }
        _mats = Kokkos::View<value_type*,HostSpaceType>("Tacho::Solver::mats", ap(nblocks));
        Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, nblocks),
                             [&](const ordinal_type k) {
                               auto &block = _HF.Value(k);
                               block.Flat().setExternalMatrix(block.NumRows(),
                                                              block.NumCols(),
                                                              -1, -1,
                                                              Kokkos::subview(_mats, range_type(ap(k), ap(k+1))));
                               block.copyToFlat();
                             } );
        // HostSpaceType::execution_space::fence();        
        // std::cout <<" file writing \n";
        // std::ofstream file("from_tacho_sym.mtx");
        // MatrixMarket::write(file, _FF);
        // std::cout <<" file writing done \n";
      }

      // ** nullify FF
      //_FF = CrsMatrixBaseHostType();
      
      // ** dense nested blocks construction
      if (_mb) {
        const auto nblocks = _HF.NumNonZeros();

        // scan ap
        Kokkos::View<size_type*,HostSpaceType> ap("ap", nblocks + 1);
        ap(0) = 0;
        for (size_type k=0;k<nblocks;++k) {
          const auto &block = _HF.Value(k);
          ordinal_type hm, hn;
          DenseMatrixTools::getDimensionOfHierMatrix(hm, hn, block.Flat(), _mb, _mb);
          ap(k+1) = ap(k) + hm*hn;
        }

        _blks = Kokkos::View<DenseTaskViewHostType*,HostSpaceType>("DenseBlocksInCrsBlocks", ap(nblocks));
        Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, nblocks),
                             [&](const ordinal_type k) {
                               auto &block = _HF.Value(k);
                               ordinal_type hm, hn;
                               DenseMatrixTools::getDimensionOfHierMatrix(hm, hn, block.Flat(), _mb, _mb);
                               block.Hier().setExternalMatrix(hm, hn,
                                                              -1, -1,
                                                              Kokkos::subview(_blks, range_type(ap(k), ap(k+1))));
                               DenseMatrixTools::getHierMatrix(block.Hier(),
                                                               block.Flat(),
                                                               _mb, _mb);
                             } );
      }

      return 0;
    }

    int factorize() {
      {
        const auto prob_set = _prob_set & ProblemSetMask::Matrix;
        TACHO_TEST_FOR_EXCEPTION(prob_set != ProblemSetMask::Matrix, std::runtime_error, 
                                 "Matrix (A) is not properly set");
      }

      CrsTaskHierViewHostType TF(_HF);
      future_type future;

      if (_mb) // call nested block version
        future = _policy.proc_create_team(Chol<Uplo::Upper,AlgoChol::ByBlocks,Variant::Three>
                                          ::createTaskFunctor(_policy, TF), 0);

      else     // call plain block version
        future = _policy.proc_create_team(Chol<Uplo::Upper,AlgoChol::ByBlocks,Variant::Two>
                                          ::createTaskFunctor(_policy, TF), 0);
      _policy.spawn(future);
      Kokkos::Experimental::wait(_policy);
      TACHO_TEST_FOR_EXCEPTION(future.get(), std::runtime_error, 
                               "Fail to perform CholeskySuperNodesByBlocks");
      


      // // ** dense block to view
      // {
      //   const auto nblocks = _HF.NumNonZeros();
      //   Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, nblocks),
      //                        [&](const ordinal_type k) {
      //                          auto &block = _HF.Value(k);
      //                          block.copyToView();
      //                        } );
      //   // HostSpaceType::execution_space::fence();        
      //   // std::cout <<" file writing \n";
      //   // std::ofstream file("from_tacho_factors.mtx");
      //   // MatrixMarket::write(file, _FF);
      //   // std::cout <<" file writing done \n";
      // }

      return 0;
    }
    
    int solve() {
      {
        const auto prob_set = _prob_set & ProblemSetMask::SetProblem;
        TACHO_TEST_FOR_EXCEPTION(prob_set != ProblemSetMask::SetProblem, std::runtime_error, 
                                 "Problem is not properly set  for AX = B");
      }
      
      const auto m = _BB.NumRows();
      const auto n = _BB.NumCols();
      if (m && n) {
        // create YY if necessary
        {
          const auto work_alloc = _prob_set & ProblemSetMask::WorkspaceAllocated;
          if (!work_alloc) {
            _YY.setLabel("YY");
            _YY.createConfTo(_XX);
            
            DenseMatrixTools::createHierMatrix(_HY, _YY,
                                               _nblks,
                                               _range,
                                               _nb);
            
          }
        }

        // copy BB to XX with permutation
        DenseMatrixTools::applyRowPermutation(_YY, _BB, _perm);
        // for (auto i=0;i<m;++i)
        //   for (auto j=0;j<n;++j)
        //     _YY.Value(_perm(i), j) = _BB.Value(i, j);

        CrsTaskHierViewHostType TF(_HF);
        DenseTaskHierViewHostType TY(_HY);

        auto future_forward_solve
          = _policy.proc_create_team(TriSolve<Uplo::Upper,Trans::ConjTranspose,
                                     AlgoTriSolve::ByBlocks,Variant::Two>
                                     ::createTaskFunctor(_policy,
                                                         Diag::NonUnit, TF, TY), 0);
        _policy.spawn(future_forward_solve);

        auto future_backward_solve
          = _policy.proc_create_team(TriSolve<Uplo::Upper,Trans::NoTranspose,
                                     AlgoTriSolve::ByBlocks,Variant::Two>
                                     ::createTaskFunctor(_policy,
                                                         Diag::NonUnit, TF, TY), 1);

        _policy.add_dependence(future_backward_solve, future_forward_solve);
        _policy.spawn(future_backward_solve);
        
        Kokkos::Experimental::wait(_policy);
        
        TACHO_TEST_FOR_ABORT(future_forward_solve.get(),  "Fail to perform TriSolveSuperNodesByBlocks (forward)");
        TACHO_TEST_FOR_ABORT(future_backward_solve.get(), "Fail to perform TriSolveSuperNodesByBlocks (backward)");

        // copy YY to XX with permutation
        DenseMatrixTools::applyRowPermutation(_XX, _YY, _peri);
        // for (auto i=0;i<m;++i)
        //   for (auto j=0;j<n;++j)
        //     _XX.Value(_peri(i), j) = _YY.Value(i, j);
      }
      
      return 0;
    }

    int check(double &norm, double &err) {

      norm = 0.0;
      err  = 0.0;
      
      const auto m = _BB.NumRows();
      const auto n = _BB.NumCols();
      Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, m),
                           [&](const ordinal_type i) {
                             const auto nnz  = _AA.NumNonZerosInRow(i);
                             const auto cols = _AA.ColsInRow(i);
                             const auto vals = _AA.ValuesInRow(i);
                             
                             for (auto rhs=0;rhs<n;++rhs) {
                               value_type tmp = 0;
                               for (ordinal_type j=0;j<nnz;++j)
                                 tmp += vals(j)*_XX.Value(cols(j), rhs);
                               
                               _YY.Value(i, rhs) = tmp - _BB.Value(i, rhs);
                             }
                           } );
      HostSpaceType::execution_space::fence();

      for (auto j=0;j<n;++j)
        for (auto i=0;i<m;++i) {
          const auto  val = std::abs(_BB.Value(i, j));
          const auto diff = std::abs(_YY.Value(i, j));

          norm += val*val;
          err  += diff*diff;
        }
      
      return 0;
    }
    

  };

}

#endif
