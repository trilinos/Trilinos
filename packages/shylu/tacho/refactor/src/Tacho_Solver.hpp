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
    
    // Policy on device
    typedef Kokkos::TaskPolicy<DeviceSpaceType> PolicyType;
    
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
    typedef Kokkos::Future<int,HostSpaceType> future_type;

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

    // cross over parameter
    ordinal_type _algo_hier_minsize; //, _algo_flat_maxsize;
    
    // input matrices
    CrsMatrixBaseHostType _AA;
    DenseMatrixBaseHostType _BB, _XX, _YY;

    // extract permutation and blocking information
    CrsMatrixBaseHostType _PAA;    
    ordinal_type _nblks, _maxrange;
    Kokkos::View<ordinal_type*,HostSpaceType> _perm, _peri, _range, _tree;
    
    // block hierachical matrix factored
    CrsHierBaseHostType _HF;
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

  public:
    Solver(const std::string label = "Tacho::Solver")
      : _label(label)
    {}

    Solver(const Solver &b) = default;

    CrsMatrixBaseHostType   getA() const { return _AA; }
    DenseMatrixBaseHostType getB() const { return _BB; }
    DenseMatrixBaseHostType getX() const { return _XX; }

    CrsHierBaseHostType getHierU() const { return _HF; }
    ordinal_type getMaxRangeSize() const { return _maxrange; }

    void setPolicy(const size_type max_concurrency = 25000, 
                   const size_type memory_pool_grain_size = 16) {
      const size_type max_task_size = (3*sizeof(CrsTaskViewHostType)+sizeof(PolicyType)+128);

      _policy = PolicyType(typename PolicyType::memory_space(),
                           max_task_size*max_concurrency,
                           memory_pool_grain_size );
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

    void setBlocksize(const ordinal_type mb, 
                      const ordinal_type nb) {
      _mb = mb;
      _nb = nb;
    }

    void setCrossOverSize(const ordinal_type algo_hier_minsize) {
      _algo_hier_minsize = algo_hier_minsize;
    }

    int reorder(int treecut = 4, int prunecut = 0) {
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
                     |    SCOTCH_STRATLEVELMAX
                     |    SCOTCH_STRATLEVELMIN
                     |    SCOTCH_STRATLEAFSIMPLE
                     |    SCOTCH_STRATSEPASIMPLE
                     );

      S.computeOrdering(treecut);
      S.pruneTree(prunecut);

      CrsMatrixBaseHostType AA_scotch("AA_scotch");
      AA_scotch.createConfTo(_AA);

      CrsMatrixTools::copy(AA_scotch,
                           S.PermVector(),
                           S.InvPermVector(),
                           _AA);

#define USE_CAMD
#ifdef USE_CAMD
      GraphToolsHostType::getGraph(rptr, cidx, AA_scotch);
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

        _maxrange = 0;        
        for (auto i=0;i<_nblks;++i)
          _maxrange = Util::max(_maxrange, _range(i+1)-_range(i)); 

#ifdef USE_CAMD
        _PAA = AA_camd;
#else
        _PAA = AA_scotch;    
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
      
      CrsHierBaseHostType PHA("PHA");
      { 
        // find fills in block matrix
        typename CrsHierBaseHostType::size_type_array ap;
        typename CrsHierBaseHostType::ordinal_type_array aj;
        {
          // to apply symbolic factorization on the matrix of blocks, we need a graph of
          // the block matrix (full matrix).
          const bool full = true;
          CrsMatrixTools::createHierMatrix(PHA, _PAA, _nblks, _range, _tree, full);
          
          {
            const size_type nnz = PHA.NumNonZeros();
            Kokkos::View<size_type*,HostSpaceType> offs("offs", nnz + 1);
            offs(0) = 0;
            for (size_type k=0;k<nnz;++k) {
              const auto &block = PHA.Value(k);
              offs(k+1) = offs(k) + block.NumRows();
            }
            _rows = Kokkos::View<typename CrsMatrixViewHostType::row_view_type*,HostSpaceType>("RowViewsInBlocks", offs(nnz));
            Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, nnz),
                                 [&](const size_type k) {
                                   auto &block = PHA.Value(k);
                                   block.setRowViewArray(Kokkos::subview(_rows, range_type(offs(k), offs(k+1))));
                                 } );
          }
          CrsMatrixTools::filterEmptyBlocks(PHA);
          CrsMatrixTools::createSymbolicStructure(ap, aj, PHA);
        }
        
        // fill block information
        typename CrsHierBaseHostType::value_type_array ax("ax", aj.dimension(0));
        {
          Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, _nblks),
                               [&](const ordinal_type i) {
                                 const auto beg = ap(i);
                                 const auto end = ap(i+1);
                                 for (auto idx=beg;idx<end;++idx) {
                                   const auto j = aj(idx);
                                   ax(idx).setView(_PAA, _range(i), (_range(i+1) - _range(i)),
                                                   /**/  _range(j), (_range(j+1) - _range(j)));
                                 }
                               } );
        }
        
        // construct hierachical matrix
        {
          const size_type nnz = aj.dimension(0);
          const auto ap_begin = Kokkos::subview(ap, range_type(0,_nblks  ));
          const auto ap_end   = Kokkos::subview(ap, range_type(1,_nblks+1));
          _HF = CrsHierBaseHostType("HF", _nblks, _nblks, nnz, ap_begin, ap_end, aj, ax);
          _HF.setNumNonZeros();
        }
      }

      // copy row view structure to block symbolic factors
      {
        Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, _nblks),
                             [&](const ordinal_type i) {
                               const auto cols_a = PHA.ColsInRow(i);
                               const auto vals_a = PHA.ValuesInRow(i);
                                 
                               const auto cols_f = _HF.ColsInRow(i);
                               const auto vals_f = _HF.ValuesInRow(i);

                               const size_type nnz_in_a = PHA.NumNonZerosInRow(i);
                               const size_type nnz_in_f = _HF.NumNonZerosInRow(i);

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

      // allocate super nodal structure
      {
        const size_type nnz = _HF.NumNonZeros();
        {
          Kokkos::View<size_type*,HostSpaceType> offs("offs", nnz + 1);
          offs(0) = 0;
          for (size_type k=0;k<nnz;++k) {
            const auto &block = _HF.Value(k);
            offs(k+1) = offs(k) + block.NumRows()*block.NumCols();
          }
          
          _mats = Kokkos::View<value_type*,HostSpaceType>("MatsInnz", offs(nnz));
          Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, nnz),
                               [&](const ordinal_type k) {
                                 auto &block = _HF.Value(k);
                                 block.Flat().setExternalMatrix(block.NumRows(),
                                                                block.NumCols(),
                                                                -1, -1,
                                                                Kokkos::subview(_mats, range_type(offs(k), offs(k+1))));
                               } );
        }
        if (_mb) {
          Kokkos::View<size_type*,HostSpaceType> offs("offs", nnz + 1);
          offs(0) = 0;
          for (size_type k=0;k<nnz;++k) {
            const auto &block = _HF.Value(k);
            ordinal_type hm, hn;
            DenseMatrixTools::getDimensionOfHierMatrix(hm, hn, block.Flat(), _mb, _mb);
            offs(k+1) = offs(k) + hm*hn;
          }
          
          _blks = Kokkos::View<DenseTaskViewHostType*,HostSpaceType>("DenseBlocksInCrsBlocks", offs(nnz));
          Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, nnz),
                               [&](const ordinal_type k) {
                                 auto &block = _HF.Value(k);
                                 ordinal_type hm, hn;
                                 DenseMatrixTools::getDimensionOfHierMatrix(hm, hn, block.Flat(), _mb, _mb);
                                 block.Hier().setExternalMatrix(hm, hn,
                                                                -1, -1,
                                                                Kokkos::subview(_blks, range_type(offs(k), offs(k+1))));
                                 Impl::DenseMatrixTools::Serial
                                   ::getHierMatrix(block.Hier(),
                                                   block.Flat(),
                                                   _mb, _mb);
                               } );
        }
      }

      return 0;
    }

    int getFactorizationAlgorithmVariant() const {
      int variant = 0;
      if (_mb > 0) {
        const ordinal_type m = _PAA.NumRows();
        const ordinal_type hm = _maxrange/_mb + 1;

        variant = (hm >= _algo_hier_minsize ? Variant::Three : Variant::Two);
      } else {
        variant = Variant::Two;
      }
      return variant;
    }
    
    int factorize(const int var = 0) {
      {
        const auto prob_set = _prob_set & ProblemSetMask::Matrix;
        TACHO_TEST_FOR_EXCEPTION(prob_set != ProblemSetMask::Matrix, std::runtime_error, 
                                 "Matrix (A) is not properly set");
      }

      // copy the matrix 
      {
        Kokkos::deep_copy(_mats, 0);
        
        const size_type nnz = _HF.NumNonZeros();
        Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, nnz),
                             [&](const ordinal_type k) {
                               auto &block = _HF.Value(k);
                               block.copyToFlat();
                             } );
      }
      
      // factorize
      {
        int variant = getFactorizationAlgorithmVariant();
        if (var) variant = var;

        CrsTaskHierViewHostType TF(_HF);
        future_type future;

        switch (variant) {
        case Variant::Three: 
          future = _policy.host_spawn(Chol<Uplo::Upper,AlgoChol::ByBlocks,Variant::Three>
                                      ::createTaskFunctor(_policy, TF));
          break;
        case Variant::Two: 
          future = _policy.host_spawn(Chol<Uplo::Upper,AlgoChol::ByBlocks,Variant::Two>
                                      ::createTaskFunctor(_policy, TF));
          break;
        }
        TACHO_TEST_FOR_EXCEPTION(future.is_null(), std::runtime_error,
                                 ">> host_spawn returns a null future");
        
        Kokkos::wait(_policy);
        TACHO_TEST_FOR_EXCEPTION(future.get(), std::runtime_error, 
                                 "Fail to perform CholeskySuperNodesByBlocks");
      }
      
      return 0;
    }
    
    int solve() {
      {
        const auto prob_set = _prob_set & ProblemSetMask::SetProblem;
        TACHO_TEST_FOR_EXCEPTION(prob_set != ProblemSetMask::SetProblem, std::runtime_error, 
                                 "Problem is not properly set  for AX = B");
      }
      
      const ordinal_type m = _BB.NumRows();
      const ordinal_type n = _BB.NumCols();
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
          = _policy.host_spawn(TriSolve<Uplo::Upper,Trans::ConjTranspose,
                               AlgoTriSolve::ByBlocks,Variant::Two>
                               ::createTaskFunctor(_policy,
                                                   Diag::NonUnit, TF, TY));
        TACHO_TEST_FOR_EXCEPTION(future_forward_solve.is_null(), std::runtime_error,
                                 ">> host_spawn returns a null future");
        
        auto future_backward_solve
          = _policy.host_spawn(TriSolve<Uplo::Upper,Trans::NoTranspose,
                               AlgoTriSolve::ByBlocks,Variant::Two>
                               ::createTaskFunctor(_policy,
                                                   Diag::NonUnit, TF, TY),
                               future_forward_solve);
        TACHO_TEST_FOR_EXCEPTION(future_backward_solve.is_null(), std::runtime_error,
                                 ">> host_spawn returns a null future");

        Kokkos::wait(_policy);
        
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
