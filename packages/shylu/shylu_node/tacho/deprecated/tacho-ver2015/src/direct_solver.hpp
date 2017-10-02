#pragma once
#ifndef __DIRECT_SOLVER_HPP__
#define __DIRECT_SOLVER_HPP__

/// \file direct_solver.hpp
/// \brief Wrapper for direct solver.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "dense_matrix_base.hpp"
#include "dense_matrix_view.hpp"

#include "graph_helper_scotch.hpp"
#include "symbolic_factor_helper.hpp"
#include "crs_matrix_helper.hpp"
#include "dense_matrix_helper.hpp"

#include "task_view.hpp"
#include "task_factory.hpp"

#include "chol.hpp"
#include "tri_solve.hpp"


namespace Tacho {

  using namespace std;

  template<typename ValueType,
           typename OrdinalType,
           typename SizeType = OrdinalType,
           typename SpaceType = void,
           typename MemoryTraits = void>
  class DirectSolver : public Disp {
  public:
    typedef ValueType    value_type;
    typedef OrdinalType  ordinal_type;
    typedef SpaceType    space_type;
    typedef SizeType     size_type;
    typedef MemoryTraits memory_traits;

    typedef TaskFactory<Kokkos::Experimental::TaskPolicy<space_type>,
                        Kokkos::Experimental::Future<int,space_type> > TaskFactoryType;

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,space_type,memory_traits> CrsMatrixBaseType;
    typedef GraphHelper_Scotch<CrsMatrixBaseType> GraphHelperType;
    typedef SymbolicFactorHelper<CrsMatrixBaseType> SymbolicFactorHelperType;

    typedef CrsMatrixView<CrsMatrixBaseType> CrsMatrixViewType;
    typedef TaskView<CrsMatrixViewType,TaskFactoryType> CrsTaskViewType;

    typedef CrsMatrixBase<CrsTaskViewType,ordinal_type,size_type,SpaceType,MemoryTraits> CrsHierMatrixBaseType;

    typedef CrsMatrixView<CrsHierMatrixBaseType> CrsHierMatrixViewType;
    typedef TaskView<CrsHierMatrixViewType,TaskFactoryType> CrsHierTaskViewType;

    typedef DenseMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> DenseMatrixBaseType;

    typedef DenseMatrixView<DenseMatrixBaseType> DenseMatrixViewType;
    typedef TaskView<DenseMatrixViewType,TaskFactoryType> DenseTaskViewType;

    typedef DenseMatrixBase<DenseTaskViewType,ordinal_type,size_type,SpaceType,MemoryTraits> DenseHierMatrixBaseType;

    typedef DenseMatrixView<DenseHierMatrixBaseType> DenseHierMatrixViewType;
    typedef TaskView<DenseHierMatrixViewType,TaskFactoryType> DenseHierTaskViewType;

    typedef typename TaskFactoryType::policy_type policy_type;

  private:
    string _label;
    policy_type _policy;

    CrsMatrixBaseType _AA, _PA, _UU;
    DenseMatrixBaseType _BB;

    GraphHelperType _S;

    size_type _max_concurrency, _max_task_size, _max_task_dependence, _team_size;
    bool _team_interface;

    ordinal_type _prunecut;
    size_type _league_size;

    size_type _mb, _nb;

    CrsHierMatrixBaseType _HU;
    DenseHierMatrixBaseType _HB;

  public:
    DirectSolver() 
      : _label("DirectSolver"),
        _AA(),
        _BB(),
        _S()
    { }

    DirectSolver(const string label,
                 const CrsMatrixBaseType AA)
      : _label(label),
        _AA(AA),
        _BB(),
        _S(AA, 0)
    { }

    DirectSolver(const DirectSolver &b) = default;

    KOKKOS_INLINE_FUNCTION
    void setB(const DenseMatrixBaseType BB) { _BB = BB; }

    KOKKOS_INLINE_FUNCTION
    void permuteB(const bool is_inverse = false) { 
      DenseMatrixBaseType tmp;
      tmp.copy(_BB);
      auto p = (is_inverse ? _S.PermVector() : _S.InvPermVector());
      for (ordinal_type j=0;j<_BB.NumCols();++j)
        for (ordinal_type i=0;i<_BB.NumRows();++i) 
          _BB.Value(p[i], j) = tmp.Value(i, j); 
    }

    KOKKOS_INLINE_FUNCTION
    void setLabel(const string label) { _label = label; }

    KOKKOS_INLINE_FUNCTION
    string Label() const { return _label; }

    KOKKOS_INLINE_FUNCTION
    void setDefaultParameters(const ordinal_type nthreads) {
      // task policy
      _max_concurrency = 1024;
      _max_task_size = 3*sizeof(CrsTaskViewType)+128;
      _max_task_dependence = 3;
      _team_size = 1;

      _team_interface = true;

      // graph ordering
      _prunecut = 0;

      _league_size = nthreads;

      _nb = 128;
    }

    KOKKOS_INLINE_FUNCTION
    int init() {
      _policy = policy_type(_max_concurrency,
                            _max_task_size,
                            _max_task_dependence,
                            _team_size);
      
      TaskFactoryType::setUseTeamInterface(_team_interface);
      TaskFactoryType::setMaxTaskDependence(_max_task_dependence);
      TaskFactoryType::setPolicy(&_policy);

      return 0;
    }

    KOKKOS_INLINE_FUNCTION
    int reorder() {
      _S.computeOrdering();
      _S.pruneTree(_prunecut);
      _PA.copy(_S.PermVector(), _S.InvPermVector(), _AA);

      return 0;
    }

    KOKKOS_INLINE_FUNCTION
    int analyze() {
      SymbolicFactorHelperType F(_PA, 1);//_league_size);
      F.createNonZeroPattern(Uplo::Upper, _UU);
      CrsMatrixHelper::flat2hier(Uplo::Upper, _UU, _HU,
                                 _S.NumBlocks(),
                                 _S.RangeVector(),
                                 _S.TreeVector());

      for (ordinal_type k=0;k<_HU.NumNonZeros();++k)
        _HU.Value(k).fillRowViewArray();

      return 0;
    }

    KOKKOS_INLINE_FUNCTION
    int factorize() {
      CrsHierTaskViewType TU(&_HU);
      auto future = TaskFactoryType::Policy().create_team
        (Chol<Uplo::Upper,AlgoChol::ByBlocks>::
         TaskFunctor<CrsHierTaskViewType>(TU), 0);
      TaskFactoryType::Policy().spawn(future);
      Kokkos::Experimental::wait(TaskFactoryType::Policy());
      return 0;
    }

    int solve() {
      DenseMatrixHelper::flat2hier(_BB, _HB,
                                   _S.NumBlocks(),
                                   _S.RangeVector(),
                                   _nb);

      CrsHierTaskViewType TU(&_HU);
      DenseHierTaskViewType TB(&_HB);

      auto future_forward_solve = TaskFactoryType::Policy().create_team
        (TriSolve<Uplo::Upper,Trans::ConjTranspose,AlgoTriSolve::ByBlocks>
         ::TaskFunctor<CrsHierTaskViewType,DenseHierTaskViewType>
         (Diag::NonUnit, TU, TB), 0);

      auto future_backward_solve = TaskFactoryType::Policy().create_team
        (TriSolve<Uplo::Upper,Trans::NoTranspose,AlgoTriSolve::ByBlocks>
         ::TaskFunctor<CrsHierTaskViewType,DenseHierTaskViewType>
         (Diag::NonUnit, TU, TB), 1);

      TaskFactoryType::Policy().spawn(future_forward_solve);
      TaskFactoryType::Policy().add_dependence(future_backward_solve, future_forward_solve);
      TaskFactoryType::Policy().spawn(future_backward_solve);

      Kokkos::Experimental::wait(TaskFactoryType::Policy());

      return 0;
    }



    ostream& showMe(ostream &os) const {
      return os;
    }

  };

}

#endif
