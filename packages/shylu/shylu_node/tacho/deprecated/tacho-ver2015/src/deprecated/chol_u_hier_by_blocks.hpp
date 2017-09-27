#pragma once
#ifndef __CHOL_U_HIER_BY_BLOCKS_HPP__
#define __CHOL_U_HIER_BY_BLOCKS_HPP__

/// \file chol_u_hier_by_blocks.hpp
/// \brief Hierarchical Cholesky factorization by-blocks
/// \author Kyungjoo Kim (kyukim@sandia.gov)

/// This perform Cholesky factorization with subblocks defined in a given block (hierarchical view)

// basic utils
#include "util.hpp"
#include "control.hpp"
#include "partition.hpp"

namespace Tacho { 

  using namespace std;

  // specialization for different task generation in right looking by-blocks algorithm
  // =================================================================================
  template<int ArgVariant, template<int,int> class ControlType>
  class Chol<Uplo::Upper,AlgoChol::HierByBlocks,ArgVariant,ControlType> {
  public:

    // function interface
    // ==================
    template<typename ExecViewType>
    KOKKOS_INLINE_FUNCTION
    static int sync(typename ExecViewType::policy_type &policy,
                    const typename ExecViewType::policy_type::member_type &member,
                    ExecViewType &A) {
      // each algorithm and its variants should be specialized
      ERROR(MSG_INVALID_TEMPLATE_ARGS);
      return -1;
    }

    // task-data parallel interface
    // ============================
    template<typename ExecViewType>
    class TaskFunctor {
    public:
      typedef typename ExecViewType::policy_type policy_type;
      typedef typename policy_type::member_type member_type;

      typedef typename ExecViewType::task_factory_type task_factory_type;
      typedef typename ExecViewType::future_type future_type;

      typedef int value_type;
      
    private:
      ExecViewType _A;
      
      policy_type &_policy;
      
    public:
      TaskFunctor(const ExecViewType A)
        : _A(A),
          _policy(task_factory_type::Policy())
      { } 
      
      string Label() const { return "Chol"; }

      // task-data execution
      void apply(const member_type &member, value_type &r_val) {
        // first check if it has subblocks
        if (_A.hasSubBlockView()) {
          typedef typename ExecViewType::subblock_type subblock_type;

          if (!_respawn) {
            future_type f = task_factory_type::create(policy,
                                                      typename Chol<Uplo::Upper,
                                                      CtrlDetail(ControlType,AlgoChol::HierByBlocks,ArgVariant,CholByBlocks)>
                                                      template TaskFunctor<subblock_type>(_A.SubBlockView()));
            task_factory_type::spawn(policy, f);
            
            _respawn = true;
            task_factory_type::clearDependence(this);
            task_factory_type::addDependence(policy, this, f);
            task_factory_type::respawn(policy, this);
          } else {
            future_type f = task_factory_type::getDependence(policy, this, 0);
            if (f.get() == 0) { // recursive tasks are successfully generated
              _respawn = false;

              // get subblock future to make sure they are executed
              // and perform additional repacking operation  
              sync(policy, member, _A.SubBlockView());
            } else {
              // error handling
            }
          }

        } else {
          // work on A
          r_val = Chol<Uplo::Upper,CtrlDetail(ControlType,AlgoChol::HierByBlocs,ArgVariant,CholScalar)>
            ::invoke(_policy, member, _A);
        }
      }
      
      // task execution
      void apply(value_type &r_val) {
        apply(_policy.member_single(), r_val);
      }

    };

  };
}

#endif
