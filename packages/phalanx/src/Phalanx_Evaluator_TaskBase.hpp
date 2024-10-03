// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHALANX_EVALUATOR_TASK_BASE_HPP
#define PHALANX_EVALUATOR_TASK_BASE_HPP

#ifdef PHX_ENABLE_KOKKOS_AMT

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
// #include "Teuchos_TypeNameTraits.hpp"
// #include <typeinfo>

namespace PHX {

  
  //! Generic task wrapper to use a PHX::Evaluator with Kokkos AMT
  template<typename Space,typename Functor>
  struct TaskWrap {
    
    typedef void value_type;
    typedef Kokkos::TaskScheduler<Space> policy_type;
    
    const int work_size;
    const Functor functor;

    TaskWrap(const int ws,const Functor& f) : work_size(ws),functor(f) {}
    
    void operator() (typename policy_type::member_type & member)
    {
      // std::cout << "In functor: "
      //           << Teuchos::demangleName(typeid(*this).name())
      //           << "\n  team_rank = " << member.team_rank()
      //           << ", team_size = " << member.team_size()
      //           << std::endl;
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member,work_size),functor);
    }
  };

  
  template<typename Traits,typename Derived>
  struct TaskBase : public PHX::EvaluatorWithBaseImpl<Traits> {

    virtual Kokkos::Future<void,PHX::exec_space>
    createTask(Kokkos::TaskScheduler<PHX::exec_space>& policy,
	       const int& work_size,
               const std::vector<Kokkos::Future<void,PHX::exec_space>>& dependent_futures,
	       typename Traits::EvalData ) override
    {
      if (dependent_futures.size() == 0)
        return policy.host_spawn(PHX::TaskWrap<PHX::exec_space,Derived>(work_size,*dynamic_cast<Derived*>(this)),Kokkos::TaskTeam);
      else if (dependent_futures.size() == 1)
        return policy.host_spawn(PHX::TaskWrap<PHX::exec_space,Derived>(work_size,*dynamic_cast<Derived*>(this)),Kokkos::TaskTeam,dependent_futures[0]);

      auto aggregate_future = policy.when_all(dependent_futures.size(),dependent_futures.data());
      return policy.host_spawn(PHX::TaskWrap<PHX::exec_space,Derived>(work_size,*dynamic_cast<Derived*>(this)),Kokkos::TaskTeam,aggregate_future);
    }

    //! Returns the size of the task functor in bytes
    unsigned taskSize() const override
    {
      return sizeof(PHX::TaskWrap<PHX::exec_space,Derived>);
    }
  };

}

#endif

#endif
