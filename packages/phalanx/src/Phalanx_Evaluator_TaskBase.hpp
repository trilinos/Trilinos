#ifndef PHALANX_EVALAUTOR_TASK_BASE_HPP
#define PHALANX_EVALAUTOR_TASK_BASE_HPP

#ifdef PHX_ENABLE_KOKKOS_AMT

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"

namespace PHX_example {

  //! Generic task wrapper to use a PHX::Evaluator with Kokkos AMT
  template<typename Space,typename Functor>
  struct TaskWrap {
    
    typedef void value_type;
    typedef Kokkos::TaskPolicy<Space> policy_type;
    
    const int work_size;
    const Functor functor;

    TaskWrap(const int ws,const Functor& f) : work_size(ws),functor(f) {}
    
    void operator() (typename policy_type::member_type & member)
    {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member,work_size),functor);
    }
  };

  
  template<typename Traits,typename Derived>
  struct TaskBase : public PHX::EvaluatorWithBaseImpl<Traits> {

    virtual Kokkos::Future<void,PHX::Device::execution_space>
    createTask(Kokkos::TaskPolicy<PHX::Device::execution_space>& policy,
	       const int& work_size,
               const std::vector<Kokkos::Future<void,PHX::Device::execution_space>>& dependent_futures,
	       typename Traits::EvalData ) override
    {
      if (dependent_futures.size() == 0)
        return policy.host_spawn(PHX_example::TaskWrap<PHX::Device::execution_space,Derived>(work_size,*dynamic_cast<Derived*>(this)),Kokkos::TaskTeam);
      else if (dependent_futures.size() == 1)
        return policy.host_spawn(PHX_example::TaskWrap<PHX::Device::execution_space,Derived>(work_size,*dynamic_cast<Derived*>(this)),Kokkos::TaskTeam,dependent_futures[0]);

      auto aggregate_future = policy.when_all(dependent_futures.size(),dependent_futures.data());
      return policy.host_spawn(PHX_example::TaskWrap<PHX::Device::execution_space,Derived>(work_size,*dynamic_cast<Derived*>(this)),Kokkos::TaskTeam,aggregate_future);
    }

    //! Returns the size of the task functor in bytes
    unsigned taskSize() const override
    {
      return sizeof(PHX_example::TaskWrap<PHX::Device::execution_space,Derived>);
    }
  };

}

#endif

#endif
