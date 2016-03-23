#ifndef __TACHO_TASK_FACTORY_HPP__
#define __TACHO_TASK_FACTORY_HPP__

/// \file Tacho_TaskFactory.hpp
/// \brief Front interface that wrap Kokkos tasking with safe guards.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho { 

  /// \class TaskFactory
  /// \brief Front interface to Kokkos tasking.
  class TaskFactory {
  public:
    template<typename FutureType,
             typename PolicyType,
             typename TaskFunctorType>
    static KOKKOS_INLINE_FUNCTION
    FutureType create(PolicyType &policy, 
                      const TaskFunctorType &func, 
                      const int max_task_dependence) {
      FutureType f = policy.task_create_team(func, max_task_dependence);
      TACHO_TEST_FOR_ABORT( f.is_null(), ">> Creating team task is failed, out of memory" );

      return f;
    }
    
    template<typename PolicyType,
             typename FutureType>
    static KOKKOS_INLINE_FUNCTION
    void spawn(PolicyType &policy, 
               const FutureType &obj, 
               const bool priority = false) {
      policy.spawn(obj, priority);
    }
    
    template<typename PolicyType, 
             typename TaskFunctorType>
    static KOKKOS_INLINE_FUNCTION
    void respawn(PolicyType &policy, 
                 TaskFunctorType *func,
                 const bool priority = false) {
      policy.respawn(func, priority);
    }
    
    template<typename PolicyType, 
             typename FutureTypeAfter,
             typename FutureTypeBefore>
    static KOKKOS_INLINE_FUNCTION
    void depend(PolicyType &policy, 
                const FutureTypeAfter  &after, 
                const FutureTypeBefore &before) {
      policy.add_dependence(after, before);
    }
    
    template<typename PolicyType, 
             typename TaskFunctorType,
             typename FutureType>
    static KOKKOS_INLINE_FUNCTION
    void depend(PolicyType &policy, 
                TaskFunctorType *after, 
                const FutureType &before) {
      policy.add_dependence(after, before);
    }

    template<typename PolicyType, 
             typename TaskFunctorType>
    static KOKKOS_INLINE_FUNCTION
    void clear(PolicyType &policy, 
               TaskFunctorType *func) {
      policy.clear_dependence(func);
    }

  };
}

#endif
