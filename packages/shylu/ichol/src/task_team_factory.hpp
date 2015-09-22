#pragma once
#ifndef __TASK_TEAM_FACTORY_HPP__
#define __TASK_TEAM_FACTORY_HPP__

/// \file task_factory.hpp
/// \brief A wrapper for task policy and future with a provided space type.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example {

  using namespace std;

  /// \class TaskFactory
  /// \brief Minimal interface to Kokkos tasking.
  ///
  /// TaskFactory is attached to blocks as a template argument in order to
  /// create and manage tasking future objects. Note that policy (shared
  /// pointer to the task generator) is not a member object in this class.
  /// This class includes minimum interface for tasking with type decralation
  /// of the task policy and template alias of future so that future objects
  /// generated in this class will match to their policy and its execution space.
  ///
  template<typename PolicyType,
           typename FutureType,
           template<typename,typename> class ThreadLoopRegionType>
  class TaskTeamFactory : public TaskFactory<PolicyType,FutureType>,
                          public TeamFactory<PolicyType,ThreadLoopRegionType> {
  public:
    typedef PolicyType policy_type;
    typedef FutureType future_type;
  };

}

#endif
