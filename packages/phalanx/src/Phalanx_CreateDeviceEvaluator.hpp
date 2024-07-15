// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHALANX_CREATE_DEVICE_EVALUATOR_HPP
#define PHALANX_CREATE_DEVICE_EVALUATOR_HPP

#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_DeviceEvaluator.hpp"
#include <memory>

namespace PHX {

  /** \brief Function to allow evaluators to allocate a
      DeviceEvaluator on the device. The object is lambda captured on
      device and the copy constructor is used to allocate on device
      with Kokkos::malloc so that the vtable is allocated on the
      device and can be called on the device. EvaluatorType must be
      the concrete underlying type, not a base class.
  */
  template<typename DeviceEvaluatorType, typename Traits, typename ExecutionSpace, typename MemorySpace, typename ...CtorArgs>
  PHX::DeviceEvaluator<Traits>* createDeviceEvaluator(CtorArgs&&... args)
  {
    DeviceEvaluatorType evaluator_to_clone_on_device(args...);
    DeviceEvaluator<Traits>* e = static_cast<DeviceEvaluatorType*>(Kokkos::kokkos_malloc<MemorySpace>(sizeof(DeviceEvaluatorType)));
    Kokkos::parallel_for(Kokkos::RangePolicy<ExecutionSpace>(0,1),
                         KOKKOS_LAMBDA (const int i)
                         {
                           new (e) DeviceEvaluatorType(evaluator_to_clone_on_device);
                         });
    ExecutionSpace().fence();
    return e;
  }

  /** \brief Kokkos functor that wraps a users evaluator functor (Decorator design pattern).
      This eliminates the user having to store/assign workset data.
  */
  template<typename DeviceEvaluatorType>
  struct DevEvalWrapper {
    mutable DeviceEvaluatorType device_evaluator_;
    typename std::remove_reference<typename DeviceEvaluatorType::traits::EvalData>::type data_;

    DevEvalWrapper(const DeviceEvaluatorType& de, typename DeviceEvaluatorType::traits::EvalData data)
      : device_evaluator_(de), data_(data) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) const
    {
      device_evaluator_.prepareForRecompute(team, data_);
      device_evaluator_.evaluate(team, data_);
    }
  };

  /** \brief Stand alone function for creating a DevEvalWrapper.
      Eliminates user needing to define template parameter.
  */
  template<typename DeviceEvaluatorType>
  PHX::DevEvalWrapper<DeviceEvaluatorType> 
  make_dev_eval(const DeviceEvaluatorType& de, typename DeviceEvaluatorType::traits::EvalData data)
  {
    return PHX::DevEvalWrapper<DeviceEvaluatorType>(de,data);
  }

}

#endif
