// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHALANX_DEVICE_EVALUATOR_HPP
#define PHALANX_DEVICE_EVALUATOR_HPP

#include "Kokkos_Core.hpp"

namespace PHX {
  
  //! Pure virtual interface for instantiating an evaluator on device 
  template<typename Traits>
  struct DeviceEvaluator {
    using team_policy = Kokkos::TeamPolicy<PHX::exec_space>;
    using member_type = team_policy::member_type;
    using traits = Traits;
    
    KOKKOS_DEFAULTED_FUNCTION virtual ~DeviceEvaluator() = default;

    //! Used to bind EvalData objects to functor
    KOKKOS_FUNCTION virtual void
    prepareForRecompute(const member_type& , typename Traits::EvalData ) {}

    //! Performs the evaluation
    KOKKOS_FUNCTION virtual void
    evaluate(const member_type& team, typename Traits::EvalData d) = 0;
  };

  //! Struct for holding evaluator pointers in a Kokkos::View
  template<typename Traits>
  struct DeviceEvaluatorPtr {
    DeviceEvaluator<Traits>* ptr;
  };

}

#endif
