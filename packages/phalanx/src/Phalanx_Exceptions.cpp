// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Phalanx_Exceptions.hpp"

namespace PHX {

  circular_dag_exception::circular_dag_exception(const std::string message) : 
    std::runtime_error(message) {}

  missing_evaluator_exception::
  missing_evaluator_exception(const std::string message) : 
    std::runtime_error(message) {}

  multiple_evaluator_for_field_exception::
  multiple_evaluator_for_field_exception(const std::string message) : 
    std::runtime_error(message) {}

}
