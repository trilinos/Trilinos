// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_ETI_HPP
#define PHX_ETI_HPP

#include "MyTraits.hpp"

// Macros for explicit template instatiation

#define PHX_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL(name) \
  template class name<PHX::MyTraits::Residual, PHX::MyTraits>; 

#define PHX_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN(name) \
  template class name<PHX::MyTraits::Jacobian, PHX::MyTraits>; 

/*
#define PHX_INSTANTIATE_TEMPLATE_CLASS_JV(name) \
  template class name<PHX::MyTraits::Jv, PHX::MyTraits>; 

#define PHX_INSTANTIATE_TEMPLATE_CLASS(name) \
  PHX_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL(name) \
  PHX_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN(name) \
  PHX_INSTANTIATE_TEMPLATE_CLASS_JV(name)
*/

#define PHX_INSTANTIATE_TEMPLATE_CLASS(name)    \
  PHX_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL(name) \
  PHX_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN(name)

#endif
