// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EXPLICIT_TEMPLATE_INSTANTIATION_HPP
#define PANZER_EXPLICIT_TEMPLATE_INSTANTIATION_HPP

#include "Panzer_Traits.hpp"

// ONE template argument 
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_ONE_T(name) \
  template class name<panzer::Traits::Residual>; 

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_TANGENT_ONE_T(name) \
  template class name<panzer::Traits::Tangent>; 

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_ONE_T(name) \
  template class name<panzer::Traits::Jacobian>; 

#ifdef Panzer_BUILD_HESSIAN_SUPPORT 
  #define PANZER_INSTANTIATE_TEMPLATE_CLASS_HESSIAN_ONE_T(name) \
    template class name<panzer::Traits::Hessian>;
#else 
  #define PANZER_INSTANTIATE_TEMPLATE_CLASS_HESSIAN_ONE_T(name)
#endif

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_ONE_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_ONE_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_TANGENT_ONE_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_ONE_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_HESSIAN_ONE_T(name)

// TWO template arguments
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_TWO_T(name) \
  template class name<panzer::Traits::Residual, panzer::Traits>; 

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_TANGENT_TWO_T(name) \
  template class name<panzer::Traits::Tangent, panzer::Traits>; 

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_TWO_T(name) \
  template class name<panzer::Traits::Jacobian, panzer::Traits>; 

#ifdef Panzer_BUILD_HESSIAN_SUPPORT 
  #define PANZER_INSTANTIATE_TEMPLATE_CLASS_HESSIAN_TWO_T(name) \
    template class name<panzer::Traits::Hessian, panzer::Traits>;
#else
  #define PANZER_INSTANTIATE_TEMPLATE_CLASS_HESSIAN_TWO_T(name) 
#endif

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_TWO_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_TANGENT_TWO_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_TWO_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_HESSIAN_TWO_T(name)

// THREE (one user defined) template arguments
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_THREE_T(name,ExtraT) \
  template class name<panzer::Traits::Residual, panzer::Traits,ExtraT>; 

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_TANGENT_THREE_T(name,ExtraT) \
  template class name<panzer::Traits::Tangent, panzer::Traits,ExtraT>; 

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_THREE_T(name,ExtraT) \
  template class name<panzer::Traits::Jacobian, panzer::Traits,ExtraT>; 

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
  #define PANZER_INSTANTIATE_TEMPLATE_CLASS_HESSIAN_THREE_T(name,ExtraT) \
    template class name<panzer::Traits::Hessian, panzer::Traits,ExtraT>;
#else
  #define PANZER_INSTANTIATE_TEMPLATE_CLASS_HESSIAN_THREE_T(name,ExtraT) 
#endif

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_THREE_T(name,ExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_THREE_T(name,ExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_TANGENT_THREE_T(name,ExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_THREE_T(name,ExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_HESSIAN_THREE_T(name,ExtraT)

// THREE (two user defined) template arguments
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_THREE_2U_T(name,FirstExtraT,SecondExtraT) \
  template class name<panzer::Traits::Residual,FirstExtraT,SecondExtraT>; 

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_TANGENT_THREE_2U_T(name,FirstExtraT,SecondExtraT) \
  template class name<panzer::Traits::Tangent,FirstExtraT,SecondExtraT>; 

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_THREE_2U_T(name,FirstExtraT,SecondExtraT) \
  template class name<panzer::Traits::Jacobian,FirstExtraT,SecondExtraT>; 

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
  #define PANZER_INSTANTIATE_TEMPLATE_CLASS_HESSIAN_THREE_2U_T(name,FirstExtraT,SecondExtraT) \
    template class name<panzer::Traits::Hessian,FirstExtraT,SecondExtraT>;
#else
  #define PANZER_INSTANTIATE_TEMPLATE_CLASS_HESSIAN_THREE_2U_T(name,FirstExtraT,SecondExtraT)
#endif

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_THREE_2U_T(name,FirstExtraT,SecondExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_THREE_2U_T(name,FirstExtraT,SecondExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_TANGENT_THREE_2U_T(name,FirstExtraT,SecondExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_THREE_2U_T(name,FirstExtraT,SecondExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_HESSIAN_THREE_2U_T(name,FirstExtraT,SecondExtraT)

// FOUR (two user defined) template arguments
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_FOUR_T(name,FirstExtraT,SecondExtraT) \
  template class name<panzer::Traits::Residual, panzer::Traits,FirstExtraT,SecondExtraT>; 

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_TANGENT_FOUR_T(name,FirstExtraT,SecondExtraT) \
  template class name<panzer::Traits::Tangent, panzer::Traits,FirstExtraT,SecondExtraT>; 

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_FOUR_T(name,FirstExtraT,SecondExtraT) \
  template class name<panzer::Traits::Jacobian, panzer::Traits,FirstExtraT,SecondExtraT>; 

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
  #define PANZER_INSTANTIATE_TEMPLATE_CLASS_HESSIAN_FOUR_T(name,FirstExtraT,SecondExtraT) \
    template class name<panzer::Traits::Hessian, panzer::Traits,FirstExtraT,SecondExtraT>;
#else 
  #define PANZER_INSTANTIATE_TEMPLATE_CLASS_HESSIAN_FOUR_T(name,FirstExtraT,SecondExtraT)
#endif

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_FOUR_T(name,FirstExtraT,SecondExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_FOUR_T(name,FirstExtraT,SecondExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_TANGENT_FOUR_T(name,FirstExtraT,SecondExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_FOUR_T(name,FirstExtraT,SecondExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_HESSIAN_FOUR_T(name,FirstExtraT,SecondExtraT)

#endif
