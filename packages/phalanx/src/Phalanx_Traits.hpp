// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_TRAITS_HPP
#define PHX_TRAITS_HPP

#include "Kokkos_View_Fad.hpp"
#include "Kokkos_DynRankView_Fad.hpp"
#include "Phalanx_config.hpp"

namespace PHX {
  
  struct TraitsBase {};

  template <typename mpl_vector> struct eval_scalar_types;

  // template <typename T> 
  // struct eval_scalar_types
  // { PHALANX_ERROR_MissingTraitsSpecializationFor_eval_scalar_types(); };

}
   
#endif
