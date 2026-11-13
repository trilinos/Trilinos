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

#include <Kokkos_Core_fwd.hpp>
#if !defined(KOKKOS_ENABLE_IMPL_VIEW_LEGACY)
#include "Sacado.hpp"
#else
#include "Kokkos_View_Fad.hpp"
#include "Kokkos_DynRankView_Fad.hpp"
#endif
#include "Phalanx_config.hpp"
#include "Phalanx_TypeList.hpp"

namespace PHX {
  
  struct TraitsBase {};

  /** \brief Specialize for each evaluation type, providing a nested
      \c type that is a PHX::mpl::type_list of the scalar types fields
      of that evaluation type may use. */
  template <typename EvalT> struct eval_scalar_types;

  // template <typename T> 
  // struct eval_scalar_types
  // { PHALANX_ERROR_MissingTraitsSpecializationFor_eval_scalar_types(); };

}
   
#endif
