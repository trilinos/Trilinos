// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHALANX_MDFIELD_EXTENT_TRAITS_HPP
#define PHALANX_MDFIELD_EXTENT_TRAITS_HPP

#include "Phalanx_Print.hpp"
#include <type_traits>

namespace PHX {
  /// Identifies that a user defined struct is a dimension template parameter for MDFields. Users must specialize for true_types.
  template<typename Extent> struct is_extent : std::false_type {};
}

/** \brief  Macro for implementing the body of a simple ArrayDimTag
 *  \param ADT  name of the tag.
 */
#define PHX_IS_EXTENT( ADT ) \
  namespace PHX { \
    template<> struct is_extent<ADT> : std::true_type {}; \
  }

/** \brief  Macro for implementing the body of a simple ArrayDimTag
 *  \param ADT  name of the tag.
 */
#define PHX_EXTENT( ADT ) \
  struct ADT {}; \
  PHX_IS_EXTENT(ADT)

#endif
