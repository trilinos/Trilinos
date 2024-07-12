// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHALANX_PRINT_HPP
#define PHALANX_PRINT_HPP

#include "Teuchos_TypeNameTraits.hpp"
#include <string>

namespace PHX {

  /// Default print function for a user defined object, typically used for user defined MDField extents. Users can specialize this.
  template<typename T> std::string print()
  {return Teuchos::demangleName(typeid(T).name());}

}

#endif
