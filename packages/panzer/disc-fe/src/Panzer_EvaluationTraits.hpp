// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATION_TRAITS_HPP
#define PANZER_EVALUATION_TRAITS_HPP

namespace panzer {
  
  struct EvaluationTraits {
    template <class T> struct apply {
      typedef typename T::ScalarT type;
    };
  };

}

#endif
