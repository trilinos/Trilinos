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

  /** \brief Traits class used by Sacado::ScalarParameterLibrary (see
    * panzer::ParamLib) to map an evaluation type to its scalar type.
    */
  struct EvaluationTraits {
    /// \brief Maps evaluation type T to its scalar type, as T::ScalarT.
    template <class T> struct apply {
      typedef typename T::ScalarT type;
    };
  };

}

#endif
