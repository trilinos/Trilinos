// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_TESTING_MY_TRAITS_HPP
#define PHX_TESTING_MY_TRAITS_HPP

#include "Phalanx_config.hpp" // for std::vector
#include "Phalanx_Traits.hpp"
#include "Sacado_mpl_vector.hpp"

namespace PHX {

  /*! \brief Traits class for testing. */
  struct MyTraits {

    // ******************************************************************
    // *** Evaluation Types
    // ******************************************************************
    struct Residual { typedef double ScalarT; };
    struct Jacobian { typedef double ScalarT; };
    typedef Sacado::mpl::vector<Residual,Jacobian> EvalTypes;

    // ******************************************************************
    // *** User Defined Object Passed in for Evaluation Method
    // ******************************************************************
    typedef int SetupData;
    typedef int EvalData;
    typedef int PreEvalData;
    typedef int PostEvalData;

  };

}

namespace PHX {
  template<>
  struct eval_scalar_types<PHX::MyTraits::Residual> 
  { typedef Sacado::mpl::vector<double> type; };
}

namespace PHX {
  template<>
  struct eval_scalar_types<PHX::MyTraits::Jacobian> 
  { typedef Sacado::mpl::vector<double> type; };
}

#endif
