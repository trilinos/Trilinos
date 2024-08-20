// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_DFAD_TRAITS_HPP
#define PHX_DFAD_TRAITS_HPP

#include "Sacado_mpl_vector.hpp"
#include "Sacado_mpl_find.hpp"
#include "Phalanx_Traits.hpp"

// Include User Data Types
#include "Phalanx_config.hpp"
#include "Sacado.hpp"
#include "Dimension.hpp"
#include "Cell.hpp"
#include "Workset.hpp"

namespace PHX {

  struct MyTraits {
    
    // ******************************************************************
    // *** Scalar Types
    // ******************************************************************
    
    // Scalar types we plan to use
    typedef double RealType;
    typedef Sacado::Fad::DFad<double> FadType;
    
    // ******************************************************************
    // *** Evaluation Types
    // ******************************************************************
    struct Residual { typedef RealType ScalarT; };
    struct Jacobian { typedef FadType ScalarT;  };
    typedef Sacado::mpl::vector<Residual, Jacobian> EvalTypes;

    // ******************************************************************
    // *** User Defined Object Passed in for Evaluation Method
    // ******************************************************************
    typedef void* SetupData;
    typedef const MyWorkset& EvalData;
    typedef void* PreEvalData;
    typedef void* PostEvalData;

  };

  template<>
  struct eval_scalar_types<PHX::MyTraits::Residual> 
  { typedef Sacado::mpl::vector<PHX::MyTraits::RealType> type; };

  template<>
  struct eval_scalar_types<PHX::MyTraits::Jacobian> 
  { typedef Sacado::mpl::vector<PHX::MyTraits::RealType,PHX::MyTraits::FadType> type; };

}

#endif
