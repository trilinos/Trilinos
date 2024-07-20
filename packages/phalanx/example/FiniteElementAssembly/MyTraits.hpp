// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_MY_TRAITS_HPP
#define PHX_MY_TRAITS_HPP

// mpl (Meta Programming Library) templates
#include "Sacado_mpl_vector.hpp"
#include "Sacado_mpl_find.hpp"

// Include User Data Types
#include "Phalanx_config.hpp"
#include "Phalanx_Traits.hpp"
#include "Sacado.hpp" // FAD type
#include "Workset.hpp"

namespace PHX {

  struct MyTraits {
    
    // ******************************************************************
    // *** Scalar Types
    // ******************************************************************
    
    // Scalar types we plan to use
    typedef double RealType;
    typedef Sacado::Fad::DFad<double> FadType;
    typedef Sacado::Fad::SFad<double,1> JvFadType;
    
    // ******************************************************************
    // *** Evaluation Types
    // ******************************************************************
    struct Residual { typedef RealType ScalarT; };
    struct Jacobian { typedef FadType ScalarT;  };
    struct Jv { typedef JvFadType ScalarT;  };
    //typedef Sacado::mpl::vector<Residual, Jacobian, Jv> EvalTypes;
    typedef Sacado::mpl::vector<Residual, Jacobian> EvalTypes;

    // ******************************************************************
    // *** User Defined Object Passed in for Evaluation Method
    // ******************************************************************
    typedef void* SetupData;
    typedef const Workset& EvalData;
    typedef void* PreEvalData;
    typedef void* PostEvalData;

  };

  // Data Types for each evaluation type
  template<>
  struct eval_scalar_types<PHX::MyTraits::Residual> 
  { typedef Sacado::mpl::vector<PHX::MyTraits::RealType> type; };

  template<>
  struct eval_scalar_types<PHX::MyTraits::Jacobian> 
  { typedef Sacado::mpl::vector<PHX::MyTraits::FadType> type; };

  template<>
  struct eval_scalar_types<PHX::MyTraits::Jv> 
  { typedef Sacado::mpl::vector<PHX::MyTraits::JvFadType> type; };
}

#endif
