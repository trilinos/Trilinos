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

#include "Phalanx_Traits.hpp"

// Include User Data Types
#include "Phalanx_config.hpp" // for std::vector
#include "CellData.hpp"
#include "Sacado.hpp"

namespace PHX {

  /*! \brief Traits class for testing.
    
  */
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
    typedef PHX::mpl::type_list<Residual, Jacobian> EvalTypes;

    // ******************************************************************
    // *** User Defined Object Passed in for Evaluation Method
    // ******************************************************************
    typedef void* SetupData;
    typedef std::vector<CellData>& EvalData;
    typedef void* PreEvalData;
    typedef void* PostEvalData;

  };


}

namespace PHX {

  template<>
  struct eval_scalar_types<PHX::MyTraits::Residual> 
  { typedef PHX::mpl::type_list<PHX::MyTraits::RealType> type; };

  template<>
  struct eval_scalar_types<PHX::MyTraits::Jacobian> 
  { typedef PHX::mpl::type_list<PHX::MyTraits::RealType,PHX::MyTraits::FadType> type; };

}

#endif
