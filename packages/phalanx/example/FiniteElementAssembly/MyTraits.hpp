// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
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
