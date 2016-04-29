// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef PANZER_TRAITS_HPP
#define PANZER_TRAITS_HPP

#include "PanzerDiscFE_config.hpp"

// Teuchos includes
#include "Teuchos_RCP.hpp"

// mpl (Meta Programming Library) templates
#include "Sacado_mpl_vector.hpp"
#include "Sacado_mpl_find.hpp"

// Scalar types
#include "Sacado.hpp"
//#include "Sacado_CacheFad_DFad.hpp"
//#include "Sacado_ELRFad_DFad.hpp"
//#include "Sacado_ELRCacheFad_DFad.hpp"

#include "Phalanx_Traits.hpp"

// Include User Data Types
//#include "Phalanx_Allocator_Contiguous.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_GlobalEvaluationDataContainer.hpp"

// Debugging information
//#include "Phalanx_TypeStrings.hpp"

namespace panzer {
  
  class LinearObjContainer;

  struct Traits {

    // ******************************************************************
    // *** Scalar Types
    // ******************************************************************
    
    // Scalar types we plan to use
    typedef double RealType;
    // typedef Sacado::Fad::DFad<double> FadType;
    // typedef Sacado::CacheFad::DFad<double> FadType;
    // typedef Sacado::ELRFad::DFad<double> FadType;
    // typedef Sacado::ELRCacheFad::DFad<double> FadType;
    // typedef Sacado::Fad::SLFad<double,8> FadType;
    typedef PANZER_FADTYPE FadType;

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
    // typedef Sacado::Fad::SFad<FadType,1> HessianType;
    typedef Sacado::Fad::DFad<Sacado::Fad::SFad<RealType,1> > HessianType;
#endif
    
    // ******************************************************************
    // *** Evaluation Types
    // ******************************************************************
    struct Residual { typedef RealType ScalarT; };
    struct Jacobian { typedef FadType ScalarT;  };
    struct Tangent { typedef FadType ScalarT;  };

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
    struct Hessian { typedef HessianType ScalarT;  };
#endif

    typedef Sacado::mpl::vector< Residual
                               , Jacobian 
                               , Tangent
#ifdef Panzer_BUILD_HESSIAN_SUPPORT
                               , Hessian
#endif
                                > EvalTypes;

    // ******************************************************************
    // *** User Defined Object Passed in for Evaluation Method
    // ******************************************************************

    struct SD { 
      Teuchos::RCP< std::vector<panzer::Workset> > worksets_;
    };
    typedef SD SetupData;

    typedef panzer::Workset& EvalData;

    typedef struct {
      GlobalEvaluationDataContainer gedc;
      std::string sensitivities_name;
    } PreEvalData;
    // typedef GlobalEvaluationDataContainer& PreEvalData;

    typedef void* PostEvalData;

  };
 
}

namespace PHX {

  template<>
  struct eval_scalar_types<panzer::Traits::Residual> 
  { typedef Sacado::mpl::vector<panzer::Traits::RealType,bool> type; };

  template<>
  struct eval_scalar_types<panzer::Traits::Jacobian> 
  { typedef Sacado::mpl::vector<panzer::Traits::FadType,panzer::Traits::RealType,bool> type; };

  template<>
  struct eval_scalar_types<panzer::Traits::Tangent> 
  { typedef Sacado::mpl::vector<panzer::Traits::FadType,bool> type; };

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
  template<>
  struct eval_scalar_types<panzer::Traits::Hessian> 
  { typedef Sacado::mpl::vector<panzer::Traits::HessianType,bool> type; };
#endif

}

#endif
