// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
//#include "Panzer_GlobalEvaluationDataContainer.hpp"

// Debugging information
//#include "Phalanx_Print.hpp"

// forward declaration
namespace Intrepid2 {
class Orientation;
}

namespace panzer {

  class GlobalEvaluationDataContainer;
  
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
      Teuchos::RCP<const std::vector<panzer::Workset>> worksets_;
      Teuchos::RCP<const std::vector<Intrepid2::Orientation>> orientations_;
    };
    using SetupData = const SD&;

    using EvalData = const panzer::Workset&;

    struct PED {
      PED();
      Teuchos::RCP<GlobalEvaluationDataContainer> gedc;
      std::string first_sensitivities_name;
      std::string second_sensitivities_name;
    };
    using PreEvalData = const PED&;

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
  { typedef Sacado::mpl::vector<panzer::Traits::FadType,panzer::Traits::RealType,bool> type; };

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
  template<>
  struct eval_scalar_types<panzer::Traits::Hessian> 
  { typedef Sacado::mpl::vector<panzer::Traits::HessianType,bool> type; };
#endif

}

#endif
