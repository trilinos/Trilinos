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

  /**
   * \struct Traits
   *
   * \brief Panzer's specialization of the Phalanx traits class.
   *
   * This struct binds together all of the type information that the
   * Phalanx `DAG Manager` and `EvaluationContainer` machinery needs in
   * order to build and run field evaluators for panzer: the scalar
   * types used for each evaluation type (residual, Jacobian, etc.), the
   * `panzer::Traits::EvalTypes` mpl::vector enumerating those evaluation
   * types, and the four user-defined data types (`SetupData`, `EvalData`,
   * `PreEvalData`, `PostEvalData`) that Phalanx passes through to
   * `PHX::Evaluator::postRegistrationSetup()`, `evaluateFields()`,
   * `preEvaluate()`, and `postEvaluate()` respectively.
   */
  struct Traits {

    // ******************************************************************
    // *** Scalar Types
    // ******************************************************************

    /// Scalar type used for the Residual evaluation type (plain double).
    typedef double RealType;
    // typedef Sacado::Fad::DFad<double> FadType;
    // typedef Sacado::CacheFad::DFad<double> FadType;
    // typedef Sacado::ELRFad::DFad<double> FadType;
    // typedef Sacado::ELRCacheFad::DFad<double> FadType;
    // typedef Sacado::Fad::SLFad<double,8> FadType;
    /// Sacado forward-mode AD scalar type used for the Jacobian and Tangent evaluation types.
    typedef PANZER_FADTYPE FadType;

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
    // typedef Sacado::Fad::SFad<FadType,1> HessianType;
    /// Nested Sacado AD scalar type (FAD-over-FAD) used for the Hessian evaluation type. Only defined when built with Hessian support.
    typedef Sacado::Fad::DFad<Sacado::Fad::SFad<RealType,1> > HessianType;
#endif

    // ******************************************************************
    // *** Evaluation Types
    // ******************************************************************

    /// Evaluation type for computing the residual only, using #RealType as the scalar type.
    struct Residual { typedef RealType ScalarT; };
    /// Evaluation type for computing the residual and its Jacobian, using #FadType as the scalar type.
    struct Jacobian { typedef FadType ScalarT;  };
    /// Evaluation type for computing directional derivatives (e.g. parameter sensitivities), using #FadType as the scalar type.
    struct Tangent { typedef FadType ScalarT;  };

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
    /// Evaluation type for computing second derivatives, using #HessianType as the scalar type. Only defined when built with Hessian support.
    struct Hessian { typedef HessianType ScalarT;  };
#endif

    /// Sacado::mpl::vector enumerating all evaluation types supported by this build (Residual, Jacobian, Tangent, and optionally Hessian).
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

    /**
     * \struct SD
     * \brief User-defined data passed to `PHX::Evaluator::postRegistrationSetup()` during the one-time setup phase.
     *
     * Gives evaluators access to the full set of worksets and cell
     * orientations up front, so any data needed at evaluation time that
     * is expensive to look up (e.g. field pointers) can be cached during
     * setup rather than on every call to `evaluateFields()`.
     */
    struct SD {
      /// All worksets that will be evaluated over the lifetime of this evaluator.
      Teuchos::RCP<const std::vector<panzer::Workset>> worksets_;
      /// Cell orientations corresponding to the worksets above.
      Teuchos::RCP<const std::vector<Intrepid2::Orientation>> orientations_;
    };
    /// Data type passed to `PHX::Evaluator::postRegistrationSetup()`. \sa SD
    using SetupData = const SD&;

    /// Data type passed to `PHX::Evaluator::evaluateFields()`: a single workset holding the cell data for one evaluation call.
    using EvalData = const panzer::Workset&;

    /**
     * \struct PED
     * \brief User-defined data passed to `PHX::Evaluator::preEvaluate()`, called once before each residual/Jacobian fill.
     */
    struct PED {
      /// Default constructor. Allocates #gedc.
      PED();
      /// Container of global (non-cell-local) data, e.g. distributed vectors needed by gather/scatter evaluators.
      Teuchos::RCP<GlobalEvaluationDataContainer> gedc;
      /// Name under which first-order sensitivity data is registered in #gedc, if any.
      std::string first_sensitivities_name;
      /// Name under which second-order sensitivity data is registered in #gedc, if any.
      std::string second_sensitivities_name;
    };
    /// Data type passed to `PHX::Evaluator::preEvaluate()`. \sa PED
    using PreEvalData = const PED&;

    /// Data type passed to `PHX::Evaluator::postEvaluate()`, called once after each residual/Jacobian fill. Unused by panzer (always null).
    typedef void* PostEvalData;

  };

}

namespace PHX {

  /**
   * \brief Declares which scalar types `PHX::MDField`s may hold for the Residual evaluation type.
   *
   * Phalanx uses this specialization to determine, for a given
   * evaluation type, which concrete scalar types are legal for fields
   * (in addition to `Traits::Residual::ScalarT` itself) — here plain
   * `bool` fields alongside #panzer::Traits::RealType.
   */
  template<>
  struct eval_scalar_types<panzer::Traits::Residual>
  { typedef Sacado::mpl::vector<panzer::Traits::RealType,bool> type; };

  /// \brief Declares which scalar types `PHX::MDField`s may hold for the Jacobian evaluation type. \sa eval_scalar_types<panzer::Traits::Residual>
  template<>
  struct eval_scalar_types<panzer::Traits::Jacobian>
  { typedef Sacado::mpl::vector<panzer::Traits::FadType,panzer::Traits::RealType,bool> type; };

  /// \brief Declares which scalar types `PHX::MDField`s may hold for the Tangent evaluation type. \sa eval_scalar_types<panzer::Traits::Residual>
  template<>
  struct eval_scalar_types<panzer::Traits::Tangent>
  { typedef Sacado::mpl::vector<panzer::Traits::FadType,panzer::Traits::RealType,bool> type; };

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
  /// \brief Declares which scalar types `PHX::MDField`s may hold for the Hessian evaluation type. Only defined when built with Hessian support. \sa eval_scalar_types<panzer::Traits::Residual>
  template<>
  struct eval_scalar_types<panzer::Traits::Hessian>
  { typedef Sacado::mpl::vector<panzer::Traits::HessianType,bool> type; };
#endif

}

#endif
