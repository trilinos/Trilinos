// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   PANZER_EVALUATOR_GRADBASISTIMESSCALAR_DECL_HPP
#define   PANZER_EVALUATOR_GRADBASISTIMESSCALAR_DECL_HPP

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// C++
#include <string>

// Kokkos
#include "Kokkos_DynRankView.hpp"

// Panzer
#include "Panzer_EvaluatorStyle.hpp"
#include "Panzer_Evaluator_WithBaseImpl.hpp"

// Phalanx
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"

/**
 *  \brief Computes \f$ Ma(x)b(x)\cdots\int s(x)\nabla\phi(x)\,dx \f$.
 *
 *  Evaluates the integral
 *  \f[
 *    Ma(x)b(x)\cdots\int s(x)\nabla\phi(x)\,dx,
 *  \f]
 *  where \f$ M \f$ is some constant, \f$ a(x) \f$, \f$ b(x) \f$, etc., are
 *  some fields that depend on position, \f$ s \f$ is some scalar function,
 *  and \f$ \phi \f$ is some vector basis.
 */
namespace panzer
{
  template<typename EvalT, typename Traits>
  class Integrator_GradBasisTimesScalar
    :
    public panzer::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<EvalT, Traits>
  {
    public:

      /**
       *  \brief Main Constructor.
       *
       *  Creates an `Evaluator` to evaluate the integral
       *  \f[
       *    Ma(x)b(x)\cdots\int s(x)\nabla\phi(x)\,dx,
       *  \f]
       *  where \f$ M \f$ is some constant, \f$ a(x) \f$, \f$ b(x) \f$, etc.,
       *  are some fields that depend on position, \f$ s \f$ is some scalar
       *  function, and \f$ \phi \f$ is some vector basis.
       *
       *  \param[in] evalStyle  An `enum` declaring the behavior of this
       *                        `Evaluator`, which is to either:
       *                        - compute and contribute (`CONTRIBUTES`), or
       *                        - compute and store (`EVALUATES`).
       *  \param[in] resNames   The names of either the contributed or
       *                        evaluated fields, depending on `evalStyle`.
       *  \param[in] scalar     The name of the scalar function being
       *                        integrated (\f$ s \f$).
       *  \param[in] basis      The basis that you'd like to use (\f$ \phi
       *                        \f$).
       *  \param[in] ir         The integration rule that you'd like to use.
       *  \param[in] multiplier The scalar multiplier out in front of the
       *                        integral you're computing (\f$ M \f$).  If not
       *                        specified, this defaults to 1.
       *  \param[in] fmNames    A list of names of fields that are multipliers
       *                        out in front of the integral you're computing
       *                        (\f$ a(x) \f$, \f$ b(x) \f$, etc.).  If not
       *                        specified, this defaults to an empty `vector`.
       *
       *  \throws std::invalid_argument If any of the inputs are invalid.
       *  \throws std::logic_error      If the `basis` supplied does not
       *                                support the gradient operator, or if
       *                                the dimension of the space doesn't
       *                                match the size of `resNames`.
       */
      Integrator_GradBasisTimesScalar(
        const panzer::EvaluatorStyle&   evalStyle,
        const std::vector<std::string>& resNames,
        const std::string&              scalar,
        const panzer::BasisIRLayout&    basis,
        const panzer::IntegrationRule&  ir,
        const double&                   multiplier = 1,
        const std::vector<std::string>& fmNames    =
          std::vector<std::string>());

      /**
       *  \brief `ParameterList` Constructor.
       *
       *  Creates an `Evaluator` to evaluate the integral
       *  \f[
       *    Ma(x)b(x)\cdots\int s(x)\nabla\phi(x)\,dx,
       *  \f]
       *  where \f$ M \f$ is some constant, \f$ a(x) \f$, \f$ b(x) \f$, etc.,
       *  are some fields that depend on position, \f$ s \f$ is some scalar
       *  function, and \f$ \phi \f$ is some vector basis.
       *
       *  \note This constructor exists to preserve the older way of creating
       *        an `Evaluator` with a `ParameterList`; however, it is
       *        _strongly_ advised that you _not_ use this `ParameterList`
       *        Constructor, but rather that you favor the main constructor
       *        with its compile-time argument checking instead.
       *
       *  \param[in] p A `ParameterList` of the form
                       \code{.xml}
                       <ParameterList>
                         <Parameter name = "Residual Names"     type = "std::vector<std::string>"            value = (required)    />
                         <Parameter name = "Scalar Name"        type = "std::string"                         value = (required)    />
                         <Parameter name = "Basis"              type = "RCP<panzer::BasisIRLayout>"          value = (required)    />
                         <Parameter name = "IR"                 type = "RCP<panzer::IntegrationRule>"        value = (required)    />
                         <Parameter name = "Multiplier"         type = "double"                              value = (required)    />
                         <Parameter name = "Field Multipliers"  type = "RCP<const std::vector<std::string>>" value = null (default)/>
                       </ParameterList>
                       \endcode
       *               where
       *               - "Residual Names" are the names for the terms this
       *                 `Evaluator` is evaluating,
       *               - "Scalar Name" is the name of the scalar function being
       *                 integrated (\f$ s \f$),
       *               - "Basis" is the basis that you'd like to use (\f$ \phi
       *                 \f$),
       *               - "IR" is the integration rule that you'd like to use,
       *               - "Multiplier" is the scalar multiplier out in front of
       *                 the integral you're computing (\f$ M \f$), and
       *               - "Field Multipliers" is an optional list of names of
       *                 fields that are multipliers out in front of the
       *                 integral you're computing (\f$ a(x) \f$, \f$ b(x) \f$,
       *                 etc.).
       */
      Integrator_GradBasisTimesScalar(
        const Teuchos::ParameterList& p);

      /**
       *  \brief Post-Registration Setup.
       *
       *  Sets the basis index, and gets the PHX::View versions of the field
       *  multipliers, if there are any.
       *
       *  \param[in] sd Essentially a list of `Workset`s, which are collections
       *                of cells (elements) that all live on a single process.
       *  \param[in] fm This is an unused part of the `Evaluator` interface.
       */
      void
      postRegistrationSetup(
        typename Traits::SetupData d,
        PHX::FieldManager<Traits>& fm);

      /**
       *  \brief Evaluate Fields.
       *
       *  This actually performs the integration by calling `operator()()` in a
       *  `Kokkos::parallel_for` over the cells in the `Workset`.
       *
       *  \param[in] workset The `Workset` on which you're going to do the
       *                     integration.
       */
      void
      evaluateFields(
        typename Traits::EvalData d);

      /**
       *  \brief This empty struct allows us to optimize `operator()()`
       *         depending on the number of field multipliers.
       */
      template<int NUM_FIELD_MULT>
      struct FieldMultTag
      {
      }; // end of struct FieldMultTag

      /**
       *  \brief Perform the integration.
       *
       *  Generally speaking, for a given cell in the `Workset`, this routine
       *  loops over quadrature points, bases, and dimensions to perform the
       *  integration, scaling the vector field to be integrated by the
       *  multiplier (\f$ M \f$) and any field multipliers (\f$ a(x) \f$, \f$
       *  b(x) \f$, etc.).
       *
       *  \note Optimizations are made for the cases in which we have no field
       *        multipliers or only a single one.
       *
       *  \param[in] tag  An indication of the number of field multipliers we
       *                  have; either 0, 1, or something else.
       *  \param[in] cell The cell in the `Workset` over which to integrate.
       */
      template<int NUM_FIELD_MULT>
      KOKKOS_INLINE_FUNCTION
      void
      operator()(
        const FieldMultTag<NUM_FIELD_MULT>& tag,
        const std::size_t&                  cell) const;

    private:

      /**
       *  \brief Get Valid Parameters.
       *
       *  Get all the parameters that we support such that the `ParameterList`
       *  Constructor can do some validation of the input `ParameterList`.
       *
       *  \returns A `ParameterList` with all the valid parameters (keys) in
       *           it.  The values tied to those keys are meaningless default
       *           values.
       */
      Teuchos::RCP<Teuchos::ParameterList>
      getValidParameters() const;

      /**
       *  \brief The scalar type.
       */
      using ScalarT = typename EvalT::ScalarT;

      /**
       *  \brief An `enum` determining the behavior of this `Evaluator`.
       *
       *  This `Evaluator` will compute the result of its integration and then:
       *  - CONTRIBUTES:  contribute it to a specified residual, not saving
       *                  anything; or
       *  - EVALUATES:    save it under a specified name for future use.
       */
      const panzer::EvaluatorStyle evalStyle_;

      /**
       *  \brief The fields to which we'll contribute, or in which we'll store,
       *         the result of computing this integral.
       */
      std::vector<PHX::MDField<ScalarT,Cell,BASIS>> fields_host_;
      using InnerView = PHX::UnmanagedView<ScalarT**>;
      using OuterView = PHX::View<InnerView*>;
      OuterView fields_;

      /**
       *  \brief A field representing the scalar function we're integrating
       *         (\f$ s \f$).
       */
      PHX::MDField<const ScalarT, Cell, IP> scalar_;

      /**
       *  \brief The scalar multiplier out in front of the integral (\f$ M
       *         \f$).
       */
      ScalarT multiplier_;

      /**
       *  \brief The scalar multiplier out in front of the integral (\f$ M
       *         \f$).
       */
      std::vector<PHX::MDField<const ScalarT, Cell, IP>> fieldMults_;

      /**
       *  \brief The `PHX::View` representation of the (possibly empty) list
       *         of fields that are multipliers out in front of the integral
       *         (\f$ a(x) \f$, \f$ b(x) \f$, etc.).
       */
      PHX::View<PHX::UnmanagedView<const ScalarT**>*> kokkosFieldMults_;

      /**
       *  \brief The number of dimensions associated with the gradient.
       */
      int numDims_;

      /**
       *  \brief The name of the basis we're using.
       */
      std::string basisName_;

      /**
       *  \brief The index in the `Workset` bases for our particular
       *         `BasisIRLayout` name.
       */
      std::size_t basisIndex_;

      /**
       *  \brief The gradient vector basis information necessary for
       *         integration.
       */
      PHX::MDField<double, panzer::Cell, panzer::BASIS, panzer::IP,
        panzer::Dim> basis_;

  }; // end of class Integrator_GradBasisTimesScalar

} // end of namespace panzer

#endif // PANZER_EVALUATOR_GRADBASISTIMESSCALAR_DECL_HPP
