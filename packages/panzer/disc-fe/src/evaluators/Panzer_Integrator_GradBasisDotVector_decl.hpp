// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   __Panzer_Integrator_GradBasisDotVector_decl_hpp__
#define   __Panzer_Integrator_GradBasisDotVector_decl_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// C++
#include <string>

// Panzer
#include "Panzer_EvaluatorStyle.hpp"
#include "Panzer_Evaluator_WithBaseImpl.hpp"

// Phalanx
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"

namespace panzer
{
  /**
   *  \brief Computes \f$ Ma(x)b(x)\cdots\int\vec{s}(x)\cdot\nabla\phi(x)\,dx
   *         \f$.
   *
   *  Evaluates the integral
   *  \f[
   *    Ma(x)b(x)\cdots\int\vec{s}(x)\cdot\nabla\phi(x)\,dx,
   *  \f]
   *  where \f$ M \f$ is some constant, \f$ a(x) \f$, \f$ b(x) \f$, etc., are
   *  some fields that depend on position, \f$ \vec{s} \f$ is some
   *  vector-valued function, and \f$ \phi \f$ is some basis.
   */
  template<typename EvalT, typename Traits>
  class Integrator_GradBasisDotVector
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
       *    Ma(x)b(x)\cdots\int\vec{s}(x)\cdot\nabla\phi(x)\,dx,
       *  \f]
       *  where \f$ M \f$ is some constant, \f$ a(x) \f$, \f$ b(x) \f$, etc.,
       *  are some fields that depend on position, \f$ \vec{s} \f$ is some
       *  vector-valued function, and \f$ \phi \f$ is some basis.
       *
       *  \param[in] evalStyle  An `enum` declaring the behavior of this
       *                        `Evaluator`, which is to either:
       *                        - compute and contribute (`CONTRIBUTES`), or
       *                        - compute and store (`EVALUATES`).
       *  \param[in] resName    The name of either the contributed or evaluated
       *                        field, depending on `evalStyle`.
       *  \param[in] fluxName   The name of the vector-valued function being
       *                        integrated (\f$ \vec{s} \f$).
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
       *  \param[in] vecDL      The vector data layout that you'd like to use.
       *                        If not specified, this defaults to
       *                        `Teuchos::null` and the vector data layout from
       *                        the given `ir` is used.
       *
       *  \throws std::invalid_argument If any of the inputs are invalid.
       *  \throws std::logic_error      If the `basis` supplied does not
       *                                support the gradient operator, or if
       *                                the dimension of the space exceeds the
       *                                dimension of the vector data layout.
       */
      Integrator_GradBasisDotVector(
        const panzer::EvaluatorStyle&        evalStyle,
        const std::string&                   resName,
        const std::string&                   fluxName,
        const panzer::BasisIRLayout&         basis,
        const panzer::IntegrationRule&       ir,
        const double&                        multiplier = 1,
        const std::vector<std::string>&      fmNames    =
          std::vector<std::string>(),
        const Teuchos::RCP<PHX::DataLayout>& vecDL      = Teuchos::null);

      /**
       *  \brief `ParameterList` Constructor.
       *
       *  Creates an `Evaluator` to evaluate the integral
       *  \f[
       *    Ma(x)b(x)\cdots\int\vec{s}(x)\cdot\nabla\phi(x)\,dx,
       *  \f]
       *  where \f$ M \f$ is some constant, \f$ a(x) \f$, \f$ b(x) \f$, etc.,
       *  are some fields that depend on position, \f$ \vec{s} \f$ is some
       *  vector-valued function, and \f$ \phi \f$ is some basis.
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
                         <Parameter name = "Residual Name"      type = "std::string"                         value = (required)    />
                         <Parameter name = "Flux Name"          type = "std::string"                         value = (required)    />
                         <Parameter name = "Basis"              type = "RCP<panzer::BasisIRLayout>"          value = (required)    />
                         <Parameter name = "IR"                 type = "RCP<panzer::IntegrationRule>"        value = (required)    />
                         <Parameter name = "Multiplier"         type = "double"                              value = (required)    />
                         <Parameter name = "Field Multipliers"  type = "RCP<const std::vector<std::string>>" value = null (default)/>
                         <Parameter name = "Vector Data Layout" type = "RCP<PHX::DataLayout>"                value = null (default)/>
                       </ParameterList>
                       \endcode
       *               where
       *               - "Residual Name" is the name for the term this
       *                 `Evaluator` is evaluating,
       *               - "Flux Name" is the name of the vector-valued function
       *                 being integrated (\f$ \vec{s} \f$),
       *               - "Basis" is the basis that you'd like to use (\f$ \phi
       *                 \f$),
       *               - "IR" is the integration rule that you'd like to use,
       *               - "Multiplier" is the scalar multiplier out in front of
       *                 the integral you're computing (\f$ M \f$), and
       *               - "Field Multipliers" is an optional list of names of
       *                 fields that are multipliers out in front of the
       *                 integral you're computing (\f$ a(x) \f$, \f$ b(x) \f$,
       *                 etc.).
       *               - "Vector Data Layout" is the vector data layout that
       *                 you'd like to use.
       */
      Integrator_GradBasisDotVector(
        const Teuchos::ParameterList& p);

      /**
       *  \brief Post-Registration Setup.
       *
       *  Sets the basis index, and gets the PHX::View versions of the field
       *  multiplier, if there are any.
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
       *  \brief This empty struct allows us to optimize
       *         `operator()()` depending on the number of field
       *         multipliers. This is the version that does not use
       *         shared memory.
       */
      template<int NUM_FIELD_MULT>
      struct FieldMultTag {};

      /**
       *  \brief This empty struct allows us to optimize
       *         `operator()()` depending on the number of field
       *         multipliers. This is the shared memory version.
       */
      template<int NUM_FIELD_MULT>
      struct SharedFieldMultTag {};

      /**
       *  \brief Perform the integration.
       *
       *  Generally speaking, for a given cell in the `Workset`, this routine
       *  loops over quadrature points, vector dimensions, and bases to perform
       *  the integration, scaling the vector field to be integrated by the
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
        const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) const;

      /**
       *  \brief Perform the integration.
       *
       *  Generally speaking, for a given cell in the `Workset`, this routine
       *  loops over quadrature points, vector dimensions, and bases to perform
       *  the integration, scaling the vector field to be integrated by the
       *  multiplier (\f$ M \f$) and any field multipliers (\f$ a(x) \f$, \f$
       *  b(x) \f$, etc.). Uses Shared memory.
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
        const SharedFieldMultTag<NUM_FIELD_MULT>& tag,
        const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) const;

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

      /// Type for shared memory
      using scratch_view = Kokkos::View<ScalarT* ,typename PHX::DevLayout<ScalarT>::type,typename PHX::exec_space::scratch_memory_space,Kokkos::MemoryUnmanaged>;

      /**
       *  \brief An `enum` determining the behavior of this `Evaluator`.
       *
       *  This `Evaluator` will compute the result of its integration and then:
       *  - CONTRIBUTES:                contribute it to a specified residual,
       *                                not saving anything; or
       *  - EVALUATES:                  save it under a specified name for
       *                                future use.
       */
      const panzer::EvaluatorStyle evalStyle_;

      /**
       *  \brief A field to which we'll contribute, or in which we'll store,
       *         the result of computing this integral.
       */
      PHX::MDField<ScalarT, panzer::Cell, panzer::BASIS> field_;

      /**
       *  \brief A field representing the vector-valued function we're
       *         integrating (\f$ \vec{s} \f$).
       */
      PHX::MDField<const ScalarT, panzer::Cell, panzer::IP, panzer::Dim>
      vector_;

      /**
       *  \brief The scalar multiplier out in front of the integral (\f$ M
       *         \f$).
       */
      double multiplier_;

      /**
       *  \brief The (possibly empty) list of fields that are multipliers out
       *         in front of the integral (\f$ a(x) \f$, \f$ b(x) \f$, etc.).
       */
      std::vector<PHX::MDField<const ScalarT, panzer::Cell, panzer::IP>>
      fieldMults_;

      /**
       *  \brief The `PHX::View` representation of the (possibly empty) list
       *         of fields that are multipliers out in front of the integral
       *         (\f$ a(x) \f$, \f$ b(x) \f$, etc.).
       */
    PHX::View<PHX::UnmanagedView<const ScalarT**>* > kokkosFieldMults_;

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

    /// Temporary used when shared memory is disabled
    PHX::View<ScalarT*> tmp_;

  }; // end of class Integrator_GradBasisDotVector

} // end of namespace panzer

#endif // __Panzer_Integrator_GradBasisDotVector_decl_hpp__
