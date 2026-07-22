// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   Panzer_Integrator_CurlBasisDotVector_hpp
#define   Panzer_Integrator_CurlBasisDotVector_hpp

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
#include "Panzer_Evaluator_WithBaseImpl.hpp"
#include "Panzer_EvaluatorStyle.hpp"

// Phalanx
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"

namespace panzer
{
  /**
   *  \brief Computes \f$
   *         Ma(x)b(x)\cdots\int\nabla\times\vec{\phi}\cdot\vec{v}\, dx \f$.
   *
   *  Evaluates the integral
   *  \f[
   *    Ma(x)b(x)\cdots\int\nabla\times\vec{\phi}\cdot\vec{v}\, dx,
   *  \f]
   *  where \f$ M \f$ is some constant, \f$ a(x) \f$, \f$ b(x) \f$, etc., are
   *  some fields that depend on position, \f$ \vec{v} \f$ is some
   *  vector-valued function, and \f$ \vec{\phi} \f$ is some vector basis.
   *
   *  \note The name can be misleading.  In contrast to 3-D, the curl of a
   *        vector in 2-D is simply a scalar.  This `Evaluator` handles both
   *        cases.
   */
  template<typename EvalT, typename Traits>
  class Integrator_CurlBasisDotVector
    :
    public panzer::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<EvalT, Traits>
  {
    public:

      /**
       *  \brief Main Constructor
       *
       *  Creates an `Evaluator` to evaluate the integral
       *  \f[
       *    Ma(x)b(x)\cdots\int\nabla\times\vec{\phi}\cdot\vec{v}\, dx,
       *  \f]
       *  where \f$ M \f$ is some constant, \f$ a(x) \f$, \f$ b(x) \f$, etc.,
       *  are some fields that depend on position, \f$ \vec{v} \f$ is some
       *  vector-valued function, and \f$ \vec{\phi} \f$ is some vector basis.
       *
       *  \param[in] evalStyle  An `enum` declaring the behavior of this
       *                        `Evaluator`, which is to either:
       *                        - compute and contribute (`CONTRIBUTES`), or
       *                        - compute and store (`EVALUATES`).
       *  \param[in] resName    The name of either the contributed or evaluated
       *                        field, depending on `evalStyle`.
       *  \param[in] valName    The name of the vector value being integrated
       *                        (\f$ \vec{s} \f$).
       *  \param[in] basis      The vector basis that you'd like to use (\f$
                                \vec{\phi} \f$).
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
       *  \throws std::logic_error      If the `basis` supplied is not a vector
       *                                basis, or if it doesn't require
       *                                orientations.
       */
      Integrator_CurlBasisDotVector(
        const panzer::EvaluatorStyle&   evalStyle,
        const std::string&              resName,
        const std::string&              valName,
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
       *    Ma(x)b(x)\cdots\int\nabla\times\vec{\phi}\cdot\vec{v}\, dx,
       *  \f]
       *  where \f$ M \f$ is some constant, \f$ a(x) \f$, \f$ b(x) \f$, etc.,
       *  are some fields that depend on position, \f$ \vec{s} \f$ is some
       *  vector-valued function, and \f$ \vec{\phi} \f$ is some vector basis.
       *
       *  \note This constructor exists to preserve the older way of creating
       *        an `Evaluator` with a `ParameterList`; however, it is
       *        _strongly_ advised that you _not_ use this `ParameterList`
       *        Constructor, but rather that you favor the Main Constructor
       *        with its compile-time argument checking instead.
       *
       *  \param[in] p A `ParameterList` of the form
                       \code{.xml}
                       <ParameterList>
                         <Parameter name = "Residual Name"     type = "std::string"                         value = (required)    />
                         <Parameter name = "Value Name"        type = "std::string"                         value = (required)    />
                         <Parameter name = "Basis"             type = "RCP<panzer::BasisIRLayout>"          value = (required)    />
                         <Parameter name = "IR"                type = "RCP<panzer::IntegrationRule>"        value = (required)    />
                         <Parameter name = "Multiplier"        type = "double"                              value = (required)    />
                         <Parameter name = "Field Multipliers" type = "RCP<const std::vector<std::string>>" value = null (default)/>
                       </ParameterList>
                       \endcode
       *               where
       *               - "Residual Name" is the name for the term this
       *                 `Evaluator` is evaluating,
       *               - "Value Name" is the name of the vector value being
       *                 integrated (\f$ \vec{s} \f$),
       *               - "Basis" is the vector basis that you'd like to use
       *                 (\f$ \vec{\phi} \f$),
       *               - "IR" is the integration rule that you'd like to use,
       *               - "Multiplier" is the scalar multiplier out in front of
       *                 the integral you're computing (\f$ M \f$), and
       *               - "Field Multipliers" is an optional list of names of
       *                 fields that are multipliers out in front of the
       *                 integral you're computing (\f$ a(x) \f$, \f$ b(x) \f$,
       *                 etc.).
       */
      Integrator_CurlBasisDotVector(
        const Teuchos::ParameterList& p);

      /**
       *  \brief `FieldTag` Constructor.
       *
       *  Creates an `Evaluator` to evaluate the integral
       *  \f[
       *    Ma(x)b(x)\cdots\int\nabla\times\vec{\phi}\cdot\vec{v}\, dx,
       *  \f]
       *  where \f$ M \f$ is some constant, \f$ a(x) \f$, \f$ b(x) \f$, etc.,
       *  are some fields that depend on position, \f$ \vec{s} \f$ is some
       *  vector-valued function, and \f$ \vec{\phi} \f$ is some vector basis.
       *
       *  \param[in] evalStyle   An `enum` declaring the behavior of this
       *                         `Evaluator`, which is to either:
       *                         - compute and contribute (`CONTRIBUTES`), or
       *                         - compute and store (`EVALUATES`).
       *  \param[in] resTag      The tag of either the contributed or evaluated
       *                         field, depending on `evalStyle`.
       *  \param[in] valTag      The tag of the vector value being integrated
       *                         (\f$ \vec{s} \f$).
       *  \param[in] bd          The vector basis descriptor that you'd like to
       *                         use (\f$ \vec{\phi} \f$).
       *  \param[in] id          The integration descriptor that you'd like to
       *                         use.
       *  \param[in] spaceDim    The spatial dimensionality of the problem.  If
       *                         not specified, this defaults to 3.
       *  \param[in] multiplier  The scalar multiplier out in front of the
       *                         integral you're computing (\f$ M \f$).  If not
       *                         specified, this defaults to 1.
       *  \param[in] multipliers A list of tags of fields that are multipliers
       *                         out in front of the integral you're computing
       *                         (\f$ a(x) \f$, \f$ b(x) \f$, etc.).  If not
       *                         specified, this defaults to an empty `vector`.
       *
       *  \throws std::invalid_argument If any of the inputs are invalid.
       *  \throws std::logic_error      If the `basis` supplied is not a vector
       *                                basis, or if it doesn't require
       *                                orientations.
       */
      Integrator_CurlBasisDotVector(
        const panzer::EvaluatorStyle&        evalStyle,
        const PHX::FieldTag&                 resTag,
        const PHX::FieldTag&                 valTag,
        const panzer::BasisDescriptor&       bd,
        const panzer::IntegrationDescriptor& id,
        const int&                           spaceDim    = 3,
        const double&                        multiplier  = 1,
        const std::vector<PHX::FieldTag>&    multipliers =
          std::vector<PHX::FieldTag>());

      /**
       *  \brief Post-Registration Setup.
       *
       *  Sets the `PHX::View`s for all of the field multipliers, sets the
       *  basis index, and sets up the field that will be used to build up the
       *  result of the integration.
       *
       *  \param[in] sd Essentially a list of `Workset`s, which are collections
       *                of cells (elements) that all live on a single process.
       *  \param[in] fm The `FieldManager` used to create the field to build up
       *                the result.
       */
      void
      postRegistrationSetup(
        typename Traits::SetupData sd,
        PHX::FieldManager<Traits>& fm);

      /**
       *  \brief Evaluate Fields.
       *
       *  This actually performs the integration using a handful of functors in
       *  `Kokkos::parallel_for`s, looping over the cells in the `Workset`.
       *
       *  \param[in] workset The `Workset` on which you're going to do the
       *                     integration.
       */
      void
      evaluateFields(
        typename Traits::EvalData d);

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
      Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

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
       *  \brief A flag indicating whether or not to use the descriptor
       *         interface.
       */
      bool useDescriptors_;

      /**
       *  \brief The `BasisDescriptor` for the basis to use.
       */
      panzer::BasisDescriptor bd_;

      /**
       *  \brief The `IntegrationDescriptor` for the quadrature to use.
       */
      panzer::IntegrationDescriptor id_;

      /**
       *  \brief A field to which we'll contribute, or in which we'll store,
       *         the result of computing this integral.
       */
      PHX::MDField<ScalarT, panzer::Cell, panzer::BASIS> field_;

      /**
       *  \brief A field representing the vector-valued function we're
       *         integrating (\f$ \vec{s} \f$) in a two-dimensional problem.
       */
      PHX::MDField<const ScalarT, panzer::Cell, panzer::IP> vector2D_;

      /**
       *  \brief A field representing the vector-valued function we're
       *         integrating (\f$ \vec{s} \f$) in a three-dimensional problem.
       */
      PHX::MDField<const ScalarT, panzer::Cell, panzer::IP, panzer::Dim>
      vector3D_;

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
    PHX::View<Kokkos::View<const ScalarT**, typename PHX::DevLayout<ScalarT>::type, Kokkos::MemoryUnmanaged>*> kokkosFieldMults_;

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
       *  \brief The spatial dimension of the vector-valued function we're
       *         integrating, either 2 or 3.
       */
      int spaceDim_;

      /**
       *  \brief A field used to build up the result of this integral when
       *         working on a two-dimensional vector field.
       */
      PHX::MDField<ScalarT, panzer::Cell, panzer::IP> result2D_;

      /**
       *  \brief A field used to build up the result of this integral when
       *         working on a three-dimensional vector field.
       */
      PHX::MDField<ScalarT, panzer::Cell, panzer::IP, panzer::Dim> result3D_;
  }; // end of class Integrator_CurlBasisDotVector

} // end of namespace panzer

#endif // Panzer_Integrator_CurlBasisDotVector_hpp
