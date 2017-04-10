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

#ifndef   __Panzer_Integrator_BasisTimesScalar_decl_hpp__
#define   __Panzer_Integrator_BasisTimesScalar_decl_hpp__

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
   *  \brief Computes \f$ Ma(x)b(x)\cdots\int s(x)\phi(x)\,dx \f$.
   *
   *  Evaluates the integral
   *  \f[
        Ma(x)b(x)\cdots\int s(x)\phi(x)\,dx,
      \f]
   *  where \f$ M \f$ is some constant, \f$ a(x) \f$, \f$ b(x) \f$, etc., are
   *  some fields that depend on position, \f$ s \f$ is some scalar function,
   *  and \f$ \phi \f$ is some basis.
   */
  template<typename EvalT, typename Traits>
  class Integrator_BasisTimesScalar
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
            Ma(x)b(x)\cdots\int s(x)\phi(x)\,dx,
          \f]
       *  where \f$ M \f$ is some constant, \f$ a(x) \f$, \f$ b(x) \f$, etc.,
       *  are some fields that depend on position, \f$ s \f$ is some scalar
       *  function, and \f$ \phi \f$ is some basis.
       *
       *  \param[in] evalStyle  An `enum` declaring the behavior of this
       *                        `Evaluator`, which is to either:
       *                        - compute and contribute (`CONTRIBUTES`), or
       *                        - compute and store (`EVALUATES`).
       *  \param[in] resName    The name of either the contributed or evaluated
       *                        field, depending on `evalStyle`.
       *  \param[in] valName    The name of the scalar value being integrated
       *                        (\f$ s \f$).
       *  \param[in] basis      The basis that you'd like to use (\f$ \phi
                                \f$).
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
       *  \throws std::logic_error      If the `basis` supplied is not a scalar
       *                                basis.
       */
      Integrator_BasisTimesScalar(
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
            Ma(x)b(x)\cdots\int s(x)\phi(x)\,dx,
          \f]
       *  where \f$ M \f$ is some constant, \f$ a(x) \f$, \f$ b(x) \f$, etc.,
       *  are some fields that depend on position, \f$ s \f$ is some scalar
       *  function, and \f$ \phi \f$ is some basis.
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
       *               - "Value Name" is the name of the scalar value being
       *                 integrated (\f$ s \f$),
       *               - "Basis" is the basis that you'd like to use (\f$ \phi
                         \f$),
       *               - "IR" is the integration rule that you'd like to use,
       *               - "Multiplier" is the scalar multiplier out in front of
       *                 the integral you're computing (\f$ M \f$), and
       *               - "Field Multipliers" is an optional list of names of
       *                 fields that are multipliers out in front of the
       *                 integral you're computing (\f$ a(x) \f$, \f$ b(x) \f$,
       *                 etc.).
       */
      Integrator_BasisTimesScalar(
        const Teuchos::ParameterList& p);

      /**
       *  \brief Post-Registration Setup.
       *
       *  Sets the number of nodes and quadrature points, sets the basis index,
       *  sets the contributed field data in the field manager (if applicable),
       *  and then creates the `tmp_` `Kokkos::View`.
       *
       *  \param[in] sd Essentially a list of `Workset`s, which are collections
       *                of cells (elements) that all live on a single process.
       *  \param[in] fm The field manager, used in setting the field data for
       *                the contributed field, if there is one.
       */
      void
      postRegistrationSetup(
        typename Traits::SetupData sd,
        PHX::FieldManager<Traits>& fm);

      /**
       *  \brief Evaluate Fields.
       *
       *  This actually performs the integration by looping over cells in the
       *  `Workset`, and then over bases and integration points on the cell.
       *
       *  \param[in] workset The `Workset` on which you're going to do the
       *                     integration.
       */
      void
      evaluateFields(
        typename Traits::EvalData workset);

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
      typedef typename EvalT::ScalarT ScalarT;

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
       *  \brief A field representing the scalar function we're integrating
       *         (\f$ s \f$).
       */
      PHX::MDField<const ScalarT, panzer::Cell, panzer::IP> scalar_;

      /**
       *  \brief The scalar multiplier out in front of the integral (\f$ M
                 \f$).
       */
      double multiplier_;

      /**
       *  \brief The (possibly empty) list of fields that are multipliers out
       *         in front of the integral (\f$ a(x) \f$, \f$ b(x) \f$, etc.).
       */
      std::vector<PHX::MDField<const ScalarT, panzer::Cell, panzer::IP>>
      fieldMults_;

      /**
       *  \brief The number of nodes for each cell.
       */
      int numNodes_;

      /**
       *  \brief The number of quadrature points for each cell.
       */
      int numQP_;

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
       *  \brief A temporary `Kokkos::View` that we'll use in computing the
       *         result  of this integral.
       */
      Kokkos::DynRankView<ScalarT, PHX::Device> tmp_;

  }; // end of class Integrator_BasisTimesScalar

} // end of namespace panzer

#endif // __Panzer_Integrator_BasisTimesScalar_decl_hpp__
