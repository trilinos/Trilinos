// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_TEST_SCATTER_DECL_HPP
#define PANZER_EVALUATOR_TEST_SCATTER_DECL_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Panzer_Evaluator_Macros.hpp"

namespace panzer {

/** A test-only evaluator that computes each cell's node average of an
  * input field and scatters a scaled copy of that average back out to
  * every node of an output field. Used for exercising the scatter
  * machinery in unit tests, not for production physics.
  */
template<typename EvalT, typename Traits>
class TestScatter
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    /// \brief Construct from a ParameterList specifying the input/output field names and shared data layout.
    TestScatter(
      const Teuchos::ParameterList& p);

    /// \brief Looks up the number of nodes per cell needed by evaluateFields(). Called once before the first evaluation.
    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

    /// \brief Computes the output field for this workset per the class description.
    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;
  PHX::MDField<ScalarT,Cell,NODE> scatter_value;
  PHX::MDField<const ScalarT,Cell,NODE> value;
  int localOffset;

  std::size_t num_nodes;
  static int offset;
}; // end of class TestScatter


}

#endif
