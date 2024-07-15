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

template<typename EvalT, typename Traits>
class TestScatter
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    TestScatter(
      const Teuchos::ParameterList& p);

    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

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
