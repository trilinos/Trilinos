#ifndef PANZER_EVALUATOR_TEST_SCATTER_HPP
#define PANZER_EVALUATOR_TEST_SCATTER_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

namespace panzer {

PHX_EVALUATOR_CLASS(TestScatter)
  PHX::MDField<ScalarT,Cell,NODE> scatter_value;
  PHX::MDField<ScalarT,Cell,NODE> value;
  int localOffset;

  std::size_t num_nodes;
  static int offset;
PHX_EVALUATOR_CLASS_END

}

#include "Panzer_TestScatterT.hpp"

#endif
