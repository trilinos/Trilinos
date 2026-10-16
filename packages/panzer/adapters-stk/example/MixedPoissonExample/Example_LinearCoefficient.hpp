// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Example_LinearCoefficient_hpp__
#define __Example_LinearCoefficient_hpp__

#include "PanzerAdaptersSTK_config.hpp"

#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_FieldLibrary.hpp"

#include <string>

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace Example {

  using panzer::Cell;
  using panzer::Point;

/** A spatially varying coefficient, 1 + x_d, at the integration points,
  * where x_d is the coordinate in direction d.
  */
template<typename EvalT, typename Traits>
class LinearCoefficient : public panzer::EvaluatorWithBaseImpl<Traits>,
                          public PHX::EvaluatorDerived<EvalT, Traits>  {

public:
    LinearCoefficient(const std::string & name,
                      const panzer::IntegrationRule & ir,
                      const int direction);

    void postRegistrationSetup(typename Traits::SetupData d,
                               PHX::FieldManager<Traits>& fm);

    void evaluateFields(typename Traits::EvalData d);


private:
  typedef typename EvalT::ScalarT ScalarT;

  PHX::MDField<ScalarT,Cell,Point> coeff;
  int ir_degree, ir_index, direction;
};

}

#include "Example_LinearCoefficient_impl.hpp"

#endif
