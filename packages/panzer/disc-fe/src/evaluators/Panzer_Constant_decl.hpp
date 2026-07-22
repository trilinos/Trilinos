// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_CONSTANT_DECL_HPP
#define PANZER_EVALUATOR_CONSTANT_DECL_HPP

#include "PanzerDiscFE_config.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Panzer_Evaluator_Macros.hpp"

namespace panzer {
    
/** Build a constant Phalanx field using a specified data layout
  
    \verbatim
    <ParameterList>
       <Parameter name="Name" type="string" value=(required)/>
       <Parameter name="Value" type="double" value="(required)"/>
       <Parameter name="Data Layout" type="RCP<PHX::DataLayout>" value=(required)/>
    </ParameterList>
    \endverbatim
  */
template<typename EvalT, typename Traits>
class Constant
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Constant(const Teuchos::ParameterList& p);

    void
    postRegistrationSetup(typename Traits::SetupData d,
			  PHX::FieldManager<Traits>& fm);

    void
    evaluateFields(typename Traits::EvalData d);

  private:

  using ScalarT = typename EvalT::ScalarT;
  
  double value;
  
  PHX::MDField<ScalarT> constant;
  
}; // end of class Constant


}

#endif
