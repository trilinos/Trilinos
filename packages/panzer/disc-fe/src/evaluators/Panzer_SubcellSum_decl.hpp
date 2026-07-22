// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_SubcellSum_decl_hpp__
#define __Panzer_SubcellSum_decl_hpp__

#include <string>

#include "PanzerDiscFE_config.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_FieldPattern.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Panzer_Evaluator_Macros.hpp"

namespace panzer {
    
/** This performs a sum over all the fields limited to the subcell
  * specified in the workset. It is useful for computing high-order
  * surface integrals as responses. 
  *
  * The field specified with "Sum Name" will be dimensioned as the number
  * of cells in the workset. The "Field Name" object is dimension as the number
  * of cells by the number of basis functions specified by the "Basis" object.
  * The "Evaluate On Closure" indicates if the subcells are to use the closure
  * index (i.e. all subcells of lesser dimension contained within a subcell) or
  * simply sum on those fields on the subcell proper.

  \verbatim
    <ParameterList>
      <Parameter name="Sum Name" type="string" value="<Name to give to the summed field>"/>
      <Parameter name="Field Name" type="string" value="<Name of field to sum>"/>
      <Parameter name="Basis" type="RCP<const PureBasis>" value="<user specified PureBasis object>"/>
      <Parameter name="Multiplier" type="double" value="<Scaling factor, default=1>"/>
      <Parameter name="Evaluate On Closure" type="bool" value="false"/>
    </ParameterList>
  \endverbatim
  */
template<typename EvalT, typename Traits>
class SubcellSum
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    SubcellSum(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;
  
  PHX::MDField<ScalarT,Cell> outField;  // result
    
  PHX::MDField<const ScalarT,Cell,BASIS> inField; // function to be integrated

  double multiplier;

public:

  const PHX::FieldTag & getFieldTag() const 
  { return outField.fieldTag(); }

private:
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;
 
  // This is used to lookup closure indices (local Ids that live on a subcell)
  Teuchos::RCP<const panzer::FieldPattern> fieldPattern_;
  
  // evalaute on the "closure" of the indicated sub-cells
  bool evaluateOnClosure_;

}; // end of class SubcellSum


}

#endif
