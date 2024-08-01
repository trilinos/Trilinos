// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_COPY_HPP
#define PANZER_EVALUATOR_COPY_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Panzer_Evaluator_Macros.hpp"

namespace panzer {
    
/** Copies the contents of one field to anouther with a different
  * name. This basically allows easy renaming of fields. The constructor
  * takes a parameter list of the form
    \verbatim
    <ParameterList>
      <Parameter name="Source Name" type="string" value="<Input Field Name>"/>
      <Parameter name="Destination Name" type="string" value="<Output Field Name>"/>
      <Parameter name="Data Layout" type="RCP<PHX::DataLayout>" value="<Pointer to data layout describing input and output fields>"/>
    <ParameterList/>
    \endverbatim
  */
template<typename EvalT, typename Traits>
class Copy
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Copy(
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
  
  PHX::MDField<const ScalarT> input;
  PHX::MDField<ScalarT> output;

}; // end of class Copy


}

#endif
