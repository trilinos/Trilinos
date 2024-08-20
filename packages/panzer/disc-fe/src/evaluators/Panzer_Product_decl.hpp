// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_PRODUCT_HPP
#define PANZER_PRODUCT_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Panzer_Evaluator_Macros.hpp"

namespace panzer {
    
/** Product of entries on a single data layout
 
    \verbatim
    <ParameterList>
      <ParameterList name="Product Name" type="string" value="<destination field name>"/>
      <ParameterList name="Values Names" type="Teuchos::RCP<std::vector<std::string> >" value="<Source field names>"/>
      <ParameterList name="Data Layout" type="Teuchos::RCP<PHX::DataLayout>" value="<data layout of all associated fields>"/>
      <ParameterList name="Scaling" type="double" value="<data of the scaling>/> <!-- Optional -->
    </ParameterList>
    \endverbatim
  */
template<typename EvalT, typename Traits>
class Product
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Product(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

  using ScalarT = typename EvalT::ScalarT;

  double scaling;
  PHX::MDField<ScalarT> product;
  std::vector< PHX::MDField<const ScalarT> > values;

}; // end of class Product


}

#endif
