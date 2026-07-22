// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_CellAverage_hpp__
#define __Panzer_CellAverage_hpp__

#include <string>
#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Kokkos_DynRankView.hpp"

#include "Panzer_Evaluator_Macros.hpp"

namespace panzer {
    
/** This integrates a scalar quanity over each cell.
  * It is useful for comptuing integral responses.

  \verbatim
    <ParameterList>
      <Parameter name="Average Name" type="string" value="<Name to give to the average field>"/>
      <Parameter name="Field Name" type="string" value="<Name of field to find average of>"/>
      <Parameter name="IR" type="RCP<IntegrationRule>" value="<user specified IntegrationRule>"/>
      <Parameter name="Multiplier" type="double" value="<Scaling factor, default=1>"/>
      <Parameter name="Field Multipliers" type="RCP<const vector<string> >" value="<Other scalar multiplier fields>"/>
    </ParameterList>
  \endverbatim
  */
template<typename EvalT, typename Traits>
class CellAverage
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    CellAverage(
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
  
  PHX::MDField<ScalarT,Cell> average;  // result
    
  PHX::MDField<const ScalarT,Cell,IP> scalar; // function to be integrated

  std::vector<PHX::MDField<const ScalarT,Cell,IP> > field_multipliers;
  double multiplier;

  std::size_t num_qp;
  std::size_t quad_index;
  int quad_order;
 
public:
  // for testing purposes
  const PHX::FieldTag & getFieldTag() const 
  { return average.fieldTag(); }

private:
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class CellAverage


/** This is a function constructor for an evaluator
  * that builds scalars from a single vector field. The user specifies
  * the layouts (assumed compatible) and then uses a postfix for each
  * of the scalar fields.
  * 
  * \param[in] vectorName Name of the vector 
  * \param[in] postfix Vector specifying the postfix to use when naming
  *                    each scalar field
  * \param[in] vectorLayout Data layout for the vector field
  * \param[in] scalarLayout Data layout for the scalars
  */
template <typename EvalT,typename Traits>
Teuchos::RCP<PHX::Evaluator<Traits> > cellAverageEvaluator(const std::string & averageName,
                                                           const std::string & fieldName,
                                                           const Teuchos::RCP<const panzer::IntegrationRule> & ir)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;

  Teuchos::ParameterList input;
  input.set("Average Name",averageName);
  input.set("Field Name",fieldName);
  input.set("IR",rcp_const_cast<panzer::IntegrationRule>(ir));

  return rcp(new CellAverage<EvalT,Traits>(input));
}

}

#endif
