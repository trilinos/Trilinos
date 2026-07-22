// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_DotProduct_DECL_HPP
#define PANZER_EVALUATOR_DotProduct_DECL_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Panzer_Evaluator_Macros.hpp"

namespace panzer {
    
/** \brief Evaluates dot product at a set of points

    v_a \cdot v_b

  <Parameter name="Result Name" type="string" value="<Name to give to dot product field>"/>
  <Parameter name="Point Rule" type="RCP<const PointRule>" value="<user specified point rule>"/>
  <Parameter name="Vector A Name" type="string" value="<vector a name>"/>
  <Parameter name="Vector B Name" type="string" value="<vector b name>"/>
  <Parameter name="Multiplier" type="double" value="Multiplier value"/>
  <Parameter name="Field Multiplier" type="string" value="Multiplier name"/>
*/
template<typename EvalT, typename Traits>
class DotProduct
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    DotProduct(
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
  
  PHX::MDField<ScalarT> vec_a_dot_vec_b;
  PHX::MDField<const ScalarT> vec_a, vec_b;
  PHX::MDField<const ScalarT> multiplier_field;

  int num_pts;
  int num_dim;

  bool multiplier_field_on;
  double multiplier_value;
}; // end of class DotProduct


/** \brief Build a dot product evaluator. Evaluates dot product at a set of points

    mv * fm * v_a \cdot v_b
   
  * \param[in] resultName Destination scalar field sized by <code>pr.dl_scalar</code>
  *                       containing the dot product
  * \param[in] pr Point rule class that defines the size of the fields
  * \param[in] vecA The a vector
  * \param[in] vecB The b vector
  * \param[in] multiplier Constant multiplier (mv above)
  * \param[in] fieldMultiplier Field to multiply by (fm above)
  */
template <typename EvalT,typename TraitsT>
Teuchos::RCP<DotProduct<EvalT,TraitsT> > buildEvaluator_DotProduct(const std::string & resultName,
                                                                     const panzer::PointRule & pr,
                                                                     const std::string & vecA,
                                                                     const std::string & vecB,
                                                                     double multiplier=1,
                                                                     const std::string & fieldMultiplier="");

}

#endif
