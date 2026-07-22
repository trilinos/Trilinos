// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_CrossProduct_HPP
#define PANZER_EVALUATOR_CrossProduct_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Panzer_Evaluator_Macros.hpp"

namespace panzer {
    
  /** \brief Evaluates cross product at a set of points

      v_a \times v_b

    <Parameter name="Result Name" type="string" value="<Name to give to cross product field>"/>
    <Parameter name="Point Rule" type="RCP<const PointRule>" value="<user specified point rule>"/>
    <Parameter name="Vector A Name" type="string" value="<vector a name>"/>
    <Parameter name="Vector B Name" type="string" value="<vector b name>"/>
  */
template<typename EvalT, typename Traits>
class CrossProduct
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    CrossProduct(
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
  
  PHX::MDField<ScalarT> vec_a_cross_vec_b;
  PHX::MDField<const ScalarT> vec_a, vec_b;

  bool useScalarField;

  int num_pts;
  int num_dim;

}; // end of class CrossProduct


}

#endif
