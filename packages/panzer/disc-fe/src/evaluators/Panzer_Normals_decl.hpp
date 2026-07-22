// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_NORMALS_DECL_HPP
#define PANZER_EVALUATOR_NORMALS_DECL_HPP

#include <string>
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Panzer_Evaluator_Macros.hpp"

namespace panzer {
    
/** Compute normals on a particular side of an element.
  * By default the normals are normalized. A second option
  * would be for the normals to be unormalized values.
  
    <ParameterList name="Name" type="string" value="<Name to give to the normals field>"/>
    <ParameterList name="Side Id" type="int" value="<side id to use for computing normals>"/>
    <ParameterList name="IR" type="RCP<IntegrationRule>" value="<user specified IntegrationRule>"/>
    <ParameterList name="Normalize" type="bool" value="true"/>
  
  * The Name used to define the normals field is specified by "Name"
  * and the data layout is defined by the dl_vector field in the IntegrationRule.
  * The side ID must be legitimate for this topology and will be used to
  * construct an outward facing normal. The normals will be normalized by
  * default.  However, if Normalize=false then the resulting normals will have
  * the determinant of the side jacobian built in thus they correspond to a
  * differential on the side.
  */
template<typename EvalT, typename Traits>
class Normals
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Normals(
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

  int side_id;
  int quad_order, quad_index;

  std::size_t num_qp, num_dim;

  PHX::MDField<ScalarT,Cell,Point,Dim> normals;
  PHX::View<ScalarT*> scratch;
  bool normalize;

public:
  // for testing purposes
  const PHX::FieldTag & getFieldTag() const 
  { return normals.fieldTag(); }

}; // end of class Normals


}

#endif
