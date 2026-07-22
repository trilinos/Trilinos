// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_PROJECT_TO_FACES_DECL_HPP
#define PANZER_EVALUATOR_PROJECT_TO_FACES_DECL_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_CloneableEvaluator.hpp"

namespace panzer {

/** \brief Given a function stored as a vector and the tangents at each edge, project the vector onto the edge basis
*/
template<typename EvalT, typename Traits> 
class ProjectToFaces
  : public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<EvalT, Traits>,
    public CloneableEvaluator  {
   
public:
  
  ProjectToFaces(const Teuchos::ParameterList& p);
  
  void postRegistrationSetup(typename Traits::SetupData d,
			     PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData d);

  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  {return Teuchos::rcp(new ProjectToFaces<EvalT,Traits>(pl));}
  
private:

  typedef typename EvalT::ScalarT ScalarT;

  std::string dof_name_;
  Teuchos::RCP<const PureBasis> basis_;
  int num_faces_;
  int num_dim_;

  PHX::MDField<const ScalarT,Cell,BASIS,Dim> normals_;
  PHX::MDField<const ScalarT,Cell,BASIS,Dim> vector_values_;
  PHX::MDField<ScalarT,Cell,BASIS> result_;
};

}

// **************************************************************
#endif
