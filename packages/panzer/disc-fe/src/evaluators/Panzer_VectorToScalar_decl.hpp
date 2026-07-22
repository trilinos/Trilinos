// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_VECTOR_TO_SCALAR_DECL_HPP
#define PANZER_EVALUATOR_VECTOR_TO_SCALAR_DECL_HPP

#include <vector>
#include <string>
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Panzer_Evaluator_Macros.hpp"

namespace panzer {
    
template<typename EvalT, typename Traits>
class VectorToScalar
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    VectorToScalar(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;
  
  std::vector< PHX::MDField<ScalarT,Cell,Point> > scalar_fields;
  PHX::MDField<const ScalarT,Cell,Point,Dim> vector_field;

public:
 
  /**
   * \brief Tag only constructor for this class.
   */
  VectorToScalar(const PHX::FieldTag & input,
                 const std::vector<PHX::Tag<ScalarT>> & output);

}; // end of class VectorToScalar


/** This is a function constructor for an evaluator
  * that builds scalars from a single vector field. The user specifies
  * the layouts (assumed compatible) and then uses a postfix for each
  * of the scalar fields.
  * 
  * \param[in] vectorName Name of the vector 
  * \param[in] scalarPrefix Name to prefix to the scalar values
  * \param[in] postfix Vector specifying the postfix to use when naming
  *                    each scalar field
  * \param[in] vectorLayout Data layout for the vector field
  * \param[in] scalarLayout Data layout for the scalars
  */
template <typename EvalT,typename Traits>
Teuchos::RCP<PHX::Evaluator<Traits> > vectorToScalarEvaluator(const std::string & vectorName,
                                                              const std::string & scalarPrefix,
                                                              const std::vector<std::string> & postfix,
                                                              const Teuchos::RCP<const PHX::DataLayout> & vectorLayout,
                                                              const Teuchos::RCP<const PHX::DataLayout> & scalarLayout)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;

  RCP<std::vector<std::string> > scalarNames = rcp(new std::vector<std::string>);
  for(std::size_t i=0;i<postfix.size();i++)
    scalarNames->push_back(scalarPrefix+postfix[i]);

  Teuchos::ParameterList input;
  input.set("Vector Name",vectorName);
  input.set("Scalar Names",scalarNames.getConst());
  input.set("Data Layout Vector",rcp_const_cast<PHX::DataLayout>(vectorLayout));
  input.set("Data Layout Scalar",rcp_const_cast<PHX::DataLayout>(scalarLayout));

  return rcp(new VectorToScalar<EvalT,Traits>(input));
}

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
Teuchos::RCP<PHX::Evaluator<Traits> > vectorToScalarEvaluator(const std::string & vectorName,
                                                              const std::vector<std::string> & postfix,
                                                              const Teuchos::RCP<const PHX::DataLayout> & vectorLayout,
                                                              const Teuchos::RCP<const PHX::DataLayout> & scalarLayout)
{
  return vectorToScalarEvaluator<EvalT,Traits>(vectorName,vectorName,postfix,vectorLayout,scalarLayout);
}

}

#endif
