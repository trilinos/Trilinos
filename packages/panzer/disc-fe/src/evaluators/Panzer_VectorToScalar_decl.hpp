// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
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
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

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
