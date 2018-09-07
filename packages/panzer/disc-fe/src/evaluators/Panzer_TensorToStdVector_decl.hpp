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

#ifndef PANZER_EVALUATOR_TENSOR_TO_STD_VECTOR_DECL_HPP
#define PANZER_EVALUATOR_TENSOR_TO_STD_VECTOR_DECL_HPP

#include <vector>
#include <string>

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

namespace panzer {
    
/*! \brief Transform at Tensor to a std::vector of PHX-vectors
 *
 *  Since Phalanx/Panzer heavily relies on componentwise compuations,
 *  a tensor PHX::MDField<ScalarT,Cell,Point,Dim,Dim> is often represented by
 *  a std::vector<PHX::MDField<ScalarT,Cell,Point,Dim> >.
 *
 *  This class transforms a tensor to a std::vector representation.
 *
 *  \author mayr.mt \date 09/2015
 */
template<typename EvalT, typename Traits>
class TensorToStdVector
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    TensorToStdVector(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

    //! Tensor (to be distributed to vector)
    PHX::MDField<const ScalarT,Cell,Point,Dim,Dim> tensor_field;

    //! Vector (to be filled)
    std::vector<PHX::MDField<ScalarT,Cell,Point,Dim> > vector_fields;

}; // end of class TensorToStdVector


/** This is a function constructor for an evaluator
  * that builds vectors from a single tensor field. The user specifies
  * the layouts (assumed compatible) and then uses a postfix for each
  * of the vector fields.
  * 
  * \param[in] tensorName Name of the tensor
  * \param[in] vectorPrefix Name to prefix to the vector values
  * \param[in] postfix Vector specifying the postfix to use when naming
  *                    each vector field
  * \param[in] tensorLayout Data layout for the tensor field
  * \param[in] vectorLayout Data layout for the vectors
  */
template <typename EvalT,typename Traits>
Teuchos::RCP<PHX::Evaluator<Traits> > tensorToStdVectorEvaluator(const std::string & tensorName,
                                                              const std::string & vectorPrefix,
                                                              const std::vector<std::string> & postfix,
                                                              const Teuchos::RCP<const PHX::DataLayout> & tensorLayout,
                                                              const Teuchos::RCP<const PHX::DataLayout> & vectorLayout)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;

  RCP<std::vector<std::string> > vectorNames = rcp(new std::vector<std::string>);
  for(std::size_t i=0;i<postfix.size();i++)
    vectorNames->push_back(vectorPrefix+postfix[i]);

  Teuchos::ParameterList input;
  input.set("Tensor Name", tensorName);
  input.set("Vector Names", vectorNames.getConst());
  input.set("Data Layout Tensor",rcp_const_cast<PHX::DataLayout>(tensorLayout));
  input.set("Data Layout Vector",rcp_const_cast<PHX::DataLayout>(vectorLayout));

  return rcp(new TensorToStdVector<EvalT,Traits>(input));
}

/** This is a function constructor for an evaluator
  * that builds vectors from a single tensor field. The user specifies
  * the layouts (assumed compatible) and then uses a postfix for each
  * of the vector fields.
  * 
  * \param[in] tensorName Name of the tensor
  * \param[in] postfix Vector specifying the postfix to use when naming
  *                    each vector field
  * \param[in] tensorLayout Data layout for the tensor field
  * \param[in] vectorLayout Data layout for the vectors
  */
template <typename EvalT,typename Traits>
Teuchos::RCP<PHX::Evaluator<Traits> > tensorToStdVectorEvaluator(const std::string & tensorName,
                                                              const std::vector<std::string> & postfix,
                                                              const Teuchos::RCP<const PHX::DataLayout> & tensorLayout,
                                                              const Teuchos::RCP<const PHX::DataLayout> & vectorLayout)
{
  return tensorToStdVectorEvaluator<EvalT,Traits>(tensorName,tensorName,postfix,tensorLayout,vectorLayout);
}

}

#endif
