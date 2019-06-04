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
