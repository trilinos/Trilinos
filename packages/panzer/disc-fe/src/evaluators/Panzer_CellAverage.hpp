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
