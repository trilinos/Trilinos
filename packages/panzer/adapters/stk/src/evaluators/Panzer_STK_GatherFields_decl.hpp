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

#ifndef __PANZER_STK_GatherFields_decl_HPP__
#define __PANZER_STK_GatherFields_decl_HPP__

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_STK_Interface.hpp"

namespace panzer_stk {

/** This class is a gather operation from the mesh. It
  * takes a set of field names and basis objects and
  * then reads them in from the mesh object.
  *
  * The constructor takes a STK_Interface RCP and parameter list
  * that is required to contain the following two fields
  * "Field Names" of type <code>Teuchos::RCP<std::vector<std::string> ></code>
  * and "Basis" of type <code>Teuchos::RCP<panzer::Basis></code>.
  */
template<typename EvalT, typename Traits> 
class GatherFields
  : public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<panzer::Traits::Residual, Traits> {

public:
  GatherFields(const Teuchos::RCP<const panzer_stk::STK_Interface> & mesh,const Teuchos::ParameterList & p) {}
  
  void postRegistrationSetup(typename Traits::SetupData d,
			     PHX::FieldManager<Traits>& vm) {}
  
  void evaluateFields(typename Traits::EvalData d) { TEUCHOS_ASSERT(false); }
public:

  GatherFields();
};

// **************************************************************
// **************************************************************
// * Specializations
// **************************************************************
// **************************************************************


// **************************************************************
// Residual 
// **************************************************************
template<typename Traits>
class GatherFields<panzer::Traits::Residual,Traits>
  : public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<panzer::Traits::Residual, Traits> {
   
  
public:
  
  GatherFields(const Teuchos::RCP<const panzer_stk::STK_Interface> & mesh,const Teuchos::ParameterList & p);
  
  void postRegistrationSetup(typename Traits::SetupData d,
			     PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData d);

private:
  typedef typename panzer::Traits::Residual::ScalarT ScalarT;
  typedef panzer_stk::STK_Interface::SolutionFieldType VariableField;

  std::vector< PHX::MDField<ScalarT,panzer::Cell,panzer::NODE> > gatherFields_;
  std::vector<VariableField*> stkFields_;
 
  Teuchos::RCP<const STK_Interface> mesh_;

  GatherFields();
};

}

// **************************************************************
#endif
