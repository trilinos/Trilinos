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

#ifndef __PANZER_STK_ScatterFields_decl_HPP__
#define __PANZER_STK_ScatterFields_decl_HPP__

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Teuchos_ParameterList.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_STK_Interface.hpp"

namespace panzer_stk {

/** This class is a scatter operation to the mesh. It
  * takes a set of field names and basis objects and
  * then writes them to the mesh object. Note that <code>scaling</code> vector
  * must be the same length as the <code>names</code> vector. The scaling
  * is applied to each field.
  */
template <typename EvalT,typename TraitsT>
class ScatterFields : public PHX::EvaluatorWithBaseImpl<TraitsT>,
                      public PHX::EvaluatorDerived<EvalT, TraitsT>  { 
  typedef typename EvalT::ScalarT ScalarT;
  typedef panzer_stk::STK_Interface::SolutionFieldType VariableField;

  std::vector< PHX::MDField<ScalarT,panzer::Cell,panzer::NODE> > scatterFields_;
  Teuchos::RCP<STK_Interface> mesh_;

  std::vector<VariableField*> stkFields_;

  std::vector<double> scaling_;

  void initialize(const std::string & scatterName,
                  const Teuchos::RCP<STK_Interface> mesh,
                  const Teuchos::RCP<const panzer::PureBasis> & basis,
                  const std::vector<std::string> & names,
                  const std::vector<double> & scaling);

public:
  
  ScatterFields(const std::string & scatterName,
                const Teuchos::RCP<STK_Interface> mesh,
                const Teuchos::RCP<const panzer::PureBasis> & basis,
                const std::vector<std::string> & names);

  ScatterFields(const std::string & scatterName,
                const Teuchos::RCP<STK_Interface> mesh,
                const Teuchos::RCP<const panzer::PureBasis> & basis,
                const std::vector<std::string> & names,
                const std::vector<double> & scaling);

  void postRegistrationSetup(typename TraitsT::SetupData d, 
                             PHX::FieldManager<TraitsT>& fm);

  void evaluateFields(typename TraitsT::EvalData d);
}; 

}

// **************************************************************
#endif
