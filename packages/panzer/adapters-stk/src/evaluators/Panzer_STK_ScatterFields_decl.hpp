// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __PANZER_STK_ScatterFields_decl_HPP__
#define __PANZER_STK_ScatterFields_decl_HPP__

#include "Phalanx_config.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Teuchos_ParameterList.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_STK_Interface.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer_stk {

/** This class is a scatter operation to the mesh. It
  * takes a set of field names and basis objects and
  * then writes them to the mesh object. Note that <code>scaling</code> vector
  * must be the same length as the <code>names</code> vector. The scaling
  * is applied to each field.
  */
template <typename EvalT,typename TraitsT>
class ScatterFields : public panzer::EvaluatorWithBaseImpl<TraitsT>,
                      public PHX::EvaluatorDerived<EvalT, TraitsT>  { 
  typedef typename EvalT::ScalarT ScalarT;
  typedef panzer_stk::STK_Interface::SolutionFieldType VariableField;

  std::vector< PHX::MDField<const ScalarT,panzer::Cell,panzer::NODE> > scatterFields_;
  Teuchos::RCP<STK_Interface> mesh_;

  std::vector<double> scaling_;

  bool cellFields_;

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
