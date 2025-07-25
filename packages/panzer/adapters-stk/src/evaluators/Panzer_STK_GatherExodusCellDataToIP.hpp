// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __PANZER_STK_GatherExodusCellDataToIP_decl_HPP__
#define __PANZER_STK_GatherExodusCellDataToIP_decl_HPP__

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_STK_Interface.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer_stk {

/** This class is a gathers a cell-based field from a STK Mesh 
  */
template<typename EvalT, typename Traits> 
class GatherExodusCellDataToIP
  : public panzer::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<EvalT, Traits> {
  
public:

  /**
     Loads a set of fields from an exodus file into the field
     manager. The list of exodus names are mapped to the field names.
     The time index from the exodus file for these fields corresponds
     to the restart value in the exodus mesh reader.

     @param mesh STK mesh to read value from
     @param fieldNames Names of all fields in the field manager
     @param exodusNames Names of the fields in the exodus file 
   */
  GatherExodusCellDataToIP(const Teuchos::RCP<const panzer_stk::STK_Interface> & mesh,
                           const std::vector<std::string>& fieldNames,
                           const std::vector<std::string>& exodusNames,
                           const Teuchos::RCP<panzer::IntegrationRule>& integrationRule);
  
  void postRegistrationSetup(typename Traits::SetupData d,
			     PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData d);

private:
  typedef typename EvalT::ScalarT ScalarT;
  typedef panzer_stk::STK_Interface::SolutionFieldType VariableField;
  
  const Teuchos::RCP<const STK_Interface> mesh_;
  const std::vector<std::string> exodusNames_;
  std::vector<PHX::MDField<ScalarT,panzer::Cell,panzer::IP>> gatherFields_;
  std::vector<VariableField*> stkFields_;
};

}

// **************************************************************
#endif
