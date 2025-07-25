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

  // Strategy to pick the time index for the solution from the input file
  enum class IndexChoice {
    /// Specify an index and use it for the whole run
    Fixed,
    /// Interpolate using time values
    Interpolate,
    /// Pick the closest time value
    Closest
  };
  
  /**
     Loads a set of fields from an exodus file into the field
     manager. The list of exodus names are mapped to the field names.

     @param mesh STK mesh to read value from
     @param fieldNames Names of all fields in the field manager
     @param exodusNames Names of the fields in the exodus file 
     @param indexChoice determines how to pick the index in the exodus file.
     @param indexValue if indexChoice is set to "Fixed", this is the value to use.
   */
  GatherExodusCellDataToIP(const Teuchos::RCP<const panzer_stk::STK_Interface> & mesh,
                           const std::vector<std::string>& fieldNames,
                           const std::vector<std::string>& exodusNames,
                           const Teuchos::RCP<panzer::IntegrationRule>& integrationRule,
                           const GatherExodusCellDataToIP::IndexChoice& indexChoice,
                           const int indexValue = 0);
  
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
  const GatherExodusCellDataToIP::IndexChoice indexChoice_;
  const int indexValue_;
};

}

// **************************************************************
#endif
