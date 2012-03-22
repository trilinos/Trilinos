#ifndef __PANZER_STK_ScatterCellAvgQuantity_decl_HPP__
#define __PANZER_STK_ScatterCellAvgQuantity_decl_HPP__

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"

#include "Panzer_STK_Interface.hpp"

namespace panzer_stk {

/** This class is a scatter operation to the mesh. It
  * takes a set of field names and an integration rule 
  * those quantities are then averaged over the cell and
  * written to the mesh.
  *
  * The constructor takes a STK_Interface RCP and parameter list
  * that is required to contain the following two fields
  * "Scatter Name" string specifying the name of this evaulator
  * "Field Names" of type this is a comma seperated list of strings,
  * "IR" of type <code>Teuchos::RCP<panzer::IntegrationRule></code> and
  * "Mesh" of type <code>Teuchos::RCP<const panzer_stk::STK_Interface></code>.
  */
PHX_EVALUATOR_CLASS(ScatterCellAvgQuantity)
  typedef panzer_stk::STK_Interface::SolutionFieldType VariableField; // this is weird, but the correct thing

  std::size_t numValues_;
 
  std::vector< PHX::MDField<ScalarT,panzer::Cell,panzer::Point> > scatterFields_;
  Teuchos::RCP<STK_Interface> mesh_;

  std::vector<VariableField*> stkFields_;
 
PHX_EVALUATOR_CLASS_END

}

// **************************************************************
#endif
