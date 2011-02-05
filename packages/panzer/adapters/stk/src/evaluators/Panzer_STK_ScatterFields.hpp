#ifndef __PANZER_STK_ScatterFields_HPP__
#define __PANZER_STK_ScatterFields_HPP__

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"

namespace panzer_stk {

/** This class is a scatter operation to the mesh. It
  * takes a set of field names and basis objects and
  * then writes them to the mesh object.
  *
  * The constructor takes a STK_Interface RCP and parameter list
  * that is required to contain the following two fields
  * "Field Names" of type <code>Teuchos::RCP<std::vector<std::string> ></code>,
  * "Basis" of type <code>Teuchos::RCP<panzer::Basis></code> and
  * "Mesh" of type <code>Teuchos::RCP<const panzer_stk::STK_Interface></code>.
  */
PHX_EVALUATOR_CLASS(ScatterFields)
  typedef panzer_stk::STK_Interface::SolutionFieldType VariableField;

  std::vector< PHX::MDField<ScalarT,panzer::Cell,panzer::NODE> > scatterFields_;
  Teuchos::RCP<STK_Interface> mesh_;

  std::vector<VariableField*> stkFields_;
 
PHX_EVALUATOR_CLASS_END

}

// **************************************************************

#include "Panzer_STK_ScatterFieldsT.hpp"

// **************************************************************
#endif
