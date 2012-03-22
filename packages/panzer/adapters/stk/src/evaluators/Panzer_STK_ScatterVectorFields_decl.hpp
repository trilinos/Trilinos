#ifndef __PANZER_STK_ScatterVectorFields_decl_HPP__
#define __PANZER_STK_ScatterVectorFields_decl_HPP__

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_PointRule.hpp"
#include "Panzer_STK_Interface.hpp"

namespace panzer_stk {

/** This class is a scatter operation to the mesh. It
  * takes a set of field names and basis objects and
  * then writes them to the mesh object.
  *
  * The constructor takes a STK_Interface RCP and parameter list
  * that is required to contain the following two fields
  * "Field Names" of type <code>Teuchos::RCP<std::vector<std::string> ></code>,
  * "Basis" of type <code>Teuchos::RCP<panzer::BasisIRLayout></code> and
  * "Mesh" of type <code>Teuchos::RCP<const panzer_stk::STK_Interface></code>.
  */
PHX_EVALUATOR_CLASS(ScatterVectorFields)
  typedef panzer_stk::STK_Interface::SolutionFieldType VariableField;

  std::vector<std::string> names_;
  std::vector< PHX::MDField<ScalarT,panzer::Cell,panzer::IP,panzer::Dim> > scatterFields_;
  PHX::MDField<ScalarT,panzer::Cell,panzer::IP,panzer::Dim> pointField_;
  Teuchos::RCP<STK_Interface> mesh_;

  int spatialDimension_;

public:
  
  ScatterVectorFields(const std::string & scatterName,
                      const Teuchos::RCP<STK_Interface> mesh,
                      const Teuchos::RCP<const panzer::PointRule> & pointRule,
                      const std::vector<std::string> & names);
 
PHX_EVALUATOR_CLASS_END

}

// **************************************************************
#endif
