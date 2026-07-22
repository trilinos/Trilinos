// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_STK_SCATTER_CELL_AVG_VECTOR_DECL_HPP
#define PANZER_STK_SCATTER_CELL_AVG_VECTOR_DECL_HPP

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"

#include "Panzer_STK_Interface.hpp"

#include "Panzer_Evaluator_Macros.hpp"

namespace panzer_stk {

/** This class is a scatter operation to the mesh. It
  * takes a set of field names and an integration rule. 
  * Those quantities are components of vector fields. They are averaged over 
  * the cell and written to the mesh.
  *
  * The constructor takes a STK_Interface RCP and parameter list
  * that is required to contain the following fields:
  * "Scatter Name" string specifying the name of this evaulator
  * "Field Names" of type this is a comma seperated list of strings,
  * "IR" of type <code>Teuchos::RCP<panzer::IntegrationRule></code> and
  * "Mesh" of type <code>Teuchos::RCP<const panzer_stk::STK_Interface></code>.
  */
template<typename EvalT, typename Traits>
class ScatterCellAvgVector
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    ScatterCellAvgVector(
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

  // typedef panzer_stk::STK_Interface::SolutionFieldType VariableField; // this is weird, but the correct thing
  typedef panzer_stk::STK_Interface::VectorFieldType VariableField;
  
  std::size_t numValues_;
 
  // map of variable-name to scale-factors to be applied upon output. if
  // this is empty then no variable scaling will be performed. this
  // should be passed in as an object via the teuchos parameter list in
  // the ctor with the parameter name "Variable Scale Factors Map".
  Teuchos::RCP<std::map<std::string,double>> varScaleFactors_;

  std::vector< PHX::MDField<const ScalarT,panzer::Cell,panzer::Point,panzer::Dim> > scatterFields_;
  Teuchos::RCP<STK_Interface> mesh_;

  std::vector<VariableField*> stkFields_;
 
}; // end of class ScatterCellAvgVector


}

// **************************************************************
#endif
