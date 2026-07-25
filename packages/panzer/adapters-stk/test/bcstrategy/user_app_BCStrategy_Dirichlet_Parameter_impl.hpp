// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Panzer_PhysicsBlock.hpp"

#include "Panzer_BasisIRLayout.hpp"

// Evaluators
#include "Panzer_Parameter.hpp"

#include "Phalanx_MDField.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

// ***********************************************************************
template <typename EvalT>
user_app::BCStrategy_Dirichlet_Parameter<EvalT>::
BCStrategy_Dirichlet_Parameter(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data) :
  panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>(bc,global_data),
  global_data_(global_data)
{
  TEUCHOS_ASSERT( (this->m_bc.strategy() == "Parameter") );
}

// ***********************************************************************
template <typename EvalT>
void user_app::BCStrategy_Dirichlet_Parameter<EvalT>::
setup(const panzer::PhysicsBlock& side_pb,
      const Teuchos::ParameterList& /* user_data */)
{
  using Teuchos::RCP;
  using std::vector;
  using std::string;
  using std::pair;

  // add in gather DOF
  this->addDOF(this->m_bc.equationSetName());              // DOF Name

  // add in the target
  this->addTarget("Parameter_"+this->m_bc.equationSetName(),  // Target Name
                  this->m_bc.equationSetName(),              // DOF Name
                  "Residual_"+this->m_bc.identifier());      // Residual Name

  // find the basis for this dof 
  const vector<pair<string,RCP<panzer::PureBasis> > >& dofs = side_pb.getProvidedDOFs();

  for (vector<pair<string,RCP<panzer::PureBasis> > >::const_iterator dof_it = 
	 dofs.begin(); dof_it != dofs.end(); ++dof_it) {
    if (dof_it->first == this->m_bc.equationSetName())
      this->basis = dof_it->second;
  }



  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(this->basis), std::runtime_error,
		     "Error the name \"" << this->m_bc.equationSetName()
		     << "\" is not a valid DOF for the boundary condition:\n"
		     << this->m_bc << "\n");

}

// ***********************************************************************
template <typename EvalT>
void user_app::BCStrategy_Dirichlet_Parameter<EvalT>::
buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
			   const panzer::PhysicsBlock& side_pb,
			   const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& /* factory */,
			   const Teuchos::ParameterList& /* models */,
			   const Teuchos::ParameterList& user_data) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  const std::map<int,RCP< panzer::IntegrationRule > >& ir_map = side_pb.getIntegrationRules();
  TEUCHOS_ASSERT(ir_map.size() == 1); 
  const auto ir = ir_map.begin()->second;

  // provide a constant target value to map into residual
  {
    const std::string param_name = this->m_bc.params()->template get<std::string>("Parameter Name");
    const std::string field_name = "Parameter_" + this->m_bc.equationSetName();
    
    RCP< PHX::Evaluator<panzer::Traits> > op_ip = 
      rcp(new panzer::Parameter<EvalT,panzer::Traits>(param_name, field_name, ir->dl_scalar, *global_data_->pl));
    
    RCP<const panzer::BasisIRLayout> basis_irl = basisIRLayout(basis,*ir);
    RCP< PHX::Evaluator<panzer::Traits> > op_basis =
        rcp(new panzer::Parameter<EvalT,panzer::Traits>(param_name,field_name,basis->functional,*global_data_->pl));
 
    this->template registerEvaluator<EvalT>(fm, op_ip);
    this->template registerEvaluator<EvalT>(fm, op_basis);
  }

}
