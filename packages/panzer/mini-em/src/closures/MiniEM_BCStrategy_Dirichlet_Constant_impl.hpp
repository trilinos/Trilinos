// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _MiniEM_BCStrategy_Dirichlet_Constant_impl_hpp_
#define _MiniEM_BCStrategy_Dirichlet_Constant_impl_hpp_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Panzer_PhysicsBlock.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"

// Evaluators
#include "Panzer_Constant.hpp"
#include "Panzer_ConstantVector.hpp"

#include "Phalanx_MDField.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

// ***********************************************************************
template <typename EvalT>
mini_em::BCStrategy_Dirichlet_Constant<EvalT>::
BCStrategy_Dirichlet_Constant(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data) :
  panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>(bc,global_data)
{
  TEUCHOS_ASSERT(this->m_bc.strategy() == "Constant");
}

// ***********************************************************************
template <typename EvalT>
void mini_em::BCStrategy_Dirichlet_Constant<EvalT>::
setup(const panzer::PhysicsBlock& side_pb,
      const Teuchos::ParameterList& /* user_data */)
{
  using Teuchos::RCP;
  using std::vector;
  using std::string;
  using std::pair;

/*
  // need the dof value to form the residual
  this->required_dof_names.push_back(this->m_bc.equationSetName());

  // unique residual name
  this->residual_name = "Residual_" + this->m_bc.identifier();

  // map residual to dof 
  this->residual_to_dof_names_map[residual_name] = this->m_bc.equationSetName();

  // map residual to target field
  this->residual_to_target_field_map[residual_name] = "Constant_" + this->m_bc.equationSetName();
*/

  // gather the DOF
  this->addDOF(this->m_bc.equationSetName());              // DOF Name

  // add in the targert
  this->addTarget("Constant_"+this->m_bc.equationSetName(),  // Target Name
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
void mini_em::BCStrategy_Dirichlet_Constant<EvalT>::
buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
			   const panzer::PhysicsBlock& /* pb */,
			   const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& /* factory */,
			   const Teuchos::ParameterList& /* models */,
			   const Teuchos::ParameterList& /* user_data */) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // provide a constant target value to map into residual
  if(basis->isScalarBasis()) 
  {
    // build the scalar
    {
      ParameterList p("BC Constant Dirichlet");
      p.set("Name", "Constant_" + this->m_bc.equationSetName());
      p.set("Data Layout", basis->functional);
      p.set("Value", this->m_bc.params()->template get<double>("Value"));

      RCP< PHX::Evaluator<panzer::Traits> > op =
          rcp(new panzer::Constant<EvalT,panzer::Traits>(p));

      this->template registerEvaluator<EvalT>(fm, op);
    }

    {
      ParameterList p("BC Constant Dirichlet");
      p.set("Name", "Constant_" + this->m_bc.equationSetName());
      p.set("Data Layout", basis->functional_grad);
      p.set("Value", this->m_bc.params()->template get<double>("Value"));
    
      RCP< PHX::Evaluator<panzer::Traits> > op =
          rcp(new panzer::Constant<EvalT,panzer::Traits>(p));
    
      this->template registerEvaluator<EvalT>(fm, op);
    }
  }
  else if(basis->isVectorBasis())
  {
    ParameterList p("BC Constant Vector Dirichlet");
    p.set("Name", "Constant_" + this->m_bc.equationSetName());
    p.set("Data Layout", basis->functional_grad);

    p.set("Value X", this->m_bc.params()->template get<double>("Value X"));
    if(basis->dimension()>1)
      p.set("Value Y", this->m_bc.params()->template get<double>("Value Y"));
    if(basis->dimension()>2)
      p.set("Value Z", this->m_bc.params()->template get<double>("Value Z"));
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::ConstantVector<EvalT,panzer::Traits>(p));
    
    fm.template registerEvaluator<EvalT>(op);
  }
  else
    TEUCHOS_ASSERT(false);

}
#endif
