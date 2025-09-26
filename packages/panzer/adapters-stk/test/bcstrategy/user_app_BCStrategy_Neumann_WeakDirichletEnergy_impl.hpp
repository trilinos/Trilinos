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
#include "user_app_Evaluator_EnergyFlux.hpp"
#include "Panzer_DOF.hpp"
#include "Panzer_DOFGradient.hpp"
#include "Panzer_ScalarToVector.hpp"
#include "Panzer_Integrator_BasisTimesScalar.hpp"
#include "Panzer_Integrator_GradBasisDotVector.hpp"

#include "Phalanx_MDField.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

// ***********************************************************************
template <typename EvalT>
user_app::BCStrategy_Neumann_WeakDirichletEnergy<EvalT>::
BCStrategy_Neumann_WeakDirichletEnergy(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data) :
  panzer::BCStrategy_Neumann_DefaultImpl<EvalT>(bc,global_data)
{

}

// ***********************************************************************
template <typename EvalT>
void user_app::BCStrategy_Neumann_WeakDirichletEnergy<EvalT>::
setup(const panzer::PhysicsBlock& side_pb,
      const Teuchos::ParameterList& /* user_data */)
{
  // need the dof value to form the residual
  this->requireDOFGather(this->m_bc.equationSetName());

  const std::string residual_name = "Residual_" + this->m_bc.identifier();
  const std::string dof_name = this->m_bc.equationSetName();
  const std::string flux_name = "WeakDirichletEnergy_" + this->m_bc.equationSetName();
  const int integration_order = this->m_bc.params()->template get<int>("Integration Order");

  this->addResidualContribution(residual_name,dof_name,flux_name,integration_order,side_pb);
}

// ***********************************************************************
template <typename EvalT>
void user_app::BCStrategy_Neumann_WeakDirichletEnergy<EvalT>::
buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
			   const panzer::PhysicsBlock& side_pb,
			   const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
			   const Teuchos::ParameterList& models,
			   const Teuchos::ParameterList& user_data) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::make_rcp;

  const std::vector<std::tuple<std::string,std::string,std::string,int,Teuchos::RCP<panzer::PureBasis>,Teuchos::RCP<panzer::IntegrationRule>>> data =
    this->getResidualContributionData();

  std::string residual_name = std::get<0>(data[0]);
  std::string dof_name = std::get<1>(data[0]);
  std::string flux_name = std::get<2>(data[0]);
  auto basis = std::get<4>(data[0]);
  Teuchos::RCP<panzer::IntegrationRule> ir = std::get<5>(data[0]);
  auto basis_ir_layout = panzer::basisIRLayout(basis,*ir);

  // Enable the weak dirichlet, otherwise use as outflow condition
  bool apply_weak_dirichlet = true;
  double value = 0.0;
  {
    Teuchos::ParameterList bc_params = *(this->bc().params());
    Teuchos::ParameterList validParams;
    validParams.set<int>("Integration Order",2,"Integration order of accuracy.");
    validParams.set<bool>("Apply Weak Dirichlet", true, "If set to true, use the Velocity Values in this parameter list for the Dirichlet BC. If false, then the DOF values are used and it acts like an outflow bc.");
    validParams.set<double>("Value",0.0);
    bc_params.validateParametersAndSetDefaults(validParams);
    apply_weak_dirichlet = bc_params.template get<bool>("Apply Weak Dirichlet");
    value = bc_params.template get<double>("Value");
  }

  // Add in the closure models
  {
    side_pb.buildAndRegisterClosureModelEvaluatorsForType<EvalT>(fm,factory,models,user_data);
    // side_pb.buildAndRegisterDOFProjectionsToIPEvaluators(fm,Teuchos::null,user_data);
  }

  // weak form of energy flux on boundary
  {
    // DOF names are currently hard coded for a particular test problem
    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new user_app::EnergyFluxEvaluator<EvalT,panzer::Traits>(flux_name,
                                                                  dof_name,
                                                                  basis_ir_layout->name(),
                                                                  "GRAD_OTHER_TEMPERATURE",
                                                                  "OTHER_U",
                                                                  "OTHER_DENSITY",
                                                                  "OTHER_HEAT_CAPACITY",
                                                                  "OTHER_THERMAL_CONDUCTIVITY",
                                                                  ir->dl_vector,
                                                                  ir->dl_scalar,
                                                                  apply_weak_dirichlet,
                                                                  value));

    this->template registerEvaluator<EvalT>(fm, op);
  }

  // Combine into velocity vector
  {
    Teuchos::ParameterList params;
    Teuchos::RCP<std::vector<std::string>> scalar_names = Teuchos::make_rcp<std::vector<std::string>>();
    scalar_names->push_back("OTHER_UX");
    scalar_names->push_back("OTHER_UY");
    params.set<Teuchos::RCP<const std::vector<std::string>>>("Scalar Names",scalar_names);
    params.set("Vector Name", "OTHER_U");
    params.set("Data Layout Scalar",ir->dl_scalar);
    params.set("Data Layout Vector",ir->dl_vector);
    auto op = Teuchos::make_rcp<panzer::ScalarToVector<EvalT,panzer::Traits>>(params);
    this->template registerEvaluator<EvalT>(fm, op);
  }

  // Compute DOF and GRAD_DOF at integration points
  {
    // New path requires partitioned worksets
    bool use_new_workset_ctor_path=false;
    if (use_new_workset_ctor_path) {
      // DOF
      panzer::BasisDescriptor basis_descriptor(basis->order(),basis->type());
      auto layout_basis = basis->functional;
      auto layout_ip = ir->dl_scalar;
      PHX::Tag<typename EvalT::ScalarT> dof(dof_name,layout_basis);
      PHX::Tag<typename EvalT::ScalarT> dof_at_ip(dof_name,layout_ip);
      auto op_dof = Teuchos::make_rcp<panzer::DOF<EvalT,panzer::Traits>>(dof,dof_at_ip,basis_descriptor,*ir);
      this->template registerEvaluator<EvalT>(fm, op_dof);

      // GRAD_DOF
      auto layout_basis_grad = basis->functional_grad;
      auto layout_ip_grad = ir->dl_vector;
      PHX::Tag<typename EvalT::ScalarT> dof_grad(dof_name,layout_basis);
      PHX::Tag<typename EvalT::ScalarT> dof_grad_at_ip("GRAD_"+dof_name,layout_ip_grad);
      auto op_dof_grad = Teuchos::make_rcp<panzer::DOFGradient<EvalT,panzer::Traits>>(dof_grad,dof_grad_at_ip,basis_descriptor,*ir);
      this->template registerEvaluator<EvalT>(fm, op_dof_grad);
    }
    else {
      { // DOF
        Teuchos::ParameterList params;
        params.set("Name",dof_name);
        params.set("Basis",basis_ir_layout);
        params.set("IR",ir);
        auto op_dof = Teuchos::make_rcp<panzer::DOF<EvalT,panzer::Traits>>(params);
        this->template registerEvaluator<EvalT>(fm, op_dof);
      }
      { // GRAD_DOF
        Teuchos::ParameterList params;
        params.set("Name",dof_name);
        params.set("Gradient Name",std::string("GRAD_")+dof_name);
        params.set("Basis",basis_ir_layout);
        params.set("IR",ir);
        auto op_dof_grad = Teuchos::make_rcp<panzer::DOFGradient<EvalT,panzer::Traits>>(params);
        this->template registerEvaluator<EvalT>(fm, op_dof_grad);
      }
    }
  }

  // When applying the weak dirichlet, add in the extra terms for the
  // symmetric interior penalty method as shown in "Weak imposition of
  // Dirichlet boundary conditions in fluid mechanics," Bazilev and
  // Hughes, Computers and Fluids 36, 2007:
  if (apply_weak_dirichlet) {
    {
      Teuchos::ParameterList p;
      p.set("Residual Name",residual_name);
      p.set("Value Name",std::string("WEAK_DBC_PENALTY_TERM_")+dof_name);
      p.set("Basis",basis_ir_layout);
      p.set("IR",ir);
      p.set("Multiplier",1.0);
      auto evalType = panzer::EvaluatorStyle::CONTRIBUTES;
      auto op = Teuchos::make_rcp<panzer::Integrator_BasisTimesScalar<EvalT,panzer::Traits>>(p,evalType);
      this->template registerEvaluator<EvalT>(fm, op);
    }
    {
      Teuchos::ParameterList p;
      p.set("Residual Name",residual_name);
      p.set("Flux Name",std::string("WEAK_DBC_FLUX_TERM_")+dof_name);
      p.set("Basis",basis_ir_layout);
      p.set("IR",ir);
      p.set("Multiplier",-1.0);
      auto evalType = panzer::EvaluatorStyle::CONTRIBUTES;
      auto op = Teuchos::make_rcp<panzer::Integrator_GradBasisDotVector<EvalT,panzer::Traits>>(p,evalType);
      this->template registerEvaluator<EvalT>(fm, op);
    }
  }

}

// ***********************************************************************
template <typename EvalT>
void user_app::BCStrategy_Neumann_WeakDirichletEnergy<EvalT>::
postRegistrationSetup(typename panzer::Traits::SetupData /* d */,
		      PHX::FieldManager<panzer::Traits>& /* vm */)
{

}


// ***********************************************************************
template <typename EvalT>
void user_app::BCStrategy_Neumann_WeakDirichletEnergy<EvalT>::
evaluateFields(typename panzer::Traits::EvalData /* d */)
{

}
