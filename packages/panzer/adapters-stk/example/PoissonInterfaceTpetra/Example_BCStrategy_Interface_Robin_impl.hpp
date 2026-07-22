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

#include "Panzer_PureBasis.hpp"

// Evaluators
#include "Panzer_Integrator_BasisTimesScalar.hpp"
#include "Panzer_Constant.hpp"
#include "Panzer_DOF.hpp"
#include "Panzer_Sum.hpp"
#include "Panzer_Product.hpp"
#include "Panzer_DOFGradient.hpp"
#include "Panzer_DotProduct.hpp"

#include "Phalanx_MDField.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

template <typename EvalT>
Example::BCStrategy_Interface_Robin<EvalT>::
BCStrategy_Interface_Robin(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data)
  : panzer::BCStrategy_Interface_DefaultImpl<EvalT>(bc, global_data)
{
  using std::string;

  TEUCHOS_ASSERT(this->m_bc.strategy() == "Robin Interface");

  const char* coeff_names[] = {"a", "b", "c"};
  {
    Teuchos::ParameterList vp;
    vp.set<string>("Coupling DOF Name", "INVALID_DOF", "Field to which to couple at interface");
    for (int i = 0; i < 3; ++i)
      vp.set<double>(coeff_names[i], 1.0, "Coefficient");
    vp.set<bool>("Nonlinear", false);
    bc.nonconstParams()->validateParametersAndSetDefaults(vp);
  }
  coupling_dof_name_ = bc.params()->get<string>("Coupling DOF Name");
  TEUCHOS_TEST_FOR_EXCEPTION(coupling_dof_name_ == "INVALID_DOF", std::logic_error,
                             "Must provide 'Coupling DOF Name'");
  for (int i = 0; i < 3; ++i) {
    // Negative so that a = 1, b = c = 0 means f matches phi's normal
    // derivative.
    coeffs_[i] = -bc.params()->get<double>(coeff_names[i]);
  }

  nonlinear_ = bc.params()->get<bool>("Nonlinear");
}

template <typename EvalT>
void Example::BCStrategy_Interface_Robin<EvalT>::
setup(const panzer::PhysicsBlock& side_pb,
      const Teuchos::ParameterList& /* user_data */)
{
  using Teuchos::RCP;
  using std::vector;
  using std::string;
  using std::pair;

  const int di = this->getDetailsIndex();

  dof_name_ = di == 0 ? this->m_bc.equationSetName() : this->m_bc.equationSetName2();
  other_dof_name_ = di == 1 ? this->m_bc.equationSetName() : this->m_bc.equationSetName2();

  this->requireDOFGather(dof_name_);
  this->requireDOFGather(coupling_dof_name_);

  // Unique residual name.
  const string residual_name = "Residual_" + this->m_bc.equationSetName();

  const std::map<int,RCP< panzer::IntegrationRule > >& ir = side_pb.getIntegrationRules();
  TEUCHOS_ASSERT(ir.size() == 1); 
  
  const int integration_order = ir.begin()->second->order();

  this->addResidualContribution(residual_name, dof_name_, "", integration_order, side_pb);
}

template <typename EvalT>
void Example::BCStrategy_Interface_Robin<EvalT>::
setCombineValues(Teuchos::ParameterList& p,
                 const std::string value_name1, const double scalar1,
                 const std::string value_name2, const double scalar2,
                 const std::string value_name3, const double scalar3)
{
  std::vector<std::string> values_names(2);
  values_names[0] = value_name1;
  values_names[1] = value_name2;
  if (value_name3 != "") values_names.push_back(value_name3);
  p.set< Teuchos::RCP<std::vector<std::string> > >(
    "Values Names", Teuchos::rcp(new std::vector<std::string>(values_names)));
  std::vector<double> scalars(2);
  scalars[0] = scalar1;
  scalars[1] = scalar2;
  if (values_names.size() > 2) scalars.push_back(scalar3);
  p.set< Teuchos::RCP<const std::vector<double> > >(
    "Scalars", Teuchos::rcp(new std::vector<double>(scalars)));
}

template <typename EvalT>
void Example::BCStrategy_Interface_Robin<EvalT>::
buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
			   const panzer::PhysicsBlock& pb,
			   const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& /* factory */,
			   const Teuchos::ParameterList& /* models */,
			   const Teuchos::ParameterList& /* user_data */) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::string; 

  const std::vector<std::tuple<string,string,string,int,Teuchos::RCP<panzer::PureBasis>,
    Teuchos::RCP<panzer::IntegrationRule> > > data = this->getResidualContributionData();

  const string residual_name = std::get<0>(data[0]);

  RCP<panzer::IntegrationRule> ir = std::get<5>(data[0]);
  RCP<const panzer::FieldLayoutLibrary> fll = pb.getFieldLibrary()->buildFieldLayoutLibrary(*ir);
  RCP<panzer::BasisIRLayout> basis = fll->lookupLayout(dof_name_);

  // Implement the interface condition by substituting the rhs of the condition
  // into the flux term for f. The flux term is the natural BC on this side of
  // the element block, so it doesn't appear in the following evaluators, just
  // like in the implementation of the constant Neumann boundary condition.
  if (this->getDetailsIndex() == 0) {
    const string
      my_normal_name = "my_normal",
      my_dof_grad_name = dof_name_ + "_gradient",
      coupling_dof_grad_name = coupling_dof_name_ + "_gradient",
      grad_sum_name = "grad_sum",
      normal_dot_coupling_grad_name = "normal_dot_coupling_grad",
      product_name = "product",
      sum_contributions_name = "sum_contributions";
    {
      ParameterList p("My DOF");
      p.set("Name", dof_name_);
      p.set("Basis", basis); 
      p.set("IR", ir);
      const RCP< PHX::Evaluator<panzer::Traits> >
        op = rcp(new panzer::DOF<EvalT,panzer::Traits>(p));
      this->template registerEvaluator<EvalT>(fm, op);
    }
    {
      ParameterList p("My Side Normal");
      p.set("Name", my_normal_name);
      p.set("Side ID", pb.cellData().side());
      p.set("IR", ir);
      p.set("Normalize", true);
      const RCP< PHX::Evaluator<panzer::Traits> >
        op = rcp(new panzer::Normals<EvalT,panzer::Traits>(p));
      this->template registerEvaluator<EvalT>(fm, op);
    }
    {
      ParameterList p("Coupling DOF Gradient");
      p.set("Name", coupling_dof_name_);
      p.set("Gradient Name", coupling_dof_grad_name);
      p.set("Basis", basis);
      p.set("IR", ir);
      const RCP< PHX::Evaluator<panzer::Traits> >
        op = rcp(new panzer::DOFGradient<EvalT,panzer::Traits>(p));
      this->template registerEvaluator<EvalT>(fm, op);
    }
    {
      ParameterList p("dot(grad phi, normal)");
      p.set("Result Name", normal_dot_coupling_grad_name);
      p.set("Vector A Name", coupling_dof_grad_name);
      p.set("Vector B Name", my_normal_name);
      p.set("Point Rule", Teuchos::rcp_dynamic_cast<const panzer::PointRule>(ir));
      const RCP< PHX::Evaluator<panzer::Traits> >
        op = rcp(new panzer::DotProduct<EvalT,panzer::Traits>(p));
      this->template registerEvaluator<EvalT>(fm, op);
    }
    if (nonlinear_) {
      ParameterList p("f_me dot(grad phi, normal)");
      p.set("Product Name", product_name);
      setCombineValues(p,
                       normal_dot_coupling_grad_name, 1,
                       dof_name_, 1);
      p.set("Data Layout", ir->dl_scalar);
      const RCP< PHX::Evaluator<panzer::Traits> >
        op = rcp(new panzer::Product<EvalT,panzer::Traits>(p));
      this->template registerEvaluator<EvalT>(fm, op);      
    }
    {
      ParameterList p("a dot(grad phi, normal) [f_me if nonlinear] + b f_me + c f_other");
      p.set("Sum Name", sum_contributions_name);
      setCombineValues(p,
                       nonlinear_ ? product_name : normal_dot_coupling_grad_name, coeffs_[0],
                       dof_name_, coeffs_[1],
                       other_dof_name_, coeffs_[2]);
      p.set("Data Layout", ir->dl_scalar);
      const RCP< PHX::Evaluator<panzer::Traits> >
        op = rcp(new panzer::Sum<EvalT,panzer::Traits>(p));
      this->template registerEvaluator<EvalT>(fm, op);
    }
    {
      ParameterList p("a dot(grad phi, normal) [f_me if nonlinear] + b f_me + c f_other residual");
      p.set("Residual Name", residual_name);
      p.set("Value Name", sum_contributions_name);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", 1.0);
      const RCP< PHX::Evaluator<panzer::Traits> >
        op = rcp(new panzer::Integrator_BasisTimesScalar<EvalT,panzer::Traits>(p));
      this->template registerEvaluator<EvalT>(fm, op);
    }
  } else {
    { // Get values on other side.
      ParameterList p("Other DOF");
      p.set("Name", dof_name_);
      p.set("Basis", basis); 
      p.set("IR", ir);
      const RCP< PHX::Evaluator<panzer::Traits> >
        op = rcp(new panzer::DOF<EvalT,panzer::Traits>(p));
      this->template registerEvaluator<EvalT>(fm, op);
    }
  }
}

template <typename EvalT>
void Example::BCStrategy_Interface_Robin<EvalT>::
buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& fm,
			                       const panzer::PhysicsBlock& pb,
				               const panzer::LinearObjFactory<panzer::Traits> & lof,
				               const Teuchos::ParameterList& user_data) const
{
  pb.buildAndRegisterGatherAndOrientationEvaluators(fm, lof, user_data);
}

template <typename EvalT>
void Example::BCStrategy_Interface_Robin<EvalT>::
postRegistrationSetup(typename panzer::Traits::SetupData /* d */,
		      PHX::FieldManager<panzer::Traits>& /* vm */)
{ 
}

template <typename EvalT>
void Example::BCStrategy_Interface_Robin<EvalT>::
evaluateFields(typename panzer::Traits::EvalData /* d */)
{  
}
