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
#include "Panzer_DOFGradient.hpp"
#include "Panzer_DotProduct.hpp"

#include "Phalanx_MDField.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

// ***********************************************************************
template <typename EvalT>
Example::BCStrategy_Interface_WeakDirichletMatch<EvalT>::
BCStrategy_Interface_WeakDirichletMatch(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data) :
  panzer::BCStrategy_Interface_DefaultImpl<EvalT>(bc, global_data)
{
  TEUCHOS_ASSERT(this->m_bc.strategy() == "Weak Dirichlet Match Interface");
}

// ***********************************************************************
template <typename EvalT>
void Example::BCStrategy_Interface_WeakDirichletMatch<EvalT>::
setup(const panzer::PhysicsBlock& side_pb,
      const Teuchos::ParameterList& /* user_data */)
{
  using Teuchos::RCP;
  using std::vector;
  using std::string;
  using std::pair;

  const int di = this->getDetailsIndex();

  // obtain the dof name
  const string dof_name = di == 0 ? this->m_bc.equationSetName() : this->m_bc.equationSetName2();
  other_dof_name = di == 1 ? this->m_bc.equationSetName() : this->m_bc.equationSetName2();

  // need the dof value to form the residual
  this->requireDOFGather(dof_name);

  // unique residual name
  const string residual_name = "Residual_" + this->m_bc.equationSetName();
  const string diff_name = "Difference";

  const std::map<int,RCP< panzer::IntegrationRule > >& ir = side_pb.getIntegrationRules();
  TEUCHOS_ASSERT(ir.size() == 1); 
  
  const int integration_order = ir.begin()->second->order();

  this->addResidualContribution(residual_name,dof_name,diff_name,integration_order,side_pb);
}

// ***********************************************************************
template <typename EvalT>
void Example::BCStrategy_Interface_WeakDirichletMatch<EvalT>::
setSumValues(Teuchos::ParameterList& p,
             const std::string value_name1, const double scalar1,
             const std::string value_name2, const double scalar2)
{
  std::vector<std::string> values_names(2);
  values_names[0] = value_name1;
  values_names[1] = value_name2;
  p.set< Teuchos::RCP<std::vector<std::string> > >(
    "Values Names", Teuchos::rcp(new std::vector<std::string>(values_names)));
  std::vector<double> scalars(2);
  scalars[0] = scalar1;
  scalars[1] = scalar2;
  p.set< Teuchos::RCP<const std::vector<double> > >(
    "Scalars", Teuchos::rcp(new std::vector<double>(scalars)));
}

template <typename EvalT>
void Example::BCStrategy_Interface_WeakDirichletMatch<EvalT>::
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

  string residual_name = std::get<0>(data[0]);
  string dof_name = std::get<1>(data[0]);
  string diff_name = std::get<2>(data[0]);

  RCP<panzer::IntegrationRule> ir = std::get<5>(data[0]);
  RCP<const panzer::FieldLayoutLibrary> fll = pb.getFieldLibrary()->buildFieldLayoutLibrary(*ir);
  RCP<panzer::BasisIRLayout> basis = fll->lookupLayout(dof_name);

  if (this->getDetailsIndex() == 0) {
    const std::string
      dof_grad_name = dof_name + "_gradient",
      cancel_natural_name = dof_name + "_cancel",
      my_normal_name = "My_Normal",
      sum_contributions_name = "Sum_Contributions";
    // Weak Dirichlet match.
    {
      { // Get values on my side.
        ParameterList p("My DOF");
        p.set("Name", dof_name);
        p.set("Basis", basis); 
        p.set("IR", ir);
        const RCP< PHX::Evaluator<panzer::Traits> >
          op = rcp(new panzer::DOF<EvalT,panzer::Traits>(p));
        this->template registerEvaluator<EvalT>(fm, op);
      }
      { // Other DOF - my DOF.
        ParameterList p("other DOF - my DOF");
        p.set("Sum Name", diff_name);
        setSumValues(p, other_dof_name, 1, dof_name, -1);
        p.set("Data Layout", ir->dl_scalar);
        const RCP< PHX::Evaluator<panzer::Traits> >
          op = rcp(new panzer::Sum<EvalT,panzer::Traits>(p));
        this->template registerEvaluator<EvalT>(fm, op);
      }
    }
    // Cancel my natural (Neumann) BC.
    {
      { // Normal.
        ParameterList p("My Side Normal");
        p.set("Name", my_normal_name);
        p.set("Side ID", pb.cellData().side());
        p.set("IR", ir);
        p.set("Normalize", true);
        const RCP< PHX::Evaluator<panzer::Traits> >
          op = rcp(new panzer::Normals<EvalT,panzer::Traits>(p));
        this->template registerEvaluator<EvalT>(fm, op);
      }
      { // Gradient.
        ParameterList p("My DOF gradient");
        p.set("Name", dof_name);
        p.set("Gradient Name", dof_grad_name);
        p.set("Basis", basis); 
        p.set("IR", ir);
        const RCP< PHX::Evaluator<panzer::Traits> >
          op = rcp(new panzer::DOFGradient<EvalT,panzer::Traits>(p));
        this->template registerEvaluator<EvalT>(fm, op);
      }
      { // dot(DOF gradient, normal).
        ParameterList p("dot(my DOF gradient, my normal)");
        p.set("Result Name", cancel_natural_name);
        p.set("Vector A Name", dof_grad_name);
        p.set("Vector B Name", my_normal_name);
        p.set("Point Rule", Teuchos::rcp_dynamic_cast<const panzer::PointRule>(ir));
        const RCP< PHX::Evaluator<panzer::Traits> >
          op = rcp(new panzer::DotProduct<EvalT,panzer::Traits>(p));
        this->template registerEvaluator<EvalT>(fm, op);
      }
    }
    // Add contributions to the residual.
    {
      { // Weak Dirichlet Match + Cancel Neumann
        ParameterList p("Weak Dirichlet Match + Cancel Neumann");
        p.set("Sum Name", sum_contributions_name);
        setSumValues(p, diff_name, 1e5, cancel_natural_name, -1);
        p.set("Data Layout", ir->dl_scalar);
        const RCP< PHX::Evaluator<panzer::Traits> >
          op = rcp(new panzer::Sum<EvalT,panzer::Traits>(p));
        this->template registerEvaluator<EvalT>(fm, op);
      }
      {
        using panzer::EvaluatorStyle;
        using panzer::Integrator_BasisTimesScalar;
        using panzer::Traits;
        using PHX::Evaluator;
        double multiplier(1);
        const RCP<Evaluator<Traits>> op = rcp(new
          Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::EVALUATES,
          residual_name, sum_contributions_name, *basis, *ir, multiplier));
        this->template registerEvaluator<EvalT>(fm, op);
      }
    }
  } else {
    { // Get values on other side.
      ParameterList p("Other DOF");
      p.set("Name", dof_name);
      p.set("Basis", basis); 
      p.set("IR", ir);
      const RCP< PHX::Evaluator<panzer::Traits> >
        op = rcp(new panzer::DOF<EvalT,panzer::Traits>(p));
      this->template registerEvaluator<EvalT>(fm, op);
    }
  }
}

// ***********************************************************************
template <typename EvalT>
void Example::BCStrategy_Interface_WeakDirichletMatch<EvalT>::
buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& fm,
			                       const panzer::PhysicsBlock& pb,
				               const panzer::LinearObjFactory<panzer::Traits> & lof,
				               const Teuchos::ParameterList& user_data) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::vector;
  using std::map;
  using std::string;
  using std::pair;

  // Gather
  pb.buildAndRegisterGatherAndOrientationEvaluators(fm,lof,user_data);
}

// ***********************************************************************
template <typename EvalT>
void Example::BCStrategy_Interface_WeakDirichletMatch<EvalT>::
postRegistrationSetup(typename panzer::Traits::SetupData /* d */,
		      PHX::FieldManager<panzer::Traits>& /* vm */)
{
  
}


// ***********************************************************************
template <typename EvalT>
void Example::BCStrategy_Interface_WeakDirichletMatch<EvalT>::
evaluateFields(typename panzer::Traits::EvalData /* d */)
{
  
}
