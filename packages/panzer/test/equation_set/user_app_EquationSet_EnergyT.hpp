#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Basis.hpp"

// include evaluators here
#include "Panzer_Integrator_BasisTimesScalar.hpp"
#include "Panzer_Integrator_GradBasisDotVector.hpp"
#include "Panzer_ScalarToVector.hpp"
#include "Panzer_Sum.hpp"
#include "Panzer_Constant.hpp"
#include "user_app_Convection.hpp"

// ***********************************************************************
template <typename EvalT>
user_app::EquationSet_Energy<EvalT>::
EquationSet_Energy(const panzer::InputEquationSet& ies,
		   const panzer::CellData& cell_data) :
  panzer::EquationSet_DefaultImpl<EvalT>(ies, cell_data)
{
  this->m_eqset_prefix = ies.prefix;

  // ********************
  // Assemble DOF names and Residual names
  // ********************
  this->m_dof_names = Teuchos::rcp(new std::vector<std::string>);
  this->m_dof_names->push_back(this->m_eqset_prefix+"TEMPERATURE");

  this->m_dof_gradient_names = Teuchos::rcp(new std::vector<std::string>);
  this->m_dof_gradient_names->push_back(this->m_eqset_prefix+"GRAD_TEMPERATURE");

  this->m_residual_names = Teuchos::rcp(new std::vector<std::string>);
  this->m_residual_names->push_back(this->m_eqset_prefix+"RESIDUAL_TEMPERATURE");

  this->m_scatter_name = "Scatter_"+this->m_eqset_prefix+"RESIDUAL_TEMPERATURE";

  // ********************
  // Build Basis Functions and Integration Rules
  // ********************
  
  this->m_int_rule = 
    Teuchos::rcp(new panzer::IntegrationRule(ies.integration_order,cell_data));
  
  this->m_basis = Teuchos::rcp(new panzer::Basis(ies.basis, 
						 *(this->m_int_rule)));

  this->m_eval_plist = Teuchos::rcp(new Teuchos::ParameterList);
  this->m_eval_plist->set("IR", this->m_int_rule);
  this->m_eval_plist->set("Basis", this->m_basis);
  this->m_eval_plist->set("Equation Dimension", cell_data.baseCellDimension());

  this->m_provided_dofs.push_back(std::make_pair((*(this->m_dof_names))[0],
						 this->m_basis));
}

// ***********************************************************************
template <typename EvalT>
void user_app::EquationSet_Energy<EvalT>::
buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				      const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > > & dofs) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  
  const std::string& prefix = this->m_eqset_prefix;

  // ********************
  // Energy Equation
  // ********************

  // Transient Operator
  {
    //! \todo Transients
  }

  // Diffusion Operator
  {
    double diffusion_coefficient = 1.0;

    ParameterList p("Diffusion Residual");
    p.set("Residual Name", prefix+"RESIDUAL_TEMPERATURE_DIFFUSION_OP");
    p.set("Flux Name", prefix+"GRAD_TEMPERATURE");
    p.set("Basis", this->m_basis);
    p.set("IR", this->m_int_rule);
    p.set("Multiplier", -diffusion_coefficient);
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::Integrator_GradBasisDotVector<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }
  
  // Convection Operator
  {
    // Combine scalar velocities into a velocity vector
    {
      ParameterList p("Velocity: ScalarToVector");
      RCP<std::vector<std::string> > scalar_names = rcp(new std::vector<std::string>);
      scalar_names->push_back(prefix+"UX");
      scalar_names->push_back(prefix+"UY");
      p.set<RCP<const std::vector<std::string> > >("Scalar Names", scalar_names);
      p.set("Vector Name", prefix+"U");
      p.set("Data Layout Scalar",this->m_int_rule->dl_scalar);
      p.set("Data Layout Vector",this->m_int_rule->dl_vector);

      RCP< PHX::Evaluator<panzer::Traits> > op = 
	rcp(new panzer::ScalarToVector<EvalT,panzer::Traits>(p));
      
      fm.template registerEvaluator<EvalT>(op);
    }
    
    // Constant value for convection 
    {
      // UX
      {
	ParameterList p("Constant UX");
	p.set("Value", 0.0);
	p.set("Name", prefix+"UX");
	p.set("Data Layout", this->m_int_rule->dl_scalar);
	
	RCP< PHX::Evaluator<panzer::Traits> > op = 
	  rcp(new panzer::Constant<EvalT,panzer::Traits>(p));
	
	fm.template registerEvaluator<EvalT>(op);
      }

      // UY
      {
	ParameterList p("Constant UY");
	p.set("Value", 0.0);
	p.set("Name", prefix+"UY");
	p.set("Data Layout", this->m_int_rule->dl_scalar);
	
	RCP< PHX::Evaluator<panzer::Traits> > op = 
	  rcp(new panzer::Constant<EvalT,panzer::Traits>(p));
	
	fm.template registerEvaluator<EvalT>(op);
      }
      
    }

    // Evaluator to assemble convection term
    {
      ParameterList p("Convection Operator");
      p.set("IR", this->m_int_rule);
      p.set("Operator Name", prefix+"TEMPERATURE_CONVECTION_OP");
      p.set("A Name", prefix+"U");
      p.set("Gradient Name", prefix+"GRAD_TEMPERATURE");
      p.set("Multiplier", 1.0);

      RCP< PHX::Evaluator<panzer::Traits> > op = 
	rcp(new user_app::Convection<EvalT,panzer::Traits>(p));
      
      fm.template registerEvaluator<EvalT>(op);
    }

    // Integration operator (could sum this into source for efficiency)
    {
      ParameterList p("Convection Residual");
      p.set("Residual Name",prefix+"RESIDUAL_TEMPERATURE_CONVECTION_OP");
      p.set("Value Name", prefix+"TEMPERATURE_CONVECTION_OP");
      p.set("Basis", this->m_basis);
      p.set("IR", this->m_int_rule);
      p.set("Multiplier", 1.0);
      
      RCP< PHX::Evaluator<panzer::Traits> > op = 
	rcp(new panzer::Integrator_BasisTimesScalar<EvalT,panzer::Traits>(p));
      
      fm.template registerEvaluator<EvalT>(op);
    }
  }

  // Source Operator
  {    
    {
      ParameterList p("Source Operator");
      p.set("Value", 1.0);
      p.set("Name", prefix+"SOURCE_TEMPERATURE");
      p.set("Data Layout", this->m_int_rule->dl_scalar);
      
      RCP< PHX::Evaluator<panzer::Traits> > op = 
	rcp(new panzer::Constant<EvalT,panzer::Traits>(p));
      
      fm.template registerEvaluator<EvalT>(op);
    }

    {
      ParameterList p("Source Residual");
      p.set("Residual Name", prefix+"RESIDUAL_TEMPERATURE_SOURCE_OP");
      p.set("Value Name", prefix+"SOURCE_TEMPERATURE");
      p.set("Basis", this->m_basis);
      p.set("IR", this->m_int_rule);
      p.set("Multiplier", -1.0);
      
      RCP< PHX::Evaluator<panzer::Traits> > op = 
	rcp(new panzer::Integrator_BasisTimesScalar<EvalT,panzer::Traits>(p));
      
      fm.template registerEvaluator<EvalT>(op);
    }

  }

  // Use a sum operator to form the overall residual for the equation
  // - this way we avoid loading each operator separately into the
  // global residual and Jacobian
  {
    ParameterList p;
    p.set("Sum Name", prefix+"RESIDUAL_TEMPERATURE");

    RCP<std::vector<std::string> > sum_names = 
      rcp(new std::vector<std::string>);
    sum_names->push_back(prefix+"RESIDUAL_TEMPERATURE_DIFFUSION_OP");
    sum_names->push_back(prefix+"RESIDUAL_TEMPERATURE_SOURCE_OP");
    sum_names->push_back(prefix+"RESIDUAL_TEMPERATURE_CONVECTION_OP");
    p.set("Values Names", sum_names);
    p.set("Data Layout", this->m_basis->functional);

    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::Sum<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }

}

// ***********************************************************************
