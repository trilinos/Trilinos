#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Basis.hpp"

// include evaluators here

// ***********************************************************************
template <typename EvalT>
user_app::EquationSet_Energy<EvalT>::
EquationSet_Energy(const panzer::InputEquationSet& ies,
		   const panzer::CellData& cell_data) :
  panzer::EquationSet_DefaultImpl<EvalT>(ies, cell_data)
{
  eqSetPrefix = ies.prefix;

  // ********************
  // Assemble DOF names
  // ********************
  m_dof_names = Teuchos::rcp(new std::vector<std::string>);
  m_dof_names->push_back(eqSetPrefix+"TEMPERATURE");


  // ********************
  // Build Basis Functions and Integration Rules
  // ********************
  
  m_int_rule = Teuchos::rcp(new panzer::IntegrationRule(ies.integration_order,
							cell_data));
  
  m_basis = Teuchos::rcp(new panzer::Basis(ies.basis, *m_int_rule));

  m_eval_plist = Teuchos::rcp(new Teuchos::ParameterList);
  m_eval_plist->set("IR", m_int_rule);
  m_eval_plist->set("Basis", m_basis);
  m_eval_plist->set("Equation Dimension", cell_data.baseCellDimension());

  m_provided_dofs.push_back(std::make_pair((*m_dof_names)[0],m_basis));
}

// ***********************************************************************
template <typename EvalT>
void user_app::EquationSet_Energy<EvalT>::
buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				      const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > > & dofs) const
{

  /*

  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  
  // ********************
  // Energy Equation
  // ********************

  // Transient Operator
  bool do_transient = false;
  if (do_transient) {
    ParameterList p("Transient Term");
    p.set("Residual Name", m_names.res.T + m_names.op.src);
    p.set("Value Name", m_names.dof.src_prefix+m_names.dof.T);
    //p.set("Data Layout Basis", m_node_scalar);
    //p.set("Data Layout IP", m_qp_scalar);
    p.set("Multiplier", -1.0);
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::Integrator_BasisTimesScalar<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }

  // Diffusion Operator
  {
    ParameterList p("Diffusion Residual");
    p.set("Residual Name", m_names.res.T + m_names.op.diff);
    p.set("Flux Name", m_names.field.q);
    p.set("Basis", m_basis);
    p.set("IR", m_int_rule);
    p.set("Multiplier", -1.0);
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::Integrator_GradBasisDotVector<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }

  // Convection Operator
  bool do_convection = true;
  if (do_convection) {
    // Combine scalar velocities into a velocity vector
    {
      ParameterList p("Velocity: ScalarToVector");
      RCP<const std::vector<std::string> > scalar_names = rcp(&(m_names.dof.u),false);
      p.set("Scalar Names", scalar_names);
      p.set("Vector Name", m_names.field.u);
      p.set("Data Layout Scalar",m_int_rule->dl_scalar);
      p.set("Data Layout Vector",m_int_rule->dl_vector);

      RCP< PHX::Evaluator<panzer::Traits> > op = 
	rcp(new panzer::ScalarToVector<EvalT,panzer::Traits>(p));
      
      fm.template registerEvaluator<EvalT>(op);
    }

    // Evaluator to assemble convection term
    {
      ParameterList p("Convection Operator");
      p.set("IR", m_int_rule);
      p.set("Operator Name", m_names.dof.T + m_names.op.conv);
      p.set("A Name", m_names.field.u);
      p.set("Gradient Name", m_names.dof.grad_prefix + m_names.dof.T);
      p.set("Multiplier", 1.0);

      RCP< PHX::Evaluator<panzer::Traits> > op = 
	rcp(new panzer::rf::Convection<EvalT,panzer::Traits>(p));
      
      fm.template registerEvaluator<EvalT>(op);
    }

    // Integration operator (could sum this into source for efficiency)
    {
      ParameterList p("Convection Residual");
      p.set("Residual Name", m_names.res.T + m_names.op.conv);
      p.set("Value Name", m_names.dof.T + m_names.op.conv);
      p.set("Basis", m_basis);
      p.set("IR", m_int_rule);
      p.set("Multiplier", 1.0);
      
      RCP< PHX::Evaluator<panzer::Traits> > op = 
	rcp(new panzer::Integrator_BasisTimesScalar<EvalT,panzer::Traits>(p));
      
      fm.template registerEvaluator<EvalT>(op);
    }
  }

  // Source Operator
  {
    ParameterList p("Source Residual");
    p.set("Residual Name", m_names.res.T + m_names.op.src);
    p.set("Value Name", m_names.dof.src_prefix+m_names.dof.T);
    p.set("Basis", m_basis);
    p.set("IR", m_int_rule);
    p.set("Multiplier", -1.0);
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::Integrator_BasisTimesScalar<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }

  // Check for magnetics equation and of found, add joule heating to
  // energe equation
  bool add_joule_heating = false;
  { 
    std::string bx_name = m_names.dof.B[0];
    typedef typename std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > >::const_iterator dof_it_type;
    for (dof_it_type i = dofs.begin(); i != dofs.end(); ++i) {
      if (i->first == bx_name)
	add_joule_heating = true;
    }
  }
  if (add_joule_heating) {     

    // Evaluator for lorentz force
    {      
      ParameterList p("Joule Heating Evaluator");
      const RCP<const panzer::rf::NameConventions> names_ptr = 
	rcp(&m_names, false);
      p.set("Names", names_ptr);
      p.set("IR", m_int_rule);

      RCP< PHX::Evaluator<panzer::Traits> > op = 
	rcp(new panzer::rf::JouleHeating<EvalT,panzer::Traits>(p));
      
      fm.template registerEvaluator<EvalT>(op);
    }

    {
      ParameterList p("Source Residual");
      p.set("Residual Name", m_names.res.T + m_names.op.joule_heating);
      p.set("Value Name", m_names.field.joule_heating);
      p.set("Basis", m_basis);
      p.set("IR", m_int_rule);
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
    p.set("Sum Name", m_names.res.T);

    RCP<std::vector<std::string> > sum_names = 
      rcp(new std::vector<std::string>);
    sum_names->push_back(m_names.res.T + m_names.op.diff);
    sum_names->push_back(m_names.res.T + m_names.op.src);
    if (do_convection) 
      sum_names->push_back(m_names.res.T + m_names.op.conv);
    if (add_joule_heating) 
      sum_names->push_back(m_names.res.T + m_names.op.joule_heating);
    p.set("Values Names", sum_names);
    p.set("Data Layout", m_basis->functional);

    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::Sum<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }

  // DOFs: Scalar value @ basis --> Scalar value @ IP 
  for (std::size_t i = 0; i < m_dof_names->size(); ++i) {
    std::string grad_name = m_names.dof.grad_prefix + (*m_dof_names)[i];
    ParameterList p;
    p.set("Name", (*m_dof_names)[i]);
    p.set("Basis", m_basis); 
    p.set("IR", m_int_rule);

    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::DOF<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }

  // Gradients of DOFs: Scalar value @ basis --> Vector value @ IP 
  for (std::size_t i = 0; i < m_dof_names->size(); ++i) {
    std::string grad_name = m_names.dof.grad_prefix + (*m_dof_names)[i];
    ParameterList p;
    p.set("Name", (*m_dof_names)[i]);
    p.set("Gradient Name", grad_name);
    p.set("Basis", m_basis); 
    p.set("IR", m_int_rule);

    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::DOFGradient<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }

  */

}

// ***********************************************************************
template <typename EvalT>
void user_app::EquationSet_Energy<EvalT>::
buildAndRegisterGatherScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > > & dofs) const
{

  /*

  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // ********************
  // DOFs (unknowns)
  // ********************

  // Gather
  {
    ParameterList p("Gather");
    p.set("Basis", m_basis);
    p.set("DOF Names", m_dof_names);
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::rf::GatherSolution<EvalT,panzer::Traits>(p));
    
    fm.template registerEvaluator<EvalT>(op);
  }
  
  // Scatter
  RCP<std::map<std::string,std::string> > names_map =
    rcp(new std::map<std::string,std::string>);
  RCP< std::vector<std::string> > residual_names = 
    rcp(new std::vector<std::string>);
  for (std::vector<std::string>::iterator i=m_dof_names->begin();
       i != m_dof_names->end(); ++i) {
    residual_names->push_back(m_names.res.residual_prefix + *i);
    names_map->insert(std::make_pair(m_names.res.residual_prefix + *i,*i));
  }

  {
    ParameterList p("Scatter");
    p.set("Scatter Name", eqSetPrefix+"Scatter_EnergyRawEvaluators");
    //p.set("Dummy Data Layout", m_dummy);
    p.set("Basis", m_basis);
    p.set("Evaluated Names", residual_names);
    p.set("Evaluated Map", names_map);

    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::ScatterResidual<EvalT,panzer::Traits>(p));
    
    fm.template registerEvaluator<EvalT>(op);
  }

  // Require variables
  {
    std::string reqField = eqSetPrefix+"Scatter_EnergyRawEvaluators";

    PHX::Tag<typename EvalT::ScalarT> tag(reqField, 
			        Teuchos::rcp(new PHX::MDALayout<Dummy>(0)));
    fm.template requireField<EvalT>(tag);
  }
 
  */
 
}

// ***********************************************************************
template <typename EvalT>
const Teuchos::RCP<Teuchos::ParameterList> 
user_app::EquationSet_Energy<EvalT>::
getEvaluatorParameterList() const
{
  return m_eval_plist;
}

// ***********************************************************************
template <typename EvalT>
const std::vector<std::string> & 
user_app::EquationSet_Energy<EvalT>::
getDOFNames() const
{
   return *m_dof_names;
}
// ***********************************************************************
template <typename EvalT>
const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > > & 
user_app::EquationSet_Energy<EvalT>::
getProvidedDOFs() const
{
   return m_provided_dofs;
}
// ***********************************************************************
