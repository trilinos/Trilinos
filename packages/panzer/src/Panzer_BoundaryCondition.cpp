#include "Panzer_BoundaryCondition.hpp"
#include "Teuchos_TestForException.hpp"
#include <sstream>

//=======================================================================
// "Constant" bcs
//=======================================================================
panzer::BoundaryCondition::BoundaryCondition(int bc_id,
					     BCType bc_type,
					     int sideset_id,
					     int element_block_id,
					     std::string dof_name,
					     double value) :
  m_bcf_type(panzer::BoundaryCondition::BCFT_Constant),
  m_constant_value(value)
{
  this->setup(bc_id, bc_type, sideset_id, element_block_id, dof_name);
  
  if (this->bcType() == BCT_Dirichlet)
    m_strategy = "Dirichlet: Constant";
  else if (this->bcType() == BCT_Neumann)
    m_strategy = "Neumann: Constant";
  else
    TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Error: A constant bc was input, but the bc_type was neither dirichlet or neumann. Please fix your input.  The offending boundary condition is:\n" << *this << std::endl);
  
}

//=======================================================================
// "Strategy" bcs - Single variable provider
//=======================================================================
panzer::BoundaryCondition::BoundaryCondition(int bc_id,
					     BCType bc_type,
					     int sideset_id,
					     int element_block_id,
					     std::string equation_set_name,
					     std::string strategy) :
  m_bcf_type(panzer::BoundaryCondition::BCFT_Strategy),
  m_strategy(strategy)
{
  this->setup(bc_id, bc_type, sideset_id, element_block_id, equation_set_name);
}

//=======================================================================
// "Method Function" bcs - method function supplies extra input data
//=======================================================================
panzer::BoundaryCondition::BoundaryCondition(int bc_id,
					     BCType bc_type,
					     int sideset_id,
					     int element_block_id,
					     std::string equation_set_name,
					     std::string strategy,
					     int method_function_id) :
  m_bcf_type(panzer::BoundaryCondition::BCFT_MethodFunction),
  m_strategy(strategy),
  m_method_function_id(method_function_id)
{
  this->setup(bc_id, bc_type, sideset_id, element_block_id, equation_set_name);
}

//=======================================================================
//=======================================================================
panzer::BoundaryCondition::~BoundaryCondition()
{ }

//=======================================================================
//=======================================================================
void panzer::BoundaryCondition::setup(int bc_id,
				      BCType bc_type,
				      int sideset_id,
				      int element_block_id,
				      std::string equation_set_name)
{
  m_bc_id = bc_id;
  m_bc_type = bc_type;
  m_sideset_id = sideset_id;
  m_element_block_id = element_block_id;
  m_equation_set_name = equation_set_name;
}

//=======================================================================
//=======================================================================
int panzer::BoundaryCondition::bcID() const
{
  return m_bc_id;
}

//=======================================================================
//=======================================================================
panzer::BCType panzer::BoundaryCondition::bcType() const
{
  return m_bc_type;
}

//=======================================================================
//=======================================================================
int panzer::BoundaryCondition::sidesetID() const
{
  return m_sideset_id;
}

//=======================================================================
//=======================================================================
int panzer::BoundaryCondition::elementBlockID() const
{
  return m_element_block_id;
}

//=======================================================================
//=======================================================================
std::string panzer::BoundaryCondition::equationSetName() const
{
  return m_equation_set_name;
}

//=======================================================================
//=======================================================================
double panzer::BoundaryCondition::constantValue() const
{
  return  m_constant_value;
}

//=======================================================================
//=======================================================================
std::string panzer::BoundaryCondition::strategy() const
{
  return m_strategy;
}

//=======================================================================
//=======================================================================
int panzer::BoundaryCondition::methodFunctionID() const
{
  return m_method_function_id;
}

//=======================================================================
//=======================================================================
bool panzer::BoundaryCondition::isConstant() const
{
  return m_bcf_type == panzer::BoundaryCondition::BCFT_Constant;
}

//=======================================================================
//=======================================================================
bool panzer::BoundaryCondition::isStrategy() const
{
  return m_bcf_type == panzer::BoundaryCondition::BCFT_Strategy;
}

//=======================================================================
//=======================================================================
bool panzer::BoundaryCondition::isMethodFunction() const
{
  return m_bcf_type == panzer::BoundaryCondition::BCFT_MethodFunction;
}

//=======================================================================
//=======================================================================
std::string panzer::BoundaryCondition::identifier() const
{
  std::ostringstream os;
  os << "BC(" <<  bcID() << ")";
  return os.str();
}

//=======================================================================
//=======================================================================
void panzer::BoundaryCondition::print(std::ostream& os) const
{
  using std::endl;

  os << "panzer::BoundaryCondition" << endl;

  os << "  BoundaryCondition ID =" << m_bc_id << endl;

  std::string type;
  if (m_bc_type == BCT_Dirichlet)
    type = "Dirichlet";
  else
    type = "Neumann";

  os << "  Type = " << type << endl;
  
  os << "  Side Set ID = " << m_sideset_id << endl;

  os << "  Element Block ID =" << m_element_block_id << endl;

  os << "  Function Type = ";
  if (m_bcf_type == panzer::BoundaryCondition::BCFT_Constant) {
    os << "Constant" << endl;
    os << "  Variable Name = " << m_equation_set_name << endl;
    os << "  Value = " << m_constant_value;
  }
  else if (m_bcf_type == panzer::BoundaryCondition::BCFT_Strategy) {
    os << "Strategy " << endl;
    os << "  Variable Name(s) = " << m_equation_set_name << endl;
    os << "  Strategy Name = " << m_strategy;
  }
  else if (m_bcf_type == panzer::BoundaryCondition::BCFT_MethodFunction) {
    os << "Method Function" << endl;
    os << "  Variable Name(s) = " << m_equation_set_name << endl;
    os << "  Strategy Name = " << m_strategy << endl;
    os << "  Method Function ID = " << m_method_function_id;
  }
}

//=======================================================================
//=======================================================================
std::ostream& 
panzer::operator<<(std::ostream & os, const panzer::BoundaryCondition& bc)
{
  bc.print(os);
  return os;
}
