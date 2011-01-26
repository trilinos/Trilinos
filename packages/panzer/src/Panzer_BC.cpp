#include "Panzer_BC.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_ParameterList.hpp"

//=======================================================================
//=======================================================================
panzer::BC::BC(std::size_t bc_id,
	       BCType bc_type,
	       std::string sideset_id,
	       std::string element_block_id,
	       std::string eq_set_name,
	       std::string strategy) :
  m_bc_id(bc_id),
  m_bc_type(bc_type),
  m_sideset_id(sideset_id),
  m_element_block_id(element_block_id),
  m_equation_set_name(eq_set_name),
  m_strategy(strategy)
{
}

//=======================================================================
//=======================================================================
panzer::BC::BC(std::size_t bc_id,
	       BCType bc_type,
	       std::string sideset_id,
	       std::string element_block_id,
	       std::string eq_set_name,
	       std::string strategy,
	       const Teuchos::ParameterList& p) :
  m_bc_id(bc_id),
  m_bc_type(bc_type),
  m_sideset_id(sideset_id),
  m_element_block_id(element_block_id),
  m_equation_set_name(eq_set_name),
  m_strategy(strategy)
{
  m_params = Teuchos::rcp(new Teuchos::ParameterList);
  *m_params = p;
}

//=======================================================================
//=======================================================================
panzer::BC::BC(const Teuchos::ParameterList& p)
{
  this->validateParameters(p);

  m_bc_id = p.get<std::size_t>("ID");
  std::string type = p.get<std::string>("Type");
  if (type == "Dirichlet")
    m_bc_type = BCT_Dirichlet;
  else if (type == "Neumann")
    m_bc_type = BCT_Neumann;

  m_sideset_id = p.get<std::string>("Sideset ID");
  m_element_block_id = p.get<std::string>("Element Block ID");
  m_equation_set_name = p.get<std::string>("Equation Set Name");
  m_strategy = p.get<std::string>("Strategy");
  m_params = Teuchos::rcp(new Teuchos::ParameterList);
  *m_params = p.sublist("Data");
}

//=======================================================================
//=======================================================================
panzer::BC::~BC()
{ }

//=======================================================================
//=======================================================================
std::size_t panzer::BC::bcID() const
{
  return m_bc_id;
}

//=======================================================================
//=======================================================================
panzer::BCType panzer::BC::bcType() const
{
  return m_bc_type;
}

//=======================================================================
//=======================================================================
std::string panzer::BC::sidesetID() const
{
  return m_sideset_id;
}

//=======================================================================
//=======================================================================
std::string panzer::BC::elementBlockID() const
{
  return m_element_block_id;
}

//=======================================================================
//=======================================================================
std::string panzer::BC::equationSetName() const
{
  return m_equation_set_name;
}

//=======================================================================
//=======================================================================
std::string panzer::BC::strategy() const
{
  return m_strategy;
}

//=======================================================================
//=======================================================================
Teuchos::RCP<const Teuchos::ParameterList> panzer::BC::params() const
{
  return m_params;
}

//=======================================================================
//=======================================================================
std::string panzer::BC::identifier() const
{
  std::ostringstream os;
  os << "BC(" <<  bcID() << ")";
  return os.str();
}

//=======================================================================
//=======================================================================
void panzer::BC::print(std::ostream& os) const
{
  using std::endl;

  os << "panzer::BC" << endl;

  os << "  BC ID =" << m_bc_id << endl;

  std::string type;
  if (m_bc_type == BCT_Dirichlet)
    type = "Dirichlet";
  else
    type = "Neumann";

  os << "  Type = " << type << endl;
  os << "  Side Set ID = " << m_sideset_id << endl;
  os << "  Element Block ID = " << m_element_block_id << endl;
  os << "  Strategy = " << m_strategy << endl;
  os << "  Variable Name(s) = " << m_equation_set_name << endl;
  os << "  Strategy Name = " << m_strategy;
  
  if (!Teuchos::is_null(m_params))
    os << endl << m_params;

}

//=======================================================================
//=======================================================================
void panzer::BC::validateParameters(const Teuchos::ParameterList& p) const
{
  Teuchos::ParameterList valid_params;
  
  valid_params.set<std::size_t>("ID", 0);
  valid_params.set<std::string>("Type", "Dirichlet");
  valid_params.set<std::string>("Sideset ID", "???");
  valid_params.set<std::string>("Element Block ID", "???");
  valid_params.set<std::string>("Equation Set Name", "???");
  valid_params.set<std::string>("Strategy", "???");
  valid_params.sublist("Data");

  p.validateParameters(valid_params, 0);
}

//=======================================================================
//=======================================================================
std::ostream& 
panzer::operator<<(std::ostream & os, const panzer::BC& bc)
{
  bc.print(os);
  return os;
}
