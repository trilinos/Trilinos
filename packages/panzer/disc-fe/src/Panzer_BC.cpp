// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_BC.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Panzer_WorksetDescriptor.hpp"

//=======================================================================
//=======================================================================
void 
panzer::buildBCs(std::vector<panzer::BC>& bcs,const Teuchos::ParameterList& p, const Teuchos::RCP<panzer::GlobalData> global_data)
{
  using Teuchos::ParameterList;
  
  bcs.clear();
  
  // Check for non-backward compatible change
  TEUCHOS_TEST_FOR_EXCEPTION(p.isParameter("Number of Boundary Conditions"),
			     std::logic_error,
			     "Error - the parameter \"Number of Boundary Conditions\" is no longer valid for the boundary condition sublist.  Please remove this from your input file!");
  
  std::size_t bc_index = 0;
  for (ParameterList::ConstIterator bc_pl=p.begin(); bc_pl != p.end(); ++bc_pl,++bc_index) {
    TEUCHOS_TEST_FOR_EXCEPTION( !(bc_pl->second.isList()), std::logic_error,
				"Error - All objects in the boundary condition sublist must be BC sublists!" );
    ParameterList& sublist = Teuchos::getValue<Teuchos::ParameterList>(bc_pl->second);
    
    panzer::BC bc(bc_index,sublist,global_data);
    bcs.push_back(bc);
  }
  
}

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
  m_params = Teuchos::rcp(new Teuchos::ParameterList(""));
  *m_params = p;
}

//=======================================================================
//=======================================================================
panzer::BC::BC(std::size_t bc_id,const Teuchos::ParameterList& p)
{
  Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::parameterList();
  *params = p;

  this->validateParameters(*params);

  m_bc_id = bc_id;
  std::string type = params->get<std::string>("Type");
  if (type == "Dirichlet")
    m_bc_type = BCT_Dirichlet;
  else if (type == "Neumann")
    m_bc_type = BCT_Neumann;
  else if (type == "Interface")
    m_bc_type = BCT_Interface;

  m_sideset_id = params->get<std::string>("Sideset ID");
  m_element_block_id = params->get<std::string>("Element Block ID");
  m_equation_set_name = params->get<std::string>("Equation Set Name");
  m_strategy = params->get<std::string>("Strategy");
  m_params = Teuchos::sublist(params,"Data");
  if (type == "Interface") {
     m_element_block_id2 = params->get<std::string>("Element Block ID2");
     m_equation_set_name2 = params->get<std::string>("Equation Set Name2");
  }
}

//=======================================================================
//=======================================================================
panzer::BC::BC(std::size_t bc_id,const Teuchos::ParameterList& p, const Teuchos::RCP<panzer::GlobalData> gd)
{
  Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::parameterList();
  *params = p;

  m_gd = gd;

  this->validateParameters(*params);

  m_bc_id = bc_id;
  std::string type = params->get<std::string>("Type");
  if (type == "Dirichlet")
    m_bc_type = BCT_Dirichlet;
  else if (type == "Neumann")
    m_bc_type = BCT_Neumann;
  else if (type == "Interface")
    m_bc_type = BCT_Interface;

  m_sideset_id = params->get<std::string>("Sideset ID");
  m_element_block_id = params->get<std::string>("Element Block ID");
  m_equation_set_name = params->get<std::string>("Equation Set Name");
  m_strategy = params->get<std::string>("Strategy");
  m_params = Teuchos::sublist(params,"Data");
  if (type == "Interface") {
    m_element_block_id2 = params->get<std::string>("Element Block ID2");
    m_equation_set_name2 = params->get<std::string>("Equation Set Name2");
  }
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
std::string panzer::BC::elementBlockID2() const
{
  return m_element_block_id2;
}

//=======================================================================
//=======================================================================
std::string panzer::BC::equationSetName() const
{
  return m_equation_set_name;
}

//=======================================================================
//=======================================================================
std::string panzer::BC::equationSetName2() const
{
  return m_equation_set_name2;
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
Teuchos::RCP<panzer::GlobalData> panzer::BC::global_data() const
{
  return m_gd;
}

//=======================================================================
//=======================================================================
Teuchos::RCP<Teuchos::ParameterList> 
panzer::BC::nonconstParams() const
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

  os << "  BC ID = " << m_bc_id << endl;

  std::string type;
  if (m_bc_type == BCT_Dirichlet)
    type = "Dirichlet";
  else if (m_bc_type == BCT_Neumann)
    type = "Neumann";
  else if (m_bc_type == BCT_Interface)
    type = "Interface";
  else
    type = "Neumann";

  os << "  Type = " << type << endl;
  os << "  Identifier = " << identifier() << endl;
  os << "  Side Set ID = " << m_sideset_id << endl;
  os << "  Element Block ID = " << m_element_block_id << endl;
  if (m_bc_type == BCT_Interface)
    os << "  Second Element Block ID = " << m_element_block_id2 << endl;
  os << "  Strategy = " << m_strategy << endl;
  os << "  Variable Name(s) = " << m_equation_set_name << endl;
  if (m_bc_type == BCT_Interface)
    os << "  Second Variable Name(s) = " << m_equation_set_name2 << endl;
  os << "  Strategy Name = " << m_strategy;
  
  if (!Teuchos::is_null(m_params))
    os << endl << m_params;

}

//=======================================================================
//=======================================================================
void panzer::BC::validateParameters(Teuchos::ParameterList& p) const
{
  Teuchos::ParameterList valid_params;
  
  valid_params.set<std::string>("Type", "Dirichlet");
  valid_params.set<std::string>("Sideset ID", "???");
  valid_params.set<std::string>("Element Block ID", "???");
  valid_params.set<std::string>("Element Block ID2", "???");
  valid_params.set<std::string>("Equation Set Name", "???");
  valid_params.set<std::string>("Equation Set Name2", "???");
  valid_params.set<std::string>("Strategy", "???");
  valid_params.sublist("Data").disableRecursiveValidation();

  p.validateParametersAndSetDefaults(valid_params);
}

//=======================================================================
//=======================================================================
std::ostream& 
panzer::operator<<(std::ostream & os, const panzer::BC& bc)
{
  bc.print(os);
  return os;
}

//=======================================================================
//=======================================================================

panzer::WorksetDescriptor 
panzer::bcDescriptor(const panzer::BC & bc)
{
  if(bc.bcType()==BCT_Interface) {
    WorksetDescriptor desc(bc.elementBlockID(),bc.elementBlockID2(),bc.sidesetID(),bc.sidesetID());
 
    return desc;
  }
  else {
    WorksetDescriptor desc(bc.elementBlockID(),bc.sidesetID());

    return desc;
  }
}
