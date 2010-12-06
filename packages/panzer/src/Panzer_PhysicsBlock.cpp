#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"

#include "Phalanx_FieldManager.hpp"
#include "Phalanx_Evaluator_Factory.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_InputPhysicsBlock.hpp"
#include "Panzer_EquationSet_Factory.hpp"
#include "Shards_CellTopology.hpp"

// *******************************************************************
panzer::PhysicsBlock::
PhysicsBlock(const panzer::InputPhysicsBlock& ipb,
             const std::string & element_block_id,
	     const panzer::CellData cell_data,
	     const panzer::EquationSetFactory& factory) :
  m_physics_id(ipb.physics_block_id),
  m_element_block_id(element_block_id),
  m_cell_data(cell_data)
{
  using Teuchos::RCP;
  
  const std::vector<panzer::InputEquationSet>& input_eq_sets = ipb.eq_sets;

  TEST_FOR_EXCEPTION(input_eq_sets.size() < 1, std::runtime_error,
		     "There are no equation sets in the input file.  In order to use the phalanx assembly routines, you must add equation sets to a physics block!");

  m_equation_sets.clear();
  for (std::size_t i=0; i < input_eq_sets.size(); ++i) {
 
    RCP<panzer::EquationSet_TemplateManager<panzer::Traits> > eq_set
      = factory.buildEquationSet(input_eq_sets[i], m_cell_data);
 
    // add this equation set in
    m_equation_sets.push_back(eq_set);

    // figure out the DOF names from equation set
    const std::vector<std::string> & eqNames = eq_set->begin()->getDOFNames();
    for(std::size_t j=0;j<eqNames.size();j++)
       m_dof_names.push_back(eqNames[j]);

    // figure out the DOF names from equation set
    const std::vector<StrBasisPair> & sbNames = eq_set->begin()->getProvidedDOFs();
    for(std::size_t j=0;j<sbNames.size();j++)
       m_provided_dofs.push_back(sbNames[j]);

    // Get unique list of bases
    for(std::size_t j=0;j<sbNames.size();j++)
      m_bases[sbNames[j].second->name()] = sbNames[j].second;

  }
}

// *******************************************************************
void panzer::PhysicsBlock::
buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm) const
{
  using namespace std;
  using namespace panzer;
  using namespace Teuchos;

  // Loop over equation set template managers
  vector< RCP<EquationSet_TemplateManager<panzer::Traits> > >::const_iterator 
    eq_set = m_equation_sets.begin();
  for (;eq_set != m_equation_sets.end(); ++eq_set) {

    // Loop over evaluation types
    EquationSet_TemplateManager<panzer::Traits> eqstm = *(*eq_set);
    EquationSet_TemplateManager<panzer::Traits>::iterator eval_type =
      eqstm.begin();
    for (; eval_type != eqstm.end(); ++eval_type) {
      eval_type->buildAndRegisterEquationSetEvaluators(fm, m_provided_dofs);
    }
  }
}

// *******************************************************************
void panzer::PhysicsBlock::
buildAndRegisterGatherScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm) const
{
  using namespace std;
  using namespace panzer;
  using namespace Teuchos;

  // Loop over equation set template managers
  vector< RCP<EquationSet_TemplateManager<panzer::Traits> > >::const_iterator 
    eq_set = m_equation_sets.begin();
  for (;eq_set != m_equation_sets.end(); ++eq_set) {

    // Loop over evaluation types
    EquationSet_TemplateManager<panzer::Traits> eqstm = *(*eq_set);
    EquationSet_TemplateManager<panzer::Traits>::iterator eval_type =
      eqstm.begin();
    for (; eval_type != eqstm.end(); ++eval_type) {
      eval_type->buildAndRegisterGatherScatterEvaluators(fm, m_provided_dofs);
    }
  }
}

// *******************************************************************
void panzer::PhysicsBlock::
buildAndRegisterModelEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				const std::map<std::string,Teuchos::RCP<panzer::ModelFactory_TemplateManager<panzer::Traits> > >& factories,
				const std::vector<Teuchos::ParameterList>& models) const
{
  using namespace std;
  using namespace panzer;
  using namespace Teuchos;

  // Loop over equation set template managers
  vector< RCP<EquationSet_TemplateManager<panzer::Traits> > >::const_iterator 
    eq_set = m_equation_sets.begin();
  for (;eq_set != m_equation_sets.end(); ++eq_set) {

    // Loop over evaluation types
    EquationSet_TemplateManager<panzer::Traits> eqstm = *(*eq_set);
    EquationSet_TemplateManager<panzer::Traits>::iterator eval_type =
      eqstm.begin();
    for (; eval_type != eqstm.end(); ++eval_type) {
      eval_type->buildAndRegisterModelEvaluators(fm, m_provided_dofs, factories, models);
    }
  }
}

// *******************************************************************
const std::vector<std::string>& panzer::PhysicsBlock::getDOFNames() const
{
  return m_dof_names;
}

// *******************************************************************
const std::vector<panzer::PhysicsBlock::StrBasisPair>& panzer::PhysicsBlock::getProvidedDOFs() const
{
  return m_provided_dofs;
}

// *******************************************************************
const std::map<std::string,Teuchos::RCP<panzer::Basis> >& 
panzer::PhysicsBlock::getBases() const
{
  return m_bases;
}

// *******************************************************************
const shards::CellTopology panzer::PhysicsBlock::getBaseCellTopology() const
{
  TEST_FOR_EXCEPTION(m_bases.size() == 0, std::runtime_error,
		     "Cannot return a basis since none exist in this physics block.");
  return m_bases.begin()->second->getIntrepidBasis()->getBaseCellTopology();
}

// *******************************************************************
std::string panzer::PhysicsBlock::physicsBlockID() const
{
  return m_physics_id;
}

// *******************************************************************
std::string panzer::PhysicsBlock::elementBlockID() const
{
  return m_element_block_id;
}

// *******************************************************************
const panzer::CellData panzer::PhysicsBlock::cellData() const
{
  return m_cell_data;
}

// *******************************************************************
