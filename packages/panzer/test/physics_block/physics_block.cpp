#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "Phalanx_FieldManager.hpp"

#include "Panzer_Traits.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_InputEquationSet.hpp"
#include "Panzer_PhysicsBlock.hpp"

#include "user_app_EquationSetFactory.hpp"
#include "user_app_ModelFactory.hpp"
#include "Panzer_ModelFactory_TemplateManager.hpp"
#include "user_app_ModelFactory_TemplateBuilder.hpp"

namespace panzer_test_utils {

  Teuchos::RCP<panzer::PhysicsBlock> createPhysicsBlock()
  {
    
    panzer::InputEquationSet ies_1;
    {
      ies_1.name = "Energy";
      ies_1.basis = "Q2";
      ies_1.integration_order = 1;
      ies_1.model_id = 6;
      ies_1.model_factory = "rf";
      ies_1.prefix = "";
      ies_1.params.set<int>("junk", 1);
    }
    
    panzer::InputEquationSet ies_2;
    {
      ies_2.name = "Energy";
      ies_2.basis = "Q1";
      ies_2.integration_order = 1;
      ies_2.model_id = 6;
      ies_2.model_factory = "rf";
      ies_2.prefix = "ION_";
      ies_2.params.set<int>("junk", 1);
    }
    
    std::size_t num_cells = 20;
    std::size_t base_cell_dimension = 3;
    panzer::CellData cd(num_cells, base_cell_dimension);

    panzer::InputPhysicsBlock ipb;
    {
      ipb.physics_block_id = "4";
      ipb.eq_sets.push_back(ies_1);
      ipb.eq_sets.push_back(ies_2);
    }
    
    user_app::MyFactory eqs_factory;
    
    std::string element_block_id = "eblock_id";
    Teuchos::RCP<panzer::PhysicsBlock> physics_block = 
      Teuchos::rcp(new panzer::PhysicsBlock(ipb,element_block_id,cd,eqs_factory));
    
    return physics_block;
  }

  Teuchos::RCP<panzer::ModelFactory_TemplateManager<panzer::Traits> >
  buildModelFactory() 
  {
    user_app::MyModelFactory_TemplateBuilder builder;

    Teuchos::RCP<panzer::ModelFactory_TemplateManager<panzer::Traits> > 
      model_factory = 
      Teuchos::rcp(new panzer::ModelFactory_TemplateManager<panzer::Traits>);
    
    model_factory->buildObjects(builder);

    return model_factory;
  }

  std::vector<Teuchos::ParameterList> buildModelDescriptors()
  {    
    Teuchos::ParameterList p;
    {
      p.set("Object Type","Constant");
      p.set("Name","UX");
      p.set("Value",4.5);
    }

    std::vector<Teuchos::ParameterList> evaluators_to_build;
    evaluators_to_build.push_back(p);
    evaluators_to_build.push_back(p);
    evaluators_to_build.push_back(p);

    return evaluators_to_build;
  }

}

namespace panzer {


  TEUCHOS_UNIT_TEST(physics_block, getDOFNames)
  {
    Teuchos::RCP<panzer::PhysicsBlock> physics_block = 
      panzer_test_utils::createPhysicsBlock();

    const std::vector<std::string>& dof_names = physics_block->getDOFNames();

    TEST_EQUALITY(dof_names.size(), 2);
    TEST_EQUALITY(dof_names[0], "TEMPERATURE");
    TEST_EQUALITY(dof_names[1], "ION_TEMPERATURE");
  }

  TEUCHOS_UNIT_TEST(physics_block, getProvidedDOFs)
  {
    Teuchos::RCP<panzer::PhysicsBlock> physics_block = 
      panzer_test_utils::createPhysicsBlock();

    const std::vector<panzer::PhysicsBlock::StrBasisPair>& basis = 
      physics_block->getProvidedDOFs();

    TEST_EQUALITY(basis.size(), 2);
    TEST_EQUALITY(basis[0].first, "TEMPERATURE");
    TEST_EQUALITY(basis[1].first, "ION_TEMPERATURE");
    TEST_EQUALITY(basis[0].second->name(), "Q2");
    TEST_EQUALITY(basis[1].second->name(), "Q1");
    TEST_EQUALITY(basis[0].second->getCardinality(), 27);
    TEST_EQUALITY(basis[1].second->getCardinality(), 8);
  }

  TEUCHOS_UNIT_TEST(physics_block, getBases)
  {
    Teuchos::RCP<panzer::PhysicsBlock> physics_block = 
      panzer_test_utils::createPhysicsBlock();

    const std::map<std::string,Teuchos::RCP<panzer::Basis> >& unique_basis = 
      physics_block->getBases();

    TEST_EQUALITY(unique_basis.size(), 2);
    TEST_ASSERT(unique_basis.find("Q2") != unique_basis.end());
    TEST_ASSERT(unique_basis.find("Q1") != unique_basis.end());
    TEST_EQUALITY(unique_basis.find("Q2")->second->getCardinality(), 27);
    TEST_EQUALITY(unique_basis.find("Q1")->second->getCardinality(), 8);
  }

  TEUCHOS_UNIT_TEST(physics_block, getBaseCellTopology)
  {
    Teuchos::RCP<panzer::PhysicsBlock> physics_block = 
      panzer_test_utils::createPhysicsBlock();

    TEST_EQUALITY(physics_block->getBaseCellTopology().getDimension(), 3);
  }

  TEUCHOS_UNIT_TEST(physics_block, physicsBlockID)
  {
    Teuchos::RCP<panzer::PhysicsBlock> physics_block = 
      panzer_test_utils::createPhysicsBlock();

    TEST_EQUALITY(physics_block->physicsBlockID(), "4");
  }

  TEUCHOS_UNIT_TEST(physics_block, getCellData)
  {
    Teuchos::RCP<panzer::PhysicsBlock> physics_block = 
      panzer_test_utils::createPhysicsBlock();

    TEST_EQUALITY(physics_block->cellData().numCells(), 20);
    TEST_EQUALITY(physics_block->cellData().isSide(), false);
  }

  TEUCHOS_UNIT_TEST(physics_block, nontemplate_evaluator_builders)
  {

    Teuchos::RCP<panzer::PhysicsBlock> physics_block = 
      panzer_test_utils::createPhysicsBlock();

    PHX::FieldManager<panzer::Traits> fm;

    physics_block->buildAndRegisterEquationSetEvaluators(fm);
    physics_block->buildAndRegisterGatherScatterEvaluators(fm);

    Teuchos::RCP<panzer::ModelFactory_TemplateManager<panzer::Traits> > factory =
      panzer_test_utils::buildModelFactory(); 
    std::map<std::string,Teuchos::RCP<panzer::ModelFactory_TemplateManager<panzer::Traits> > > my_factories;
    my_factories["rf"] = factory;

    std::vector<Teuchos::ParameterList> models = panzer_test_utils::buildModelDescriptors();

    physics_block->buildAndRegisterModelEvaluators(fm,my_factories,models);
  }

  TEUCHOS_UNIT_TEST(physics_block, elementBlockID)
  {

    Teuchos::RCP<panzer::PhysicsBlock> physics_block = 
      panzer_test_utils::createPhysicsBlock();
   

    TEST_EQUALITY(physics_block->elementBlockID(),"eblock_id");
  }

  TEUCHOS_UNIT_TEST(physics_block, templated_evaluator_builders)
  {

    Teuchos::RCP<panzer::PhysicsBlock> physics_block = 
      panzer_test_utils::createPhysicsBlock();

    PHX::FieldManager<panzer::Traits> fm;

    physics_block->buildAndRegisterEquationSetEvaluatorsForType<panzer::Traits::Residual>(fm);
    physics_block->buildAndRegisterGatherScatterEvaluatorsForType<panzer::Traits::Residual>(fm);
    physics_block->buildAndRegisterEquationSetEvaluatorsForType<panzer::Traits::Jacobian>(fm);
    physics_block->buildAndRegisterGatherScatterEvaluatorsForType<panzer::Traits::Jacobian>(fm);

    Teuchos::RCP<panzer::ModelFactory_TemplateManager<panzer::Traits> > factory =
      panzer_test_utils::buildModelFactory(); 
    std::map<std::string,Teuchos::RCP<panzer::ModelFactory_TemplateManager<panzer::Traits> > > my_factories;
    my_factories["rf"] = factory;

    std::vector<Teuchos::ParameterList> models = panzer_test_utils::buildModelDescriptors();

    physics_block->buildAndRegisterModelEvaluatorsForType<panzer::Traits::Residual>(fm, my_factories, models);
    physics_block->buildAndRegisterModelEvaluatorsForType<panzer::Traits::Jacobian>(fm, my_factories, models);
  }


}
