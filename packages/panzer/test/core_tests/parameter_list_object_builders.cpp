#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "Panzer_ParameterList_ObjectBuilders.hpp"
#include "Panzer_InputPhysicsBlock.hpp"
#include "Panzer_BC.hpp"

namespace panzer {

  TEUCHOS_UNIT_TEST(object_builders, input_physics_block)
  {
    Teuchos::ParameterList p("Test");
    p.set("Physics Blocks", "fluid,solid");
    Teuchos::ParameterList& fluid = p.sublist("fluid");
    fluid.set("Number of Equation Sets", 2);
    Teuchos::ParameterList& fluid_eqs0 = fluid.sublist("EQ 0");
    fluid_eqs0.set("Name", "Continuity");
    fluid_eqs0.set("Basis", "Q1");
    fluid_eqs0.set("Integration Order", 1);
    fluid_eqs0.set("Model ID", "fluid");
    fluid_eqs0.set("Prefix", "ION_");
    fluid_eqs0.sublist("Options").set("Tau_C", "Codina");
    Teuchos::ParameterList& fluid_eqs1 = fluid.sublist("EQ 1");
    fluid_eqs1.set("Name", "Momentum");
    fluid_eqs1.set("Basis", "Q2");
    fluid_eqs1.set("Integration Order", 2);
    fluid_eqs1.set("Model ID", "fluid");
    fluid_eqs1.set("Prefix", "ION_");
    fluid_eqs1.sublist("Options").set("Tau_M", "Codina");
    
    Teuchos::ParameterList& solid = p.sublist("solid");
    solid.set("Number of Equation Sets", 1);
    Teuchos::ParameterList& solid_eqs0 = solid.sublist("EQ 0");
    solid_eqs0.set("Name", "Energy");
    solid_eqs0.set("Basis", "Q1");
    solid_eqs0.set("Integration Order", 1);
    solid_eqs0.set("Model ID", "solid");
    solid_eqs0.set("Prefix", "ION_");
    solid_eqs0.sublist("Options").set("junk", 1);

    std::map<std::string,panzer::InputPhysicsBlock> ipb;

    panzer::buildInputPhysicsBlocks(ipb, p);

    TEST_EQUALITY(ipb["fluid"].eq_sets.size(), 2);
    TEST_EQUALITY(ipb["solid"].eq_sets.size(), 1);
    // check non-existent case
    TEST_EQUALITY(ipb["garbage"].eq_sets.size(), 0);
    
    TEST_EQUALITY(ipb["fluid"].physics_block_id, "fluid");
    TEST_EQUALITY(ipb["solid"].physics_block_id, "solid");

    // Implicit assumption on ordering of equation sets.  Only needed
    // for testing.
    TEST_EQUALITY(ipb["fluid"].eq_sets[0].name, "Continuity");
    TEST_EQUALITY(ipb["fluid"].eq_sets[0].basis, "Q1");
  }

  TEUCHOS_UNIT_TEST(object_builders, bc)
  {
    Teuchos::ParameterList bc_params;

    std::vector<panzer::BC> bcs;
    bc_params.set<int>("Number of Boundary Conditions", 2);
    Teuchos::ParameterList& bc_0 = bc_params.sublist("BC 0");
    bc_0.set<std::size_t>("ID", 0);
    bc_0.set("Type", "Dirichlet");
    bc_0.set("Sideset ID", "4");
    bc_0.set("Element Block ID", "fluid");
    bc_0.set("Equation Set Name", "UX");
    bc_0.set("Strategy", "Constant");
    bc_0.sublist("Data").set("Value",1.0);
    Teuchos::ParameterList& bc_1 = bc_params.sublist("BC 1");
    bc_1.set<std::size_t>("ID", 1);
    bc_1.set("Type", "Neumann");
    bc_1.set("Sideset ID", "4");
    bc_1.set("Element Block ID", "fluid");
    bc_1.set("Equation Set Name", "UX");
    bc_1.set("Strategy", "Constant");
    bc_1.sublist("Data").set("Value",1.0);

    panzer::buildBCs(bcs, bc_params);

    TEST_EQUALITY(bcs.size(), 2);
    TEST_EQUALITY(bcs[0].bcID(), 0);
    TEST_EQUALITY(bcs[1].bcID(), 1);
    TEST_EQUALITY(bcs[0].bcType(), panzer::BCT_Dirichlet);
    TEST_EQUALITY(bcs[1].bcType(), panzer::BCT_Neumann);
  }

  TEUCHOS_UNIT_TEST(object_builders, block_id_to_physics_id)
  {
    Teuchos::ParameterList p;
    p.set("eblock-0_0", "fluid");
    p.set("eblock-0_1", "fluid");
    p.set("eblock-1_0", "solid");
    p.set("eblock-1_1", "solid");
    
    std::map<std::string,std::string> b_to_p;
    
    panzer::buildBlockIdToPhysicsIdMap(b_to_p, p);

    TEST_EQUALITY(b_to_p["eblock-0_0"], "fluid");
    TEST_EQUALITY(b_to_p["eblock-0_1"], "fluid");
    TEST_EQUALITY(b_to_p["eblock-1_0"], "solid");
    TEST_EQUALITY(b_to_p["eblock-1_1"], "solid");
  }

  TEUCHOS_UNIT_TEST(object_builders, StringTokenizerMultiple)
  {
    std::string names =  "fluid,solid";
    std::vector<std::string> tokens;
    panzer::StringTokenizer(tokens, names);
    
    TEST_EQUALITY(tokens.size(), 2);
    TEST_EQUALITY(tokens[0], "fluid");
    TEST_EQUALITY(tokens[1], "solid");
  }
  TEUCHOS_UNIT_TEST(object_builders, StringTokenizerSingle)
  {
    std::string names =  "fluid";
    std::vector<std::string> tokens;
    panzer::StringTokenizer(tokens, names);
    
    TEST_EQUALITY(tokens.size(), 1);
    TEST_EQUALITY(tokens[0], "fluid");
  }

}
