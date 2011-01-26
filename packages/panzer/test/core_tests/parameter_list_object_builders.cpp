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
    fluid_eqs0.set("Model ID", 6);
    fluid_eqs0.set("Model Factory", "rf");
    fluid_eqs0.set("Prefix", "ION_");
    fluid_eqs0.sublist("Options").set("Tau_C", "Codina");
    Teuchos::ParameterList& fluid_eqs1 = fluid.sublist("EQ 1");
    fluid_eqs1.set("Name", "Momentum");
    fluid_eqs1.set("Basis", "Q2");
    fluid_eqs1.set("Integration Order", 2);
    fluid_eqs1.set("Model ID", 6);
    fluid_eqs1.set("Model Factory", "rf");
    fluid_eqs1.set("Prefix", "ION_");
    fluid_eqs1.sublist("Options").set("Tau_M", "Codina");
    
    Teuchos::ParameterList& solid = p.sublist("solid");
    solid.set("Number of Equation Sets", 1);
    Teuchos::ParameterList& solid_eqs0 = solid.sublist("EQ 0");
    solid_eqs0.set("Name", "Energy");
    solid_eqs0.set("Basis", "Q1");
    solid_eqs0.set("Integration Order", 1);
    solid_eqs0.set("Model ID", 6);
    solid_eqs0.set("Model Factory", "rf");
    solid_eqs0.set("Prefix", "ION_");
    solid_eqs0.sublist("Options").set("junk", 1);

    std::vector<panzer::InputPhysicsBlock> ipb;

    panzer::buildInputPhysicsBlocks(ipb, p);
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
    bc_1.set<std::size_t>("ID", 0);
    bc_1.set("Type", "Dirichlet");
    bc_1.set("Sideset ID", "4");
    bc_1.set("Element Block ID", "fluid");
    bc_1.set("Equation Set Name", "UX");
    bc_1.set("Strategy", "Constant");
    bc_1.sublist("Data").set("Value",1.0);

    panzer::buildBCs(bcs, bc_params);
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
