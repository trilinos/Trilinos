// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "Phalanx_FieldManager.hpp"

#include "Panzer_Traits.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"

#include "user_app_EquationSetFactory.hpp"
#include "Panzer_ClosureModel_Factory_TemplateManager.hpp"
#include "user_app_ClosureModel_Factory_TemplateBuilder.hpp"
#include "user_app_ClosureModel_Factory.hpp"
#include "UnitTest_UniqueGlobalIndexer.hpp"

#include "Epetra_MpiComm.h"

namespace panzer_test_utils {

  Teuchos::RCP<panzer::PhysicsBlock> createPhysicsBlock()
  {
    
    Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::parameterList("4");

    {
      Teuchos::ParameterList& p = plist->sublist("a");
      p.set("Type","Energy");
      p.set("Prefix","");
      p.set("Model ID","solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",2);
    }

    {
      Teuchos::ParameterList& p = plist->sublist("b");
      p.set("Type","Energy");
      p.set("Prefix","ION_");
      p.set("Model ID","ion solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",1);
    }

    std::size_t num_cells = 20;
    const int default_int_order = 1;
    panzer::CellData cd(num_cells,Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Hexahedron<8> >())));

    Teuchos::RCP<user_app::MyFactory> eqs_factory = Teuchos::rcp(new user_app::MyFactory);
    
    std::string element_block_id = "eblock_id";

    Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();

    Teuchos::RCP<panzer::PhysicsBlock> physics_block = 
      Teuchos::rcp(new panzer::PhysicsBlock(plist,element_block_id,default_int_order,cd,eqs_factory,gd,false));
    
    return physics_block;
  }

  Teuchos::RCP<panzer::ClosureModelFactory_TemplateManager<panzer::Traits> >
  buildModelFactory() 
  {
    user_app::MyModelFactory_TemplateBuilder builder;

    Teuchos::RCP<panzer::ClosureModelFactory_TemplateManager<panzer::Traits> > 
      model_factory = 
      Teuchos::rcp(new panzer::ClosureModelFactory_TemplateManager<panzer::Traits>);
    
    model_factory->buildObjects(builder);

    return model_factory;
  }

  Teuchos::RCP<Teuchos::ParameterList> buildModelDescriptors()
  {    
    Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList("Closure Models"));
    {
      p->sublist("solid").sublist("SOURCE_TEMPERATURE").set<double>("Value",1.0);
      p->sublist("ion solid").sublist("SOURCE_ION_TEMPERATURE").set<double>("Value",1.0);
    }

    return p;
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

    const std::vector<panzer::StrPureBasisPair>& basis = 
      physics_block->getProvidedDOFs();

    TEST_EQUALITY(basis.size(), 2);
    TEST_EQUALITY(basis[0].first, "TEMPERATURE");
    TEST_EQUALITY(basis[1].first, "ION_TEMPERATURE");
    TEST_EQUALITY(basis[0].second->name(), "HGrad:2");
    TEST_EQUALITY(basis[1].second->name(), "HGrad:1");
    TEST_EQUALITY(basis[0].second->cardinality(), 27);
    TEST_EQUALITY(basis[1].second->cardinality(), 8);
  }

  TEUCHOS_UNIT_TEST(physics_block, getBases)
  {
    Teuchos::RCP<panzer::PhysicsBlock> physics_block = 
      panzer_test_utils::createPhysicsBlock();

    const std::map<std::string,Teuchos::RCP<panzer::PureBasis> >& unique_basis = 
      physics_block->getBases();

    TEST_EQUALITY(unique_basis.size(), 2);
    TEST_ASSERT(unique_basis.find("HGrad:2") != unique_basis.end());
    TEST_ASSERT(unique_basis.find("HGrad:1") != unique_basis.end());
    TEST_EQUALITY(unique_basis.find("HGrad:2")->second->cardinality(), 27);
    TEST_EQUALITY(unique_basis.find("HGrad:1")->second->cardinality(), 8);
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
    Teuchos::RCP<panzer::UniqueGlobalIndexer<short,int> > ugi 
          = Teuchos::rcp(new panzer::unit_test::UniqueGlobalIndexer(0,1));
    Teuchos::RCP<const Epetra_Comm> comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
    panzer::EpetraLinearObjFactory<panzer::Traits,short> elof(comm,ugi);

    Teuchos::RCP<panzer::PhysicsBlock> physics_block = 
      panzer_test_utils::createPhysicsBlock();

    PHX::FieldManager<panzer::Traits> fm;

    Teuchos::ParameterList user_data("User Data");

    physics_block->buildAndRegisterEquationSetEvaluators(fm, user_data);
    physics_block->buildAndRegisterGatherAndOrientationEvaluators(fm, elof, user_data);
    physics_block->buildAndRegisterDOFProjectionsToIPEvaluators(fm, user_data);
    physics_block->buildAndRegisterScatterEvaluators(fm, elof, user_data);

    Teuchos::RCP<panzer::ClosureModelFactory_TemplateManager<panzer::Traits> > factory =
      panzer_test_utils::buildModelFactory(); 

    Teuchos::RCP<Teuchos::ParameterList> models = panzer_test_utils::buildModelDescriptors();

    physics_block->buildAndRegisterClosureModelEvaluators(fm,*factory,*models, user_data);
  }

  TEUCHOS_UNIT_TEST(physics_block, elementBlockID)
  {

    Teuchos::RCP<panzer::PhysicsBlock> physics_block = 
      panzer_test_utils::createPhysicsBlock();
   

    TEST_EQUALITY(physics_block->elementBlockID(),"eblock_id");
  }

  TEUCHOS_UNIT_TEST(physics_block, templated_evaluator_builders)
  {
    Teuchos::RCP<panzer::UniqueGlobalIndexer<short,int> > ugi 
          = Teuchos::rcp(new panzer::unit_test::UniqueGlobalIndexer(0,1));
    Teuchos::RCP<const Epetra_Comm> comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
    panzer::EpetraLinearObjFactory<panzer::Traits,short> elof(comm,ugi);

    Teuchos::RCP<panzer::PhysicsBlock> physics_block = 
      panzer_test_utils::createPhysicsBlock();

    Teuchos::ParameterList user_data("User Data");

    PHX::FieldManager<panzer::Traits> fm;

    physics_block->buildAndRegisterEquationSetEvaluatorsForType<panzer::Traits::Residual>(fm, user_data);
    physics_block->buildAndRegisterEquationSetEvaluatorsForType<panzer::Traits::Jacobian>(fm, user_data);

    physics_block->buildAndRegisterGatherAndOrientationEvaluatorsForType<panzer::Traits::Residual>(fm, elof, user_data);
    physics_block->buildAndRegisterGatherAndOrientationEvaluatorsForType<panzer::Traits::Jacobian>(fm, elof, user_data);

    physics_block->buildAndRegisterDOFProjectionsToIPEvaluatorsForType<panzer::Traits::Residual>(fm, user_data);
    physics_block->buildAndRegisterDOFProjectionsToIPEvaluatorsForType<panzer::Traits::Jacobian>(fm, user_data);

    physics_block->buildAndRegisterScatterEvaluatorsForType<panzer::Traits::Residual>(fm, elof, user_data);
    physics_block->buildAndRegisterScatterEvaluatorsForType<panzer::Traits::Jacobian>(fm, elof, user_data);

    Teuchos::RCP<panzer::ClosureModelFactory_TemplateManager<panzer::Traits> > factory =
      panzer_test_utils::buildModelFactory(); 

    Teuchos::RCP<Teuchos::ParameterList> models = panzer_test_utils::buildModelDescriptors();

    physics_block->buildAndRegisterClosureModelEvaluatorsForType<panzer::Traits::Residual>(fm, *factory, *models, user_data);
    physics_block->buildAndRegisterClosureModelEvaluatorsForType<panzer::Traits::Jacobian>(fm, *factory, *models, user_data);
  }


}
