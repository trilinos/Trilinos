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

using Teuchos::RCP;
using Teuchos::rcp;

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Panzer_config.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_AssemblyEngine.hpp"
#include "Panzer_AssemblyEngine_TemplateManager.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Panzer_SGEpetraLinearObjFactory.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_PauseToAttach.hpp"
#include "Panzer_WorksetContainer.hpp"
#include "user_app_EquationSetFactory.hpp"
#include "user_app_ClosureModel_Factory_TemplateBuilder.hpp"
#include "user_app_BCStrategy_Factory.hpp"

#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

#include "Epetra_MpiComm.h"

#include "Stokhos_HermiteBasis.hpp"
#include "Stokhos_CompletePolynomialBasis.hpp"
#include "Stokhos_QuadOrthogPolyExpansion.hpp"
#include "Stokhos_TensorProductQuadrature.hpp"

#include "UnitTest_UniqueGlobalIndexer.hpp"

namespace panzer {

/** Pure virtual base class used to construct 
  * worksets on volumes and side sets.
  */
class TestWorksetFactory : public panzer::WorksetFactoryBase {
    mutable std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > > volume_worksets;
    mutable std::map<panzer::BC,Teuchos::RCP<std::map<unsigned,panzer::Workset> >,panzer::LessBC> bc_worksets;

public:
   TestWorksetFactory(const Teuchos::RCP<const shards::CellTopology> & topo,
                      std::vector<std::size_t> & cellIds,
                      std::vector<std::size_t> & sideIds,
                      Intrepid::FieldContainer<double> coords,
                      const panzer::PhysicsBlock& pb,
                      int workset_size,
                      std::vector<panzer::BC> & bcs,
                      int myRank) {
     volume_worksets["block_0"] = panzer::buildWorksets(pb,cellIds,coords,workset_size);
     bc_worksets[bcs[myRank]] = panzer::buildBCWorkset(bcs[myRank],pb,cellIds,sideIds,coords);
   }
   virtual ~TestWorksetFactory() {}

   
   Teuchos::RCP<std::vector<panzer::Workset> >
   getWorksets(const WorksetDescriptor & wd,
               const panzer::PhysicsBlock & pb) const
   { return volume_worksets[wd.getElementBlock()]; } // lazy

   Teuchos::RCP<std::vector<panzer::Workset> >
   getVolumeWorksets(const std::string & eBlock,
                     const panzer::PhysicsBlock & pb) const
   { return volume_worksets[eBlock]; }

   Teuchos::RCP<std::map<unsigned,panzer::Workset> > 
   getSideWorksets(const panzer::BC & bc,
		   const panzer::PhysicsBlock & pb) const
   { return bc_worksets[bc]; }
};


void testInitialzation(const Teuchos::RCP<Teuchos::ParameterList>& ipb,
		       std::vector<panzer::BC>& bcs);

Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > buildExpansion(int numDim,int order)
{
   Teuchos::Array<Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(numDim);
   for(int i=0;i<numDim;i++)
      bases[i] = Teuchos::rcp(new Stokhos::HermiteBasis<int,double>(order));
   Teuchos::RCP<Stokhos::ProductBasis<int,double> > basis = Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));

   Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > Cijk = basis->computeTripleProductTensor();
   Teuchos::RCP<Stokhos::Quadrature<int,double> > quadrature = Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));

   return Teuchos::rcp(new Stokhos::QuadOrthogPolyExpansion<int,double>(basis,Cijk,quadrature));
}

TEUCHOS_UNIT_TEST(field_manager_builder, basic)
{
  // build global (or serial communicator)
  #ifdef HAVE_MPI
     Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
  #else
     Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_SerialComm());
  #endif
 
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  int myRank = eComm->MyPID();
  int numProc = eComm->NumProc();

  Teuchos::RCP<panzer::GlobalData> global_data = panzer::createGlobalData();

  RCP<Stokhos::OrthogPolyExpansion<int,double> > sgExpansion = buildExpansion(3,5);
  RCP<unit_test::UniqueGlobalIndexer> indexer
        = rcp(new unit_test::UniqueGlobalIndexer(myRank,numProc));

  Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList("Physics Blocks");
  std::vector<panzer::BC> bcs;
  testInitialzation(ipb, bcs);

  Teuchos::RCP<panzer::FieldManagerBuilder> fmb = 
    Teuchos::rcp(new panzer::FieldManagerBuilder);

  // build physics blocks
  //////////////////////////////////////////////////////////////
  const std::size_t workset_size = 20;
  Teuchos::RCP<user_app::MyFactory> eqset_factory = Teuchos::rcp(new user_app::MyFactory);
  user_app::BCFactory bc_factory;
  std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;

  Teuchos::RCP<const shards::CellTopology> topo = 
      Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));
  {
    std::map<std::string,std::string> block_ids_to_physics_ids;
    block_ids_to_physics_ids["block_0"] = "test physics";

    std::map<std::string,Teuchos::RCP<const shards::CellTopology> > block_ids_to_cell_topo;
    block_ids_to_cell_topo["block_0"] = topo;
    
    int default_integration_order = 1;
      
    panzer::buildPhysicsBlocks(block_ids_to_physics_ids,
                               block_ids_to_cell_topo,
			       ipb,
			       default_integration_order,
			       workset_size,
                               eqset_factory,
                               global_data,
			       false,
                               physicsBlocks);
  }

  // build worksets
  //////////////////////////////////////////////////////////////
  Intrepid::FieldContainer<double> coords;
  std::vector<std::size_t> cellIds(1,0);
  std::vector<std::size_t> sideIds(1,0);
  sideIds[0] = (myRank==0) ? 0 : 1;
  indexer->getCoordinates(cellIds[0],coords);

   Teuchos::RCP<panzer::WorksetFactoryBase> wkstFactory
     = Teuchos::rcp(new TestWorksetFactory(topo,cellIds,sideIds,coords,*physicsBlocks[0],workset_size,bcs,myRank));
   Teuchos::RCP<panzer::WorksetContainer> wkstContainer
     = Teuchos::rcp(new panzer::WorksetContainer(wkstFactory,physicsBlocks,workset_size));

  // build DOF Manager
  /////////////////////////////////////////////////////////////

  Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,short> > eLinObjFactory
        = Teuchos::rcp(new panzer::EpetraLinearObjFactory<panzer::Traits,short>(eComm.getConst(),indexer));
  Teuchos::RCP<panzer::SGEpetraLinearObjFactory<panzer::Traits,short> > sgeLinObjFactory
        = Teuchos::rcp(new panzer::SGEpetraLinearObjFactory<panzer::Traits,short>(eLinObjFactory,sgExpansion,Teuchos::null));
  Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory = sgeLinObjFactory;

  // setup field manager build
  /////////////////////////////////////////////////////////////
 
  // Add in the application specific closure model factory
  user_app::MyModelFactory_TemplateBuilder cm_builder;
  panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory;
  cm_factory.buildObjects(cm_builder);

  Teuchos::ParameterList closure_models("Closure Models");
  closure_models.sublist("solid").sublist("SOURCE_TEMPERATURE").set<double>("Value",1.0);
  closure_models.sublist("ion solid").sublist("SOURCE_ION_TEMPERATURE").set<double>("Value",1.0);

  Teuchos::ParameterList user_data("User Data");

  fmb->setWorksetContainer(wkstContainer);
  fmb->setupVolumeFieldManagers(physicsBlocks,cm_factory,closure_models,*linObjFactory,user_data);
  fmb->setupBCFieldManagers(bcs,physicsBlocks,*eqset_factory,cm_factory,bc_factory,closure_models,*linObjFactory,user_data);

  panzer::AssemblyEngine_TemplateManager<panzer::Traits> ae_tm;
  panzer::AssemblyEngine_TemplateBuilder builder(fmb,linObjFactory);
  ae_tm.buildObjects(builder);

  // setup deterministic linear object containers
  RCP<panzer::LinearObjContainer> ghosted = eLinObjFactory->buildGhostedLinearObjContainer();
  RCP<panzer::LinearObjContainer> global = eLinObjFactory->buildLinearObjContainer();
  RCP<panzer::EpetraLinearObjContainer> e_ghosted = rcp_dynamic_cast<panzer::EpetraLinearObjContainer>(ghosted);
  RCP<panzer::EpetraLinearObjContainer> e_global = rcp_dynamic_cast<panzer::EpetraLinearObjContainer>(global);

  eLinObjFactory->initializeGhostedContainer(panzer::EpetraLinearObjContainer::X |
                                             panzer::EpetraLinearObjContainer::DxDt |
                                             panzer::EpetraLinearObjContainer::F |
                                             panzer::EpetraLinearObjContainer::Mat,*ghosted);
  eLinObjFactory->initializeContainer(panzer::EpetraLinearObjContainer::X |
                                      panzer::EpetraLinearObjContainer::DxDt |
                                      panzer::EpetraLinearObjContainer::F |
                                      panzer::EpetraLinearObjContainer::Mat,*global);

  panzer::AssemblyEngineInArgs input(ghosted,global);
  input.alpha = 0;
  input.beta = 1;

  // setup stochastic linear object containers
  RCP<panzer::LinearObjContainer> sg_ghosted = linObjFactory->buildGhostedLinearObjContainer();
  RCP<panzer::LinearObjContainer> sg_global = linObjFactory->buildLinearObjContainer();
  RCP<panzer::SGEpetraLinearObjContainer> sgeGhosted = rcp_dynamic_cast<panzer::SGEpetraLinearObjContainer>(sg_ghosted);
  RCP<panzer::SGEpetraLinearObjContainer> sgeGlobal = rcp_dynamic_cast<panzer::SGEpetraLinearObjContainer>(sg_global);

  sgeLinObjFactory->initializeGhostedContainer(panzer::EpetraLinearObjContainer::X |
                                             panzer::EpetraLinearObjContainer::DxDt |
                                             panzer::EpetraLinearObjContainer::F |
                                             panzer::EpetraLinearObjContainer::Mat,*sg_ghosted);
  sgeLinObjFactory->initializeContainer(panzer::EpetraLinearObjContainer::X |
                                      panzer::EpetraLinearObjContainer::DxDt |
                                      panzer::EpetraLinearObjContainer::F |
                                      panzer::EpetraLinearObjContainer::Mat,*sg_global);

  panzer::AssemblyEngineInArgs sg_input(sg_ghosted,sg_global);
  sg_input.alpha = 0;
  sg_input.beta = 1;

  // evaluate both deterministic and stochastic problem
  ae_tm.getAsObject<panzer::Traits::Jacobian>()->evaluate(input);
  ae_tm.getAsObject<panzer::Traits::SGJacobian>()->evaluate(sg_input);

  out << "Determ Solution" << std::endl;
  e_global->get_A()->Print(out);
  e_global->get_f()->Print(out);

  out << "SG Solution" << std::endl;
  (*sgeGlobal->begin())->get_A()->Print(out);
  (*sgeGlobal->begin())->get_f()->Print(out);

  double diff,exact;

  // test residuals
  {
     (*sgeGlobal->begin())->get_f()->Update(-1.0,*e_global->get_f(),1.0);
     e_global->get_f()->Norm2(&exact);
     (*sgeGlobal->begin())->get_f()->Norm2(&diff);
     TEST_ASSERT(diff/exact < 1e-15);
  }

  // test jacobian
  {
     Epetra_Vector x(e_global->get_A()->RowMap());
     Epetra_Vector y1(e_global->get_A()->RowMap());
     Epetra_Vector y2(e_global->get_A()->RowMap());

     // test a bunch of vectors
     for(int i=0;i<10;i++) {
        x.Random();
        y1.PutScalar(0);
        y2.PutScalar(0);

        e_global->get_A()->Apply(x,y1);
        (*sgeGlobal->begin())->get_A()->Apply(x,y2);

        y2.Update(-1.0,y1,1.0);
        y1.Norm2(&exact);
        y2.Norm2(&diff);

        TEST_ASSERT(diff/exact < 1e-15);
     }
  }
}
  void testInitialzation(const Teuchos::RCP<Teuchos::ParameterList>& ipb,
			 std::vector<panzer::BC>& bcs)
  {
    // Physics block
    Teuchos::ParameterList& physics_block = ipb->sublist("test physics");
    {
      Teuchos::ParameterList& p = physics_block.sublist("a");
      p.set("Type","Energy");
      p.set("Prefix","");
      p.set("Model ID","solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",1);
    }
    {
      Teuchos::ParameterList& p = physics_block.sublist("b");
      p.set("Type","Energy");
      p.set("Prefix","ION_");
      p.set("Model ID","ion solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",1);
    }

  {
    std::size_t bc_id = 0;
    panzer::BCType neumann = BCT_Dirichlet;
    std::string sideset_id = "left";
    std::string element_block_id = "block_0";
    std::string dof_name = "TEMPERATURE";
    std::string strategy = "Constant";
    double value = 5.0;
    Teuchos::ParameterList p;
    p.set("Value",value);
    panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name, 
		    strategy, p);
    bcs.push_back(bc);
  }    
  {
    std::size_t bc_id = 1;
    panzer::BCType neumann = BCT_Dirichlet;
    std::string sideset_id = "right";
    std::string element_block_id = "block_0";
    std::string dof_name = "ION_TEMPERATURE";
    std::string strategy = "Constant";
    double value = 5.0;
    Teuchos::ParameterList p;
    p.set("Value",value);
    panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name, 
		    strategy, p);
    bcs.push_back(bc);
  }   
}

}
