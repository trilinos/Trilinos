#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Panzer_config.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_AssemblyEngine.hpp"
#include "Panzer_AssemblyEngine_TemplateManager.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Panzer_SGEpetraLinearObjFactory.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_PauseToAttach.hpp"
#include "user_app_EquationSetFactory.hpp"
#include "user_app_ClosureModel_Factory_TemplateBuilder.hpp"
#include "user_app_BCStrategy_Factory.hpp"

#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

#include "Stokhos_HermiteBasis.hpp"
#include "Stokhos_CompletePolynomialBasis.hpp"
#include "Stokhos_QuadOrthogPolyExpansion.hpp"
#include "Stokhos_TensorProductQuadrature.hpp"

#include "UnitTest_UniqueGlobalIndexer.hpp"

namespace panzer {

void testInitialzation(panzer::InputPhysicsBlock& ipb,
			 std::vector<panzer::BC>& bcs);

Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > buildExpansion(int numDim,int order)
{
   Teuchos::Array<Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(numDim);
   for(int i=0;i<numDim;i++)
      bases[i] = Teuchos::rcp(new Stokhos::HermiteBasis<int,double>(order));
   Teuchos::RCP<Stokhos::ProductBasis<int,double> > basis = Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));

   // build Cijk and "expansion"
   int kExpOrder = basis->size();
   // if(!fullExpansion)
   //    kExpOrder = numDim+1;
   Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > Cijk = basis->computeTripleProductTensor(kExpOrder);
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

  // panzer::pauseToAttach();

  RCP<Stokhos::OrthogPolyExpansion<int,double> > sgExpansion = buildExpansion(3,5);
  RCP<unit_test::UniqueGlobalIndexer> indexer
        = rcp(new unit_test::UniqueGlobalIndexer(myRank,numProc));

  panzer::InputPhysicsBlock ipb;
  std::vector<panzer::BC> bcs;
  testInitialzation(ipb, bcs);

  Teuchos::RCP<panzer::FieldManagerBuilder<short,int> > fmb = 
    Teuchos::rcp(new panzer::FieldManagerBuilder<short,int>);

  // build physics blocks
  //////////////////////////////////////////////////////////////
  const std::size_t workset_size = 20;
  user_app::MyFactory eqset_factory;
  user_app::BCFactory bc_factory;
  std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;

  {
    std::map<std::string,std::string> block_ids_to_physics_ids;
    block_ids_to_physics_ids["block_0"] = "test physics";
    
    std::map<std::string,panzer::InputPhysicsBlock> 
      physics_id_to_input_physics_blocks;
    physics_id_to_input_physics_blocks["test physics"] = ipb;

    fmb->buildPhysicsBlocks(block_ids_to_physics_ids,
                            physics_id_to_input_physics_blocks,
                            Teuchos::as<int>(2), workset_size,
                            eqset_factory,
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

  std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > > volume_worksets;
  volume_worksets["block_0"] = panzer::buildWorksets("block_0",cellIds,coords, ipb, workset_size,2);

  std::map<panzer::BC,Teuchos::RCP<std::map<unsigned,panzer::Workset> >,panzer::LessBC> bc_worksets;
  bc_worksets[bcs[myRank]] = panzer::buildBCWorkset(bcs[myRank],cellIds,sideIds,coords, ipb,2);

  // build DOF Manager
  /////////////////////////////////////////////////////////////

  Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,short> > eLinObjFactory
        = Teuchos::rcp(new panzer::EpetraLinearObjFactory<panzer::Traits,short>(eComm.getConst(),indexer));
  Teuchos::RCP<panzer::SGEpetraLinearObjFactory<panzer::Traits,short> > sgeLinObjFactory
        = Teuchos::rcp(new panzer::SGEpetraLinearObjFactory<panzer::Traits,short>(eLinObjFactory,sgExpansion));
  Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory = sgeLinObjFactory;

  // setup field manager build
  /////////////////////////////////////////////////////////////
 
  // Add in the application specific closure model factory
  Teuchos::RCP<const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> > cm_factory = 
    Teuchos::rcp(new panzer::ClosureModelFactory_TemplateManager<panzer::Traits>);
  user_app::MyModelFactory_TemplateBuilder cm_builder;
  (Teuchos::rcp_const_cast<panzer::ClosureModelFactory_TemplateManager<panzer::Traits> >(cm_factory))->buildObjects(cm_builder);

  Teuchos::ParameterList closure_models("Closure Models");
  closure_models.sublist("solid").sublist("SOURCE_TEMPERATURE").set<double>("Value",1.0);
  closure_models.sublist("ion solid").sublist("SOURCE_ION_TEMPERATURE").set<double>("Value",1.0);

  Teuchos::ParameterList user_data("User Data");

  fmb->setupVolumeFieldManagers(volume_worksets,physicsBlocks,*cm_factory,closure_models,indexer,*linObjFactory,user_data);

  // fmb->setupBCFieldManagers(bc_worksets,physicsBlocks,eqset_factory,*cm_factory,bc_factory,closure_models,*linObjFactory,user_data);

  panzer::AssemblyEngine_TemplateManager<panzer::Traits,short,int> ae_tm;
  panzer::AssemblyEngine_TemplateBuilder<short,int> builder(fmb,linObjFactory);
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
  e_global->A->Print(out);
  e_global->f->Print(out);

  out << "SG Solution" << std::endl;
  (*sgeGlobal->begin())->A->Print(out);
  (*sgeGlobal->begin())->f->Print(out);

  double diff,exact;

  // test residuals
  {
     (*sgeGlobal->begin())->f->Update(-1.0,*e_global->f,1.0);
     e_global->f->Norm2(&exact);
     (*sgeGlobal->begin())->f->Norm2(&diff);
     TEST_ASSERT(diff/exact < 1e-15);
  }

  // test jacobian
  {
     Epetra_Vector x(e_global->A->RowMap());
     Epetra_Vector y1(e_global->A->RowMap());
     Epetra_Vector y2(e_global->A->RowMap());

     // test a bunch of vectors
     for(int i=0;i<10;i++) {
        x.Random();
        y1.PutScalar(0);
        y2.PutScalar(0);

        e_global->A->Apply(x,y1);
        (*sgeGlobal->begin())->A->Apply(x,y2);

        y2.Update(-1.0,y1,1.0);
        y1.Norm2(&exact);
        y2.Norm2(&diff);

        TEST_ASSERT(diff/exact < 1e-15);
     }
  }
}

void testInitialzation(panzer::InputPhysicsBlock& ipb,
			 std::vector<panzer::BC>& bcs)
{
  panzer::InputEquationSet ies_1;
  ies_1.name = "Energy";
  ies_1.basis = "Q1";
  ies_1.integration_order = 1;
  ies_1.model_id = "solid";
  ies_1.prefix = "";

  panzer::InputEquationSet ies_2;
  ies_2.name = "Energy";
  ies_2.basis = "Q1";
  ies_2.integration_order = 1;
  ies_2.model_id = "ion solid";
  ies_2.prefix = "ION_";

  ipb.physics_block_id = "4";
  ipb.eq_sets.push_back(ies_1);
  ipb.eq_sets.push_back(ies_2);


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
