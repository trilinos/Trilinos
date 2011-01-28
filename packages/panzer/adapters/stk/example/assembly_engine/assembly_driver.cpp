#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_FancyOStream.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

#include "Thyra_get_Epetra_Operator.hpp"

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_AssemblyEngine.hpp"
#include "Panzer_AssemblyEngine_TemplateManager.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "user_app_EquationSetFactory.hpp"
#include "user_app_ModelFactory_TemplateBuilder.hpp"
#include "user_app_BCStrategy_Factory.hpp"

#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_MultiVectorOut.h"

#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_EpetraLinearOp.hpp"

#include "write_solution_data.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

void pause_to_attach()
{
   MPI_Comm mpicomm = MPI_COMM_WORLD;
   Teuchos::RCP<Teuchos::Comm<int> > comm = Teuchos::createMpiComm<int>(
         Teuchos::rcp(new Teuchos::OpaqueWrapper<MPI_Comm>(mpicomm)));
   Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
   out.setShowProcRank(true);
   out.setOutputToRootOnly(-1);

   out << "PID = " << getpid();

   if (comm->getRank() == 0)
      getchar();
   comm->barrier();
}


void testInitialzation(panzer::InputPhysicsBlock& ipb,
		       std::vector<panzer::BC>& bcs);

// calls MPI_Init and MPI_Finalize
int main(int argc,char * argv[])
{
   using Teuchos::RCP;
   using panzer::StrBasisPair;
   using panzer::StrBasisComp;

   Teuchos::GlobalMPISession mpiSession(&argc,&argv);
   RCP<Epetra_Comm> Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
   Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
   out.setOutputToRootOnly(0);
   out.setShowProcRank(true);

   // pause_to_attach();

   // variable declarations
   ////////////////////////////////////////////////////

   // factory definitions
   user_app::MyFactory eqset_factory;
   panzer_stk::SquareQuadMeshFactory mesh_factory;
   user_app::BCFactory bc_factory;

   // other declarations
   const std::size_t workset_size = 20;
   Teuchos::RCP<panzer::FieldManagerBuilder<int,int> > fmb = 
         Teuchos::rcp(new panzer::FieldManagerBuilder<int,int>);

   RCP<panzer_stk::STK_Interface> mesh;
   int base_cell_dimension = -1;

   // construction of uncommitted (no elements) mesh 
   ////////////////////////////////////////////////////////

   // set mesh factory parameters
   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",2);
   pl->set("Y Blocks",1);
   pl->set("X Elements",10);
   pl->set("Y Elements",10);
   mesh_factory.setParameterList(pl);
   mesh = mesh_factory.buildUncommitedMesh(MPI_COMM_WORLD);
   base_cell_dimension = Teuchos::as<int>(mesh->getDimension());

   // construct input physics and physics block
   ////////////////////////////////////////////////////////

   out << "BUILD PHYSICS" << std::endl;
   panzer::InputPhysicsBlock ipb;
   std::vector<panzer::BC> bcs;
   std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;
   {
      std::map<std::string,panzer::InputPhysicsBlock> physics_id_to_input_physics_blocks;
      std::map<std::string,std::string> block_ids_to_physics_ids;

      testInitialzation(ipb, bcs);

      block_ids_to_physics_ids["eblock-0_0"] = "test physics";
      block_ids_to_physics_ids["eblock-1_0"] = "test physics";
   
      physics_id_to_input_physics_blocks["test physics"] = ipb; // copying

      // build physicsBlocks map
      fmb->buildPhysicsBlocks(block_ids_to_physics_ids,
                             physics_id_to_input_physics_blocks,
                             base_cell_dimension, workset_size,
	                     eqset_factory,
                             physicsBlocks);
   }

   // finish building mesh, set required field variables and mesh bulk data
   ////////////////////////////////////////////////////////////////////////

   {
      std::vector<Teuchos::RCP<panzer::PhysicsBlock> >::const_iterator physIter;
      for(physIter=physicsBlocks.begin();physIter!=physicsBlocks.end();++physIter) {
         Teuchos::RCP<const panzer::PhysicsBlock> pb = *physIter;
         const std::vector<StrBasisPair> & blockFields = pb->getProvidedDOFs();

         // insert all fields into a set
         std::set<StrBasisPair,StrBasisComp> fieldNames;
         fieldNames.insert(blockFields.begin(),blockFields.end());

         // add basis to DOF manager: block specific
         std::set<StrBasisPair>::const_iterator fieldItr;
         for (fieldItr=fieldNames.begin();fieldItr!=fieldNames.end();++fieldItr) {

            mesh->addSolutionField(fieldItr->first,pb->elementBlockID());
         }
      }

      mesh_factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD);
   }


   // build worksets
   ////////////////////////////////////////////////////////

   // build worksets
   out << "BUILD WORKSETS" << std::endl;
   std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > > 
         volume_worksets = panzer_stk::buildWorksets(*mesh, ipb, workset_size);

   out << "block count = " << volume_worksets.size() << std::endl;
   out << "workset count = " << volume_worksets["eblock-0_0"]->size() << std::endl;
   
   out << "BUILD BC WORKSETS" << std::endl;
   std::map<panzer::BC,Teuchos::RCP<std::map<unsigned,panzer::Workset> >,panzer::LessBC> 
         bc_worksets = panzer_stk::buildBCWorksets(*mesh,ipb,bcs);

   // build DOF Manager
   /////////////////////////////////////////////////////////////
 
   out << "BUILD CONN MANAGER" << std::endl;
   // build the connection manager 
   const Teuchos::RCP<panzer::ConnManager<int,int> > 
     conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

   panzer::DOFManagerFactory<int,int> globalIndexerFactory;
   RCP<panzer::UniqueGlobalIndexer<int,int> > dofManager 
         = globalIndexerFactory.buildUniqueGlobalIndexer(MPI_COMM_WORLD,physicsBlocks,conn_manager);

   // construct some linear algebra object, build object to pass to evaluators
   Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory
         = Teuchos::rcp(new panzer::EpetraLinearObjFactory<panzer::Traits,int>(Comm,dofManager));

   // setup field manager build
   /////////////////////////////////////////////////////////////
 
   out << "SETUP FMB" << std::endl;
   fmb->setupVolumeFieldManagers(volume_worksets,physicsBlocks,dofManager,*linObjFactory);
   fmb->setupBCFieldManagers(bc_worksets,physicsBlocks,eqset_factory,bc_factory,*linObjFactory);

   // setup assembly engine
   /////////////////////////////////////////////////////////////

   // build assembly engine
   panzer::AssemblyEngine_TemplateManager<panzer::Traits,int,int> ae_tm;
   panzer::AssemblyEngine_TemplateBuilder<int,int> builder(fmb,linObjFactory);
   ae_tm.buildObjects(builder);

   // setup linear algebra and solve 
   /////////////////////////////////////////////////////////////

   // build ghosted variables
   out << "BUILD LA" << std::endl;
   RCP<panzer::EpetraLinearObjContainer> ghostCont 
         = Teuchos::rcp_dynamic_cast<panzer::EpetraLinearObjContainer>(linObjFactory->buildGhostedLinearObjContainer());
   RCP<panzer::EpetraLinearObjContainer> container 
         = Teuchos::rcp_dynamic_cast<panzer::EpetraLinearObjContainer>(linObjFactory->buildLinearObjContainer());

   panzer::AssemblyEngineInArgs input(ghostCont,container);

   // evaluate physics
   out << "EVALUTE" << std::endl;
   ae_tm.getAsObject<panzer::Traits::Residual>()->evaluate(input);
   ae_tm.getAsObject<panzer::Traits::Jacobian>()->evaluate(input);

   out << "RAN SUCCESSFULLY!" << std::endl;

   out << "SOLVE" << std::endl;

   // redistribute vectors and matrices so that we can parallel solve
   linObjFactory->ghostToGlobalContainer(*ghostCont,*container);

   Teuchos::RCP<const Thyra::LinearOpBase<double> > th_A = Thyra::epetraLinearOp(container->A);
   Teuchos::RCP<const Thyra::VectorSpaceBase<double> > range  = th_A->range();
   Teuchos::RCP<const Thyra::VectorSpaceBase<double> > domain = th_A->domain();

   Teuchos::RCP<Thyra::VectorBase<double> > th_x = Thyra::create_Vector(container->x,domain);
   Teuchos::RCP<Thyra::VectorBase<double> > th_f = Thyra::create_Vector(container->f,range);

   // solve with amesos
   Stratimikos::DefaultLinearSolverBuilder solverBuilder;
   Teuchos::RCP<Teuchos::ParameterList> validList = Teuchos::rcp(new Teuchos::ParameterList(*solverBuilder.getValidParameters()));
   solverBuilder.setParameterList(validList);
   
   RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory = solverBuilder.createLinearSolveStrategy("Amesos");
   RCP<Thyra::LinearOpWithSolveBase<double> > lows = Thyra::linearOpWithSolve(*lowsFactory, th_A.getConst());
   Thyra::solve<double>(*lows, Thyra::NOTRANS, *th_f, th_x.ptr());

   {
      EpetraExt::RowMatrixToMatrixMarketFile("a_op.mm",*container->A);
      EpetraExt::VectorToMatrixMarketFile("x_vec.mm",*container->x);
      EpetraExt::VectorToMatrixMarketFile("b_vec.mm",*container->f);
   }

   // redistribute solution vector
   linObjFactory->globalToGhostContainer(*container,*ghostCont);

   out << "WRITE" << std::endl;

   write_solution_data(*Teuchos::rcp_dynamic_cast<panzer::DOFManager<int,int> >(dofManager),*mesh,*ghostCont->x);
   mesh->writeToExodus("output.exo");

   return 0;
}

void testInitialzation(panzer::InputPhysicsBlock& ipb,
		       std::vector<panzer::BC>& bcs)
{
     panzer::InputEquationSet ies_1;
     ies_1.name = "Energy";
     ies_1.basis = "Q2";
     ies_1.integration_order = 2;
     ies_1.model_id = 6;
     ies_1.model_factory = "rf";
     ies_1.prefix = "";
     ies_1.params.set<int>("junk", 1);
   
     panzer::InputEquationSet ies_2;
     ies_2.name = "Energy";
     ies_2.basis = "Q2";
     ies_2.integration_order = 2;
     ies_2.model_id = 6;
     ies_2.model_factory = "rf";
     ies_2.prefix = "ION_";
     ies_2.params.set<int>("junk", 1);
   
     ipb.physics_block_id = "4";
     ipb.eq_sets.push_back(ies_1);
     ipb.eq_sets.push_back(ies_2);
   
     {
       std::size_t bc_id = 0;
       panzer::BCType neumann = panzer::BCT_Dirichlet;
       std::string sideset_id = "left";
       std::string element_block_id = "eblock-0_0";
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
       panzer::BCType neumann = panzer::BCT_Dirichlet;
       std::string sideset_id = "top";
       std::string element_block_id = "eblock-0_0";
       std::string dof_name = "TEMPERATURE";
       std::string strategy = "Constant";
       double value = -5.0;
       Teuchos::ParameterList p;
       p.set("Value",value);
       panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name, 
   		    strategy, p);
       bcs.push_back(bc);
     }    

     {
       std::size_t bc_id = 2;
       panzer::BCType neumann = panzer::BCT_Dirichlet;
       std::string sideset_id = "top";
       std::string element_block_id = "eblock-1_0";
       std::string dof_name = "ION_TEMPERATURE";
       std::string strategy = "Constant";
       double value = 20.0;
       Teuchos::ParameterList p;
       p.set("Value",value);
       panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name, 
   		    strategy, p);
       bcs.push_back(bc);
     }   

     {
       std::size_t bc_id = 3;
       panzer::BCType neumann = panzer::BCT_Dirichlet;
       std::string sideset_id = "right";
       std::string element_block_id = "eblock-1_0";
       std::string dof_name = "ION_TEMPERATURE";
       std::string strategy = "Constant";
       double value = -10.0;
       Teuchos::ParameterList p;
       p.set("Value",value);
       panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name, 
   		    strategy, p);
       bcs.push_back(bc);
     }
}
