#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_FancyOStream.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

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
#include "Panzer_DOFManager.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "user_app_EquationSetFactory.hpp"
#include "user_app_ModelFactory_TemplateBuilder.hpp"
#include "user_app_BCStrategy_Factory.hpp"

#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_VectorOut.h"

#include "write_solution_data.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

void testInitialzation(const panzer_stk::STK_Interface& mesh,
                       panzer::InputPhysicsBlock& ipb,
		       std::vector<panzer::BC>& bcs);

void prepForParallelSolve(Epetra_Export & exporter,const Epetra_CrsMatrix & inJac,const Epetra_Vector & inX,
                                                   Epetra_CrsMatrix & outJac,Epetra_Vector & outX);
void redistributeSolution(Epetra_Import & importer,const Epetra_Vector & inX,Epetra_Vector & outX);

   // calls MPI_Init and MPI_Finalize
int main(int argc,char * argv[])
{
   using Teuchos::RCP;

   Teuchos::GlobalMPISession mpiSession(&argc,&argv);

   Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
   out.setOutputToRootOnly(0);
   out.setShowProcRank(true);

   // build mesh
   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",2);
   pl->set("Y Blocks",1);
   pl->set("X Elements",10);
   pl->set("Y Elements",10);

   out << "BUILD MESH" << std::endl;
   panzer_stk::SquareQuadMeshFactory factory;
   factory.setParameterList(pl);

   // two part build f
   RCP<panzer_stk::STK_Interface> mesh = factory.buildUncommitedMesh(MPI_COMM_WORLD);
   mesh->addSolutionField("TEMPERATURE","eblock-0_0");
   mesh->addSolutionField("TEMPERATURE","eblock-1_0");
   mesh->addSolutionField("ION_TEMPERATURE","eblock-0_0");
   mesh->addSolutionField("ION_TEMPERATURE","eblock-1_0");
   factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD);
 
   // construct input physics
   out << "BUILD PHYSICS" << std::endl;
   panzer::InputPhysicsBlock ipb;
   std::vector<panzer::BC> bcs;
   testInitialzation(*mesh, ipb, bcs);

   std::map<std::string,std::string> block_ids_to_physics_ids;
   block_ids_to_physics_ids["eblock-0_0"] = "test physics";
   block_ids_to_physics_ids["eblock-1_0"] = "test physics";
   
   std::map<std::string,panzer::InputPhysicsBlock> 
     physics_id_to_input_physics_blocks;
   physics_id_to_input_physics_blocks["test physics"] = ipb;
 
   // build worksets: note that buildWorksets and buildBCWorksets both call testInitialization
   out << "BUILD WORKSETS" << std::endl;
   const std::size_t workset_size = 20;
   std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > > 
         volume_worksets = panzer_stk::buildWorksets(*mesh, ipb, workset_size);

   out << "block count = " << volume_worksets.size() << std::endl;
   out << "workset count = " << volume_worksets["eblock-0_0"]->size() << std::endl;
   
   out << "BUILD BC WORKSETS" << std::endl;
   std::map<panzer::BC,Teuchos::RCP<std::map<unsigned,panzer::Workset> >,panzer::LessBC> 
         bc_worksets = panzer_stk::buildBCWorksets(*mesh,ipb,bcs);
 
   out << "BUILD CONN MANAGER" << std::endl;
   // build the connection manager 
   const Teuchos::RCP<panzer::ConnManager<int,int> > 
     conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));
   
   // construct and setup the field manager builder
   out << "BUILD FMB" << std::endl;
   user_app::MyFactory eqset_factory;
 				  
   Teuchos::RCP<panzer::FieldManagerBuilder<int,int> > fmb = 
         Teuchos::rcp(new panzer::FieldManagerBuilder<int,int>);
 
   user_app::BCFactory bc_factory;
 
   out << "SETUP" << std::endl;
   fmb->setup(conn_manager,
 	      MPI_COMM_WORLD,
 	      block_ids_to_physics_ids,
 	      physics_id_to_input_physics_blocks,
 	      volume_worksets,
 	      bc_worksets,
 	      Teuchos::as<int>(mesh->getDimension()),
 	      eqset_factory,
 	      bc_factory,
 	      workset_size);

   // build assembly engine
   panzer::AssemblyEngine_TemplateManager<panzer::Traits,int,int> ae_tm;
   panzer::AssemblyEngine_TemplateBuilder<int,int> builder(fmb);
   ae_tm.buildObjects(builder);

   out << "BUILD LA" << std::endl;
 
   // construct some linear algebra object, build object to pass to evaluators
   RCP<panzer::DOFManager<int,int> > dofManager = 
     ae_tm.getAsObject<panzer::Traits::Residual>()->getManagerBuilder()->getDOFManager();
   RCP<Epetra_Map> ghosted_map = dofManager->getOverlapMap();
   RCP<Epetra_CrsGraph> ghosted_graph = dofManager->getOverlapGraph();
   RCP<Epetra_Export> exporter = dofManager->getOverlapExport();
   RCP<Epetra_Import> importer = dofManager->getOverlapImport();
   RCP<Epetra_Map> map = dofManager->getMap();
   RCP<Epetra_CrsGraph> graph = dofManager->getGraph();
   
   panzer::AssemblyEngineInArgs input;
   input.x = rcp(new Epetra_Vector(*ghosted_map));
   input.dxdt = rcp(new Epetra_Vector(*ghosted_map));
   input.f = rcp(new Epetra_Vector(*ghosted_map));
   input.j = rcp(new Epetra_CrsMatrix(Copy, *ghosted_graph));

   // evaluate physics
   out << "EVALUTE" << std::endl;
   ae_tm.getAsObject<panzer::Traits::Residual>()->evaluate(input);
   ae_tm.getAsObject<panzer::Traits::Jacobian>()->evaluate(input);

   out << "RAN SUCCESSFULLY!" << std::endl;

   out << "SOLVE" << std::endl;
   RCP<Epetra_CrsMatrix> epetra_A = rcp(new Epetra_CrsMatrix(Copy, *graph));
   RCP<Epetra_Vector> epetra_x = rcp(new Epetra_Vector(*map));
   RCP<Epetra_Vector> epetra_b = rcp(new Epetra_Vector(*map));

   // redistribute vectors and matrices so that we can parallel solve
   prepForParallelSolve(*exporter,*input.j,*input.f,*epetra_A,*epetra_b);

   // setup thyra stuff
   RCP<const Thyra::LinearOpBase<double> > A = Thyra::epetraLinearOp( epetra_A );
   RCP<Thyra::VectorBase<double> > x = Thyra::create_Vector( epetra_x, A->domain() );
   RCP<const Thyra::VectorBase<double> > b = Thyra::create_Vector( epetra_b, A->range());

   // solve with amesos
   Stratimikos::DefaultLinearSolverBuilder solverBuilder;
   Teuchos::RCP<Teuchos::ParameterList> validList = Teuchos::rcp(new Teuchos::ParameterList(*solverBuilder.getValidParameters()));
   solverBuilder.setParameterList(validList);
   
   RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory = solverBuilder.createLinearSolveStrategy("Amesos");
   RCP<Thyra::LinearOpWithSolveBase<double> > lows = Thyra::linearOpWithSolve(*lowsFactory, A);
   Thyra::solve<double>(*lows, Thyra::NOTRANS, *b, x.ptr());

   // redistribute solution vector
   redistributeSolution(*importer,*epetra_x,*input.x);

   out << "WRITE" << std::endl;
   write_solution_data(*dofManager,*mesh,*input.x);
   mesh->writeToExodus("output.exo");

   return 0;
}

void testInitialzation(const panzer_stk::STK_Interface& mesh,
			 panzer::InputPhysicsBlock& ipb,
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
       double value = 10.0;
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

void prepForParallelSolve(Epetra_Export & exporter,const Epetra_CrsMatrix & inJac,const Epetra_Vector & inX,Epetra_CrsMatrix & outJac,Epetra_Vector & outX)
{
    outJac.PutScalar(0.0);
    outJac.Export( inJac,exporter,Add);

    outX.PutScalar(0.0);
    outX.Export( inX,exporter,Add);
}

void redistributeSolution(Epetra_Import & importer,const Epetra_Vector & inX,Epetra_Vector & outX)
{
    outX.PutScalar(0.0);
    outX.Import(inX,importer,Insert);
}
