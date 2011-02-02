#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Panzer_Workset_Builder.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_AssemblyEngine.hpp"
#include "Panzer_AssemblyEngine_TemplateManager.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_Basis.hpp"
#include "Panzer_AuxiliaryEvaluator_TemplateManager.hpp"

#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STK_GatherFields.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_STK_AuxiliaryVariables.hpp"

#include "user_app_EquationSetFactory.hpp"
#include "user_app_ModelFactory_TemplateBuilder.hpp"
#include "user_app_BCStrategy_Factory.hpp"

#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

#include <cstdio> // for get char
#include <vector>
#include <string>

namespace panzer {

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

  class DogBuilder {
     Teuchos::RCP<panzer_stk::STK_Interface> mesh_;
     Teuchos::RCP<panzer::Basis> basis_;
  public:
     DogBuilder(const Teuchos::RCP<panzer_stk::STK_Interface> & mesh,
                const Teuchos::RCP<panzer::Basis> & basis) 
        : mesh_(mesh), basis_(basis) {}
 
 
     template <typename EvalT>
     Teuchos::RCP<panzer::Base> build() const
     {
        Teuchos::RCP<std::vector<std::string> > fieldNames = Teuchos::rcp(new std::vector<std::string>);
        fieldNames->push_back("dog");

        Teuchos::ParameterList pl;
        pl.set<Teuchos::RCP<std::vector<std::string> > >("Field Names",fieldNames);
        pl.set("Basis",basis_);
  
        return Teuchos::rcp(new panzer_stk::AuxiliaryVariables<EvalT>(mesh_,pl));
     }
  };

  Teuchos::RCP<panzer::Basis> buildLinearBasis(std::size_t worksetSize);

  void testInitialzation(panzer::InputPhysicsBlock& ipb,
			 std::vector<panzer::BC>& bcs);

  Teuchos::RCP<panzer_stk::STK_Interface> buildMesh(int elemX,int elemY);

  TEUCHOS_UNIT_TEST(gs_evaluators, gather_constr)
  {
    const std::size_t workset_size = 20;
    Teuchos::RCP<panzer::Basis> linBasis = buildLinearBasis(workset_size);

    Teuchos::RCP<std::vector<std::string> > fieldNames
        = Teuchos::rcp(new std::vector<std::string>);
    fieldNames->push_back("dog");

    Teuchos::ParameterList pl;
    pl.set("Basis",linBasis);
    pl.set("Field Names",fieldNames);

    Teuchos::RCP<panzer_stk::STK_Interface> mesh = buildMesh(2,2);
    Teuchos::RCP<PHX::Evaluator<panzer::Traits> > gatherEval
        = Teuchos::rcp(new panzer_stk::GatherFields<panzer::Traits::Residual,panzer::Traits>(mesh,pl));

    RCP<Epetra_Comm> Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));

    panzer::InputPhysicsBlock ipb;
    std::vector<panzer::BC> bcs;
    testInitialzation(ipb, bcs);

    Teuchos::RCP<panzer::FieldManagerBuilder<int,int> > fmb = 
      Teuchos::rcp(new panzer::FieldManagerBuilder<int,int>);

    // build worksets
    //////////////////////////////////////////////////////////////
    std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > > 
      volume_worksets = panzer_stk::buildWorksets(*mesh,ipb, workset_size);

    const std::map<panzer::BC,Teuchos::RCP<std::map<unsigned,panzer::Workset> >,panzer::LessBC> bc_worksets 
       = panzer_stk::buildBCWorksets(*mesh,ipb,bcs);

    // build physics blocks
    //////////////////////////////////////////////////////////////
    user_app::MyFactory eqset_factory;
    user_app::BCFactory bc_factory;
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;

    {
      std::map<std::string,std::string> block_ids_to_physics_ids;
      block_ids_to_physics_ids["eblock-0_0"] = "test physics";
      block_ids_to_physics_ids["eblock-1_0"] = "test physics";
      
      std::map<std::string,panzer::InputPhysicsBlock> 
        physics_id_to_input_physics_blocks;
      physics_id_to_input_physics_blocks["test physics"] = ipb;
  
      fmb->buildPhysicsBlocks(block_ids_to_physics_ids,
                              physics_id_to_input_physics_blocks,
                              Teuchos::as<int>(mesh->getDimension()), workset_size,
                              eqset_factory,
                              physicsBlocks);
    }

    // build DOF Manager
    /////////////////////////////////////////////////////////////

    // build the connection manager 
    const Teuchos::RCP<panzer::ConnManager<int,int> > 
      conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    panzer::DOFManagerFactory<int,int> globalIndexerFactory;
    RCP<panzer::UniqueGlobalIndexer<int,int> > dofManager 
         = globalIndexerFactory.buildUniqueGlobalIndexer(MPI_COMM_WORLD,physicsBlocks,conn_manager);
 
    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory
          = Teuchos::rcp(new panzer::EpetraLinearObjFactory<panzer::Traits,int>(Comm,dofManager));

    // setup field manager build
    /////////////////////////////////////////////////////////////
 
    const DogBuilder db(mesh,linBasis);
    std::map<std::string,Teuchos::RCP<panzer::AuxiliaryEvaluator_TemplateManager<panzer::Traits> > > auxEval;
    auxEval["eblock-0_0"] = Teuchos::rcp(new panzer::AuxiliaryEvaluator_TemplateManager<panzer::Traits>);
    auxEval["eblock-0_0"]->buildAndPushBackObjects(db);
 
    fmb->setupVolumeFieldManagers(volume_worksets,physicsBlocks,dofManager,*linObjFactory,auxEval);

    fmb->setupBCFieldManagers(bc_worksets,physicsBlocks,eqset_factory,bc_factory,*linObjFactory);

    panzer::AssemblyEngine_TemplateManager<panzer::Traits,int,int> ae_tm;
    panzer::AssemblyEngine_TemplateBuilder<int,int> builder(fmb,linObjFactory);
    ae_tm.buildObjects(builder);

    panzer::AssemblyEngineInArgs input(
             linObjFactory->buildGhostedLinearObjContainer(),
             linObjFactory->buildLinearObjContainer());
    ae_tm.getAsObject<panzer::Traits::Residual>()->evaluate(input);
    ae_tm.getAsObject<panzer::Traits::Jacobian>()->evaluate(input);
  }

  Teuchos::RCP<panzer::Basis> buildLinearBasis(std::size_t worksetSize)
  {
     panzer::CellData cellData(worksetSize,2);
     panzer::IntegrationRule intRule(1,cellData);

     return Teuchos::rcp(new panzer::Basis("Q1",intRule)); 
  }

  Teuchos::RCP<panzer_stk::STK_Interface> buildMesh(int elemX,int elemY)
  {
    typedef panzer_stk::STK_Interface::SolutionFieldType VariableField;
    typedef panzer_stk::STK_Interface::VectorFieldType CoordinateField;

    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",2);
    pl->set("Y Blocks",1);
    pl->set("X Elements",elemX);
    pl->set("Y Elements",elemY);
    
    panzer_stk::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildUncommitedMesh(MPI_COMM_WORLD);
 
    // add in some fields
    mesh->addSolutionField("dog","eblock-0_0");
    mesh->addSolutionField("dog","eblock-1_0");

    factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD); 

    VariableField * field = mesh->getMetaData()->get_field<VariableField>("dog");
    CoordinateField * cField = mesh->getMetaData()->get_field<CoordinateField>("coordinates");
    TEUCHOS_ASSERT(field!=0);
    TEUCHOS_ASSERT(cField!=0);

    // fill the fields with data
    const std::vector<stk::mesh::Bucket*> nodeData 
        = mesh->getBulkData()->buckets(mesh->getNodeRank());
    for(std::size_t b=0;b<nodeData.size();++b) {
       stk::mesh::Bucket * bucket = nodeData[b];

       // build all buckets
       for(stk::mesh::Bucket::iterator itr=bucket->begin();
           itr!=bucket->end();++itr) {

          stk::mesh::EntityArray<CoordinateField> coordinates(*cField,*itr);
          stk::mesh::EntityArray<VariableField> dog_array(*field,*itr);

          double x = coordinates(0);
          double y = coordinates(1);

          dog_array() = 4.0*x*x+y;
       }
    }
    
    if(mesh->isWritable())
       mesh->writeToExodus("output.exo");

    return mesh;
  }

  void testInitialzation(panzer::InputPhysicsBlock& ipb,
			 std::vector<panzer::BC>& bcs)
  {
    panzer::InputEquationSet ies_1;
    ies_1.name = "Energy";
    ies_1.basis = "Q2";
    ies_1.integration_order = 1;
    ies_1.model_id = 6;
    ies_1.model_factory = "rf";
    ies_1.prefix = "";
    ies_1.params.set<int>("junk", 1);

    panzer::InputEquationSet ies_2;
    ies_2.name = "Energy";
    ies_2.basis = "Q1";
    ies_2.integration_order = 1;
    ies_2.model_id = 6;
    ies_2.model_factory = "rf";
    ies_2.prefix = "ION_";
    ies_2.params.set<int>("junk", 1);

    ipb.physics_block_id = "4";
    ipb.eq_sets.push_back(ies_1);
    ipb.eq_sets.push_back(ies_2);


    {
      std::size_t bc_id = 0;
      panzer::BCType neumann = BCT_Dirichlet;
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
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "right";
      std::string element_block_id = "eblock-1_0";
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
      std::size_t bc_id = 2;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "top";
      std::string element_block_id = "eblock-1_0";
      std::string dof_name = "TEMPERATURE";
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
