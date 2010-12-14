#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_AssemblyEngine.hpp"
#include "Panzer_AssemblyEngine_TemplateManager.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Panzer_DOFManager.hpp"
#include "user_app_EquationSetFactory.hpp"
#include "user_app_ModelFactory_TemplateBuilder.hpp"
#include "user_app_BCStrategy_Factory.hpp"

#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

#include <cstdio> // for get char

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


  std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > > 
  buildWorksets(const RCP<panzer_stk::STK_Interface>& mesh,
		const std::size_t workset_size);

  const std::map<panzer::BC,Teuchos::RCP<std::map<unsigned,panzer::Workset> >,panzer::LessBC> buildBCWorksets(const RCP<panzer_stk::STK_Interface>& mesh);

  template<typename Array>
  void getIdsAndVertices(const panzer_stk::STK_Interface& si,
			 std::string blockId,
			 std::vector<std::size_t>& localIds,
			 Array& vertices);

  void getNodeIds(const stk::mesh::Entity * element,
		  std::vector<stk::mesh::EntityId> & nodeIds);

  void getSideElements(const panzer_stk::STK_Interface & mesh,
		       const std::string & blockId, 
		       const std::vector<stk::mesh::Entity*> & sides,
		       std::vector<std::size_t> & localSideIds, 
		       std::vector<stk::mesh::Entity*> & elements);

  void testInitialzation(const panzer_stk::STK_Interface& mesh,
			 panzer::InputPhysicsBlock& ipb,
			 std::vector<panzer::BC>& bcs);



  TEUCHOS_UNIT_TEST(field_manager_builder, basic)
  {
    using Teuchos::RCP;
  
    // pause_to_attach();

    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",2);
    pl->set("Y Blocks",1);
    pl->set("X Elements",6);
    pl->set("Y Elements",4);
    
    panzer_stk::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
    RCP<Epetra_Comm> Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));

    panzer::InputPhysicsBlock ipb;
    std::vector<panzer::BC> bcs;
    testInitialzation(*mesh, ipb, bcs);

    const std::size_t workset_size = 20;
    std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > > 
      volume_worksets = panzer_stk::buildWorksets(*mesh,ipb, workset_size);

    const std::map<panzer::BC,Teuchos::RCP<std::map<unsigned,panzer::Workset> >,panzer::LessBC> bc_worksets 
       = panzer_stk::buildBCWorksets(*mesh,ipb,bcs);

    std::map<std::string,std::string> block_ids_to_physics_ids;
    block_ids_to_physics_ids["eblock-0_0"] = "test physics";
    block_ids_to_physics_ids["eblock-1_0"] = "test physics";
    
    std::map<std::string,panzer::InputPhysicsBlock> 
      physics_id_to_input_physics_blocks;
    physics_id_to_input_physics_blocks["test physics"] = ipb;

    const Teuchos::RCP<panzer::ConnManager<int,int> > 
      conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));
    
    user_app::MyFactory eqset_factory;
				  
    Teuchos::RCP<panzer::FieldManagerBuilder<int,int> > fmb = 
      Teuchos::rcp(new panzer::FieldManagerBuilder<int,int>);

    user_app::BCFactory bc_factory;

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

    panzer::AssemblyEngine_TemplateManager<panzer::Traits,int,int> ae_tm;
    panzer::AssemblyEngine_TemplateBuilder<int,int> builder(fmb);
    ae_tm.buildObjects(builder);


    RCP<DOFManager<int,int> > dofManager = 
      ae_tm.getAsObject<panzer::Traits::Residual>()->getManagerBuilder()->getDOFManager();
    panzer::EpetraLinearObjFactory<int> linObjFactory(Comm,dofManager);
    RCP<Epetra_Map> ghosted_map = linObjFactory.getGhostedMap();
    RCP<Epetra_CrsGraph> ghosted_graph = linObjFactory.getGhostedGraph();
    
    panzer::AssemblyEngineInArgs input;
    input.x = rcp(new Epetra_Vector(*ghosted_map));
    input.dxdt = rcp(new Epetra_Vector(*ghosted_map));
    input.f = rcp(new Epetra_Vector(*ghosted_map));
    input.j = rcp(new Epetra_CrsMatrix(Copy, *ghosted_graph));

    ae_tm.getAsObject<panzer::Traits::Residual>()->evaluate(input);
    ae_tm.getAsObject<panzer::Traits::Jacobian>()->evaluate(input);

    //input.f->Print(std::cout);
    //input.j->Print(std::cout);
  }

  std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > > 
  buildWorksets(const RCP<panzer_stk::STK_Interface>& mesh,
		const std::size_t workset_size)
  {
    std::vector<std::string> element_blocks;
    mesh->getElementBlockNames(element_blocks);
    int base_cell_dimension = 2;

    panzer::InputPhysicsBlock ipb;
    std::vector<panzer::BC> bcs;
    testInitialzation(*mesh, ipb, bcs);


    std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > > worksets;

    for (std::vector<std::string>::size_type i=0; i < element_blocks.size(); 
	 ++i) {

      std::vector<std::size_t> local_cell_ids;
      Intrepid::FieldContainer<double> cell_vertex_coordinates;

      panzer::getIdsAndVertices(*mesh, element_blocks[i], local_cell_ids, 
				cell_vertex_coordinates);

      worksets.insert(std::make_pair(element_blocks[i],panzer::buildWorksets(element_blocks[i],
							      local_cell_ids,
							      cell_vertex_coordinates,
							      ipb,
							      workset_size,
							      base_cell_dimension)));
    }
    
    return worksets;
  }

  const std::map<panzer::BC,Teuchos::RCP<std::map<unsigned,panzer::Workset> >,panzer::LessBC>
  buildBCWorksets(const RCP<panzer_stk::STK_Interface>& mesh)
  {
    using Teuchos::RCP;
    
    unsigned dim = mesh->getDimension();

    int base_cell_dimension = 2;

    panzer::InputPhysicsBlock ipb;
    std::vector<panzer::BC> bcs;
    testInitialzation(*mesh, ipb, bcs);

    std::vector<std::string> sideSets; 
    std::vector<std::string> elementBlocks; 
    mesh->getSidesetNames(sideSets);
    mesh->getElementBlockNames(elementBlocks);

    
    std::map<panzer::BC,Teuchos::RCP<std::map<unsigned,panzer::Workset> >,panzer::LessBC> bc_worksets;
    
    for (std::vector<panzer::BC>::const_iterator bc = bcs.begin();
	 bc != bcs.end(); ++bc) {
      
      std::vector<stk::mesh::Entity*> sideEntities; 
      mesh->getMySides(bc->sidesetID(),bc->elementBlockID(),sideEntities);
   
      
      std::vector<stk::mesh::Entity*> elements;
      std::vector<std::size_t> local_cell_ids;
      std::vector<std::size_t> local_side_ids;
      getSideElements(*mesh, bc->elementBlockID(),
		      sideEntities,local_side_ids,elements);

      Intrepid::FieldContainer<double> vertices;
      vertices.resize(elements.size(),4,dim);  
      
      // loop over elements of this block
      for(std::size_t elm=0;elm<elements.size();++elm) {
	std::vector<stk::mesh::EntityId> nodes;
	stk::mesh::Entity * element = elements[elm];
	
	local_cell_ids.push_back(mesh->elementLocalId(element));
	getNodeIds(element,nodes);
	
	TEUCHOS_ASSERT(nodes.size()==4);
	
	for(std::size_t v=0;v<nodes.size();++v) {
	  const double * coord = mesh->getNodeCoordinates(nodes[v]);
          
	  for(unsigned d=0;d<dim;++d) 
	    vertices(elm,v,d) = coord[d]; 
	}
      }
      
      Teuchos::RCP<std::map<unsigned,panzer::Workset> > workset = 
	buildBCWorkset(*bc, local_cell_ids, local_side_ids,
		       vertices, ipb, base_cell_dimension);
      

      bc_worksets[*bc] = workset;
    }
    
    return bc_worksets;
  }

  
// *******************************************
// functions
// *******************************************
 
  template<typename Array>
  void getIdsAndVertices(const panzer_stk::STK_Interface& si,
			 std::string blockId,
			 std::vector<std::size_t>& localIds,
			 Array& vertices) {
    
    std::vector<stk::mesh::Entity*> elements;
    si.getMyElements(blockId,elements);
    
    unsigned dim = si.getDimension();
    
    vertices.resize(elements.size(),4,dim);
    
    // loop over elements of this block
    for(std::size_t elm=0;elm<elements.size();++elm) {
      std::vector<stk::mesh::EntityId> nodes;
      stk::mesh::Entity * element = elements[elm];
      
      localIds.push_back(si.elementLocalId(element));
      getNodeIds(element,nodes);
      
      TEUCHOS_ASSERT(nodes.size()==4);
      
      for(std::size_t v=0;v<nodes.size();++v) {
	const double * coord = si.getNodeCoordinates(nodes[v]);
	
	for(unsigned d=0;d<dim;++d) 
	  vertices(elm,v,d) = coord[d]; 
      }
    }
  }


  void getNodeIds(const stk::mesh::Entity * element,
		  std::vector<stk::mesh::EntityId> & nodeIds)
  {
    stk::mesh::PairIterRelation nodeRel = element->relations(stk::mesh::Node);
    
    stk::mesh::PairIterRelation::iterator itr;
    for(itr=nodeRel.begin();itr!=nodeRel.end();++itr) 
      nodeIds.push_back(itr->entity()->identifier());
  }
  
  /** This function loops over the passed in set of "Sides" and looks
   * at there related elements. It is then determined which elements
   * belong in the requested element block, and what the local ID of 
   * the side is.
   *
   * \param[in] mesh STK mesh interface
   * \param[in] blockId Requested element block identifier
   * \param[in] sides Set of sides (entities of dimension-1) where
   *                  there is assumed part membership (induced or not)
   *                  in the requested element block.
   * \param[out] localSideIds On output this will contain the local side ids. 
   *             Assumed that on input <code>sides.size()==0</code>
   * \param[out] elements On output this will contain the elements associated
   *             with each side in the requested block. Assumed that on input
   *             <code>elements.size()==0</code>
   *
   * \note Some elements may be repeated in the lists, however the
   *       local side ID should be distinct for each of those.
   */
  void getSideElements(const panzer_stk::STK_Interface & mesh,
		       const std::string & blockId, 
		       const std::vector<stk::mesh::Entity*> & sides,
		       std::vector<std::size_t> & localSideIds, 
		       std::vector<stk::mesh::Entity*> & elements) 
  {
    // for verifying that an element is in specified block
    stk::mesh::Part * blockPart = mesh.getElementBlockPart(blockId);
    
    // loop over each side extracting elements and local side ID that
    // are containted in specified block.
    std::vector<stk::mesh::Entity*>::const_iterator sideItr;
    for(sideItr=sides.begin();sideItr!=sides.end();++sideItr) {
      stk::mesh::Entity * side = *sideItr;
      
      stk::mesh::PairIterRelation relations = 
	side->relations(stk::mesh::Element);

      for(std::size_t e=0;e<relations.size();++e) {
	stk::mesh::Entity * element = relations[e].entity();
	std::size_t sideId = relations[e].identifier();
	
         // is this element in requested block
	bool inBlock = element->bucket().member(*blockPart);
	if(inBlock) {
	  // add element and Side ID to output vectors
	  elements.push_back(element);
	  localSideIds.push_back(sideId);
	}
      }
    }
  }

  void testInitialzation(const panzer_stk::STK_Interface& mesh,
			 panzer::InputPhysicsBlock& ipb,
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
