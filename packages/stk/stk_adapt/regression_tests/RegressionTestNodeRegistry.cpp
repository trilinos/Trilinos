/*--------------------------------------------------------------------*/
/*    Copyright 2009, 2011 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>

#include <math.h>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>

#include <Teuchos_ScalarTraits.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_percept/Util.hpp>
#include <stk_percept/ExceptionWatch.hpp>

#include <stk_percept/function/StringFunction.hpp>
#include <stk_percept/function/FieldFunction.hpp>
#include <stk_percept/function/ConstantFunction.hpp>
#include <stk_percept/PerceptMesh.hpp>

#include <stk_adapt/NodeRegistry.hpp>
#include <stk_percept/fixtures/Fixture.hpp>

namespace stk
{
  namespace adapt
  {
    namespace regression_tests
    {

#include "RegressionTestFileLoc.hpp"


#define EXTRA_PRINT 0

      //======================================================================================================================      
      //======================================================================================================================      
      //======================================================================================================================      
#if 1
      STKUNIT_UNIT_TEST(nodeRegistry_regr, test_parallel_0)
      {
        using namespace mesh;
        using namespace fem;

        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        // start_demo_nodeRegistry_test_parallel_1
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if ( p_size == 3)
          {

            // Set up meta and bulk data
            const unsigned spatial_dim = 2;
            FEMMetaData meta_data(spatial_dim);
            meta_data.commit();
            BulkData mesh(FEMMetaData::get_meta_data(meta_data), pm);
            //unsigned p_rank = mesh.parallel_rank();
            //unsigned p_size = mesh.parallel_size();
            const EntityRank elem_rank = meta_data.element_rank();

            // Begin modification cycle so we can create the entities and relations
            mesh.modification_begin();

            // We're just going to add everything to the universal part
            stk::mesh::PartVector empty_parts;
            stk::mesh::Part *  owned_part = & meta_data.locally_owned_part();
            stk::mesh::Part *  shared_part = & meta_data.globally_shared_part();
            stk::mesh::PartVector owned_parts(1,owned_part);
            stk::mesh::PartVector shared_parts(1,shared_part);

            // Create nodes
            Entity *node=0;
            node = &mesh.declare_entity(0,
                                        p_rank+1,
                                        empty_parts);


            Entity& elem = mesh.declare_entity(elem_rank,
                                               p_rank+1,
                                               empty_parts);


            if (p_rank==0 || p_rank==2)
              {
                node = &mesh.declare_entity(0,
                                            2,
                                            empty_parts);
                mesh.declare_relation(elem, *node, 0);
              }
            else
              {
                mesh.declare_relation(elem, *node, 0);
              }



            mesh.modification_end();
          }

        //exit(123);
      }
#endif
      //======================================================================================================================      
      //======================================================================================================================      
      //======================================================================================================================      
#if 1
      STKUNIT_UNIT_TEST(nodeRegistry_regr, test_parallel_1)
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        // start_demo_nodeRegistry_test_parallel_1
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {

            percept::PerceptMesh eMesh(2u);
            eMesh.open(input_files_loc+"break_test._.quad._.square._.square_quad4.e");
            eMesh.commit();

            //eMesh.print_info("square_quad4 before dist");
            //eMesh.save_as("./cube_hex9-orig.e");

            //unsigned p_size = eMesh.get_parallel_size();
            unsigned p_rank = eMesh.get_rank();
            Util::setRank(eMesh.get_rank());
        
            if (p_size != 3) // FIXME
              return;

            bool useGhosting = true;
            NodeRegistry nodeRegistry(eMesh, useGhosting);
            nodeRegistry.initialize();

            if (p_size == 3)
              {
                //if (p_rank != 0)
                {
                  // pick an element on the processor boundary
                  /* P[1] element_local = 1 Elem: 5 nodes: 9 11 12 10 
                   * P[1] element_ghost = 0 Elem: 11 nodes: 12 25 26 10 
                   * P[2] element_local = 1 Elem: 11 nodes: 12 25 26 10 
                   * P[2] element_ghost = 0 Elem: 5 nodes: 9 11 12 10 
                   */

                  // for proc 1
                  unsigned elem_num_local = 11;  // edge #3
                  unsigned elem_num_ghost = 5;  // edge #2
                  unsigned elem_20 = 20;

                  stk::mesh::Entity* element_local_p = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_num_local);
                  stk::mesh::Entity* element_ghost_p = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_num_ghost);
                  if (p_rank == 2)
                    {
                      element_local_p = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_num_ghost);
                      element_ghost_p = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_num_local);
                    }
                  if (p_rank == 0)
                    {
                      element_local_p = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_20);
                      element_ghost_p = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_20);
                    }

                  dw() << "P["<<p_rank<<"] elem_num_local = " << elem_num_local << DWENDL;
                  dw() << "P["<<p_rank<<"] elem_num_ghost = " << elem_num_ghost << DWENDL;

                  stk::mesh::Entity& element_local = *element_local_p;
                  stk::mesh::Entity& element_ghost = *element_ghost_p;

                  std::cout << "P["<<p_rank<<"] element_local isGhost = " << eMesh.isGhostElement(element_local) << " " << element_local << std::endl;
                  std::cout << "P["<<p_rank<<"] element_ghost isGhost = " << eMesh.isGhostElement(element_ghost) << " " << element_ghost << std::endl;

                  // choose edges to be used for new node locations (i.e., this would model a serendipity-like element with only edge Lagrange nodes)
                  NeededEntityType needed_entity_rank(eMesh.edge_rank(), 1u);
                  std::vector<NeededEntityType> needed_entity_ranks(1, needed_entity_rank);
            
                  /*
                   * 1st of three steps to create and associate new nodes - register need for new nodes, then check if node is remote, then get
                   *   from remote proc if necessary; finally, the local node database is ready to be queried 
                   *
                   * The pattern is to begin the step, loop over all elements (including ghosts) and invoke the local operation
                   * The method doForAllSubEntities is a utility for performing the operation on all the sub entities.
                   * If more granularity is desired, the member functions can be invoked directly for a particular sub-entity.
                   */
                  nodeRegistry.beginRegistration();
                  if (p_rank) 
                    {
                      nodeRegistry.doForAllSubEntities(&NodeRegistry::registerNeedNewNode, element_local, needed_entity_ranks);
                      nodeRegistry.doForAllSubEntities(&NodeRegistry::registerNeedNewNode, element_ghost, needed_entity_ranks);
                    }
                  nodeRegistry.endRegistration();

                  std::cout << "P["<<p_rank<<"] nodeRegistry size  = " << nodeRegistry.total_size() << std::endl;
                  std::cout << "P["<<p_rank<<"] nodeRegistry lsize = " << nodeRegistry.local_size() << std::endl;

                  dw() << "P["<<p_rank<<"] nodeRegistry size       = " << nodeRegistry.total_size() << DWENDL;
                  dw() << "P["<<p_rank<<"] nodeRegistry lsize      = " << nodeRegistry.local_size() << DWENDL;

                  // could do local create of elements here
                  nodeRegistry.beginLocalMeshMods();
                  nodeRegistry.endLocalMeshMods();

                  // check if the newly requested nodes are local or remote
                  nodeRegistry.beginCheckForRemote();  
                  if (p_rank)
                    {
                      nodeRegistry.doForAllSubEntities(&NodeRegistry::checkForRemote, element_local, needed_entity_ranks);
                      nodeRegistry.doForAllSubEntities(&NodeRegistry::checkForRemote, element_ghost, needed_entity_ranks);
                    }
                  nodeRegistry.endCheckForRemote();

                  // get the new nodes from other procs if they are nonlocal
                  nodeRegistry.beginGetFromRemote();
                  if (p_rank)
                    {
                      nodeRegistry.doForAllSubEntities(&NodeRegistry::getFromRemote, element_local, needed_entity_ranks);
                      nodeRegistry.doForAllSubEntities(&NodeRegistry::getFromRemote, element_ghost, needed_entity_ranks);
                    }
                  nodeRegistry.endGetFromRemote();

                  if (0)
                    {
                      MPI_Barrier( MPI_COMM_WORLD );
                      std::cout << "P["<< eMesh.get_rank() 
                                <<"] ========================================================================================================================" << std::endl;
                      nodeRegistry.checkDB();
                      MPI_Barrier( MPI_COMM_WORLD );
                      std::cout << "P["<< eMesh.get_rank() 

                                <<"] ========================================================================================================================" << std::endl;
                    }

                  // now we can get the new node's id and entity
                  if (p_rank)
                    {
                      unsigned iSubDimOrd = 3u;
                      if (p_rank == 2) 
                        {
                          iSubDimOrd = 2u;
                        }
                      NodeIdsOnSubDimEntityType& nodeIds_onSE = *(nodeRegistry.getNewNodesOnSubDimEntity(element_local, needed_entity_rank.first, iSubDimOrd));

                      //if (!nodeIds_onSE[0])
                      //  throw std::logic_error("nodeRegistry_regr.parallel_2 logic err3");

                      stk::mesh::Entity*  node   = eMesh.get_bulk_data()->get_entity(stk::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE.m_entity_id_vector[0]);

                      //EXPECT_EQ(nodeIds_onSE.m_entity_id_vector[0], 42u);
                      EXPECT_EQ(nodeIds_onSE.m_entity_id_vector[0], 41u);
                      // should be the same node on each proc
                      std::cout << "P[" << p_rank << "] nodeId = " << nodeIds_onSE << " node= " << node << std::endl;
                    }

                  // end_demo

                }
                //std::cout << "P[" << p_rank << "] exiting " << std::endl;
                //Util::pause(true);
                //eMesh.save_as("./cube_hex9.e");
                if (0)
                  {
                    MPI_Barrier( MPI_COMM_WORLD );
                    exit(1);
                  }
              }
          }
        //exit(123);

      }
#endif

      //======================================================================================================================      
      //======================================================================================================================      
      //======================================================================================================================      
#if 1
      STKUNIT_UNIT_TEST(nodeRegistry_regr, test_parallel_2)
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        // start_demo_nodeRegistry_test_parallel_2
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        std::cout << "nodeRegistry_regr.test_parallel_2: p_size = " << p_size << std::endl;

        if (p_size == 1 || p_size == 3)
          {

            percept::PerceptMesh eMesh(2u);
            eMesh.open(input_files_loc+"break_test._.quad._.square._.square_quad4.e");

            eMesh.commit();

            //eMesh.print_info("square_quad4 before dist");
            //eMesh.save_as("./cube_hex9-orig.e");

            //unsigned p_size = eMesh.get_parallel_size();
            unsigned p_rank = eMesh.get_rank();
            Util::setRank(eMesh.get_rank());
        
            if (p_size != 3) // FIXME
              return;

            bool useGhosting = true;
            NodeRegistry nodeRegistry(eMesh, useGhosting);
            nodeRegistry.initialize();

            if (p_size == 3)
              {
                //if (p_rank != 0)
                {
                  // pick an element on the processor boundary
                  /* P[1] element_local = 1 Elem: 5 nodes: 9 11 12 10 
                   * P[1] element_ghost = 0 Elem: 11 nodes: 12 25 26 10 
                   * P[2] element_local = 1 Elem: 11 nodes: 12 25 26 10 
                   * P[2] element_ghost = 0 Elem: 5 nodes: 9 11 12 10 
                   */

                  // for proc 1
                  unsigned elem_num_local = 11;  // edge #3
                  unsigned elem_num_ghost = 5;  // edge #2
                  unsigned elem_20 = 20;

                  unsigned elem_num_local_proc_0 = elem_20;
                  unsigned elem_num_local_proc_1 = elem_num_local;
                  unsigned elem_num_local_proc_2 = elem_num_ghost;

                  stk::mesh::Entity* element_local_p = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_num_local_proc_1);
                  stk::mesh::Entity* element_ghost_p = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_num_local_proc_2);
                  if (p_rank == 2)
                    {
                      element_local_p = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_num_local_proc_2);
                      element_ghost_p = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_num_local_proc_1);
                    }
                  if (p_rank == 0)
                    {
                      element_local_p = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_num_local_proc_0);
                      element_ghost_p = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_num_local_proc_0);
                    }

                  dw() << "P["<<p_rank<<"] elem_num_local = " << elem_num_local << DWENDL;
                  dw() << "P["<<p_rank<<"] elem_num_ghost = " << elem_num_ghost << DWENDL;

                  stk::mesh::Entity& element_local = *element_local_p;
                  stk::mesh::Entity& element_ghost = *element_ghost_p;

                  std::cout << "P["<< p_rank <<"] element_local isGhost = " << eMesh.isGhostElement(element_local) << " " << element_local << std::endl;
                  std::cout << "P["<< p_rank <<"] element_ghost isGhost = " << eMesh.isGhostElement(element_ghost) << " " << element_ghost << std::endl;

                  // choose edges to be used for new node locations (i.e., this would model a serendipity-like element with only edge Lagrange nodes)
                  std::vector<NeededEntityType> needed_entity_ranks(2);
                  needed_entity_ranks[0] = NeededEntityType(eMesh.edge_rank(), 1u);
                  needed_entity_ranks[1] = NeededEntityType(eMesh.element_rank(), 1u);
            
                  /*
                   * 1st of three steps to create and associate new nodes - register need for new nodes, then check if node is remote, then get
                   *   from remote proc if necessary; finally, the local node database is ready to be queried 
                   *
                   * The pattern is to begin the step, loop over all elements (including ghosts) and invoke the local operation
                   * The method doForAllSubEntities is a utility for performing the operation on all the sub entities.
                   * If more granularity is desired, the member functions can be invoked directly for a particular sub-entity.
                   */
                  nodeRegistry.beginRegistration();
                  if (p_rank) 
                    {
                      nodeRegistry.doForAllSubEntities(&NodeRegistry::registerNeedNewNode, element_local, needed_entity_ranks);
                      nodeRegistry.doForAllSubEntities(&NodeRegistry::registerNeedNewNode, element_ghost, needed_entity_ranks);
                    }
                  nodeRegistry.endRegistration();

                  std::cout << "P["<<p_rank<<"] nodeRegistry size  = " << nodeRegistry.total_size() << std::endl;
                  std::cout << "P["<<p_rank<<"] nodeRegistry lsize = " << nodeRegistry.local_size() << std::endl;

                  dw() << "P["<<p_rank<<"] nodeRegistry size       = " << nodeRegistry.total_size() << DWENDL;
                  dw() << "P["<<p_rank<<"] nodeRegistry lsize      = " << nodeRegistry.local_size() << DWENDL;

                  // could do local create of elements here
                  nodeRegistry.beginLocalMeshMods();
                  nodeRegistry.endLocalMeshMods();

                  // check if the newly requested nodes are local or remote
                  nodeRegistry.beginCheckForRemote();  
                  if (p_rank)
                    {
                      nodeRegistry.doForAllSubEntities(&NodeRegistry::checkForRemote, element_local, needed_entity_ranks);
                      nodeRegistry.doForAllSubEntities(&NodeRegistry::checkForRemote, element_ghost, needed_entity_ranks);
                    }
                  nodeRegistry.endCheckForRemote();

                  // get the new nodes from other procs if they are nonlocal
                  nodeRegistry.beginGetFromRemote();
                  if (p_rank)
                    {
                      nodeRegistry.doForAllSubEntities(&NodeRegistry::getFromRemote, element_local, needed_entity_ranks);
                      nodeRegistry.doForAllSubEntities(&NodeRegistry::getFromRemote, element_ghost, needed_entity_ranks);
                    }
                  nodeRegistry.endGetFromRemote();

                  if (0)
                    {
                      MPI_Barrier( MPI_COMM_WORLD );
                      std::cout << "P["<< eMesh.get_rank() 
                                <<"] ========================================================================================================================" << std::endl;
                      nodeRegistry.checkDB();
                      MPI_Barrier( MPI_COMM_WORLD );
                      std::cout << "P["<< eMesh.get_rank() 

                                <<"] ========================================================================================================================" << std::endl;
                    }

                  // now we can get the new node's id and entity
                  if (p_rank)
                    {
                      unsigned iSubDimOrd = 3u;
                      if (p_rank == 2) 
                        {
                          iSubDimOrd = 2u;
                        }
                      NodeIdsOnSubDimEntityType& nodeIds_onSE = *(nodeRegistry.getNewNodesOnSubDimEntity(element_local, needed_entity_ranks[0].first, iSubDimOrd));
                      if (!nodeIds_onSE[0])
                        throw std::logic_error("nodeRegistry_regr.parallel_2 logic err1");
                      stk::mesh::Entity*  node   = eMesh.get_bulk_data()->get_entity(stk::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE[0]->identifier());

                      //EXPECT_EQ(nodeId, 42u);
                      // should be the same node on each proc
                      std::cout << "P[" << p_rank << "] nodeId = " << nodeIds_onSE << " node= " << node << std::endl;
                    }

                  if (1)
                    {
                      if (p_rank)
                        {
                          NodeIdsOnSubDimEntityType& nodeIds_onSE_1 = *(nodeRegistry.getNewNodesOnSubDimEntity(element_local, needed_entity_ranks[1].first, 0u));
                          if (!nodeIds_onSE_1[0])
                            throw std::logic_error("nodeRegistry_regr.parallel_2 logic err2");

                          stk::mesh::Entity*  node_1   = eMesh.get_bulk_data()->get_entity(stk::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE_1[0]->identifier());

                          std::cout << "P[" << p_rank << "] nodeId_1 = " << nodeIds_onSE_1 << " node_1= " << node_1 << std::endl;


                          unsigned expectedId= 37u;
#if NODE_REGISTRY_MAP_TYPE_BOOST
#endif
#if NODE_REGISTRY_MAP_TYPE_STD
                          expectedId= 38u;
#endif
#if NODE_REGISTRY_MAP_TYPE_GOOGLE
                          expectedId= 39u;
#endif
#if SDS_ENTITY_TYPE_ID
                          expectedId= 37u;
#else
                          expectedId= 37u;
#endif
                          

                          unsigned expectedId_p2= 41u;
#if SDS_ENTITY_TYPE_ID
                          expectedId_p2= 41u;
#else
                          expectedId_p2= 41u;
#endif

                          if (p_rank==1) std::cout << "P["<<p_rank<<"] nodeIds_onSE_1[0]= " << nodeIds_onSE_1.m_entity_id_vector[0] << "should be " << expectedId    << std::endl;
                          if (p_rank==2) std::cout << "P["<<p_rank<<"] nodeIds_onSE_1[0]= " << nodeIds_onSE_1.m_entity_id_vector[0] << "should be " << expectedId_p2 << std::endl;
                          
                          if (p_rank==1) EXPECT_EQ(nodeIds_onSE_1[0]->identifier(), expectedId); 
                          if (p_rank==2) EXPECT_EQ(nodeIds_onSE_1[0]->identifier(), expectedId_p2);
                        }

                    }

                  // end_demo
                }
                //std::cout << "P[" << p_rank << "] exiting " << std::endl;
                //Util::pause(true);
                //eMesh.save_as("./cube_hex9.e");
                if (0)
                  {
                    MPI_Barrier( MPI_COMM_WORLD );
                    exit(1);
                  }
              }
          }
      }
#endif

    }//    namespace unit_tests
  }//  namespace adapt
}// namespace stk

