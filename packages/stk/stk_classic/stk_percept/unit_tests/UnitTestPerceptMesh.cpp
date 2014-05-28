/*--------------------------------------------------------------------*/
/*    Copyright 2009, 2011 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stk_percept/function/StringFunction.hpp>
#include <stk_percept/function/FieldFunction.hpp>
#include <stk_percept/function/ConstantFunction.hpp>
#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/Util.hpp>
#include <stk_percept/ExceptionWatch.hpp>
#include <stk_percept/fixtures/Fixture.hpp>
#include <stk_percept/fixtures/Fixture.hpp>
#include <stk_percept/fixtures/QuadFixture.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <Teuchos_ScalarTraits.hpp>

#include <stk_percept/fixtures/Fixture.hpp>
#include <stk_percept/fixtures/BeamFixture.hpp>
#include <stk_percept/fixtures/HeterogeneousFixture.hpp>
#include <stk_percept/fixtures/QuadFixture.hpp>
#include <stk_percept/fixtures/WedgeFixture.hpp>

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>
#include <math.h>

namespace stk 
{
  namespace percept 
  {
    namespace unit_tests 
    {

#define EXTRA_PRINT 0
      static int print_infoLevel = 0;

      //the following defines where to put the input and output files created by this set of functions
#if 0
      static const std::string input_files_loc="./input_files/";
      static const std::string output_files_loc="./output_files/";
#else
      static const std::string input_files_loc="./input_files_";
      static const std::string output_files_loc="./output_files_";
#endif


      //=============================================================================
      //=============================================================================
      //=============================================================================

      /// create a mesh of hex elements for use in other tests below
      //STKUNIT_UNIT_TEST(unit_perceptMesh, build_meshes)
      static void fixture_setup_0()
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        // start_demoperceptMesh_hex8_build
        {
          percept::PerceptMesh eMesh(3u);

          unsigned p_size = eMesh.get_parallel_size();

          // generate a 4x4x(4*p_size) mesh
          std::string gmesh_spec = std::string("4x4x")+toString(4*p_size)+std::string("|bbox:0,0,0,1,1,1");
          // NLM eMesh.new_mesh(percept::GMeshSpec(gmesh_spec));
          eMesh.new_mesh(percept::GMeshSpec(gmesh_spec));
          eMesh.commit();
          eMesh.save_as(input_files_loc+"hex_fixture.e");

          // end_demo
        }

      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      /// using the QuadFixture, generate meshes with and without sidesets
      // STKUNIT_UNIT_TEST(unit_perceptMesh, quad4_quad4_4_test_1)
      static void fixture_setup_1()
      {

        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        MPI_Barrier( MPI_COMM_WORLD );

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        //if (p_size == 1 || p_size == 3)
        if (p_size <= 2)
          {
            //const unsigned p_rank = stk::parallel_machine_rank( pm );
            //const unsigned p_size = stk::parallel_machine_size( pm );

            const unsigned n = 12;
            //const unsigned nx = n , ny = n , nz = p_size*n ;
            const unsigned nx = n , ny = n;

            percept::QuadFixture<double> fixture( pm , nx , ny, true);
            fixture.meta_data.commit();
            fixture.generate_mesh();

            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data);
            eMesh.print_info("quad fixture",  print_infoLevel);
            eMesh.save_as(input_files_loc+"quad_fixture.e");
          }

        if (p_size <= 2)
          {
            //const unsigned p_rank = stk::parallel_machine_rank( pm );
            //const unsigned p_size = stk::parallel_machine_size( pm );

            const unsigned n = 12;
            //const unsigned nx = n , ny = n , nz = p_size*n ;
            const unsigned nx = n , ny = n;

            percept::QuadFixture<double> fixture( pm , nx , ny, false);
            fixture.meta_data.commit();
            fixture.generate_mesh();

            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data);
            eMesh.print_info("quad fixture no sidesets",  print_infoLevel);
            eMesh.save_as(input_files_loc+"quad_fixture_no_sidesets.e");
          }
      }
 
      static void fixture_setup()
      {
        static bool is_setup = false;
        if (is_setup) return;
        fixture_setup_0();
        fixture_setup_1();
        is_setup = true;
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      /// generate a mesh with wedge elements 

      STKUNIT_UNIT_TEST(unit_perceptMesh, wedge6_1)
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        // start_demo_perceptMesh_wedge6_1

        percept::PerceptMesh eMesh(3u);

        unsigned p_size = eMesh.get_parallel_size();
        if (p_size == 1)
          {

            percept::WedgeFixture wedgeFixture;
            wedgeFixture.createMesh(MPI_COMM_WORLD,
                                    4, 3, 2,
                                    0, 1,
                                    0, 1,
                                    0, 1,
                                    std::string("swept-wedge_0.e") );
          }
      }


      //=============================================================================
      //=============================================================================
      //=============================================================================

      /// a tutorial on some of the innards of stk_mesh database

      STKUNIT_UNIT_TEST(perceptMesh, walk_nodes)
      {
        fixture_setup();
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        MPI_Barrier( MPI_COMM_WORLD );

        const unsigned p_size = stk::parallel_machine_size( pm );
        const unsigned p_rank = stk::parallel_machine_rank( pm );
        if (p_size <= 2)
          {

            // create a 12x12 quad mesh with sidesets
            const unsigned n = 12;
            //const unsigned nx = n , ny = n , nz = p_size*n ;
            const unsigned nx = n , ny = n;

            bool sidesets_on = true;
            percept::QuadFixture<double> fixture( pm , nx , ny, sidesets_on);
            fixture.meta_data.commit();
            fixture.generate_mesh();

            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data);
            eMesh.print_info("quad fixture",  print_infoLevel);
            //eMesh.save_as("./output_files/quad_fixture.e");

            stk::mesh::fem::FEMMetaData& metaData = *eMesh.get_fem_meta_data();

            const std::vector< stk::mesh::Part * > & parts = metaData.get_parts();

            unsigned nparts = parts.size();
            if (1) std::cout << "Number of parts = " << nparts << std::endl;

            int surface_id = 2;
            std::string surface_name = "surface_"+toString(surface_id);
            stk::mesh::Part *part = eMesh.get_non_const_part(surface_name);
            stk::mesh::Selector in_surface_selector(*part);
            stk::mesh::BulkData& bulkData = *eMesh.get_bulk_data();
            VectorFieldType* coordField = eMesh.get_coordinates_field();

            const std::vector<stk::mesh::Bucket*> & buckets = bulkData.buckets( (eMesh.get_spatial_dim() == 2 ? eMesh.edge_rank() : eMesh.face_rank() ) );  // Note
            double sum = 0.0;

            for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
              {
                if (in_surface_selector(**k)) 
                  {
                    stk::mesh::Bucket & bucket = **k ;

                    // in case the cell topology is needed
                    const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(bucket);
                    shards::CellTopology cell_topo(cell_topo_data);

                    const unsigned num_elements_in_bucket = bucket.size();
                
                    for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                      {
                        stk::mesh::Entity& element = bucket[iElement];

                        const stk::mesh::PairIterRelation& elem_nodes = element.relations( stk::mesh::fem::FEMMetaData::NODE_RANK );  

                        unsigned num_node = elem_nodes.size(); 
                        for (unsigned inode=0; inode < num_node; inode++)
                          {
                            stk::mesh::Entity & node = *elem_nodes[ inode ].entity();
                            //stk::mesh::EntityId nid = node.identifier();

                            double * const coord = stk::mesh::field_data( *coordField , node );
                            // do something with coord's
                            sum += coord[0]*coord[0] + coord[1]*coord[1];
                          }
                      }
                  }
              }
            std::cout << "P[" << p_rank << ":" << p_size << "] sum = " << sum << std::endl;
          }
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      /// Test the mesh_difference capability of PerceptMesh and the interface to stk_io
      ///   1. read (and write and read back in) meshes generated above (quad_fixture)
      ///   2. invoke PerceptMesh::print_info(ostringstream...) to create a string representation of the mesh
      ///   3. compare the string with the saved, gold value of the string
      ///   4. invoke mesh_difference to ensure it behaves as expected (two meshes are shown as identical)
      ///   5. modify one mesh and ensure mesh_difference shows the meshes as being different

      STKUNIT_UNIT_TEST(perceptMesh, test_mesh_diff)
      {
        fixture_setup();
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        std::string expected_serialized_mesh_string = 
	  "P[0] ======================================================== P[0] ========================================================P[0] ========================================================P[0] PerceptMesh::print_info: quad fixtureP[0] Uses { Node = 169 Edge = 48 Face = 144 Elem = 0 }P[0] info>    Number of parts = 26 P[0] info>    Part subset info:  P[0] info>     Part[0]= {UNIVERSAL} topology = null primary_entity_rank = 4294967295 subsets = {{OWNS} , {SHARES} , {FEM_ROOT_CELL_TOPOLOGY_PART_Node} , {FEM_ROOT_CELL_TOPOLOGY_PART_Line_2} , {FEM_ROOT_CELL_TOPOLOGY_PART_Line_3} , {FEM_ROOT_CELL_TOPOLOGY_PART_Particle} , {FEM_ROOT_CELL_TOPOLOGY_PART_Triangle_3} , {FEM_ROOT_CELL_TOPOLOGY_PART_Triangle_6} , {FEM_ROOT_CELL_TOPOLOGY_PART_Triangle_4} , {FEM_ROOT_CELL_TOPOLOGY_PART_Quadrilateral_4} , {FEM_ROOT_CELL_TOPOLOGY_PART_Quadrilateral_8} , {FEM_ROOT_CELL_TOPOLOGY_PART_Quadrilateral_9} , {FEM_ROOT_CELL_TOPOLOGY_PART_Beam_2} , {FEM_ROOT_CELL_TOPOLOGY_PART_Beam_3} , {FEM_ROOT_CELL_TOPOLOGY_PART_ShellLine_2} , {FEM_ROOT_CELL_TOPOLOGY_PART_ShellLine_3} , block_1 , surface_1 , surface_2 , surface_3 , surface_4 , surface_quad4_edge2_1 , surface_quad4_edge2_2 , surface_quad4_edge2_3 , surface_quad4_edge2_4}P[0] info>     Part[1]= {OWNS} topology = null primary_entity_rank = 4294967295 subsets = {}P[0] info>     Part[2]= {SHARES} topology = null primary_entity_rank = 4294967295 subsets = {}P[0] info>     Part[3]= {FEM_ROOT_CELL_TOPOLOGY_PART_Node} topology = Node primary_entity_rank = 0 subsets = {}P[0] info>     Part[4]= {FEM_ROOT_CELL_TOPOLOGY_PART_Line_2} topology = Line_2 primary_entity_rank = 1 subsets = {surface_quad4_edge2_1 , surface_quad4_edge2_2 , surface_quad4_edge2_3 , surface_quad4_edge2_4}P[0] info>     Part[5]= {FEM_ROOT_CELL_TOPOLOGY_PART_Line_3} topology = Line_3 primary_entity_rank = 1 subsets = {}P[0] info>     Part[6]= {FEM_ROOT_CELL_TOPOLOGY_PART_Particle} topology = Particle primary_entity_rank = 2 subsets = {}P[0] info>     Part[7]= {FEM_ROOT_CELL_TOPOLOGY_PART_Triangle_3} topology = Triangle_3 primary_entity_rank = 2 subsets = {}P[0] info>     Part[8]= {FEM_ROOT_CELL_TOPOLOGY_PART_Triangle_6} topology = Triangle_6 primary_entity_rank = 2 subsets = {}P[0] info>     Part[9]= {FEM_ROOT_CELL_TOPOLOGY_PART_Triangle_4} topology = Triangle_4 primary_entity_rank = 2 subsets = {}P[0] info>     Part[10]= {FEM_ROOT_CELL_TOPOLOGY_PART_Quadrilateral_4} topology = Quadrilateral_4 primary_entity_rank = 2 subsets = {block_1}P[0] info>     Part[11]= {FEM_ROOT_CELL_TOPOLOGY_PART_Quadrilateral_8} topology = Quadrilateral_8 primary_entity_rank = 2 subsets = {}P[0] info>     Part[12]= {FEM_ROOT_CELL_TOPOLOGY_PART_Quadrilateral_9} topology = Quadrilateral_9 primary_entity_rank = 2 subsets = {}P[0] info>     Part[13]= {FEM_ROOT_CELL_TOPOLOGY_PART_Beam_2} topology = Beam_2 primary_entity_rank = 2 subsets = {}P[0] info>     Part[14]= {FEM_ROOT_CELL_TOPOLOGY_PART_Beam_3} topology = Beam_3 primary_entity_rank = 2 subsets = {}P[0] info>     Part[15]= {FEM_ROOT_CELL_TOPOLOGY_PART_ShellLine_2} topology = ShellLine_2 primary_entity_rank = 2 subsets = {}P[0] info>     Part[16]= {FEM_ROOT_CELL_TOPOLOGY_PART_ShellLine_3} topology = ShellLine_3 primary_entity_rank = 2 subsets = {}P[0] info>     Part[17]= block_1 topology = Quadrilateral_4 primary_entity_rank = 2 subsets = {}P[0] info>     Part[18]= surface_1 topology = null primary_entity_rank = 1 subsets = {surface_quad4_edge2_1}P[0] info>     Part[19]= surface_2 topology = null primary_entity_rank = 1 subsets = {surface_quad4_edge2_2}P[0] info>     Part[20]= surface_3 topology = null primary_entity_rank = 1 subsets = {surface_quad4_edge2_3}P[0] info>     Part[21]= surface_4 topology = null primary_entity_rank = 1 subsets = {surface_quad4_edge2_4}P[0] info>     Part[22]= surface_quad4_edge2_1 topology = Line_2 primary_entity_rank = 1 subsets = {}P[0] info>     Part[23]= surface_quad4_edge2_2 topology = Line_2 primary_entity_rank = 1 subsets = {}P[0] info>     Part[24]= surface_quad4_edge2_3 topology = Line_2 primary_entity_rank = 1 subsets = {}P[0] info>     Part[25]= surface_quad4_edge2_4 topology = Line_2 primary_entity_rank = 1 subsets = {} P[0] info>     Part Uses information:  P[0] info>     Part[0]= {UNIVERSAL} : Uses { Node = 169 Edge = 48 Face = 144 Elem = 0 }P[0] info>     Part[1]= {OWNS} : Uses { Node = 169 Edge = 48 Face = 144 Elem = 0 }P[0] info>     Part[2]= {SHARES} : Uses { Node = 0 Edge = 0 Face = 0 Elem = 0 }P[0] info>     Part[3]= {FEM_ROOT_CELL_TOPOLOGY_PART_Node} : Uses { Node = 0 Edge = 0 Face = 0 Elem = 0 }P[0] info>     Part[4]= {FEM_ROOT_CELL_TOPOLOGY_PART_Line_2} : Uses { Node = 48 Edge = 48 Face = 0 Elem = 0 }P[0] info>     Part[5]= {FEM_ROOT_CELL_TOPOLOGY_PART_Line_3} : Uses { Node = 0 Edge = 0 Face = 0 Elem = 0 }P[0] info>     Part[6]= {FEM_ROOT_CELL_TOPOLOGY_PART_Particle} : Uses { Node = 0 Edge = 0 Face = 0 Elem = 0 }P[0] info>     Part[7]= {FEM_ROOT_CELL_TOPOLOGY_PART_Triangle_3} : Uses { Node = 0 Edge = 0 Face = 0 Elem = 0 }P[0] info>     Part[8]= {FEM_ROOT_CELL_TOPOLOGY_PART_Triangle_6} : Uses { Node = 0 Edge = 0 Face = 0 Elem = 0 }P[0] info>     Part[9]= {FEM_ROOT_CELL_TOPOLOGY_PART_Triangle_4} : Uses { Node = 0 Edge = 0 Face = 0 Elem = 0 }P[0] info>     Part[10]= {FEM_ROOT_CELL_TOPOLOGY_PART_Quadrilateral_4} : Uses { Node = 169 Edge = 48 Face = 144 Elem = 0 }P[0] info>     Part[11]= {FEM_ROOT_CELL_TOPOLOGY_PART_Quadrilateral_8} : Uses { Node = 0 Edge = 0 Face = 0 Elem = 0 }P[0] info>     Part[12]= {FEM_ROOT_CELL_TOPOLOGY_PART_Quadrilateral_9} : Uses { Node = 0 Edge = 0 Face = 0 Elem = 0 }P[0] info>     Part[13]= {FEM_ROOT_CELL_TOPOLOGY_PART_Beam_2} : Uses { Node = 0 Edge = 0 Face = 0 Elem = 0 }P[0] info>     Part[14]= {FEM_ROOT_CELL_TOPOLOGY_PART_Beam_3} : Uses { Node = 0 Edge = 0 Face = 0 Elem = 0 }P[0] info>     Part[15]= {FEM_ROOT_CELL_TOPOLOGY_PART_ShellLine_2} : Uses { Node = 0 Edge = 0 Face = 0 Elem = 0 }P[0] info>     Part[16]= {FEM_ROOT_CELL_TOPOLOGY_PART_ShellLine_3} : Uses { Node = 0 Edge = 0 Face = 0 Elem = 0 }P[0] info>     Part[17]= block_1 : Uses { Node = 169 Edge = 48 Face = 144 Elem = 0 }P[0] info>     Part[18]= surface_1 : Uses { Node = 13 Edge = 12 Face = 0 Elem = 0 }P[0] info>     Part[19]= surface_2 : Uses { Node = 13 Edge = 12 Face = 0 Elem = 0 }P[0] info>     Part[20]= surface_3 : Uses { Node = 13 Edge = 12 Face = 0 Elem = 0 }P[0] info>     Part[21]= surface_4 : Uses { Node = 13 Edge = 12 Face = 0 Elem = 0 }P[0] info>     Part[22]= surface_quad4_edge2_1 : Uses { Node = 13 Edge = 12 Face = 0 Elem = 0 }P[0] info>     Part[23]= surface_quad4_edge2_2 : Uses { Node = 13 Edge = 12 Face = 0 Elem = 0 }P[0] info>     Part[24]= surface_quad4_edge2_3 : Uses { Node = 13 Edge = 12 Face = 0 Elem = 0 }P[0] info>     Part[25]= surface_quad4_edge2_4 : Uses { Node = 13 Edge = 12 Face = 0 Elem = 0 }P[0] info>    Number of fields = 6P[0] info>    Field[0]= coordinates rank= 1P[0] info>    number of field restrictions= 1P[0] info>    field restriction 0 stride[0] = 2 type= 0 ord= 0 which corresponds to Part= {UNIVERSAL}P[0] info>    Field[1]= surface_1_df rank= 1P[0] info>    number of field restrictions= 1P[0] info>    field restriction 0 stride[0] = 2 type= 1 ord= 22 which corresponds to Part= surface_quad4_edge2_1P[0] info>    Field[2]= surface_2_df rank= 1P[0] info>    number of field restrictions= 1P[0] info>    field restriction 0 stride[0] = 2 type= 1 ord= 23 which corresponds to Part= surface_quad4_edge2_2P[0] info>    Field[3]= surface_3_df rank= 1P[0] info>    number of field restrictions= 1P[0] info>    field restriction 0 stride[0] = 2 type= 1 ord= 24 which corresponds to Part= surface_quad4_edge2_3P[0] info>    Field[4]= surface_4_df rank= 1P[0] info>    number of field restrictions= 1P[0] info>    field restriction 0 stride[0] = 2 type= 1 ord= 25 which corresponds to Part= surface_quad4_edge2_4P[0] info>    Field[5]= distribution_factors rank= 0P[0] info>    number of field restrictions= 0 P[0] ======================================================== P[0] ========================================================P[0] ========================================================";


        const unsigned p_size = stk::parallel_machine_size( pm );
        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        if (p_size <= 2)
          {

            percept::PerceptMesh eMesh_0(2u);
            eMesh_0.open_read_only(input_files_loc+"quad_fixture.e");
            eMesh_0.save_as(input_files_loc+"quad_fixture_readwrite.e");

            MPI_Barrier( MPI_COMM_WORLD );

            percept::PerceptMesh eMesh_1(2u);
            percept::PerceptMesh eMesh_2(2u);
            eMesh_1.open_read_only(input_files_loc+"quad_fixture_readwrite.e");
            eMesh_2.open_read_only(input_files_loc+"quad_fixture.e");

            if (p_size == 1)
              {
                std::ostringstream serialized_mesh_1;
                std::ostringstream serialized_mesh_2;
                bool add_newlines = false;
                eMesh_1.print_info(serialized_mesh_1, "quad fixture", 2, add_newlines);
                eMesh_2.print_info(serialized_mesh_2, "quad fixture", 2, add_newlines);
                std::string serialized_mesh_string_1 = serialized_mesh_1.str();
                std::string serialized_mesh_string_2 = serialized_mesh_2.str();
                std::cout << "expected_serialized_mesh_string.size()= " << expected_serialized_mesh_string.size() << std::endl;
                std::cout << "serialized_mesh_1.size()= " << serialized_mesh_string_1.size() << std::endl;
                std::cout << "serialized_mesh_1=\n" << serialized_mesh_string_1 << std::endl;
                std::cout << "...serialized_mesh_1" << std::endl;
                //std::cout << "serialized_mesh_2= " << serialized_mesh_string_2 << std::endl;
                STKUNIT_EXPECT_TRUE(expected_serialized_mesh_string == serialized_mesh_string_2);
                STKUNIT_EXPECT_TRUE(serialized_mesh_string_1 == serialized_mesh_string_2);
              }

            MPI_Barrier( MPI_COMM_WORLD );

            {
              std::string diff_msg = "diff report: \n";
              bool diff = PerceptMesh::mesh_difference(eMesh_1, eMesh_2, diff_msg, true);
              STKUNIT_EXPECT_TRUE(!diff);
            }

            stk::mesh::fem::FEMMetaData& metaData_1 = *eMesh_1.get_fem_meta_data();
            stk::mesh::fem::FEMMetaData& metaData_2 = *eMesh_2.get_fem_meta_data();

            stk::mesh::BulkData& bulkData_1 = *eMesh_1.get_bulk_data();
            VectorFieldType* coordField_1 = eMesh_1.get_coordinates_field();
            stk::mesh::BulkData& bulkData_2 = *eMesh_2.get_bulk_data();
            //VectorFieldType* coordField_2 = eMesh_2.get_coordinates_field();

            MPI_Barrier( MPI_COMM_WORLD );

            {
              std::string diff_msg = "diff report 1: \n";
              bool diff = PerceptMesh::mesh_difference(metaData_1, metaData_2, bulkData_1, bulkData_2, diff_msg, true);
              STKUNIT_EXPECT_TRUE(!diff);
            }

            const std::vector<stk::mesh::Bucket*> & buckets = bulkData_1.buckets( stk::mesh::fem::FEMMetaData::NODE_RANK );

            for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
              {
                //if (in_surface_selector(**k)) 
                {
                  stk::mesh::Bucket & bucket = **k ;

                  const unsigned num_elements_in_bucket = bucket.size();
                
                  for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                    {
                      stk::mesh::Entity& entity = bucket[iEntity];

                      double * const coord = stk::mesh::field_data( *coordField_1 , entity );

                      coord[0] += 0.01;
                    }
                }
              }
            MPI_Barrier( MPI_COMM_WORLD );

            {
              std::string diff_msg = "diff report after mod: \n";
              bool diff = PerceptMesh::mesh_difference(eMesh_1, eMesh_2, diff_msg, true);
              STKUNIT_EXPECT_TRUE(diff);
            }


          }
        MPI_Barrier( MPI_COMM_WORLD );

      }


      STKUNIT_UNIT_TEST(perceptMesh, create_skewed_mesh)
      {
        bool notActive = false;
        if (notActive) return;

        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        MPI_Barrier( MPI_COMM_WORLD );

        const unsigned p_size = stk::parallel_machine_size( pm );
        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        if (p_size > 1) return;

        int thetas[] = {0,10,20,30,40,45,50,60,70,80};
        for (unsigned itheta=0; itheta < sizeof(thetas)/sizeof(thetas[0]); itheta++)
          {
            double theta = thetas[itheta];
            theta *= M_PI/180.0;
            // create a nxn quad mesh with sidesets
            const unsigned n = 20;
            //const unsigned nx = n , ny = n , nz = p_size*n ;
            const unsigned nx = n , ny = n;

            bool sidesets_on = true;
            percept::QuadFixture<double> fixture( pm , nx , ny, sidesets_on);
            fixture.set_bounding_box(0,1,0,1);
            fixture.meta_data.commit();
            fixture.generate_mesh();

            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data);
            //stk::mesh::fem::FEMMetaData& metaData = *eMesh.get_fem_meta_data();
            //const std::vector< stk::mesh::Part * > & parts = metaData.get_parts();

            if (0 == itheta) eMesh.save_as("2d_duct.e");
            stk::mesh::BulkData& bulkData = *eMesh.get_bulk_data();
            VectorFieldType* coordField = eMesh.get_coordinates_field();

            const std::vector<stk::mesh::Bucket*> & buckets = bulkData.buckets( stk::mesh::fem::FEMMetaData::NODE_RANK );  
            // right-shear, theta=angle from vertical, dxdy=tan(theta)
            // up-shear, theta=angle from horizontal, dydx=tan(theta)
            // choose one or other, set other to 0.0
            double dxdy = std::tan(theta);
            double dydx = 0.0;

            for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
              {
                //if (in_surface_selector(**k)) 
                {
                  stk::mesh::Bucket & bucket = **k ;

                  const unsigned num_nodes_in_bucket = bucket.size();
                
                  for (unsigned iNode = 0; iNode < num_nodes_in_bucket; iNode++)
                    {
                      stk::mesh::Entity& node = bucket[iNode];
                      //stk::mesh::EntityId nid = node.identifier();

                      double * const coord = stk::mesh::field_data( *coordField , node );
                  
                      coord[0] += dxdy*coord[1];
                      coord[1] += dydx*coord[0];
                    }
                }
              }
            if (dydx != 0)
              eMesh.save_as(std::string("slantThetaDYDX")+boost::lexical_cast<std::string>(thetas[itheta])+".g");
            else
              eMesh.save_as(std::string("slantTheta")+boost::lexical_cast<std::string>(thetas[itheta])+".g");

          }
      }

      STKUNIT_UNIT_TEST(perceptMesh, create_quad_streaming_mesh)
      {
        bool notActive = false;
        if (notActive) return;

        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        MPI_Barrier( MPI_COMM_WORLD );

        const unsigned p_size = stk::parallel_machine_size( pm );
        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        if (p_size != 2) return;

        const unsigned n = 2;
        //const unsigned nx = n , ny = n , nz = p_size*n ;
        const unsigned nx = n , ny = n;

        bool sidesets_on = false;
        percept::QuadFixture<double> fixture( pm , nx , ny, sidesets_on);
        fixture.set_bounding_box(0,1,0,1);
        fixture.meta_data.commit();
        fixture.generate_mesh();

        percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data);

        eMesh.save_as("quad_streaming.e");
      }

      STKUNIT_UNIT_TEST(perceptMesh, test_states)
      {
        fixture_setup();
        PerceptMesh eMesh(0);
        eMesh.open(input_files_loc+"hex_fixture.e");
        eMesh.add_coordinate_state_fields();
        eMesh.commit();
        //std::cout << "eMesh.get_coordinates_field()->number_of_states() = "  << eMesh.get_coordinates_field()->number_of_states() << std::endl;
        // field, dst, src
        stk::mesh::FieldBase * coordinates_N = eMesh.get_field("coordinates_N");
        //stk::mesh::FieldBase * coordinates_NM1 = eMesh.get_field("coordinates_NM1");
        stk::mesh::FieldBase * coordinates_None = eMesh.get_coordinates_field();

        // dst,src
        eMesh.copy_field(coordinates_N, coordinates_None);
        Math::Matrix rm = Math::rotationMatrix(0, 30);
        eMesh.transform_mesh(rm);
        eMesh.save_as(output_files_loc+"hex_copy_fields_rot.e");
        eMesh.copy_field(coordinates_None, coordinates_N);
        eMesh.save_as(output_files_loc+"hex_copy_fields.e");
      }

    }
  }
}
