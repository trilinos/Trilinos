/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/util/PrintTable.hpp>

#include <Teuchos_ScalarTraits.hpp>

#include <stk_percept/Util.hpp>
#include <stk_percept/ExceptionWatch.hpp>
#include <stk_percept/function/StringFunction.hpp>
#include <stk_percept/function/FieldFunction.hpp>
#include <stk_percept/function/ConstantFunction.hpp>
#include <stk_percept/PerceptMesh.hpp>

#include <stk_adapt/UniformRefinerPattern.hpp>
#include <stk_adapt/UniformRefiner.hpp>
#include <unit_tests/TestLocalRefinerTri.hpp>
#include <unit_tests/TestLocalRefinerTri1.hpp>
#include <unit_tests/TestLocalRefinerTri2.hpp>
#include <unit_tests/TestLocalRefinerTri_N.hpp>
#include <unit_tests/TestLocalRefinerTri_N_1.hpp>

#include <unit_tests/TestLocalRefinerTet_N_1.hpp>
#include <unit_tests/TestLocalRefinerTet_N_2.hpp>


#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <unit_tests/UnitTestSupport.hpp>

#include <stk_io/IossBridge.hpp>

#include <boost/lexical_cast.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/util/PrintTable.hpp>

#include <stk_percept/Percept.hpp>
#include <stk_percept/Util.hpp>
#include <stk_percept/ExceptionWatch.hpp>

#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>
#include <stk_percept/RunEnvironment.hpp>

#include <stk_percept/fixtures/Fixture.hpp>
#include <stk_percept/fixtures/SingleTetFixture.hpp>
#include <stk_percept/fixtures/BeamFixture.hpp>
#include <stk_percept/fixtures/HeterogeneousFixture.hpp>
#include <stk_percept/fixtures/QuadFixture.hpp>
#include <stk_percept/fixtures/WedgeFixture.hpp>

#include <use_cases/UseCase_3.hpp>

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <typeinfo>

#include <math.h>
#include <stk_util/parallel/Parallel.hpp>

namespace stk {
  namespace adapt {
    namespace unit_tests {

      static int printInfoLevel = 0;

      /// configuration: you can choose where to put the generated Exodus files (see variables input_files_loc, output_files_loc)
      /// The following defines where to put the input and output files created by this set of functions

#if 1
      const std::string input_files_loc="./input_files_";
      const std::string output_files_loc="./output_files_";
#else
      const std::string input_files_loc="./input_files/";
      const std::string output_files_loc="./output_files/";
#endif

#define EXTRA_PRINT 0

      /// This function either writes the given mesh to a file in Exodus format (option 0)
      ///   or, under option 1, checks if the file already exists, and if so, treats that
      ///   file as the "gold" copy and does a regression difference check.

      static void save_or_diff(PerceptMesh& eMesh, std::string filename, int option = 0)
      {
        return UnitTestSupport::save_or_diff(eMesh, filename, option);
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================
      /// check triangulate_tet

      STKUNIT_UNIT_TEST(unit_localRefiner, triangulate_tet)
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size(pm);

        if (p_size <= 1)
          {
            {
              // create the mesh

              stk::percept::SingleTetFixture mesh(pm, false);
              stk::io::put_io_part_attribute(  mesh.m_block_tet );
              mesh.m_metaData.commit();
              mesh.populate();

              std::cout << "here" << std::endl;
              bool isCommitted = true;
              percept::PerceptMesh eMesh(&mesh.m_metaData, &mesh.m_bulkData, isCommitted);
              eMesh.saveAs(input_files_loc+"local_tet_0.e");
            }

            {
              PerceptMesh eMesh;
              eMesh.open(input_files_loc+"local_tet_0.e");
              Local_Tet4_Tet4_N break_tet(eMesh);
              eMesh.commit();

              TestLocalRefinerTet_N_1 breaker(eMesh, break_tet, 0);
              breaker.setRemoveOldElements(false);
              breaker.doBreak();

              save_or_diff(eMesh, output_files_loc+"local_tet_N_1_1.e");
            }

            {
              PerceptMesh eMesh;
              eMesh.open(input_files_loc+"local_tet_0.e");
              Local_Tet4_Tet4_N break_tet(eMesh);
              eMesh.commit();

              TestLocalRefinerTet_N_2 breaker(eMesh, break_tet, 0);
              breaker.setRemoveOldElements(true);
              breaker.doBreak();

              save_or_diff(eMesh, output_files_loc+"local_tet_N_2_1.e");
              exit(123);
            }

          }
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      /// Create a triangle mesh using the QuadFixture with the option of breaking the quads into triangles
      /// Refine the triangle mesh, write the results.

      /// Refine a triangle mesh

      STKUNIT_UNIT_TEST(unit_localRefiner, break_tri_to_tri)
      {
        //fixture_setup();
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size <= 3)
          {
            // start_demo_local_refiner_break_tri_to_tri

            const unsigned n = 1;
            const unsigned nx = n , ny = n;

            bool createEdgeSets = false;
            percept::QuadFixture<double, shards::Triangle<3> > fixture( pm , nx , ny, createEdgeSets);

            bool isCommitted = false;
            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, isCommitted);

            Local_Tri3_Tri3_2 break_tri_to_tri_2(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
            eMesh.addField("proc_rank_edge", eMesh.edge_rank(), scalarDimension);
            eMesh.commit();

            fixture.generate_mesh();

            eMesh.printInfo("local tri mesh", printInfoLevel);
            save_or_diff(eMesh, output_files_loc+"local_tri_0.e");

            TestLocalRefinerTri breaker(eMesh, break_tri_to_tri_2, proc_rank_field);
            breaker.setRemoveOldElements(false);
            breaker.doBreak();

            eMesh.printInfo("local tri mesh refined", printInfoLevel);
            save_or_diff(eMesh, output_files_loc+"local_tri_1.e");

            // end_demo
          }

      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      /// Refine a triangle mesh by trying to mark only one edge per triangle, in a random-ish way


      STKUNIT_UNIT_TEST(unit_localRefiner, break_tri_to_tri_1)
      {
        //fixture_setup();
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size <= 3)
          {
            // start_demo_local_refiner_break_tri_to_tri_1

            const unsigned n = 4;
            const unsigned nx = n , ny = n;

            bool createEdgeSets = false;
            percept::QuadFixture<double, shards::Triangle<3> > fixture( pm , nx , ny, createEdgeSets);

            bool isCommitted = false;
            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, isCommitted);

            Local_Tri3_Tri3_2 break_tri_to_tri_2(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
            eMesh.addField("proc_rank_edge", eMesh.edge_rank(), scalarDimension);
            eMesh.commit();

            fixture.generate_mesh();

            eMesh.printInfo("local tri mesh", printInfoLevel);
            save_or_diff(eMesh, output_files_loc+"local_tri_1_0.e");

            bool diagonals=true;
            TestLocalRefinerTri2 breaker(eMesh, break_tri_to_tri_2, proc_rank_field, diagonals);
            //breaker.setRemoveOldElements(false);
            breaker.doBreak();

            eMesh.printInfo("local tri mesh refined", printInfoLevel);
            save_or_diff(eMesh, output_files_loc+"local_tri_1_1.e");

            // end_demo
          }

      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      /// Refine a triangle mesh by trying to mark only one edge per triangle, in a random-ish way


      STKUNIT_UNIT_TEST(unit_localRefiner, break_tri_to_tri_2)
      {
        //fixture_setup();
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size <= 3)
          {
            // start_demo_local_refiner_break_tri_to_tri_2

            const unsigned n = 4;
            const unsigned nx = n , ny = n;

            bool createEdgeSets = false;
            percept::QuadFixture<double, shards::Triangle<3> > fixture( pm , nx , ny, createEdgeSets);

            bool isCommitted = false;
            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, isCommitted);

            Local_Tri3_Tri3_2 break_tri_to_tri_2(eMesh);
            //Local_Tri3_Tri3_N break_tri_to_tri_2(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
            eMesh.addField("proc_rank_edge", eMesh.edge_rank(), scalarDimension);
            eMesh.commit();

            fixture.generate_mesh();

            eMesh.printInfo("local tri mesh", printInfoLevel);
            save_or_diff(eMesh, output_files_loc+"local_tri_2_0.e");

            bool diagonals=false;
            TestLocalRefinerTri2 breaker(eMesh, break_tri_to_tri_2, proc_rank_field, diagonals);
            //breaker.setRemoveOldElements(false);
            breaker.doBreak();

            eMesh.printInfo("local tri mesh refined",  printInfoLevel);
            save_or_diff(eMesh, output_files_loc+"local_tri_2_1.e");

            // end_demo
          }

      }

      //=============================================================================
      //=============================================================================
      //=============================================================================


#if 0
      STKUNIT_UNIT_TEST(unit_localRefiner, break_tri_to_tri_N)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size <= 3)
          {
            // start_demo_local_refiner_break_tri_to_tri_1

            const unsigned n = 4;
            const unsigned nx = n , ny = n;

            bool createEdgeSets = false;
            percept::QuadFixture<double, shards::Triangle<3> > fixture( pm , nx , ny, createEdgeSets);

            bool isCommitted = false;
            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, isCommitted);

            Local_Tri3_Tri3_N break_tri_to_tri_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
            eMesh.addField("proc_rank_edge", eMesh.edge_rank(), scalarDimension);
            eMesh.commit();

            fixture.generate_mesh();

            eMesh.printInfo("local tri mesh", printInfoLevel);
            save_or_diff(eMesh, output_files_loc+"local_tri_N_0.e");

            TestLocalRefinerTri_N breaker(eMesh, break_tri_to_tri_N, proc_rank_field);
            breaker.setRemoveOldElements(false);
            breaker.doBreak();

            eMesh.printInfo("local tri mesh refined",  printInfoLevel);
            //eMesh.dumpElements();
            save_or_diff(eMesh, output_files_loc+"local_tri_N_1.e");

            //breaker.unrefineAll();
            ElementUnrefineCollection elements_to_unref = breaker.buildTestUnrefList();
            breaker.unrefineTheseElements(elements_to_unref);

            save_or_diff(eMesh, output_files_loc+"local_tri_N_1_unref.e");

            // end_demo
          }

      }
#endif

      //=============================================================================
      //=============================================================================
      //=============================================================================


#if 0
      STKUNIT_UNIT_TEST(unit_localRefiner, break_tri_to_tri_N_1)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size <= 3)
          {
            // start_demo_local_refiner_break_tri_to_tri_1

            const unsigned n = 4;
            const unsigned nx = n , ny = n;

            bool createEdgeSets = false;
            percept::QuadFixture<double, shards::Triangle<3> > fixture( pm , nx , ny, createEdgeSets);

            bool isCommitted = false;
            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, isCommitted);

            Local_Tri3_Tri3_N break_tri_to_tri_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
            eMesh.addField("proc_rank_edge", eMesh.edge_rank(), scalarDimension);
            eMesh.commit();

            fixture.generate_mesh();

            eMesh.printInfo("local tri mesh", printInfoLevel);
            save_or_diff(eMesh, output_files_loc+"local_tri_N_1_0.e");

            TestLocalRefinerTri_N_1 breaker(eMesh, break_tri_to_tri_N, proc_rank_field);
            breaker.setRemoveOldElements(false);
            breaker.doBreak();

            eMesh.printInfo("local tri mesh refined",  printInfoLevel);
            //eMesh.dumpElements();
            save_or_diff(eMesh, output_files_loc+"local_tri_N_1_1.e");

            //breaker.unrefineAll();
            ElementUnrefineCollection elements_to_unref = breaker.buildTestUnrefList();
            breaker.unrefineTheseElements(elements_to_unref);

            save_or_diff(eMesh, output_files_loc+"local_tri_N_1_1_unref.e");

            // end_demo
          }

      }
#endif



      //=============================================================================
      //=============================================================================
      //=============================================================================

      struct SingleTriangleFixture
      {
        static PerceptMesh * create()
        {
            const unsigned n = 1;
            const unsigned nx = n , ny = n;

            stk::ParallelMachine pm = MPI_COMM_WORLD ;
            bool createEdgeSets = false;
            percept::QuadFixture<double, shards::Triangle<3> > *fixture = new percept::QuadFixture<double, shards::Triangle<3> >( pm , nx , ny, createEdgeSets);

            bool isCommitted = false;
            percept::PerceptMesh * eMesh = new percept::PerceptMesh(&fixture->meta_data, &fixture->bulk_data, isCommitted);

            eMesh->commit();

            fixture->generate_mesh();

            // delete the first element
            eMesh->getBulkData()->modification_begin();
            stk::mesh::Entity* element = &( (**(eMesh->getBulkData()->buckets(eMesh->element_rank()).begin()))[0]);
            if ( ! eMesh->getBulkData()->destroy_entity( element ) )
              {
                throw std::logic_error("failed in deleting element");
              }
            eMesh->getBulkData()->modification_end();

            // single element left
            //stk::mesh::Entity& element = (**(eMesh->getBulkData()->buckets(eMesh->element_rank()).begin()))[0];

            return eMesh;
        }
      };


      static void set_node_coords(percept::PerceptMesh& eMesh, mesh::PairIterRelation& elem_nodes, double tri_coords[3][3])
      {
        for (unsigned inode=0; inode < elem_nodes.size(); inode++)
          {
            stk::mesh::Entity *node = elem_nodes[inode].entity();
            double *fdata = stk::mesh::field_data( *eMesh.getCoordinatesField() , *node );
            for (int dim=0; dim < eMesh.getSpatialDim(); dim++)
              {
                fdata[dim] = tri_coords[inode][dim];
              }
          }
      }

      std::vector<int> convert_tuple(tri_tuple_type_local& tuple)
      {
        std::vector<int> cv(3);
        cv[0] = tuple.get<0>();
        cv[1] = tuple.get<1>();
        cv[2] = tuple.get<2>();
        return cv;
      }

      static bool in_set(tri_tuple_type_local& expected, vector<tri_tuple_type_local>& base, bool reverse=false)
      {
        std::vector<int> cv_expected = convert_tuple(expected);

        for (unsigned ie = 0; ie < base.size(); ie++)
          {
            std::vector<int> cv_base = convert_tuple(base[ie]);
            for (int i = 0; i < 3; i++)
              {
                bool found = true;
                if (reverse)
                  {
                    int k=0;
                    for (int j = 2; j >= 0; --j)
                      {
                        if (cv_expected[k++] != cv_base[(i+j)%3])
                          {
                            found = false;
                            break;
                          }
                      }
                  }
                else
                  {
                    for (int j = 0; j < 3; j++)
                      {
                        if (cv_expected[j] != cv_base[(i+j)%3])
                          {
                            found = false;
                            break;
                          }
                      }
                  }

                if (found)
                  {
                    return true;
                  }
              }
          }
        return false;
      }

      /// Create a single triangle mesh and mark the edges, call RefinerPattern_Tri3_Tri3_N::triangulate_face
      ///   and check properties of the result - reverse the triangle polarity and check for consistency

      STKUNIT_UNIT_TEST(unit_localRefiner, check_triangulate_face)
      {
        //fixture_setup();
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size <= 1)
          {
            percept::PerceptMesh& eMesh = *SingleTriangleFixture::create();

            eMesh.saveAs(output_files_loc+"tri_face_0.e");

            // single element left
            stk::mesh::Entity& element = (**(eMesh.getBulkData()->buckets(eMesh.element_rank()).begin()))[0];
            std::cout << "element = " << element << std::endl;

            mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);


            stk::mesh::Entity *elem_nodes_vector[3];
            for (unsigned inode=0; inode < elem_nodes.size(); inode++)
              {
                stk::mesh::Entity *node = elem_nodes[inode].entity();
                elem_nodes_vector[inode] = node;
              }
            vector<tri_tuple_type_local> elems_local;

            // test 1
            {
              unsigned edge_marks[3] = {1,1,0};
              double tri_coords[3][3] = {{0,0,0}, {1,0,0}, {0,1,0}};
              set_node_coords(eMesh, elem_nodes, tri_coords);
              Local_Tri3_Tri3_N::triangulate_face(eMesh, elem_nodes_vector, edge_marks, elems_local);

              // expected:
              vector<tri_tuple_type_local> elems_local_expected(3);
              elems_local_expected[0] = tri_tuple_type_local(0,3,4);
              elems_local_expected[1] = tri_tuple_type_local(3,1,4);
              elems_local_expected[2] = tri_tuple_type_local(0,4,2);
              
              std::cout << "test1: elems_local_expected= " << elems_local_expected << std::endl;
              std::cout << "test1: elems_local= " << elems_local << std::endl;

              STKUNIT_EXPECT_TRUE(in_set(elems_local_expected[0], elems_local));
              STKUNIT_EXPECT_TRUE(in_set(elems_local_expected[1], elems_local));
              STKUNIT_EXPECT_TRUE(in_set(elems_local_expected[2], elems_local));

              STKUNIT_EXPECT_TRUE(!in_set(elems_local_expected[0], elems_local, true));
              STKUNIT_EXPECT_TRUE(!in_set(elems_local_expected[1], elems_local, true));
              STKUNIT_EXPECT_TRUE(!in_set(elems_local_expected[2], elems_local, true));
            }

            // test2: same as test 1 but mirror image (emulating a face shared between two tets)
            {
              stk::mesh::Entity* node1 = elem_nodes_vector[1];
              elem_nodes_vector[1] = elem_nodes_vector[2];
              elem_nodes_vector[2] = node1;

              unsigned edge_marks[3] = {0,1,1};
              Local_Tri3_Tri3_N::triangulate_face(eMesh, elem_nodes_vector, edge_marks, elems_local);

              // expected:
              vector<tri_tuple_type_local> elems_local_expected(3);
              elems_local_expected[0] = tri_tuple_type_local(0,1,4);
              elems_local_expected[1] = tri_tuple_type_local(0,4,5);
              elems_local_expected[2] = tri_tuple_type_local(2,5,4);
              
              std::cout << "test2: elems_local_expected= " << elems_local_expected << std::endl;
              std::cout << "test2: elems_local= " << elems_local << std::endl;

              STKUNIT_EXPECT_TRUE(in_set(elems_local_expected[0], elems_local));
              STKUNIT_EXPECT_TRUE(in_set(elems_local_expected[1], elems_local));
              STKUNIT_EXPECT_TRUE(in_set(elems_local_expected[2], elems_local));

              STKUNIT_EXPECT_TRUE(!in_set(elems_local_expected[0], elems_local, true));
              STKUNIT_EXPECT_TRUE(!in_set(elems_local_expected[1], elems_local, true));
              STKUNIT_EXPECT_TRUE(!in_set(elems_local_expected[2], elems_local, true));
            }

            
        
            // end_demo
          }

      }



    } // namespace unit_tests
  } // namespace adapt
} // namespace stk


