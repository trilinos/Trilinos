
/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#include <stk_percept/Percept.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>

#include <stk_percept/Util.hpp>
#include <stk_percept/ExceptionWatch.hpp>
#include <stk_percept/GeometryVerifier.hpp>
#include <stk_percept/function/StringFunction.hpp>
#include <stk_percept/function/FieldFunction.hpp>
#include <stk_percept/function/ConstantFunction.hpp>
#include <stk_percept/PerceptMesh.hpp>

#include <stk_percept/mesh/mod/smoother/SpacingFieldUtil.hpp>

#include <stk_adapt/UniformRefinerPattern.hpp>
#include <stk_adapt/UniformRefiner.hpp>

#include <stk_percept/fixtures/SingleTetFixture.hpp>

#include <stk_mesh/fixtures/HexFixture.hpp>


#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <unit_tests/UnitTestSupport.hpp>

#include <stk_io/IossBridge.hpp>

#include <boost/lexical_cast.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>

#include <stk_percept/Percept.hpp>
#include <stk_percept/Util.hpp>
#include <stk_percept/ExceptionWatch.hpp>

#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>
#include <stk_percept/RunEnvironment.hpp>

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
#include <fstream>
#include <string>
#include <typeinfo>

#include <math.h>
#include <stk_util/parallel/Parallel.hpp>

namespace stk {
  namespace adapt {
    namespace regression_tests {

      static int print_infoLevel = 0;

      /// configuration: you can choose where to put the generated Exodus files (see variables input_files_loc, output_files_loc)
      /// The following defines where to put the input and output files created by this set of functions

#include "RegressionTestFileLoc.hpp"

#define EXTRA_PRINT 0

      static void save(PerceptMesh& eMesh, std::string filename)
      {
        eMesh.save_as(filename);
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      /// Creates meshes for use in later tests
      /// 1. Create hex mesh from a fixture and write it in Exodus format for use later.
      /// 2. Read the hex mesh and convert it to tet elements using stk_adapt/UniformRefiner, write it in Exodus format

      //STKUNIT_UNIT_TEST(unit_uniformRefiner, build_meshes)
      static void fixture_setup_0()
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );
        // start_demo_uniformRefiner_hex8_build
        {
          percept::PerceptMesh eMesh(3u);

          unsigned p_size = eMesh.get_parallel_size();

          // generate a 4x4x(4*p_size) mesh
          std::string gmesh_spec = std::string("4x4x")+toString(4*p_size)+std::string("|bbox:0,0,0,1,1,1");
          eMesh.new_mesh(percept::GMeshSpec(gmesh_spec));
          eMesh.commit();
          save(eMesh, input_files_loc+"hex_spc_fixture.e");

          // end_demo
        }

        // start_demo_uniformRefiner_hex8_build_1
        {
          percept::PerceptMesh eMesh(3u);

          //unsigned p_size = eMesh.get_parallel_size();
          eMesh.open(input_files_loc+"hex_spc_fixture.e");

          Hex8_Tet4_24 break_hex_to_tet(eMesh);

          int scalarDimension = 0; // a scalar
          stk::mesh::FieldBase* proc_rank_field = eMesh.add_field("proc_rank", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);

          eMesh.commit();

          UniformRefiner breaker(eMesh, break_hex_to_tet, proc_rank_field);
          breaker.doBreak();
          save(eMesh, input_files_loc+"tet_spc_fixture.e");
          save(eMesh, input_files_loc+"tet_from_hex_spc_fixture.e");
          // end_demo
        }
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      /// Creates meshes for use in later tests - quad meshes with and without sidesets

      //STKUNIT_UNIT_TEST(unit_uniformRefiner, quad4_quad4_4_test_1)
      static void fixture_setup_1()
      {
        EXCEPTWATCH;

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        //if (p_size == 1 || p_size == 3)
        if (p_size <= 3)
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
            eMesh.print_info("quad fixture", print_infoLevel);
            save(eMesh, input_files_loc+"quad_spc_fixture.e");
          }

        if (p_size <= 3)
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
            eMesh.print_info("quad fixture no sidesets", print_infoLevel);
            save(eMesh, input_files_loc+"quad_spc_fixture_no_sidesets.e");
          }
      }

      static void fixture_setup()
      {
        //std::cout << "tmp fixture_setup" << std::endl;
        static bool is_setup = false;
        if (is_setup) return;
        fixture_setup_0();
        fixture_setup_1();
        is_setup = true;
      }

#if defined ( STK_PERCEPT_HAS_GEOMETRY )

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Refine a triangle mesh with spacing
      STKUNIT_UNIT_TEST(regr_uniformRefiner, break_quad_to_quad_spc)
      {
        fixture_setup();
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size <= 3)
          {
            const unsigned n = 12;
            //const unsigned nx = n , ny = n , nz = p_size*n ;
            const unsigned nx = n , ny = n;

            bool createEdgeSets = true;
            percept::QuadFixture<double > fixture( pm , nx , ny, createEdgeSets);
            fixture.set_bounding_box(0,1,0,1);

            bool isCommitted = false;
            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, isCommitted);

            Quad4_Quad4_4 break_quad_to_quad_4(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.add_field("proc_rank", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);

            //stk::mesh::FieldBase* proc_rank_field_edge =
            eMesh.add_field("proc_rank_edge", eMesh.edge_rank(), scalarDimension);

            //             std::cout << "proc_rank_field rank= " << proc_rank_field->rank() << std::endl;
            //             std::cout << "proc_rank_field_edge rank= " << proc_rank_field_edge->rank() << std::endl;

            eMesh.add_spacing_fields();
            eMesh.set_respect_spacing(true);
            //fixture.meta_data.commit();
            eMesh.commit();

            fixture.generate_mesh();

            const std::vector<stk::mesh::Bucket*> & buckets = eMesh.get_bulk_data()->buckets( stk::mesh::MetaData::NODE_RANK );

            for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
              {
                  {
                    stk::mesh::Bucket & bucket = **k ;

                    const unsigned num_elements_in_bucket = bucket.size();
                
                    for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                      {
                        stk::mesh::Entity node = bucket[iEntity];

                        double * data = eMesh.field_data( *eMesh.get_coordinates_field() , node );
                        double iy = data[1]; // /double(nele);
                        iy = iy*iy;
                        data[1] = iy; // *double(nele);
                        
                        data[0] = data[0]*data[0]*data[0];
                      }
                  }
              }
#if 0
            Math::Matrix rmz = Math::rotationMatrix(2, 30);
            eMesh.transform_mesh(rmz);
#endif

            SpacingFieldUtil sfu(eMesh);
            sfu.compute_spacing_field();

            //eMesh.print_info("quad mesh", 5);
            eMesh.print_info("quad mesh");
            save(eMesh, output_files_loc+"quad_spc_fixture_quad3_spc_0.e");

            UniformRefiner breaker(eMesh, break_quad_to_quad_4, proc_rank_field);
            breaker.doBreak();

            //eMesh.print_info("quad mesh refined", 5);
            eMesh.print_info("quad mesh refined");
            save(eMesh, output_files_loc+"quad_spc_fixture_quad3_spc_1.e");

            // end_demo
          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Refine a triangle mesh with spacing
      STKUNIT_UNIT_TEST(regr_uniformRefiner, break_tri_to_tri_spc)
      {
        fixture_setup();
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size <= 3)
          {
            const unsigned n = 12;
            //const unsigned nx = n , ny = n , nz = p_size*n ;
            const unsigned nx = n , ny = n;

            bool createEdgeSets = true;
            percept::QuadFixture<double, shards::Triangle<3> > fixture( pm , nx , ny, createEdgeSets);
            fixture.set_bounding_box(0,1,0,1);

            bool isCommitted = false;
            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, isCommitted);

            Tri3_Tri3_4 break_tri_to_tri_4(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.add_field("proc_rank", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);

            //stk::mesh::FieldBase* proc_rank_field_edge =
            eMesh.add_field("proc_rank_edge", eMesh.edge_rank(), scalarDimension);

            //             std::cout << "proc_rank_field rank= " << proc_rank_field->rank() << std::endl;
            //             std::cout << "proc_rank_field_edge rank= " << proc_rank_field_edge->rank() << std::endl;

            eMesh.add_spacing_fields();
            eMesh.set_respect_spacing(true);
            //fixture.meta_data.commit();
            eMesh.commit();

            fixture.generate_mesh();

            const std::vector<stk::mesh::Bucket*> & buckets = eMesh.get_bulk_data()->buckets( stk::mesh::MetaData::NODE_RANK );

            for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
              {
                  {
                    stk::mesh::Bucket & bucket = **k ;

                    const unsigned num_elements_in_bucket = bucket.size();
                
                    for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                      {
                        stk::mesh::Entity node = bucket[iEntity];

                        double * data = eMesh.field_data( *eMesh.get_coordinates_field() , node );
                        double iy = data[1]; // /double(nele);
                        iy = iy*iy;
                        data[1] = iy; // *double(nele);
                      }
                  }
              }
#if 1
            Math::Matrix rmz = Math::rotationMatrix(2, 30);
            eMesh.transform_mesh(rmz);
#endif

            SpacingFieldUtil sfu(eMesh);
            sfu.compute_spacing_field();

            //eMesh.print_info("tri mesh", 5);
            eMesh.print_info("tri mesh");
            save(eMesh, output_files_loc+"quad_spc_fixture_tri3_spc_0.e");

            UniformRefiner breaker(eMesh, break_tri_to_tri_4, proc_rank_field);
            breaker.doBreak();

            //eMesh.print_info("tri mesh refined", 5);
            eMesh.print_info("tri mesh refined");
            save(eMesh, output_files_loc+"quad_spc_fixture_tri3_spc_1.e");

            // end_demo
          }

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      // A cube with an indented bump on the bottom

      STKUNIT_UNIT_TEST(regr_uniformRefiner, hex_4)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        MPI_Barrier( MPI_COMM_WORLD );

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size <= 3)
          {
            unsigned n = 12;
            std::string gmesh_spec = toString(n)+"x"+toString(n)+"x"+toString(n)+std::string("|bbox:0,0,0,1,1,1|sideset:xXyYzZ");
            PerceptMesh eMesh(3);
            eMesh.new_mesh(percept::GMeshSpec(gmesh_spec));
            Hex8_Hex8_8 break_hex_to_hex(eMesh);
            eMesh.add_spacing_fields();
            eMesh.set_respect_spacing(true);
            eMesh.commit();

            const std::vector<stk::mesh::Bucket*> & buckets = eMesh.get_bulk_data()->buckets( stk::mesh::MetaData::NODE_RANK );

            // cluster the mesh towards the bump
            for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
              {
                  {
                    stk::mesh::Bucket & bucket = **k ;

                    const unsigned num_elements_in_bucket = bucket.size();
                
                    for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                      {
                        stk::mesh::Entity entity = bucket[iEntity];

                        double * data = eMesh.field_data( *eMesh.get_coordinates_field() , entity );
                        data[2] = data[2]*data[2];
                      }
                  }
              }
            SpacingFieldUtil sfu(eMesh);
            sfu.compute_spacing_field();

            save(eMesh, output_files_loc+"hex_4_spc.0.e");

            UniformRefiner breaker(eMesh, break_hex_to_hex, 0);
            breaker.doBreak();

            save(eMesh, output_files_loc+"hex_4_spc.1.e");

          }
      }

#endif

    } // namespace regression_tests
  } // namespace adapt
} // namespace stk

