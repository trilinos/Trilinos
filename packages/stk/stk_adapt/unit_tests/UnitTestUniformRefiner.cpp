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
#include <unit_tests/TestLocalRefinerTri_N.hpp>
#include <unit_tests/TestLocalRefinerTri_N_1.hpp>

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

          unsigned p_size = eMesh.getParallelSize();

          // generate a 4x4x(4*p_size) mesh
          std::string gmesh_spec = std::string("4x4x")+toString(4*p_size)+std::string("|bbox:0,0,0,1,1,1");
          eMesh.newMesh(percept::PerceptMesh::GMeshSpec(gmesh_spec));
          eMesh.commit();
          save_or_diff(eMesh, input_files_loc+"hex_fixture.e");

          // end_demo
        }

        // start_demo_uniformRefiner_hex8_build_1
        {
          percept::PerceptMesh eMesh(3u);

          //unsigned p_size = eMesh.getParallelSize();
          eMesh.open(input_files_loc+"hex_fixture.e");

          Hex8_Tet4_24 break_hex_to_tet(eMesh);

          int scalarDimension = 0; // a scalar
          stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

          eMesh.commit();

          UniformRefiner breaker(eMesh, break_hex_to_tet, proc_rank_field);
          breaker.doBreak();
          save_or_diff(eMesh, input_files_loc+"tet_fixture.e");
          save_or_diff(eMesh, input_files_loc+"tet_from_hex_fixture.e");
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
            eMesh.printInfo("quad fixture", 2);
            save_or_diff(eMesh, input_files_loc+"quad_fixture.e");
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
            eMesh.printInfo("quad fixture no sidesets", 2);
            save_or_diff(eMesh, input_files_loc+"quad_fixture_no_sidesets.e");
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

      //=====================================================================================================================================================================================================
      //=====================================================================================================================================================================================================
      //=====================================================================================================================================================================================================

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

#if 1

      /// Refine a quad mesh

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, break_quad_to_quad_sierra_1_test)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size <= 2)
          {
            const unsigned n = 2;
            //const unsigned nx = n , ny = n , nz = p_size*n ;
            const unsigned nx = 1 , ny = n;

            bool createEdgeSets = false;
            percept::QuadFixture<double > fixture( pm , nx , ny, createEdgeSets);

            bool isCommitted = false;
            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, isCommitted);

            Quad4_Quad4_4 break_quad_to_quad_4(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            //stk::mesh::FieldBase* proc_rank_field_edge =
            eMesh.addField("proc_rank_edge", eMesh.edge_rank(), scalarDimension);

            //             std::cout << "proc_rank_field rank= " << proc_rank_field->rank() << std::endl;
            //             std::cout << "proc_rank_field_edge rank= " << proc_rank_field_edge->rank() << std::endl;

            //fixture.meta_data.commit();
            eMesh.commit();

            fixture.generate_mesh();

            UniformRefiner breaker(eMesh, break_quad_to_quad_4, proc_rank_field);
            breaker.setRemoveOldElements(false);
            breaker.doBreak();
            MPI_Barrier( MPI_COMM_WORLD );

            eMesh.dumpElementsCompact();

            MPI_Barrier( MPI_COMM_WORLD );
            //exit(123);

            
            // end_demo
          }

      }

      /// Refine a triangle mesh

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, break_tri_to_tri_sierra_1_test)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size <= 2)
          {
            const unsigned n = 2;
            //const unsigned nx = n , ny = n , nz = p_size*n ;
            const unsigned nx = n , ny = n;

            bool createEdgeSets = false;
            percept::QuadFixture<double, shards::Triangle<3> > fixture( pm , nx , ny, createEdgeSets);

            bool isCommitted = false;
            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, isCommitted);

            Tri3_Tri3_4 break_tri_to_tri_4(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            //stk::mesh::FieldBase* proc_rank_field_edge =
            eMesh.addField("proc_rank_edge", eMesh.edge_rank(), scalarDimension);

            //             std::cout << "proc_rank_field rank= " << proc_rank_field->rank() << std::endl;
            //             std::cout << "proc_rank_field_edge rank= " << proc_rank_field_edge->rank() << std::endl;

            //fixture.meta_data.commit();
            eMesh.commit();

            fixture.generate_mesh();

            UniformRefiner breaker(eMesh, break_tri_to_tri_4, proc_rank_field);
            breaker.setRemoveOldElements(false);
            breaker.doBreak();
            eMesh.dumpElementsCompact();

            MPI_Barrier( MPI_COMM_WORLD );
            //exit(123);

            
            // end_demo
          }

      }

#endif

#if 1
      //=============================================================================
      //=============================================================================
      //=============================================================================

      STKUNIT_UNIT_TEST(unit_localRefiner, break_tri_to_tri_N_1)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size <= 3)
          {
            std::string post_fix[4] = {"np0", "np1", "np2", "np3"};

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

            //eMesh.printInfo("local tri mesh",2);
            save_or_diff(eMesh, output_files_loc+"local_tri_N_1_0_"+post_fix[p_size]+".e");

            TestLocalRefinerTri_N_1 breaker(eMesh, break_tri_to_tri_N, proc_rank_field);
            breaker.setRemoveOldElements(false);
            breaker.doBreak();

            eMesh.dumpElementsCompact();

            //eMesh.printInfo("local tri mesh refined", 2);
            save_or_diff(eMesh, output_files_loc+"local_tri_N_1_1_"+post_fix[p_size]+".e");

            MPI_Barrier( MPI_COMM_WORLD );
            //exit(123);

            //breaker.unrefineAll();
            ElementUnrefineCollection elements_to_unref = breaker.buildTestUnrefList();
            breaker.unrefineTheseElements(elements_to_unref);

            // FIXME
            eMesh.saveAs( output_files_loc+"local_tri_N_1_1_unref_"+post_fix[p_size]+".e");
            //save_or_diff(eMesh, output_files_loc+"local_tri_N_1_1_unref_"+post_fix[p_size]+".e");

            // end_demo
          }

      }
#endif

#if 0
      //=============================================================================
      //=============================================================================
      //=============================================================================


      STKUNIT_UNIT_TEST(unit_tmp, break_tri_to_tri_N)
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

            eMesh.printInfo("local tri mesh",2);
            save_or_diff(eMesh, output_files_loc+"local_tri_N_0.e");

            TestLocalRefinerTri_N breaker(eMesh, break_tri_to_tri_N, proc_rank_field);
            breaker.setRemoveOldElements(false);
            breaker.doBreak();

            eMesh.printInfo("local tri mesh refined", 2);
            //eMesh.dumpElements();
            save_or_diff(eMesh, output_files_loc+"local_tri_N_1.e");

            //breaker.unrefineAll();
            ElementUnrefineCollection elements_to_unref = breaker.buildTestUnrefList();
            breaker.unrefineTheseElements(elements_to_unref);

            // FIXME
            eMesh.saveAs( output_files_loc+"local_tri_N_1_unref.e");
            //save_or_diff(eMesh, output_files_loc+"local_tri_N_1_unref.e");

            //exit(123);
            // end_demo
          }

      }
#endif

      //=============================================================================
      //=============================================================================
      //=============================================================================

      /// Refine quad elements
      /// uses the Sierra-ported tables from framework/{element,mesh_modification}
      STKUNIT_UNIT_TEST(unit_uniformRefiner, break_quad_to_quad_sierra)
      {
        fixture_setup();
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

            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, false);

            UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4, SierraPort > break_quad_to_quad_4(eMesh);
            // FIXME
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            //fixture.meta_data.commit();
            eMesh.commit();

            fixture.generate_mesh();

            eMesh.printInfo("quad mesh");

            UniformRefiner breaker(eMesh, break_quad_to_quad_4, proc_rank_field);
            //breaker.setRemoveOldElements(false);
            breaker.doBreak();

            //!!eMesh, "./square_quad4_ref_sierra_out.e");
            // end_demo
          }

      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      /// Create a triangle mesh using the QuadFixture with the option of breaking the quads into triangles
      /// Refine the triangle mesh, write the results.

      STKUNIT_UNIT_TEST(unit_uniformRefiner, break_tri_to_tri_sierra)
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

            percept::QuadFixture<double, shards::Triangle<3> > fixture( pm , nx , ny, true);

            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, false);

            //             UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4, SierraPort > break_tri_to_tri_4(eMesh);
            //             // FIXME
            //             int scalarDimension = 0; // a scalar
            //             FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            //fixture.meta_data.commit();
            eMesh.commit();

            fixture.generate_mesh();

            eMesh.printInfo("tri mesh");

            //             UniformRefiner breaker(eMesh, break_tri_to_tri_4, proc_rank_field);
            //             breaker.doBreak();

            save_or_diff(eMesh, input_files_loc+"quad_fixture_tri3.e");
            // end_demo
          }

      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      /// Refine a hex8 mesh
      STKUNIT_UNIT_TEST(unit_uniformRefiner, hex8_hex8_8_1)
      {
        fixture_setup();
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        // start_demo_uniformRefiner_hex8_hex8_8_1

        percept::PerceptMesh eMesh(3u);

        unsigned p_size = eMesh.getParallelSize();

        // generate a 4x4x(4*p_size) mesh
        std::string gmesh_spec = std::string("4x4x")+toString(4*p_size)+std::string("|bbox:0,0,0,1,1,1");
        eMesh.newMesh(percept::PerceptMesh::GMeshSpec(gmesh_spec));
        //eMesh.commit();
        //eMesh.reopen();

        Hex8_Hex8_8 break_hex_to_hex(eMesh);

        int scalarDimension = 0; // a scalar
        stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

        eMesh.commit();

        UniformRefiner breaker(eMesh, break_hex_to_hex, proc_rank_field);
        breaker.doBreak();
        // end_demo
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      /// Create and write a wedge mesh using the WedgeFixture

      STKUNIT_UNIT_TEST(unit_uniformRefiner, wedge6_1)
      {
        fixture_setup();
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        // start_demo_uniformRefiner_wedge6_1

        percept::PerceptMesh eMesh(3u);

        unsigned p_size = eMesh.getParallelSize();
        if (p_size == 1)  // this fixture only works in serial mode
          {
            percept::WedgeFixture wedgeFixture;
            wedgeFixture.createMesh(MPI_COMM_WORLD,
                                    4, 3, 2,
                                    0, 1,
                                    0, 1,
                                    0, 1,
                                    std::string(input_files_loc+"swept-wedge_0.e") );
          }
      }


      //=============================================================================================================================================================
      //=============================================================================================================================================================
      //=============================================================================================================================================================

      /// Create a Beam mesh and enrich it

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, beam_enrich)
      {
        fixture_setup();
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size(pm);

        if (p_size <= 1)
          {
            // create the mesh
            {

              stk::percept::BeamFixture mesh(pm, false);
              stk::io::put_io_part_attribute(  mesh.m_block_beam );
              mesh.m_metaData.commit();
              mesh.populate();

              bool isCommitted = true;
              percept::PerceptMesh em1(&mesh.m_metaData, &mesh.m_bulkData, isCommitted);
              save_or_diff(em1, input_files_loc+"beam_enrich_0.e");

            }

            // enrich
            {
              stk::percept::PerceptMesh eMesh(3u);
              eMesh.open(input_files_loc+"beam_enrich_0.e");
              //URP_Heterogeneous_3D break_pattern(eMesh);
              Beam2_Beam3_1 break_pattern(eMesh);
              int scalarDimension = 0; // a scalar
              stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
              eMesh.commit();

              save_or_diff(eMesh, output_files_loc+"beam_enrich_0.e");

              eMesh.printInfo("beam", 2);

              UniformRefiner breaker(eMesh, break_pattern, proc_rank_field);
              //breaker.setRemoveOldElements(false);
              breaker.setIgnoreSideSets(true);
              breaker.doBreak();

              save_or_diff(eMesh, output_files_loc+"beam_enrich_1.e");

            }
          }
      }
      //=============================================================================
      //=============================================================================
      //=============================================================================

      /// Create a beam mesh and refine it

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, beam_refine)
      {
        fixture_setup();
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size(pm);

        if (p_size <= 1)
          {
            // create the mesh
            {

              stk::percept::BeamFixture mesh(pm, false);
              stk::io::put_io_part_attribute(  mesh.m_block_beam );
              mesh.m_metaData.commit();
              mesh.populate();

              bool isCommitted = true;
              percept::PerceptMesh em1(&mesh.m_metaData, &mesh.m_bulkData, isCommitted);
              save_or_diff(em1, input_files_loc+"beam_0.e");

            }

            // refine
            {
              stk::percept::PerceptMesh eMesh(3u);
              eMesh.open(input_files_loc+"beam_0.e");
              //URP_Heterogeneous_3D break_pattern(eMesh);
              Beam2_Beam2_2 break_pattern(eMesh);
              int scalarDimension = 0; // a scalar
              stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
              eMesh.commit();

              save_or_diff(eMesh, output_files_loc+"beam_0.e");

              eMesh.printInfo("beam", 2);

              UniformRefiner breaker(eMesh, break_pattern, proc_rank_field);
              //breaker.setRemoveOldElements(false);
              breaker.setIgnoreSideSets(true);
              breaker.doBreak();

              save_or_diff(eMesh, output_files_loc+"beam_1.e");

            }
          }
      }


      //======================================================================================================================
      //======================================================================================================================
      //===================== Table generation
      //======================================================================================================================

      /// This code generates C++ tables used by stk_adapt - tables contain node numbering, parametric coordinates consistent
      ///    with Intrepid, and related information needed by UniformRefiner.  The generated code should be compared with
      ///    and merged into <stk_adapt/sierra_element/GeneratedRefinementTable.hpp> as appropriate if there is a change
      ///    in the tables in that package, or an additional element type is added. 
      /** This comment is from the generated code and tells how to bootstrap this process.
       * Bootstrapping this file: to create this file, run the regression test RegressionTestUniformRefiner.cpp :: generate_tables after putting in
       *   a dummy entry in ./sierra_element/GeneratedRefinementTable.hpp.  The run will produce a local file, generated_refinement_tables.hpp 
       *   which can be checked against the gold copy of GeneratedRefinementTable.hpp, then copied over it.  Add a call below to generate the 
       *   actual new table data. 
       */


      STKUNIT_UNIT_TEST(unit1_uniformRefiner, generate_tables)
      {
        fixture_setup();
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );
        Elem::StdMeshObjTopologies::bootstrap();

        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1)
          {
            std::ofstream file("./generated_refinement_tables.hpp");

            file << "#ifndef STK_ADAPT_GENERATED_REFINEMENT_TABLES_HPP" << std::endl;
            file << "#define STK_ADAPT_GENERATED_REFINEMENT_TABLES_HPP" << std::endl;

            file <<
              "/**  New ref topo info \n"
              "*  ------------------\n"
              "*\n"
              "*  {Ord, Rnk-assoc, Ord-rnk-assoc, Ord-node-on-subcell, num-rnk-assoc, param-coord}\n"
              "*\n"
              "*   struct RefinementTopologyExtraEntry\n"
              "*   {\n"
              "*     unsigned ordinal_of_node;               // ordinal of node in the total list of nodes - corresponds to the shards node ordinal\n"
              "*     unsigned rank_of_subcell;               // rank of the subcell this node is associated with                                   \n"
              "*     unsigned ordinal_of_subcell;            // ordinal of the subcell in the shards numbering (e.g. edge # 3)\n"
              "*     unsigned ordinal_of_node_on_subcell;    // ordinal of the node on the subcell (whcih node it is on a subcell that has multiple nodes)\n"
              "*     unsigned num_nodes_on_subcell;          // how many nodes exist on the subcell                                                       \n"
              "*     double parametric_coordinates[3];\n"
              "*   };\n"
              "*       \n"
              "* Bootstrapping this file: to create this file, run the regression test RegressionTestUniformRefiner.cpp :: generate_tables after putting in\n"
              "*   a dummy entry in ./sierra_element/GeneratedRefinementTable.hpp.  The run will produce a local file, generated_refinement_tables.hpp \n"
              "*   which can be checked against the gold copy of GeneratedRefinementTable.hpp, then copied over it.  Add a call below to generate the \n"
              "*   actual new table data. \n"
              "*/\n\n"
                 << std::endl;

            // FIXME
#if !(defined(__PGI) && defined(USE_PGI_7_1_COMPILER_BUG_WORKAROUND))

            Line2_Line2_2            :: printRefinementTopoX_Table(file);

            Beam2_Beam2_2            :: printRefinementTopoX_Table(file);

            ShellLine2_ShellLine2_2  :: printRefinementTopoX_Table(file);
            ShellLine3_ShellLine3_2  :: printRefinementTopoX_Table(file);
            Quad4_Quad4_4            :: printRefinementTopoX_Table(file);
            Tri3_Tri3_4              :: printRefinementTopoX_Table(file);
            ShellTri3_ShellTri3_4    :: printRefinementTopoX_Table(file);
            ShellTri6_ShellTri6_4    :: printRefinementTopoX_Table(file);
            ShellQuad4_ShellQuad4_4  :: printRefinementTopoX_Table(file);
            ShellQuad8_ShellQuad8_4  :: printRefinementTopoX_Table(file);
            Tet4_Tet4_8              :: printRefinementTopoX_Table(file);
            Hex8_Hex8_8              :: printRefinementTopoX_Table(file);
            Wedge6_Wedge6_8          :: printRefinementTopoX_Table(file);
            Wedge15_Wedge15_8        :: printRefinementTopoX_Table(file);

            // Not supported by Sierra
            // Wedge18_Wedge18_8        :: printRefinementTopoX_Table(file);

            Line3_Line3_2            :: printRefinementTopoX_Table(file);
            Beam3_Beam3_2            :: printRefinementTopoX_Table(file);

            Tri6_Tri6_4              :: printRefinementTopoX_Table(file);
            Quad8_Quad8_4            :: printRefinementTopoX_Table(file);
            Quad9_Quad9_4            :: printRefinementTopoX_Table(file);
            Hex27_Hex27_8            :: printRefinementTopoX_Table(file);
            Hex20_Hex20_8            :: printRefinementTopoX_Table(file);
            Tet10_Tet10_8            :: printRefinementTopoX_Table(file);
#endif
            file << "#endif" << std::endl;
          }

      }


      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Code to generate Dot/Graphviz files representing the topology of element refinement based on the internal tables
      
      static void output_draw(std::string filename, std::string toFile)
      {
        std::ofstream file(filename.c_str());
        file << toFile;
      }

      STKUNIT_UNIT_TEST(unit_uniformRefiner, draw1)
      {
        fixture_setup();
        //std::cout << Quad4_Quad4_4::draw() << std::endl;
        std::string dir = "./";
        output_draw(dir+"quad4.dot", Quad4_Quad4_4::draw(true) );
        output_draw(dir+"tet4.dot",  Tet4_Tet4_8::draw() );
        output_draw(dir+"hex8.dot",  Hex8_Hex8_8::draw(true) );
        output_draw(dir+"hex27.dot",  Hex27_Hex27_8::draw(true, true) );
        output_draw(dir+"hex20.dot",  Hex20_Hex20_8::draw(true, true) );
        output_draw(dir+"wedge6.dot",  Wedge6_Wedge6_8::draw() );

        output_draw(dir+"quad9.dot", Quad9_Quad9_4::draw(true, true));
      }

      STKUNIT_UNIT_TEST(unit_uniformRefiner, draw)
      {
        fixture_setup();
        //std::cout << Quad4_Quad4_4::draw() << std::endl;
        std::string dir = "./";
        output_draw(dir+"quad4.dot", Quad4_Quad4_4::draw(true) );
        output_draw(dir+"tet4.dot",  Tet4_Tet4_8::draw() );
        output_draw(dir+"hex8.dot",  Hex8_Hex8_8::draw(true) );
        output_draw(dir+"wedge6.dot",  Wedge6_Wedge6_8::draw() );

        output_draw(dir+"quad9.dot", Quad9_Quad9_4::draw(true));

        // refine
#if 0
        std::cout << Line2_Line2_2::draw() << std::endl;
        std::cout << Tri3_Tri3_4::draw() << std::endl;
        std::cout << Tet4_Tet4_8::draw() << std::endl;
        std::cout << Hex8_Hex8_8::draw() << std::endl;
#endif

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Convert a quad mesh to triangles with 6 triangles per quad

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, break_quad_to_tri_6)
      {
        fixture_setup();
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            // start_demo_uniformRefiner_break_quad_to_tri_6
            percept::PerceptMesh eMesh(2u);
            eMesh.open(input_files_loc+"quad_fixture.e");

            typedef  UniformRefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>, 6 > Quad4_Tri3_6;

            UniformRefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>, 6 > break_quad_to_tri_6(eMesh);

            int scalarDimension = 0; // a scalar
            //         int vectorDimension = 3;

            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();

            eMesh.printInfo("quad mesh");

            //UniformRefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>, 6 > break_quad_to_tri_6;
            UniformRefiner breaker(eMesh, break_quad_to_tri_6, proc_rank_field);
            breaker.setIgnoreSideSets(true);
            breaker.setRemoveOldElements(false);
            breaker.doBreak();

            save_or_diff(eMesh, output_files_loc+"square_quad4_out.e");
            // end_demo
          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Convert a quad mesh to triangles with 4 triangles per quad

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, break_quad_to_tri_4)
      {
        fixture_setup();
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {

            // start_demo_uniformRefiner_break_quad_to_tri_4
            percept::PerceptMesh eMesh(2u);
            eMesh.open(input_files_loc+"quad_fixture.e");

            UniformRefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>, 4, Specialization > break_quad_to_tri_4(eMesh);

            int scalarDimension = 0; // a scalar

            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();

            eMesh.printInfo("quad mesh");

            //UniformRefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>, 6 > break_quad_to_tri_6;
            UniformRefiner breaker(eMesh, break_quad_to_tri_4, proc_rank_field);
            breaker.setRemoveOldElements(false);
            breaker.setIgnoreSideSets(true);
            breaker.doBreak();

            save_or_diff(eMesh, output_files_loc+"square_quad4_tri3_4_out.e");
            // end_demo
          }
      }


      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Refine a quad mesh using the "standalone" refinement pattern:
      ///      UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4 >
      /// This pattern is an example (like the convert-type patterns) showing how to write a new pattern with no dependencies
      //     on other (say tabular) data/info.

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, break_quad_to_quad)
      {
        fixture_setup();
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            // start_demo_uniformRefiner_break_quad_to_quad
            percept::PerceptMesh eMesh(2u);
            eMesh.open(input_files_loc+"quad_fixture_no_sidesets.e");

            UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4 > break_quad_to_quad_4(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
            eMesh.commit();

            eMesh.printInfo("quad mesh");

            UniformRefiner breaker(eMesh, break_quad_to_quad_4, proc_rank_field);
            //breaker.setRemoveOldElements(false);
            breaker.setIgnoreSideSets(true);
            breaker.doBreak();

            save_or_diff(eMesh, output_files_loc+"square_quad4_ref_out.e");
            // end_demo
          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Refine a quad mesh
      /// uses the Sierra-ported tables from framework/{element,mesh_modification}

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, break_quad_to_quad_sierra)
      {
        fixture_setup();
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            // start_demo_uniformRefiner_break_quad_to_quad_sierra
            percept::PerceptMesh eMesh(2u);
            eMesh.open(input_files_loc+"quad_fixture.e");

            UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4, SierraPort > break_quad_to_quad_4(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
            eMesh.commit();

            eMesh.printInfo("quad mesh");

            UniformRefiner breaker(eMesh, break_quad_to_quad_4, proc_rank_field);
            //breaker.setRemoveOldElements(false);
            breaker.doBreak();

            save_or_diff(eMesh, output_files_loc+"square_quad4_ref_sierra_out.e");
            // end_demo
          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Refine a quad mesh with sidesets
      /// uses the Sierra-ported tables from framework/{element,mesh_modification}

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, break_quad_to_quad_sierra_sidesets)
      {
        fixture_setup();
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );

        if (p_size == 1 || p_size == 2)
          {
            // start_demo_uniformRefiner_break_quad_to_quad_sierra_sidesets
            percept::PerceptMesh eMesh(2u);
            eMesh.open(input_files_loc+"quad_fixture.e");

            UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4, SierraPort > break_quad_to_quad_4(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
            eMesh.commit();

            eMesh.printInfo("quad mesh");

            UniformRefiner breaker(eMesh, break_quad_to_quad_4, proc_rank_field);
            //breaker.setRemoveOldElements(false);
            breaker.doBreak();

            eMesh.printInfo("after refinement break_quad_to_quad_sierra_sidesets");

            save_or_diff(eMesh, output_files_loc+"quad_sidesets_sierra_out.e");
          }
        // end_demo
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Convert a hex mesh to tets using 24 tets per hex

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, hex8_tet4_24_1)
      {
        fixture_setup();
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        // start_demo_uniformRefiner_hex8_tet4_24_1

        percept::PerceptMesh eMesh(3u);

        unsigned p_size = eMesh.getParallelSize();
        //unsigned p_rank = eMesh.getRank();
        Util::setRank(eMesh.getRank());

        std::string gmesh_spec = std::string("1x1x")+toString(p_size)+std::string("|bbox:0,0,0,1,1,"+toString(p_size) );
        eMesh.newMesh(percept::PerceptMesh::GMeshSpec(gmesh_spec));

        UniformRefinerPattern<shards::Hexahedron<8>, shards::Tetrahedron<4>, 24 > break_hex_to_tet(eMesh);

        int scalarDimension = 0; // a scalar
        //         int vectorDimension = 3;

        stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
        //         eMesh.addField("velocity", mesh::Node, vectorDimension);
        //         eMesh.addField("element_volume", eMesh.element_rank(), scalarDimension);

        eMesh.commit();
        eMesh.printInfo();
        save_or_diff(eMesh, std::string(output_files_loc+"")+std::string("hex_tet_24_cube1x1x")+toString(p_size)+std::string("-orig.e"));

        UniformRefiner breaker(eMesh, break_hex_to_tet, proc_rank_field);
        breaker.setRemoveOldElements(true);

        breaker.doBreak();

        save_or_diff(eMesh, std::string(output_files_loc+"")+std::string("hex_tet_24_cube1x1x")+toString(p_size)+std::string(".e"));

        // end_demo

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Convert a hex mesh using 6 tets per hex

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, hex8_tet4_6_12_1)
      {
        fixture_setup();
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        // start_demo_uniformRefiner_hex8_tet4_6_12_1

        percept::PerceptMesh eMesh(3u);

        unsigned p_size = eMesh.getParallelSize();

        std::string gmesh_spec = std::string("1x1x")+toString(p_size)+std::string("|bbox:0,0,0,1,1,"+toString(p_size) );
        eMesh.newMesh(percept::PerceptMesh::GMeshSpec(gmesh_spec));

        Hex8_Tet4_6_12 break_hex_to_tet(eMesh);

        int scalarDimension = 0; // a scalar
        //         int vectorDimension = 3;

        stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
        //         eMesh.addField("velocity", mesh::Node, vectorDimension);
        //         eMesh.addField("element_volume", eMesh.element_rank(), scalarDimension);

        eMesh.commit();
        eMesh.printInfo();
        save_or_diff(eMesh, std::string(output_files_loc+"")+std::string("hex_tet_6_12_cube1x1x")+toString(p_size)+std::string("-orig.e"));

        UniformRefiner breaker(eMesh, break_hex_to_tet, proc_rank_field);
        breaker.setRemoveOldElements(true);
        //breaker.setIgnoreSideSets(true);

        breaker.doBreak();

        save_or_diff(eMesh, std::string(output_files_loc+"")+std::string("hex_tet_6_12_cube1x1x")+toString(p_size)+std::string(".e"));

        // end_demo
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Convert a hex mesh using 6 tets per hex

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, hex8_tet4_6_12_2)
      {
        fixture_setup();

        EXCEPTWATCH;

        // start_demo_uniformRefiner_hex8_tet4_6_12_2
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            percept::PerceptMesh eMesh(3u);

            eMesh.open(input_files_loc+"hex_fixture.e");

            Hex8_Tet4_6_12 break_hex_to_tet(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();
            //eMesh.printInfo("test",2);
            eMesh.printInfo();
            save_or_diff(eMesh, output_files_loc+"hex8_tet4_6_12_0.e");

            UniformRefiner breaker(eMesh, break_hex_to_tet, proc_rank_field);
            //breaker.setIgnoreSideSets(true);
            breaker.doBreak();

            save_or_diff(eMesh, output_files_loc+"hex8_tet4_6_12_1.e");

            // end_demo
          }

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Refine a quad mesh

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, quad4_quad4_4_test_1)
      {
        fixture_setup();
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
            eMesh.printInfo("quad fixture");
            save_or_diff(eMesh, output_files_loc+"quad_fixture_test_1.e");
          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Refine a quad mesh; test the multiple refinement feature

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, break_quad_to_quad_sierra_1)
      {
        fixture_setup();
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        bool doGenSideSets = true;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        //if (p_size == 1 || p_size == 3)
        if (p_size <= 3)
          {
            const unsigned n = 12;
            //const unsigned nx = n , ny = n , nz = p_size*n ;
            const unsigned nx = n , ny = n;

            percept::QuadFixture<double> fixture( pm , nx , ny, doGenSideSets);

            // Adopt the meta/bulk data
            bool isCommitted = false;
            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, isCommitted);

            eMesh.commit();

            fixture.generate_mesh();

            //eMesh.printInfo("quad mesh");

            save_or_diff(eMesh, output_files_loc+"quad_fixture_0.e");
            save_or_diff(eMesh, input_files_loc+"quad_fixture_0.e");
            eMesh.close();

            for (int iBreak = 0; iBreak < 2; iBreak++)
              {
                std::cout << "\n\n\n ================ tmp Refine Pass = " << iBreak << std::endl;

                percept::PerceptMesh eMesh1(2);
                std::string fileName = std::string(input_files_loc+"quad_fixture_")+toString(iBreak)+std::string(".e");
                eMesh1.open(fileName);
                UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4, SierraPort > break_quad_to_quad_4(eMesh1);
                int scalarDimension = 0; // a scalar
                stk::mesh::FieldBase* proc_rank_field = eMesh1.addField("proc_rank", eMesh.element_rank(), scalarDimension);
                eMesh1.commit();

                //                 if (iBreak != 0)
                //                   proc_rank_field = eMesh1.getField("proc_rank");

                UniformRefiner breaker(eMesh1, break_quad_to_quad_4, proc_rank_field);

                //breaker.setRemoveOldElements(false);
                breaker.doBreak();
                std::string fileName1 = std::string(output_files_loc+"quad_fixture_")+toString(iBreak+1)+std::string(".e");
                std::string fileName2 = std::string(input_files_loc+"quad_fixture_")+toString(iBreak+1)+std::string(".e");

                //eMesh1.printInfo("quad_fixture_1.e");

                save_or_diff(eMesh1, fileName2);
                save_or_diff(eMesh1, fileName1);
                eMesh1.close();

                if (0 && iBreak==0)
                  {
                    percept::PerceptMesh e1(2);
                    std::cout << "\n\n\n ================ tmp eMesh1.openReadOnly(quad_fixture_1.e) \n\n\n " << std::endl;
                    e1.openReadOnly(input_files_loc+"quad_fixture_1.e");
                    e1.printInfo("quad_fixture_1_read.e");
                    e1.close();
                  }
              }

            // end_demo
          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Refine a quad mesh; test the multiple refinement feature

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, break_quad_to_quad_sierra_2)
      {
        fixture_setup();
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        //if (p_size == 1 || p_size == 3)
        if (p_size <= 3)
          {
            const unsigned n = 12;
            //const unsigned nx = n , ny = n , nz = p_size*n ;
            const unsigned nx = n , ny = n;

            bool doGenSideSets = true;
            percept::QuadFixture<double> fixture( pm , nx , ny, doGenSideSets);

            // Adopt the meta/bulk data
            bool isCommitted = false;
            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, isCommitted);

            eMesh.commit();

            fixture.generate_mesh();

            //eMesh.printInfo("quad mesh");

            save_or_diff(eMesh, output_files_loc+"quad_fixture_mbreak_0.e");
            save_or_diff(eMesh, input_files_loc+"quad_fixture_mbreak_0.e");
            eMesh.close();


            percept::PerceptMesh eMesh1(2);
            std::string fileName = std::string(input_files_loc+"quad_fixture_mbreak_0.e");
            eMesh1.open(fileName);
            UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4, SierraPort > break_quad_to_quad_4(eMesh1);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh1.addField("proc_rank", eMesh.element_rank(), scalarDimension);
            eMesh1.commit();

            UniformRefiner breaker(eMesh1, break_quad_to_quad_4, proc_rank_field);

            for (int iBreak = 0; iBreak < 2; iBreak++)
              {
                std::cout << "\n\n\n ================ tmp Refine Pass = " << iBreak << std::endl;

                breaker.doBreak();
                std::string fileName1 = std::string(output_files_loc+"quad_fixture_mbreak_")+toString(iBreak+1)+std::string(".e");
                save_or_diff(eMesh1, fileName1);
              }

            // end_demo
          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Enrich a quad mesh (convert linear Quad4 elements to quadratic Quad9)

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, break_quad4_to_quad9)
      {
        fixture_setup();
        EXCEPTWATCH;


        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        //if (p_size == 1 || p_size == 3)
        if (p_size <= 3)
          {
            const unsigned n = 12;
            //const unsigned nx = n , ny = n , nz = p_size*n ;
            const unsigned nx = n , ny = n;

            bool doGenSideSets = true;
            percept::QuadFixture<double> fixture( pm , nx , ny, doGenSideSets);

            // Adopt the meta/bulk data
            bool isCommitted = false;
            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, isCommitted);

            Quad4_Quad9_1 break_quad4_to_quad9_1(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();

            fixture.generate_mesh();

            save_or_diff(eMesh, output_files_loc+"quad_fixture_quad9_0.e");

            UniformRefiner breaker(eMesh, break_quad4_to_quad9_1, proc_rank_field);

            breaker.doBreak();
            save_or_diff(eMesh, output_files_loc+"quad_fixture_quad9_1.e");

            // end_demo
          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Enrich a quad mesh (convert linear Quad4 elements to serendepity Quad8)

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, break_quad4_to_quad8)
      {
        fixture_setup();
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size <= 3)
          {
            const unsigned n = 12;
            const unsigned nx = n , ny = n;

            bool doGenSideSets = true;
            percept::QuadFixture<double> fixture( pm , nx , ny, doGenSideSets);

            // Adopt the meta/bulk data
            bool isCommitted = false;
            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, isCommitted);

            Quad4_Quad8_1 break_quad4_to_quad8_1(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();

            fixture.generate_mesh();

            save_or_diff(eMesh, output_files_loc+"quad_fixture_quad8_0.e");

            UniformRefiner breaker(eMesh, break_quad4_to_quad8_1, proc_rank_field);

            breaker.doBreak();
            save_or_diff(eMesh, output_files_loc+"quad_fixture_quad8_1.e");
            save_or_diff(eMesh, input_files_loc+"quad_fixture_quad8_quad8_0.e");

            // end_demo
          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Refine a quad8/serendepity mesh

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, break_quad8_to_quad8)
      {
        fixture_setup();
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        //if (p_size == 1 || p_size == 3)
        if (p_size <= 3)
          {
            percept::PerceptMesh eMesh(2u);
            eMesh.open(input_files_loc+"quad_fixture_quad8_quad8_0.e");

            Quad8_Quad8_4 break_quad8_to_quad8_4(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();

            UniformRefiner breaker(eMesh, break_quad8_to_quad8_4, proc_rank_field);
            breaker.setIgnoreSideSets(false);

            breaker.doBreak();
            save_or_diff(eMesh, output_files_loc+"quad_fixture_quad8_quad8_1.e");

          }
      }


      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Enrich a quad4 mesh to quad9 then refine it

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, break_quad4_to_quad9_to_quad9_0)
      {
        fixture_setup();
        EXCEPTWATCH;

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        bool doGenSideSets = false;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        //if (p_size == 1 || p_size == 3)
        if (p_size <= 1)
          {
            {
              const unsigned n = 1;
              //const unsigned nx = n , ny = n , nz = p_size*n ;
              const unsigned nx = n , ny = n;

              percept::QuadFixture<double> fixture( pm , nx , ny, doGenSideSets);

              // Adopt the meta/bulk data
              bool isCommitted = false;
              percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, isCommitted);

              Quad4_Quad9_1 break_quad4_to_quad9_1(eMesh);
              int scalarDimension = 0; // a scalar
              stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

              eMesh.commit();

              fixture.generate_mesh();


              UniformRefiner breaker(eMesh, break_quad4_to_quad9_1, proc_rank_field);

              breaker.doBreak();
              save_or_diff(eMesh, input_files_loc+"quad_1x1_quad9_quad9_0.e");
            }

            {
              percept::PerceptMesh em1(2);
              em1.open(input_files_loc+"quad_1x1_quad9_quad9_0.e");
              Quad9_Quad9_4 break_q9_q9(em1);
              int scalarDimension = 0; // a scalar
              stk::mesh::FieldBase* proc_rank_field = em1.addField("proc_rank", em1.element_rank(), scalarDimension);

              em1.commit();


              UniformRefiner breaker(em1, break_q9_q9, proc_rank_field);
              breaker.setIgnoreSideSets(!doGenSideSets);

              breaker.doBreak();
              save_or_diff(em1, output_files_loc+"quad_1x1_quad9_quad9_1.e");

            }
            // end_demo
          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Enrich a quad4 mesh to quad9 then refine it

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, break_quad4_to_quad9_to_quad9)
      {
        fixture_setup();
        EXCEPTWATCH;

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        bool doGenSideSets = true;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        //if (p_size == 1 || p_size == 3)
        if (p_size <= 3)
          {
            {
              //FIXME const unsigned n = 12;
              const unsigned n = 2;
              //const unsigned nx = n , ny = n , nz = p_size*n ;
              const unsigned nx = n , ny = n;

              percept::QuadFixture<double> fixture( pm , nx , ny, doGenSideSets);

              // Adopt the meta/bulk data
              bool isCommitted = false;
              percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, isCommitted);

              Quad4_Quad9_1 break_quad4_to_quad9_1(eMesh);
              int scalarDimension = 0; // a scalar
              stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

              eMesh.commit();

              fixture.generate_mesh();


              UniformRefiner breaker(eMesh, break_quad4_to_quad9_1, proc_rank_field);
              std::cout << "break_quad4_to_quad9_1.fixSurfaceAndEdgeSetNamesMap().size()= "
                        << break_quad4_to_quad9_1.fixSurfaceAndEdgeSetNamesMap().size() << std::endl;

              breaker.doBreak();
              save_or_diff(eMesh, input_files_loc+"quad_fixture_quad9_quad9_0.e");
              //eMesh.printInfo("quad_fixture_quad9_quad9_0.e", 2);
            }

            {
              percept::PerceptMesh em1(2);
              em1.open(input_files_loc+"quad_fixture_quad9_quad9_0.e");
              Quad9_Quad9_4 break_q9_q9(em1);
              int scalarDimension = 0; // a scalar
              stk::mesh::FieldBase* proc_rank_field = em1.addField("proc_rank", em1.element_rank(), scalarDimension);

              em1.commit();


              UniformRefiner breaker(em1, break_q9_q9, proc_rank_field);
              breaker.setIgnoreSideSets(!doGenSideSets);

              breaker.doBreak();
              save_or_diff(em1, output_files_loc+"quad_fixture_quad9_quad9_1.e");

            }
            // end_demo
          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Refine a triangle mesh

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, break_tri_to_tri_sierra_0)
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

            percept::QuadFixture<double, shards::Triangle<3> > fixture( pm , nx , ny, true);

            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, false);

            eMesh.commit();

            fixture.generate_mesh();

            eMesh.printInfo("tri mesh");

            save_or_diff(eMesh, output_files_loc+"quad_fixture_tri3.e");
            // end_demo
          }

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Refine a triangle mesh

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, break_tri_to_tri_sierra_1)
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

            bool isCommitted = false;
            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, isCommitted);

            Tri3_Tri3_4 break_tri_to_tri_4(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            //stk::mesh::FieldBase* proc_rank_field_edge =
            eMesh.addField("proc_rank_edge", eMesh.edge_rank(), scalarDimension);

            //             std::cout << "proc_rank_field rank= " << proc_rank_field->rank() << std::endl;
            //             std::cout << "proc_rank_field_edge rank= " << proc_rank_field_edge->rank() << std::endl;

            //fixture.meta_data.commit();
            eMesh.commit();

            fixture.generate_mesh();

            //eMesh.printInfo("tri mesh", 5);
            eMesh.printInfo("tri mesh");
            save_or_diff(eMesh, output_files_loc+"quad_fixture_tri3_0.e");

            UniformRefiner breaker(eMesh, break_tri_to_tri_4, proc_rank_field);
            breaker.doBreak();

            //eMesh.printInfo("tri mesh refined", 5);
            eMesh.printInfo("tri mesh refined");
            save_or_diff(eMesh, output_files_loc+"quad_fixture_tri3_1.e");

            if (0)
              {
                percept::PerceptMesh e1(2);
                e1.openReadOnly(input_files_loc+"quad_fixture_tri3_1.e");
                e1.printInfo("after read", 3);
              }
            // end_demo
          }

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Enrich a triangle mesh

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, break_tri3_to_tri6_sierra)
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

            bool isCommitted = false;
            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, isCommitted);

            Tri3_Tri6_1 break_tri3_to_tri6(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            //stk::mesh::FieldBase* proc_rank_field_edge =
            eMesh.addField("proc_rank_edge", eMesh.edge_rank(), scalarDimension);

            //             std::cout << "proc_rank_field rank= " << proc_rank_field->rank() << std::endl;
            //             std::cout << "proc_rank_field_edge rank= " << proc_rank_field_edge->rank() << std::endl;

            //fixture.meta_data.commit();
            eMesh.commit();

            fixture.generate_mesh();

            //eMesh.printInfo("tri mesh", 5);
            eMesh.printInfo("tri mesh tri6");
            save_or_diff(eMesh, output_files_loc+"quad_fixture_tri3_tri6_0.e");

            UniformRefiner breaker(eMesh, break_tri3_to_tri6, proc_rank_field);
            breaker.doBreak();

            //eMesh.printInfo("tri mesh refined", 5);
            eMesh.printInfo("tri mesh enriched");
            save_or_diff(eMesh, output_files_loc+"quad_fixture_tri3_tri6_1.e");
            save_or_diff(eMesh, input_files_loc+"quad_fixture_tri6_tri6_0.e");

            if (0)
              {
                percept::PerceptMesh e1(2);
                e1.openReadOnly(input_files_loc+"quad_fixture_tri3_1.e");
                e1.printInfo("after read", 3);
              }
            // end_demo
          }

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Enrich a triangle mesh then refine it

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, break_tri3_to_tri6_to_tri6_sierra)
      {
        fixture_setup();
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size <= 3)
          {
            // start_demo_break_tri3_to_tri6_to_tri6_sierra

            percept::PerceptMesh eMesh(2u);
            eMesh.open(input_files_loc+"quad_fixture_tri6_tri6_0.e");

            Tri6_Tri6_4 break_tri6_to_tri6(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            //stk::mesh::FieldBase* proc_rank_field_edge =
            eMesh.addField("proc_rank_edge", eMesh.edge_rank(), scalarDimension);
            eMesh.commit();

            eMesh.printInfo("tri mesh tri6");
            save_or_diff(eMesh, output_files_loc+"quad_fixture_tri6_tri6_0.e");

            UniformRefiner breaker(eMesh, break_tri6_to_tri6, proc_rank_field);
            breaker.doBreak();

            eMesh.printInfo("tri mesh refined");
            save_or_diff(eMesh, output_files_loc+"quad_fixture_tri6_tri6_1.e");
            // end_demo
          }

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Refine a linear tet mesh

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, break_tet4_tet4_0)
      {
        fixture_setup();
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            // start_demo_uniformRefiner_break_tet4_tet4_0
            percept::PerceptMesh eMesh(3u);
            eMesh.openReadOnly(input_files_loc+"tet_from_hex_fixture.e");
            save_or_diff(eMesh, input_files_loc+"tet_from_hex_fixture_0.e");
            // end_demo

          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Refine a linear tet mesh

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, break_tet4_tet4_1)
      {
        fixture_setup();
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            // start_demo_uniformRefiner_break_tet4_tet4_1
            percept::PerceptMesh eMesh(3u);
            eMesh.open(input_files_loc+"tet_from_hex_fixture_0.e");

            Tet4_Tet4_8 break_tet_tet(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();
            eMesh.printInfo("tet mesh");

            UniformRefiner breaker(eMesh, break_tet_tet, proc_rank_field);
            //breaker.setRemoveOldElements(false);
            //breaker.setIgnoreSideSets(true);
            breaker.doBreak();
            save_or_diff(eMesh, output_files_loc+"tet4_refined_1.e");

            breaker.doBreak();
            save_or_diff(eMesh, output_files_loc+"tet4_refined_2.e");
            // end_demo

          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Enrich a linear tet mesh

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, break_tet4_tet10_1)
      {
        fixture_setup();
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            // start_demo_uniformRefiner_break_tet4_tet10_1
            percept::PerceptMesh eMesh(3u);
            eMesh.open(input_files_loc+"tet_from_hex_fixture_0.e");

            Tet4_Tet10_1 break_tet_tet(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();
            eMesh.printInfo("tet mesh");

            UniformRefiner breaker(eMesh, break_tet_tet, proc_rank_field);
            //breaker.setRemoveOldElements(false);
            //breaker.setIgnoreSideSets(true);
            breaker.doBreak();
            save_or_diff(eMesh, output_files_loc+"tet10_1.e");
            // end_demo


          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Enrich a linear tet mesh then refine it

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, break_tet4_tet10_tet10_1)
      {
        fixture_setup();
        // FIXME
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            // start_demo_uniformRefiner_break_tet4_tet10_tet10_1
            percept::PerceptMesh eMesh(3u);
            eMesh.open(input_files_loc+"tet_from_hex_fixture_0.e");

            Tet4_Tet10_1 break_tet_tet(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();
            eMesh.printInfo("tet mesh");

            UniformRefiner breaker(eMesh, break_tet_tet, proc_rank_field);
            //breaker.setRemoveOldElements(false);
            //breaker.setIgnoreSideSets(true);
            breaker.doBreak();
            save_or_diff(eMesh, input_files_loc+"tet10_1.e");
            eMesh.printInfo("tet10_1");
            // end_demo

          }

        if (p_size == 1 || p_size == 3)
          {
            // start_demo_uniformRefiner_break_tet4_tet10_tet10_2
            percept::PerceptMesh eMesh(3u);
            eMesh.open(input_files_loc+"tet10_1.e");

            Tet10_Tet10_8 break_tet_tet(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();
            //eMesh.printInfo("tet mesh");

            UniformRefiner breaker(eMesh, break_tet_tet, proc_rank_field);
            //breaker.setRemoveOldElements(false);
            //breaker.setIgnoreSideSets(true);
            //breaker.doBreak();

            unsigned numRefines = 1;
            for (unsigned iBreak = 0; iBreak < numRefines; iBreak++)
              {
                breaker.doBreak();
              }
            
            save_or_diff(eMesh, output_files_loc+"tet10_tet10_"+toString(numRefines)+".e");
            // end_demo


          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Refine a linear hex mesh

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, hex8_hex8_8_1)
      {
        fixture_setup();
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        // start_demo_uniformRefiner_hex8_hex8_8_1

        percept::PerceptMesh eMesh(3u);

        unsigned p_size = eMesh.getParallelSize();

        std::string gmesh_spec = std::string("1x1x")+toString(p_size)+std::string("|bbox:0,0,0,1,1,"+toString(p_size) );
        eMesh.newMesh(percept::PerceptMesh::GMeshSpec(gmesh_spec));

        Hex8_Hex8_8 break_hex_to_hex(eMesh);

        int scalarDimension = 0; // a scalar
        //         int vectorDimension = 3;

        stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
        //         eMesh.addField("velocity", mesh::Node, vectorDimension);
        //         eMesh.addField("element_volume", eMesh.element_rank(), scalarDimension);

        eMesh.commit();
        eMesh.printInfo();
        save_or_diff(eMesh, std::string(output_files_loc+"")+std::string("hex_hex_cube1x1x")+toString(p_size)+std::string("-orig.e"));

        UniformRefiner breaker(eMesh, break_hex_to_hex, proc_rank_field);
        breaker.setRemoveOldElements(true);

        breaker.doBreak();

        save_or_diff(eMesh, std::string(output_files_loc+"")+std::string("hex_hex_cube1x1x")+toString(p_size)+std::string(".e"));

        // end_demo

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Refine a linear hex mesh

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, hex8_hex8_8_2)
      {
        fixture_setup();
        EXCEPTWATCH;

        // start_demo_uniformRefiner_hex8_hex8_8_2
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            percept::PerceptMesh eMesh(3u);

            eMesh.open(input_files_loc+"hex_fixture.e");

            Hex8_Hex8_8 break_hex_to_hex(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();
            eMesh.printInfo();
            save_or_diff(eMesh, output_files_loc+"hex8_0.e");

            UniformRefiner breaker(eMesh, break_hex_to_hex, proc_rank_field);
            breaker.doBreak();

            save_or_diff(eMesh, output_files_loc+"hex8_1.e");

            breaker.doBreak();
            save_or_diff(eMesh, output_files_loc+"hex8_2.e");

            // end_demo
          }

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Enrich a linear hex mesh

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, hex8_hex27_1_1)
      {
        fixture_setup();
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        // start_demo_uniformRefiner_hex8_hex27_1_1

        percept::PerceptMesh eMesh(3u);

        unsigned p_size = eMesh.getParallelSize();

        std::string gmesh_spec = std::string("1x1x")+toString(p_size)+std::string("|bbox:0,0,0,1,1,"+toString(p_size) );
        eMesh.newMesh(percept::PerceptMesh::GMeshSpec(gmesh_spec));

        Hex8_Hex27_1 break_hex_to_hex(eMesh);

        int scalarDimension = 0; // a scalar
        stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

        eMesh.commit();
        eMesh.printInfo();
        save_or_diff(eMesh, std::string(output_files_loc+"")+std::string("hex8_hex27_cube1x1x")+toString(p_size)+std::string("-orig.e"));

        UniformRefiner breaker(eMesh, break_hex_to_hex, proc_rank_field);
        breaker.setRemoveOldElements(true);

        breaker.doBreak();

        save_or_diff(eMesh, std::string(output_files_loc+"")+std::string("hex8_hex27_cube1x1x")+toString(p_size)+std::string(".e"));

        // end_demo

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Enrich a linear hex mesh

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, hex8_hex27_1_2)
      {
        fixture_setup();
        EXCEPTWATCH;

        // start_demo_uniformRefiner_hex8_hex27_1_2
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            percept::PerceptMesh eMesh(3u);

            eMesh.open(input_files_loc+"hex_fixture.e");

            Hex8_Hex27_1 break_hex_to_hex(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();
            eMesh.printInfo();
            save_or_diff(eMesh, output_files_loc+"hex27_0.e");

            UniformRefiner breaker(eMesh, break_hex_to_hex, proc_rank_field);
            breaker.doBreak();

            save_or_diff(eMesh, output_files_loc+"hex27_1.e");


            // end_demo
          }

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Enrich a linear hex mesh to serendepity hex20 elements

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, hex8_hex20_1_1)
      {
        fixture_setup();
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        // start_demo_uniformRefiner_hex8_hex20_1_1

        percept::PerceptMesh eMesh(3u);

        unsigned p_size = eMesh.getParallelSize();

        std::string gmesh_spec = std::string("1x1x")+toString(p_size)+std::string("|bbox:0,0,0,1,1,"+toString(p_size) );
        eMesh.newMesh(percept::PerceptMesh::GMeshSpec(gmesh_spec));

        Hex8_Hex20_1 break_hex_to_hex(eMesh);

        int scalarDimension = 0; // a scalar
        stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

        eMesh.commit();
        eMesh.printInfo();
        save_or_diff(eMesh, std::string(output_files_loc+"")+std::string("hex8_hex20_cube1x1x")+toString(p_size)+std::string("-orig.e"));

        UniformRefiner breaker(eMesh, break_hex_to_hex, proc_rank_field);
        breaker.setRemoveOldElements(true);

        breaker.doBreak();

        save_or_diff(eMesh, std::string(output_files_loc+"")+std::string("hex8_hex20_cube1x1x")+toString(p_size)+std::string(".e"));
        save_or_diff(eMesh, std::string(input_files_loc+"")+std::string("hex20_hex20_cube1x1x")+toString(p_size)+std::string("_0.e"));


        // end_demo

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Enrich a linear hex mesh to serendepity hex20 elements

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, hex8_hex20_1_2)
      {
        fixture_setup();
        EXCEPTWATCH;

        // start_demo_uniformRefiner_hex8_hex20_1_2
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            percept::PerceptMesh eMesh(3u);

            eMesh.open(input_files_loc+"hex_fixture.e");

            Hex8_Hex20_1 break_hex_to_hex(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();
            eMesh.printInfo();
            save_or_diff(eMesh, output_files_loc+"hex20_0.e");

            UniformRefiner breaker(eMesh, break_hex_to_hex, proc_rank_field);
            breaker.doBreak();

            save_or_diff(eMesh, output_files_loc+"hex20_1.e");
            save_or_diff(eMesh, input_files_loc+"hex20_hex20_0.e");


            // end_demo
          }

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Refine a serendepity hex20 mesh

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, hex20_hex20_1)
      {
        fixture_setup();
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        // start_demo_uniformRefiner_hex20_hex20_1

        percept::PerceptMesh eMesh(3u);

        unsigned p_size = eMesh.getParallelSize();

        if (p_size <= 3)
          {
            eMesh.open(std::string(input_files_loc+"")+std::string("hex20_hex20_cube1x1x")+toString(p_size)+std::string("_0.e"));

            Hex20_Hex20_8 break_hex_to_hex(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();
            save_or_diff(eMesh, std::string(output_files_loc+"")+std::string("hex20_hex20_cube1x1x")+toString(p_size)+std::string("_0.e"));

            UniformRefiner breaker(eMesh, break_hex_to_hex, proc_rank_field);
            breaker.setRemoveOldElements(true);

            breaker.doBreak();

            save_or_diff(eMesh, std::string(output_files_loc+"")+std::string("hex20_hex20_cube1x1x")+toString(p_size)+std::string("_1.e"));
          }

        // end_demo

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Refine a serendepity hex20 mesh

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, hex20_hex20_1_2)
      {
        fixture_setup();
        EXCEPTWATCH;

        // start_demo_uniformRefiner_hex20_hex20_1_2
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            percept::PerceptMesh eMesh(3u);

            eMesh.open(input_files_loc+"hex20_hex20_0.e");

            Hex20_Hex20_8 break_hex_to_hex(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();
            save_or_diff(eMesh, output_files_loc+"hex20_hex20_0.e");

            UniformRefiner breaker(eMesh, break_hex_to_hex, proc_rank_field);
            //breaker.setIgnoreSideSets(true);
            breaker.doBreak();

            save_or_diff(eMesh, output_files_loc+"hex20_hex20_1.e");


            // end_demo
          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Refine a quadratic hex27 mesh

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, hex27_hex27_0)
      {
        fixture_setup();
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        // start_demo_uniformRefiner_hex27_hex27_0

        int scalarDimension = 0; // a scalar
        {
          percept::PerceptMesh eMesh(3u);

          unsigned p_size = eMesh.getParallelSize();

          std::string gmesh_spec = std::string("1x1x")+toString(p_size)+std::string("|bbox:0,0,0,1,1,"+toString(p_size) );
          eMesh.newMesh(percept::PerceptMesh::GMeshSpec(gmesh_spec));

          Hex8_Hex27_1 break_hex_to_hex(eMesh);

          stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

          eMesh.commit();
          eMesh.printInfo();
          save_or_diff(eMesh, std::string(output_files_loc+"")+std::string("hex27_hex27_cube1x1x")+toString(p_size)+std::string("-orig.e"));

          UniformRefiner breaker(eMesh, break_hex_to_hex, proc_rank_field);
          breaker.setRemoveOldElements(true);

          breaker.doBreak();

          save_or_diff(eMesh, std::string(input_files_loc+"")+std::string("hex27_hex27_cube1x1x")+toString(p_size)+std::string("_0.e"));
        }


        {
          percept::PerceptMesh eMesh(3u);

          unsigned p_size = eMesh.getParallelSize();
          eMesh.open(std::string(input_files_loc+"")+std::string("hex27_hex27_cube1x1x")+toString(p_size)+std::string("_0.e"));

          Hex27_Hex27_8 break_hex_to_hex(eMesh);

          //stk::mesh::FieldBase* proc_rank_field = eMesh.getField("proc_rank");
          stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

          eMesh.commit();
          //eMesh.printInfo();

          UniformRefiner breaker(eMesh, break_hex_to_hex, proc_rank_field);
          // FIXME
          breaker.setIgnoreSideSets(true);
          breaker.setRemoveOldElements(true);

          breaker.doBreak();

          save_or_diff(eMesh, std::string(output_files_loc+"")+std::string("hex27_hex27_cube1x1x")+toString(p_size)+std::string("_1.e"));
        }

        // end_demo

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Refine a quadratic hex27 mesh

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, hex8_hex27_hex27_1)
      {
        fixture_setup();
        EXCEPTWATCH;

        // start_demo_uniformRefiner_hex8_hex27_hex27_1
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            {
              percept::PerceptMesh eMesh(3u);

              eMesh.open(input_files_loc+"hex_fixture.e");

              Hex8_Hex27_1 break_hex_to_hex(eMesh);

              int scalarDimension = 0; // a scalar
              stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

              eMesh.commit();
              eMesh.printInfo();
              save_or_diff(eMesh, output_files_loc+"hex8_hex27_0.e");

              UniformRefiner breaker(eMesh, break_hex_to_hex, proc_rank_field);
              breaker.doBreak();

              save_or_diff(eMesh, input_files_loc+"hex8_hex27_1.e");
            }

            {
              percept::PerceptMesh eMesh(3u);

              eMesh.open(input_files_loc+"hex8_hex27_1.e");

              Hex27_Hex27_8 break_hex_to_hex(eMesh);

              //stk::mesh::FieldBase* proc_rank_field = eMesh.getField("proc_rank");
              int scalarDimension = 0;
              stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

              eMesh.commit();

              UniformRefiner breaker(eMesh, break_hex_to_hex, proc_rank_field);
              //FIXME breaker.setIgnoreSideSets(false);
              breaker.setRemoveOldElements(true);

              breaker.doBreak();

              save_or_diff(eMesh, output_files_loc+"hex8_hex27_hex27_1.e");

            }

            // end_demo
          }

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Refine a linear wedge mesh

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, wedge6_2)
      {
        fixture_setup();
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        // start_demo_unit1_uniformRefiner_wedge6_2

        percept::PerceptMesh eMesh(3u);

        unsigned p_size = eMesh.getParallelSize();
        if (p_size == 1)
          {

            //         void createMesh(stk::ParallelMachine parallel_machine,
            //                         unsigned n_nodes_x, unsigned n_nodes_y, unsigned n_nodes_z,
            //                         double xmin, double xmax,
            //                         double ymin, double ymax,
            //                         double zmin, double zmax,
            //                         std::string output_filename
            //                         )
            percept::WedgeFixture wedgeFixture;

            wedgeFixture.createMesh(MPI_COMM_WORLD,
                                    4, 3, 2,
                                    0, 1,
                                    0, 1,
                                    0, 1,
                                    std::string(input_files_loc+"swept-wedge_0.e") );

            eMesh.open(input_files_loc+"swept-wedge_0.e");

            Wedge6_Wedge6_8 break_wedge(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();

            UniformRefiner breaker(eMesh, break_wedge, proc_rank_field);
            breaker.doBreak();

            save_or_diff(eMesh, output_files_loc+"swept-wedge_1.e");

          }
        // end_demo
        
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Enrich a linear wedge mesh to serendepity Wedge15

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, wedge6_enrich_1)
      {
        fixture_setup();
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        // start_demo_unit1_uniformRefiner_wedge6_enrich_1

        percept::PerceptMesh eMesh(3u);

        unsigned p_size = eMesh.getParallelSize();
        if (p_size == 1)
          {
            //         void createMesh(stk::ParallelMachine parallel_machine,
            //                         unsigned n_nodes_x, unsigned n_nodes_y, unsigned n_nodes_z,
            //                         double xmin, double xmax,
            //                         double ymin, double ymax,
            //                         double zmin, double zmax,
            //                         std::string output_filename
            //                         )
            percept::WedgeFixture wedgeFixture;

            wedgeFixture.createMesh(MPI_COMM_WORLD,
                                    4, 3, 2,
                                    0, 1,
                                    0, 1,
                                    0, 1,
                                    std::string(input_files_loc+"swept-wedge_enrich_0.e") );

            eMesh.open(input_files_loc+"swept-wedge_enrich_0.e");

            Wedge6_Wedge15_1 break_wedge(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();

            UniformRefiner breaker(eMesh, break_wedge, proc_rank_field);
            breaker.doBreak();

            save_or_diff(eMesh, output_files_loc+"swept-wedge_enrich_1.e");
            save_or_diff(eMesh, input_files_loc+"swept-wedge_enrich_refine_0.e");

          }
        // end_demo
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Enrich a linear wedge mesh to serendepity Wedge15 then refine it

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, wedge6_enrich_refine)
      {
        fixture_setup();
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        // start_demo_unit1_uniformRefiner_wedge6_enrich_refine

        const unsigned p_size = stk::parallel_machine_size( MPI_COMM_WORLD );

        if (p_size == 1)
          {
            PerceptMesh eMesh(3u);
            percept::WedgeFixture wedgeFixture;

            wedgeFixture.createMesh(MPI_COMM_WORLD,
                                    4, 2, 2,
                                    0, 1,
                                    0, 1,
                                    0, 1,
                                    std::string(input_files_loc+"tmp-swept-wedge_enrich_0.e") );

            eMesh.open(input_files_loc+"tmp-swept-wedge_enrich_0.e");

            Wedge6_Wedge15_1 break_wedge(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
            eMesh.commit();
            UniformRefiner breaker(eMesh, break_wedge, proc_rank_field);
            breaker.doBreak();
            save_or_diff(eMesh, input_files_loc+"swept-wedge_2_enrich_refine_0.e");
          }

        percept::PerceptMesh eMesh(3u);

        if (p_size == 1)
          {
            eMesh.open(input_files_loc+"swept-wedge_2_enrich_refine_0.e");

            Wedge15_Wedge15_8 break_wedge(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();

            UniformRefiner breaker(eMesh, break_wedge, proc_rank_field);
            breaker.setIgnoreSideSets(true);
            breaker.doBreak();

            save_or_diff(eMesh, output_files_loc+"swept-wedge_2_enrich_refine_1.e");

          }
        // end_demo

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Generate a heterogeneous mesh (tet, hex, wedge elements) then refine it

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, heterogeneous_mesh)
      {
        fixture_setup();
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        // start_demo_heterogeneous_mesh

        //const unsigned p_rank = stk::parallel_machine_rank( MPI_COMM_WORLD);
        const unsigned p_size = stk::parallel_machine_size( MPI_COMM_WORLD);

        if (p_size <= 1)
          {
            // create the mesh
            {
              stk::percept::HeterogeneousFixture mesh(MPI_COMM_WORLD, false);
              stk::io::put_io_part_attribute(  mesh.m_block_hex );
              stk::io::put_io_part_attribute(  mesh.m_block_wedge );
              stk::io::put_io_part_attribute(  mesh.m_block_tet );

#if HET_FIX_INCLUDE_EXTRA_ELEM_TYPES
              stk::io::put_io_part_attribute(  mesh.m_block_pyramid );
              stk::io::put_io_part_attribute(  mesh.m_block_quad_shell );
              stk::io::put_io_part_attribute(  mesh.m_block_tri_shell );
#endif

              mesh.m_metaData.commit();
              mesh.populate();

              bool isCommitted = true;
              percept::PerceptMesh em1(&mesh.m_metaData, &mesh.m_bulkData, isCommitted);

              //em1.printInfo("heterogeneous", 4);

              save_or_diff(em1, input_files_loc+"heterogeneous_0.e");
              em1.close();
            }

            std::string input_mesh = input_files_loc+"heterogeneous_0.e";
            if (p_size > 1)
              {
                RunEnvironment::doLoadBalance(pm, input_mesh);
              }

            // refine the mesh
            if (1)
              {
                percept::PerceptMesh eMesh1(3);

                eMesh1.open(input_files_loc+"heterogeneous_0.e");

                URP_Heterogeneous_3D break_pattern(eMesh1);
                int scalarDimension = 0; // a scalar
                stk::mesh::FieldBase* proc_rank_field = eMesh1.addField("proc_rank", eMesh1.element_rank(), scalarDimension);
                eMesh1.commit();

                UniformRefiner breaker(eMesh1, break_pattern, proc_rank_field);

                //breaker.setRemoveOldElements(false);
                breaker.setIgnoreSideSets(true);
                breaker.doBreak();

                RefinementInfoByType::printTable(std::cout, breaker.getRefinementInfoByType(), 0, true );
                RefinementInfoByType::printTable(std::cout, breaker.getRefinementInfoByType(), 0, false );

                save_or_diff(eMesh1, output_files_loc+"heterogeneous_1.e");
                eMesh1.close();
              }
          }
        // end_demo
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Enrich a heterogeneous mesh

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, heterogeneous_mesh_enrich)
      {
        fixture_setup();
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        // start_demo_heterogeneous_mesh_enrich

        //const unsigned p_rank = stk::parallel_machine_rank( MPI_COMM_WORLD);
        const unsigned p_size = stk::parallel_machine_size( MPI_COMM_WORLD);

        if (p_size <= 1)
          {
            // create the mesh
            {
              stk::percept::HeterogeneousFixture mesh(MPI_COMM_WORLD, false);
              stk::io::put_io_part_attribute(  mesh.m_block_hex );
              stk::io::put_io_part_attribute(  mesh.m_block_wedge );
              stk::io::put_io_part_attribute(  mesh.m_block_tet );

#if HET_FIX_INCLUDE_EXTRA_ELEM_TYPES
              stk::io::put_io_part_attribute(  mesh.m_block_pyramid );
              stk::io::put_io_part_attribute(  mesh.m_block_quad_shell );
              stk::io::put_io_part_attribute(  mesh.m_block_tri_shell );
#endif

              mesh.m_metaData.commit();
              mesh.populate();

              bool isCommitted = true;
              percept::PerceptMesh em1(&mesh.m_metaData, &mesh.m_bulkData, isCommitted);
              
              //em1.printInfo("hetero_enrich", 4);


              save_or_diff(em1, input_files_loc+"heterogeneous_enrich_0.e");
              em1.close();
            }

            std::string input_mesh = input_files_loc+"heterogeneous_enrich_0.e";
            if (p_size > 1)
              {
                RunEnvironment::doLoadBalance(pm, input_mesh);
              }

            // enrich the mesh
            if (1)
              {
                percept::PerceptMesh eMesh1(3);

                eMesh1.open(input_files_loc+"heterogeneous_enrich_0.e");
                //eMesh1.printInfo("hetero_enrich_2", 4);

                URP_Heterogeneous_Enrich_3D break_pattern(eMesh1);
                //int scalarDimension = 0; // a scalar
                stk::mesh::FieldBase* proc_rank_field = 0;      //eMesh1.addField("proc_rank", eMesh.element_rank(), scalarDimension);
                eMesh1.commit();
                //eMesh1.printInfo("hetero_enrich_2", 4);

                UniformRefiner breaker(eMesh1, break_pattern, proc_rank_field);

                //breaker.setRemoveOldElements(false);
                breaker.setIgnoreSideSets(true);
                breaker.doBreak();

                save_or_diff(eMesh1, output_files_loc+"heterogeneous_enrich_1.e");
                save_or_diff(eMesh1, input_files_loc+"heterogeneous_quadratic_refine_0.e");
                eMesh1.close();
              }
          }
        // end_demo
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Refine the enriched heterogeneous mesh

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, heterogeneous_quadratic_refine)
      {
        fixture_setup();
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        // start_demo_heterogeneous_quadratic_refine

        //const unsigned p_rank = stk::parallel_machine_rank( MPI_COMM_WORLD);
        const unsigned p_size = stk::parallel_machine_size(pm);

        if (p_size <= 1)
          {
            // refine the mesh
            if (1)
              {
                percept::PerceptMesh eMesh1(3);

                eMesh1.open(input_files_loc+"heterogeneous_quadratic_refine_0.e");

                URP_Heterogeneous_QuadraticRefine_3D break_pattern(eMesh1);
                int scalarDimension = 0; // a scalar
                stk::mesh::FieldBase* proc_rank_field = eMesh1.addField("proc_rank", eMesh1.element_rank(), scalarDimension);
                eMesh1.commit();

                UniformRefiner breaker(eMesh1, break_pattern, proc_rank_field);

                //breaker.setRemoveOldElements(false);
                breaker.setIgnoreSideSets(true);
                breaker.doBreak();

                save_or_diff(eMesh1, output_files_loc+"heterogeneous_quadratic_refine_1.e");
                eMesh1.close();
              }
          }
        // end_demo
      }

      //here

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// Enrich a wedge6 mesh to wedge18

      STKUNIT_UNIT_TEST(unit1_uniformRefiner, wedge6_wedge18_enrich)
      {
        fixture_setup();
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        // start_demo_unit1_uniformRefiner_wedge6_wedge18_enrich

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( MPI_COMM_WORLD);
        const unsigned p_size = stk::parallel_machine_size(pm);


        //unsigned p_size = eMesh.getParallelSize();
        if (p_size == 1)
          {
            //         void createMesh(stk::ParallelMachine parallel_machine,
            //                         unsigned n_nodes_x, unsigned n_nodes_y, unsigned n_nodes_z,
            //                         double xmin, double xmax,
            //                         double ymin, double ymax,
            //                         double zmin, double zmax,
            //                         std::string output_filename
            //                         )
            percept::WedgeFixture wedgeFixture;

            mesh::BulkData *bulk = 
              wedgeFixture.createMesh(MPI_COMM_WORLD,
                                      4, 3, 2,
                                      0, 1,
                                      0, 1,
                                      0, 1,
                                      std::string(""));
            //std::string("swept-wedge6_18_enrich_0.e") );

            percept::PerceptMesh eMesh(wedgeFixture.getMetaData(), bulk, false);
            //percept::PerceptMesh eMesh;
            //eMesh.open("swept-wedge6_18_enrich_0.e");

            Wedge6_Wedge18_1 break_wedge(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();
            //eMesh.printInfo();

            wedgeFixture.createBulkAfterMetaCommit(MPI_COMM_WORLD);

            UniformRefiner breaker(eMesh, break_wedge, proc_rank_field);
            breaker.doBreak();

          }
        // end_demo
      }




    } // namespace unit_tests
  } // namespace adapt
} // namespace stk

