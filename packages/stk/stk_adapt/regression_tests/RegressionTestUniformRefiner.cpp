/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
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
#include <fstream>
#include <string>
#include <typeinfo>

#include <math.h>
#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/util/PrintTable.hpp>

#include <Teuchos_ScalarTraits.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <boost/lexical_cast.hpp>
#include <stk_io/IossBridge.hpp>

#include <stk_percept/Percept.hpp>
#include <stk_percept/Util.hpp>
#include <stk_percept/ExceptionWatch.hpp>

#include <stk_percept/function/StringFunction.hpp>
#include <stk_percept/function/FieldFunction.hpp>
#include <stk_percept/function/ConstantFunction.hpp>
#include <stk_percept/PerceptMesh.hpp>

#include <stk_adapt/UniformRefinerPattern.hpp>
#include <stk_adapt/UniformRefiner.hpp>
#include <unit_tests/TestLocalRefiner.hpp>
#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>
#include <stk_percept/RunEnvironment.hpp>

#include <stk_percept/fixtures/Fixture.hpp>
#include <stk_percept/fixtures/BeamFixture.hpp>
#include <stk_percept/fixtures/HeterogeneousFixture.hpp>
#include <stk_percept/fixtures/QuadFixture.hpp>
#include <stk_percept/fixtures/WedgeFixture.hpp>

#include <use_cases/UseCase_3.hpp>

// this is for testing the local-refine refactoring 
#define UNIFORM_REFINER UniformRefiner
//#define UNIFORM_REFINER TestLocalRefiner

namespace stk
{
  namespace adapt
  {
    namespace regression_tests
    {

#define EXTRA_PRINT 0

      static void output_draw(std::string filename, std::string toFile)
      {
        std::ofstream file(filename.c_str());
        file << toFile;
      }


      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================
      //========= AREA for tests in progress of being debugged
      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================


#if 0
      STKUNIT_UNIT_TEST(regr_uniformRefiner, break_quad4_to_quad9_to_quad9_shell_1)
      {
        EXCEPTWATCH;

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );

        // this case can't be load balanced (I presume there are too few elements)

        if (p_size <= 1)
          {
            // start_demo_break_quad4_to_quad9_to_quad9_shell
            std::string input_mesh = "./input_files/shell-tests/freshell_quad4.g";
            if (p_size > 1)
              {
                RunEnvironment::doLoadBalance(pm, input_mesh);
              }

            percept::PerceptMesh eMesh(3u);
            eMesh.open(input_mesh);

            ShellQuad4_ShellQuad9_1 break_quad4_to_quad9_1(eMesh);

            int scalarDimension = 0; // a scalar
            //         int vectorDimension = 3;

            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();

            eMesh.printInfo("quad mesh");
            eMesh.saveAs("./output_files/freshell_quad4_quad9_0.g");

            UNIFORM_REFINER breaker(eMesh, break_quad4_to_quad9_1, proc_rank_field);
            //breaker.setIgnoreSideSets(true);
            breaker.doBreak();

            //eMesh.printInfo("quad mesh refined", 5);
            eMesh.printInfo("quad shell mesh enriched");
            eMesh.saveAs("./output_files/freshell_quad4_quad9_1.g");
            eMesh.saveAs("./input_files/freshell_quad9_quad9_0.g");

          }
      }
#endif


      STKUNIT_UNIT_TEST(regr_uniformRefiner, beam_enrich)
      {
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
              em1.saveAs("./input_files/beam_enrich_0.e");

            }

            // enrich
            {
              stk::percept::PerceptMesh eMesh(3u);
              eMesh.open("./input_files/beam_enrich_0.e");
              //URP_Heterogeneous_3D break_pattern(eMesh);
              Beam2_Beam3_1 break_pattern(eMesh);
              int scalarDimension = 0; // a scalar
              stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
              eMesh.commit();

              eMesh.saveAs("./output_files/beam_enrich_0.e");

              eMesh.printInfo("beam", 2);

              UNIFORM_REFINER breaker(eMesh, break_pattern, proc_rank_field);
              //breaker.setRemoveOldElements(false);
              breaker.setIgnoreSideSets(true);
              breaker.doBreak();

              eMesh.saveAs("./output_files/beam_enrich_1.e");

            }
          }
      }

      STKUNIT_UNIT_TEST(regr_uniformRefiner, beam_refine)
      {
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
              em1.saveAs("./input_files/beam_0.e");

            }

            // refine
            {
              stk::percept::PerceptMesh eMesh(3u);
              eMesh.open("./input_files/beam_0.e");
              //URP_Heterogeneous_3D break_pattern(eMesh);
              Beam2_Beam2_2 break_pattern(eMesh);
              int scalarDimension = 0; // a scalar
              stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
              eMesh.commit();

              eMesh.saveAs("./output_files/beam_0.e");

              eMesh.printInfo("beam", 2);

              UNIFORM_REFINER breaker(eMesh, break_pattern, proc_rank_field);
              //breaker.setRemoveOldElements(false);
              breaker.setIgnoreSideSets(true);
              breaker.doBreak();

              eMesh.saveAs("./output_files/beam_1.e");

            }
          }
      }


      //======================================================================================================================
      //======================================================================================================================
      //===================== Table generation
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, generate_tables)
      {
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
      //===================== Shell elements testing
      //======================================================================================================================
      //======================================================================================================================

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================


      STKUNIT_UNIT_TEST(regr_uniformRefiner, break_quad4_to_quad8_to_quad8_shell)
      {
        EXCEPTWATCH;

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );

        // this case can't be load balanced (I presume there are too few elements)

        if (p_size <= 1)
          {
            // start_demo_break_quad4_to_quad8_to_quad8_shell
            std::string input_mesh = "./input_files/shell-tests/freshell_quad4.g";
            if (p_size > 1)
              {
                RunEnvironment::doLoadBalance(pm, input_mesh);
              }

            percept::PerceptMesh eMesh(3u);
            eMesh.open(input_mesh);

            ShellQuad4_ShellQuad8_1 break_quad4_to_quad8_1(eMesh);

            int scalarDimension = 0; // a scalar
            //         int vectorDimension = 3;

            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();

            eMesh.printInfo("quad mesh");
            eMesh.saveAs("./output_files/freshell_quad4_quad8_0.g");

            UNIFORM_REFINER breaker(eMesh, break_quad4_to_quad8_1, proc_rank_field);
            //breaker.setIgnoreSideSets(true);
            breaker.doBreak();

            //eMesh.printInfo("quad mesh refined", 5);
            eMesh.printInfo("quad shell mesh enriched");
            eMesh.saveAs("./output_files/freshell_quad4_quad8_1.g");
            eMesh.saveAs("./input_files/freshell_quad8_quad8_0.g");

          }

        if (1 && p_size <= 1)
          {

            percept::PerceptMesh eMesh(3u);
            eMesh.open("./input_files/freshell_quad8_quad8_0.g");

            ShellQuad8_ShellQuad8_4 break_quad8_to_quad_8(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();

            //eMesh.printInfo("quad mesh");
            //eMesh.saveAs("./output_files/freshell_quad4_quad8_0.g");

            UNIFORM_REFINER breaker(eMesh, break_quad8_to_quad_8, proc_rank_field);
            //breaker.setIgnoreSideSets(true);
            breaker.doBreak();

            //eMesh.printInfo("quad mesh refined", 5);
            eMesh.printInfo("quad shell mesh enriched and refined");
            eMesh.saveAs("./output_files/freshell_quad8_quad8_1.g");
            // end_demo

          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================


#if 0
      STKUNIT_UNIT_TEST(regr_uniformRefiner, break_quad4_to_quad9_to_quad9_shell)
      {
        EXCEPTWATCH;

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );

        // this case can't be load balanced (I presume there are too few elements)

        if (p_size <= 1)
          {
            // start_demo_break_quad4_to_quad9_to_quad9_shell
            std::string input_mesh = "./input_files/shell-tests/freshell_quad4.g";
            if (p_size > 1)
              {
                RunEnvironment::doLoadBalance(pm, input_mesh);
              }

            percept::PerceptMesh eMesh(3u);
            eMesh.open(input_mesh);

            ShellQuad4_ShellQuad9_1 break_quad4_to_quad9_1(eMesh);

            int scalarDimension = 0; // a scalar
            //         int vectorDimension = 3;

            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();

            eMesh.printInfo("quad mesh");
            eMesh.saveAs("./output_files/freshell_quad4_quad9_0.g");

            UNIFORM_REFINER breaker(eMesh, break_quad4_to_quad9_1, proc_rank_field);
            //breaker.setIgnoreSideSets(true);
            breaker.doBreak();

            //eMesh.printInfo("quad mesh refined", 5);
            eMesh.printInfo("quad shell mesh enriched");
            eMesh.saveAs("./output_files/freshell_quad4_quad9_1.g");
            eMesh.saveAs("./input_files/freshell_quad9_quad9_0.g");

          }

#if 0
        if (1 && p_size <= 1)
          {

            percept::PerceptMesh eMesh(3u);
            eMesh.open("./input_files/freshell_quad9_quad9_0.g");

            ShellQuad9_ShellQuad9_4 break_quad9_to_quad_9(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();

            //eMesh.printInfo("quad mesh");
            //eMesh.saveAs("./output_files/freshell_quad4_quad9_0.g");

            UNIFORM_REFINER breaker(eMesh, break_quad9_to_quad_9, proc_rank_field);
            //breaker.setIgnoreSideSets(true);
            breaker.doBreak();

            //eMesh.printInfo("quad mesh refined", 5);
            eMesh.printInfo("quad shell mesh enriched and refined");
            eMesh.saveAs("./output_files/freshell_quad9_quad9_1.g");
            // end_demo

          }
#endif
      }
#endif

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================


      STKUNIT_UNIT_TEST(regr_uniformRefiner, break_quad_to_quad_shell)
      {
        EXCEPTWATCH;

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );

        // this case can't be load balanced (I presume there are too few elements)

        if (p_size <= 1)
          {
            // start_demo_break_quad_to_quad_shell
            std::string input_mesh = "./input_files/shell-tests/freshell_quad4.g";
            if (p_size > 1)
              {
                RunEnvironment::doLoadBalance(pm, input_mesh);
              }

            percept::PerceptMesh eMesh(3u);
            eMesh.open(input_mesh);

            ShellQuad4_ShellQuad4_4 break_quad_to_quad_4(eMesh);

            int scalarDimension = 0; // a scalar
            //         int vectorDimension = 3;

            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();

            eMesh.printInfo("quad mesh");
            eMesh.saveAs("./output_files/freshell_quad4_0.g");

            UNIFORM_REFINER breaker(eMesh, break_quad_to_quad_4, proc_rank_field);
            //breaker.setIgnoreSideSets(true);
            breaker.doBreak();

            //eMesh.printInfo("quad mesh refined", 5);
            eMesh.printInfo("quad shell mesh refined");
            eMesh.saveAs("./output_files/freshell_quad4_1.g");
            // end_demo

          }
      }

      STKUNIT_UNIT_TEST(regr_uniformRefiner, break_tri_to_tri_shell)
      {
        EXCEPTWATCH;

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        // this case can't be load balanced (I presume there are too few elements)

        if (p_size <= 1)
          {
            // start_demo_break_tri_to_tri_shell

            std::string input_mesh = "./input_files/shell-tests/freshell_tri3.g";
            if (p_size > 1)
              {
                RunEnvironment::doLoadBalance(pm, input_mesh);
              }

            percept::PerceptMesh eMesh(3u);
            eMesh.open(input_mesh);


            ShellTri3_ShellTri3_4 break_tri_to_tri_4(eMesh);

            int scalarDimension = 0; // a scalar
            //         int vectorDimension = 3;

            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();

            eMesh.printInfo("tri mesh");
            eMesh.saveAs("./output_files/freshell_tri3_0.g");

            UNIFORM_REFINER breaker(eMesh, break_tri_to_tri_4, proc_rank_field);
            //breaker.setIgnoreSideSets(true);
            breaker.doBreak();

            //eMesh.printInfo("tri mesh refined", 5);
            eMesh.printInfo("tri shell mesh refined");
            eMesh.saveAs("./output_files/freshell_tri3_1.g");
            // end_demo

          }

      }



      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, draw1)
      {
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

      STKUNIT_UNIT_TEST(regr_uniformRefiner, draw)
      {
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
        // enrich
        //     typedef  UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<9>, 1, SierraPort >            Quad4_Quad9_1;
        //     typedef  UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<8>, 1, SierraPort >            Quad4_Quad8_1;
        //     typedef  UniformRefinerPattern<shards::Triangle<3>,      shards::Triangle<6>,      1, SierraPort >            Tri3_Tri6_1;
        //     typedef  UniformRefinerPattern<shards::Tetrahedron<4>,   shards::Tetrahedron<10>,  1, SierraPort >            Tet4_Tet10_1;
        //     typedef  UniformRefinerPattern<shards::Hexahedron<8>,    shards::Hexahedron<27>,   1, SierraPort >            Hex8_Hex27_1;
        //     typedef  UniformRefinerPattern<shards::Hexahedron<8>,    shards::Hexahedron<20>,   1, SierraPort >            Hex8_Hex20_1;

        //     // convert
        //     typedef  UniformRefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>,      6 >                        Quad4_Tri3_6;
        //     typedef  UniformRefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>,      4, Specialization >        Quad4_Tri3_4;
        //     typedef  UniformRefinerPattern<shards::Hexahedron<8>,    shards::Tetrahedron<4>,  24 >                        Hex8_Tet4_24;

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, break_quad_to_tri_6)
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            // start_demo_uniformRefiner_break_quad_to_tri
            percept::PerceptMesh eMesh(2u);
            eMesh.open("./input_files/break_test/quad/square/square_quad4.e");

            typedef  UniformRefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>, 6 > Quad4_Tri3_6;

            UniformRefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>, 6 > break_quad_to_tri_6(eMesh);

            int scalarDimension = 0; // a scalar
            //         int vectorDimension = 3;

            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();

            eMesh.printInfo("quad mesh");

            //UniformRefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>, 6 > break_quad_to_tri_6;
            UNIFORM_REFINER breaker(eMesh, break_quad_to_tri_6, proc_rank_field);
            breaker.setRemoveOldElements(false);
            breaker.setIgnoreSideSets(true);
            breaker.doBreak();

            //eMesh.saveAs("./output_files/break_test/quad/square/square_quad4_out.e");
            eMesh.saveAs("./output_files/square_quad4_out.e");
            // end_demo
          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, break_quad_to_tri_4)
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {

            // start_demo_uniformRefiner_break_quad_to_tri
            percept::PerceptMesh eMesh(2u);
            eMesh.open("./input_files/break_test/quad/square/square_quad4.e");

            UniformRefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>, 4, Specialization > break_quad_to_tri_4(eMesh);

            int scalarDimension = 0; // a scalar

            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();

            eMesh.printInfo("quad mesh");

            //UniformRefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>, 6 > break_quad_to_tri_6;
            UNIFORM_REFINER breaker(eMesh, break_quad_to_tri_4, proc_rank_field);
            breaker.setRemoveOldElements(false);
            breaker.setIgnoreSideSets(true);
            breaker.doBreak();

            //eMesh.saveAs("./input_files/break_test/quad/square/square_quad4_out.e");
            eMesh.saveAs("./output_files/square_quad4_tri3_4_out.e");
            // end_demo
          }
      }


      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, break_quad_to_quad)
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            // start_demo_uniformRefiner_break_quad_to_quad
            percept::PerceptMesh eMesh(2u);
            eMesh.open("./input_files/break_test/quad/square/square_quad4.e");

            UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4 > break_quad_to_quad_4(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
            eMesh.commit();

            eMesh.printInfo("quad mesh");

            UNIFORM_REFINER breaker(eMesh, break_quad_to_quad_4, proc_rank_field);
            //breaker.setRemoveOldElements(false);
            breaker.setIgnoreSideSets(true);
            breaker.doBreak();

            eMesh.saveAs("./output_files/square_quad4_ref_out.e");
            // end_demo
          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// uses the Sierra-ported tables from framework/{element,mesh_modification}

      STKUNIT_UNIT_TEST(regr_uniformRefiner, break_quad_to_quad_sierra)
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            // start_demo_uniformRefiner_break_quad_to_quad
            percept::PerceptMesh eMesh(2u);
            eMesh.open("./input_files/break_test/quad/square/square_quad4.e");

            UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4, SierraPort > break_quad_to_quad_4(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
            eMesh.commit();

            eMesh.printInfo("quad mesh");

            UNIFORM_REFINER breaker(eMesh, break_quad_to_quad_4, proc_rank_field);
            //breaker.setRemoveOldElements(false);
            breaker.doBreak();

            eMesh.saveAs("./output_files/square_quad4_ref_sierra_out.e");
            // end_demo
          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      /// uses the Sierra-ported tables from framework/{element,mesh_modification}

      STKUNIT_UNIT_TEST(regr_uniformRefiner, break_quad_to_quad_sierra_sidesets)
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );

        //FIXME
        if (0)
          {
            percept::PerceptMesh eMesh(2u);
            eMesh.open("./input_files/break_test/quad/sidesets/quad_sidesets.e");
            eMesh.commit();
            eMesh.printInfo("quad mesh");
          }

        if (p_size == 1 || p_size == 2)
          {
            // start_demo_uniformRefiner_break_quad_to_quad
            percept::PerceptMesh eMesh(2u);
            eMesh.open("./input_files/break_test/quad/sidesets/quad_sidesets.e");

            UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4, SierraPort > break_quad_to_quad_4(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
            eMesh.commit();

            eMesh.printInfo("quad mesh");

            UNIFORM_REFINER breaker(eMesh, break_quad_to_quad_4, proc_rank_field);
            //breaker.setRemoveOldElements(false);
            breaker.doBreak();

            eMesh.printInfo("after refinement break_quad_to_quad_sierra_sidesets");

            eMesh.saveAs("./output_files/quad_sidesets_sierra_out.e");
          }
        // end_demo
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================
      // FIXME - move and/or copy to unit tests
      STKUNIT_UNIT_TEST(regr_uniformRefiner, hex8_tet4_24_1)
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        // start_demo_uniformRefiner_hex8_tet4_24_1

        percept::PerceptMesh eMesh(3u);

        unsigned p_size = eMesh.getParallelSize();
        //unsigned p_rank = eMesh.getRank();
        Util::setRank(eMesh.getRank());

        std::string gmesh_spec = std::string("1x1x")+toString(p_size)+std::string("|bbox:0,0,0,1,1,"+toString(p_size) );
        eMesh.newMesh(percept::GMeshSpec(gmesh_spec));

        UniformRefinerPattern<shards::Hexahedron<8>, shards::Tetrahedron<4>, 24 > break_hex_to_tet(eMesh);

        int scalarDimension = 0; // a scalar
        //         int vectorDimension = 3;

        stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
        //         eMesh.addField("velocity", mesh::Node, vectorDimension);
        //         eMesh.addField("element_volume", eMesh.element_rank(), scalarDimension);

        eMesh.commit();
        eMesh.printInfo();
        eMesh.saveAs(std::string("./output_files/")+std::string("hex_tet_24_cube1x1x")+toString(p_size)+std::string("-orig.e"));

        UNIFORM_REFINER breaker(eMesh, break_hex_to_tet, proc_rank_field);
        breaker.setRemoveOldElements(true);

        breaker.doBreak();

        //eMesh.saveAs("./input_files/break_test/quad/square/square_quad4_out.e");
        eMesh.saveAs(std::string("./output_files/")+std::string("hex_tet_24_cube1x1x")+toString(p_size)+std::string(".e"));

        // end_demo

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, hex8_tet4_6_12_1)
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        // start_demo_uniformRefiner_hex8_tet4_6_12_1

        percept::PerceptMesh eMesh(3u);

        unsigned p_size = eMesh.getParallelSize();

        std::string gmesh_spec = std::string("1x1x")+toString(p_size)+std::string("|bbox:0,0,0,1,1,"+toString(p_size) );
        eMesh.newMesh(percept::GMeshSpec(gmesh_spec));

        Hex8_Tet4_6_12 break_hex_to_tet(eMesh);

        int scalarDimension = 0; // a scalar
        //         int vectorDimension = 3;

        stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
        //         eMesh.addField("velocity", mesh::Node, vectorDimension);
        //         eMesh.addField("element_volume", eMesh.element_rank(), scalarDimension);

        eMesh.commit();
        eMesh.printInfo();
        eMesh.saveAs(std::string("./output_files/")+std::string("hex_tet_6_12_cube1x1x")+toString(p_size)+std::string("-orig.e"));

        UNIFORM_REFINER breaker(eMesh, break_hex_to_tet, proc_rank_field);
        breaker.setRemoveOldElements(true);
        //breaker.setIgnoreSideSets(true);

        breaker.doBreak();

        eMesh.saveAs(std::string("./output_files/")+std::string("hex_tet_6_12_cube1x1x")+toString(p_size)+std::string(".e"));

        // end_demo
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, hex8_tet4_6_12_2)
      {

        EXCEPTWATCH;

        // start_demo_uniformRefiner_hex8_tet4_6_12_2
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            percept::PerceptMesh eMesh(3u);

            eMesh.open("./input_files/break_test/hex/cylinder/cylinder_hex8.e");

            Hex8_Tet4_6_12 break_hex_to_tet(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();
            //eMesh.printInfo("test",2);
            eMesh.printInfo();
            eMesh.saveAs("./output_files/cylinder_hex8_tet4_6_12_0.e");

            UNIFORM_REFINER breaker(eMesh, break_hex_to_tet, proc_rank_field);
            //breaker.setIgnoreSideSets(true);
            breaker.doBreak();

            eMesh.saveAs("./output_files/cylinder_hex8_tet4_6_12_1.e");

            // end_demo
          }

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, quad4_quad4_4_test_1)
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
            eMesh.printInfo("quad fixture");
            eMesh.saveAs("./output_files/quad_fixture.e");
          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, break_quad_to_quad_sierra_1)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        bool doGenSideSets = true;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        //if (p_size == 1 || p_size == 3)
        if (p_size == 1 || p_size == 3)
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

            eMesh.saveAs("./output_files/quad_fixture_0.e");
            eMesh.close();

            for (int iBreak = 0; iBreak < 2; iBreak++)
              {
                std::cout << "\n\n\n ================ tmp Refine Pass = " << iBreak << std::endl;

                percept::PerceptMesh eMesh1(2);
                std::string fileName = std::string("./input_files/quad_fixture_")+toString(iBreak)+std::string(".e");
                eMesh1.open(fileName);
                UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4, SierraPort > break_quad_to_quad_4(eMesh1);
                int scalarDimension = 0; // a scalar
                stk::mesh::FieldBase* proc_rank_field = eMesh1.addField("proc_rank", eMesh.element_rank(), scalarDimension);
                eMesh1.commit();

                //                 if (iBreak != 0)
                //                   proc_rank_field = eMesh1.getField("proc_rank");

                UNIFORM_REFINER breaker(eMesh1, break_quad_to_quad_4, proc_rank_field);

                //breaker.setRemoveOldElements(false);
                breaker.doBreak();
                std::string fileName1 = std::string("./output_files/quad_fixture_")+toString(iBreak+1)+std::string(".e");
                //eMesh1.saveAs(fileName+"_ref.e");
                //eMesh1.printInfo("quad_fixture_1.e");

                eMesh1.saveAs(fileName1);
                eMesh1.close();

                if (0 && iBreak==0)
                  {
                    percept::PerceptMesh e1(2);
                    std::cout << "\n\n\n ================ tmp eMesh1.openReadOnly(quad_fixture_1.e) \n\n\n " << std::endl;
                    e1.openReadOnly("./input_files/quad_fixture_1.e");
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

      STKUNIT_UNIT_TEST(regr_uniformRefiner, break_quad_to_quad_sierra_2)
      {
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

            eMesh.saveAs("./output_files/quad_fixture_mbreak_0.e");
            eMesh.close();


            percept::PerceptMesh eMesh1(2);
            std::string fileName = std::string("./input_files/quad_fixture_mbreak_0.e");
            eMesh1.open(fileName);
            UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4, SierraPort > break_quad_to_quad_4(eMesh1);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh1.addField("proc_rank", eMesh.element_rank(), scalarDimension);
            eMesh1.commit();

            UNIFORM_REFINER breaker(eMesh1, break_quad_to_quad_4, proc_rank_field);

            for (int iBreak = 0; iBreak < 2; iBreak++)
              {
                std::cout << "\n\n\n ================ tmp Refine Pass = " << iBreak << std::endl;

                breaker.doBreak();
                std::string fileName1 = std::string("./output_files/quad_fixture_mbreak_")+toString(iBreak+1)+std::string(".e");
                eMesh1.saveAs(fileName1);
              }

            // end_demo
          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, break_quad4_to_quad9)
      {
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

            eMesh.saveAs("./output_files/quad_fixture_quad9_0.e");

            UNIFORM_REFINER breaker(eMesh, break_quad4_to_quad9_1, proc_rank_field);

            breaker.doBreak();
            eMesh.saveAs("./output_files/quad_fixture_quad9_1.e");

            // end_demo
          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, break_quad4_to_quad8)
      {
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

            Quad4_Quad8_1 break_quad4_to_quad8_1(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();

            fixture.generate_mesh();

            eMesh.saveAs("./output_files/quad_fixture_quad8_0.e");

            UNIFORM_REFINER breaker(eMesh, break_quad4_to_quad8_1, proc_rank_field);

            breaker.doBreak();
            eMesh.saveAs("./output_files/quad_fixture_quad8_1.e");
            eMesh.saveAs("./input_files/quad_fixture_quad8_quad8_0.e");

            // end_demo
          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, break_quad8_to_quad8)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        //if (p_size == 1 || p_size == 3)
        if (p_size <= 3)
          {
            percept::PerceptMesh eMesh(2u);
            eMesh.open("./input_files/quad_fixture_quad8_quad8_0.e");

            Quad8_Quad8_4 break_quad8_to_quad8_4(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();

            UNIFORM_REFINER breaker(eMesh, break_quad8_to_quad8_4, proc_rank_field);
            breaker.setIgnoreSideSets(false);

            breaker.doBreak();
            eMesh.saveAs("./output_files/quad_fixture_quad8_quad8_1.e");

          }
      }


      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================
      STKUNIT_UNIT_TEST(regr_uniformRefiner, break_quad4_to_quad9_to_quad9_0)
      {
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

              //eMesh.saveAs("./output_files/quad_fixture_quad9_0.e");

              UNIFORM_REFINER breaker(eMesh, break_quad4_to_quad9_1, proc_rank_field);

              breaker.doBreak();
              eMesh.saveAs("./input_files/quad_1x1_quad9_quad9_0.e");
            }

            {
              percept::PerceptMesh em1(2);
              em1.open("./input_files/quad_1x1_quad9_quad9_0.e");
              Quad9_Quad9_4 break_q9_q9(em1);
              int scalarDimension = 0; // a scalar
              stk::mesh::FieldBase* proc_rank_field = em1.addField("proc_rank", em1.element_rank(), scalarDimension);

              em1.commit();

              //em1.saveAs("./output_files/quad_1x1_quad9_0.e");

              UNIFORM_REFINER breaker(em1, break_q9_q9, proc_rank_field);
              breaker.setIgnoreSideSets(!doGenSideSets);

              breaker.doBreak();
              em1.saveAs("./output_files/quad_1x1_quad9_quad9_1.e");

            }
            // end_demo
          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================
      STKUNIT_UNIT_TEST(regr_uniformRefiner, break_quad4_to_quad9_to_quad9)
      {
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

              //eMesh.saveAs("./output_files/quad_fixture_quad9_0.e");

              UNIFORM_REFINER breaker(eMesh, break_quad4_to_quad9_1, proc_rank_field);
              std::cout << "break_quad4_to_quad9_1.fixSurfaceAndEdgeSetNamesMap().size()= "
                        << break_quad4_to_quad9_1.fixSurfaceAndEdgeSetNamesMap().size() << std::endl;

              breaker.doBreak();
              eMesh.saveAs("./input_files/quad_fixture_quad9_quad9_0.e");
              //eMesh.printInfo("quad_fixture_quad9_quad9_0.e", 2);
            }

            {
              percept::PerceptMesh em1(2);
              em1.open("./input_files/quad_fixture_quad9_quad9_0.e");
              Quad9_Quad9_4 break_q9_q9(em1);
              int scalarDimension = 0; // a scalar
              stk::mesh::FieldBase* proc_rank_field = em1.addField("proc_rank", em1.element_rank(), scalarDimension);

              em1.commit();

              //em1.saveAs("./output_files/quad_fixture_quad9_0.e");

              UNIFORM_REFINER breaker(em1, break_q9_q9, proc_rank_field);
              breaker.setIgnoreSideSets(!doGenSideSets);

              breaker.doBreak();
              em1.saveAs("./output_files/quad_fixture_quad9_quad9_1.e");

            }
            // end_demo
          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, break_tri_to_tri_sierra_0)
      {
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

            eMesh.saveAs("./output_files/quad_fixture_tri3.e");
            // end_demo
          }

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, break_tri_to_tri_sierra_1)
      {
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
            eMesh.saveAs("./output_files/quad_fixture_tri3_0.e");

            UNIFORM_REFINER breaker(eMesh, break_tri_to_tri_4, proc_rank_field);
            breaker.doBreak();

            //eMesh.printInfo("tri mesh refined", 5);
            eMesh.printInfo("tri mesh refined");
            eMesh.saveAs("./output_files/quad_fixture_tri3_1.e");

            if (0)
              {
                percept::PerceptMesh e1(2);
                e1.openReadOnly("./input_files/quad_fixture_tri3_1.e");
                e1.printInfo("after read", 3);
              }
            // end_demo
          }

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, break_tri3_to_tri6_sierra)
      {
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
            eMesh.saveAs("./output_files/quad_fixture_tri3_tri6_0.e");

            UNIFORM_REFINER breaker(eMesh, break_tri3_to_tri6, proc_rank_field);
            breaker.doBreak();

            //eMesh.printInfo("tri mesh refined", 5);
            eMesh.printInfo("tri mesh enriched");
            eMesh.saveAs("./output_files/quad_fixture_tri3_tri6_1.e");
            eMesh.saveAs("./input_files/quad_fixture_tri6_tri6_0.e");

            if (0)
              {
                percept::PerceptMesh e1(2);
                e1.openReadOnly("./input_files/quad_fixture_tri3_1.e");
                e1.printInfo("after read", 3);
              }
            // end_demo
          }

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, break_tri3_to_tri6_to_tri6_sierra)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size <= 3)
          {
            // start_demo_break_tri3_to_tri6_to_tri6_sierra

            percept::PerceptMesh eMesh(2u);
            eMesh.open("./input_files/quad_fixture_tri6_tri6_0.e");

            Tri6_Tri6_4 break_tri6_to_tri6(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            //stk::mesh::FieldBase* proc_rank_field_edge =
            eMesh.addField("proc_rank_edge", eMesh.edge_rank(), scalarDimension);
            eMesh.commit();

            eMesh.printInfo("tri mesh tri6");
            eMesh.saveAs("./output_files/quad_fixture_tri6_tri6_0.e");

            UNIFORM_REFINER breaker(eMesh, break_tri6_to_tri6, proc_rank_field);
            breaker.doBreak();

            eMesh.printInfo("tri mesh refined");
            eMesh.saveAs("./output_files/quad_fixture_tri6_tri6_1.e");
            // end_demo
          }

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, break_tet4_tet4_0)
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            // start_demo_uniformRefiner_break_tet4_tet4
            percept::PerceptMesh eMesh(3u);
            eMesh.openReadOnly("./input_files/break_test/tet/cylinder-from-hex/cylinder_tet4.e");

            eMesh.saveAs("./output_files/cylinder_tet4_0.e");
            // end_demo

          }
      }

#if 1
      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, break_tet4_tet4_1)
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            // start_demo_uniformRefiner_break_tet4_tet4_1
            percept::PerceptMesh eMesh(3u);
            eMesh.open("./input_files/break_test/tet/cylinder-from-hex/cylinder_tet4_0.e");

            Tet4_Tet4_8 break_tet_tet(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();
            eMesh.printInfo("tet mesh");

            UNIFORM_REFINER breaker(eMesh, break_tet_tet, proc_rank_field);
            //breaker.setRemoveOldElements(false);
            //breaker.setIgnoreSideSets(true);
            breaker.doBreak();
            eMesh.saveAs("./output_files/cylinder_tet4_1.e");

            breaker.doBreak();
            eMesh.saveAs("./output_files/cylinder_tet4_2.e");
            // end_demo

          }
      }
#endif
      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, break_tet4_tet10_1)
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            // start_demo_uniformRefiner_break_tet4_tet4_1
            percept::PerceptMesh eMesh(3u);
            eMesh.open("./input_files/break_test/tet/cylinder-from-hex/cylinder_tet4_0.e");

            Tet4_Tet10_1 break_tet_tet(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();
            eMesh.printInfo("tet mesh");

            UNIFORM_REFINER breaker(eMesh, break_tet_tet, proc_rank_field);
            //breaker.setRemoveOldElements(false);
            //breaker.setIgnoreSideSets(true);
            breaker.doBreak();
            eMesh.saveAs("./output_files/cylinder_tet10_1.e");
            // end_demo


          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, break_tet4_tet10_tet10_1)
      {
        // FIXME
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            // start_demo_uniformRefiner_break_tet4_tet4_1
            percept::PerceptMesh eMesh(3u);
            eMesh.open("./input_files/break_test/tet/cylinder-from-hex/cylinder_tet4_0.e");

            Tet4_Tet10_1 break_tet_tet(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();
            eMesh.printInfo("tet mesh");

            UNIFORM_REFINER breaker(eMesh, break_tet_tet, proc_rank_field);
            //breaker.setRemoveOldElements(false);
            //breaker.setIgnoreSideSets(true);
            breaker.doBreak();
            eMesh.saveAs("./input_files/cylinder_tet10_1.e");
            eMesh.printInfo("cylinder_tet10_1");
            // end_demo

          }

        if (p_size == 1 || p_size == 3)
          {
            // start_demo_uniformRefiner_break_tet4_tet4_1
            percept::PerceptMesh eMesh(3u);
            eMesh.open("./input_files/cylinder_tet10_1.e");

            Tet10_Tet10_8 break_tet_tet(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();
            //eMesh.printInfo("tet mesh");

            UNIFORM_REFINER breaker(eMesh, break_tet_tet, proc_rank_field);
            //breaker.setRemoveOldElements(false);
            //breaker.setIgnoreSideSets(true);
            //breaker.doBreak();

            unsigned numRefines = 1;
            for (unsigned iBreak = 0; iBreak < numRefines; iBreak++)
              {
                breaker.doBreak();
              }
            
            eMesh.saveAs("./output_files/cylinder_tet10_tet10_"+toString(numRefines)+".e");
            // end_demo


          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, hex8_hex8_8_1)
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        // start_demo_uniformRefiner_hex8_hex8_8_1

        percept::PerceptMesh eMesh(3u);

        unsigned p_size = eMesh.getParallelSize();

        std::string gmesh_spec = std::string("1x1x")+toString(p_size)+std::string("|bbox:0,0,0,1,1,"+toString(p_size) );
        eMesh.newMesh(percept::GMeshSpec(gmesh_spec));

        Hex8_Hex8_8 break_hex_to_hex(eMesh);

        int scalarDimension = 0; // a scalar
        //         int vectorDimension = 3;

        stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
        //         eMesh.addField("velocity", mesh::Node, vectorDimension);
        //         eMesh.addField("element_volume", eMesh.element_rank(), scalarDimension);

        eMesh.commit();
        eMesh.printInfo();
        eMesh.saveAs(std::string("./output_files/")+std::string("hex_hex_cube1x1x")+toString(p_size)+std::string("-orig.e"));

        UNIFORM_REFINER breaker(eMesh, break_hex_to_hex, proc_rank_field);
        breaker.setRemoveOldElements(true);

        breaker.doBreak();

        eMesh.saveAs(std::string("./output_files/")+std::string("hex_hex_cube1x1x")+toString(p_size)+std::string(".e"));

        // end_demo

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, hex8_hex8_8_2)
      {
        EXCEPTWATCH;

        // start_demo_uniformRefiner_hex8_hex8_8_2
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            percept::PerceptMesh eMesh(3u);

            eMesh.open("./input_files/break_test/hex/cylinder/cylinder_hex8.e");

            Hex8_Hex8_8 break_hex_to_hex(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();
            eMesh.printInfo();
            eMesh.saveAs("./output_files/cylinder_hex8_0.e");

            UNIFORM_REFINER breaker(eMesh, break_hex_to_hex, proc_rank_field);
            breaker.doBreak();

            eMesh.saveAs("./output_files/cylinder_hex8_1.e");

            breaker.doBreak();
            eMesh.saveAs("./output_files/cylinder_hex8_2.e");

            // end_demo
          }

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, hex8_hex27_1_1)
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        // start_demo_uniformRefiner_hex8_hex27_1_1

        percept::PerceptMesh eMesh(3u);

        unsigned p_size = eMesh.getParallelSize();

        std::string gmesh_spec = std::string("1x1x")+toString(p_size)+std::string("|bbox:0,0,0,1,1,"+toString(p_size) );
        eMesh.newMesh(percept::GMeshSpec(gmesh_spec));

        Hex8_Hex27_1 break_hex_to_hex(eMesh);

        int scalarDimension = 0; // a scalar
        stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

        eMesh.commit();
        eMesh.printInfo();
        eMesh.saveAs(std::string("./output_files/")+std::string("hex8_hex27_cube1x1x")+toString(p_size)+std::string("-orig.e"));

        UNIFORM_REFINER breaker(eMesh, break_hex_to_hex, proc_rank_field);
        breaker.setRemoveOldElements(true);

        breaker.doBreak();

        eMesh.saveAs(std::string("./output_files/")+std::string("hex8_hex27_cube1x1x")+toString(p_size)+std::string(".e"));

        // end_demo

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, hex8_hex27_1_2)
      {
        EXCEPTWATCH;

        // start_demo_uniformRefiner_hex8_hex27_1_2
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            percept::PerceptMesh eMesh(3u);

            eMesh.open("./input_files/break_test/hex/cylinder/cylinder_hex8.e");

            Hex8_Hex27_1 break_hex_to_hex(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();
            eMesh.printInfo();
            eMesh.saveAs("./output_files/cylinder_hex27_0.e");

            UNIFORM_REFINER breaker(eMesh, break_hex_to_hex, proc_rank_field);
            breaker.doBreak();

            eMesh.saveAs("./output_files/cylinder_hex27_1.e");


            // end_demo
          }

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, hex8_hex20_1_1)
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        // start_demo_uniformRefiner_hex8_hex20_1_1

        percept::PerceptMesh eMesh(3u);

        unsigned p_size = eMesh.getParallelSize();

        std::string gmesh_spec = std::string("1x1x")+toString(p_size)+std::string("|bbox:0,0,0,1,1,"+toString(p_size) );
        eMesh.newMesh(percept::GMeshSpec(gmesh_spec));

        Hex8_Hex20_1 break_hex_to_hex(eMesh);

        int scalarDimension = 0; // a scalar
        stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

        eMesh.commit();
        eMesh.printInfo();
        eMesh.saveAs(std::string("./output_files/")+std::string("hex8_hex20_cube1x1x")+toString(p_size)+std::string("-orig.e"));

        UNIFORM_REFINER breaker(eMesh, break_hex_to_hex, proc_rank_field);
        breaker.setRemoveOldElements(true);

        breaker.doBreak();

        eMesh.saveAs(std::string("./output_files/")+std::string("hex8_hex20_cube1x1x")+toString(p_size)+std::string(".e"));
        eMesh.saveAs(std::string("./input_files/")+std::string("hex20_hex20_cube1x1x")+toString(p_size)+std::string("_0.e"));


        // end_demo

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, hex8_hex20_1_2)
      {
        EXCEPTWATCH;

        // start_demo_uniformRefiner_hex8_hex20_1_2
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            percept::PerceptMesh eMesh(3u);

            eMesh.open("./input_files/break_test/hex/cylinder/cylinder_hex8.e");

            Hex8_Hex20_1 break_hex_to_hex(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();
            eMesh.printInfo();
            eMesh.saveAs("./output_files/cylinder_hex20_0.e");

            UNIFORM_REFINER breaker(eMesh, break_hex_to_hex, proc_rank_field);
            breaker.doBreak();

            eMesh.saveAs("./output_files/cylinder_hex20_1.e");
            eMesh.saveAs("./input_files/cylinder_hex20_hex20_0.e");


            // end_demo
          }

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, hex20_hex20_1)
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        // start_demo_uniformRefiner_hex20_hex20_1

        percept::PerceptMesh eMesh(3u);

        unsigned p_size = eMesh.getParallelSize();

        if (p_size <= 3)
          {
            eMesh.open(std::string("./input_files/")+std::string("hex20_hex20_cube1x1x")+toString(p_size)+std::string("_0.e"));

            Hex20_Hex20_8 break_hex_to_hex(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();
            eMesh.saveAs(std::string("./output_files/")+std::string("hex20_hex20_cube1x1x")+toString(p_size)+std::string("_0.e"));

            UNIFORM_REFINER breaker(eMesh, break_hex_to_hex, proc_rank_field);
            breaker.setRemoveOldElements(true);

            breaker.doBreak();

            eMesh.saveAs(std::string("./output_files/")+std::string("hex20_hex20_cube1x1x")+toString(p_size)+std::string("_1.e"));
          }

        // end_demo

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, hex20_hex20_1_2)
      {
        EXCEPTWATCH;

        // start_demo_uniformRefiner_hex20_hex20_1_2
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            percept::PerceptMesh eMesh(3u);

            eMesh.open("./input_files/cylinder_hex20_hex20_0.e");

            Hex20_Hex20_8 break_hex_to_hex(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();
            eMesh.saveAs("./output_files/cylinder_hex20_hex20_0.e");

            UNIFORM_REFINER breaker(eMesh, break_hex_to_hex, proc_rank_field);
            //breaker.setIgnoreSideSets(true);
            breaker.doBreak();

            eMesh.saveAs("./output_files/cylinder_hex20_hex20_1.e");


            // end_demo
          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, hex27_hex27_0)
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        // start_demo_uniformRefiner_hex27_hex27_0

        int scalarDimension = 0; // a scalar
        {
          percept::PerceptMesh eMesh(3u);

          unsigned p_size = eMesh.getParallelSize();

          std::string gmesh_spec = std::string("1x1x")+toString(p_size)+std::string("|bbox:0,0,0,1,1,"+toString(p_size) );
          eMesh.newMesh(percept::GMeshSpec(gmesh_spec));

          Hex8_Hex27_1 break_hex_to_hex(eMesh);

          stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

          eMesh.commit();
          eMesh.printInfo();
          eMesh.saveAs(std::string("./output_files/")+std::string("hex27_hex27_cube1x1x")+toString(p_size)+std::string("-orig.e"));

          UNIFORM_REFINER breaker(eMesh, break_hex_to_hex, proc_rank_field);
          breaker.setRemoveOldElements(true);

          breaker.doBreak();

          eMesh.saveAs(std::string("./input_files/")+std::string("hex27_hex27_cube1x1x")+toString(p_size)+std::string("_0.e"));
        }


        {
          percept::PerceptMesh eMesh(3u);

          unsigned p_size = eMesh.getParallelSize();
          eMesh.open(std::string("./input_files/")+std::string("hex27_hex27_cube1x1x")+toString(p_size)+std::string("_0.e"));

          Hex27_Hex27_8 break_hex_to_hex(eMesh);

          //stk::mesh::FieldBase* proc_rank_field = eMesh.getField("proc_rank");
          stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

          eMesh.commit();
          //eMesh.printInfo();
          //eMesh.saveAs(std::string("./output_files/")+std::string("hex27_hex27_cube1x1x")+toString(p_size)+std::string("-orig.e"));

          UNIFORM_REFINER breaker(eMesh, break_hex_to_hex, proc_rank_field);
          // FIXME
          breaker.setIgnoreSideSets(true);
          breaker.setRemoveOldElements(true);

          breaker.doBreak();

          eMesh.saveAs(std::string("./output_files/")+std::string("hex27_hex27_cube1x1x")+toString(p_size)+std::string("_1.e"));
        }

        // end_demo

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, hex8_hex27_hex27_1)
      {
        EXCEPTWATCH;

        // start_demo_uniformRefiner_hex8_hex27_1_2
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            {
              percept::PerceptMesh eMesh(3u);

              eMesh.open("./input_files/break_test/hex/cylinder/cylinder_hex8.e");

              Hex8_Hex27_1 break_hex_to_hex(eMesh);

              int scalarDimension = 0; // a scalar
              stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

              eMesh.commit();
              eMesh.printInfo();
              eMesh.saveAs("./output_files/cylinder_hex8_hex27_0.e");

              UNIFORM_REFINER breaker(eMesh, break_hex_to_hex, proc_rank_field);
              breaker.doBreak();

              eMesh.saveAs("./input_files/cylinder_hex8_hex27_1.e");
            }

            {
              percept::PerceptMesh eMesh(3u);

              eMesh.open("./input_files/cylinder_hex8_hex27_1.e");

              Hex27_Hex27_8 break_hex_to_hex(eMesh);

              //stk::mesh::FieldBase* proc_rank_field = eMesh.getField("proc_rank");
              int scalarDimension = 0;
              stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

              eMesh.commit();

              UNIFORM_REFINER breaker(eMesh, break_hex_to_hex, proc_rank_field);
              //FIXME breaker.setIgnoreSideSets(false);
              breaker.setRemoveOldElements(true);

              breaker.doBreak();

              eMesh.saveAs("./output_files/cylinder_hex8_hex27_hex27_1.e");

            }

            // end_demo
          }

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, wedge6_2)
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        // start_demo_regr_uniformRefiner_wedge6_1

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
                                    std::string("swept-wedge_0.e") );

            eMesh.open("swept-wedge_0.e");

            Wedge6_Wedge6_8 break_wedge(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();

            UNIFORM_REFINER breaker(eMesh, break_wedge, proc_rank_field);
            breaker.doBreak();

            eMesh.saveAs("./output_files/swept-wedge_1.e");

          }
        // end_demo
        
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, wedge6_enrich_1)
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        // start_demo_regr_uniformRefiner_wedge6_enrich_1

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
                                    std::string("swept-wedge_enrich_0.e") );

            eMesh.open("swept-wedge_enrich_0.e");

            Wedge6_Wedge15_1 break_wedge(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();

            UNIFORM_REFINER breaker(eMesh, break_wedge, proc_rank_field);
            breaker.doBreak();

            eMesh.saveAs("./output_files/swept-wedge_enrich_1.e");
            eMesh.saveAs("./input_files/swept-wedge_enrich_refine_0.e");

          }
        // end_demo
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, wedge6_enrich_refine)
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        // start_demo_regr_uniformRefiner_wedge6_enrich_refine

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
                                    std::string("tmp-swept-wedge_enrich_0.e") );

            eMesh.open("tmp-swept-wedge_enrich_0.e");

            Wedge6_Wedge15_1 break_wedge(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
            eMesh.commit();
            UNIFORM_REFINER breaker(eMesh, break_wedge, proc_rank_field);
            breaker.doBreak();
            eMesh.saveAs("./input_files/swept-wedge_enrich_refine_0.e");
          }

        percept::PerceptMesh eMesh(3u);

        if (p_size == 1)
          {
            eMesh.open("./input_files/swept-wedge_enrich_refine_0.e");

            Wedge15_Wedge15_8 break_wedge(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();

            UNIFORM_REFINER breaker(eMesh, break_wedge, proc_rank_field);
            breaker.setIgnoreSideSets(true);
            breaker.doBreak();

            eMesh.saveAs("./output_files/swept-wedge_enrich_refine_1.e");

          }
        // end_demo

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, heterogeneous_mesh)
      {
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

              em1.saveAs("./input_files/heterogeneous_0.e");
              em1.close();
            }

            std::string input_mesh = "./input_files/heterogeneous_0.e";
            if (p_size > 1)
              {
                RunEnvironment::doLoadBalance(pm, input_mesh);
              }

            // refine the mesh
            if (1)
              {
                percept::PerceptMesh eMesh1(3);

                eMesh1.open("./input_files/heterogeneous_0.e");

                URP_Heterogeneous_3D break_pattern(eMesh1);
                int scalarDimension = 0; // a scalar
                stk::mesh::FieldBase* proc_rank_field = eMesh1.addField("proc_rank", eMesh1.element_rank(), scalarDimension);
                eMesh1.commit();

                UNIFORM_REFINER breaker(eMesh1, break_pattern, proc_rank_field);

                //breaker.setRemoveOldElements(false);
                breaker.setIgnoreSideSets(true);
                breaker.doBreak();

                eMesh1.saveAs("./output_files/heterogeneous_1.e");
                eMesh1.close();
              }
          }
        // end_demo
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================


#if 0
      STKUNIT_UNIT_TEST(regr_uniformRefiner, beam_refine)
      {
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
              em1.saveAs("./input_files/beam_0.e");

            }

            // refine
            {
              stk::percept::PerceptMesh eMesh(3u);
              eMesh.open("./input_files/beam_0.e");
              //URP_Heterogeneous_3D break_pattern(eMesh);
              Beam2_Beam2_2 break_pattern(eMesh);
              int scalarDimension = 0; // a scalar
              stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
              eMesh.commit();

              eMesh.saveAs("./output_files/beam_0.e");

              eMesh.printInfo("beam", 2);

              UNIFORM_REFINER breaker(eMesh, break_pattern, proc_rank_field);
              //breaker.setRemoveOldElements(false);
              breaker.setIgnoreSideSets(true);
              breaker.doBreak();

              eMesh.saveAs("./output_files/beam_1.e");

            }
          }
      }

#endif
      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================
#if 1

      STKUNIT_UNIT_TEST(regr_uniformRefiner, biplane_refine)
      {
        EXCEPTWATCH;

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );

        // this case can't be load balanced?

        if (p_size <= 1)
          {
            // start_demo_biplane_refine
            std::string input_mesh = "./input_files/salinas/biplane.e";

            if (p_size > 1)
              {
                RunEnvironment::doLoadBalance(pm, input_mesh);
              }

            percept::PerceptMesh eMesh(3u);
            eMesh.open(input_mesh);

            URP_Heterogeneous_3D break_pattern(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
            eMesh.commit();

            eMesh.printInfo("biplane", 2);
            eMesh.saveAs("./output_files/biplane_0.e");

            UNIFORM_REFINER breaker(eMesh, break_pattern, proc_rank_field);
            //breaker.setRemoveOldElements(false);
            breaker.setIgnoreSideSets(false);
            breaker.doBreak();

            //eMesh.printInfo("biplane", 2);
            eMesh.saveAs("./output_files/biplane_1.e");

          }
      }

#endif
      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, break_tet_shell3_tet)
      {
        EXCEPTWATCH;

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        std::cout << "p_size= " << p_size << std::endl;
        // this case can't be load balanced (I presume there are too few elements)
        if (p_size <= 1)
          {
            // start_demo_break_tet_shell3_tet

            std::string input_mesh = "./input_files/shell-tests/tet_shell3_tet.g";
            if (p_size > 1)
              {
                RunEnvironment::doLoadBalance(pm, input_mesh);
              }

            percept::PerceptMesh eMesh(3u);
            eMesh.open(input_mesh);

            URP_Heterogeneous_3D break_pattern(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
            eMesh.commit();

            UNIFORM_REFINER breaker(eMesh, break_pattern, proc_rank_field);

            //breaker.setRemoveOldElements(false);
            //breaker.setIgnoreSideSets(true);
            unsigned numRefines = 2;
            for (unsigned iBreak = 0; iBreak < numRefines; iBreak++)
              {
                breaker.doBreak();
              }

            eMesh.saveAs("./output_files/tet_shell3_tet_"+toString(numRefines)+".g");
            eMesh.close();

            // end_demo

          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, break_hex_shell4_hex)
      {
        EXCEPTWATCH;

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        // this case can't be load balanced (I presume there are too few elements)
        //if (p_size <= 3)
        if (p_size == 1)
          {
            // start_demo_break_hex_shell4_hex
            std::string input_mesh = "./input_files/shell-tests/hex_shell4_hex.g";
            if (p_size > 1)
              {
                RunEnvironment::doLoadBalance(pm, input_mesh);
              }

            percept::PerceptMesh eMesh(3u);
            eMesh.open(input_mesh);

            URP_Heterogeneous_3D break_pattern(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
            eMesh.commit();

            UNIFORM_REFINER breaker(eMesh, break_pattern, proc_rank_field);

            //breaker.setRemoveOldElements(false);
            //breaker.setIgnoreSideSets(true);
            breaker.doBreak();
            eMesh.saveAs("./output_files/hex_shell4_hex_1.g");
            eMesh.close();

            // end_demo

          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, heterogeneous_mesh_enrich)
      {
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


              em1.saveAs("./input_files/heterogeneous_enrich_0.e");
              em1.close();
            }

            std::string input_mesh = "./input_files/heterogeneous_enrich_0.e";
            if (p_size > 1)
              {
                RunEnvironment::doLoadBalance(pm, input_mesh);
              }

            // enrich the mesh
            if (1)
              {
                percept::PerceptMesh eMesh1(3);

                eMesh1.open("./input_files/heterogeneous_enrich_0.e");
                //eMesh1.printInfo("hetero_enrich_2", 4);

                URP_Heterogeneous_Enrich_3D break_pattern(eMesh1);
                //int scalarDimension = 0; // a scalar
                stk::mesh::FieldBase* proc_rank_field = 0;      //eMesh1.addField("proc_rank", eMesh.element_rank(), scalarDimension);
                eMesh1.commit();
                //eMesh1.printInfo("hetero_enrich_2", 4);

                UNIFORM_REFINER breaker(eMesh1, break_pattern, proc_rank_field);

                //breaker.setRemoveOldElements(false);
                breaker.setIgnoreSideSets(true);
                breaker.doBreak();

                eMesh1.saveAs("./output_files/heterogeneous_enrich_1.e");
                eMesh1.saveAs("./input_files/heterogeneous_quadratic_refine_0.e");
                eMesh1.close();
              }
          }
        // end_demo
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, heterogeneous_quadratic_refine)
      {
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

                eMesh1.open("./input_files/heterogeneous_quadratic_refine_0.e");

                URP_Heterogeneous_QuadraticRefine_3D break_pattern(eMesh1);
                int scalarDimension = 0; // a scalar
                stk::mesh::FieldBase* proc_rank_field = eMesh1.addField("proc_rank", eMesh1.element_rank(), scalarDimension);
                eMesh1.commit();

                UNIFORM_REFINER breaker(eMesh1, break_pattern, proc_rank_field);

                //breaker.setRemoveOldElements(false);
                breaker.setIgnoreSideSets(true);
                breaker.doBreak();

                eMesh1.saveAs("./output_files/heterogeneous_quadratic_refine_1.e");
                eMesh1.close();
              }
          }
        // end_demo
      }

      //here

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      STKUNIT_UNIT_TEST(regr_uniformRefiner, wedge6_wedge18_enrich)
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        // start_demo_regr_uniformRefiner_wedge6_wedge18_enrich

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

            UNIFORM_REFINER breaker(eMesh, break_wedge, proc_rank_field);
            breaker.doBreak();

            //eMesh.saveAs("./output_files/swept-wedge6_18_enrich_0.e");
            //eMesh.saveAs("./input_files/swept-wedge6_18_enrich_refine_0.e");
          }
        // end_demo
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

#if 0

      STKUNIT_UNIT_TEST(regr_uniformRefiner, wedge6_18_enrich_refine)
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        // start_demo_regr_uniformRefiner_wedge6_18_enrich_refine

        percept::PerceptMesh eMesh(3u);

        unsigned p_size = eMesh.getParallelSize();
        if (p_size == 1)
          {
            eMesh.open("./input_files/swept-wedge_enrich_refine_0.e");

            Wedge15_Wedge15_8 break_wedge(eMesh);

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

            eMesh.commit();

            UNIFORM_REFINER breaker(eMesh, break_wedge, proc_rank_field);
            breaker.setIgnoreSideSets(true);
            breaker.doBreak();

            eMesh.saveAs("./output_files/swept-wedge_enrich_refine_1.e");
          }
        // end_demo
      }
#endif

    }//    namespace regression_tests
  }//  namespace adapt
}// namespace stk

