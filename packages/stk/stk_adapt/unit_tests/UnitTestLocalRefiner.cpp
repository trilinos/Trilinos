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

            eMesh.printInfo("local tri mesh",2);
            save_or_diff(eMesh, output_files_loc+"local_tri_0.e");

            TestLocalRefinerTri breaker(eMesh, break_tri_to_tri_2, proc_rank_field);
            breaker.setRemoveOldElements(false);
            breaker.doBreak();

            eMesh.printInfo("local tri mesh refined", 2);
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

            eMesh.printInfo("local tri mesh",2);
            save_or_diff(eMesh, output_files_loc+"local_tri_1_0.e");

            bool diagonals=true;
            TestLocalRefinerTri2 breaker(eMesh, break_tri_to_tri_2, proc_rank_field, diagonals);
            //breaker.setRemoveOldElements(false);
            breaker.doBreak();

            eMesh.printInfo("local tri mesh refined", 2);
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

            eMesh.printInfo("local tri mesh",2);
            save_or_diff(eMesh, output_files_loc+"local_tri_2_0.e");

            bool diagonals=false;
            TestLocalRefinerTri2 breaker(eMesh, break_tri_to_tri_2, proc_rank_field, diagonals);
            //breaker.setRemoveOldElements(false);
            breaker.doBreak();

            eMesh.printInfo("local tri mesh refined", 2);
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

            eMesh.printInfo("local tri mesh",2);
            save_or_diff(eMesh, output_files_loc+"local_tri_N_1_0.e");

            TestLocalRefinerTri_N_1 breaker(eMesh, break_tri_to_tri_N, proc_rank_field);
            breaker.setRemoveOldElements(false);
            breaker.doBreak();

            eMesh.printInfo("local tri mesh refined", 2);
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



    } // namespace unit_tests
  } // namespace adapt
} // namespace stk


