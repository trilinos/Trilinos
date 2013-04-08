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
#include <stk_util/diag/PrintTable.hpp>

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
#include <stk_adapt/Refiner.hpp>

#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>
#include <stk_percept/RunEnvironment.hpp>

#include <stk_percept/fixtures/Fixture.hpp>
#include <stk_percept/fixtures/BeamFixture.hpp>
#include <stk_percept/fixtures/HeterogeneousFixture.hpp>
#include <stk_percept/fixtures/QuadFixture.hpp>
#include <stk_percept/fixtures/WedgeFixture.hpp>

#include <stk_adapt/IEdgeBasedAdapterPredicate.hpp>
#include <stk_adapt/IElementBasedAdapterPredicate.hpp>
#include <stk_adapt/PredicateBasedElementAdapter.hpp>
#include <stk_adapt/PredicateBasedEdgeAdapter.hpp>

#include <stk_percept/function/ElementOp.hpp>

#include <regression_tests/RegressionTestLocalRefiner.hpp>

//#include "RegressionTestLocalRefiner.hpp"

namespace stk
{
  namespace adapt
  {
    namespace regression_tests
    {

#if !defined(__IBMCPP__)
      bool DO_TESTS=true;
#else
      // maybe it takes too long on dawn so we turn off these long-running tests...
      bool DO_TESTS=false;
#endif

      bool LARGE_TEST_ONLY=false;

#include "RegressionTestFileLoc.hpp"

//       const std::string input_files_loc="./input_files/";
//       const std::string output_files_loc="./output_files/";

      static std::string post_fix[9] = {"np0", "np1", "np2", "np3", "np4", "np5", "np6", "np7", "np8"};

#define EXTRA_PRINT 0

      STKUNIT_UNIT_TEST(regr_localRefiner, break_tet_to_tet_N_5_ElementBased)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        if (LARGE_TEST_ONLY || !DO_TESTS) return;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            // start_demo_local_refiner_break_tet_to_tet_2

            percept::PerceptMesh eMesh;
            eMesh.open(input_files_loc+"cylinder_with_5_holes.e");

            Local_Tet4_Tet4_N break_tet_to_tet_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.add_field("proc_rank", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);
            stk::mesh::FieldBase* refine_field = eMesh.add_field("refine_field", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);

            eMesh.commit();

            SetRefineField set_ref_field(eMesh);
            eMesh.elementOpLoop(set_ref_field, refine_field);

            SetUnrefineField set_unref_field(eMesh);
            //eMesh.elementOpLoop(set_ref_field, refine_field);

            eMesh.save_as( output_files_loc+"local_tet_N_5_ElementBased_0_"+post_fix[p_size]+".e");

            ElementRefinePredicate erp(0, refine_field, 0.0);

            PredicateBasedElementAdapter<ElementRefinePredicate>
              breaker(erp,
                      eMesh, break_tet_to_tet_N, proc_rank_field);

            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);
            for (int ipass=0; ipass < 3; ipass++)
              {
                eMesh.elementOpLoop(set_ref_field, refine_field);

                std::cout << "P[" << eMesh.get_rank() << "] ipass= " << ipass << std::endl;
                breaker.doBreak();
                std::cout << "P[" << eMesh.get_rank() << "] done... ipass= " << ipass << std::endl;
                eMesh.save_as(output_files_loc+"local_tet_N_5_ElementBased_1_ipass_"+toString(ipass)+"_"+post_fix[p_size]+".e");
              }

            breaker.deleteParentElements();
            eMesh.save_as(output_files_loc+"local_tet_N_5_ElementBased_1_"+post_fix[p_size]+".e");

#if 0
            for (int iunref_pass=0; iunref_pass < 4; iunref_pass++)
              {
                eMesh.elementOpLoop(set_unref_field, refine_field);
                std::cout << "P[" << eMesh.get_rank() << "] iunref_pass= " << iunref_pass << std::endl;
                ElementUnrefineCollection elements_to_unref = breaker.buildUnrefineList();
                breaker.unrefineTheseElements(elements_to_unref);
                eMesh.save_as(output_files_loc+"local_tet_N_5_ElementBased_1_unref_ipass_"+toString(iunref_pass)+"_"+post_fix[p_size]+".e");
              }

            eMesh.save_as( output_files_loc+"local_tet_N_5_ElementBased_1_unref_"+post_fix[p_size]+".e");
#endif
            // end_demo
          }

      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      STKUNIT_UNIT_TEST(regr_localRefiner, break_tet_to_tet_N_5_EdgeBased)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        if (LARGE_TEST_ONLY || !DO_TESTS) return;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            // start_demo_local_refiner_break_tet_to_tet_2

            percept::PerceptMesh eMesh;
            eMesh.open(input_files_loc+"cylinder_with_5_holes.e");

            Local_Tet4_Tet4_N break_tet_to_tet_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.add_field("proc_rank", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);
            stk::mesh::FieldBase* refine_field = eMesh.add_field("refine_field", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);
            eMesh.commit();

            if (0)
              {
                SetRefineField set_ref_field(eMesh);
                eMesh.elementOpLoop(set_ref_field, refine_field);
              }

            eMesh.save_as( output_files_loc+"local_tet_N_5_EdgeBased_0_"+post_fix[p_size]+".e");

            MyEdgeBasedRefinePredicate mrp(0, refine_field, 0.0);

            PredicateBasedEdgeAdapter<MyEdgeBasedRefinePredicate>
              breaker(mrp,
                      eMesh, break_tet_to_tet_N, proc_rank_field);

            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);
            for (int ipass=0; ipass < 5; ipass++)
              {
                //eMesh.elementOpLoop(set_ref_field, refine_field);

                std::cout << "P[" << eMesh.get_rank() << "] ipass= " << ipass << std::endl;
                breaker.doBreak();
                std::cout << "P[" << eMesh.get_rank() << "] done... ipass= " << ipass << std::endl;
                eMesh.save_as(output_files_loc+"local_tet_N_5_EdgeBased_1_ipass_"+toString(ipass)+"_"+post_fix[p_size]+".e");
              }

            //breaker.deleteParentElements();
            eMesh.save_as(output_files_loc+"local_tet_N_5_EdgeBased_1_"+post_fix[p_size]+".e");

#if 1
            for (int iunref_pass=0; iunref_pass < 5; iunref_pass++)
              {
                std::cout << "P[" << eMesh.get_rank() << "] iunref_pass= " << iunref_pass << std::endl;
                ElementUnrefineCollection elements_to_unref = breaker.buildUnrefineList();
                breaker.unrefineTheseElements(elements_to_unref);
                eMesh.save_as(output_files_loc+"local_tet_N_5_EdgeBased_1_unref_ipass_"+toString(iunref_pass)+"_"+post_fix[p_size]+".e");
              }

            eMesh.save_as( output_files_loc+"local_tet_N_5_EdgeBased_1_unref_"+post_fix[p_size]+".e");
#endif
            // end_demo
          }

      }

      STKUNIT_UNIT_TEST(regr_localRefiner, break_tet_to_tet_N_5_EdgeBased_shock)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        if (LARGE_TEST_ONLY || !DO_TESTS) return;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            // start_demo_local_refiner_break_tet_to_tet_2

            percept::PerceptMesh eMesh;
            eMesh.open(input_files_loc+"cylinder_with_5_holes.e");

            Local_Tet4_Tet4_N break_tet_to_tet_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.add_field("proc_rank", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);
            stk::mesh::FieldBase* refine_field = eMesh.add_field("refine_field", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);
            stk::mesh::FieldBase* nodal_refine_field = eMesh.add_field("nodal_refine_field", eMesh.node_rank(), scalarDimension);
            eMesh.commit();

            if (0)
              {
                SetRefineField set_ref_field(eMesh);
                eMesh.elementOpLoop(set_ref_field, refine_field);
              }

            eMesh.save_as( output_files_loc+"local_tet_N_5_EdgeBased_shock_0_"+post_fix[p_size]+".e");

            PlaneShock shock;

            ShockBasedRefinePredicate srp(nodal_refine_field, eMesh, 0, refine_field, 0.0, shock, 0.0);

            PredicateBasedEdgeAdapter<ShockBasedRefinePredicate>
              breaker(srp,
                      eMesh, break_tet_to_tet_N, proc_rank_field);

            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);

            for (int ipass=0; ipass < 3; ipass++)
              {
                //eMesh.elementOpLoop(set_ref_field, refine_field);

                std::cout << "P[" << eMesh.get_rank() << "] ipass= " << ipass << std::endl;
                breaker.doBreak();
                std::cout << "P[" << eMesh.get_rank() << "] done... ipass= " << ipass << std::endl;
                eMesh.save_as(output_files_loc+"local_tet_N_5_EdgeBased_shock_1_ipass_"+toString(ipass)+"_"+post_fix[p_size]+".e");
              }

            breaker.deleteParentElements();
            eMesh.save_as(output_files_loc+"local_tet_N_5_EdgeBased_shock_1_"+post_fix[p_size]+".e");

#if 0
            for (int iunref_pass=0; iunref_pass < 4; iunref_pass++)
              {
                std::cout << "P[" << eMesh.get_rank() << "] iunref_pass= " << iunref_pass << std::endl;
                ElementUnrefineCollection elements_to_unref = breaker.buildUnrefineList();
                breaker.unrefineTheseElements(elements_to_unref);
                eMesh.save_as(output_files_loc+"local_tet_N_5_EdgeBased_shock_1_unref_ipass_"+toString(iunref_pass)+"_"+post_fix[p_size]+".e");
              }

            eMesh.save_as( output_files_loc+"local_tet_N_5_EdgeBased_shock_1_unref_"+post_fix[p_size]+".e");
#endif
            // end_demo
          }

      }

      //=============================================================================
      //=============================================================================
      //=============================================================================


      static void do_moving_shock_test(int num_time_steps, bool save_intermediate=false, bool delete_parents=false)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        shock_width = 1./5.0;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            // start_demo_local_refiner_break_tet_to_tet_2

            percept::PerceptMesh eMesh;
            eMesh.open(input_files_loc+"cylinder_with_5_holes.e");

            Local_Tet4_Tet4_N break_tet_to_tet_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.add_field("proc_rank", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);
            stk::mesh::FieldBase* refine_field = eMesh.add_field("refine_field", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);
            stk::mesh::FieldBase* nodal_refine_field = eMesh.add_field("nodal_refine_field", eMesh.node_rank(), scalarDimension);
            eMesh.commit();

            std::cout << "moving_shock initial number elements= " << eMesh.get_number_elements() << std::endl;

            eMesh.save_as( output_files_loc+"moving_shock_"+post_fix[p_size]+".e.0");

            PlaneShock shock;
            shock.plane_point_init[0] = 2.0;
            shock.plane_point_init[1] = 0.0;
            shock.plane_point_init[2] = 0.0;
            shock.plane_normal[0] = 1;
            shock.plane_normal[1] = 0;
            shock.plane_normal[2] = 0;

            ShockBasedRefinePredicate srp(nodal_refine_field, eMesh, 0, refine_field, 0.0, shock, 0.0, 0.4);

            PredicateBasedEdgeAdapter<ShockBasedRefinePredicate>
              breaker(srp,
                      eMesh, break_tet_to_tet_N, proc_rank_field);

            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);

            double delta_shock_displacement = 0.2;
            double shock_displacement = -2.0;
            int num_ref_passes = 2;
            int num_unref_passes = 3;

            for (int istep = 0; istep < num_time_steps; istep++)
              {
                std::cout << "P[" << eMesh.get_rank() << "] istep= " << istep << std::endl;

                breaker.getRefinePredicate().m_shock_displacement = shock_displacement;

                for (int ipass=0; ipass < num_ref_passes; ipass++)
                  {
                    std::cout << "P[" << eMesh.get_rank() << "] ipass= " << ipass << std::endl;
                    breaker.doBreak();
                    std::cout << "P[" << eMesh.get_rank() << "] done... ipass= " << ipass << std::endl;
                    if (save_intermediate)
                      eMesh.save_as(output_files_loc+"tmp_moving_shock_ref_istep_ipass_"+toString(istep)+"_"+toString(ipass)+"_"+post_fix[p_size]+".e");
                  }

                //breaker.getNodeRegistry().init_entity_repo();
                for (int iunref_pass=0; iunref_pass < num_unref_passes; iunref_pass++)
                  {
                    ElementUnrefineCollection elements_to_unref = breaker.buildUnrefineList();
                    std::cout << "P[" << eMesh.get_rank() << "] iunref_pass= " << iunref_pass << " unref list size= " << elements_to_unref.size() << std::endl;
                    breaker.unrefineTheseElements(elements_to_unref);
                    if (save_intermediate)
                      eMesh.save_as(output_files_loc+"tmp_moving_shock_unref_istep_ipass_"+toString(istep)+"_"+toString(iunref_pass)+"_"+post_fix[p_size]+".e");
                  }

                if (delete_parents && istep == num_time_steps-1)
                  {
                    breaker.deleteParentElements();
                  }
                if (istep == num_time_steps-1 || save_intermediate)
                  eMesh.save_as(output_files_loc+"moving_shock_"+post_fix[p_size]+".e."+toString(istep+1) );

                shock_displacement += delta_shock_displacement;

              }

            eMesh.save_as(output_files_loc+"final_moving_shock_"+post_fix[p_size]+".e."+toString(num_time_steps) );
            for (int iunref=0; iunref < 10; iunref++)
              {
                breaker.unrefineAll();
              }
            breaker.deleteParentElements();
            std::cout << "moving_shock final number elements= " << eMesh.get_number_elements() << std::endl;
            eMesh.save_as(output_files_loc+"final_unrefed_moving_shock_"+post_fix[p_size]+".e."+toString(num_time_steps) );

            // end_demo
          }

      }

      STKUNIT_UNIT_TEST(regr_localRefiner, break_tet_to_tet_N_5_EdgeBased_moving_shock)
      {
        //if (1) return;
        const bool do_full_demo = false;
        if (LARGE_TEST_ONLY || !DO_TESTS) return;
        if (do_full_demo)
          {
            int num_time_steps = 10;  // 10 for stress testing
            for (int istep=1; istep <= num_time_steps; istep++)
              do_moving_shock_test(istep, false, true);
          }
        else
          // normal regression testing
          {
            int num_time_steps = 3;  // 10 for stress testing
            bool save_intermediate=true;
            do_moving_shock_test(num_time_steps, save_intermediate);
          }
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================


      static void do_moving_shock_test_cyl_sidesets(int num_time_steps, bool save_intermediate=false, bool delete_parents=false)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        shock_width = 1./25.0;
        double shock_diff_criterion = 0.04;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size==1 || p_size == 3)
          {
            percept::PerceptMesh eMesh;
            eMesh.open(input_files_loc+"cylinder_tet4_0.e");

            Local_Tet4_Tet4_N break_tet_to_tet_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field    = eMesh.add_field("proc_rank", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);
            stk::mesh::FieldBase* refine_field       = eMesh.add_field("refine_field", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);
            stk::mesh::FieldBase* nodal_refine_field = eMesh.add_field("nodal_refine_field", eMesh.node_rank(), scalarDimension);
            eMesh.commit();
            //eMesh.delete_side_sets();

            std::cout << "moving_shock initial number elements= " << eMesh.get_number_elements() << std::endl;

            eMesh.save_as( output_files_loc+"cyl_sidesets_moving_shock_"+post_fix[p_size]+".e.0");

            PlaneShock shock;
            shock.plane_point_init[0] = 0.0;
            shock.plane_point_init[1] = 0.0;
            shock.plane_point_init[2] = 0.0;
            shock.plane_normal[0] = 1;
            shock.plane_normal[1] = 0;
            shock.plane_normal[2] = 0;

            ShockBasedRefinePredicate srp(nodal_refine_field, eMesh, 0, refine_field, 0.0, shock, 0.0, shock_diff_criterion);

            PredicateBasedEdgeAdapter<ShockBasedRefinePredicate>
              breaker(srp,
                      eMesh, break_tet_to_tet_N, proc_rank_field);

            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);

            // x,y,z: [-.08,.08],[-.08,.08],[-.5,-.2]
            double delta_shock_displacement = 0.2;
            double shock_displacement = 0.0;
            int num_ref_passes = 3;
            int num_unref_passes = 3;

            for (int istep = 0; istep < num_time_steps; istep++)
              {
                std::cout << "P[" << eMesh.get_rank() << "] istep= " << istep << std::endl;

                breaker.getRefinePredicate().m_shock_displacement = shock_displacement;

                for (int ipass=0; ipass < num_ref_passes; ipass++)
                  {
                    std::cout << "P[" << eMesh.get_rank() << "] ipass= " << ipass <<  std::endl;
                    breaker.doBreak();
                    std::cout << "P[" << eMesh.get_rank() << "] done... ipass= " << ipass << " moving_shock number elements= " << eMesh.get_number_elements() << std::endl;
                  }

                //breaker.deleteParentElements();
                eMesh.save_as(output_files_loc+"tmp_cyl_sidesets_moving_shock_"+post_fix[p_size]+".e");
                //exit(123);

                //breaker.getNodeRegistry().init_entity_repo();
                for (int iunref_pass=0; iunref_pass < num_unref_passes; iunref_pass++)
                  {
                    std::cout << "P[" << eMesh.get_rank() << "] iunref_pass= " << iunref_pass <<  std::endl;
                    ElementUnrefineCollection elements_to_unref = breaker.buildUnrefineList();
                    breaker.unrefineTheseElements(elements_to_unref);
                    std::cout << "P[" << eMesh.get_rank() << "] done... iunref_pass= " << iunref_pass << " moving_shock number elements= " << eMesh.get_number_elements() << std::endl;
                  }

//                 breaker.deleteParentElements();
                 eMesh.save_as(output_files_loc+"tmp_cyl_sidesets_moving_shock_"+post_fix[p_size]+".e");
                 //exit(123);

                if (delete_parents && istep == num_time_steps-1)
                  {
                    breaker.deleteParentElements();
                  }
                if (istep == num_time_steps-1 || save_intermediate)
                  eMesh.save_as(output_files_loc+"cyl_sidesets_moving_shock_"+post_fix[p_size]+".e."+toString(istep+1) );

                shock_displacement += delta_shock_displacement;

              }

            eMesh.save_as(output_files_loc+"cyl_sidesets_final_moving_shock_"+post_fix[p_size]+".e."+toString(num_time_steps) );
            for (int iunref=0; iunref < 10; iunref++)
              {
                std::cout << "P[" << eMesh.get_rank() << "] iunrefAll_pass= " << iunref <<  std::endl;
                breaker.unrefineAll();
                std::cout << "P[" << eMesh.get_rank() << "] done... iunrefAll_pass= " << iunref << " moving_shock number elements= " << eMesh.get_number_elements() << std::endl;
                eMesh.save_as(output_files_loc+"cyl_sidesets_final_moving_shock_"+post_fix[p_size]+"_unrefAll_pass_"+toString(iunref)+".e."+toString(num_time_steps) );
              }

            if (0)
              breaker.deleteParentElements();
            std::cout << "moving_shock final number elements= " << eMesh.get_number_elements() << std::endl;
            eMesh.save_as(output_files_loc+"cyl_sidesets_final_unrefed_moving_shock_"+post_fix[p_size]+".e."+toString(num_time_steps) );
            //exit(123);

            // end_demo
          }
      }

      STKUNIT_UNIT_TEST(regr_localRefiner, break_tet_to_tet_N_5_EdgeBased_moving_shock_cyl_sidesets)
      {
        const bool do_full_demo = false;
        if (LARGE_TEST_ONLY || !DO_TESTS) return;
        if (do_full_demo)
          {
            int num_time_steps = 10;  // 10 for stress testing
            for (int istep=1; istep <= num_time_steps; istep++)
              do_moving_shock_test_cyl_sidesets(istep, false, true);
          }
        else
          // normal regression testing
          {
            int num_time_steps = 1;  // 10 for stress testing
            do_moving_shock_test_cyl_sidesets(num_time_steps);
          }
      }


      static void do_moving_shock_test_square_sidesets(PerceptMesh& eMesh, int num_time_steps, bool save_intermediate=false, bool delete_parents=false)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        shock_width = 1./25.0;
        double shock_diff_criterion = 0.04;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size==1 || p_size == 3)
          {
            Local_Tri3_Tri3_N break_tri_to_tri_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field    = eMesh.add_field("proc_rank", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);
            stk::mesh::FieldBase* refine_field       = eMesh.add_field("refine_field", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);
            stk::mesh::FieldBase* nodal_refine_field = eMesh.add_field("nodal_refine_field", eMesh.node_rank(), scalarDimension);
            eMesh.commit();
            //eMesh.delete_side_sets();

            std::cout << "moving_shock initial number elements= " << eMesh.get_number_elements() << std::endl;

            eMesh.save_as( output_files_loc+"square_sidesets_moving_shock_"+post_fix[p_size]+".e.0");

            PlaneShock shock;
            shock.plane_point_init[0] = 0.0;
            shock.plane_point_init[1] = 0.0;
            shock.plane_point_init[2] = 0.0;
            shock.plane_normal[0] = 1;
            shock.plane_normal[1] = 0;
            shock.plane_normal[2] = 0;

            ShockBasedRefinePredicate srp(nodal_refine_field, eMesh, 0, refine_field, 0.0, shock, 0.0, shock_diff_criterion);

            PredicateBasedEdgeAdapter<ShockBasedRefinePredicate>
              breaker(srp,
                      eMesh, break_tri_to_tri_N, proc_rank_field);

            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);

            // x,y,z: [-.08,.08],[-.08,.08],[-.5,-.2]
            double delta_shock_displacement = 0.2;
            double shock_displacement = 0.0;
            int num_ref_passes = 3;
            int num_unref_passes = 3;

            for (int istep = 0; istep < num_time_steps; istep++)
              {
                std::cout << "P[" << eMesh.get_rank() << "] istep= " << istep << std::endl;

                breaker.getRefinePredicate().m_shock_displacement = shock_displacement;

                for (int ipass=0; ipass < num_ref_passes; ipass++)
                  {
                    std::cout << "P[" << eMesh.get_rank() << "] ipass= " << ipass <<  std::endl;
                    breaker.doBreak();
                    std::cout << "P[" << eMesh.get_rank() << "] done... ipass= " << ipass << " moving_shock number elements= " << eMesh.get_number_elements() << std::endl;
                  }

                //breaker.deleteParentElements();
                eMesh.save_as(output_files_loc+"tmp_square_sidesets_moving_shock_ref_"+post_fix[p_size]+".e.s-"+toString(istep+1));
                //eMesh.save_as("square_anim."+toString(istep+1)+".e");
                char buf[1000];
                sprintf(buf, "%04d", istep);
                if (istep==0)
                  eMesh.save_as("square_anim.e");
                else
                  eMesh.save_as("square_anim.e-s"+std::string(buf));
                //eMesh.save_as("square_anim.e-s"+toString(istep));
                //eMesh.dump_vtk("square_anim."+toString(istep+1)+".vtk", false);

                //exit(123);

                //breaker.getNodeRegistry().init_entity_repo();
                for (int iunref_pass=0; iunref_pass < num_unref_passes; iunref_pass++)
                  {
                    std::cout << "P[" << eMesh.get_rank() << "] iunref_pass= " << iunref_pass <<  std::endl;
                    ElementUnrefineCollection elements_to_unref = breaker.buildUnrefineList();
                    breaker.unrefineTheseElements(elements_to_unref);
                    std::cout << "P[" << eMesh.get_rank() << "] done... iunref_pass= " << iunref_pass << " moving_shock number elements= " << eMesh.get_number_elements() << std::endl;
                  }

                //                 breaker.deleteParentElements();
                eMesh.save_as(output_files_loc+"tmp_square_sidesets_moving_shock_unref_"+post_fix[p_size]+".e");

                if (delete_parents && istep == num_time_steps-1)
                  {
                    breaker.deleteParentElements();
                  }
                if (istep == num_time_steps-1 || save_intermediate)
                  eMesh.save_as(output_files_loc+"square_sidesets_moving_shock_"+post_fix[p_size]+".e.s-"+toString(istep+1) );

                shock_displacement += delta_shock_displacement;

              }

            eMesh.save_as(output_files_loc+"square_sidesets_final_moving_shock_"+post_fix[p_size]+".e.s-"+toString(num_time_steps) );
            for (int iunref=0; iunref < 10; iunref++)
              {
                std::cout << "P[" << eMesh.get_rank() << "] iunrefAll_pass= " << iunref <<  std::endl;
                breaker.unrefineAll();
                std::cout << "P[" << eMesh.get_rank() << "] done... iunrefAll_pass= " << iunref << " moving_shock number elements= " << eMesh.get_number_elements() << std::endl;
                eMesh.save_as(output_files_loc+"square_sidesets_final_moving_shock_"+post_fix[p_size]+"_unrefAll_pass_"+toString(iunref)+".e."+toString(num_time_steps) );
              }

            if (1)
              breaker.deleteParentElements();
            std::cout << "moving_shock final number elements= " << eMesh.get_number_elements() << std::endl;
            eMesh.save_as(output_files_loc+"square_sidesets_final_unrefed_moving_shock_"+post_fix[p_size]+".e."+toString(num_time_steps) );
            //exit(123);

            // end_demo
          }
      }

      STKUNIT_UNIT_TEST(regr_localRefiner, break_tri_to_tri_N_5_EdgeBased_moving_shock_square_sidesets)
      {
        bool do_test=false;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        if (do_test) {

          // if this is a fresh installation, set to true to generate the initial meshes needed for this test (once only)
          bool do_bootstrap_mesh = false;
          if (do_bootstrap_mesh)
          {
            const unsigned n = 5;
            const unsigned nx = n , ny = n;

            bool createEdgeSets = true;
            percept::QuadFixture<double, shards::Triangle<3> > fixture( pm , nx , ny, createEdgeSets);
            fixture.set_bounding_box(-1,1, -1, 1);

            bool isCommitted = false;
            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, isCommitted);

            eMesh.commit();

            fixture.generate_mesh();
            eMesh.save_as(input_files_loc+"square_tri3_0.e");
            //eMesh.print_info("test1",2);
            //eMesh.dump_vtk("sqtri3.vtk",false);
          }
          //if (1) return;

          {
            PerceptMesh eMesh;
            eMesh.open(input_files_loc+"square_tri3_0.e");

            //int num_time_steps = 10;  // 10 for stress testing
            //for (int istep=1; istep <= num_time_steps; istep++)
            do_moving_shock_test_square_sidesets(eMesh, 10, false, true);
          }

          if (do_bootstrap_mesh)
          {
            PerceptMesh eMesh1;
            eMesh1.open(input_files_loc+"square_tri3_uns.e");
            eMesh1.commit();
            MDArray xform(3,3);
            xform(0,0) = 1./5.;
            xform(1,1) = 1./5.;
            xform(2,2) = 1.0;
            eMesh1.transform_mesh(xform);
            eMesh1.save_as(input_files_loc+"square_tri3_uns_xformed.e");
          }

          {
            PerceptMesh eMesh2;
            eMesh2.open(input_files_loc+"square_tri3_uns_xformed.e");
            //int num_time_steps = 10;  // 10 for stress testing
            //for (int istep=1; istep <= num_time_steps; istep++)
            do_moving_shock_test_square_sidesets(eMesh2, 10, false, true);
          }
        }


      }

    }
  }
}
