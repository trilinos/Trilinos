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
#include <stk_percept/fixtures/SingleTetFixture.hpp>

#include <stk_adapt/IEdgeBasedAdapterPredicate.hpp>
#include <stk_adapt/ElementRefinePredicate.hpp>
#include <stk_adapt/PredicateBasedElementAdapter.hpp>
#include <stk_adapt/PredicateBasedEdgeAdapter.hpp>

#include <stk_percept/function/ElementOp.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/CreateEdges.hpp>
#include <stk_mesh/base/Comm.hpp>

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

      class SetElementRefineFieldValue : public percept::ElementOp
      {
        percept::PerceptMesh& m_eMesh;
        double m_val;
      public:
        SetElementRefineFieldValue(percept::PerceptMesh& eMesh, double val=0.0) : m_eMesh(eMesh), m_val(val) {
        }

        virtual bool operator()(const stk::mesh::Entity element, stk::mesh::FieldBase *field,  const mesh::BulkData& bulkData)
        {
          double *f_data = m_eMesh.field_data(field, element);
          f_data[0] = m_val;

          return false;  // don't terminate the loop
        }
        virtual void init_elementOp() {}
        virtual void fini_elementOp() {}

      };



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

            ElementRefinePredicate erp(eMesh, 0, refine_field, 0.0);

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
#ifdef NDEBUG
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

            MyEdgeBasedRefinePredicate mrp(eMesh, 0, refine_field, 0.0);

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

#endif
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
            eMesh.set_ioss_write_options("large");

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
                std::cout << "moving shock unrefineAll pass= " << iunref << std::endl;
                breaker.unrefineAll();
              }
            std::cout << "moving shock deleteParentElements = "  << std::endl;
            breaker.deleteParentElements();
            std::cout << "moving_shock final number elements= " << eMesh.get_number_elements() << std::endl;
            eMesh.save_as(output_files_loc+"final_unrefed_moving_shock_"+post_fix[p_size]+".e."+toString(num_time_steps) );

            // end_demo
          }

      }

      STKUNIT_UNIT_TEST(regr_localRefiner, break_tet_to_tet_N_5_EdgeBased_moving_shock)
      {
#ifndef NDEBUG
        bool no_test=true;
        if (no_test) return;
#endif
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


      static void do_moving_shock_test_cyl_sidesets(std::string filename, int num_time_steps, bool save_intermediate=false, bool delete_parents=false)
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
            //eMesh.open(input_files_loc+"cylinder_tet4_0.e");
            eMesh.open(input_files_loc+filename);

            Local_Tet4_Tet4_N break_tet_to_tet_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field       = eMesh.add_field("proc_rank", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);
            stk::mesh::FieldBase* refine_field          = eMesh.add_field("refine_field", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);
            stk::mesh::FieldBase* refine_field_filtered = eMesh.add_field("refine_field_filtered", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);
            (void)refine_field_filtered;
            stk::mesh::FieldBase* nodal_refine_field    = eMesh.add_field("nodal_refine_field", eMesh.node_rank(), scalarDimension);

            stk::mesh::FieldBase* refine_level_d     = eMesh.add_field("refine_level_d", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);
            ScalarIntFieldType& refine_level       = eMesh.get_fem_meta_data()->declare_field<ScalarIntFieldType>("refine_level");
            stk::mesh::put_field( refine_level , stk::mesh::MetaData::ELEMENT_RANK , eMesh.get_fem_meta_data()->universal_part());
            stk::io::set_field_role(refine_level, Ioss::Field::TRANSIENT);

            eMesh.commit();
            eMesh.delete_side_sets();
            eMesh.output_active_children_only(true);
            if (1)
              {
                SetElementRefineFieldValue set_ref_field(eMesh);
                eMesh.elementOpLoop(set_ref_field, refine_level_d);
              }

            std::cout << "moving_shock initial number elements= " << eMesh.get_number_elements() << std::endl;

            eMesh.save_as( output_files_loc+"cyl_sidesets_moving_shock_"+post_fix[p_size]+".e.0");

            PlaneShock shock;
            shock.plane_point_init[0] = -0.8;
            shock.plane_point_init[1] = 0.0;
            shock.plane_point_init[2] = 0.0;
            shock.plane_normal[0] = 1;
            shock.plane_normal[1] = 0;
            shock.plane_normal[2] = 0;

            ShockBasedRefinePredicate1 srp(nodal_refine_field, eMesh, 0, refine_field, 0.0, shock, 0.0, shock_diff_criterion);

            PredicateBasedEdgeAdapter<ShockBasedRefinePredicate1>
              breaker(srp,
                      eMesh, break_tet_to_tet_N, proc_rank_field);

            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);
            breaker.setDoLevelBasedUnrefinement(true);

            // x,y,z: [-.08,.08],[-.08,.08],[-.5,-.2]
            double delta_shock_displacement = 0.02;
            double shock_displacement = 0.0;
            int num_ref_passes = 1;
            int num_unref_passes = 1;

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

                if (1)
                  {
                    //SetRefineField set_ref_field(eMesh);
                    //eMesh.elementOpLoop(srp, refine_field);
                    SetRefineField1 set_ref_field(breaker);
                    eMesh.elementOpLoop(set_ref_field, refine_field);
                  }

                eMesh.save_as(output_files_loc+"tmp_cyl_sidesets_moving_shock_ref_"+post_fix[p_size]+"_step_"+toString(istep)+".e");

                for (int iunref_pass=0; iunref_pass < num_unref_passes; iunref_pass++)
                  {
                    std::cout << "P[" << eMesh.get_rank() << "] iunref_pass= " << iunref_pass << " number elements before unrefine= " << eMesh.get_number_elements() <<  std::endl;
                    ElementUnrefineCollection elements_to_unref = breaker.buildUnrefineList();
                    size_t unref_sz = elements_to_unref.size();
                    //breaker.unrefineTheseElements(elements_to_unref);
                    breaker.unrefinePass2(elements_to_unref);

                    std::cout << "P[" << eMesh.get_rank() << "] done... iunref_pass= " << iunref_pass << " unref list size = " << unref_sz << " moving_shock number elements= " << eMesh.get_number_elements() << std::endl;
                  }

                eMesh.save_as(output_files_loc+"tmp_cyl_sidesets_moving_shock_unref_"+post_fix[p_size]+"_step_"+toString(istep)+".e");

                if (delete_parents && istep == num_time_steps-1)
                  {
                    breaker.deleteParentElements();
                  }
                if (istep == num_time_steps-1 || save_intermediate)
                  eMesh.save_as(output_files_loc+"cyl_sidesets_moving_shock_"+post_fix[p_size]+".e."+toString(istep+1) );

                if (1)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", istep);
                    if (istep==0)
                      eMesh.save_as("cyl_tet_anim.e");
                    else
                      eMesh.save_as("cyl_tet_anim.e-s"+std::string(buf));
                  }

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

            // end_demo
          }
      }

      STKUNIT_UNIT_TEST(regr_localRefiner, break_tet_to_tet_N_5_EdgeBased_moving_shock_cyl_sidesets)
      {
        const bool do_full_demo = false;
#ifndef NDEBUG
        bool no_test=true;
        if (no_test) return;
#endif
        if (LARGE_TEST_ONLY || !DO_TESTS) return;
        std::string filename = "cylinder_tet4_0.e";
        if (do_full_demo)
          {
            int num_time_steps = 10;  // 10 for stress testing
            for (int istep=1; istep <= num_time_steps; istep++)
              do_moving_shock_test_cyl_sidesets(filename, istep, false, true);
          }
        else
          // normal regression testing
          {
            int num_time_steps = 10;  // 10 for stress testing
            do_moving_shock_test_cyl_sidesets(filename, num_time_steps);
          }
      }

      class SetElementField : public percept::ElementOp
      {
        percept::PerceptMesh& m_eMesh;
        IEdgeAdapter *m_iea;
      public:
        SetElementField(percept::PerceptMesh& eMesh, IEdgeAdapter* iea) : m_eMesh(eMesh), m_iea(iea) {
        }

        virtual bool operator()(const stk::mesh::Entity element, stk::mesh::FieldBase *field,  const mesh::BulkData& bulkData)
        {
          double *f_data = m_eMesh.field_data(field, element);
          int unref = m_iea->markUnrefine(element);
          int ref_count = m_iea->markCountRefinedEdges(element);
          f_data[0] = (double)((unref & DO_UNREFINE ? -1 : 0));
          f_data[0] = (double)(ref_count?1.0:f_data[0]);


          return false;  // don't terminate the loop
        }
        virtual void init_elementOp() {}
        virtual void fini_elementOp() {}

      };



      template< class Breaker >
      static void do_moving_shock_test_square_sidesets(PerceptMesh& eMesh, int num_time_steps, bool save_intermediate=false, bool delete_parents=false,
                                                       std::string prefix="")
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        shock_width = 1./25.0;
        double shock_diff_criterion = 0.04;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size==1 || p_size == 3)
          {
            Breaker break_tri_to_tri_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field    = eMesh.add_field("proc_rank", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);

            // just some debugging fields, not necessary
            stk::mesh::FieldBase* refine_field       = eMesh.add_field("refine_field", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);
            stk::mesh::FieldBase* nodal_refine_field = eMesh.add_field("nodal_refine_field", eMesh.node_rank(), scalarDimension);
            eMesh.add_field("normal_kept_deleted", eMesh.node_rank(), scalarDimension);
            eMesh.add_field("refine_field_filtered", eMesh.element_rank(), scalarDimension);

            stk::mesh::FieldBase* refine_level_d     = eMesh.add_field("refine_level_d", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);
            ScalarIntFieldType& refine_level       = eMesh.get_fem_meta_data()->declare_field<ScalarIntFieldType>("refine_level");
            stk::mesh::put_field( refine_level , stk::mesh::MetaData::ELEMENT_RANK , eMesh.get_fem_meta_data()->universal_part());
            stk::io::set_field_role(refine_level, Ioss::Field::TRANSIENT);

            eMesh.commit();
            eMesh.output_active_children_only(true);
            if (1)
              {
                SetElementRefineFieldValue set_ref_field(eMesh);
                eMesh.elementOpLoop(set_ref_field, refine_level_d);
              }

            //eMesh.delete_side_sets();
            std::cout << "moving_shock initial number elements= " << eMesh.get_number_elements() << std::endl;

            eMesh.save_as( output_files_loc+prefix+"square_sidesets_moving_shock_"+post_fix[p_size]+".e.0");

            PlaneShock shock;
            shock.plane_point_init[0] = -0.8;
            shock.plane_point_init[1] = 0.0;
            shock.plane_point_init[2] = 0.0;
            shock.plane_normal[0] = 1;
            shock.plane_normal[1] = 0;
            shock.plane_normal[2] = 0;

            ShockBasedRefinePredicate1 srp(nodal_refine_field, eMesh, 0, refine_field, 0.0, shock, 0.0, shock_diff_criterion);

#define USE_EDGE_BASED 1

#if USE_EDGE_BASED
            PredicateBasedEdgeAdapter<ShockBasedRefinePredicate1>
              breaker(srp,
                      eMesh, break_tri_to_tri_N, proc_rank_field);
#else
            PredicateBasedElementAdapter<ShockBasedRefinePredicate>
              breaker(srp,
                      eMesh, break_tri_to_tri_N, proc_rank_field);
#endif

            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);
            //breaker.setDoLevelBasedUnrefinement(true);

            // x,y,z: [-.08,.08],[-.08,.08],[-.5,-.2]
            double delta_shock_displacement = 0.02;
            double shock_displacement = 0.0;
            int num_ref_passes = 3;
            int num_unref_passes = 1;
            static int cnt1=0;

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

#if USE_EDGE_BASED
                if (0)
                {
                  SetElementField set_ref_field(eMesh, &breaker);
                  eMesh.elementOpLoop(set_ref_field, refine_field);
                }
#endif

                eMesh.save_as(output_files_loc+prefix+"tmp_square_sidesets_moving_shock_ref_"+post_fix[p_size]+".e.s-"+toString(istep+1));

                if (1)
                  {
                    eMesh.nodalOpLoop(srp, nodal_refine_field);
                  }

                if (0)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", cnt1);
                    if (cnt1==0)
                      eMesh.save_as(prefix+"ref.e");
                    else
                      eMesh.save_as(prefix+"ref.e-s"+std::string(buf));
                    ++cnt1;
                  }

                if (0)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", istep);
                    if (istep==0)
                      eMesh.save_as(prefix+"square_anim.e");
                    else
                      eMesh.save_as(prefix+"square_anim.e-s"+std::string(buf));
                  }

                for (int iunref_pass=0; iunref_pass < num_unref_passes; iunref_pass++)
                  {
                    std::cout << "P[" << eMesh.get_rank() << "] iunref_pass= " << iunref_pass <<  std::endl;

                    if (istep == 12)
                      {
                        Util::setFlag(1256,true);
                      }
                    if (istep==0)
                      {
                        Util::setFlag(1257,true);
                      }

                    if (0)
                      {
                        ElementUnrefineCollection elements_to_unref = breaker.buildUnrefineList();
                        breaker.unrefineTheseElements(elements_to_unref);
                      }

                    if (1)
                      {
                        ElementUnrefineCollection elements_to_unref = breaker.buildUnrefineList();
                        breaker.unrefinePass2(elements_to_unref);
                      }

                    if (istep==11)
                      {
                        char buf1[1000];
                        sprintf(buf1, "%04d", iunref_pass);
                        eMesh.save_as("test1.e"+(iunref_pass==0?"":"-s"+std::string(buf1)));
                        if (iunref_pass==15) {
                          std::cout << "shock_displacement= " << shock_displacement << std::endl;
                          //exit(1);
                        }
                      }

                    std::cout << "P[" << eMesh.get_rank() << "] done... iunref_pass= " << iunref_pass << " moving_shock number elements= " << eMesh.get_number_elements() << std::endl;
                  }
                if (1)
                  {
                    eMesh.nodalOpLoop(srp, nodal_refine_field);
                  }


                if (0)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", cnt1);
                    if (cnt1==0)
                      eMesh.save_as(prefix+"ref.e");
                    else
                      eMesh.save_as(prefix+"ref.e-s"+std::string(buf));
                    ++cnt1;
                  }

                if (1)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", istep);
                    if (istep==0)
                      eMesh.save_as(prefix+"square_anim.e");
                    else
                      eMesh.save_as(prefix+"square_anim.e-s"+std::string(buf));
                  }

                eMesh.save_as(output_files_loc+prefix+"tmp_square_sidesets_moving_shock_unref_"+post_fix[p_size]+".e");

                if (delete_parents && istep == num_time_steps-1)
                  {
                    breaker.deleteParentElements();
                  }
                if (istep == num_time_steps-1 || save_intermediate)
                  eMesh.save_as(output_files_loc+prefix+"square_sidesets_moving_shock_"+post_fix[p_size]+".e.s-"+toString(istep+1) );

                shock_displacement += delta_shock_displacement;

              }

            eMesh.save_as(output_files_loc+prefix+"square_sidesets_final_moving_shock_"+post_fix[p_size]+".e.s-"+toString(num_time_steps) );
            for (int iunref=0; iunref < 10; iunref++)
              {
                std::cout << "P[" << eMesh.get_rank() << "] iunrefAll_pass= " << iunref <<  std::endl;
                eMesh.save_as(output_files_loc+prefix+"square_sidesets_final_moving_shock_"+post_fix[p_size]+"_unrefAll_pass_"+toString(iunref)+".e."+toString(num_time_steps) );
                breaker.unrefineAll();
                std::cout << "P[" << eMesh.get_rank() << "] done... iunrefAll_pass= " << iunref << " moving_shock number elements= " << eMesh.get_number_elements() << std::endl;
              }

            if (0)
              breaker.deleteParentElements();
            std::cout << "moving_shock final number elements= " << eMesh.get_number_elements() << std::endl;
            eMesh.save_as(output_files_loc+prefix+"square_sidesets_final_unrefed_moving_shock_"+post_fix[p_size]+".e."+toString(num_time_steps) );

            // end_demo
          }
      }

      STKUNIT_UNIT_TEST(regr_localRefiner, break_tri_to_tri_N_5_EdgeBased_moving_shock_square_sidesets)
      {
        bool do_test=true;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        bool do_demo=false; // does more timesteps

        int num_time_steps = (do_demo?80:10);

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
          }

          // structured mesh
          if (1)
          {
            PerceptMesh eMesh;
            eMesh.open(input_files_loc+"square_tri3_0.e");
            do_moving_shock_test_square_sidesets<Local_Tri3_Tri3_N>(eMesh, num_time_steps, false, true, "str-");
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

          // unstructured mesh
          if (1)
          {
            PerceptMesh eMesh2;
            eMesh2.open(input_files_loc+"square_tri3_uns_xformed.e");
            do_moving_shock_test_square_sidesets<Local_Tri3_Tri3_N>(eMesh2, num_time_steps, false, false, "uns-");
          }
        }
      }

      static void do_tet_edge_test(unsigned ntets, bool do_part_for_edges)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        if (LARGE_TEST_ONLY || !DO_TESTS) return;
        std::string ntetstr = toString(ntets);

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1)
          {
            {
              unsigned npts=5;
              //unsigned ntets=2;
              if (ntets == 1) npts=4;
              static  SingleTetFixture::Point node_coord_data[  ] = {
                { 0 , 0 , 0 } , { 1 , 0 , 0 } , { 0 , 1 , 0 } , { 0 , 0 , 1 }, {1, 1, 1} };

              // Hard coded tetra node ids for all the tetra nodes in the entire mesh
              static  SingleTetFixture::TetIds tetra_node_ids[] = {
                { 1, 2, 3, 4}, {2, 3, 4, 5} };

              stk::percept::SingleTetFixture mesh(pm, false, npts, node_coord_data, ntets, tetra_node_ids);

              stk::io::put_io_part_attribute(  mesh.m_block_tet );
              mesh.m_metaData.commit();
              mesh.populate();

              bool isCommitted = true;
              percept::PerceptMesh eMesh(&mesh.m_metaData, &mesh.m_bulkData, isCommitted);
              eMesh.save_as("tmp.e");
            }

            {
              percept::PerceptMesh eMesh;
              eMesh.open("tmp.e");

              int scalarDimension = 0; // a scalar
              stk::mesh::FieldBase* proc_rank_field = eMesh.add_field("proc_rank", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);
              stk::mesh::FieldBase* refine_field = eMesh.add_field("refine_field", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);

              const std::string part_for_edges_name = "tmp_part_for_edges";
              stk::mesh::Part * part_for_edges = 0;
              if (do_part_for_edges) {
                part_for_edges = &eMesh.get_fem_meta_data()->declare_part(part_for_edges_name, shards::CellTopology( shards::getCellTopologyData<shards::Line<2> >() ) );
              }

              if (part_for_edges) stk::io::put_io_part_attribute(*part_for_edges);
              Local_Tet4_Tet4_N break_tet_to_tet_N(eMesh);

              eMesh.commit();

              if (part_for_edges)
                {
                  stk::mesh::create_edges(*eMesh.get_bulk_data());

                  // sanity on mesh counts; overall time
                  std::vector<size_t> counts;
                  stk::mesh::comm_mesh_counts(*eMesh.get_bulk_data() , counts);

                  std::cout << "Mesh has  "
                            << counts[0] << " nodes, "
                            << counts[1] << " edges, "
                            << counts[2] << " faces, "
                            << counts[3] << " elements" << std::endl;
                  if (1)
                    {
                      stk::mesh::PartVector add_parts(1, part_for_edges), remove_parts;
                      std::vector<stk::mesh::Entity> edges;
                      const std::vector<stk::mesh::Bucket*> & buckets = eMesh.get_bulk_data()->buckets( eMesh.edge_rank() );
                      for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
                        {
                          stk::mesh::Bucket & bucket = **k ;
                          const unsigned num_edges_in_bucket = bucket.size();
                          for (unsigned iEdge = 0; iEdge < num_edges_in_bucket; iEdge++)
                            {
                              stk::mesh::Entity edge = bucket[iEdge];
                              edges.push_back(edge);
                            }
                        }

                      eMesh.get_bulk_data()->modification_begin();
                      for (unsigned iEdge = 0; iEdge < edges.size(); iEdge++)
                        {
                          eMesh.get_bulk_data()->change_entity_parts( edges[iEdge], add_parts, remove_parts );
                        }
                      eMesh.get_bulk_data()->modification_end();
                    }
                  //eMesh.print_info("after create edges", 2);
                }

              std::cout << "nele= " << eMesh.get_number_elements() << std::endl;
              std::cout << "nnode= " << eMesh.get_number_nodes() << std::endl;
              stk::mesh::Entity element_0 = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), 1);
              double *fd = eMesh.field_data(refine_field, element_0);
              fd[0] = 1.0;

              eMesh.save_as( output_files_loc+ntetstr+"_tet_with_edges_0.e");

              ElementRefinePredicate erp(eMesh, 0, refine_field, 0.0);

              PredicateBasedElementAdapter<ElementRefinePredicate>
                breaker(erp, eMesh, break_tet_to_tet_N, proc_rank_field);

              breaker.setRemoveOldElements(false);
              breaker.setAlwaysInitializeNodeRegistry(false);

              breaker.doBreak();
              eMesh.save_as(output_files_loc+ntetstr+"_tet_with_edges_1.e");


              breaker.deleteParentElements();
              eMesh.save_as(output_files_loc+ntetstr+"_tet_with_edges_2.e");
              std::cout << "nele= " << eMesh.get_number_elements() << std::endl;
              std::cout << "nnode= " << eMesh.get_number_nodes() << std::endl;

              eMesh.dump_vtk(ntetstr+"-tet.vtk",true);

              if (1)
                {
                  // sanity on mesh counts; overall time
                  std::vector<size_t> counts;
                  stk::mesh::comm_mesh_counts(*eMesh.get_bulk_data() , counts);

                  std::cout << "Mesh has  "
                            << counts[0] << " nodes, "
                            << counts[1] << " edges, "
                            << counts[2] << " faces, "
                            << counts[3] << " elements" << std::endl;
                }
            }
          }
      }

      STKUNIT_UNIT_TEST(regr_localRefiner, tet_with_edges)
      {
        //do_tet_edge_test(1, true);
        do_tet_edge_test(2, false);

      }

      //  ====================================================================================================================================
      //  ====================================================================================================================================
      //  Local quad hanging-node refinement
      //  ====================================================================================================================================
      //  ====================================================================================================================================

      class SetElementFieldQuadCorner : public percept::ElementOp
      {
        percept::PerceptMesh& m_eMesh;
        stk::mesh::FieldBase *m_field;
      public:
        SetElementFieldQuadCorner(percept::PerceptMesh& eMesh, stk::mesh::FieldBase *field) : m_eMesh(eMesh), m_field(field) {
        }

        virtual bool operator()(const stk::mesh::Entity element, stk::mesh::FieldBase *field,  const mesh::BulkData& bulkData)
        {
          double *f_data = m_eMesh.field_data(field, element);
          stk::mesh::FieldBase* coord_field = m_eMesh.get_coordinates_field();

          const MyPairIterRelation elem_nodes(bulkData, element, stk::mesh::MetaData::NODE_RANK );

          unsigned num_node = elem_nodes.size();
          double c[2] = {0,0};
          for (unsigned inode=0; inode < num_node; inode++)
            {
              mesh::Entity node = elem_nodes[ inode ].entity();
              double *c_data = m_eMesh.field_data(coord_field, node);
              c[0] += c_data[0]/double(num_node);
              c[1] += c_data[1]/double(num_node);
            }
          if ((-0.1 < c[0] && c[0] < 0.1) &&
              (-0.1 < c[1] && c[1] < 0.1) )
            {
              f_data[0] = 1.0;
            }
          else
            {
              f_data[0] = -1.0;
            }
          // FIXME tmp
          //f_data[0] = 1.0;
          return false;  // don't terminate the loop
        }
        virtual void init_elementOp()
        {
          std::vector< const stk::mesh::FieldBase *> fields(1,m_field);
          //stk::mesh::copy_owned_to_shared( *m_eMesh.get_bulk_data(), fields);
          stk::mesh::communicate_field_data(m_eMesh.get_bulk_data()->shared_aura(), fields);
        }
        virtual void fini_elementOp() {
          std::vector< const stk::mesh::FieldBase *> fields(1,m_field);
          //stk::mesh::copy_owned_to_shared( *m_eMesh.get_bulk_data(), fields);
          stk::mesh::communicate_field_data(m_eMesh.get_bulk_data()->shared_aura(), fields);
        }

      };


      static void do_quad_local_corner_refine_sidesets(PerceptMesh& eMesh, int num_time_steps, bool save_intermediate=false, bool delete_parents=false)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size==1 || p_size == 3)
          {
            Local_Quad4_Quad4_N break_quad_to_quad_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field    = eMesh.add_field("proc_rank", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);
            stk::mesh::FieldBase* refine_field       = eMesh.add_field("refine_field", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);

            // for plotting, use doubles, for internal use, use int
            stk::mesh::FieldBase* refine_level_d     = eMesh.add_field("refine_level_d", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);
            ScalarIntFieldType& refine_level       = eMesh.get_fem_meta_data()->declare_field<ScalarIntFieldType>("refine_level");
            stk::mesh::put_field( refine_level , stk::mesh::MetaData::ELEMENT_RANK , eMesh.get_fem_meta_data()->universal_part());
            stk::io::set_field_role(refine_level, Ioss::Field::TRANSIENT);

            eMesh.commit();
            if (1)
              {
                SetElementRefineFieldValue set_ref_field(eMesh);
                eMesh.elementOpLoop(set_ref_field, refine_level_d);
              }

            std::cout << "qual_local initial number elements= " << eMesh.get_number_elements() << std::endl;

            eMesh.save_as( output_files_loc+"quad_local_"+post_fix[p_size]+".e.0");

            stk::mesh::Selector univ_selector(eMesh.get_fem_meta_data()->universal_part());

            ElementRefinePredicate erp(eMesh, &univ_selector, refine_field, 0.0);
            PredicateBasedElementAdapter<ElementRefinePredicate>
              breaker(erp, eMesh, break_quad_to_quad_N, proc_rank_field);

            // special for local adaptivity
            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);
            breaker.setNeedsRemesh(false); // special for quad/hex hanging node

            SetElementFieldQuadCorner set_ref_field(eMesh, refine_field);

            int num_ref_passes = 6;
            int num_unref_passes = 4;
            int iplot=0;
            if (1)
              {
                char buf[1000];
                sprintf(buf, "%04d", iplot);
                if (iplot == 0)
                  eMesh.save_as("quad_square_anim.e");
                else
                  eMesh.save_as("quad_square_anim.e-s"+std::string(buf));
                ++iplot;
              }
            for (int ipass=0; ipass < num_ref_passes; ipass++)
              {

                eMesh.elementOpLoop(set_ref_field, refine_field);
                eMesh.save_as(output_files_loc+"quad_anim_set_field_"+post_fix[p_size]+".e.s-"+toString(ipass+1));

                bool enforce_what[3] = {false, false, true};
                erp.refine(breaker, enforce_what);
                //breaker.doBreak();

                MPI_Barrier( MPI_COMM_WORLD );
                //VERIFY_OP_ON(true,==,false,"here");
                bool check_what[3] = {false, false, true};
                bool is_valid_2_to_1 = erp.check_two_to_one(check_what);
                bool check_what_1[3] = {true, false, false};
                bool is_valid_2_to_1_1 = erp.check_two_to_one(check_what_1);
                std::cout << "P[" << eMesh.get_rank() << "] done... ipass= " << ipass << " quad_local number elements= "
                          << eMesh.get_number_elements() << " check_two_to_one= " << is_valid_2_to_1
                          << " node check_two_to_one= " << is_valid_2_to_1_1 << std::endl;

                //breaker.deleteParentElements();
                //eMesh.save_as("square_anim."+toString(ipass+1)+".e");
                if (1)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", iplot);
                    if (iplot == 0)
                      eMesh.save_as("quad_square_anim.e");
                    else
                      eMesh.save_as("quad_square_anim.e-s"+std::string(buf));
                    ++iplot;
                  }

              }

            for (int iunref_pass=0; iunref_pass < num_unref_passes; iunref_pass++)
              {
                eMesh.elementOpLoop(set_ref_field, refine_field);

                bool enforce_what[3] = {false, false, true};
                erp.unrefine(breaker, enforce_what);
                std::cout << "P[" << eMesh.get_rank() << "] done... iunref_pass= " << iunref_pass << " quad_local number elements= " << eMesh.get_number_elements() << std::endl;
                if (1)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", iplot);
                    if (iplot == 0)
                      eMesh.save_as("quad_square_anim.e");
                    else
                      eMesh.save_as("quad_square_anim.e-s"+std::string(buf));
                    ++iplot;
                  }
              }

            SetElementRefineFieldValue set_ref_field_val_unref_all(eMesh, -1.0);

            eMesh.save_as(output_files_loc+"quad_tmp_square_sidesets_quad_local_unref_"+post_fix[p_size]+".e");

            eMesh.save_as(output_files_loc+"quad_square_sidesets_final_quad_local_"+post_fix[p_size]+".e.s-"+toString(num_time_steps) );

            for (int iunref=0; iunref < 10; iunref++)
              {
                eMesh.elementOpLoop(set_ref_field_val_unref_all, refine_field);
                std::cout << "P[" << eMesh.get_rank() << "] iunrefAll_pass= " << iunref <<  std::endl;
                bool enforce_what[3] = {false, false, true};
                erp.unrefine(breaker, enforce_what);

                if (1)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", iplot);
                    if (iplot == 0)
                      eMesh.save_as("quad_square_anim.e");
                    else
                      eMesh.save_as("quad_square_anim.e-s"+std::string(buf));
                    ++iplot;
                  }

                std::cout << "P[" << eMesh.get_rank() << "] done... iunrefAll_pass= " << iunref << " quad_local number elements= " << eMesh.get_number_elements() << std::endl;
              }

            if (delete_parents)
              breaker.deleteParentElements();
            std::cout << "quad_local final number elements= " << eMesh.get_number_elements() << std::endl;
            eMesh.save_as(output_files_loc+"quad_square_sidesets_final_unrefed_quad_local_"+post_fix[p_size]+".e."+toString(num_time_steps) );
            //exit(123);

            // end_demo
          }
      }

      STKUNIT_UNIT_TEST(regr_localRefiner, break_quad_to_quad_N_5_ElementBased_quad_local_square_sidesets)
      {
        bool do_test = true;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        if (do_test) {

          // if this is a fresh installation, set to true to generate the initial meshes needed for this test (once only)
          bool do_bootstrap_mesh = true;
          if (do_bootstrap_mesh)
          {
            const unsigned n = 5;
            const unsigned nx = n , ny = n;

            bool createEdgeSets = true;
            percept::QuadFixture<double> fixture( pm , nx , ny, createEdgeSets);
            fixture.set_bounding_box(-1,1, -1, 1);

            bool isCommitted = false;
            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, isCommitted);

            eMesh.commit();

            fixture.generate_mesh();
            eMesh.save_as(input_files_loc+"quad_square_quad4_0.e");
            //eMesh.print_info("test1",2);
            //eMesh.dump_vtk("sqquad3.vtk",false);
          }
          //if (1) return;

          {
            PerceptMesh eMesh;
            eMesh.open(input_files_loc+"quad_square_quad4_0.e");

            eMesh.output_active_children_only(true);

            //int num_time_steps = 10;  // 10 for stress testing
            //for (int istep=1; istep <= num_time_steps; istep++)
            do_quad_local_corner_refine_sidesets(eMesh, 10, false, true);
          }
        }
      }

      //  ====================================================================================================================================
      //  ====================================================================================================================================
      //  Local hex hanging-node refinement
      //  ====================================================================================================================================
      //  ====================================================================================================================================

      class SetElementFieldHexCorner : public percept::ElementOp
      {
        percept::PerceptMesh& m_eMesh;
        stk::mesh::FieldBase *m_field;
      public:
        SetElementFieldHexCorner(percept::PerceptMesh& eMesh, stk::mesh::FieldBase *field) : m_eMesh(eMesh), m_field(field) {
        }

        virtual bool operator()(const stk::mesh::Entity element, stk::mesh::FieldBase *field,  const mesh::BulkData& bulkData)
        {
          double *f_data = m_eMesh.field_data(field, element);
          stk::mesh::FieldBase* coord_field = m_eMesh.get_coordinates_field();

          const MyPairIterRelation elem_nodes(bulkData, element, stk::mesh::MetaData::NODE_RANK );

          unsigned num_node = elem_nodes.size();
          double c[] = {0,0,0};
          for (unsigned inode=0; inode < num_node; inode++)
            {
              mesh::Entity node = elem_nodes[ inode ].entity();
              double *c_data = m_eMesh.field_data(coord_field, node);
              c[0] += c_data[0]/double(num_node);
              c[1] += c_data[1]/double(num_node);
              c[2] += c_data[2]/double(num_node);
            }
          if ((-0.1 < c[0] && c[0] < 0.1) &&
              (-0.1 < c[1] && c[1] < 0.1) &&
              (-0.1 < c[2] && c[2] < 0.1) )
            {
              f_data[0] = 1.0;
            }
          else
            {
              f_data[0] = -1.0;
            }
          // FIXME tmp
          //f_data[0] = 1.0;
          return false;  // don't terminate the loop
        }
        virtual void init_elementOp()
        {
          std::vector< const stk::mesh::FieldBase *> fields(1,m_field);
          //stk::mesh::copy_owned_to_shared( *m_eMesh.get_bulk_data(), fields);
          stk::mesh::communicate_field_data(m_eMesh.get_bulk_data()->shared_aura(), fields);
        }
        virtual void fini_elementOp() {
          std::vector< const stk::mesh::FieldBase *> fields(1,m_field);
          //stk::mesh::copy_owned_to_shared( *m_eMesh.get_bulk_data(), fields);
          stk::mesh::communicate_field_data(m_eMesh.get_bulk_data()->shared_aura(), fields);
        }

      };


      static void do_hex_local_corner_refine_sidesets(PerceptMesh& eMesh, int num_time_steps, bool save_intermediate=false, bool delete_parents=false)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size==1 || p_size == 3)
          {
            Local_Hex8_Hex8_N break_hex_to_hex_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field    = eMesh.add_field("proc_rank", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);
            stk::mesh::FieldBase* refine_field       = eMesh.add_field("refine_field", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);

            // for plotting, use doubles, for internal use, use int
            stk::mesh::FieldBase* refine_level_d     = eMesh.add_field("refine_level_d", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);
            ScalarIntFieldType& refine_level       = eMesh.get_fem_meta_data()->declare_field<ScalarIntFieldType>("refine_level");
            stk::mesh::put_field( refine_level , stk::mesh::MetaData::ELEMENT_RANK , eMesh.get_fem_meta_data()->universal_part());
            stk::io::set_field_role(refine_level, Ioss::Field::TRANSIENT);

            eMesh.commit();
            if (1)
              {
                SetElementRefineFieldValue sr(eMesh);
                eMesh.elementOpLoop(sr, refine_level_d);
              }

            std::cout << "hex_local initial number elements= " << eMesh.get_number_elements() << std::endl;

            eMesh.save_as( output_files_loc+"hex_local_"+post_fix[p_size]+".e.0");

            stk::mesh::Selector univ_selector(eMesh.get_fem_meta_data()->universal_part());

            ElementRefinePredicate erp(eMesh, &univ_selector, refine_field, 0.0);
            PredicateBasedElementAdapter<ElementRefinePredicate>
              breaker(erp, eMesh, break_hex_to_hex_N, proc_rank_field);

            // special for local adaptivity
            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);
            breaker.setNeedsRemesh(false); // special for quad/hex hanging node

            SetElementFieldHexCorner set_ref_field(eMesh, refine_field);

            int num_ref_passes = 6;
            int num_unref_passes = 4;
            int iplot=0;
            if (1)
              {
                char buf[1000];
                sprintf(buf, "%04d", iplot);
                if (iplot == 0)
                  eMesh.save_as("hex_cube_anim.e");
                else
                  eMesh.save_as("hex_cube_anim.e-s"+std::string(buf));
                ++iplot;
              }
            for (int ipass=0; ipass < num_ref_passes; ipass++)
              {

                eMesh.elementOpLoop(set_ref_field, refine_field);
                //eMesh.save_as(output_files_loc+"hex_anim_set_field_"+post_fix[p_size]+".e.s-"+toString(ipass+1));
                if (1)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", ipass);
                    if (ipass == 0)
                      eMesh.save_as("hex_set_field.e");
                    else
                      eMesh.save_as("hex_set_field.e-s"+std::string(buf));
                  }

                // {node-, edge-, face-neighors}
                bool enforce_what[3] = {false, true, true};
                erp.refine(breaker, enforce_what);
                //breaker.doBreak();

                MPI_Barrier( MPI_COMM_WORLD );
                //VERIFY_OP_ON(true,==,false,"here");
                bool check_what[3] = {false, false, true};
                bool is_valid_2_to_1 = erp.check_two_to_one(check_what);
                bool check_what_1[3] = {true, false, false};
                bool is_valid_2_to_1_1 = erp.check_two_to_one(check_what_1);
                std::cout << "P[" << eMesh.get_rank() << "] done... ipass= " << ipass << " hex_local number elements= "
                          << eMesh.get_number_elements() << " check_two_to_one= " << is_valid_2_to_1
                          << " node check_two_to_one= " << is_valid_2_to_1_1 << std::endl;
                STKUNIT_EXPECT_TRUE(is_valid_2_to_1);
                STKUNIT_EXPECT_TRUE(is_valid_2_to_1_1);
                //breaker.deleteParentElements();
                //eMesh.save_as("square_anim."+toString(ipass+1)+".e");
                if (1)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", iplot);
                    if (iplot == 0)
                      eMesh.save_as("hex_cube_anim.e");
                    else
                      eMesh.save_as("hex_cube_anim.e-s"+std::string(buf));
                    ++iplot;
                  }

              }

            for (int iunref_pass=0; iunref_pass < num_unref_passes; iunref_pass++)
              {
                eMesh.elementOpLoop(set_ref_field, refine_field);

                bool enforce_what[3] = {false, false, true};
                erp.unrefine(breaker, enforce_what);
                std::cout << "P[" << eMesh.get_rank() << "] done... iunref_pass= " << iunref_pass << " hex_local number elements= " << eMesh.get_number_elements() << std::endl;
                if (1)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", iplot);
                    if (iplot == 0)
                      eMesh.save_as("hex_cube_anim.e");
                    else
                      eMesh.save_as("hex_cube_anim.e-s"+std::string(buf));
                    ++iplot;
                  }
              }

            SetElementRefineFieldValue set_ref_field_val_unref_all(eMesh, -1.0);

            eMesh.save_as(output_files_loc+"hex_tmp_square_sidesets_hex_local_unref_"+post_fix[p_size]+".e");

            eMesh.save_as(output_files_loc+"hex_cube_sidesets_final_hex_local_"+post_fix[p_size]+".e.s-"+toString(num_time_steps) );

            for (int iunref=0; iunref < 10; iunref++)
              {
                eMesh.elementOpLoop(set_ref_field_val_unref_all, refine_field);
                std::cout << "P[" << eMesh.get_rank() << "] iunrefAll_pass= " << iunref <<  std::endl;
                bool enforce_what[3] = {false, false, true};
                erp.unrefine(breaker, enforce_what);

                if (1)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", iplot);
                    if (iplot == 0)
                      eMesh.save_as("hex_cube_anim.e");
                    else
                      eMesh.save_as("hex_cube_anim.e-s"+std::string(buf));
                    ++iplot;
                  }

                std::cout << "P[" << eMesh.get_rank() << "] done... iunrefAll_pass= " << iunref << " hex_local number elements= " << eMesh.get_number_elements() << std::endl;
              }

            if (delete_parents)
              breaker.deleteParentElements();
            std::cout << "hex_local final number elements= " << eMesh.get_number_elements() << std::endl;
            eMesh.save_as(output_files_loc+"hex_cube_sidesets_final_unrefed_hex_local_"+post_fix[p_size]+".e."+toString(num_time_steps) );
            //exit(123);

            // end_demo
          }
      }

      STKUNIT_UNIT_TEST(regr_localRefiner, break_hex_to_hex_N_5_ElementBased_hex_local_square_sidesets)
      {
        bool do_test = true;
        //stk::ParallelMachine pm = MPI_COMM_WORLD ;

        if (do_test) {

          // if this is a fresh installation, set to true to generate the initial meshes needed for this test (once only)
          bool do_bootstrap_mesh = true;
          if (do_bootstrap_mesh)
          {
            const unsigned n = 5;
            const unsigned nx = n , ny = n, nz = n;

            percept::PerceptMesh eMesh(3u);
            std::string gmesh_spec = toString(nx)+"x"+toString(ny)+"x"+toString(nz)+"|bbox:-1,-1,-1,1,1,1";
            eMesh.new_mesh(percept::GMeshSpec(gmesh_spec));
            eMesh.commit();

            eMesh.save_as(input_files_loc+"hex_cube_hex4_0.e");
            //eMesh.print_info("test1",2);
            //eMesh.dump_vtk("sqhex3.vtk",false);
          }
          //if (1) return;

          {
            PerceptMesh eMesh;
            eMesh.open(input_files_loc+"hex_cube_hex4_0.e");

            eMesh.output_active_children_only(true);

            //int num_time_steps = 10;  // 10 for stress testing
            //for (int istep=1; istep <= num_time_steps; istep++)
            do_hex_local_corner_refine_sidesets(eMesh, 10, false, true);
          }
        }
      }

      //=============================================================================
      // Encore/Percept mimic
      //=============================================================================
      //=============================================================================

      class Mimic_Encr_Entity_Marker {
      public:
        Mimic_Encr_Entity_Marker(stk::percept::PerceptMesh & pMesh, const stk::mesh::FieldBase * marker_field, const stk::mesh::EntityRank entity_rank)
        : my_pMesh(pMesh), my_marker_field(marker_field), my_entity_rank(entity_rank) {}
        void update_markers();
        Int get_marker(stk::mesh::Entity entity);

      protected:
        stk::percept::PerceptMesh & my_pMesh;
        const stk::mesh::FieldBase * my_marker_field;
        const stk::mesh::EntityRank my_entity_rank;
        std::vector<stk::mesh::Entity> my_entities;
        std::vector<Int> my_markers;
      };

      class Mimic_Encr_Percept_Edge_Adapter : public stk::adapt::IEdgeAdapter
      {
      public:
        Mimic_Encr_Percept_Edge_Adapter(Mimic_Encr_Entity_Marker & element_marker, stk::percept::PerceptMesh & pMesh, stk::adapt::UniformRefinerPatternBase & bp, stk::mesh::FieldBase * proc_rank_field=0)
        : stk::adapt::IEdgeAdapter(pMesh, bp, proc_rank_field), my_element_marker(element_marker) {}
        virtual int markEdge(const stk::mesh::Entity element, unsigned which_edge, stk::mesh::Entity node0, stk::mesh::Entity node1,
          double *coord0, double *coord1, std::vector<int>* existing_edge_marks);
      protected:
        Mimic_Encr_Entity_Marker & my_element_marker;
      };

      int
      Mimic_Encr_Percept_Edge_Adapter::markEdge(
          const stk::mesh::Entity element,
          unsigned which_edge,
          stk::mesh::Entity node0,
          stk::mesh::Entity node1,
          double *coord0,
          double *coord1,
          std::vector<int>* existing_edge_marks)
      {
        std::vector<stk::mesh::Entity> nodes;
        nodes.push_back(node0);
        nodes.push_back(node1);
        std::vector<stk::mesh::Entity> edge_elems;
        get_entities_through_relations(*m_eMesh.get_bulk_data(), nodes, stk::mesh::MetaData::ELEMENT_RANK, edge_elems);
        ThrowRequire(!edge_elems.empty());

        bool all_marked_for_refinement = true;
        bool all_marked_for_unrefinement = true;
        for ( UInt ie = 0 ; ie < edge_elems.size(); ++ie )
        {
          stk::mesh::Entity edge_elem = edge_elems[ie];
          ThrowRequire(m_eMesh.is_valid(edge_elem));

          Int element_marker = my_element_marker.get_marker(edge_elem);
          if (element_marker <= 0) all_marked_for_refinement = false;
          if (element_marker >= 0) all_marked_for_unrefinement = false;
        }
        ThrowRequire(!(all_marked_for_refinement && all_marked_for_unrefinement));

        int mark=stk::adapt::DO_NOTHING;
        if (all_marked_for_refinement)
        {
          mark=stk::adapt::DO_REFINE;
        }
        else if (all_marked_for_unrefinement)
        {
          mark= stk::adapt::DO_UNREFINE;
        }
        return mark;
      }

      void
      Mimic_Encr_Entity_Marker::update_markers()
      {
        my_entities.clear();
        stk::mesh::get_entities( *my_pMesh.get_bulk_data(), my_entity_rank, my_entities );
        my_markers.resize(my_entities.size());

        for (unsigned i=0; i<my_entities.size(); ++i)
        {
          mesh::Entity elem = my_entities[i];

          // if (my_pMesh.isParentElement(elem,false))
          // {
          //   my_markers[i] = stk::adapt::DO_NOTHING;
          // }
          // else
          {
            my_markers[i] = *((Int *)my_pMesh.field_data( *my_marker_field, my_entities[i] ));
          }
          //std::cout << "Storing element marker for element " << m_eMesh.identifier(elem) << " = " << my_markers[i] << std::endl;
        }
      }

      Int
      Mimic_Encr_Entity_Marker::get_marker(stk::mesh::Entity entity)
      {
        std::vector<stk::mesh::Entity>::iterator it = std::find(my_entities.begin(), my_entities.end(), entity);
        if (it == my_entities.end())
        {
          return stk::adapt::DO_NOTHING;
        }
        else
        {
          ThrowRequire(my_pMesh.is_valid(*it));
          const UInt index = std::distance(my_entities.begin(), it);
          return my_markers[index];
        }
      }

      static void mimic_encr_function_and_element_marker(stk::percept::PerceptMesh & pMesh, double time, ScalarFieldType & function_field, mesh::Field<int> & marker_field)
      {
        std::vector<stk::mesh::Entity> entities;
        stk::mesh::get_entities( *pMesh.get_bulk_data(), stk::mesh::MetaData::NODE_RANK, entities );

        VectorFieldType* coordField = pMesh.get_coordinates_field();

        for (unsigned i=0; i<entities.size(); ++i)
        {
          mesh::Entity node = entities[i];
          double *coord_data = pMesh.field_data(coordField, node);
          double *function_data = pMesh.field_data(&function_field, node);
          *function_data = (2.0*time-10.0-coord_data[0]-coord_data[1])*(2.0*time-10.0-coord_data[0]-coord_data[1]);
        }

        // Now calculate element field marker

        const double refine_lower = 0.0;
        const double refine_upper = 1.0;
        double coarsen_lower = 4.0;
        double coarsen_upper = 1000.0;

        entities.clear();
        stk::mesh::get_entities( *pMesh.get_bulk_data(), stk::mesh::MetaData::ELEMENT_RANK, entities );

        for (unsigned i=0; i<entities.size(); ++i)
        {
          mesh::Entity elem = entities[i];

          //if (pMesh.isParentElement(elem,false)) continue;

          int elem_mark = 0;
          const percept::MyPairIterRelation elem_nodes (pMesh, elem, stk::mesh::MetaData::NODE_RANK );
          unsigned num_node = elem_nodes.size();
          for (unsigned inode=0; inode < num_node; inode++)
          {
            mesh::Entity node = elem_nodes[ inode ].entity();
            double cur_value = *(pMesh.field_data(&function_field, node));

            //If one of the values is in the refine region, immediately mark and return
            if(refine_upper>=cur_value && cur_value>=refine_lower)
            {
              elem_mark = 1;
              break;
            }
            else if(coarsen_lower<=cur_value && cur_value<=coarsen_upper)
            {
              elem_mark = -1;
            }
          }
          Int * marker_data = pMesh.get_bulk_data()->field_data( marker_field, elem );
          *marker_data = elem_mark;
        }

        // update ghost element values
        {
          std::vector< const stk::mesh::FieldBase *> fields(1, &marker_field);
          //stk::mesh::copy_owned_to_shared( *m_eMesh.get_bulk_data(), fields);
          stk::mesh::communicate_field_data(pMesh.get_bulk_data()->shared_aura(), fields);
        }
      }

      STKUNIT_UNIT_TEST(regr_localRefiner, break_tri_to_tri_N_EdgeFromElementMarker_MimicEncr)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
        {
          PerceptMesh eMesh;
          eMesh.open(input_files_loc+"square_tri3_uns.e");

          Local_Tri3_Tri3_N break_tri_to_tri_N(eMesh);
          int scalarDimension = 0; // a scalar
          stk::mesh::FieldBase* proc_rank_field = eMesh.add_field("proc_rank", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);
          ScalarFieldType * function_field =
            (ScalarFieldType *) eMesh.add_field("s_node", stk::mesh::MetaData::NODE_RANK, scalarDimension);
          mesh::Field<int> & marker_field =  eMesh.get_fem_meta_data()->declare_field< mesh::Field<int> >("marker_field_1");
          stk::io::set_field_role(marker_field, Ioss::Field::TRANSIENT);
          stk::mesh::put_field( marker_field, stk::mesh::MetaData::ELEMENT_RANK, eMesh.get_fem_meta_data()->universal_part() );
          eMesh.commit();
          eMesh.output_active_children_only(true);

          eMesh.save_as(output_files_loc+"mimic_encr_percept_square_adapt_tri_"+post_fix[p_size]+".e");

          Mimic_Encr_Entity_Marker element_marker(eMesh, &marker_field, stk::mesh::MetaData::ELEMENT_RANK);
          Mimic_Encr_Percept_Edge_Adapter breaker(element_marker, eMesh, break_tri_to_tri_N, proc_rank_field);
          breaker.setRemoveOldElements(false);
          breaker.setAlwaysInitializeNodeRegistry(false);

          for (int ipass=0; ipass < 10; ipass++)
          {
            const double time = ipass+1.0;
            mimic_encr_function_and_element_marker(eMesh, time, *function_field, marker_field);
            element_marker.update_markers();
            breaker.doBreak();

            for (int iunref_pass=0; iunref_pass < 1; iunref_pass++)
            {
              std::cout << "P[" << eMesh.get_rank() << "] iunref_pass= " << iunref_pass << std::endl;
              //mimic_encr_function_and_element_marker(eMesh, time, *function_field, marker_field);
              //element_marker.update_markers();
              ElementUnrefineCollection elements_to_unref = breaker.buildUnrefineList();
              std::cout << "tmp srk ipass= " << ipass << " elements_to_unref.size() = " << elements_to_unref.size() << std::endl;
              breaker.unrefineTheseElements(elements_to_unref);
            }

            std::stringstream fileid_ss;
            fileid_ss << std::setfill('0') << std::setw(4) << ipass+2;
            eMesh.save_as(output_files_loc+"mimic_encr_percept_square_adapt_tri_"+post_fix[p_size]+".e-s"+fileid_ss.str(), time);
          }
        }
      }
    }
  }
}
