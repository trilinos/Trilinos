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
#include <stk_adapt/ElementRefinePredicate.hpp>
#include <stk_adapt/PredicateBasedElementAdapter.hpp>
#include <stk_adapt/PredicateBasedEdgeAdapter.hpp>

#include <stk_percept/function/ElementOp.hpp>

#include <regression_tests/RegressionTestLocalRefiner.hpp>

namespace stk
{
  namespace adapt
  {
    namespace heavy_tests
    {
      using namespace regression_tests;

#if 1
      static const std::string path_sep = "._.";
      static const std::string input_files_loc="./input_files"+path_sep;
      static const std::string output_files_loc="./output_files"+path_sep;
#else
      static const std::string input_files_loc="./input_files/";
      static const std::string output_files_loc="./output_files/";
#endif

      static std::string post_fix[9] = {"np0", "np1", "np2", "np3", "np4", "np5", "np6", "np7", "np8"};

      //=============================================================================
      //=============================================================================
      //=============================================================================


      static void do_moving_shock_test_large_test(int num_time_steps, bool save_intermediate=false, bool delete_parents=false)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //         shock_width = 1./50.; //1./25.0;
        //         shock_diff_criterion = 0.1;
        shock_width = 1./500.; //1./25.0;
        double shock_diff_criterion = 0.1;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size==1 || p_size == 2 || p_size == 8)
          {
            percept::PerceptMesh eMesh;
            eMesh.open(input_files_loc+"OUO_large_testh1tet.g");

            Local_Tet4_Tet4_N break_tet_to_tet_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field    = eMesh.add_field("proc_rank", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);
            stk::mesh::FieldBase* refine_field       = eMesh.add_field("refine_field", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);
            stk::mesh::FieldBase* nodal_refine_field = eMesh.add_field("nodal_refine_field", eMesh.node_rank(), scalarDimension);
            eMesh.commit();

            //eMesh.delete_side_sets();

            std::cout << "moving_shock initial number elements= " << eMesh.get_number_elements() << std::endl;

            eMesh.save_as( output_files_loc+"OUO_large_scale_moving_shock_"+post_fix[p_size]+".e.0");

            stk::mesh::Selector univ_selector(eMesh.get_fem_meta_data()->universal_part());

            PlaneShock shock;
            shock.plane_point_init[0] = 0.0;
            shock.plane_point_init[1] = -0.07;
            shock.plane_point_init[2] = -0.35;
            shock.plane_normal[0] = 0;
            shock.plane_normal[1] = .25/std::sqrt(.25*.25+1.);
            shock.plane_normal[2] = 1./std::sqrt(.25*.25+1.);

            ShockBasedRefinePredicate srp(nodal_refine_field, eMesh, &univ_selector, refine_field, 0.0, shock, 0.0, shock_diff_criterion);

            PredicateBasedEdgeAdapter<ShockBasedRefinePredicate>
              breaker(srp,
                      eMesh, break_tet_to_tet_N, proc_rank_field);

            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);
            //breaker.setIgnoreSideSets(true);

            // x,y,z: [-.08,.08],[-.08,.08],[-.5,-.2]
            double delta_shock_displacement = 0.2*0.08*(10./((double)num_time_steps));
            double shock_displacement = 0.0;
            int num_ref_passes = 2;
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

                eMesh.save_as(output_files_loc+"OUO_large_scale_moving_shock_"+post_fix[p_size]+"_ref.e");
                //exit(123);

                breaker.getNodeRegistry().init_entity_repo();
                for (int iunref_pass=0; iunref_pass < num_unref_passes; iunref_pass++)
                  {
                    ElementUnrefineCollection elements_to_unref = breaker.buildUnrefineList();
                    std::cout << "P[" << eMesh.get_rank() << "] iunref_pass= " << iunref_pass <<  std::endl;
                    breaker.unrefineTheseElements(elements_to_unref);
                    std::cout << "P[" << eMesh.get_rank() << "] done... iunref_pass= " << iunref_pass << " moving_shock number elements= " << eMesh.get_number_elements() << std::endl;
                  }

                if (delete_parents && istep == num_time_steps-1)
                  {
                    breaker.deleteParentElements();
                  }
                if (istep == num_time_steps-1 || save_intermediate)
                  eMesh.save_as(output_files_loc+"OUO_large_scale_moving_shock_"+post_fix[p_size]+".e."+toString(istep+1) );

                shock_displacement += delta_shock_displacement;

              }

            eMesh.save_as(output_files_loc+"OUO_large_scale_final_moving_shock_"+post_fix[p_size]+"."+toString(num_time_steps)+".e" );
            eMesh.save_as(output_files_loc+"OUO_large_scale_final_moving_shock_"+post_fix[p_size]+".e" );
            for (int iunref=0; iunref < 10; iunref++)
              {
                breaker.unrefineAll();
              }
            breaker.deleteParentElements();
            std::cout << "moving_shock final number elements= " << eMesh.get_number_elements() << std::endl;
            eMesh.save_as(output_files_loc+"OUO_large_scale_final_unrefed_moving_shock_"+post_fix[p_size]+".e."+toString(num_time_steps) );

            // end_demo
          }
      }

      STKUNIT_UNIT_TEST(heavy_localRefiner, break_tet_to_tet_N_5_EdgeBased_moving_shock_large_test)
      {
        const bool do_full_demo = false;
        if (do_full_demo)
          {
            int num_time_steps = 10;  // 10 for stress testing
            for (int istep=1; istep <= num_time_steps; istep++)
              do_moving_shock_test_large_test(istep, false, true);
          }
        else
          // normal regression testing
          {
            int num_time_steps = 10;  // 30 for making animation
            do_moving_shock_test_large_test(num_time_steps, true);
            //int num_time_steps = 3;  // 10 for stress testing
            //do_moving_shock_test_large_test(num_time_steps);
          }
      }

    }
  }
}
