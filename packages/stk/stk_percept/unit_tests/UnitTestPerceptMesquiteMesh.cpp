/*--------------------------------------------------------------------*/
/*    Copyright 2009, 2011 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/Util.hpp>
#include <stk_percept/ExceptionWatch.hpp>
#include <stk_percept/fixtures/Fixture.hpp>
#include <stk_percept/fixtures/Fixture.hpp>
#include <stk_percept/fixtures/QuadFixture.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/util/PrintTable.hpp>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <Teuchos_ScalarTraits.hpp>

#include <stk_percept/fixtures/Fixture.hpp>
#include <stk_percept/fixtures/BeamFixture.hpp>
#include <stk_percept/fixtures/HeterogeneousFixture.hpp>
#include <stk_percept/fixtures/QuadFixture.hpp>
#include <stk_percept/fixtures/WedgeFixture.hpp>

#ifdef STK_BUILT_IN_SIERRA
#define STK_ADAPT_HAS_GEOMETRY
#else
#undef STK_ADAPT_HAS_GEOMETRY
#endif

#if defined( STK_ADAPT_HAS_GEOMETRY )
//#include <stk_percept/mesh/geometry/kernel/GeometryKernelOpenNURBS.hpp>
//#include <stk_percept/mesh/geometry/kernel/MeshGeometry.hpp>
//#include <stk_percept/mesh/geometry/kernel/GeometryFactory.hpp>

// place Mesquite-related headers between the StackTrace redefinitions
#define StackTraceTmp StackTrace
#undef StackTrace
#include <stk_percept/mesh/mod/mesquite-interface/PerceptMesquiteMesh.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PerceptMesquiteMeshDomain.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PMMLaplaceSmoother.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PMMLaplaceSmoother1.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PMMShapeImprover.hpp>
#include <MsqDebug.hpp>
#define StackTrace StackTraceTmp

#endif

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
      //static int printInfoLevel = 0;

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

      STKUNIT_UNIT_TEST(unit_perceptMesquite, quad_1)
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        MPI_Barrier( MPI_COMM_WORLD );

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        //if (p_size == 1 || p_size == 3)
        if (p_size <= 1)
          {
            //const unsigned p_rank = stk::parallel_machine_rank( pm );
            //const unsigned p_size = stk::parallel_machine_size( pm );

            const unsigned n = 2;
            //const unsigned nx = n , ny = n , nz = p_size*n ;
            const unsigned nx = n , ny = n;

            bool sidesets=true;
            percept::QuadFixture<double> fixture( pm , nx , ny, sidesets);
            fixture.meta_data.commit();
            fixture.generate_mesh();

            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data);
            eMesh.printInfo("quad fixture",  2);
            eMesh.saveAs(input_files_loc+"quad_smooth.0.e");

            unsigned center_node_id = 5;
            stk::mesh::Entity* node = eMesh.getBulkData()->get_entity(0, center_node_id);
            double * data = stk::mesh::field_data( *eMesh.getCoordinatesField() , *node );
            std::cout << "tmp srk  center node= " << data[0] << " " << data[1] << std::endl;
            data[0] += .2;
            data[1] += .3;

            eMesh.saveAs(output_files_loc+"quad_smooth.0_perturbed.e");

            stk::mesh::Selector boundarySelector_1(*eMesh.getNonConstPart("surface_1") );
            stk::mesh::Selector boundarySelector_2(*eMesh.getNonConstPart("surface_2") );
            stk::mesh::Selector boundarySelector_3(*eMesh.getNonConstPart("surface_3") );
            stk::mesh::Selector boundarySelector_4(*eMesh.getNonConstPart("surface_4") );
            stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4;

            PerceptMesquiteMesh pmm(&eMesh, 0, &boundarySelector);
            PerceptMesquiteMeshDomain pmd(&eMesh, 0);
            percept::PMMLaplaceSmoother ls;
            ls.run(pmm, pmd);

            STKUNIT_EXPECT_DOUBLE_EQ_APPROX(data[0], 1.0);
            STKUNIT_EXPECT_DOUBLE_EQ_APPROX(data[1], 1.0);

            eMesh.saveAs(output_files_loc+"quad_smooth.1.e");
          }
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      STKUNIT_UNIT_TEST(unit_perceptMesquite, quad_2)
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        MPI_Barrier( MPI_COMM_WORLD );

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        //if (p_size == 1 || p_size == 3)
        if (p_size <= 1)
          {
            //const unsigned p_rank = stk::parallel_machine_rank( pm );
            //const unsigned p_size = stk::parallel_machine_size( pm );

            const unsigned n = 3;
            //const unsigned nx = n , ny = n , nz = p_size*n ;
            const unsigned nx = n , ny = n;

            bool sidesets=true;
            percept::QuadFixture<double> fixture( pm , nx , ny, sidesets);
            fixture.meta_data.commit();
            fixture.generate_mesh();

            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data);
            eMesh.printInfo("quad fixture",  2);
            eMesh.saveAs(input_files_loc+"quad_2_smooth.0.e");

            unsigned center_node_id[4] = {6,7,10,11};
            for (int ii=0; ii < 4; ii++)
              {
                stk::mesh::Entity* node = eMesh.getBulkData()->get_entity(0, center_node_id[ii]);
                double * data = stk::mesh::field_data( *eMesh.getCoordinatesField() , *node );
                std::cout << "tmp srk  center node= " << data[0] << " " << data[1] << std::endl;
                data[0] += .02*(ii+1);
                data[1] += .03*(ii+1);
              }

            eMesh.saveAs(input_files_loc+"quad_2_smooth.0_perturbed.e");

            stk::mesh::Selector boundarySelector_1(*eMesh.getNonConstPart("surface_1") );
            stk::mesh::Selector boundarySelector_2(*eMesh.getNonConstPart("surface_2") );
            stk::mesh::Selector boundarySelector_3(*eMesh.getNonConstPart("surface_3") );
            stk::mesh::Selector boundarySelector_4(*eMesh.getNonConstPart("surface_4") );
            stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4;

            bool do_jacobi = false;
            PerceptMesquiteMesh pmm(&eMesh, 0, &boundarySelector);
            //pmm.setDoJacobiIterations(do_jacobi);
            PerceptMesquiteMeshDomain pmd(&eMesh, 0);
            percept::PMMLaplaceSmoother1 ls;
            if (do_jacobi) 
              ls.get_smoother().do_jacobi_optimization();
            ls.run(pmm, pmd);

            //STKUNIT_EXPECT_DOUBLE_EQ_APPROX(data[0], 1.0);
            //STKUNIT_EXPECT_DOUBLE_EQ_APPROX(data[1], 1.0);

            eMesh.saveAs(output_files_loc+"quad_2_smooth.1.e");
          }

        if (p_size <= 1)
          {
            percept::PerceptMesh eMesh(2);
            eMesh.open(input_files_loc+"quad_2_smooth.0_perturbed.e");
            eMesh.commit();

            stk::mesh::Selector boundarySelector_1(*eMesh.getNonConstPart("surface_1") );
            stk::mesh::Selector boundarySelector_2(*eMesh.getNonConstPart("surface_2") );
            stk::mesh::Selector boundarySelector_3(*eMesh.getNonConstPart("surface_3") );
            stk::mesh::Selector boundarySelector_4(*eMesh.getNonConstPart("surface_4") );
            stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4;

            //bool do_jacobi = true;
            Mesquite::MsqDebug::enable(1);
            Mesquite::MsqDebug::enable(2);
            Mesquite::MsqDebug::enable(3);
            PerceptMesquiteMesh pmm(&eMesh, 0, &boundarySelector);
            //pmm.setDoJacobiIterations(do_jacobi);
            PerceptMesquiteMeshDomain pmd(&eMesh, 0);
            percept::PMMShapeImprover si;
            si.run(pmm, pmd);
            eMesh.saveAs(output_files_loc+"quad_si_smooth.1.e");

//             percept::PMMShapeImprover si2;
//             si2.run(pmm, pmd);
//             eMesh.saveAs(output_files_loc+"quad_si2_smooth.1.e");

            //STKUNIT_EXPECT_DOUBLE_EQ_APPROX(data[0], 1.0);
            //STKUNIT_EXPECT_DOUBLE_EQ_APPROX(data[1], 1.0);

          }
      }


    }
  }
}
