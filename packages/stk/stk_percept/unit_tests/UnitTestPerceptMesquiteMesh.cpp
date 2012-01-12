/*--------------------------------------------------------------------*/
/*    Copyright 2009, 2011 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#if defined(STK_BUILT_IN_SIERRA) && !defined(__IBMCPP__)
#define STK_ADAPT_HAS_GEOMETRY
#else
#undef STK_ADAPT_HAS_GEOMETRY
#endif

#if defined( STK_ADAPT_HAS_GEOMETRY )

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
      static int s_par_size_max = 1;

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
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        MPI_Barrier( MPI_COMM_WORLD );
        unsigned par_size_max = 1;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        //if (p_size == 1 || p_size == 2)
        if (p_size <= par_size_max)
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
#if 0
            eMesh.reopen();
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
            eMesh.commit();
#endif

            eMesh.printInfo("quad fixture",  2);
            eMesh.saveAs(input_files_loc+"quad_smooth.0.e");

            unsigned center_node_id = 5;
            stk::mesh::Entity* node = eMesh.getBulkData()->get_entity(0, center_node_id);
            double *data = 0;
            if (node)
              {
                data = stk::mesh::field_data( *eMesh.getCoordinatesField() , *node );
                std::cout << "tmp srk  center node= " << data[0] << " " << data[1] << std::endl;
                data[0] += .2;
                data[1] += .3;
              }

            eMesh.saveAs(output_files_loc+"quad_smooth.0_perturbed.e");

            stk::mesh::Selector boundarySelector_1(*eMesh.getNonConstPart("surface_1") );
            stk::mesh::Selector boundarySelector_2(*eMesh.getNonConstPart("surface_2") );
            stk::mesh::Selector boundarySelector_3(*eMesh.getNonConstPart("surface_3") );
            stk::mesh::Selector boundarySelector_4(*eMesh.getNonConstPart("surface_4") );
            stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4;

            if (p_size == 1)
              {
                PerceptMesquiteMesh pmm(&eMesh, 0, &boundarySelector);
                PerceptMesquiteMeshDomain pmd(&eMesh, 0);
                percept::PMMLaplaceSmoother ls;
                ls.run(pmm, pmd);
              }
            else
              {
                PerceptMesquiteMesh pmm0(&eMesh, 0, &boundarySelector);
                PerceptMesquiteMesh::PMMParallelMesh pmm(&pmm0);
                PerceptMesquiteMeshDomain pmd(&eMesh, 0);
                percept::PMMLaplaceSmoother ls;
                ls.run(pmm, pmd);
              }

            if (node)
              {
                STKUNIT_EXPECT_DOUBLE_EQ_APPROX(data[0], 1.0);
                STKUNIT_EXPECT_DOUBLE_EQ_APPROX(data[1], 1.0);
              }

            eMesh.saveAs(output_files_loc+"quad_smooth.1.e");
          }
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      STKUNIT_UNIT_TEST(unit_perceptMesquite, quad_2)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        MPI_Barrier( MPI_COMM_WORLD );
        unsigned par_size_max = 1;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        //if (p_size == 1 || p_size == 2)
          if (p_size <= par_size_max)
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
                if (node)
                  {
                    double * data = stk::mesh::field_data( *eMesh.getCoordinatesField() , *node );
                    std::cout << "tmp srk  center node= " << data[0] << " " << data[1] << std::endl;
                    data[0] += .02*(ii+1);
                    data[1] += .03*(ii+1);
                  }
              }

            eMesh.saveAs(input_files_loc+"quad_2_smooth.0_perturbed.e");

            stk::mesh::Selector boundarySelector_1(*eMesh.getNonConstPart("surface_1") );
            stk::mesh::Selector boundarySelector_2(*eMesh.getNonConstPart("surface_2") );
            stk::mesh::Selector boundarySelector_3(*eMesh.getNonConstPart("surface_3") );
            stk::mesh::Selector boundarySelector_4(*eMesh.getNonConstPart("surface_4") );
            stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4;

            bool do_jacobi = false;

            if (p_size==1)
              {
                PerceptMesquiteMesh pmm(&eMesh, 0, &boundarySelector);
                PerceptMesquiteMeshDomain pmd(&eMesh, 0);
                percept::PMMLaplaceSmoother1 ls;
                if (do_jacobi) 
                  ls.get_smoother().do_jacobi_optimization();
                ls.run(pmm, pmd);
              }
            else
              {
                PerceptMesquiteMesh pmm0(&eMesh, 0, &boundarySelector);
                PerceptMesquiteMesh::PMMParallelMesh pmm(&pmm0);
                PerceptMesquiteMeshDomain pmd(&eMesh, 0);
                percept::PMMLaplaceSmoother1 ls;
                if (do_jacobi) 
                  ls.get_smoother().do_jacobi_optimization();
                ls.run(pmm, pmd);
              }

            //STKUNIT_EXPECT_DOUBLE_EQ_APPROX(data[0], 1.0);
            //STKUNIT_EXPECT_DOUBLE_EQ_APPROX(data[1], 1.0);

            eMesh.saveAs(output_files_loc+"quad_2_smooth.1.e");
          }

        if (p_size <= par_size_max)
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
            //Mesquite::MsqDebug::enable(3);
            if (p_size == 1)
              {
                PerceptMesquiteMesh pmm(&eMesh, 0, &boundarySelector);
                //pmm.setDoJacobiIterations(do_jacobi);
                PerceptMesquiteMeshDomain pmd(&eMesh, 0);
                percept::PMMShapeImprover si;
                si.run(pmm, pmd);
              }
            else
              {
                PerceptMesquiteMesh pmm0(&eMesh, 0, &boundarySelector);
                PerceptMesquiteMesh::PMMParallelMesh pmm(&pmm0);
                //pmm.setDoJacobiIterations(do_jacobi);
                PerceptMesquiteMeshDomain pmd(&eMesh, 0);
                percept::PMMShapeImprover si;
                si.run(pmm, pmd);
              }
            eMesh.saveAs(output_files_loc+"quad_si_smooth.1.e");

            //STKUNIT_EXPECT_DOUBLE_EQ_APPROX(data[0], 1.0);
            //STKUNIT_EXPECT_DOUBLE_EQ_APPROX(data[1], 1.0);

          }
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      STKUNIT_UNIT_TEST(unit_perceptMesquite, quad_3)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        MPI_Barrier( MPI_COMM_WORLD );
        unsigned par_size_max = 1;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        //if (p_size == 1 || p_size == 2)
          if (p_size <= par_size_max)
          {
            const unsigned n = 12;
            //const unsigned nx = n , ny = n , nz = p_size*n ;
            const unsigned nx = n , ny = n;

            bool sidesets=true;
            percept::QuadFixture<double> fixture( pm , nx , ny, sidesets);
            fixture.meta_data.commit();
            fixture.generate_mesh();

            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data);
            eMesh.printInfo("quad fixture",  2);
            eMesh.saveAs(input_files_loc+"quad_3_smooth.0.e");


            stk::mesh::Selector boundarySelector_1(*eMesh.getNonConstPart("surface_1") );
            stk::mesh::Selector boundarySelector_2(*eMesh.getNonConstPart("surface_2") );
            stk::mesh::Selector boundarySelector_3(*eMesh.getNonConstPart("surface_3") );
            stk::mesh::Selector boundarySelector_4(*eMesh.getNonConstPart("surface_4") );
            stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4;

            const std::vector<stk::mesh::Bucket*> & buckets = eMesh.getBulkData()->buckets( stk::mesh::fem::FEMMetaData::NODE_RANK );

            for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
              {
                if (boundarySelector_1(**k)) 
                  {
                    stk::mesh::Bucket & bucket = **k ;

                    const unsigned num_elements_in_bucket = bucket.size();
                
                    for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                      {
                        stk::mesh::Entity& entity = bucket[iEntity];

                        double * data = stk::mesh::field_data( *eMesh.getCoordinatesField() , entity );
                        //data[0] += .02*(ii+1);
                        double ix = data[0]/double(n);
                        data[1] += (ix)*(1.0-ix)*0.8*double(n);
                        std::cout << "tmp srk surface 1 node = " << data[0] << " " << data[1] << std::endl;
                      }
                  }
              }

            eMesh.saveAs(input_files_loc+"quad_3_smooth.0_perturbed.e");

            //bool do_jacobi = true;
            Mesquite::MsqDebug::enable(1);
            Mesquite::MsqDebug::enable(2);
            //Mesquite::MsqDebug::enable(3);
            if (p_size == 1)
              {
                PerceptMesquiteMesh pmm(&eMesh, 0, &boundarySelector);
                PerceptMesquiteMeshDomain pmd(&eMesh, 0);
                percept::PMMShapeImprover si;
                si.run(pmm, pmd);
              }
            else
              {
                PerceptMesquiteMesh pmm0(&eMesh, 0, &boundarySelector);
                PerceptMesquiteMesh::PMMParallelMesh pmm(&pmm0);
                //pmm.setDoJacobiIterations(do_jacobi);
                PerceptMesquiteMeshDomain pmd(&eMesh, 0);
                percept::PMMShapeImprover si;
                si.run(pmm, pmd);
              }
            eMesh.saveAs(output_files_loc+"quad_3_si_smooth.1.e");

            //STKUNIT_EXPECT_DOUBLE_EQ_APPROX(data[0], 1.0);
            //STKUNIT_EXPECT_DOUBLE_EQ_APPROX(data[1], 1.0);


          }
      }


      //=============================================================================
      //=============================================================================
      //=============================================================================

      STKUNIT_UNIT_TEST(unit_perceptMesquite, hex_1)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        MPI_Barrier( MPI_COMM_WORLD );
        unsigned par_size_max = s_par_size_max;

        const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        //if (p_size == 1 || p_size == 3)
        if (p_size <= par_size_max)
          {
            //const unsigned p_rank = stk::parallel_machine_rank( pm );
            //const unsigned p_size = stk::parallel_machine_size( pm );

            unsigned n = 4;
            std::cout << "P["<<p_rank<<"] " << "tmp srk doing Laplace smoothing for hex_1 case, n = " << n << std::endl;
            unsigned nn = n+1;
            //std::string gmesh_spec = toString(n)+"x"+toString(n)+"x"+toString(n*p_size)+std::string("|bbox:0,0,0,1,1,1|sideset:xXyYzZ");
            std::string gmesh_spec = toString(n)+"x"+toString(n)+"x"+toString(n*2)+std::string("|bbox:0,0,0,1,1,1|sideset:xXyYzZ");
            PerceptMesh eMesh(3);
            eMesh.newMesh(percept::GMeshSpec(gmesh_spec));
            eMesh.commit();
            eMesh.saveAs(input_files_loc+"hex_1_smooth.0.e");

            stk::mesh::Selector boundarySelector_1(*eMesh.getNonConstPart("surface_1") );
            stk::mesh::Selector boundarySelector_2(*eMesh.getNonConstPart("surface_2") );
            stk::mesh::Selector boundarySelector_3(*eMesh.getNonConstPart("surface_3") );
            stk::mesh::Selector boundarySelector_4(*eMesh.getNonConstPart("surface_4") );
            stk::mesh::Selector boundarySelector_5(*eMesh.getNonConstPart("surface_5") );
            stk::mesh::Selector boundarySelector_6(*eMesh.getNonConstPart("surface_6") );
            stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4 | boundarySelector_5 | boundarySelector_6;

            double delta_max = 0.01/(double(n));
            for (unsigned ii=1; ii <= (nn*nn*nn); ii++)
              {
                stk::mesh::Entity* node = eMesh.getBulkData()->get_entity(0, ii);
                if (node)
                  {
                    if (boundarySelector(*node)) continue;
                    double * data = stk::mesh::field_data( *eMesh.getCoordinatesField() , *node );
                    std::cout << "P["<<p_rank<<"] " << "tmp srk  center node= " << data[0] << " " << data[1] << std::endl;
                    data[0] += delta_max*double(ii)/double(n);
                    data[1] += 2*delta_max*double(ii)/double(n);
                    data[2] += 3*delta_max*double(ii)/double(n);
                  }
              }

            //std::cout << "P["<<p_rank<<"] " << "tmp srk doing Laplace smoothing for hex_1 case, n = " << n << std::endl;

            eMesh.saveAs(input_files_loc+"hex_1_smooth.0_perturbed.e");

            bool do_jacobi = true;
            int numIterMax = 20;
            if (p_size == 1)
              {
                PerceptMesquiteMesh pmm(&eMesh, 0, &boundarySelector);
                PerceptMesquiteMeshDomain pmd(&eMesh, 0);
                percept::PMMLaplaceSmoother1 ls(numIterMax);
                if (do_jacobi) 
                  ls.get_smoother().do_jacobi_optimization();
                ls.run(pmm, pmd);
              }
            else
              {
                PerceptMesquiteMesh pmm0(&eMesh, 0, &boundarySelector);
                PerceptMesquiteMesh::PMMParallelMesh pmm(&pmm0);
                PerceptMesquiteMeshDomain pmd(&eMesh, 0);
                percept::PMMLaplaceSmoother1 ls(numIterMax);
                if (do_jacobi) 
                  ls.get_smoother().do_jacobi_optimization();
                ls.run(pmm, pmd);

              }
            std::cout << "tmp srk doing Laplace smoothing for hex_1 case ... done " << std::endl;

            eMesh.saveAs(output_files_loc+"hex_1_smooth.1.e");
          }

        MPI_Barrier( MPI_COMM_WORLD );
        std::cout << "\n\n P["<<p_rank<<"] " << "tmp srk BARRIER doing Shape smoothing for hex_1 case... " <<  std::endl;
        MPI_Barrier( MPI_COMM_WORLD );

        if (p_size <= par_size_max)
          {
            percept::PerceptMesh eMesh(3);
            eMesh.open(input_files_loc+"hex_1_smooth.0_perturbed.e");
            eMesh.commit();

            std::cout << "tmp srk doing Shape smoothing for hex_1 case..." << std::endl;

            stk::mesh::Selector boundarySelector_1(*eMesh.getNonConstPart("surface_1") );
            stk::mesh::Selector boundarySelector_2(*eMesh.getNonConstPart("surface_2") );
            stk::mesh::Selector boundarySelector_3(*eMesh.getNonConstPart("surface_3") );
            stk::mesh::Selector boundarySelector_4(*eMesh.getNonConstPart("surface_4") );
            stk::mesh::Selector boundarySelector_5(*eMesh.getNonConstPart("surface_5") );
            stk::mesh::Selector boundarySelector_6(*eMesh.getNonConstPart("surface_6") );
            stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4 | boundarySelector_5 | boundarySelector_6;

            //bool do_jacobi = true;
            Mesquite::MsqDebug::enable(1);
            Mesquite::MsqDebug::enable(2);
            //Mesquite::MsqDebug::enable(3);
            int  msq_debug             = 2; // 1,2,3 for more debug info
            bool always_smooth         = true;

            if (p_size == 1)
              {
                PerceptMesquiteMesh pmm(&eMesh, 0, &boundarySelector);
                PerceptMesquiteMeshDomain pmd(&eMesh, 0);
                percept::PMMShapeImprover si;
                si.run(pmm, pmd, always_smooth, msq_debug);
              }
            else
              {
                PerceptMesquiteMesh pmm0(&eMesh, 0, &boundarySelector);
                PerceptMesquiteMesh::PMMParallelMesh pmm(&pmm0);
                PerceptMesquiteMeshDomain pmd(&eMesh, 0);
                percept::PMMShapeImprover si;
                si.run(pmm, pmd, always_smooth, msq_debug);
              }
            std::cout << "tmp srk doing Shape smoothing for hex_1 case... done " << std::endl;
            eMesh.saveAs(output_files_loc+"hex_1_si_smooth.1.e");
          }
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      STKUNIT_UNIT_TEST(unit_perceptMesquite, hex_2)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        MPI_Barrier( MPI_COMM_WORLD );
        unsigned par_size_max = s_par_size_max;

        const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        //if (p_size == 1 || p_size == 3)
        if (p_size <= par_size_max)
          {
            unsigned n = 5;
            std::cout << "P["<<p_rank<<"] " << "tmp srk doing Laplace smoothing for hex_1 case, n = " << n << std::endl;
            //unsigned nn = n+1;
            //std::string gmesh_spec = toString(n)+"x"+toString(n)+"x"+toString(n*p_size)+std::string("|bbox:0,0,0,1,1,1|sideset:xXyYzZ");
            std::string gmesh_spec = toString(n)+"x"+toString(n)+"x"+toString(n)+std::string("|bbox:0,0,0,1,1,1|sideset:xXyYzZ");
            PerceptMesh eMesh(3);
            eMesh.newMesh(percept::GMeshSpec(gmesh_spec));
            eMesh.commit();
            eMesh.saveAs(input_files_loc+"hex_2_smooth.0.e");

            stk::mesh::Selector boundarySelector_1(*eMesh.getNonConstPart("surface_1") );
            stk::mesh::Selector boundarySelector_2(*eMesh.getNonConstPart("surface_2") );
            stk::mesh::Selector boundarySelector_3(*eMesh.getNonConstPart("surface_3") );
            stk::mesh::Selector boundarySelector_4(*eMesh.getNonConstPart("surface_4") );
            stk::mesh::Selector boundarySelector_5(*eMesh.getNonConstPart("surface_5") );
            stk::mesh::Selector boundarySelector_6(*eMesh.getNonConstPart("surface_6") );
            stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4 | boundarySelector_5 | boundarySelector_6;

            const std::vector<stk::mesh::Bucket*> & buckets = eMesh.getBulkData()->buckets( stk::mesh::fem::FEMMetaData::NODE_RANK );

            for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
              {
                if (boundarySelector_5(**k)) 
                  {
                    stk::mesh::Bucket & bucket = **k ;

                    const unsigned num_elements_in_bucket = bucket.size();
                
                    for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                      {
                        stk::mesh::Entity& entity = bucket[iEntity];

                        double * data = stk::mesh::field_data( *eMesh.getCoordinatesField() , entity );
                        //data[0] += .02*(ii+1);
                        double ix = data[0];
                        double iy = data[1];
                        //data[2] += (ix)*(1.0-ix)*(iy)*(1.0-iy)*2.0*4.;
                        data[2] += (ix)*(1.0-ix)*(iy)*(1.0-iy)*2.0*.5;
                        std::cout << "tmp srk surface 1 hex_2 node = " << data[0] << " " << data[1] << " " << data[2] << std::endl;
                      }
                  }
              }
            eMesh.saveAs(input_files_loc+"hex_2_smooth.0_perturbed.e");

            std::cout << "tmp srk doing Shape smoothing for hex_2 case..." << std::endl;

            //bool do_jacobi = true;
            Mesquite::MsqDebug::enable(1);
            Mesquite::MsqDebug::enable(2);
            //Mesquite::MsqDebug::enable(3);

            int  msq_debug             = 2; // 1,2,3 for more debug info
            bool always_smooth         = true;

            if (p_size == 1) 
              {
                PerceptMesquiteMesh pmm(&eMesh, 0, &boundarySelector);
                PerceptMesquiteMeshDomain pmd(&eMesh, 0);
                percept::PMMShapeImprover si;
                si.run(pmm, pmd, always_smooth, msq_debug);
              }
            else
              {
                PerceptMesquiteMesh pmm0(&eMesh, 0, &boundarySelector);
                PerceptMesquiteMesh::PMMParallelMesh pmm(&pmm0);
                PerceptMesquiteMeshDomain pmd(&eMesh, 0);
                percept::PMMShapeImprover si;
                si.run(pmm, pmd, always_smooth, msq_debug);
              }
            MPI_Barrier( MPI_COMM_WORLD );
            std::cout << "\n\n P["<<p_rank<<"] " << "tmp srk BARRIER done Shape smoothing for hex_2 case...done " <<  std::endl;
            MPI_Barrier( MPI_COMM_WORLD );

            eMesh.saveAs(output_files_loc+"hex_2_si_smooth.1.e");
          }
      }

    }
  }
}
#endif

