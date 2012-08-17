/*--------------------------------------------------------------------*/
/*    Copyright 2009, 2011 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/


#include <stk_percept/Percept.hpp>

#if defined( STK_PERCEPT_HAS_GEOMETRY ) && defined(STK_PERCEPT_HAS_MESQUITE)

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
#include <stk_percept/mesh/mod/mesquite-interface/PMMParallelShapeImprover.hpp>
#include <MsqDebug.hpp>

#include "MeshImpl.hpp"
#include "MsqTimer.hpp"
#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "Vector3D.hpp"
#include "InstructionQueue.hpp"
#include "LaplaceWrapper.hpp"
#include "PatchData.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"

/* Mesquite includes */
#include "ParallelMeshImpl.hpp"
#include "ParallelHelper.hpp"


// algorithms
#include "Randomize.hpp"
#include "ConditionNumberQualityMetric.hpp"
#include "UntangleBetaQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "LInfTemplate.hpp"
#include "SteepestDescent.hpp"
#include "ConjugateGradient.hpp"
#include "PlanarDomain.hpp"

#define StackTrace StackTraceTmp

#include <iostream>
#include <cstdlib>
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
      Mesquite::MeshImpl *create_mesquite_mesh(PerceptMesh *eMesh, stk::mesh::Selector *boundarySelector);


#define DO_TESTS 0
#if DO_TESTS

#define EXTRA_PRINT 0
      static int s_par_size_max = 2;

      //static int print_infoLevel = 0;

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

      // run Laplace smoother on an interior-node perturbed domain
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

            eMesh.reopen();
            eMesh.addParallelInfoFields(true,true);
            eMesh.commit();

            stk::mesh::Selector boundarySelector_1(*eMesh.get_non_const_part("surface_1") );
            stk::mesh::Selector boundarySelector_2(*eMesh.get_non_const_part("surface_2") );
            stk::mesh::Selector boundarySelector_3(*eMesh.get_non_const_part("surface_3") );
            stk::mesh::Selector boundarySelector_4(*eMesh.get_non_const_part("surface_4") );
            stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4;

            eMesh.populateParallelInfoFields(true,true,&boundarySelector);

            eMesh.print_info("quad fixture",  2);
            eMesh.save_as(input_files_loc+"quad_smooth.0.e");

            unsigned center_node_id = 5;
            stk::mesh::Entity* node = eMesh.get_bulk_data()->get_entity(0, center_node_id);
            double *data = 0;
            if (node)
              {
                data = stk::mesh::field_data( *eMesh.get_coordinates_field() , *node );
                //std::cout << "tmp srk  center node= " << data[0] << " " << data[1] << std::endl;
                data[0] += .2;
                data[1] += .3;
              }

            eMesh.save_as(input_files_loc+"quad_smooth.0_perturbed.e");


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

            eMesh.save_as(output_files_loc+"quad_smooth.1.e");
          }
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      // run LaplaceSmoother1, and ShapeImprover on an interior-node perturbed domain
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
            eMesh.print_info("quad fixture",  2);
            eMesh.save_as(input_files_loc+"quad_2_smooth.0.e");

            unsigned center_node_id[4] = {6,7,10,11};
            for (int ii=0; ii < 4; ii++)
              {
                stk::mesh::Entity* node = eMesh.get_bulk_data()->get_entity(0, center_node_id[ii]);
                if (node)
                  {
                    double * data = stk::mesh::field_data( *eMesh.get_coordinates_field() , *node );
                    //std::cout << "tmp srk  center node= " << data[0] << " " << data[1] << std::endl;
                    data[0] += .02*(ii+1);
                    data[1] += .03*(ii+1);
                  }
              }

            eMesh.save_as(input_files_loc+"quad_2_smooth.0_perturbed.e");

            stk::mesh::Selector boundarySelector_1(*eMesh.get_non_const_part("surface_1") );
            stk::mesh::Selector boundarySelector_2(*eMesh.get_non_const_part("surface_2") );
            stk::mesh::Selector boundarySelector_3(*eMesh.get_non_const_part("surface_3") );
            stk::mesh::Selector boundarySelector_4(*eMesh.get_non_const_part("surface_4") );
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

            eMesh.save_as(output_files_loc+"quad_2_smooth.1.e");
          }

        if (p_size <= par_size_max)
          {
            percept::PerceptMesh eMesh(2);
            eMesh.open(input_files_loc+"quad_2_smooth.0_perturbed.e");
            eMesh.commit();

            stk::mesh::Selector boundarySelector_1(*eMesh.get_non_const_part("surface_1") );
            stk::mesh::Selector boundarySelector_2(*eMesh.get_non_const_part("surface_2") );
            stk::mesh::Selector boundarySelector_3(*eMesh.get_non_const_part("surface_3") );
            stk::mesh::Selector boundarySelector_4(*eMesh.get_non_const_part("surface_4") );
            stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4;

            //bool do_jacobi = true;
            //Mesquite::MsqDebug::enable(1);
            //Mesquite::MsqDebug::enable(2);
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
            eMesh.save_as(output_files_loc+"quad_si_smooth.1.e");

            //STKUNIT_EXPECT_DOUBLE_EQ_APPROX(data[0], 1.0);
            //STKUNIT_EXPECT_DOUBLE_EQ_APPROX(data[1], 1.0);

          }
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      // run ShapeImprover on a domain with one side perturbed
      STKUNIT_UNIT_TEST(unit_perceptMesquite, quad_3)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        MPI_Barrier( MPI_COMM_WORLD );
        unsigned par_size_max = s_par_size_max;
        
        const unsigned p_rank = stk::parallel_machine_rank( pm );
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
            eMesh.print_info("quad fixture",  2);
            eMesh.save_as(input_files_loc+"quad_3_smooth.0.e");

            eMesh.reopen();
            eMesh.addParallelInfoFields(true,true);
            eMesh.commit();

            stk::mesh::Selector boundarySelector_1(*eMesh.get_non_const_part("surface_1") );
            stk::mesh::Selector boundarySelector_2(*eMesh.get_non_const_part("surface_2") );
            stk::mesh::Selector boundarySelector_3(*eMesh.get_non_const_part("surface_3") );
            stk::mesh::Selector boundarySelector_4(*eMesh.get_non_const_part("surface_4") );
            stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4;

            eMesh.populateParallelInfoFields(true,true,&boundarySelector);

            const std::vector<stk::mesh::Bucket*> & buckets = eMesh.get_bulk_data()->buckets( stk::mesh::fem::FEMMetaData::NODE_RANK );

            for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
              {
                if (boundarySelector_1(**k)) 
                  {
                    stk::mesh::Bucket & bucket = **k ;

                    const unsigned num_elements_in_bucket = bucket.size();
                
                    for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                      {
                        stk::mesh::Entity& entity = bucket[iEntity];

                        double * data = stk::mesh::field_data( *eMesh.get_coordinates_field() , entity );
                        double ix = data[0]/double(n);
                        data[1] += (ix)*(1.0-ix)*0.8*double(n);
                        //std::cout << "tmp srk surface 1 node = " << data[0] << " " << data[1] << std::endl;
                      }
                  }
              }

            eMesh.save_as(input_files_loc+"quad_3_smooth.0_perturbed.e");

            {
              Mesquite::MsqError err;
              Mesquite::MeshImpl *msqMesh = create_mesquite_mesh(&eMesh, &boundarySelector);
              std::ostringstream vtk_file;
              vtk_file << "par_original_quad_mesh." << p_size << "." << p_rank << ".vtk";
              msqMesh->write_vtk(vtk_file.str().c_str(), err); 
              if (err) {std::cout << err << endl;  STKUNIT_EXPECT_TRUE(false);}
            }

            //bool do_jacobi = true;
            //Mesquite::MsqDebug::enable(1);
            //Mesquite::MsqDebug::enable(2);
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
            eMesh.save_as(output_files_loc+"quad_3_si_smooth.1.e");

            {
              Mesquite::MsqError err;
              Mesquite::MeshImpl *msqMesh = create_mesquite_mesh(&eMesh, &boundarySelector);
              std::ostringstream vtk_file;
              vtk_file << "par_smoothed_quad_mesh." << p_size << "." << p_rank << ".vtk";
              msqMesh->write_vtk(vtk_file.str().c_str(), err); 
              if (err) {std::cout << err << endl;  STKUNIT_EXPECT_TRUE(false);}
            }

          }
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      // run PMMParallelShapeImprover on a domain with one side perturbed
      STKUNIT_UNIT_TEST(unit_perceptMesquite, quad_4)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        MPI_Barrier( MPI_COMM_WORLD );
        unsigned par_size_max = s_par_size_max;
        
        const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        //if (p_size == 1 || p_size == 2)
        if (p_size <= par_size_max)
          {
            const unsigned nele = 12;
            //const unsigned nx = nele , ny = nele , nz = p_size*nele ;
            const unsigned nx = nele , ny = nele;

            bool sidesets=true;
            percept::QuadFixture<double> fixture( pm , nx , ny, sidesets);
            fixture.meta_data.commit();
            fixture.generate_mesh();

            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data);
            //eMesh.print_info("quad fixture",  2);
            eMesh.save_as(input_files_loc+"quad_4_smooth.0.e");

            eMesh.reopen();
            eMesh.add_coordinate_state_fields();
            stk::mesh::FieldBase *proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), 0);
            //eMesh.addParallelInfoFields(true,true);
            eMesh.commit();
            eMesh.set_proc_rank_field(proc_rank_field);

            stk::mesh::Selector boundarySelector_1(*eMesh.get_non_const_part("surface_1") );
            stk::mesh::Selector boundarySelector_2(*eMesh.get_non_const_part("surface_2") );
            stk::mesh::Selector boundarySelector_3(*eMesh.get_non_const_part("surface_3") );
            stk::mesh::Selector boundarySelector_4(*eMesh.get_non_const_part("surface_4") );
            stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4;

            //eMesh.populateParallelInfoFields(true,true,&boundarySelector);

            const std::vector<stk::mesh::Bucket*> & buckets = eMesh.get_bulk_data()->buckets( stk::mesh::fem::FEMMetaData::NODE_RANK );

            for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
              {
                  {
                    stk::mesh::Bucket & bucket = **k ;

                    const unsigned num_elements_in_bucket = bucket.size();
                
                    for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                      {
                        stk::mesh::Entity& node = bucket[iEntity];

                        double * data = stk::mesh::field_data( *eMesh.get_coordinates_field() , node );
                        double iy = data[1]/double(nele);
                        iy = iy*iy;
                        data[1] = iy*double(nele);
                      }
                  }
              }

            eMesh.save_as(input_files_loc+"quad_4_smooth.1.e");

            // save state of original mesh
            // field, dst, src: 
            eMesh.copy_field(eMesh.get_field("coordinates_NM1"), eMesh.get_coordinates_field());

            for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
              {
                if (boundarySelector_1(**k)) 
                  {
                    stk::mesh::Bucket & bucket = **k ;

                    const unsigned num_elements_in_bucket = bucket.size();
                
                    for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                      {
                        stk::mesh::Entity& node = bucket[iEntity];

                        double * data = stk::mesh::field_data( *eMesh.get_coordinates_field() , node );
                        double ix = data[0]/double(nele);
                        //double bump_size=2.8; // 0.8
                        double bump_size=2.8; // 0.8
                        data[1] += (ix)*(1.0-ix)*bump_size*double(nele);
                        //std::cout << "tmp srk surface 1 node = " << data[0] << " " << data[1] << std::endl;
                      }
                  }
              }

            // save state of projected mesh
            // field, dst, src: 
            eMesh.copy_field(eMesh.get_field("coordinates_N"), eMesh.get_coordinates_field());

            eMesh.save_as(input_files_loc+"quad_4_smooth.0_perturbed.e");

            if (0)
            {
              Mesquite::MsqError err;
              Mesquite::MeshImpl *msqMesh = create_mesquite_mesh(&eMesh, &boundarySelector);
              std::ostringstream vtk_file;
              vtk_file << "par_original_quad_4_mesh." << p_size << "." << p_rank << ".vtk";
              msqMesh->write_vtk(vtk_file.str().c_str(), err); 
              if (err) {std::cout << err << endl;  STKUNIT_EXPECT_TRUE(false);}
            }

            int  msq_debug             = 0; // 1,2,3 for more debug info
            bool always_smooth         = true;
            //bool do_jacobi = true;
            Mesquite::MsqDebug::disable(1);
            //Mesquite::MsqDebug::enable(2);
            //Mesquite::MsqDebug::enable(3);
              {
                PerceptMesquiteMesh pmm(&eMesh, 0, &boundarySelector);
                //PerceptMesquiteMeshDomain pmd(&eMesh, 0);
                PlanarDomain planar_domain(PlanarDomain::XY);
                //PMMParallelShapeImprover(int innerIter=100, double gradNorm = 1.e-8, int parallelIterations=20) : 
                if (1)
                  {
                    percept::PMMParallelShapeImprover pmmpsi(1001, 1.e-4, 1);
                    //pmmpsi.run(pmm, &pmd, always_smooth, msq_debug);
                    pmmpsi.run(pmm, &planar_domain, always_smooth, msq_debug);
                  }
                else
                  {
                    percept::PMMLaplaceSmoother1 ls;
                    //if (do_jacobi) 
                    //  ls.get_smoother().do_jacobi_optimization();
                    ls.run(pmm, planar_domain);
                  }
              }

            eMesh.save_as(output_files_loc+"quad_4_si_smooth.1.e");

            if (0)
            {
              Mesquite::MsqError err;
              Mesquite::MeshImpl *msqMesh = create_mesquite_mesh(&eMesh, &boundarySelector);
              std::ostringstream vtk_file;
              vtk_file << "par_smoothed_quad_4_mesh." << p_size << "." << p_rank << ".vtk";
              msqMesh->write_vtk(vtk_file.str().c_str(), err); 
              if (err) {std::cout << err << endl;  STKUNIT_EXPECT_TRUE(false);}
            }

          }
      }


      //=============================================================================
      //=============================================================================
      //=============================================================================

      // a cube with an internally distorted mesh

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
            std::cout << "P[" << p_rank << "] " << "tmp srk doing Laplace smoothing for hex_1 case, n = " << n << std::endl;
            unsigned nn = n+1;
            //std::string gmesh_spec = toString(n)+"x"+toString(n)+"x"+toString(n*p_size)+std::string("|bbox:0,0,0,1,1,1|sideset:xXyYzZ");
            //std::string gmesh_spec = toString(n)+"x"+toString(n)+"x"+toString(n*2)+std::string("|bbox:0,0,0,1,1,1|sideset:xXyYzZ");
            std::string gmesh_spec = toString(n)+"x"+toString(n)+"x"+toString(n)+std::string("|bbox:0,0,0,1,1,1|sideset:xXyYzZ");
            PerceptMesh eMesh(3);
            eMesh.new_mesh(percept::GMeshSpec(gmesh_spec));
            eMesh.commit();
            eMesh.save_as(input_files_loc+"hex_1_smooth.0.e");

            stk::mesh::Selector boundarySelector_1(*eMesh.get_non_const_part("surface_1") );
            stk::mesh::Selector boundarySelector_2(*eMesh.get_non_const_part("surface_2") );
            stk::mesh::Selector boundarySelector_3(*eMesh.get_non_const_part("surface_3") );
            stk::mesh::Selector boundarySelector_4(*eMesh.get_non_const_part("surface_4") );
            stk::mesh::Selector boundarySelector_5(*eMesh.get_non_const_part("surface_5") );
            stk::mesh::Selector boundarySelector_6(*eMesh.get_non_const_part("surface_6") );
            stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4 | boundarySelector_5 | boundarySelector_6;

            //double delta_max = 0.01/(double(n));
            double delta_max = 0.001/(double(n));
            for (unsigned ii=1; ii <= (nn*nn*nn); ii++)
              {
                stk::mesh::Entity* node = eMesh.get_bulk_data()->get_entity(0, ii);
                if (node)
                  {
                    if (boundarySelector(*node)) continue;
                    double * data = stk::mesh::field_data( *eMesh.get_coordinates_field() , *node );
                    //std::cout << "P["<<p_rank<<"] " << "tmp srk  center node= " << data[0] << " " << data[1] << std::endl;
                    data[0] += delta_max*double(ii)/double(n);
                    data[1] += 2*delta_max*double(ii)/double(n);
                    data[2] += 3*delta_max*double(ii)/double(n);
                  }
              }

            //std::cout << "P["<<p_rank<<"] " << "tmp srk doing Laplace smoothing for hex_1 case, n = " << n << std::endl;

            eMesh.save_as(input_files_loc+"hex_1_smooth.0_perturbed.e");
            //             Mesquite::MsqDebug::enable(1);
            //             Mesquite::MsqDebug::enable(2);

            bool do_jacobi = false;
            double max_vertex_movement=1.e-8;
            int msq_debug = 0;
            if (p_size == 1)
              {
                PerceptMesquiteMesh pmm(&eMesh, 0, &boundarySelector);
                PerceptMesquiteMeshDomain pmd(&eMesh, 0);
                percept::PMMLaplaceSmoother1 ls(0.0, max_vertex_movement, 1000, false);
                if (do_jacobi) 
                  ls.get_smoother().do_jacobi_optimization();
                ls.run(pmm, pmd, true, msq_debug);
              }
            else
              {
                PerceptMesquiteMesh pmm0(&eMesh, 0, &boundarySelector);
                PerceptMesquiteMesh::PMMParallelMesh pmm(&pmm0);
                PerceptMesquiteMeshDomain pmd(&eMesh, 0);
                percept::PMMLaplaceSmoother1 ls(0.0, max_vertex_movement, 1000, false);
                if (do_jacobi) 
                  ls.get_smoother().do_jacobi_optimization();
                ls.run(pmm, pmd, true, msq_debug);
              }
            std::cout << "tmp srk doing Laplace smoothing for hex_1 case ... done " << std::endl;

            eMesh.save_as(output_files_loc+"hex_1_smooth.1.e");
          }

        if (p_size <= par_size_max)
          {
            percept::PerceptMesh eMesh(3);
            eMesh.open(input_files_loc+"hex_1_smooth.0_perturbed.e");
            eMesh.commit();

            std::cout << "tmp srk doing Shape smoothing for hex_1 case..." << std::endl;

            stk::mesh::Selector boundarySelector_1(*eMesh.get_non_const_part("surface_1") );
            stk::mesh::Selector boundarySelector_2(*eMesh.get_non_const_part("surface_2") );
            stk::mesh::Selector boundarySelector_3(*eMesh.get_non_const_part("surface_3") );
            stk::mesh::Selector boundarySelector_4(*eMesh.get_non_const_part("surface_4") );
            stk::mesh::Selector boundarySelector_5(*eMesh.get_non_const_part("surface_5") );
            stk::mesh::Selector boundarySelector_6(*eMesh.get_non_const_part("surface_6") );
            stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4 | boundarySelector_5 | boundarySelector_6;

            //bool do_jacobi = true;
            //Mesquite::MsqDebug::enable(1);
            //Mesquite::MsqDebug::enable(2);
            //Mesquite::MsqDebug::enable(3);
            int  msq_debug             = 0; // 1,2,3 for more debug info
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
            eMesh.save_as(output_files_loc+"hex_1_si_smooth.1.e");
          }
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      // Mesquite utilities...
      Mesquite::MeshImpl *create_mesquite_mesh(PerceptMesh *eMesh, stk::mesh::Selector *boundarySelector)
      {
        using namespace Mesquite;

        std::vector<size_t> gids;
        std::vector<int> pids;
        std::vector<int> fixed;
        std::vector<int> connectivity;
        std::vector<double> coords;
        unsigned num_elem=0;
        unsigned num_node=0;
        std::map<unsigned, int> local_id;
        //unsigned rank=eMesh->get_parallel_rank();
        //unsigned psize=eMesh->get_parallel_size();

        const std::vector<stk::mesh::Bucket*> & node_buckets = eMesh->get_bulk_data()->buckets( eMesh->node_rank() );
        for ( std::vector<stk::mesh::Bucket*>::const_iterator k = node_buckets.begin() ; k != node_buckets.end() ; ++k )
          {
            //if (removePartSelector(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_entity_in_bucket = bucket.size();
              for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                {
                  stk::mesh::Entity& node = bucket[ientity];
                  bool is_fixed=false;
                  if (boundarySelector && ((*boundarySelector)(node))) is_fixed = true;
                  fixed.push_back(is_fixed);
                  double * const coord = stk::mesh::field_data( *eMesh->get_coordinates_field() , node );
                  coords.push_back(coord[0]);
                  coords.push_back(coord[1]);
                  if (eMesh->get_spatial_dim()==3) 
                    coords.push_back(coord[2]);
                  else
                    coords.push_back(0.0);

                  gids.push_back(node.identifier());
                  pids.push_back(node.owner_rank());

                  local_id[node.identifier()] = num_node++;
                }
            }
          }

        const std::vector<stk::mesh::Bucket*> & buckets = eMesh->get_bulk_data()->buckets( eMesh->element_rank() );
        for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            //if (removePartSelector(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_entity_in_bucket = bucket.size();
              for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                {
                  stk::mesh::Entity& element = bucket[ientity];
                  ++num_elem;
                  const mesh::PairIterRelation elem_nodes = element.relations( stk::mesh::fem::FEMMetaData::NODE_RANK );
                  for (unsigned j = 0; j < elem_nodes.size(); j++)
                    {
                      mesh::Entity & node = * elem_nodes[ j ].entity();
                      connectivity.push_back(local_id[node.identifier()]);
                    }
                }
            }
          }
        bool *fixed_bool = new bool[num_node];
        for (unsigned ii=0; ii < num_node; ii++) fixed_bool[ii] = fixed[ii];

        Mesquite::MsqError err;
        // FIXME - works for hexes only...
        EntityTopology topo = HEXAHEDRON;
        const CellTopologyData *const topology = stk::percept::PerceptMesh::get_cell_topology(*buckets[0]);
        switch(topology->key) 
          {
          case shards::Triangle<3>::key:
            topo = TRIANGLE;
            break;
          case shards::Quadrilateral<4>::key:
            topo = QUADRILATERAL;
            break;
          case shards::Tetrahedron<4>::key:
            topo = TETRAHEDRON;
            break;
          case shards::Hexahedron<8>::key:
            topo = HEXAHEDRON;
            break;
          case shards::Wedge<6>::key:
            topo = PRISM;
            break;

          case shards::Node::key:
          case shards::Particle::key:
          case shards::Line<2>::key:
          case shards::Line<3>::key:
          case shards::ShellLine<2>::key:
          case shards::ShellLine<3>::key:
          case shards::Beam<2>::key:
          case shards::Beam<3>::key:
      
          case shards::Triangle<4>::key:
          case shards::Triangle<6>::key:
          case shards::ShellTriangle<3>::key:
          case shards::ShellTriangle<6>::key:
      
          case shards::Quadrilateral<8>::key:
          case shards::Quadrilateral<9>::key:
          case shards::ShellQuadrilateral<4>::key:
          case shards::ShellQuadrilateral<8>::key:
          case shards::ShellQuadrilateral<9>::key:
      
          case shards::Tetrahedron<8>::key:
          case shards::Tetrahedron<10>::key:
          case shards::Tetrahedron<11>::key:
      
          case shards::Hexahedron<20>::key:
          case shards::Hexahedron<27>::key:
      
          case shards::Pyramid<5>::key:
          case shards::Pyramid<13>::key:
          case shards::Pyramid<14>::key:
      
          case shards::Wedge<15>::key:
          case shards::Wedge<18>::key:
      
          case shards::Pentagon<5>::key:
          case shards::Hexagon<6>::key:
          default:
            throw std::runtime_error("unknown/unhandled topology in create_mesquite_mesh");
            break;

          }

        Mesquite::MeshImpl *mesh = new Mesquite::MeshImpl(num_node, num_elem, topo, 
                                                          fixed_bool, &coords[0], &connectivity[0]);

        std::vector<Mesh::VertexHandle> vertices;
        mesh->get_all_vertices( vertices, err); MSQ_ERRZERO(err);

        //size_t default_gid=0;
        void *gid_tag = mesh->tag_create( "GLOBAL_ID", Mesh::HANDLE, 1, 0, err );   MSQ_ERRZERO(err);
        mesh->tag_set_vertex_data( gid_tag, num_node, &vertices[0], &gids[0], err );  MSQ_ERRZERO(err);

        //int default_pid=0;
        void *pid_tag = mesh->tag_create( "PROCESSOR_ID", Mesh::INT, 1, 0, err );  MSQ_ERRZERO(err);
        mesh->tag_set_vertex_data( pid_tag, num_node, &vertices[0], &pids[0], err );  MSQ_ERRZERO(err);

        return mesh;

      }      

      //=============================================================================
      //=============================================================================
      //=============================================================================

      // A cube with an indented bump on the bottom

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
            eMesh.new_mesh(percept::GMeshSpec(gmesh_spec));
            eMesh.addParallelInfoFields(true,true);
            eMesh.commit();

            stk::mesh::Selector boundarySelector_1(*eMesh.get_non_const_part("surface_1") );
            stk::mesh::Selector boundarySelector_2(*eMesh.get_non_const_part("surface_2") );
            stk::mesh::Selector boundarySelector_3(*eMesh.get_non_const_part("surface_3") );
            stk::mesh::Selector boundarySelector_4(*eMesh.get_non_const_part("surface_4") );
            stk::mesh::Selector boundarySelector_5(*eMesh.get_non_const_part("surface_5") );
            stk::mesh::Selector boundarySelector_6(*eMesh.get_non_const_part("surface_6") );
            stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4 | boundarySelector_5 | boundarySelector_6;

            eMesh.populateParallelInfoFields(true,true,&boundarySelector);
            eMesh.save_as(input_files_loc+"hex_2_smooth.0.e");

            const std::vector<stk::mesh::Bucket*> & buckets = eMesh.get_bulk_data()->buckets( stk::mesh::fem::FEMMetaData::NODE_RANK );

            for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
              {
                if (boundarySelector_5(**k)) 
                  {
                    stk::mesh::Bucket & bucket = **k ;

                    const unsigned num_elements_in_bucket = bucket.size();
                
                    for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                      {
                        stk::mesh::Entity& entity = bucket[iEntity];

                        double * data = stk::mesh::field_data( *eMesh.get_coordinates_field() , entity );
                        double ix = data[0];
                        double iy = data[1];
                        data[2] = (ix)*(1.0-ix)*(iy)*(1.0-iy)*2.0*.5;
                      }
                  }
              }
            eMesh.save_as(input_files_loc+"hex_2_smooth.0_perturbed_small.e");

            Mesquite::MsqError err;
            Mesquite::MeshImpl *msqMesh = create_mesquite_mesh(&eMesh, &boundarySelector);
            std::ostringstream vtk_file;
            vtk_file << "par_original_hex_mesh." << p_size << "." << p_rank << ".vtk";

            msqMesh->write_vtk(vtk_file.str().c_str(), err); 
            if (err) {std::cout << err << endl;  STKUNIT_EXPECT_TRUE(false);}

            for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
              {
                if (boundarySelector_5(**k)) 
                  {
                    stk::mesh::Bucket & bucket = **k ;

                    const unsigned num_elements_in_bucket = bucket.size();
                
                    for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                      {
                        stk::mesh::Entity& entity = bucket[iEntity];

                        double * data = stk::mesh::field_data( *eMesh.get_coordinates_field() , entity );
                        double ix = data[0];
                        double iy = data[1];
                        data[2] = (ix)*(1.0-ix)*(iy)*(1.0-iy)*2.0*4.;
                      }
                  }
              }
            eMesh.save_as(input_files_loc+"hex_2_smooth.0_perturbed.e");

            msqMesh = create_mesquite_mesh(&eMesh, &boundarySelector);
            std::ostringstream vtk_file1;
            vtk_file1 << "par_untangle_original_hex_mesh." << p_size << "." << p_rank << ".vtk";

            msqMesh->write_vtk(vtk_file1.str().c_str(), err); 
            if (err) {std::cout << err << endl;  STKUNIT_EXPECT_TRUE(false);}

            std::cout << "tmp srk doing Shape smoothing for hex_2 case..." << std::endl;

            //bool do_jacobi = true;
            //Mesquite::MsqDebug::enable(1);
            //Mesquite::MsqDebug::enable(2);
            //Mesquite::MsqDebug::enable(3);

            int  msq_debug             = 2; // 1,2,3 for more debug info
            bool always_smooth         = true;
            int innerIter = 100;

            if (p_size == 1) 
              {
                PerceptMesquiteMesh pmm(&eMesh, 0, &boundarySelector);
                PerceptMesquiteMeshDomain pmd(&eMesh, 0);
                percept::PMMShapeImprover si(innerIter);
                si.run(pmm, pmd, always_smooth, msq_debug);
              }
            else
              {
                PerceptMesquiteMesh pmm0(&eMesh, 0, &boundarySelector);
                PerceptMesquiteMesh::PMMParallelMesh pmm(&pmm0);
                PerceptMesquiteMeshDomain pmd(&eMesh, 0);
                percept::PMMShapeImprover si(innerIter);
                si.run(pmm, pmd, always_smooth, msq_debug);
              }

            eMesh.save_as(output_files_loc+"hex_2_si_smooth.1.e");

          }
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      // A cube with an indented bump on the bottom, new parallel smoother

      STKUNIT_UNIT_TEST(unit_perceptMesquite, hex_4)
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
            unsigned n = 12;
            std::cout << "P["<<p_rank<<"] " << "tmp srk doing Laplace smoothing for hex_1 case, n = " << n << std::endl;
            //unsigned nn = n+1;
            //std::string gmesh_spec = toString(n)+"x"+toString(n)+"x"+toString(n*p_size)+std::string("|bbox:0,0,0,1,1,1|sideset:xXyYzZ");
            std::string gmesh_spec = toString(n)+"x"+toString(n)+"x"+toString(n)+std::string("|bbox:0,0,0,1,1,1|sideset:xXyYzZ");
            PerceptMesh eMesh(3);
            eMesh.new_mesh(percept::GMeshSpec(gmesh_spec));
            eMesh.addParallelInfoFields(true,true);
            eMesh.add_coordinate_state_fields();
            eMesh.commit();

            stk::mesh::Selector boundarySelector_1(*eMesh.get_non_const_part("surface_1") );
            stk::mesh::Selector boundarySelector_2(*eMesh.get_non_const_part("surface_2") );
            stk::mesh::Selector boundarySelector_3(*eMesh.get_non_const_part("surface_3") );
            stk::mesh::Selector boundarySelector_4(*eMesh.get_non_const_part("surface_4") );
            stk::mesh::Selector boundarySelector_5(*eMesh.get_non_const_part("surface_5") );
            stk::mesh::Selector boundarySelector_6(*eMesh.get_non_const_part("surface_6") );
            stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4 | boundarySelector_5 | boundarySelector_6;

            eMesh.populateParallelInfoFields(true,true,&boundarySelector);

            const std::vector<stk::mesh::Bucket*> & buckets = eMesh.get_bulk_data()->buckets( stk::mesh::fem::FEMMetaData::NODE_RANK );

            // cluster the mesh towards the bump
            for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
              {
                //if (boundarySelector_5(**k)) 
                  {
                    stk::mesh::Bucket & bucket = **k ;

                    const unsigned num_elements_in_bucket = bucket.size();
                
                    for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                      {
                        stk::mesh::Entity& entity = bucket[iEntity];

                        double * data = stk::mesh::field_data( *eMesh.get_coordinates_field() , entity );
                        data[2] = data[2]*data[2];
                      }
                  }
              }
            eMesh.save_as(input_files_loc+"hex_4_smooth.0.e");

            // save state of original mesh
            // field, dst, src: 
            eMesh.copy_field(eMesh.get_field("coordinates_NM1"), eMesh.get_coordinates_field());

            Mesquite::MsqError err;
            Mesquite::MeshImpl *msqMesh = create_mesquite_mesh(&eMesh, &boundarySelector);
            std::ostringstream vtk_file;
            vtk_file << "par_original_hex_mesh." << p_size << "." << p_rank << ".vtk";

            msqMesh->write_vtk(vtk_file.str().c_str(), err); 
            if (err) {std::cout << err << endl;  STKUNIT_EXPECT_TRUE(false);}

            for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
              {
                if (boundarySelector_5(**k)) 
                  {
                    stk::mesh::Bucket & bucket = **k ;

                    const unsigned num_elements_in_bucket = bucket.size();
                
                    for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                      {
                        stk::mesh::Entity& entity = bucket[iEntity];

                        double * data = stk::mesh::field_data( *eMesh.get_coordinates_field() , entity );
                        double ix = data[0];
                        double iy = data[1];
                        data[2] = (ix)*(1.0-ix)*(iy)*(1.0-iy)*2.0*4.;
                      }
                  }
              }
            // save state of projected mesh
            // field, dst, src: 
            eMesh.copy_field(eMesh.get_field("coordinates_N"), eMesh.get_coordinates_field());

            eMesh.save_as(input_files_loc+"hex_4_smooth.0_perturbed.e");

            msqMesh = create_mesquite_mesh(&eMesh, &boundarySelector);
            std::ostringstream vtk_file1;
            vtk_file1 << "par_untangle_original_hex_mesh." << p_size << "." << p_rank << ".vtk";

            msqMesh->write_vtk(vtk_file1.str().c_str(), err); 
            if (err) {std::cout << err << endl;  STKUNIT_EXPECT_TRUE(false);}

            std::cout << "tmp srk doing Shape smoothing for hex_4 case..." << std::endl;

            //bool do_jacobi = true;
            //Mesquite::MsqDebug::enable(1);
            Mesquite::MsqDebug::disable(1);
            //Mesquite::MsqDebug::enable(2);
            //Mesquite::MsqDebug::enable(3);

            int  msq_debug             = 0; // 1,2,3 for more debug info
            bool always_smooth         = true;
            int innerIter = 1001;

            {
              PerceptMesquiteMesh pmm(&eMesh, 0, &boundarySelector);
              percept::PMMParallelShapeImprover pmmpsi(innerIter, 1.e-4, 1);
              pmmpsi.run(pmm, 0, always_smooth, msq_debug);
            }

            eMesh.save_as(output_files_loc+"hex_4_si_smooth.1.e");

          }
      }


#define MSQ_TEST 0
#if MSQ_TEST

      //=============================================================================
      //=============================================================================
      //=============================================================================

      // A test of Mesquite read/write vtk, parallel smooth

      static void run_laplace(Mesquite::Mesh *mesh, Mesquite::ParallelMesh *pmesh, Mesquite::MsqError& err)
      {
        using namespace Mesquite;
        Settings settings;
        QualityAssessor qa;

        IdealWeightInverseMeanRatio qa_metric;
        qa.add_quality_assessment( &qa_metric );
  
        LaplacianSmoother smoother;
        TerminationCriterion outer("<type:laplace_outer>"), inner("<type:laplace_inner>");
        inner.add_iteration_limit( 1 );
        outer.add_iteration_limit( 100 );
        inner.add_absolute_vertex_movement_edge_length( 1.e-3 );
        outer.add_absolute_vertex_movement_edge_length( 1.e-3 );

        smoother.set_inner_termination_criterion( &inner );
        smoother.set_outer_termination_criterion( &outer );
  
        InstructionQueue q;
        q.add_quality_assessor( &qa, err ); MSQ_ERRRTN(err);
        q.set_master_quality_improver( &smoother, err ); MSQ_ERRRTN(err);
        q.add_quality_assessor( &qa, err ); MSQ_ERRRTN(err);
        q.run_common( mesh, pmesh, 0, &settings, err ); MSQ_ERRRTN(err);
      }

      STKUNIT_UNIT_TEST(unit_perceptMesquite, msq_hex_2)
      {

        using namespace Mesquite;

#define VTK_3D_DIR "./"

        using namespace std;
  
        /* init MPI */
        int rank, nprocs;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

        /* create processor-specific file names */
        ostringstream in_name, out_name, echo_name;
        if (nprocs > 1)
          {
            in_name << VTK_3D_DIR << "par_original_hex_mesh." << nprocs << "." << rank <<  ".vtk";
            out_name << "par_smoothed_hex_mesh." << nprocs << "." << rank <<  ".vtk";
            echo_name << "par_input_echo_hex_mesh." << nprocs << "." << rank <<  ".vtk";
          }
        else
          {
            in_name << VTK_3D_DIR << "par_original_hex_mesh.1.0.vtk" ;
            out_name << "par_smoothed_hex_mesh.1.0.vtk";
            echo_name << "par_input_echo_hex_mesh.1.0.vtk";
          }

        //out_name << VTK_3D_DIR << "par_smoothed_hex_mesh." << nprocs << "." << rank <<  ".vtk";

        /* load different mesh files on each processor */
        Mesquite::MsqError err;
        Mesquite::MeshImpl mesh;
        std::cout << "P[" << rank << "] tmp srk ::msq_hex_2 read_vtk..." << std::endl;
        mesh.read_vtk(in_name.str().c_str(), err);
        std::cout << "P[" << rank << "] tmp srk ::msq_hex_2 read_vtk...done" << std::endl;
        if (err) {cerr << err << endl;  STKUNIT_EXPECT_TRUE(false);}

        std::cout << "P[" << rank << "] tmp srk ::msq_hex_2 write_vtk echo..." << std::endl;
        mesh.write_vtk(echo_name.str().c_str(),err);
        if (err) {cerr << err << endl;  STKUNIT_EXPECT_TRUE(false);}
        std::cout << "P[" << rank << "] tmp srk ::msq_hex_2 write_vtk echo...done" << std::endl;

        /* create parallel mesh instance, specifying tags 
         * containing parallel data */
        if (nprocs > 1)
          {
            Mesquite::ParallelMeshImpl parallel_mesh(&mesh, "GLOBAL_ID", "PROCESSOR_ID");
            Mesquite::ParallelHelperImpl helper;
            helper.set_communicator(MPI_COMM_WORLD);
            helper.set_parallel_mesh(&parallel_mesh);
            parallel_mesh.set_parallel_helper(&helper);

            /* do Laplacian smooth */
            //LaplaceWrapper optimizer;
            //optimizer.run_instructions(&parallel_mesh, err);
            run_laplace(&mesh, &parallel_mesh, err);

            if (err) {cerr << err << endl;  STKUNIT_EXPECT_TRUE(false);}
          }
        else
          {
            //LaplaceWrapper optimizer;
            //optimizer.run_instructions(&mesh, err);
            run_laplace(&mesh, 0, err);

            if (err) {cerr << err << endl;  STKUNIT_EXPECT_TRUE(false);}
          }

        /* write mesh */
        mesh.write_vtk(out_name.str().c_str(),err);
        if (err) {cerr << err << endl;  STKUNIT_EXPECT_TRUE(false);}

        print_timing_diagnostics(cout);

      }
#endif      

#endif

    }



  }

}
#endif

