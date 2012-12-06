/*--------------------------------------------------------------------*/
/*    Copyright 2009, 2011 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stk_percept/function/StringFunction.hpp>
#include <stk_percept/function/FieldFunction.hpp>
#include <stk_percept/function/ConstantFunction.hpp>
#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/Util.hpp>
#include <stk_percept/ExceptionWatch.hpp>
#include <stk_percept/fixtures/Fixture.hpp>
#include <stk_percept/fixtures/QuadFixture.hpp>
#include <stk_percept/mesh/mod/smoother/JacobianUtil.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <Teuchos_ScalarTraits.hpp>

#include <stk_percept/fixtures/BeamFixture.hpp>
#include <stk_percept/fixtures/HeterogeneousFixture.hpp>
#include <stk_percept/fixtures/QuadFixture.hpp>
#include <stk_percept/fixtures/WedgeFixture.hpp>

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

      //=============================================================================
      //=============================================================================
      //=============================================================================

      /// test stretch_eigens

      STKUNIT_UNIT_TEST(perceptMesh, jacobian_stretch_eigens)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        MPI_Barrier( MPI_COMM_WORLD );

        const unsigned p_size = stk::parallel_machine_size( pm );
        const unsigned p_rank = stk::parallel_machine_rank( pm );
        if (p_size <= 2)
          {
            // create a 12x12 quad mesh with sidesets
            const unsigned n = 12;
            //const unsigned nx = n , ny = n , nz = p_size*n ;
            const unsigned nx = n , ny = n;

            bool sidesets_on = true;
            percept::QuadFixture<double> fixture( pm , nx , ny, sidesets_on);
            fixture.meta_data.commit();
            fixture.generate_mesh();

            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data);
            eMesh.print_info("quad fixture",  0);
            //eMesh.save_as("./output_files/quad_fixture.e");

            stk::mesh::BulkData& bulkData = *eMesh.get_bulk_data();
            VectorFieldType* coordField = eMesh.get_coordinates_field();

            const std::vector<stk::mesh::Bucket*> & buckets = bulkData.buckets( eMesh.element_rank() );

            for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
              {
                //if (in_surface_selector(**k))
                  {
                    stk::mesh::Bucket & bucket = **k ;

                    // in case the cell topology is needed
                    const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(bucket);
                    shards::CellTopology cell_topo(cell_topo_data);

                    const unsigned num_elements_in_bucket = bucket.size();

                    for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                      {
                        stk::mesh::Entity element = bucket[iElement];
                        double eigens[3];

                        JacobianUtil jac;
                        jac.stretch_eigens(eMesh, element, eigens, coordField, cell_topo_data);
                        
                        const double tol=1.e-6;
                        STKUNIT_EXPECT_DOUBLE_EQ_APPROX_TOL(eigens[0], 1.0, tol);
                        STKUNIT_EXPECT_DOUBLE_EQ_APPROX_TOL(eigens[1], 1.0, tol);
                        STKUNIT_EXPECT_DOUBLE_EQ_APPROX_TOL(eigens[2], 1.0, tol);
                        if (0) std::cout << "P[" << p_rank << ":" << p_size << "] element= " << element.identifier() << " eigens = " << eigens[0] << " " << eigens[1] << " " << eigens[2] << std::endl;
                      }
                  }
              }

          }

        if (p_size <= 2)
          {
            // create a 12x12 quad mesh with sidesets
            const unsigned n = 12;
            //const unsigned nx = n , ny = n , nz = p_size*n ;
            const unsigned nx = n , ny = n;

            bool sidesets_on = true;
            percept::QuadFixture<double> fixture( pm , nx , ny, sidesets_on);
            fixture.meta_data.commit();
            fixture.set_bounding_box(0,1,0,2);
            fixture.generate_mesh();

            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data);
            eMesh.print_info("quad fixture",  0);
            eMesh.save_as("jac-test-1.e");

            double min_max_ave[3]={0,0,0};
            double hmesh = eMesh.hmesh_stretch_eigens(min_max_ave);
            const double tol=1.e-6;
            STKUNIT_EXPECT_DOUBLE_EQ_APPROX_TOL(hmesh, 2.0/12.0, tol);
            STKUNIT_EXPECT_DOUBLE_EQ_APPROX_TOL(min_max_ave[0], 2.0/12.0, tol);
            STKUNIT_EXPECT_DOUBLE_EQ_APPROX_TOL(min_max_ave[1], 2.0/12.0, tol);
            STKUNIT_EXPECT_DOUBLE_EQ_APPROX_TOL(min_max_ave[2], 2.0/12.0, tol);

            stk::mesh::BulkData& bulkData = *eMesh.get_bulk_data();
            VectorFieldType* coordField = eMesh.get_coordinates_field();

            const std::vector<stk::mesh::Bucket*> & buckets = bulkData.buckets( eMesh.element_rank() );

            for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
              {
                //if (in_surface_selector(**k))
                  {
                    stk::mesh::Bucket & bucket = **k ;

                    // in case the cell topology is needed
                    const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(bucket);
                    shards::CellTopology cell_topo(cell_topo_data);

                    const unsigned num_elements_in_bucket = bucket.size();

                    for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                      {
                        stk::mesh::Entity element = bucket[iElement];
                        double eigens[3];

                        JacobianUtil jac;
                        jac.stretch_eigens(eMesh, element, eigens, coordField, cell_topo_data);
                        
                        STKUNIT_EXPECT_DOUBLE_EQ_APPROX_TOL(eigens[0], 1.0, tol);
                        STKUNIT_EXPECT_DOUBLE_EQ_APPROX_TOL(eigens[1], 2.0/12.0, tol);
                        STKUNIT_EXPECT_DOUBLE_EQ_APPROX_TOL(eigens[2], 1.0/12.0, tol);
                        if (0) std::cout << "P[" << p_rank << ":" << p_size << "] element= " << element.identifier() << " eigens = " << eigens[0] << " " << eigens[1] << " " << eigens[2] << std::endl;
                      }
                  }
              }

          }

        if (0 && p_size <= 2)
          {
            // create a 12x12 quad mesh with sidesets
            const unsigned n = 12;
            //const unsigned nx = n , ny = n , nz = p_size*n ;
            const unsigned nx = n , ny = n;

            bool sidesets_on = true;
            percept::QuadFixture<double, shards::Triangle<3> > fixture( pm , nx , ny, sidesets_on);
            fixture.meta_data.commit();
            fixture.set_bounding_box(0,1,0,2);
            fixture.generate_mesh();

            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data);
            eMesh.print_info("tri fixture",  0);
            eMesh.save_as("jac-test-tri.e");

            double min_max_ave[3]={0,0,0};
            double hmesh = eMesh.hmesh_stretch_eigens(min_max_ave);

            const double tol=1.e-6;
            STKUNIT_EXPECT_DOUBLE_EQ_APPROX_TOL(hmesh, 2.0/12.0, tol);
            STKUNIT_EXPECT_DOUBLE_EQ_APPROX_TOL(min_max_ave[0], 2.0/12.0, tol);
            STKUNIT_EXPECT_DOUBLE_EQ_APPROX_TOL(min_max_ave[1], 2.0/12.0, tol);
            STKUNIT_EXPECT_DOUBLE_EQ_APPROX_TOL(min_max_ave[2], 2.0/12.0, tol);

            stk::mesh::BulkData& bulkData = *eMesh.get_bulk_data();
            VectorFieldType* coordField = eMesh.get_coordinates_field();

            const std::vector<stk::mesh::Bucket*> & buckets = bulkData.buckets( eMesh.element_rank() );

            for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
              {
                //if (in_surface_selector(**k))
                  {
                    stk::mesh::Bucket & bucket = **k ;

                    // in case the cell topology is needed
                    const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(bucket);
                    shards::CellTopology cell_topo(cell_topo_data);

                    const unsigned num_elements_in_bucket = bucket.size();

                    for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                      {
                        stk::mesh::Entity element = bucket[iElement];
                        double eigens[3];

                        JacobianUtil jac;
                        jac.stretch_eigens(eMesh, element, eigens, coordField, cell_topo_data);
                        
//                         STKUNIT_EXPECT_DOUBLE_EQ_APPROX_TOL(eigens[0], 1.0, tol);
//                         STKUNIT_EXPECT_DOUBLE_EQ_APPROX_TOL(eigens[1], 2.0/12.0, tol);
//                         STKUNIT_EXPECT_DOUBLE_EQ_APPROX_TOL(eigens[2], 1.0/12.0, tol);
                        if (1) std::cout << "P[" << p_rank << ":" << p_size << "] element= " << element.identifier() << " eigens = " << eigens[0] << " " << eigens[1] << " " << eigens[2] << std::endl;
                      }
                  }
              }

          }


      }



    }
  }
}
