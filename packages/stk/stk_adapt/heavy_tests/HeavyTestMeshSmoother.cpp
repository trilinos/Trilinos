/*--------------------------------------------------------------------*/
/*    Copyright 2009, 2011 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/


#include <stk_percept/Percept.hpp>

#if defined( STK_PERCEPT_HAS_GEOMETRY )

#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/Util.hpp>
#include <stk_percept/ExceptionWatch.hpp>
#include <stk_percept/fixtures/Fixture.hpp>
#include <stk_percept/fixtures/Fixture.hpp>
#include <stk_percept/fixtures/QuadFixture.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>
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

#include <stk_percept/mesh/mod/smoother/MeshSmoother.hpp>
#include <stk_percept/mesh/mod/smoother/ReferenceMeshSmoother1.hpp>
#include <stk_percept/mesh/mod/smoother/SpacingFieldUtil.hpp>

#include <stk_adapt/UniformRefinerPattern.hpp>
#include <stk_adapt/UniformRefiner.hpp>

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

#define DO_TESTS 1
#if DO_TESTS

namespace stk 
{
  namespace adapt
  {
    namespace heavy_tests 
    {


#define EXTRA_PRINT 0
      static int s_par_size_max = 2;

#if 1
      static const std::string path_sep = "._.";
      static const std::string input_files_loc="./input_files"+path_sep;
      static const std::string output_files_loc="./output_files"+path_sep;
#else
      static const std::string input_files_loc="./input_files/";
      static const std::string output_files_loc="./output_files/";
#endif

#define EXTRA_PRINT 0

      static std::string procs_string[9] = {"np0", "np1", "np2", "np3", "np4", "np5", "np6", "np7", "np8"};

      //=============================================================================
      //=============================================================================
      //=============================================================================

      // run PMMParallelShapeImprover on a domain with one side perturbed
      STKUNIT_UNIT_TEST(heavy_perceptMeshSmoother, quad_4)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        MPI_Barrier( MPI_COMM_WORLD );
        unsigned par_size_max = s_par_size_max;
        
        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        //if (p_size == 1 || p_size == 2)
        if (1 || p_size <= par_size_max)
          {
            const unsigned nele = 12;
            //const unsigned nx = nele , ny = nele , nz = p_size*nele ;
            const unsigned nx = nele , ny = nele;

            bool sidesets=true;
            percept::QuadFixture<double> fixture( pm , nx , ny, sidesets);
            fixture.set_bounding_box(0,1,0,1);
            fixture.meta_data.commit();
            fixture.generate_mesh();

            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data);
            //eMesh.print_info("quad fixture",  2);
            eMesh.save_as(input_files_loc+"quad_4_smooth.0.e");

            eMesh.reopen();
            eMesh.add_coordinate_state_fields();
            eMesh.add_spacing_fields();
            stk::mesh::FieldBase *proc_rank_field = eMesh.add_field("proc_rank", stk::mesh::MetaData::ELEMENT_RANK, 0);
            eMesh.commit();
            eMesh.set_proc_rank_field(proc_rank_field);

            stk::mesh::Selector boundarySelector_1(*eMesh.get_non_const_part("surface_1") );
            stk::mesh::Selector boundarySelector_2(*eMesh.get_non_const_part("surface_2") );
            stk::mesh::Selector boundarySelector_3(*eMesh.get_non_const_part("surface_3") );
            stk::mesh::Selector boundarySelector_4(*eMesh.get_non_const_part("surface_4") );
            stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4;

            const std::vector<stk::mesh::Bucket*> & buckets = eMesh.get_bulk_data()->buckets( stk::mesh::MetaData::NODE_RANK );

            for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
              {
                {
                  stk::mesh::Bucket & bucket = **k ;

                  const unsigned num_elements_in_bucket = bucket.size();
                
                  for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                    {
                      stk::mesh::Entity node = bucket[iEntity];

                      double * data = stk::mesh::field_data( *eMesh.get_coordinates_field() , node );
                      double iy = data[1]; // /double(nele);
                      iy = iy*iy;
                      data[1] = iy; // *double(nele);
                    }
                }
              }
            SpacingFieldUtil sfu(eMesh);
            sfu.compute_spacing_field();
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
                        stk::mesh::Entity node = bucket[iEntity];

                        double * data = stk::mesh::field_data( *eMesh.get_coordinates_field() , node );
                        double ix = data[0]; // /double(nele);
                        //double bump_size=2.8; // 0.8
                        double bump_size=2.8; // 0.8
                        data[1] += (ix)*(1.0-ix)*bump_size; //*double(nele);
                        //std::cout << "tmp srk surface 1 node = " << data[0] << " " << data[1] << std::endl;
                      }
                  }
              }

            // save state of projected mesh
            // field, dst, src: 
            eMesh.copy_field(eMesh.get_field("coordinates_N"), eMesh.get_coordinates_field());

            eMesh.save_as(input_files_loc+"quad_4_smooth.0_perturbed.e");

            int  msq_debug             = 0; // 1,2,3 for more debug info
            bool always_smooth         = true;
            {
              percept::ReferenceMeshSmoother1 pmmpsi(&eMesh, &boundarySelector, 0, 1001, 1.e-4, 1);
              pmmpsi.run( always_smooth, msq_debug);
            }

            eMesh.save_as(output_files_loc+"quad_4_si_smooth.1.e");


          }
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      // run new parallel smoother on an interior-node perturbed domain
      STKUNIT_UNIT_TEST(heavy_perceptMeshSmoother, quad_4_1)
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
            //eMesh.addParallelInfoFields(true,true);
            eMesh.add_coordinate_state_fields();
            stk::mesh::FieldBase *proc_rank_field = eMesh.add_field("proc_rank", stk::mesh::MetaData::ELEMENT_RANK, 0);
            //eMesh.addParallelInfoFields(true,true);
            eMesh.commit();
            eMesh.set_proc_rank_field(proc_rank_field);

            stk::mesh::Selector boundarySelector_1(*eMesh.get_non_const_part("surface_1") );
            stk::mesh::Selector boundarySelector_2(*eMesh.get_non_const_part("surface_2") );
            stk::mesh::Selector boundarySelector_3(*eMesh.get_non_const_part("surface_3") );
            stk::mesh::Selector boundarySelector_4(*eMesh.get_non_const_part("surface_4") );
            stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4;

            //eMesh.populateParallelInfoFields(true,true,&boundarySelector);

            //eMesh.print_info("quad fixture",  2);
            eMesh.save_as(input_files_loc+"quad_4_1_smooth.0.e");

            eMesh.copy_field(eMesh.get_field("coordinates_NM1"), eMesh.get_coordinates_field());

            unsigned center_node_id = 5;
            stk::mesh::Entity node = eMesh.get_bulk_data()->get_entity(0, center_node_id);
            double *data = 0;
            if (node.is_valid())
              {
                data = stk::mesh::field_data( *eMesh.get_coordinates_field() , node );
                //std::cout << "tmp srk  center node= " << data[0] << " " << data[1] << std::endl;
                if (1)
                  {
                    data[0] += .2;
                    data[1] += .3;
                  }
              }

            eMesh.save_as(input_files_loc+"quad_4_1_smooth.0_perturbed.e");
            eMesh.copy_field(eMesh.get_field("coordinates_N"), eMesh.get_coordinates_field());

            {
              percept::ReferenceMeshSmoother1 pmmpsi(&eMesh, &boundarySelector, 0, 1001, 1.e-4, 1);
              pmmpsi.run(true, 0);
            }

            if (node.is_valid())
              {
                STKUNIT_EXPECT_DOUBLE_EQ_APPROX(data[0], 1.0);
                STKUNIT_EXPECT_DOUBLE_EQ_APPROX(data[1], 1.0);
              }

            eMesh.save_as(output_files_loc+"quad_4_1_smooth.1.e");
          }
      }


      //=============================================================================
      //=============================================================================
      //=============================================================================

      // run PMMParallelShapeImprover on a domain with one side perturbed
      STKUNIT_UNIT_TEST(heavy_perceptMeshSmoother, tri_4)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        MPI_Barrier( MPI_COMM_WORLD );
        unsigned par_size_max = s_par_size_max;
        
        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        //if (p_size == 1 || p_size == 2)
        if (1 || p_size <= par_size_max)
          {
            const unsigned nele = 12;
            //const unsigned nx = nele , ny = nele , nz = p_size*nele ;
            const unsigned nx = nele , ny = nele;

            bool sidesets=true;
            percept::QuadFixture<double, Triangle<3> > fixture( pm , nx , ny, sidesets);
            fixture.set_bounding_box(0,1,0,1);
            fixture.meta_data.commit();
            fixture.generate_mesh();

            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data);
            //eMesh.print_info("tri fixture",  2);
            eMesh.save_as(input_files_loc+"tri_4_smooth.0.e");

            eMesh.reopen();
            eMesh.add_coordinate_state_fields();
            eMesh.add_spacing_fields();
            stk::mesh::FieldBase *proc_rank_field = eMesh.add_field("proc_rank", stk::mesh::MetaData::ELEMENT_RANK, 0);
            eMesh.commit();
            eMesh.set_proc_rank_field(proc_rank_field);

            stk::mesh::Selector boundarySelector_1(*eMesh.get_non_const_part("surface_1") );
            stk::mesh::Selector boundarySelector_2(*eMesh.get_non_const_part("surface_2") );
            stk::mesh::Selector boundarySelector_3(*eMesh.get_non_const_part("surface_3") );
            stk::mesh::Selector boundarySelector_4(*eMesh.get_non_const_part("surface_4") );
            stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4;

            const std::vector<stk::mesh::Bucket*> & buckets = eMesh.get_bulk_data()->buckets( stk::mesh::MetaData::NODE_RANK );

            for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
              {
                {
                  stk::mesh::Bucket & bucket = **k ;

                  const unsigned num_elements_in_bucket = bucket.size();
                
                  for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                    {
                      stk::mesh::Entity node = bucket[iEntity];

                      double * data = stk::mesh::field_data( *eMesh.get_coordinates_field() , node );
                      double iy = data[1]; // /double(nele);
                      iy = iy*iy;
                      data[1] = iy; // *double(nele);
                    }
                }
              }
            SpacingFieldUtil sfu(eMesh);
            sfu.compute_spacing_field();
            eMesh.save_as(input_files_loc+"tri_4_smooth.1.e");

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
                        stk::mesh::Entity node = bucket[iEntity];

                        double * data = stk::mesh::field_data( *eMesh.get_coordinates_field() , node );
                        double ix = data[0]; // /double(nele);
                        //double bump_size=2.8; // 0.8
                        double bump_size=2.8; // 0.8
                        data[1] += (ix)*(1.0-ix)*bump_size; //*double(nele);
                        //std::cout << "tmp srk surface 1 node = " << data[0] << " " << data[1] << std::endl;
                      }
                  }
              }

            // save state of projected mesh
            // field, dst, src: 
            eMesh.copy_field(eMesh.get_field("coordinates_N"), eMesh.get_coordinates_field());

            eMesh.save_as(input_files_loc+"tri_4_smooth.0_perturbed.e");

            int  msq_debug             = 0; // 1,2,3 for more debug info
            bool always_smooth         = true;
            {
              percept::ReferenceMeshSmoother1 pmmpsi(&eMesh, &boundarySelector, 0, 1001, 1.e-4, 1);
              pmmpsi.run( always_smooth, msq_debug);
            }

            eMesh.save_as(output_files_loc+"tri_4_si_smooth.1.e");


          }
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      // A cube with an indented bump on the bottom, new parallel smoother

      STKUNIT_UNIT_TEST(heavy_perceptMeshSmoother, hex_4)
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

            const std::vector<stk::mesh::Bucket*> & buckets = eMesh.get_bulk_data()->buckets( stk::mesh::MetaData::NODE_RANK );

            // cluster the mesh towards the bump
            for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
              {
                //if (boundarySelector_5(**k)) 
                {
                  stk::mesh::Bucket & bucket = **k ;

                  const unsigned num_elements_in_bucket = bucket.size();
                
                  for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                    {
                      stk::mesh::Entity entity = bucket[iEntity];

                      double * data = stk::mesh::field_data( *eMesh.get_coordinates_field() , entity );
                      data[2] = data[2]*data[2];
                    }
                }
              }
            eMesh.save_as(input_files_loc+"hex_4_smooth.0.e");

            // save state of original mesh
            // field, dst, src: 
            eMesh.copy_field(eMesh.get_field("coordinates_NM1"), eMesh.get_coordinates_field());

            for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
              {
                if (boundarySelector_5(**k)) 
                  {
                    stk::mesh::Bucket & bucket = **k ;

                    const unsigned num_elements_in_bucket = bucket.size();
                
                    for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                      {
                        stk::mesh::Entity entity = bucket[iEntity];

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

            std::cout << "tmp srk doing Shape smoothing for hex_4 case..." << std::endl;

            //bool do_jacobi = true;

            int  msq_debug             = 0; // 1,2,3 for more debug info
            bool always_smooth         = true;
            int innerIter = 1001;

            {
              percept::ReferenceMeshSmoother1 pmmpsi(&eMesh, &boundarySelector, 0, innerIter, 1.e-4, 1);
              pmmpsi.run( always_smooth, msq_debug);
            }

            eMesh.save_as(output_files_loc+"hex_4_si_smooth.1.e");

          }
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      // A cube with an indented bump on the bottom, new parallel smoother, convert to tet

      // heavy_perceptMeshSmoother.tet_4
      STKUNIT_UNIT_TEST(heavy_perceptMeshSmoother, tet_4)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        MPI_Barrier( MPI_COMM_WORLD );

        const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        //if (p_size == 1 || p_size == 3)
        if (p_size <= 2)
          {
            unsigned n = 12;
            std::cout << "P["<<p_rank<<"] " << "tmp srk doing Laplace smoothing for tet_1 case, n = " << n << std::endl;
            std::string gmesh_spec = toString(n)+"x"+toString(n)+"x"+toString(n)+std::string("|bbox:0,0,0,1,1,1|sideset:xXyYzZ");
            PerceptMesh eMesh(3);
            eMesh.new_mesh(percept::GMeshSpec(gmesh_spec));
            eMesh.addParallelInfoFields(true,true);
            eMesh.add_coordinate_state_fields();
            int scalarDimension=0;
            stk::mesh::FieldBase* proc_rank_field = eMesh.add_field("proc_rank", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);
            Hex8_Tet4_6_12 convert_pattern(eMesh);
            eMesh.commit();

            UniformRefiner converter(eMesh, convert_pattern, proc_rank_field);
            converter.doBreak();
            eMesh.save_as(input_files_loc+"tet_4_smooth.0.0.e");

            stk::mesh::Selector boundarySelector_1(*eMesh.get_non_const_part("surface_1") );
            stk::mesh::Selector boundarySelector_2(*eMesh.get_non_const_part("surface_2") );
            stk::mesh::Selector boundarySelector_3(*eMesh.get_non_const_part("surface_3") );
            stk::mesh::Selector boundarySelector_4(*eMesh.get_non_const_part("surface_4") );
            stk::mesh::Selector boundarySelector_5(*eMesh.get_non_const_part("surface_5") );
            stk::mesh::Selector boundarySelector_6(*eMesh.get_non_const_part("surface_6") );
            stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4 | boundarySelector_5 | boundarySelector_6;

            eMesh.populateParallelInfoFields(true,true,&boundarySelector);

            const std::vector<stk::mesh::Bucket*> & buckets = eMesh.get_bulk_data()->buckets( stk::mesh::MetaData::NODE_RANK );

            // cluster the mesh towards the bump
            for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
              {
                //if (boundarySelector_5(**k)) 
                  {
                    stk::mesh::Bucket & bucket = **k ;

                    const unsigned num_elements_in_bucket = bucket.size();
                
                    for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                      {
                        stk::mesh::Entity entity = bucket[iEntity];

                        double * data = stk::mesh::field_data( *eMesh.get_coordinates_field() , entity );
                        data[2] = data[2]*data[2];
                      }
                  }
              }
            eMesh.save_as(input_files_loc+"tet_4_smooth.0.e");

            // save state of original mesh
            // field, dst, src: 
            eMesh.copy_field(eMesh.get_field("coordinates_NM1"), eMesh.get_coordinates_field());

            for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
              {
                if (boundarySelector_5(**k)) 
                  {
                    stk::mesh::Bucket & bucket = **k ;

                    const unsigned num_elements_in_bucket = bucket.size();
                
                    for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                      {
                        stk::mesh::Entity entity = bucket[iEntity];

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

            eMesh.save_as(input_files_loc+"tet_4_smooth.0_perturbed.e");

            std::cout << "tmp srk doing Shape smoothing for tet_4 case..." << std::endl;

            //bool do_jacobi = true;

            int  msq_debug             = 0; // 1,2,3 for more debug info
            bool always_smooth         = true;
            int innerIter = 1001;

            {
              stk::percept::ReferenceMeshSmoother1 pmmpsi(&eMesh, &boundarySelector, 0, innerIter, 1.e-4, 1);
              pmmpsi.run( always_smooth, msq_debug);
            }

            eMesh.save_as(output_files_loc+"tet_4_si_smooth.1.e");

          }

      }




    }

  }
}
#endif
#endif

