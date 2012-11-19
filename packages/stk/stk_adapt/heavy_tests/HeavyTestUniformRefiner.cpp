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
#include <stk_adapt/UniformRefiner.hpp>
#include <unit_tests/TestLocalRefiner.hpp>
#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>
#include <stk_percept/RunEnvironment.hpp>

#include <stk_percept/fixtures/Fixture.hpp>
#include <stk_percept/fixtures/BeamFixture.hpp>
#include <stk_percept/fixtures/HeterogeneousFixture.hpp>
#include <stk_percept/fixtures/PyramidFixture.hpp>
#include <stk_percept/fixtures/QuadFixture.hpp>
#include <stk_percept/fixtures/WedgeFixture.hpp>

// smoothing tests
#include <stk_percept/mesh/mod/smoother/MeshSmoother.hpp>
#include <stk_percept/mesh/mod/smoother/ReferenceMeshSmoother1.hpp>


// this is for testing the local-refine refactoring
#define UNIFORM_REFINER UniformRefiner
//#define UNIFORM_REFINER TestLocalRefiner


namespace stk
{
  namespace adapt
  {
    namespace heavy_tests
    {
      //using namespace regression_tests;

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


      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================
      //========= AREA for tests in progress of being debugged
      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================


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

            UNIFORM_REFINER converter(eMesh, convert_pattern, proc_rank_field);
            converter.doBreak();
            eMesh.save_as(input_files_loc+"tet_4_smooth.0.0.e");

            stk::mesh::Selector boundarySelector_1(*eMesh.get_non_const_part("surface_1#Triangle_3#") );
            stk::mesh::Selector boundarySelector_2(*eMesh.get_non_const_part("surface_2#Triangle_3#") );
            stk::mesh::Selector boundarySelector_3(*eMesh.get_non_const_part("surface_3#Triangle_3#") );
            stk::mesh::Selector boundarySelector_4(*eMesh.get_non_const_part("surface_4#Triangle_3#") );
            stk::mesh::Selector boundarySelector_5(*eMesh.get_non_const_part("surface_5#Triangle_3#") );
            stk::mesh::Selector boundarySelector_6(*eMesh.get_non_const_part("surface_6#Triangle_3#") );
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


    }//    namespace heavy_tests
  }//  namespace adapt
}// namespace stk

