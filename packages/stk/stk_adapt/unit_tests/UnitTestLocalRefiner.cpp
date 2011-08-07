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
#include <stk_percept/GeometryVerifier.hpp>
#include <stk_percept/function/StringFunction.hpp>
#include <stk_percept/function/FieldFunction.hpp>
#include <stk_percept/function/ConstantFunction.hpp>
#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/MeshUtil.hpp>

#include <stk_adapt/UniformRefinerPattern.hpp>
#include <stk_adapt/UniformRefiner.hpp>

#include <unit_tests/TestLocalRefinerTri.hpp>
#include <unit_tests/TestLocalRefinerTri1.hpp>
#include <unit_tests/TestLocalRefinerTri2.hpp>
#include <unit_tests/TestLocalRefinerTri_N.hpp>
#include <unit_tests/TestLocalRefinerTri_N_1.hpp>
#include <unit_tests/TestLocalRefinerTri_N_2.hpp>
#include <unit_tests/TestLocalRefinerTri_N_3.hpp>

#include <unit_tests/TestLocalRefinerTri_N_3_EdgeMarker.hpp>
#include <unit_tests/TestLocalRefinerTri_N_3_ElementMarker.hpp>

#include <stk_adapt/FieldBasedMarkerPredicate.hpp>
#include <stk_adapt/PredicateBasedMarker.hpp>

#include <unit_tests/TestLocalRefinerTet_N_1.hpp>
#include <unit_tests/TestLocalRefinerTet_N_2.hpp>
#include <unit_tests/TestLocalRefinerTet_N_2_1.hpp>
#include <unit_tests/TestLocalRefinerTet_N_3.hpp>
#include <unit_tests/TestLocalRefinerTet_N_3_1.hpp>
#include <unit_tests/TestLocalRefinerTet_N_4.hpp>


#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <unit_tests/UnitTestSupport.hpp>

#include <stk_io/IossBridge.hpp>

#include <boost/lexical_cast.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/util/PrintTable.hpp>

#include <stk_percept/Percept.hpp>
#include <stk_percept/Util.hpp>
#include <stk_percept/ExceptionWatch.hpp>

#include <stk_percept/function/ElementOp.hpp>

#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>
#include <stk_percept/RunEnvironment.hpp>

#include <stk_percept/fixtures/Fixture.hpp>
#include <stk_percept/fixtures/SingleTetFixture.hpp>
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
#include <cstdlib>

#include <math.h>
#include <stk_util/parallel/Parallel.hpp>

namespace stk {
  namespace adapt {
    namespace unit_tests {

      static int printInfoLevel = 0;
      static std::string post_fix[4] = {"np0", "np1", "np2", "np3"};

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

      static double  tet_volume(SingleTetFixture::Point *node_coord_data, SingleTetFixture::TetIds& tetra_node_ids, unsigned node_id_offset=0)
      {
        double mat[3][3];
        for (int inode=0; inode < 3; inode++)
          {
            for (int ipt=0; ipt < 3; ipt++)
              {
                mat[inode][ipt] = node_coord_data[tetra_node_ids[inode+1]-node_id_offset][ipt] - node_coord_data[tetra_node_ids[0]-node_id_offset][ipt];
              }
          }
        double vol = (mat[0][0]*(mat[1][1]*mat[2][2]-mat[1][2]*mat[2][1])
                      -mat[0][1]*(mat[1][0]*mat[2][2] - mat[1][2]*mat[2][0])
                      +mat[0][2]*(mat[1][0]*mat[2][1] - mat[1][1]*mat[2][0])
                      ) / 6.0;
        return vol;
      }

      static double totalVolume(PerceptMesh& eMesh)
      {
        double totVol=0.0;

        SingleTetFixture::Point node_coord_data[4] = {{0,0,0},{0,0,0},{0,0,0},{0,0,0}};
        static  SingleTetFixture::TetIds tetra_node_ids[] = { {0, 1, 2, 3} };

        const vector<stk::mesh::Bucket*> & buckets = eMesh.getBulkData()->buckets( eMesh.element_rank() );

        for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
          {
            stk::mesh::Bucket & bucket = **k ;

            const unsigned num_elements_in_bucket = bucket.size();
            for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
              {
                stk::mesh::Entity& element = bucket[iElement];
                stk::mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);
                for (unsigned inode=0; inode < elem_nodes.size(); inode++)
                  {
                    stk::mesh::Entity *node = elem_nodes[inode].entity();
                    double *fdata = stk::mesh::field_data( *eMesh.getCoordinatesField() , *node);
                    node_coord_data[inode][0] = fdata[0];
                    node_coord_data[inode][1] = fdata[1];
                    node_coord_data[inode][2] = fdata[2];
                  }
                double vol = tet_volume(node_coord_data, tetra_node_ids[0]);
                if (vol < 0)
                  {
                    std::cout << "ERROR neg vol = " << vol << std::endl;
                  }
                totVol += vol;
              }
          }        
        return totVol;
      }


      static void fixture_setup_NxNxN_box_hex_and_tet_mesh()
      {
        EXCEPTWATCH;
        static int entered=0;
        if (!entered)
          entered = 1;
        else
          return;

        MPI_Barrier( MPI_COMM_WORLD );

        int N=4;

        // start_demo_uniformRefiner_hex8_build
        {
          percept::PerceptMesh eMesh(3u);

          //unsigned p_size = eMesh.getParallelSize();

          // generate a N x N x N mesh
          std::string gmesh_spec = 
            toString(N)+"x"+
            toString(N)+"x"+
            toString(N)+
            std::string("|bbox:0,0,0,")+
            toString(N)+","+
            toString(N)+","+
            toString(N);
            
          eMesh.newMesh(percept::PerceptMesh::GMeshSpec(gmesh_spec));
          eMesh.commit();

          eMesh.saveAs(input_files_loc+"hex_fixture_NxNxN.e");

          // end_demo
        }

        // start_demo_uniformRefiner_hex8_build_1
        {
          percept::PerceptMesh eMesh(3u);

          //unsigned p_size = eMesh.getParallelSize();
          eMesh.open(input_files_loc+"hex_fixture_NxNxN.e");

          Hex8_Tet4_24 break_hex_to_tet(eMesh);

          int scalarDimension = 0; // a scalar
          stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);

          eMesh.commit();

          UniformRefiner breaker(eMesh, break_hex_to_tet, proc_rank_field);
          breaker.doBreak();
          save_or_diff(eMesh, input_files_loc+"tet_fixture_NxNxN.e");


          if (0)
          {
            PerceptMesh em1;
            em1.openReadOnly(input_files_loc+"tet_fixture_NxNxN_tmp.e");
            em1.saveAs(input_files_loc+"tet_fixture_NxNxN.e");
            em1.printInfo("srk tmp tet_fixture_NxNxN.e after reopen", 2);
          }
          // end_demo
        }
      }


      //=============================================================================
      //=============================================================================
      //=============================================================================
      /// check triangulate_tet

      STKUNIT_UNIT_TEST(unit_localRefiner, triangulate_tet)
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size(pm);

        if (p_size <= 1)
          {
            {
              // create the mesh

              stk::percept::SingleTetFixture mesh(pm, false);
              stk::io::put_io_part_attribute(  mesh.m_block_tet );
              mesh.m_metaData.commit();
              mesh.populate();

              bool isCommitted = true;
              percept::PerceptMesh eMesh(&mesh.m_metaData, &mesh.m_bulkData, isCommitted);
              eMesh.saveAs(input_files_loc+"local_tet_0.e");
            }

            {
              PerceptMesh eMesh;
              eMesh.open(input_files_loc+"local_tet_0.e");
              Local_Tet4_Tet4_N break_tet(eMesh);
              eMesh.commit();

              unsigned edge_mark_bitcode = 4u;
              TestLocalRefinerTet_N_2_1 breaker(eMesh, break_tet, 0, edge_mark_bitcode);
              breaker.setRemoveOldElements(true);
              breaker.doBreak();

              save_or_diff(eMesh, output_files_loc+"local_tet_N_2_1_bitcode_4.e");
            }

            {
              {

                // create the mesh
                unsigned npts=4;
                unsigned ntets=1;
                static  SingleTetFixture::Point node_coord_data[  ] = {
                  { 10 , 2 , 0 } , { 11 , 2 , 0 } , { 10 , 3 , 0 } , { 10 , 2 , 1 } };

                // Hard coded tetra node ids for all the tetra nodes in the entire mesh
                static  SingleTetFixture::TetIds tetra_node_ids[] = {
                  { 1, 2, 3, 4} };

                stk::percept::SingleTetFixture mesh(pm, false, npts, node_coord_data, ntets, tetra_node_ids);
                stk::io::put_io_part_attribute(  mesh.m_block_tet );
                mesh.m_metaData.commit();
                mesh.populate();

                std::cout << "here" << std::endl;
                bool isCommitted = true;
                percept::PerceptMesh eMesh(&mesh.m_metaData, &mesh.m_bulkData, isCommitted);
                eMesh.saveAs(input_files_loc+"local_tet_41.e");
              }

              {
                for (unsigned edge_mark_bitcode = 12u; edge_mark_bitcode <= 14u; edge_mark_bitcode++)
                  {
                    PerceptMesh eMesh;
                    eMesh.open(input_files_loc+"local_tet_41.e");
                    Local_Tet4_Tet4_N break_tet(eMesh);
                    eMesh.commit();

                    double totalVol0 = totalVolume(eMesh);
                    std::cout << "tmp edge_mark_bitcode= " << edge_mark_bitcode << " totalVol0 = " << totalVol0 << std::endl;

                    TestLocalRefinerTet_N_2_1 breaker(eMesh, break_tet, 0, edge_mark_bitcode);
                    breaker.setRemoveOldElements(true);
                    breaker.doBreak();

                    double totalVol1 = totalVolume(eMesh);
                    std::cout << "tmp edge_mark_bitcode= " << edge_mark_bitcode << " totalVol1 = " << totalVol1 << std::endl;

                    eMesh.saveAs(output_files_loc+"local_tet_N_2_1_bitcode_"+toString((int)edge_mark_bitcode)+".e" );
                  }
              }
            }

            {
              PerceptMesh eMesh;
              eMesh.open(input_files_loc+"local_tet_0.e");
              Local_Tet4_Tet4_N break_tet(eMesh);
              eMesh.commit();

              TestLocalRefinerTet_N_1 breaker(eMesh, break_tet, 0);
              breaker.setRemoveOldElements(false);
              breaker.doBreak();

              save_or_diff(eMesh, output_files_loc+"local_tet_N_1_1.e");
            }


            {
              PerceptMesh eMesh;
              eMesh.open(input_files_loc+"local_tet_0.e");
              Local_Tet4_Tet4_N break_tet(eMesh);
              eMesh.commit();

              TestLocalRefinerTet_N_2 breaker(eMesh, break_tet, 0);
              breaker.setRemoveOldElements(true);
              breaker.doBreak();

              save_or_diff(eMesh, output_files_loc+"local_tet_N_2_1.e");
            }

            {
              PerceptMesh eMesh;
              eMesh.open(input_files_loc+"local_tet_0.e");
              Local_Tet4_Tet4_N break_tet(eMesh);
              eMesh.commit();

              TestLocalRefinerTet_N_2 breaker(eMesh, break_tet, 0, 2);
              breaker.setRemoveOldElements(true);
              breaker.doBreak();

              save_or_diff(eMesh, output_files_loc+"local_tet_N_2_2.e");
            }

            {
              PerceptMesh eMesh;
              eMesh.open(input_files_loc+"local_tet_0.e");
              Local_Tet4_Tet4_N break_tet(eMesh);
              eMesh.commit();

              TestLocalRefinerTet_N_2 breaker(eMesh, break_tet, 0, 3);
              breaker.setRemoveOldElements(true);
              breaker.doBreak();

              save_or_diff(eMesh, output_files_loc+"local_tet_N_2_3.e");
            }

            {
              PerceptMesh eMesh;
              eMesh.open(input_files_loc+"local_tet_0.e");
              Local_Tet4_Tet4_N break_tet(eMesh);
              eMesh.commit();

              TestLocalRefinerTet_N_2 breaker(eMesh, break_tet, 0, 4);
              breaker.setRemoveOldElements(true);
              breaker.doBreak();

              save_or_diff(eMesh, output_files_loc+"local_tet_N_2_4.e");
            }

            {
              PerceptMesh eMesh;
              eMesh.open(input_files_loc+"local_tet_0.e");
              Local_Tet4_Tet4_N break_tet(eMesh);
              eMesh.commit();

              TestLocalRefinerTet_N_2 breaker(eMesh, break_tet, 0, 5);
              breaker.setRemoveOldElements(true);
              breaker.doBreak();

              save_or_diff(eMesh, output_files_loc+"local_tet_N_2_5.e");
            }

            {
              PerceptMesh eMesh;
              eMesh.open(input_files_loc+"local_tet_0.e");
              Local_Tet4_Tet4_N break_tet(eMesh);
              eMesh.commit();

              TestLocalRefinerTet_N_2 breaker(eMesh, break_tet, 0, 6);
              breaker.setRemoveOldElements(true);
              breaker.doBreak();

              save_or_diff(eMesh, output_files_loc+"local_tet_N_2_6.e");
            }
          }
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================
      /// check triangulate_tet - two tets sharing a face

      STKUNIT_UNIT_TEST(unit_localRefiner, triangulate_tet_2)
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size(pm);

        if (p_size <= 1)
          {
            {

              // create the mesh
              unsigned npts=5;
              unsigned ntets=2;
              static  SingleTetFixture::Point node_coord_data[  ] = {
                { 0 , 0 , 0 } , { 1 , 0 , 0 } , { 0 , 1 , 0 } , { 0 , 0 , 1 }, {1, 1, 1} };

              // Hard coded tetra node ids for all the tetra nodes in the entire mesh
              static  SingleTetFixture::TetIds tetra_node_ids[] = {
                { 1, 2, 3, 4}, {2, 3, 4, 5} };

              stk::percept::SingleTetFixture mesh(pm, false, npts, node_coord_data, ntets, tetra_node_ids);
              stk::io::put_io_part_attribute(  mesh.m_block_tet );
              mesh.m_metaData.commit();
              mesh.populate();

              std::cout << "here" << std::endl;
              bool isCommitted = true;
              percept::PerceptMesh eMesh(&mesh.m_metaData, &mesh.m_bulkData, isCommitted);
              eMesh.saveAs(input_files_loc+"local_tet_2.e");
            }

            {
              PerceptMesh eMesh;
              eMesh.open(input_files_loc+"local_tet_2.e");
              Local_Tet4_Tet4_N break_tet(eMesh);
              eMesh.commit();

              TestLocalRefinerTet_N_2 breaker(eMesh, break_tet, 0, 1);
              breaker.setRemoveOldElements(false);
              breaker.doBreak();

              save_or_diff(eMesh, output_files_loc+"local_tet_N_2tet_1.e");
            }

          }
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================


      /// check triangulate_tet - two tets sharing a face, random coords

      // Pathscale is the only platform that doesn't pass this test
#ifndef __PATHSCALE__
      STKUNIT_UNIT_TEST(unit_localRefiner, triangulate_tet_2_rand)
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size(pm);

        if (p_size <= 1)
          {
            srandom(1234);
            int ncases = 100;
            int icase = 0;
            for (int jcase = 0; jcase < ncases; jcase++)
            {

              // create the mesh
              unsigned npts=5;
              unsigned ntets=2;
              SingleTetFixture::Point node_coord_data[] = {
                { 0 , 0 , 0 } , { 1 , 0 , 0 } , { 0 , 1 , 0 } , { 0 , 0 , 1 }, {1, 1, 1} };

              // Hard coded tetra node ids for all the tetra nodes in the entire mesh
              static  SingleTetFixture::TetIds tetra_node_ids[] = {
                { 1, 2, 3, 4}, {2, 3, 4, 5} };

              for (unsigned ipts = 0; ipts < npts; ipts++)
                {
                  for (int ii=0; ii < 3; ii++) node_coord_data[ipts][ii] = ((double)random())/((double)RAND_MAX);
                }

              double vol0 = tet_volume(node_coord_data, tetra_node_ids[0], 1);
              double vol1 = tet_volume(node_coord_data, tetra_node_ids[1], 1);
              std::cout << "tmp vol0 = " << vol0 << " vol1 = " << vol1 << std::endl;
              if (vol0 < 0.0 || vol1 < 0.0)
                continue;

              {
                stk::percept::SingleTetFixture mesh(pm, false, npts, node_coord_data, ntets, tetra_node_ids);
                stk::io::put_io_part_attribute(  mesh.m_block_tet );
                mesh.m_metaData.commit();
                mesh.populate();



                bool isCommitted = true;
                percept::PerceptMesh eMesh(&mesh.m_metaData, &mesh.m_bulkData, isCommitted);


                eMesh.saveAs(input_files_loc+"local_tet_2_rand.e."+toString(icase));
              }

              {
                PerceptMesh eMesh;
                eMesh.open(input_files_loc+"local_tet_2_rand.e."+toString(icase));
                Local_Tet4_Tet4_N break_tet(eMesh);
                eMesh.commit();
                double totalVol0 = totalVolume(eMesh);
                std::cout << "tmp totalVol0 = " << totalVol0 << std::endl;

                int edge_mark_bitcode = 0;
                edge_mark_bitcode = (int)(63.*((double)random())/((double)RAND_MAX));
                std::cout << "tmp edge_mark_bitcode= " << edge_mark_bitcode << std::endl;
                if (edge_mark_bitcode <= 0) edge_mark_bitcode = 1;
                if (edge_mark_bitcode >= 63) edge_mark_bitcode = 63;

                TestLocalRefinerTet_N_3_1 breaker(eMesh, break_tet, 0, edge_mark_bitcode);
                // we do this (for now) since old elements that aren't refined are being removed (FIXME in Refiner.cpp)
                breaker.setRemoveOldElements(false);
                breaker.doBreak();
                breaker.deleteParentElements();

                double totalVol1 = totalVolume(eMesh);
                std::cout << "tmp edge_mark_bitcode= " << edge_mark_bitcode << " totalVol1 = " << totalVol1 << std::endl;

                if (std::abs(totalVol0 - totalVol1) > 1.e-6)
                  {
                    eMesh.saveAs( output_files_loc+"local_tet_2_rand_error."+toString(icase)+".e" );
                    throw std::runtime_error("triangulate_tet_2_rand:: error, volumes don't match");
                  }
                save_or_diff(eMesh, output_files_loc+"local_tet_2_rand.e."+toString(icase) );
              }

              {
                PerceptMesh eMesh;
                eMesh.open(input_files_loc+"local_tet_2_rand.e."+toString(icase));
                Local_Tet4_Tet4_N break_tet(eMesh);
                eMesh.commit();

                //eMesh.printInfo("srk tmp SingleTetFixture 2a", 2);

                int edge_mark_bitcode = 0;
                edge_mark_bitcode = (int)(63.*((double)random())/((double)RAND_MAX));
                if (edge_mark_bitcode <= 0) edge_mark_bitcode = 1;
                if (edge_mark_bitcode >= 63) edge_mark_bitcode = 63;

                TestLocalRefinerTet_N_3_1 breaker(eMesh, break_tet, 0, edge_mark_bitcode);
                breaker.setRemoveOldElements(false);
                breaker.doBreak();

                //eMesh.printInfo("srk tmp SingleTetFixture 3", 2);

#if 1
                MeshUtil::m_debug = true;
                bool isConsistent = percept::MeshUtil::facesConsistent(eMesh);
                MeshUtil::m_debug = false;
                if (!isConsistent)
                  {
                    std::cout << "tmp error isConsistent= " << isConsistent << std::endl;
                    eMesh.saveAs( output_files_loc+"local_tet_2_rand_error.e."+toString(icase) );
                  }
                STKUNIT_EXPECT_TRUE(isConsistent);
#endif
              }

              ++icase;
            }
          }
      }
#endif

      //=============================================================================
      //=============================================================================
      //=============================================================================
      /// check triangulate_tet - all 64 cases for a single tet

      STKUNIT_UNIT_TEST(unit_localRefiner, triangulate_tet_64)
      {
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size(pm);

        if (p_size <= 1)
          {
            {

              // create the mesh
              static  SingleTetFixture::Point node_coord_data[  ] = {
                { 0 , 0 , 0 } , { 1 , 0 , 0 } , { 0 , 1 , 0 } , { 0 , 0 , 1 } };

              SingleTetFixture::Point pts[64*4];

              // Hard coded tetra node ids for all the tetra nodes in the entire mesh
              static  SingleTetFixture::TetIds tets[64];

              unsigned ntets = 0;
              unsigned npts = 0;
              unsigned edge_mark_bitcode = 0u;
              unsigned ipts = 0;
              int is = 0;
              int ie = 7;
              int js = 0;
              int je = 7;
#if 0
              is = 0;
              ie = 2;
              js = 2;
              je = 6;
#endif
              for (int j = js; j <= je; j++)
                {
                  for (int i = is; i <= ie; i++)
                    {
                      //std::cout << "\ntmp i = " << i << " j= " << j << "\n--------------------------------\n" << std::endl;

                      for (int k = 0; k < 4; k++)
                        {
                          //pts[ipts][0] = node_coord_data[k][0] + i * 2;
                          //pts[ipts][1] = node_coord_data[k][1] + j * 2;
                          pts[ipts][0] = node_coord_data[k][0] + i * 1.25;
                          pts[ipts][1] = node_coord_data[k][1] + j * 1.25;
                          pts[ipts][2] = node_coord_data[k][2];
                          ++npts;
                          ++ipts;
                        }
                      
                      tets[edge_mark_bitcode][0] = edge_mark_bitcode*4 + 1;
                      tets[edge_mark_bitcode][1] = edge_mark_bitcode*4 + 2;
                      tets[edge_mark_bitcode][2] = edge_mark_bitcode*4 + 3;
                      tets[edge_mark_bitcode][3] = edge_mark_bitcode*4 + 4;
                      
                      ++edge_mark_bitcode;
                      ++ntets;
                    }
                }

              stk::percept::SingleTetFixture mesh(pm, false, npts, pts, ntets, tets);
              stk::io::put_io_part_attribute(  mesh.m_block_tet );
              mesh.m_metaData.commit();
              mesh.populate();

              bool isCommitted = true;
              percept::PerceptMesh eMesh(&mesh.m_metaData, &mesh.m_bulkData, isCommitted);

              if (0)
                {
                  percept::GeometryVerifier gv(true);
                  std::cout << "tmp GeometryVerifier= " << eMesh.getBulkData() << std::endl;
                  bool igb = gv.isGeometryBad(*eMesh.getBulkData(), true);
                  std::cout << "tmp isGeometryBad= " << igb << std::endl;
                }

              double totalVol0 = totalVolume(eMesh);
              std::cout << "tmp 64 totalVol0 = " << totalVol0 << std::endl;

              eMesh.saveAs(input_files_loc+"local_tet_64.e");
            }

            {
              PerceptMesh eMesh;
              eMesh.open(input_files_loc+"local_tet_64.e");
              Local_Tet4_Tet4_N break_tet(eMesh);
              eMesh.commit();

              TestLocalRefinerTet_N_3 breaker(eMesh, break_tet, 0);
              breaker.setRemoveOldElements(true);
              breaker.doBreak();
              //breaker.deleteParentElements();

              if (0)
                {
                  percept::GeometryVerifier gv(true);
                  std::cout << "tmp GeometryVerifier= " << eMesh.getBulkData() << std::endl;
                  bool igb = gv.isGeometryBad(*eMesh.getBulkData(), true);
                  std::cout << "tmp isGeometryBad= " << igb << std::endl;
                }
              
              double totalVol1 = totalVolume(eMesh);
              std::cout << "tmp 64 totalVol1 = " << totalVol1 << std::endl;

              save_or_diff(eMesh, output_files_loc+"local_tet_N_3_64tet_1.e");
            }

          }
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================
      /// check triangulate_tet

      STKUNIT_UNIT_TEST(unit_localRefiner, triangulate_tet_planes)
      {
        EXCEPTWATCH;
        fixture_setup_NxNxN_box_hex_and_tet_mesh();
        MPI_Barrier( MPI_COMM_WORLD );

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size(pm);

        if (p_size <= 3)
          {

            {
              PerceptMesh eMesh;

              eMesh.open(input_files_loc+"tet_fixture_NxNxN.e");
              Local_Tet4_Tet4_N break_tet(eMesh);
              int scalarDimension = 0; // a scalar
              stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
              eMesh.commit();

              TestLocalRefinerTet_N_4 breaker(eMesh, break_tet, proc_rank_field);
              breaker.setRemoveOldElements(false);
              breaker.setAlwaysInitializeNodeRegistry(false);
              
              int nref = 3;
              int nunref = 5;
              for (int ipass = 0; ipass < nref; ipass++)
                {
                  breaker.doBreak();
                  eMesh.saveAs( output_files_loc+"local_tet_N_4_planes_iref_"+toString(ipass)+".e");
                }

              //save_or_diff(eMesh, output_files_loc+"local_tet_N_4_planes.e");
              eMesh.saveAs( output_files_loc+"local_tet_N_4_planes.e");

              for (int iunref_pass=0; iunref_pass < nunref; ++iunref_pass)
                {
                  ElementUnrefineCollection elements_to_unref = breaker.buildTestUnrefineList();
                  breaker.unrefineTheseElements(elements_to_unref);
                }

              //save_or_diff(eMesh, output_files_loc+"local_tet_N_4_planes_unref.e");
              eMesh.saveAs( output_files_loc+"local_tet_N_4_planes_unref.e");

              breaker.deleteParentElements();
              eMesh.saveAs( output_files_loc+"local_tet_N_4_planes_unref_noParentElements.e");

            }

          }
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

            eMesh.printInfo("local tri mesh", printInfoLevel);
            save_or_diff(eMesh, output_files_loc+"local_tri_0.e");

            TestLocalRefinerTri breaker(eMesh, break_tri_to_tri_2, proc_rank_field);
            breaker.setRemoveOldElements(false);
            breaker.doBreak();

            eMesh.printInfo("local tri mesh refined", printInfoLevel);
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

            eMesh.printInfo("local tri mesh", printInfoLevel);
            save_or_diff(eMesh, output_files_loc+"local_tri_1_0.e");

            bool diagonals=true;
            TestLocalRefinerTri2 breaker(eMesh, break_tri_to_tri_2, proc_rank_field, diagonals);
            //breaker.setRemoveOldElements(false);
            breaker.doBreak();

            eMesh.printInfo("local tri mesh refined", printInfoLevel);
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

            eMesh.printInfo("local tri mesh", printInfoLevel);
            save_or_diff(eMesh, output_files_loc+"local_tri_2_0.e");

            bool diagonals=false;
            TestLocalRefinerTri2 breaker(eMesh, break_tri_to_tri_2, proc_rank_field, diagonals);
            //breaker.setRemoveOldElements(false);
            breaker.doBreak();

            eMesh.printInfo("local tri mesh refined",  printInfoLevel);
            save_or_diff(eMesh, output_files_loc+"local_tri_2_1.e");

            // end_demo
          }

      }


      //=============================================================================
      //=============================================================================
      //=============================================================================

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

            //eMesh.printInfo("local tri mesh",2);
            save_or_diff(eMesh, output_files_loc+"local_tri_N_1_0_"+post_fix[p_size]+".e");

            TestLocalRefinerTri_N_1 breaker(eMesh, break_tri_to_tri_N, proc_rank_field);
            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);
            breaker.doBreak();

            //eMesh.dumpElementsCompact();

            //eMesh.printInfo("local tri mesh refined", 2);
            save_or_diff(eMesh, output_files_loc+"local_tri_N_1_1_"+post_fix[p_size]+".e");

            //MPI_Barrier( MPI_COMM_WORLD );

            //breaker.unrefineAll();
            ElementUnrefineCollection elements_to_unref = breaker.buildTestUnrefineList();
            breaker.unrefineTheseElements(elements_to_unref);

            // FIXME
            eMesh.saveAs( output_files_loc+"local_tri_N_1_1_unref_"+post_fix[p_size]+".e");
            //save_or_diff(eMesh, output_files_loc+"local_tri_N_1_1_unref_"+post_fix[p_size]+".e");

            // end_demo
          }

      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      STKUNIT_UNIT_TEST(unit_localRefiner, break_tri_to_tri_N_2)
      {
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

            Local_Tri3_Tri3_N break_tri_to_tri_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
            eMesh.addField("proc_rank_edge", eMesh.edge_rank(), scalarDimension);
            eMesh.commit();

            fixture.generate_mesh();

            //eMesh.printInfo("local tri mesh",2);
            save_or_diff(eMesh, output_files_loc+"local_tri_N_2_0_"+post_fix[p_size]+".e");

            TestLocalRefinerTri_N_2 breaker(eMesh, break_tri_to_tri_N, proc_rank_field);
            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);
            for (int ipass=0; ipass < 4; ipass++)
              {
                std::cout << "P[" << eMesh.getRank() << "] ipass= " << ipass << std::endl;
                breaker.doBreak();
              }

            //eMesh.dumpElementsCompact();

            //eMesh.printInfo("local tri mesh refined", 2);
            //save_or_diff(eMesh, output_files_loc+"local_tri_N_2_1_"+post_fix[p_size]+".e");
            eMesh.saveAs(output_files_loc+"local_tri_N_2_1_"+post_fix[p_size]+".e");

            //MPI_Barrier( MPI_COMM_WORLD );
#if 1           

            //breaker.unrefineAll();
            ElementUnrefineCollection elements_to_unref = breaker.buildTestUnrefineList();
            breaker.unrefineTheseElements(elements_to_unref);

            // FIXME
            eMesh.saveAs( output_files_loc+"local_tri_N_2_1_unref_"+post_fix[p_size]+".e");
            //save_or_diff(eMesh, output_files_loc+"local_tri_N_2_1_unref_"+post_fix[p_size]+".e");
#endif
            // end_demo
          }

      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      STKUNIT_UNIT_TEST(unit_localRefiner, break_tri_to_tri_N_3_1)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size <= 3)
          {

            // start_demo_local_refiner_break_tri_to_tri_2

            const unsigned n = 2;
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

            //eMesh.printInfo("local tri mesh",2);
            save_or_diff(eMesh, output_files_loc+"local_tri_N_3_1_0_"+post_fix[p_size]+".e");

            TestLocalRefinerTri_N_3 breaker(eMesh, break_tri_to_tri_N, proc_rank_field);
            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);
            for (int ipass=0; ipass < 4; ipass++)
              {
                std::cout << "P[" << eMesh.getRank() << "] ipass= " << ipass << std::endl;
                breaker.doBreak();
                eMesh.saveAs(output_files_loc+"local_tri_N_3_1_1_ipass"+toString(ipass)+"_"+post_fix[p_size]+".e");
              }

            //eMesh.dumpElementsCompact();

            //eMesh.printInfo("local tri mesh refined", 2);
            //save_or_diff(eMesh, output_files_loc+"local_tri_N_3_1_1_"+post_fix[p_size]+".e");
            eMesh.saveAs(output_files_loc+"local_tri_N_3_1_1_"+post_fix[p_size]+".e");

            //MPI_Barrier( MPI_COMM_WORLD );
#if 1

            for (int iunref_pass=0; iunref_pass < 4; iunref_pass++)
              {
                std::cout << "P[" << eMesh.getRank() << "] iunref_pass= " << iunref_pass << std::endl;
                ElementUnrefineCollection elements_to_unref = breaker.buildTestUnrefineList();
                breaker.unrefineTheseElements(elements_to_unref);
                eMesh.saveAs(output_files_loc+"local_tri_N_3_1_1_unref_ipass"+toString(iunref_pass)+"_"+post_fix[p_size]+".e");
                //breaker.unrefineAll();
              }

            // FIXME
            eMesh.saveAs( output_files_loc+"local_tri_N_3_1_1_unref_"+post_fix[p_size]+".e");
            //save_or_diff(eMesh, output_files_loc+"local_tri_N_3_1_1_unref_"+post_fix[p_size]+".e");
#endif
            // end_demo
          }

      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      STKUNIT_UNIT_TEST(unit_localRefiner, break_tri_to_tri_N_3_2)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size <= 3)
          {

            // start_demo_local_refiner_break_tri_to_tri_2

            const unsigned n = 2;
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

            //eMesh.printInfo("local tri mesh",2);
            save_or_diff(eMesh, output_files_loc+"local_tri_N_3_2_0_"+post_fix[p_size]+".e");

            TestLocalRefinerTri_N_3 breaker(eMesh, break_tri_to_tri_N, proc_rank_field);
            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);
            for (int ipass=0; ipass < 8; ipass++)
              {
                std::cout << "P[" << eMesh.getRank() << "] ipass= " << ipass << std::endl;
                breaker.doBreak();
                eMesh.saveAs(output_files_loc+"local_tri_N_3_2_1_ipass"+toString(ipass)+"_"+post_fix[p_size]+".e");
              }

            //eMesh.dumpElementsCompact();

            //eMesh.printInfo("local tri mesh refined", 2);
            //save_or_diff(eMesh, output_files_loc+"local_tri_N_3_2_1_"+post_fix[p_size]+".e");
            eMesh.saveAs(output_files_loc+"local_tri_N_3_2_1_"+post_fix[p_size]+".e");

            //MPI_Barrier( MPI_COMM_WORLD );
#if 1

            for (int iunref_pass=0; iunref_pass < 7; iunref_pass++)
              {
                std::cout << "P[" << eMesh.getRank() << "] iunref_pass= " << iunref_pass << std::endl;
                ElementUnrefineCollection elements_to_unref = breaker.buildTestUnrefineList();
                breaker.unrefineTheseElements(elements_to_unref);
                eMesh.saveAs(output_files_loc+"local_tri_N_3_2_1_unref_ipass"+toString(iunref_pass)+"_"+post_fix[p_size]+".e");
                //breaker.unrefineAll();
              }



            // FIXME
            eMesh.saveAs( output_files_loc+"local_tri_N_3_2_1_unref_"+post_fix[p_size]+".e");
            //save_or_diff(eMesh, output_files_loc+"local_tri_N_3_2_1_unref_"+post_fix[p_size]+".e");

            if (1)
              {
                for (int ipass=8; ipass < 16; ipass++)
                  {
                    std::cout << "P[" << eMesh.getRank() << "] ipass= " << ipass << std::endl;
                    breaker.doBreak();
                    eMesh.saveAs(output_files_loc+"local_tri_N_3_2_1_ipass"+toString(ipass)+"_"+post_fix[p_size]+".e");
                  }

                //eMesh.dumpElementsCompact();

                // FIXME FIXME FIXME
                breaker.deleteParentElements();

                //eMesh.printInfo("local tri mesh refined", 2);

                //save_or_diff(eMesh, output_files_loc+"local_tri_N_3_2_1_"+post_fix[p_size]+".e");
                eMesh.saveAs(output_files_loc+"local_tri_N_3_2_16_"+post_fix[p_size]+".e");
              }


#endif
            // end_demo
          }

      }


#if 0  
      //=============================================================================
      //=============================================================================
      //=============================================================================

      STKUNIT_UNIT_TEST(unit_localRefiner, break_tri_to_tri_N_3)
      {
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

            Local_Tri3_Tri3_N break_tri_to_tri_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
            eMesh.addField("proc_rank_edge", eMesh.edge_rank(), scalarDimension);
            eMesh.commit();

            fixture.generate_mesh();

            //eMesh.printInfo("local tri mesh",2);
            save_or_diff(eMesh, output_files_loc+"local_tri_N_3_0_"+post_fix[p_size]+".e");

            TestLocalRefinerTri_N_3 breaker(eMesh, break_tri_to_tri_N, proc_rank_field);
            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);
            for (int ipass=0; ipass < 4; ipass++)
              {
                std::cout << "P[" << eMesh.getRank() << "] ipass= " << ipass << std::endl;
                breaker.doBreak();
              }

            //eMesh.dumpElementsCompact();

            //eMesh.printInfo("local tri mesh refined", 2);
            //save_or_diff(eMesh, output_files_loc+"local_tri_N_3_2_"+post_fix[p_size]+".e");
            eMesh.saveAs(output_files_loc+"local_tri_N_3_2_"+post_fix[p_size]+".e");

            //MPI_Barrier( MPI_COMM_WORLD );
#if 1

            for (int iunref_pass=0; iunref_pass < 2; iunref_pass++)
              {
                std::cout << "P[" << eMesh.getRank() << "] iunref_pass= " << iunref_pass << std::endl;
                ElementUnrefineCollection elements_to_unref = breaker.buildTestUnrefineList();
                breaker.unrefineTheseElements(elements_to_unref);
                //breaker.unrefineAll();
              }

            // FIXME
            eMesh.saveAs( output_files_loc+"local_tri_N_3_2_unref_"+post_fix[p_size]+".e");
            //save_or_diff(eMesh, output_files_loc+"local_tri_N_3_2_unref_"+post_fix[p_size]+".e");
            exit(123);
#endif
            // end_demo
          }

      }
#endif

#if 0
      //=============================================================================
      //=============================================================================
      //=============================================================================


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

            eMesh.printInfo("local tri mesh", printInfoLevel);
            save_or_diff(eMesh, output_files_loc+"local_tri_N_0.e");

            TestLocalRefinerTri_N breaker(eMesh, break_tri_to_tri_N, proc_rank_field);
            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);
            breaker.doBreak();

            eMesh.printInfo("local tri mesh refined", printInfoLevel);
            //eMesh.dumpElements();
            save_or_diff(eMesh, output_files_loc+"local_tri_N_1.e");

            //breaker.unrefineAll();
            ElementUnrefineCollection elements_to_unref = breaker.buildTestUnrefineList();
            breaker.unrefineTheseElements(elements_to_unref);

            // FIXME
            eMesh.saveAs( output_files_loc+"local_tri_N_1_unref.e");
            //save_or_diff(eMesh, output_files_loc+"local_tri_N_1_unref.e");

            // end_demo
          }

      }
#endif



      //=============================================================================
      //=============================================================================
      //=============================================================================

      STKUNIT_UNIT_TEST(unit_localRefiner, break_tri_to_tri_N_3_1_EM)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size <= 3)
          {

            // start_demo_local_refiner_break_tri_to_tri_2

            const unsigned n = 2;
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

            //eMesh.printInfo("local tri mesh",2);
            save_or_diff(eMesh, output_files_loc+"local_tri_N_3_1_EdgeMarker_0_"+post_fix[p_size]+".e");

            TestLocalRefinerTri_N_3_EdgeMarker breaker(eMesh, break_tri_to_tri_N, proc_rank_field);
            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);
            for (int ipass=0; ipass < 4; ipass++)
              {
                std::cout << "P[" << eMesh.getRank() << "] ipass= " << ipass << std::endl;
                breaker.doBreak();
                eMesh.saveAs(output_files_loc+"local_tri_N_3_1_EdgeMarker_1_ipass"+toString(ipass)+"_"+post_fix[p_size]+".e");
              }

            //eMesh.dumpElementsCompact();

            //eMesh.printInfo("local tri mesh refined", 2);
            //save_or_diff(eMesh, output_files_loc+"local_tri_N_3_1_EdgeMarker_1_"+post_fix[p_size]+".e");
            eMesh.saveAs(output_files_loc+"local_tri_N_3_1_EdgeMarker_1_"+post_fix[p_size]+".e");

            //MPI_Barrier( MPI_COMM_WORLD );
#if 1

            for (int iunref_pass=0; iunref_pass < 4; iunref_pass++)
              {
                std::cout << "P[" << eMesh.getRank() << "] iunref_pass= " << iunref_pass << std::endl;
                ElementUnrefineCollection elements_to_unref = breaker.buildUnrefineList();
                breaker.unrefineTheseElements(elements_to_unref);
                eMesh.saveAs(output_files_loc+"local_tri_N_3_1_EdgeMarker_1_unref_ipass"+toString(iunref_pass)+"_"+post_fix[p_size]+".e");
                //breaker.unrefineAll();
              }

            // FIXME
            eMesh.saveAs( output_files_loc+"local_tri_N_3_1_EdgeMarker_1_unref_"+post_fix[p_size]+".e");
            //save_or_diff(eMesh, output_files_loc+"local_tri_N_3_1_EdgeMarker_1_unref_"+post_fix[p_size]+".e");
#endif
            // end_demo
          }

      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

#if 1
      STKUNIT_UNIT_TEST(unit_localRefiner, break_tri_to_tri_N_3_1_ElementMarker)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size <= 3)
          {

            // start_demo_local_refiner_break_tri_to_tri_2

            const unsigned n = 2;
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

            //eMesh.printInfo("local tri mesh",2);
            save_or_diff(eMesh, output_files_loc+"local_tri_N_3_1_ElementMarker_0_"+post_fix[p_size]+".e");

            TestLocalRefinerTri_N_3_ElementMarker breaker(eMesh, break_tri_to_tri_N, proc_rank_field);
            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);
            for (int ipass=0; ipass < 4; ipass++)
              {
                std::cout << "P[" << eMesh.getRank() << "] ipass= " << ipass << std::endl;
                breaker.doBreak();
                std::cout << "P[" << eMesh.getRank() << "] done... ipass= " << ipass << std::endl;
                eMesh.saveAs(output_files_loc+"local_tri_N_3_1_ElementMarker_1_ipass"+toString(ipass)+"_"+post_fix[p_size]+".e");
              }

            //eMesh.dumpElementsCompact();

            //eMesh.printInfo("local tri mesh refined", 2);
            //save_or_diff(eMesh, output_files_loc+"local_tri_N_3_1_ElementMarker_1_"+post_fix[p_size]+".e");
            eMesh.saveAs(output_files_loc+"local_tri_N_3_1_ElementMarker_1_"+post_fix[p_size]+".e");

            //MPI_Barrier( MPI_COMM_WORLD );
#if 1

            for (int iunref_pass=0; iunref_pass < 4; iunref_pass++)
              {
                std::cout << "P[" << eMesh.getRank() << "] iunref_pass= " << iunref_pass << std::endl;
                ElementUnrefineCollection elements_to_unref = breaker.buildUnrefineList();
                breaker.unrefineTheseElements(elements_to_unref);
                eMesh.saveAs(output_files_loc+"local_tri_N_3_1_ElementMarker_1_unref_ipass"+toString(iunref_pass)+"_"+post_fix[p_size]+".e");
                //breaker.unrefineAll();
              }

            // FIXME
            eMesh.saveAs( output_files_loc+"local_tri_N_3_1_ElementMarker_1_unref_"+post_fix[p_size]+".e");
            //save_or_diff(eMesh, output_files_loc+"local_tri_N_3_1_ElementMarker_1_unref_"+post_fix[p_size]+".e");
#endif
            // end_demo
          }

      }
#endif

      //=============================================================================
      //=============================================================================
      //=============================================================================

#if 1
      class SetRefineField : public percept::ElementOp
      {
        percept::PerceptMesh& m_eMesh;
      public:
        SetRefineField(percept::PerceptMesh& eMesh) : m_eMesh(eMesh) {}
        virtual bool operator()(const stk::mesh::Entity& element, stk::mesh::FieldBase *field,  const mesh::BulkData& bulkData)
        {
          const mesh::PairIterRelation elem_nodes = element.relations( stk::mesh::fem::FEMMetaData::NODE_RANK );
          unsigned num_node = elem_nodes.size();
          double *f_data = PerceptMesh::field_data_entity(field, element);
          VectorFieldType* coordField = m_eMesh.getCoordinatesField();
                
          bool found = true;
          for (unsigned inode=0; inode < num_node; inode++)
            {
              mesh::Entity & node = * elem_nodes[ inode ].entity();
              double *coord_data = PerceptMesh::field_data(coordField, node);

              //std::cout << "tmp coord_data= " << coord_data[0] << std::endl;

              if (coord_data[0] > 1.1)
                {
                  found=false;
                  break;
                }
            }
          if (found)
            f_data[0] = 1.0;
          else
            f_data[0] = 0.0;

          return false;  // don't terminate the loop
        }
        virtual void init_elementOp() {}
        virtual void fini_elementOp() {}
      };

      class SetUnrefineField : public percept::ElementOp
      {
        percept::PerceptMesh& m_eMesh;
      public:
        SetUnrefineField(percept::PerceptMesh& eMesh) : m_eMesh(eMesh) {}
        virtual bool operator()(const stk::mesh::Entity& element, stk::mesh::FieldBase *field,  const mesh::BulkData& bulkData)
        {
          const mesh::PairIterRelation elem_nodes = element.relations( stk::mesh::fem::FEMMetaData::NODE_RANK );
          unsigned num_node = elem_nodes.size();
          double *f_data = PerceptMesh::field_data_entity(field, element);
          VectorFieldType* coordField = m_eMesh.getCoordinatesField();
                
          bool found = true;
          for (unsigned inode=0; inode < num_node; inode++)
            {
              mesh::Entity & node = * elem_nodes[ inode ].entity();
              double *coord_data = PerceptMesh::field_data(coordField, node);

              //std::cout << "tmp coord_data= " << coord_data[0] << std::endl;

              if (coord_data[0] > 1.1 || coord_data[1] > 1.1)
                {
                  found=false;
                  break;
                }
            }
          if (found)
            f_data[0] = -1.0;
          else
            f_data[0] = 0.0;

          return false;  // don't terminate the loop
        }
        virtual void init_elementOp() {}
        virtual void fini_elementOp() {}
      };

      STKUNIT_UNIT_TEST(unit_localRefiner, break_tri_to_tri_N_5_FieldBased)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size <= 3)
          {
            // start_demo_local_refiner_break_tri_to_tri_2

            const unsigned n = 2;
            const unsigned nx = n , ny = n;

            bool createEdgeSets = false;
            percept::QuadFixture<double, shards::Triangle<3> > fixture( pm , nx , ny, createEdgeSets);

            bool isCommitted = false;
            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, isCommitted);

            Local_Tri3_Tri3_N break_tri_to_tri_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
            stk::mesh::FieldBase* refine_field = eMesh.addField("refine_field", eMesh.element_rank(), scalarDimension);
            stk::mesh::FieldBase* unrefine_field = eMesh.addField("unrefine_field", eMesh.element_rank(), scalarDimension);
            eMesh.commit();

            fixture.generate_mesh();

            SetRefineField set_ref_field(eMesh);
            eMesh.elementOpLoop(set_ref_field, refine_field);

            SetUnrefineField set_unref_field(eMesh);
            eMesh.elementOpLoop(set_unref_field, unrefine_field);
            
            save_or_diff(eMesh, output_files_loc+"local_tri_N_5_FieldBased_0_"+post_fix[p_size]+".e");

            stk::mesh::Selector univ_selector(eMesh.getFEM_meta_data()->universal_part());

            PredicateBasedMarker<ElementFieldBasedRefinePredicate, ElementFieldBasedUnrefinePredicate>
              breaker(ElementFieldBasedRefinePredicate(univ_selector, refine_field, 0.0),
                      ElementFieldBasedUnrefinePredicate(univ_selector, unrefine_field, 0.0),
                      eMesh, break_tri_to_tri_N, proc_rank_field);

            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);
            for (int ipass=0; ipass < 4; ipass++)
              {
                std::cout << "P[" << eMesh.getRank() << "] ipass= " << ipass << std::endl;
                breaker.doBreak();
                std::cout << "P[" << eMesh.getRank() << "] done... ipass= " << ipass << std::endl;
                eMesh.saveAs(output_files_loc+"local_tri_N_5_FieldBased_1_ipass"+toString(ipass)+"_"+post_fix[p_size]+".e");
              }

            eMesh.saveAs(output_files_loc+"local_tri_N_5_FieldBased_1_"+post_fix[p_size]+".e");

#if 1
            for (int iunref_pass=0; iunref_pass < 4; iunref_pass++)
              {
                std::cout << "P[" << eMesh.getRank() << "] iunref_pass= " << iunref_pass << std::endl;
                ElementUnrefineCollection elements_to_unref = breaker.buildUnrefineList();
                breaker.unrefineTheseElements(elements_to_unref);
                eMesh.saveAs(output_files_loc+"local_tri_N_5_FieldBased_1_unref_ipass"+toString(iunref_pass)+"_"+post_fix[p_size]+".e");
              }

            eMesh.saveAs( output_files_loc+"local_tri_N_5_FieldBased_1_unref_"+post_fix[p_size]+".e");
#endif
            // end_demo
          }

      }
#endif

      //=============================================================================
      //=============================================================================
      //=============================================================================

      struct SingleTriangleFixture
      {
        static PerceptMesh * create()
        {
            const unsigned n = 1;
            const unsigned nx = n , ny = n;

            stk::ParallelMachine pm = MPI_COMM_WORLD ;
            bool createEdgeSets = false;
            percept::QuadFixture<double, shards::Triangle<3> > *fixture = new percept::QuadFixture<double, shards::Triangle<3> >( pm , nx , ny, createEdgeSets);

            bool isCommitted = false;
            percept::PerceptMesh * eMesh = new percept::PerceptMesh(&fixture->meta_data, &fixture->bulk_data, isCommitted);

            eMesh->commit();

            fixture->generate_mesh();

            // delete the first element
            eMesh->getBulkData()->modification_begin();
            stk::mesh::Entity* element = &( (**(eMesh->getBulkData()->buckets(eMesh->element_rank()).begin()))[0]);
            if ( ! eMesh->getBulkData()->destroy_entity( element ) )
              {
                throw std::logic_error("failed in deleting element");
              }
            eMesh->getBulkData()->modification_end();

            // single element left
            //stk::mesh::Entity& element = (**(eMesh->getBulkData()->buckets(eMesh->element_rank()).begin()))[0];

            return eMesh;
        }
      };


      static void set_node_coords(percept::PerceptMesh& eMesh, mesh::PairIterRelation& elem_nodes, double tri_coords[3][3])
      {
        for (unsigned inode=0; inode < elem_nodes.size(); inode++)
          {
            stk::mesh::Entity *node = elem_nodes[inode].entity();
            double *fdata = stk::mesh::field_data( *eMesh.getCoordinatesField() , *node );
            for (int dim=0; dim < eMesh.getSpatialDim(); dim++)
              {
                fdata[dim] = tri_coords[inode][dim];
              }
          }
      }

      std::vector<int> convert_tuple(tri_tuple_type_local& tuple)
      {
        std::vector<int> cv(3);
        cv[0] = tuple.get<0>();
        cv[1] = tuple.get<1>();
        cv[2] = tuple.get<2>();
        return cv;
      }

      static bool in_set(tri_tuple_type_local& expected, vector<tri_tuple_type_local>& base, bool reverse=false)
      {
        std::vector<int> cv_expected = convert_tuple(expected);

        for (unsigned ie = 0; ie < base.size(); ie++)
          {
            std::vector<int> cv_base = convert_tuple(base[ie]);
            for (int i = 0; i < 3; i++)
              {
                bool found = true;
                if (reverse)
                  {
                    int k=0;
                    for (int j = 2; j >= 0; --j)
                      {
                        if (cv_expected[k++] != cv_base[(i+j)%3])
                          {
                            found = false;
                            break;
                          }
                      }
                  }
                else
                  {
                    for (int j = 0; j < 3; j++)
                      {
                        if (cv_expected[j] != cv_base[(i+j)%3])
                          {
                            found = false;
                            break;
                          }
                      }
                  }

                if (found)
                  {
                    return true;
                  }
              }
          }
        return false;
      }

      /// Create a single triangle mesh and mark the edges, call RefinerPattern_Tri3_Tri3_N::triangulate_face
      ///   and check properties of the result - reverse the triangle polarity and check for consistency

      STKUNIT_UNIT_TEST(unit_localRefiner, check_triangulate_face)
      {
        //fixture_setup();
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size <= 1)
          {
            percept::PerceptMesh& eMesh = *SingleTriangleFixture::create();

            eMesh.saveAs(output_files_loc+"tri_face_0.e");

            // single element left
            stk::mesh::Entity& element = (**(eMesh.getBulkData()->buckets(eMesh.element_rank()).begin()))[0];
            std::cout << "element = " << element << std::endl;

            mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);


            stk::mesh::Entity *elem_nodes_vector[3];
            for (unsigned inode=0; inode < elem_nodes.size(); inode++)
              {
                stk::mesh::Entity *node = elem_nodes[inode].entity();
                elem_nodes_vector[inode] = node;
              }
            vector<tri_tuple_type_local> elems_local;

            // test 1
            {
              unsigned edge_marks[3] = {1,1,0};
              double tri_coords[3][3] = {{0,0,0}, {1,0,0}, {0,1,0}};
              set_node_coords(eMesh, elem_nodes, tri_coords);
              Local_Tri3_Tri3_N::triangulate_face(eMesh, elem_nodes_vector, edge_marks, elems_local);

              // expected:
              vector<tri_tuple_type_local> elems_local_expected(3);
              elems_local_expected[0] = tri_tuple_type_local(0,3,4);
              elems_local_expected[1] = tri_tuple_type_local(3,1,4);
              elems_local_expected[2] = tri_tuple_type_local(0,4,2);
              
              std::cout << "test1: elems_local_expected= " << elems_local_expected << std::endl;
              std::cout << "test1: elems_local= " << elems_local << std::endl;

              STKUNIT_EXPECT_TRUE(in_set(elems_local_expected[0], elems_local));
              STKUNIT_EXPECT_TRUE(in_set(elems_local_expected[1], elems_local));
              STKUNIT_EXPECT_TRUE(in_set(elems_local_expected[2], elems_local));

              STKUNIT_EXPECT_TRUE(!in_set(elems_local_expected[0], elems_local, true));
              STKUNIT_EXPECT_TRUE(!in_set(elems_local_expected[1], elems_local, true));
              STKUNIT_EXPECT_TRUE(!in_set(elems_local_expected[2], elems_local, true));
            }

            // test2: same as test 1 but mirror image (emulating a face shared between two tets)
            {
              stk::mesh::Entity* node1 = elem_nodes_vector[1];
              elem_nodes_vector[1] = elem_nodes_vector[2];
              elem_nodes_vector[2] = node1;

              unsigned edge_marks[3] = {0,1,1};
              Local_Tri3_Tri3_N::triangulate_face(eMesh, elem_nodes_vector, edge_marks, elems_local);

              // expected:
              vector<tri_tuple_type_local> elems_local_expected(3);
              elems_local_expected[0] = tri_tuple_type_local(0,1,4);
              elems_local_expected[1] = tri_tuple_type_local(0,4,5);
              elems_local_expected[2] = tri_tuple_type_local(2,5,4);
              
              std::cout << "test2: elems_local_expected= " << elems_local_expected << std::endl;
              std::cout << "test2: elems_local= " << elems_local << std::endl;

              STKUNIT_EXPECT_TRUE(in_set(elems_local_expected[0], elems_local));
              STKUNIT_EXPECT_TRUE(in_set(elems_local_expected[1], elems_local));
              STKUNIT_EXPECT_TRUE(in_set(elems_local_expected[2], elems_local));

              STKUNIT_EXPECT_TRUE(!in_set(elems_local_expected[0], elems_local, true));
              STKUNIT_EXPECT_TRUE(!in_set(elems_local_expected[1], elems_local, true));
              STKUNIT_EXPECT_TRUE(!in_set(elems_local_expected[2], elems_local, true));
            }

            
        
            // end_demo
          }

      }



    } // namespace unit_tests
  } // namespace adapt
} // namespace stk


