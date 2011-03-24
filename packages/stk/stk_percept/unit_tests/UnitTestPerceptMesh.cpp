/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
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
#include <stk_percept/fixtures/Fixture.hpp>
#include <stk_percept/fixtures/QuadFixture.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/util/PrintTable.hpp>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <Teuchos_ScalarTraits.hpp>

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>
#include <mpi.h>
#include <math.h>

namespace stk 
{
  using namespace mesh;
  namespace percept 
  {
    namespace unit_tests 
    {

#define EXTRA_PRINT 0

      //=============================================================================
      //=============================================================================
      //=============================================================================

      STKUNIT_UNIT_TEST(perceptMesh, walk_nodes)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size( pm );
        const unsigned p_rank = stk::parallel_machine_rank( pm );
        if (p_size <= 3)
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
            eMesh.printInfo("quad fixture", 2);
            //eMesh.saveAs("./output_files/quad_fixture.e");

            mesh::fem::FEMMetaData& metaData = *eMesh.getFEM_meta_data();

            const std::vector< stk::mesh::Part * > & parts = metaData.get_parts();

            unsigned nparts = parts.size();
            if (1) std::cout << "Number of parts = " << nparts << std::endl;

            int surface_id = 2;
            std::string surface_name = "surface_"+toString(surface_id);
            mesh::Part *part = eMesh.getNonConstPart(surface_name);
            mesh::Selector in_surface_selector(*part);
            mesh::BulkData& bulkData = *eMesh.get_bulkData();
            VectorFieldType* coordField = eMesh.getCoordinatesField();

            const std::vector<Bucket*> & buckets = bulkData.buckets( (eMesh.getSpatialDim() == 2 ? mesh::Edge : mesh::Face) );  // Note
            double sum = 0.0;

            for ( std::vector<Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
              {
                if (in_surface_selector(**k)) 
                  {
                    Bucket & bucket = **k ;

                    // in case the cell topology is needed
                    const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(bucket);
                    shards::CellTopology cell_topo(cell_topo_data);

                    const unsigned num_elements_in_bucket = bucket.size();
                
                    for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                      {
                        Entity& element = bucket[iElement];

                        const PairIterRelation& elem_nodes = element.relations( mesh::Node );  

                        unsigned num_node = elem_nodes.size(); 
                        for (unsigned inode=0; inode < num_node; inode++)
                          {
                            Entity & node = *elem_nodes[ inode ].entity();
                            //EntityId nid = node.identifier();

                            double * const coord = stk::mesh::field_data( *coordField , node );
                            // do something with coord's
                            sum += coord[0]*coord[0] + coord[1]*coord[1];
                          }
                      }
                  }
              }
            std::cout << "P[" << p_rank << ":" << p_size << "] sum = " << sum << std::endl;
          }
      }

    }
  }
}
