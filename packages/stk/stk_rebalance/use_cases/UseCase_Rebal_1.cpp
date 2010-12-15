/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <use_cases/UseCase_Rebal_1.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/fem/DefaultFEM.hpp>
#include <stk_mesh/fixtures/HexFixture.hpp>

#include <stk_rebalance/Rebalance.hpp>
#include <stk_rebalance/Partition.hpp>
#include <stk_rebalance/ZoltanPartition.hpp>

//----------------------------------------------------------------------

using namespace stk::mesh::fixtures;

typedef stk::mesh::Field<double> ScalarField ;

namespace stk {
namespace rebalance {
namespace use_cases {

bool test_unequal_weights( stk::ParallelMachine pm )
{
  const unsigned p_size = stk::parallel_machine_size(pm);
  const unsigned p_rank = stk::parallel_machine_rank(pm);

  const unsigned ngx = p_size*(p_size+1)/2;

  unsigned nx = 0;
  if( 0 == p_rank )
    nx = ngx;
  unsigned ny = 1;
  unsigned nz = 1;

  stk::mesh::fixtures::HexFixture fixture(pm, nx, ny, nz);

  stk::mesh::MetaData & meta  = fixture.m_meta_data;
  stk::mesh::BulkData & bulk  = fixture.m_bulk_data;
  stk::mesh::DefaultFEM & fem = fixture.m_fem;

  // Put weights field on all elements
  const stk::mesh::EntityRank element_rank = stk::mesh::fem::element_rank(fem);
  ScalarField & weight_field( meta.declare_field< ScalarField >( "element_weights" ) );
  stk::mesh::put_field(weight_field , element_rank , meta.universal_part() );

  meta.commit();

  bulk.modification_begin();

  // Initially put all elements on proc 0
  std::vector<stk::mesh::EntityId> my_element_ids;
  for ( unsigned i = 0 ; i < nx*ny*nz; ++i )
    my_element_ids.push_back(i+1);

  fixture.generate_mesh(my_element_ids);

  // Assign weights so that a perfect rebalance is possible so long as the rebalancer can figure out
  // to put p_rank+1 elements on proc = p_rank based on these weights.
  unsigned nslabs = 0;
  double total_weight = 0.0;
  if( 0 == p_rank ) {
    for ( unsigned l = 1 ; l <= p_size ; ++l ) {
      for ( unsigned k = 0 ; k < nz ; ++k ) {
        for ( unsigned j = 0 ; j < ny ; ++j ) {
          for ( unsigned i = 0 ; i < l ; ++i ) {
            const stk::mesh::EntityId elem_id = 1 + nslabs + i + j*ngx + k*ngx*ny; 
            stk::mesh::Entity * elem = bulk.get_entity(element_rank, elem_id);
            double * const e_weight = stk::mesh::field_data( weight_field , *elem );
            *e_weight = double(ngx) / double(l);
            //std::cerr << "(proc " << p_rank << ") : assigning weight = " << *e_weight << " to (nx, ny, nz) = (" << i+nslabs << ", " << j << ", " << k << ") = elem_id " << elem_id << std::endl;
            total_weight += *e_weight;
          }
        }
      }
      nslabs += l;
    }
  }

  MPI_Barrier( pm );

  stk::all_reduce(pm, ReduceMax<1>(&total_weight)& ReduceSum<1>(&total_weight));
  //std::cerr << "(proc " << p_rank << ") : has total_weight = " << total_weight << std::endl;

  bulk.modification_end();

  // Use Zoltan to determine new partition
  Teuchos::ParameterList emptyList;
  stk::rebalance::Zoltan zoltan_partition(pm, fixture.m_spatial_dimension, emptyList);

  stk::mesh::Selector selector(meta.universal_part());

  // Force a rebalance by using imbalance_threshold < 1.0
  double imbalance_threshold = 0.5;
  bool do_rebal = stk::rebalance::rebalance_needed(bulk, weight_field, pm, imbalance_threshold);
  // Coordinates are passed to support geometric-based load balancing algorithms
  if( do_rebal )
    stk::rebalance::rebalance(bulk, selector, &fixture.m_coord_field, &weight_field, zoltan_partition);

  imbalance_threshold = 1.5;
  do_rebal = stk::rebalance::rebalance_needed(bulk, weight_field, pm, imbalance_threshold);

  if( 0 == p_rank )
    std::cerr << "imbalance_threshold after rebalance = " << imbalance_threshold << ", " << do_rebal << std::endl;

  stk::mesh::Selector owned_selector = meta.locally_owned_part();
  size_t num_local_elems = stk::mesh::count_selected_entities(owned_selector, bulk.buckets(element_rank));

  //std::cerr << "(proc " << p_rank << ") : has " << num_local_elems << " local elems after rebalance." << std::endl;
  //stk::mesh::EntityVector local_elems;
  //stk::mesh::get_selected_entities(owned_selector, bulk.buckets(element_rank), local_elems);
  //double weight_sum = 0.0;
  //for( size_t ie = 0; ie < local_elems.size(); ++ie )
  //{
  //  stk::mesh::Entity * elem = local_elems[ie];
  //  double * const e_weight = stk::mesh::field_data( weight_field , *elem );
  //  //std::cerr << "(proc " << p_rank << ") : elem " << elem->identifier() << " has weight = " << *e_weight << std::endl;
  //  weight_sum += *e_weight;
  //}
  //double imbalance = weight_sum / total_weight * (double)p_size;
  //std::cerr << "(proc " << p_rank << ") : has imbalance = " << imbalance << std::endl;

  // Check that we satisfy our threshhold
  bool result = true;
  if( 4 > p_size )
  {
    result = (fabs(imbalance_threshold - 1.0) < 1.e-8);
    result = result & (num_local_elems == p_rank+1);
  }
  else
  {
    // Would like to put something here, but Zoltan using its default algorithm (RCB)
    // isn't able to do an adequate job rebalancing
    result = !do_rebal;
  }

  return result;
}

} //namespace use_cases
} //namespace rebalance
} //namespace stk


