/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <use_cases/UseCase_Rebal_3.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
 
#include <stk_mesh/fem/DefaultFEM.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/CreateAdjacentEntities.hpp>

#include <stk_mesh/fixtures/HexFixture.hpp>

#include <stk_rebalance/Rebalance.hpp>
#include <stk_rebalance/Partition.hpp>
#include <stk_rebalance/ZoltanPartition.hpp>

#include <stk_rebalance_utils/RebalanceUtils.hpp>

//----------------------------------------------------------------------

using namespace stk::mesh::fixtures;

typedef stk::mesh::Field<double> ScalarField ;

namespace stk {
namespace rebalance {
namespace use_cases {

class Test_Case_3_Partition : public stk::rebalance::Zoltan {
  public :
  explicit Test_Case_3_Partition(ParallelMachine pm,
                                 const unsigned ndim,
                                 Teuchos::ParameterList & rebal_region_parameters,
                                 std::string parameters_name=default_parameters_name());
  virtual ~Test_Case_3_Partition();
  virtual bool partition_dependents_needed() const; 
};

Test_Case_3_Partition::Test_Case_3_Partition(ParallelMachine pm,
                                             const unsigned ndim,
                                             Teuchos::ParameterList & rebal_region_parameters,
                                             std::string parameters_name) :
  stk::rebalance::Zoltan(pm, ndim, rebal_region_parameters, parameters_name) {}
Test_Case_3_Partition::~Test_Case_3_Partition() {}
bool Test_Case_3_Partition::partition_dependents_needed() const 
{return false;} //Do NOT move dependent entities for this case

enum { nx = 4, ny = 4 };

bool test_contact_surfaces( stk::ParallelMachine comm )
{
  unsigned spatial_dimension = 2;
  std::vector<std::string> rank_names = stk::mesh::fem::entity_rank_names(spatial_dimension);
  const stk::mesh::EntityRank constraint_rank = rank_names.size();
  rank_names.push_back("Constraint");

  stk::mesh::MetaData meta_data( rank_names );
  stk::mesh::BulkData bulk_data( meta_data , comm , 100 );
  stk::mesh::DefaultFEM top_data( meta_data, spatial_dimension );
  const stk::mesh::EntityRank element_rank    = stk::mesh::fem::element_rank(top_data);
  const stk::mesh::EntityRank side_rank       = stk::mesh::fem::side_rank(top_data);
  const stk::mesh::EntityRank node_rank       = stk::mesh::fem::node_rank(top_data);

  stk::mesh::Part & quad_part( stk::mesh::declare_part<shards::Quadrilateral<4> >( meta_data, "quad" ) );
  stk::mesh::Part & side_part( stk::mesh::declare_part<shards::Line<2> >         ( meta_data, "line" ) );
  VectorField & coord_field( meta_data.declare_field< VectorField >( "coordinates" ) );
  ScalarField & weight_field( meta_data.declare_field< ScalarField >( "element_weights" ) );

  stk::mesh::put_field( coord_field , node_rank , meta_data.universal_part() );
  stk::mesh::put_field(weight_field , element_rank , meta_data.universal_part() );

  meta_data.commit();

  //const unsigned p_size = bulk_data.parallel_size();
  const unsigned p_rank = bulk_data.parallel_rank();

  bulk_data.modification_begin();

  if ( !p_rank ) { 

    std::vector<stk::mesh::EntityVector> quads(nx);
    for ( unsigned ix = 0 ; ix < nx ; ++ix ) quads[ix].resize(ny);

    const unsigned nnx = nx + 1 ; 
    const unsigned nny = ny + 1 ; 
    unsigned face_id   = 1;
    for ( unsigned iy = 0 ; iy < ny ; ++iy ) { 
      for ( unsigned ix = 0 ; ix < nx ; ++ix ) { 
        stk::mesh::EntityId elem = 1 + ix + iy * nx ;
        stk::mesh::EntityId nodes[4] ;
        nodes[0] = 1 + ix + iy * nnx ;
        nodes[1] = 2 + ix + iy * nnx ;
        nodes[2] = 2 + ix + ( iy + 1 ) * nnx ;
        nodes[3] = 1 + ix + ( iy + 1 ) * nnx ;

        stk::mesh::Entity &q = stk::mesh::declare_element( bulk_data , quad_part , elem , nodes );
        if (0==ix) {
          stk::mesh::declare_element_side( bulk_data, face_id, q, 0 /*local side id*/, &side_part);
          ++face_id;
        }
        quads[ix][iy] = &q; 
      }   
    }   

    for ( unsigned iy = 0 ; iy < ny ; ++iy ) { 
      for ( unsigned ix = 0 ; ix < nx ; ++ix ) { 
        stk::mesh::EntityId elem = 1 + ix + iy * nx ;
        stk::mesh::Entity * e = bulk_data.get_entity( element_rank, elem );
        double * const e_weight = stk::mesh::field_data( weight_field , *e );
        *e_weight = 1.0;
      }   
    }   
    for ( unsigned iy = 0 ; iy <= ny ; ++iy ) { 
      for ( unsigned ix = 0 ; ix <= nx ; ++ix ) { 
        stk::mesh::EntityId nid = 1 + ix + iy * nnx ;
        stk::mesh::Entity * n = bulk_data.get_entity( node_rank, nid );
        double * const coord = stk::mesh::field_data( coord_field , *n );
        coord[0] = .1*ix;
        coord[1] = .1*iy;
        coord[2] = 0;
      }
    }

    // Assign constraint relations between nodes at top and bottom of mesh
    {
      const unsigned iy_bottom  =  0;
      const unsigned iy_top = ny;
      stk::mesh::PartVector add(1, &meta_data.locally_owned_part());
      for ( unsigned ix = 0 ; ix <= nx ; ++ix ) {
        stk::mesh::EntityId nid_bottom  = 1 + ix + iy_bottom  * nnx ;
        stk::mesh::EntityId nid_top = 1 + ix + iy_top * nnx ;
        stk::mesh::Entity * n_bottom  = bulk_data.get_entity( node_rank, nid_bottom  );
        stk::mesh::Entity * n_top = bulk_data.get_entity( node_rank, nid_top );
        const stk::mesh::EntityId constraint_entity_id =  1 + ix + nny * nnx;
        stk::mesh::Entity & c = bulk_data.declare_entity( constraint_rank, constraint_entity_id, add );
        bulk_data.declare_relation( c , *n_bottom  , 0 );
        bulk_data.declare_relation( c , *n_top , 1 );
      }
    } // end snippet

  }

  bulk_data.modification_end();

  // Force a rebalance by using imbalance_threshold < 1.0
  stk::mesh::Selector selector(side_part);
  selector &=  meta_data.locally_owned_part();
  double imbalance_threshold = 0.5;
  bool do_rebal = stk::rebalance::rebalance_needed(bulk_data, NULL, imbalance_threshold, side_rank, &selector);
  // Coordinates are passed to support geometric-based load balancing algorithms
  if( do_rebal ) {
    // Zoltan partition is specialized form a virtual base class, stk::rebalance::Partition.
    // Other specializations are possible.
    Teuchos::ParameterList emptyList;
    stk::rebalance::use_cases::Test_Case_3_Partition zoltan_partition(comm, spatial_dimension, emptyList);
    stk::rebalance::rebalance(bulk_data, selector, &coord_field, NULL, zoltan_partition, side_rank);
  }

  imbalance_threshold = 1.5;
  do_rebal = stk::rebalance::rebalance_needed(bulk_data, NULL, imbalance_threshold, side_rank, &selector);

  if( !p_rank )
    std::cerr << std::endl 
     << "imbalance_threshold after rebalance = " << imbalance_threshold <<", "<<do_rebal << std::endl;

  // Check that we satisfy our threshhold
  bool result = !do_rebal ;

  return result;
}

} //namespace use_cases
} //namespace rebalance
} //namespace stk


