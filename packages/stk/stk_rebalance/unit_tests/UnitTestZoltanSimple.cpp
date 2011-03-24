/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_mesh/fem/CreateAdjacentEntities.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>

#include <stk_rebalance/Rebalance.hpp>
#include <stk_rebalance/Partition.hpp>
#include <stk_rebalance/ZoltanPartition.hpp>

#include <stk_rebalance_utils/RebalanceUtils.hpp>

using stk::mesh::fem::NODE_RANK;

typedef stk::mesh::Field<double> ScalarField ;
typedef stk::mesh::Field<double, stk::mesh::Cartesian> VectorField ;

enum { nx = 2, ny = 2 };

STKUNIT_UNIT_TEST(UnitTestZoltanSimple, testUnit)
{
#ifdef STK_HAS_MPI
  stk::ParallelMachine comm(MPI_COMM_WORLD);
#else
  stk::ParallelMachine comm(0);
#endif

  unsigned spatial_dimension = 2;
  std::vector<std::string> rank_names = stk::mesh::fem::entity_rank_names(spatial_dimension);
  const stk::mesh::EntityRank constraint_rank = rank_names.size();
  rank_names.push_back("Constraint");

  stk::mesh::fem::FEMMetaData fem_meta;
  fem_meta.FEM_initialize( spatial_dimension, rank_names );
  stk::mesh::MetaData & meta_data = stk::mesh::fem::FEMMetaData::get_meta_data(fem_meta);
  stk::mesh::BulkData bulk_data( meta_data , comm , 100 );
  const stk::mesh::EntityRank element_rank    = fem_meta.element_rank();

  stk::mesh::fem::CellTopology quad_top(shards::getCellTopologyData<shards::Quadrilateral<4> >());
  stk::mesh::Part & quad_part( fem_meta.declare_part("quad", quad_top ) );
  VectorField & coord_field( fem_meta.declare_field< VectorField >( "coordinates" ) );
  ScalarField & weight_field( fem_meta.declare_field< ScalarField >( "element_weights" ) );

  stk::mesh::put_field( coord_field , NODE_RANK , fem_meta.universal_part() );
  stk::mesh::put_field(weight_field , element_rank , fem_meta.universal_part() );

  fem_meta.commit();

  const unsigned p_size = bulk_data.parallel_size();
  const unsigned p_rank = bulk_data.parallel_rank();

  bulk_data.modification_begin();

  if ( p_rank == 0 ) {

    std::vector<std::vector<stk::mesh::Entity*> > quads(nx);
    for ( unsigned ix = 0 ; ix < nx ; ++ix ) quads[ix].resize(ny);

    const unsigned nnx = nx + 1 ;
    const unsigned nny = ny + 1 ;
    for ( unsigned iy = 0 ; iy < ny ; ++iy ) {
      for ( unsigned ix = 0 ; ix < nx ; ++ix ) {
        stk::mesh::EntityId elem = 1 + ix + iy * nx ;
        stk::mesh::EntityId nodes[4] ;
        nodes[0] = 1 + ix + iy * nnx ;
        nodes[1] = 2 + ix + iy * nnx ;
        nodes[2] = 2 + ix + ( iy + 1 ) * nnx ;
        nodes[3] = 1 + ix + ( iy + 1 ) * nnx ;

        stk::mesh::Entity &q = stk::mesh::fem::declare_element( bulk_data , quad_part , elem , nodes );
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
        stk::mesh::Entity * n = bulk_data.get_entity( NODE_RANK, nid );
        double * const coord = stk::mesh::field_data( coord_field , *n );
        coord[0] = .1*ix;
        coord[1] = .1*iy;
        coord[2] = 0;
      }
    }

    {
      const unsigned iy_left  =  0; 
      const unsigned iy_right = ny; 
      stk::mesh::PartVector add(1, &fem_meta.locally_owned_part());
      for ( unsigned ix = 0 ; ix <= nx ; ++ix ) {
        stk::mesh::EntityId nid_left  = 1 + ix + iy_left  * nnx ;
        stk::mesh::EntityId nid_right = 1 + ix + iy_right * nnx ;
        stk::mesh::Entity * n_left  = bulk_data.get_entity( NODE_RANK, nid_left  );
        stk::mesh::Entity * n_right = bulk_data.get_entity( NODE_RANK, nid_right );
        const stk::mesh::EntityId constraint_entity_id =  1 + ix + nny * nnx;
        stk::mesh::Entity & c = bulk_data.declare_entity( constraint_rank, constraint_entity_id, add );
        bulk_data.declare_relation( c , *n_left  , 0 );
        bulk_data.declare_relation( c , *n_right , 1 );
      }
    }

  }

  // Only P0 has any nodes or elements
  if ( p_rank == 0 ) {
    STKUNIT_ASSERT( ! bulk_data.buckets( NODE_RANK ).empty() );
    STKUNIT_ASSERT( ! bulk_data.buckets( element_rank ).empty() );
  }
  else {
    STKUNIT_ASSERT( bulk_data.buckets( NODE_RANK ).empty() );
    STKUNIT_ASSERT( bulk_data.buckets( element_rank ).empty() );
  }


  bulk_data.modification_end();

  // create some sides and faces to rebalance.
  stk::mesh::PartVector add_parts;
  stk::mesh::create_adjacent_entities(bulk_data, add_parts);

  // Zoltan partition is specialized fomm a virtual base class, stk::rebalance::Partition.
  // Other specializations are possible.
  Teuchos::ParameterList emptyList;
  stk::rebalance::Zoltan zoltan_partition(comm, spatial_dimension, emptyList);

  {
    stk::mesh::Selector selector(fem_meta.universal_part());

    stk::rebalance::rebalance(bulk_data, selector, &coord_field, &weight_field, zoltan_partition);
  }

  const double imbalance_threshold = stk::rebalance::check_balance(bulk_data, &weight_field, element_rank);
  const bool do_rebal = 1.5 < imbalance_threshold;

  // Check that we satisfy our threshhold
  STKUNIT_ASSERT( !do_rebal );
  if( (2 == p_size) || (4 == p_size) )
  {
    STKUNIT_ASSERT_NEAR(imbalance_threshold, 1.0, 1.e-8);
  }
  else
  {
    STKUNIT_ASSERT_LE(imbalance_threshold, 1.5);
  }

  // And verify that all dependent entities are on the same proc as their parent element
  {
    stk::mesh::EntityVector entities;
    stk::mesh::Selector selector = fem_meta.locally_owned_part();

    get_selected_entities(selector, bulk_data.buckets(NODE_RANK), entities);
    bool result = stk::rebalance::verify_dependent_ownership(element_rank, entities);
    //get_selected_entities(selector, bulk_data.buckets(constraint_rank), entities);
    //result &= stk::rebalance::verify_dependent_ownership(element_rank, entities);
    STKUNIT_ASSERT( result );
  }
}

/// \page stk_rebalance_unit_test_zoltan
///  \ingroup stk_rebalance_unit_test_module
/// \section stk_rebalance_unit_test_zoltan_description Simple Zoltan Unit Test
///
/// This unit test creates a 2D quad mesh on proc 0 with coordinates and 
/// Parts associated with edges, nodes, and constraints and then moves these entities
/// as determined by calling Zoltan's load balancing capability using
/// default settings. 
///
/// Using the following mesh of 4 quads,
///
///  Global node and element numbering
/// <pre>
///
///   7       8       9
///   +-------+-------+
///   |       |       |
///   |  e3   |  e4   |
///   |       |5      |
///  4+-------+-------+6   
///   |       |       |     Y
///   |  e1   |  e2   |     |
///   |       |       |     |
///   +-------+-------+     *--> X
///   1       2      3 
/// </pre>
///
///  Local node numbering
///
/// <pre>
///     
///   3       4
///   +-------+
///   |       |
///   |  e1   |
///   |       |
///   +-------+
///   1       2
/// </pre>
///
/// Use of Zoltan with default settings is achieved by instantiating the
/// appropriate Partition class with an empty parameter list as follows,
/// \dontinclude UnitTestZoltanSimple.cpp
/// \skip Zoltan partition is 
/// \until zoltan_partition(
///
/// An initial assessment of imbalance is made using an element weight field
/// followed by a call to actually do the rebalance as follows,
/// \skip Force a rebalance
/// \until rebalance::rebalance(
/// 
/// Perfect balancing should result using 2 or 4 procs, and
/// on 3 procs, the imbalance threshold should be below 1.5.
/// The test passes if these criteria are satisfied.
///
/// See \ref UnitTestZoltanSimple.cpp for the complete source listing.
///
