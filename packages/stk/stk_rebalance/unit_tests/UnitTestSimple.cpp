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

#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/DefaultFEM.hpp>

#include <stk_rebalance/Rebalance.hpp>
#include <stk_rebalance/Partition.hpp>

using stk::mesh::fem::NODE_RANK;

typedef stk::mesh::Field<double> ScalarField ;
typedef stk::mesh::Field<double, stk::mesh::Cartesian> VectorField ;

static const int spatial_dimension = 2;

enum { nx = 2, ny = 2 };

class MockPartition : public stk::rebalance::Partition
{
  public:

    enum BALANCE_TEST_STEP
      { FIRST,
        SECOND,
        THIRD   };

  MockPartition( stk::mesh::MetaData & md, stk::mesh::BulkData & bd ) :
      stk::rebalance::Partition(bd.parallel()),
      m_meta_data(md),
      m_bulk_data(bd),
      m_fem(stk::mesh::fem::get_fem_interface(m_meta_data)),
      m_step(FIRST)
    {  }

    ~MockPartition() { }

    void set_balance_step(BALANCE_TEST_STEP step)
    { m_step = step; }

    void set_mesh_info ( const std::vector<stk::mesh::Entity *> &mesh_entities,
                         const VectorField   * nodal_coord_ref,
                         const ScalarField   * elem_weight_ref)
    { total_number_entities_ = mesh_entities.size(); }

  unsigned num_elems() const
  { return total_number_entities_; }

    void determine_new_partition(bool &RebalancingNeeded)
    { RebalancingNeeded = (m_bulk_data.parallel_size() > 1); }

    int get_new_partition(std::vector<stk::mesh::EntityProc> &new_partition);

    bool partition_dependents_needed() const
    { return false; /* I handle both element and dependent node partitioning */ }

  private:

    unsigned total_number_entities_;
    stk::mesh::MetaData & m_meta_data;
    stk::mesh::BulkData & m_bulk_data;
    stk::mesh::fem::FEMInterface & m_fem;
    BALANCE_TEST_STEP m_step;
};

int
MockPartition::get_new_partition(std::vector<stk::mesh::EntityProc> &new_partition)
{
  const unsigned p_size = m_bulk_data.parallel_size();
  const unsigned p_rank = m_bulk_data.parallel_rank();
  const stk::mesh::EntityRank element_rank = stk::mesh::fem::element_rank(m_fem);

  new_partition.clear();

  if ( 1 < p_size ) {

    if( FIRST == m_step )
    {
      if ( p_rank == 0 ) {
        if ( p_size == 3 ) {
          const unsigned nnx = nx + 1 ;
          const unsigned nny = ny + 1 ;
          for ( unsigned iy = nny / 2 ; iy < nny ; ++iy ) {
            for ( unsigned ix = 0 ; ix < nnx ; ++ix ) {
              stk::mesh::EntityId id = 1 + ix + iy * nnx ;
              unsigned proc = ix < nx/2 ? 1 : 2;
              stk::mesh::EntityProc tmp( m_bulk_data.get_entity( NODE_RANK , id ) , proc );
              new_partition.push_back( tmp );
            }
          }
          for ( unsigned iy = ny / 2 ; iy < ny ; ++iy ) {
            for ( unsigned ix = 0 ; ix < nx ; ++ix ) {
              stk::mesh::EntityId id = 1 + ix + iy * nx ;
              unsigned proc = ix < nx/2 ? 1 : 2;
              stk::mesh::EntityProc tmp( m_bulk_data.get_entity( element_rank , id ) , proc );
              new_partition.push_back( tmp );
            }
          }
        }
        else
        {
          const unsigned nnx = nx + 1 ;
          const unsigned nny = ny + 1 ;
          for ( unsigned iy = nny / 2 ; iy < nny ; ++iy ) {
            for ( unsigned ix = 0 ; ix < nnx ; ++ix ) {
              stk::mesh::EntityId id = 1 + ix + iy * nnx ;
              stk::mesh::EntityProc tmp( m_bulk_data.get_entity( NODE_RANK , id ) , 1 );
              new_partition.push_back( tmp );
            }
          }
          for ( unsigned iy = ny / 2 ; iy < ny ; ++iy ) {
            for ( unsigned ix = 0 ; ix < nx ; ++ix ) {
              stk::mesh::EntityId id = 1 + ix + iy * nx ;
              stk::mesh::EntityProc tmp( m_bulk_data.get_entity( element_rank , id ) , 1 );
              new_partition.push_back( tmp );
            }
          }
        }
      }
    }

    else if( SECOND == m_step )
    {
      if ( p_rank == 0 ) {
        const unsigned nnx = nx + 1 ;
        const unsigned nny = ny + 1 ;
        for ( unsigned iy = 0 ; iy < nny / 2 ; ++iy ) {
          for ( unsigned ix = 0 ; ix < nnx ; ++ix ) {
            stk::mesh::EntityId id = 1 + ix + iy * nnx ;
            stk::mesh::EntityProc tmp( m_bulk_data.get_entity( NODE_RANK , id ) , 1 );
            new_partition.push_back( tmp );
          }
        }
        for ( unsigned iy = 0 ; iy < ny / 2 ; ++iy ) {
          for ( unsigned ix = 0 ; ix < nx ; ++ix ) {
            stk::mesh::EntityId id = 1 + ix + iy * nx ;
            stk::mesh::EntityProc tmp( m_bulk_data.get_entity( element_rank , id ) , 1 );
            new_partition.push_back( tmp );
          }
        }
      }
    }

    else if( THIRD == m_step )
    {
      if ( p_size == 3 ) {
        new_partition.clear();

        if ( p_rank == 2 ) {
          const unsigned nnx = nx + 1 ;
          const unsigned nny = ny + 1 ;
          for ( unsigned iy = nny / 2 ; iy < nny ; ++iy ) {
            for ( unsigned ix = nx / 2 ; ix < nnx ; ++ix ) {
              stk::mesh::EntityId id = 1 + ix + iy * nnx ;
              unsigned proc = 1;
              stk::mesh::EntityProc tmp( m_bulk_data.get_entity( NODE_RANK , id ) , proc );
              new_partition.push_back( tmp );
            }
          }
          for ( unsigned iy = ny / 2 ; iy < ny ; ++iy ) {
            for ( unsigned ix = nx / 2 ; ix < nx ; ++ix ) {
              stk::mesh::EntityId id = 1 + ix + iy * nx ;
              unsigned proc = 1;
              stk::mesh::EntityProc tmp( m_bulk_data.get_entity( element_rank , id ) , proc );
              new_partition.push_back( tmp );
            }
          }
        }
      }
    }
  }

  return 0;
}


STKUNIT_UNIT_TEST(UnitTestRebalanceSimple, testUnit)
{
#ifdef STK_HAS_MPI
  stk::ParallelMachine comm(MPI_COMM_WORLD);
#else
  stk::ParallelMachine comm(0);
#endif

  unsigned spatial_dimension = 2;
  stk::mesh::MetaData meta_data( stk::mesh::fem::entity_rank_names(spatial_dimension) );
  stk::mesh::BulkData bulk_data( meta_data , comm , 100 );
  stk::mesh::DefaultFEM top_data( meta_data, spatial_dimension );
  const stk::mesh::EntityRank element_rank = stk::mesh::fem::element_rank(top_data);
  stk::mesh::Part & quad_part( stk::mesh::declare_part<shards::Quadrilateral<4> >( meta_data, "quad" ) );
  VectorField & coord_field( meta_data.declare_field< VectorField >( "coordinates" ) );
  ScalarField & weight_field( meta_data.declare_field< ScalarField >( "element_weights" ) );

  stk::mesh::put_field( coord_field , NODE_RANK , meta_data.universal_part() );
  stk::mesh::put_field(weight_field , element_rank , meta_data.universal_part() );

  meta_data.commit();

  const unsigned p_size = bulk_data.parallel_size();
  const unsigned p_rank = bulk_data.parallel_rank();

  bulk_data.modification_begin();

  if ( p_rank == 0 ) {
    const unsigned nnx = nx + 1 ;
    for ( unsigned iy = 0 ; iy < ny ; ++iy ) {
      for ( unsigned ix = 0 ; ix < nx ; ++ix ) {
        stk::mesh::EntityId elem = 1 + ix + iy * nx ;
        stk::mesh::EntityId nodes[4] ;
        nodes[0] = 1 + ix + iy * nnx ;
        nodes[1] = 2 + ix + iy * nnx ;
        nodes[2] = 2 + ix + ( iy + 1 ) * nnx ;
        nodes[3] = 1 + ix + ( iy + 1 ) * nnx ;

        stk::mesh::declare_element( bulk_data , quad_part , elem , nodes );
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

  MockPartition partition(meta_data, bulk_data);
  stk::mesh::Selector selector(meta_data.universal_part());

  partition.set_balance_step(MockPartition::FIRST);
  // Exercise the threshhold calculation by using imblance_threshhold > 1.0
  double imblance_threshhold = 1.5;
  bool do_rebal = stk::rebalance::rebalance_needed(bulk_data, &weight_field, comm, imblance_threshhold);
  if( do_rebal )
  {
    // Pick a few values as negative to exercise a check in rebalance::rebalance(...)
    // which converts negative weights to 1.0
    if ( p_rank == 0 )
    {
      for ( unsigned iy = 0 ; iy < ny ; ++iy )
      {
        stk::mesh::EntityId elem = 1 + iy * nx ;
        stk::mesh::Entity * e = bulk_data.get_entity( element_rank, elem );
        double * const e_weight = stk::mesh::field_data( weight_field , *e );
        *e_weight = -2.0;
      }
    }
    stk::rebalance::rebalance(bulk_data, selector, NULL, &weight_field, partition);
  }

  partition.set_balance_step(MockPartition::SECOND);
  // Force a rebalance by using imblance_threshhold < 1.0
  imblance_threshhold = 0.5;
  do_rebal = stk::rebalance::rebalance_needed(bulk_data, &weight_field, comm, imblance_threshhold);
  if( do_rebal )
    stk::rebalance::rebalance(bulk_data, selector, NULL, &weight_field, partition);

  partition.set_balance_step(MockPartition::THIRD);
  // Force a rebalance by using imblance_threshhold < 1.0
  imblance_threshhold = 0.5;
  do_rebal = stk::rebalance::rebalance_needed(bulk_data, &weight_field, comm, imblance_threshhold);
  if( do_rebal )
    stk::rebalance::rebalance(bulk_data, selector, NULL, &weight_field, partition);

  if ( 1 < p_size ) {
    // Only P1 has any nodes or elements
    if ( p_rank == 1 ) {
      STKUNIT_ASSERT( ! bulk_data.buckets( NODE_RANK ).empty() );
      STKUNIT_ASSERT( ! bulk_data.buckets( element_rank ).empty() );
    }
    else {
      STKUNIT_ASSERT( bulk_data.buckets( NODE_RANK ).empty() );
      STKUNIT_ASSERT( bulk_data.buckets( element_rank ).empty() );
    }
  }
}
