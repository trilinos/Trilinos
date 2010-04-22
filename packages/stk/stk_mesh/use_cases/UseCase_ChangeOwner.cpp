
#include <iostream>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/EntityComm.hpp>

#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/FieldDeclarations.hpp>
#include <stk_mesh/fem/FieldTraits.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyTraits.hpp>

#include "UseCase_ChangeOwner.hpp"

namespace {

void print_entity_status( std::ostream & os , stk::mesh::BulkData & mesh , unsigned rank )
{
  std::vector< stk::mesh::Entity * > entities ;
  stk::mesh::get_entities( mesh , rank , entities );
  for ( std::vector< stk::mesh::Entity * >::iterator
        i = entities.begin() ; i != entities.end() ; ++i ) {
    stk::mesh::Entity & e = **i ;
    stk::mesh::print_entity_key( os , mesh.mesh_meta_data() , e.key() );
    if ( stk::mesh::in_receive_ghost( e ) ) {
      os << " IS GHOST" ;
    }
    os << std::endl ;
  }
}

}

Grid2D_Fixture::Grid2D_Fixture( stk::ParallelMachine comm )
  : m_meta_data( stk::mesh::fem_entity_rank_names() ),
    m_bulk_data( m_meta_data , comm , 100 ),
    m_quad_part( m_meta_data.declare_part( "quad" , stk::mesh::Element ) ),
    m_coord_field( m_meta_data.declare_field< stk::mesh::VectorField >( "coordinates" ) )
{
  stk::mesh::put_field( m_coord_field , stk::mesh::Node , m_meta_data.universal_part() );
  stk::mesh::set_cell_topology< shards::Quadrilateral<4> >( m_quad_part );

  m_meta_data.commit();
}

bool Grid2D_Fixture::test_change_owner( unsigned nx , unsigned ny )
{
  const unsigned p_size = m_bulk_data.parallel_size();
  const unsigned p_rank = m_bulk_data.parallel_rank();

  int result = 1 ;

  m_bulk_data.modification_begin();

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

        stk::mesh::declare_element( m_bulk_data , m_quad_part , elem , nodes );
      }
    }
  }

  // Only P0 has any nodes or elements
  if ( p_rank == 0 ) {
    result = ! m_bulk_data.buckets( 0 ).empty() &&
             ! m_bulk_data.buckets( 3 ).empty();
  }
  else {
    result = m_bulk_data.buckets( 0 ).empty() &&
             m_bulk_data.buckets( 3 ).empty();
  }

  m_bulk_data.modification_end();

  if ( 1 < p_size ) {
    std::vector< stk::mesh::EntityProc > change ;

    if ( p_rank == 0 ) {
      const unsigned nnx = nx + 1 ;
      const unsigned nny = ny + 1 ;
      for ( unsigned iy = nny / 2 ; iy < nny ; ++iy ) {
        for ( unsigned ix = 0 ; ix < nnx ; ++ix ) {
          stk::mesh::EntityId id = 1 + ix + iy * nnx ;
          stk::mesh::EntityProc tmp( m_bulk_data.get_entity( stk::mesh::Node , id ) , 1 );
          change.push_back( tmp );
        }
      }
      for ( unsigned iy = ny / 2 ; iy < ny ; ++iy ) {
        for ( unsigned ix = 0 ; ix < nx ; ++ix ) {
          stk::mesh::EntityId id = 1 + ix + iy * nx ;
          stk::mesh::EntityProc tmp( m_bulk_data.get_entity( stk::mesh::Element , id ) , 1 );
          change.push_back( tmp );
        }
      }
    }

    m_bulk_data.modification_begin();
    m_bulk_data.change_entity_owner( change );
    m_bulk_data.modification_end();

    change.clear();

    if ( p_rank == 0 ) {
      const unsigned nnx = nx + 1 ;
      const unsigned nny = ny + 1 ;
      for ( unsigned iy = 0 ; iy < nny / 2 ; ++iy ) {
        for ( unsigned ix = 0 ; ix < nnx ; ++ix ) {
          stk::mesh::EntityId id = 1 + ix + iy * nnx ;
          stk::mesh::EntityProc tmp( m_bulk_data.get_entity( stk::mesh::Node , id ) , 1 );
          change.push_back( tmp );
        }
      }
      for ( unsigned iy = 0 ; iy < ny / 2 ; ++iy ) {
        for ( unsigned ix = 0 ; ix < nx ; ++ix ) {
          stk::mesh::EntityId id = 1 + ix + iy * nx ;
          stk::mesh::EntityProc tmp( m_bulk_data.get_entity( stk::mesh::Element , id ) , 1 );
          change.push_back( tmp );
        }
      }
    }

    m_bulk_data.modification_begin();
    m_bulk_data.change_entity_owner( change );
    m_bulk_data.modification_end();

    // Only P1 has any nodes or elements
    if ( p_rank == 1 ) {
      result = ! m_bulk_data.buckets( 0 ).empty() &&
               ! m_bulk_data.buckets( 3 ).empty();
    }
    else {
      result = m_bulk_data.buckets( 0 ).empty() &&
               m_bulk_data.buckets( 3 ).empty();
    }
  }

  stk::all_reduce( m_bulk_data.parallel() , stk::ReduceMin<1>( & result ) );

  return result ;
}

