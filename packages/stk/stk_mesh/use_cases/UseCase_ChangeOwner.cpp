#include <iostream>
#include <set>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/EntityComm.hpp>

#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/DefaultFEM.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyTraits.hpp>

#include "UseCase_ChangeOwner.hpp"

namespace {

using stk::mesh::fem::NODE_RANK;

const unsigned spatial_dimension = 2;

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
  : m_meta_data( stk::mesh::fem::entity_rank_names(spatial_dimension) ),
    m_bulk_data( m_meta_data , comm , 100 ),
    m_topo_data( m_meta_data, spatial_dimension ),
    m_quad_part( stk::mesh::declare_part<shards::Quadrilateral<4> >( m_meta_data, "quad")),
    m_coord_field( m_meta_data.declare_field< VectorField >( "coordinates" ) ),
    m_elem_rank( stk::mesh::fem::element_rank(m_topo_data) )
{
  stk::mesh::put_field( m_coord_field , NODE_RANK , m_meta_data.universal_part() );

  m_meta_data.commit();
}

bool Grid2D_Fixture::test_change_owner( unsigned nx , unsigned ny )
{
  const unsigned p_size = m_bulk_data.parallel_size();
  const unsigned p_rank = m_bulk_data.parallel_rank();

  int result = 1 ;

  m_bulk_data.modification_begin();

  // First of all work out the node ids and declare element elem
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
    result = ! m_bulk_data.buckets( NODE_RANK ).empty() &&
      ! m_bulk_data.buckets( m_elem_rank ).empty();
  }
  else {
    result = m_bulk_data.buckets( NODE_RANK ).empty() &&
             m_bulk_data.buckets( m_elem_rank ).empty();
  }

  m_bulk_data.modification_end();

  if ( 1 < p_size ) {
    std::vector< stk::mesh::EntityProc > change ;

    if ( p_rank == 0 ) {
      if ( p_size == 3 ) {
        const unsigned nnx = nx + 1 ;
        const unsigned nny = ny + 1 ;
        for ( unsigned iy = nny / 2 ; iy < nny ; ++iy ) {
          for ( unsigned ix = 0 ; ix < nnx ; ++ix ) {
            stk::mesh::EntityId id = 1 + ix + iy * nnx ;
            unsigned proc = ix < nx/2 ? 1 : 2;
            stk::mesh::EntityProc tmp( m_bulk_data.get_entity( NODE_RANK , id ) , proc );
            change.push_back( tmp );
          }
        }
        for ( unsigned iy = ny / 2 ; iy < ny ; ++iy ) {
          for ( unsigned ix = 0 ; ix < nx ; ++ix ) {
            stk::mesh::EntityId id = 1 + ix + iy * nx ;
            unsigned proc = ix < nx/2 ? 1 : 2;
            stk::mesh::EntityProc tmp( m_bulk_data.get_entity( m_elem_rank , id ) , proc );
            change.push_back( tmp );
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
            change.push_back( tmp );
          }
        }
        for ( unsigned iy = ny / 2 ; iy < ny ; ++iy ) {
          for ( unsigned ix = 0 ; ix < nx ; ++ix ) {
            stk::mesh::EntityId id = 1 + ix + iy * nx ;
            stk::mesh::EntityProc tmp( m_bulk_data.get_entity( m_elem_rank , id ) , 1 );
            change.push_back( tmp );
          }
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
          stk::mesh::EntityProc tmp( m_bulk_data.get_entity( NODE_RANK , id ) , 1 );
          change.push_back( tmp );
        }
      }
      for ( unsigned iy = 0 ; iy < ny / 2 ; ++iy ) {
        for ( unsigned ix = 0 ; ix < nx ; ++ix ) {
          stk::mesh::EntityId id = 1 + ix + iy * nx ;
          stk::mesh::EntityProc tmp( m_bulk_data.get_entity( m_elem_rank , id ) , 1 );
          change.push_back( tmp );
        }
      }
    }

    m_bulk_data.modification_begin();
    m_bulk_data.change_entity_owner( change );
    m_bulk_data.modification_end();

    if ( p_size == 3 ) {
      change.clear();

      if ( p_rank == 2 ) {
        const unsigned nnx = nx + 1 ;
        const unsigned nny = ny + 1 ;
        for ( unsigned iy = nny / 2 ; iy < nny ; ++iy ) {
          for ( unsigned ix = nx / 2 ; ix < nnx ; ++ix ) {
            stk::mesh::EntityId id = 1 + ix + iy * nnx ;
            unsigned proc = 1;
            stk::mesh::EntityProc tmp( m_bulk_data.get_entity( NODE_RANK , id ) , proc );
            change.push_back( tmp );
          }
        }
        for ( unsigned iy = ny / 2 ; iy < ny ; ++iy ) {
          for ( unsigned ix = nx / 2 ; ix < nx ; ++ix ) {
            stk::mesh::EntityId id = 1 + ix + iy * nx ;
            unsigned proc = 1;
            stk::mesh::EntityProc tmp( m_bulk_data.get_entity( m_elem_rank , id ) , proc );
            change.push_back( tmp );
          }
        }
      }

      m_bulk_data.modification_begin();
      m_bulk_data.change_entity_owner( change );
      m_bulk_data.modification_end();
    }

    // Only P1 has any nodes or elements
    if ( p_rank == 1 ) {
      result = ! m_bulk_data.buckets( NODE_RANK ).empty() &&
               ! m_bulk_data.buckets( m_elem_rank ).empty();
    }
    else {
      result = m_bulk_data.buckets( NODE_RANK ).empty() &&
               m_bulk_data.buckets( m_elem_rank ).empty();
    }
  }

  stk::all_reduce( m_bulk_data.parallel() , stk::ReduceMin<1>( & result ) );

  return result ;
}

//----------------------------------------------------------------------------

bool test_change_owner_with_constraint( stk::ParallelMachine pm )
{
  bool success = true ;

  const unsigned p_rank = stk::parallel_machine_rank( pm );
  const unsigned p_size = stk::parallel_machine_size( pm );

  if ( p_size != 2 ) { return success ; }

  std::vector<std::string> rank_names = stk::mesh::fem::entity_rank_names(spatial_dimension);
  const stk::mesh::EntityRank constraint_rank = rank_names.size();
  rank_names.push_back("Constraint");
  stk::mesh::MetaData meta_data( rank_names );
  stk::mesh::DefaultFEM topo_data( meta_data, spatial_dimension );
  const stk::mesh::EntityRank element_rank = stk::mesh::fem::element_rank(topo_data);

  VectorField * coordinates_field =
    & put_field(
        meta_data.declare_field<VectorField>("coordinates"),
        NODE_RANK,
        meta_data.universal_part() ,
        3
        );

  stk::mesh::Part & owned_part = meta_data.locally_owned_part();
  stk::mesh::Part & quad_part  = stk::mesh::declare_part<shards::Quadrilateral<4> >( meta_data, "quad" );

  meta_data.commit();

  stk::mesh::BulkData bulk_data( meta_data, pm, 100 );
  bulk_data.modification_begin();

  unsigned nx = 3;
  unsigned ny = 3;

  if ( p_rank==0 )
  {
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

    for ( unsigned iy = 0 ; iy < ny+1 ; ++iy ) {
      for ( unsigned ix = 0 ; ix < nx+1 ; ++ix ) {
        stk::mesh::EntityId nid = 1 + ix + iy * nnx ;
        stk::mesh::Entity * n = bulk_data.get_entity( NODE_RANK, nid );
        double * const coord = stk::mesh::field_data( *coordinates_field , *n );
        coord[0] = .1*ix;
        coord[1] = .1*iy;
        coord[2] = 0;
      }
    }
  }

  bulk_data.modification_end();

  if ( p_size>1 )
  {
    std::vector<stk::mesh::EntityProc> ep;

    if ( p_rank==0 )
    {
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  3 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  4 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  7 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  8 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 11 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 12 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 15 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 16 ), 1 ) );

      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( element_rank, 3 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( element_rank, 6 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( element_rank, 9 ), 1 ) );
    }

    bulk_data.modification_begin();
    bulk_data.change_entity_owner( ep );
    bulk_data.modification_end();

    bulk_data.modification_begin();

    if ( p_rank==1 )
    {
      // create constraint

      stk::mesh::Entity * n10 = bulk_data.get_entity( NODE_RANK, 10 );
      stk::mesh::Entity * n11 = bulk_data.get_entity( NODE_RANK, 11 );
      stk::mesh::Entity * n12 = bulk_data.get_entity( NODE_RANK, 12 );

      stk::mesh::PartVector add;
      add.push_back( &owned_part );
      const stk::mesh::EntityId c_entity_id = 1;
      stk::mesh::Entity & c = bulk_data.declare_entity( constraint_rank, c_entity_id, add );
      bulk_data.declare_relation( c , *n10 , 0 );
      bulk_data.declare_relation( c , *n11 , 1 );
      bulk_data.declare_relation( c , *n12 , 2 );
    }

    bulk_data.modification_end();

    stk::mesh::Entity * n10 = bulk_data.get_entity( NODE_RANK, 10 );

    if ( p_rank==0 or p_rank==1 )
    {
      ThrowErrorMsgIf( !stk::mesh::in_shared( *n10 ), "NODE[10] not shared" );
    }

    bulk_data.modification_begin();

    if ( p_rank==1 )
    {
      // destroy constraint

      stk::mesh::Entity * c1 = bulk_data.get_entity( constraint_rank, 1 );

      ThrowErrorMsgIf( !bulk_data.destroy_entity( c1 ),
                       "failed to destroy constraint" );
    }

    bulk_data.modification_end();

    if ( p_rank==0 or p_rank==1 )
    {
      ThrowErrorMsgIf( stk::mesh::in_shared( *n10 ), "NODE[10] shared" );
    }
  }

  return success ;
}

//----------------------------------------------------------------------------

bool test_change_owner_2( stk::ParallelMachine pm )
{
  const unsigned p_rank = stk::parallel_machine_rank( pm );
  const unsigned p_size = stk::parallel_machine_size( pm );

  if ( p_size != 3 ) { return true ; }

  stk::mesh::MetaData meta_data( stk::mesh::fem::entity_rank_names( spatial_dimension ) );
  stk::mesh::DefaultFEM topo_data( meta_data, spatial_dimension );
  const stk::mesh::EntityRank element_rank = stk::mesh::fem::element_rank(topo_data);

  VectorField * coordinates_field =
    & put_field(
        meta_data.declare_field<VectorField>("coordinates"),
        NODE_RANK,
        meta_data.universal_part() ,
        3
        );

  stk::mesh::Part & quad_part  = stk::mesh::declare_part<shards::Quadrilateral<4> >( meta_data, "quad" );

  meta_data.commit();

  stk::mesh::BulkData bulk_data( meta_data, pm, 100 );
  bulk_data.modification_begin();

  unsigned nx = 3;
  unsigned ny = 3;

  if ( p_rank==0 )
  {
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

    for ( unsigned iy = 0 ; iy < ny+1 ; ++iy ) {
      for ( unsigned ix = 0 ; ix < nx+1 ; ++ix ) {
        stk::mesh::EntityId nid = 1 + ix + iy * nnx ;
        stk::mesh::Entity * n = bulk_data.get_entity( NODE_RANK, nid );
        double * const coord = stk::mesh::field_data( *coordinates_field , *n );
        coord[0] = .1*ix;
        coord[1] = .1*iy;
        coord[2] = 0;
      }
    }
  }

  bulk_data.modification_end();

  if ( p_size==3 )
  {
    std::vector<stk::mesh::EntityProc> ep;

    if ( p_rank==0 )
    {
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  9 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 13 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 14 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 15 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 16 ), 1 ) );

      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( element_rank, 7 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( element_rank, 8 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( element_rank, 9 ), 1 ) );

      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  3 ), 2 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  4 ), 2 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  7 ), 2 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  8 ), 2 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 11 ), 2 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 12 ), 2 ) );

      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( element_rank, 3 ), 2 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( element_rank, 6 ), 2 ) );
    }

    bulk_data.modification_begin();
    bulk_data.change_entity_owner( ep );
    bulk_data.modification_end();

    ep.clear();

    if ( p_rank==1 )
    {
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  9 ), 0 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 13 ), 0 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 14 ), 0 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 15 ), 0 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 16 ), 0 ) );

      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( element_rank, 7 ), 0 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( element_rank, 8 ), 0 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( element_rank, 9 ), 0 ) );
    }

    if ( p_rank==2 )
    {
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  3 ), 0 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  4 ), 0 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  7 ), 0 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  8 ), 0 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 11 ), 0 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 12 ), 0 ) );

      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( element_rank, 3 ), 0 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( element_rank, 6 ), 0 ) );
    }

    bulk_data.modification_begin();
    bulk_data.change_entity_owner( ep );
    bulk_data.modification_end();
  }

  return true ;
}

//----------------------------------------------------------------------------

bool test_change_owner_3( stk::ParallelMachine pm )
{
  const int p_rank = stk::parallel_machine_rank( pm );
  const unsigned p_size = stk::parallel_machine_size( pm );

  stk::mesh::MetaData meta_data( stk::mesh::fem::entity_rank_names(spatial_dimension) );
  stk::mesh::DefaultFEM topo_data( meta_data, spatial_dimension );
  const stk::mesh::EntityRank element_rank = stk::mesh::fem::element_rank(topo_data);

  VectorField * coordinates_field =
    & put_field(
        meta_data.declare_field<VectorField>("coordinates"),
        NODE_RANK,
        meta_data.universal_part() ,
        3
        );

  stk::mesh::Part & quad_part  = stk::mesh::declare_part<shards::Quadrilateral<4> >( meta_data, "quad" );

  meta_data.commit();

  stk::mesh::BulkData bulk_data( meta_data, pm, 100 );
  bulk_data.modification_begin();

  unsigned nx = 3;
  unsigned ny = 3;

  if ( p_rank==0 )
  {
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

    for ( unsigned iy = 0 ; iy < ny+1 ; ++iy ) {
      for ( unsigned ix = 0 ; ix < nx+1 ; ++ix ) {
        stk::mesh::EntityId nid = 1 + ix + iy * nnx ;
        stk::mesh::Entity * n = bulk_data.get_entity( NODE_RANK, nid );
        double * const coord = stk::mesh::field_data( *coordinates_field , *n );
        coord[0] = .1*ix;
        coord[1] = .1*iy;
        coord[2] = 0;
      }
    }
  }

  bulk_data.modification_end();

  if ( p_size>1 ) {

    std::vector<stk::mesh::EntityProc> ep;

    if ( p_rank==0 )
    {
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  2 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  3 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  4 ), 1 ) );

      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  6 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  7 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  8 ), 1 ) );

      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 10 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 11 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 12 ), 1 ) );

      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 14 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 15 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 16 ), 1 ) );

      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( element_rank, 2 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( element_rank, 5 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( element_rank, 8 ), 1 ) );

      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( element_rank, 3 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( element_rank, 6 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( element_rank, 9 ), 1 ) );
    }

    bulk_data.modification_begin();
    bulk_data.change_entity_owner( ep );
    bulk_data.modification_end();

    // output to debug

    ep.clear();

    if ( p_rank==0 )
    {
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  1 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK,  5 ), 1 ) );

      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( element_rank, 1 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( element_rank, 4 ), 1 ) );
    }
    else if ( p_rank==1 )
    {
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 10 ), 0 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 11 ), 0 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 12 ), 0 ) );

      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 14 ), 0 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 15 ), 0 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( NODE_RANK, 16 ), 0 ) );

      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( element_rank, 8 ), 0 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( element_rank, 9 ), 0 ) );
    }

    bulk_data.modification_begin();
    bulk_data.change_entity_owner( ep );
    bulk_data.modification_end();
  }

  return true ;
}


