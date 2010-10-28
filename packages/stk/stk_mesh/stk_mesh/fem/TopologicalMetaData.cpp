/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <string>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <algorithm>

#include <Shards_BasicTopologies.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <stk_mesh/fem/TopologicalMetaData.hpp>

namespace stk {
namespace mesh {

void TopologicalMetaData::verify_spatial_dimension(
  unsigned spatial_dimension , const char * method )
{
  ThrowRequireMsg( spatial_dimension >= 1 && 3 >= spatial_dimension,
      method << " ERRONEOUS spatial dimension = " << spatial_dimension );
}

std::vector<std::string>
TopologicalMetaData::entity_rank_names( unsigned spatial_dimension )
{
  static const char method[] =
    "stk::mesh::TopologicalMetaData::entity_rank_names" ;

  verify_spatial_dimension( spatial_dimension , method );

  std::vector< std::string > names ;

  names.reserve( spatial_dimension + 1 );

  names.push_back( std::string( "NODE" ) );

  if ( 1 < spatial_dimension ) { names.push_back( std::string("EDGE") ); }
  if ( 2 < spatial_dimension ) { names.push_back( std::string("FACE") ); }

  names.push_back( std::string("ELEMENT") );
  names.push_back( std::string("PATCH") );

  return names ;
}

const TopologicalMetaData & 
TopologicalMetaData::find_TopologicalMetaData( const MetaData & meta_data )
{
  const TopologicalMetaData * top_data = TopologicalMetaData::internal_get(meta_data);
  if (NULL == top_data) {
    throw std::runtime_error("stk::mesh::TopologicalMetaData::find_TopologicalMetaData:  ERROR, TopologicalMetaData attribute not found on MetaData.  This is probably because no TopologicalMetaData was constructed for this MetaData.");
  }
  return *top_data;
}

const TopologicalMetaData *
TopologicalMetaData::internal_get( const MetaData & meta_data )
{
  return meta_data.get_attribute< TopologicalMetaData >();
}

TopologicalMetaData::~TopologicalMetaData()
{
  m_meta_data.remove_attribute( this );
}

TopologicalMetaData::TopologicalMetaData(
  MetaData & arg_meta_data ,
  unsigned   arg_spatial_dimension )
  : spatial_dimension( arg_spatial_dimension )
  , node_rank( 0 )
  , edge_rank( 1 < spatial_dimension ? 1 : 0 )
  , side_rank( 2 < spatial_dimension ? 2 : edge_rank )
  , element_rank( spatial_dimension )
  , patch_rank( spatial_dimension + 1 )
  , m_meta_data( arg_meta_data )
  , m_top_rank()
  , m_part_top_map()
{
  static const char method[] = "stk::mesh::TopologicalMetaData::TopologicalMetaData" ;

  verify_spatial_dimension( spatial_dimension , method );

  // Attach to the meta data
  m_meta_data.declare_attribute_no_delete< TopologicalMetaData >( this );

  // Load up appropriate standard cell topologies.

  declare_cell_topology( shards::getCellTopologyData< shards::Node >() , 0 );

  declare_cell_topology( shards::getCellTopologyData< shards::Line<2> >() , 1 );
  declare_cell_topology( shards::getCellTopologyData< shards::Line<3> >() , 1 );

  declare_cell_topology( shards::getCellTopologyData< shards::Particle >() , spatial_dimension );

  if ( 1 < spatial_dimension ) {
    declare_cell_topology( shards::getCellTopologyData< shards::Triangle<3> >() , 2 );
    declare_cell_topology( shards::getCellTopologyData< shards::Triangle<6> >() , 2 );
    declare_cell_topology( shards::getCellTopologyData< shards::Triangle<4> >() , 2 );

    declare_cell_topology( shards::getCellTopologyData< shards::Quadrilateral<4> >() , 2 );
    declare_cell_topology( shards::getCellTopologyData< shards::Quadrilateral<8> >() , 2 );
    declare_cell_topology( shards::getCellTopologyData< shards::Quadrilateral<9> >() , 2 );

    declare_cell_topology( shards::getCellTopologyData< shards::Beam<2> >() , spatial_dimension );
    declare_cell_topology( shards::getCellTopologyData< shards::Beam<3> >() , spatial_dimension );
  }

  if ( 2 == spatial_dimension ) {
    declare_cell_topology( shards::getCellTopologyData< shards::ShellLine<2> >() , 2 );
    declare_cell_topology( shards::getCellTopologyData< shards::ShellLine<3> >() , 2 );
  }

  if ( 2 < spatial_dimension ) {
    declare_cell_topology( shards::getCellTopologyData< shards::Tetrahedron<4> >() , 3 );
    declare_cell_topology( shards::getCellTopologyData< shards::Tetrahedron<10> >() , 3 );
    declare_cell_topology( shards::getCellTopologyData< shards::Tetrahedron<8> >() , 3 );

    declare_cell_topology( shards::getCellTopologyData< shards::Pyramid<5> >() , 3 );
    declare_cell_topology( shards::getCellTopologyData< shards::Pyramid<13> >() , 3 );
    declare_cell_topology( shards::getCellTopologyData< shards::Pyramid<14> >() , 3 );

    declare_cell_topology( shards::getCellTopologyData< shards::Wedge<6> >() , 3 );
    declare_cell_topology( shards::getCellTopologyData< shards::Wedge<15> >() , 3 );
    declare_cell_topology( shards::getCellTopologyData< shards::Wedge<18> >() , 3 );

    declare_cell_topology( shards::getCellTopologyData< shards::Hexahedron<8> >() , 3 );
    declare_cell_topology( shards::getCellTopologyData< shards::Hexahedron<20> >() , 3 );
    declare_cell_topology( shards::getCellTopologyData< shards::Hexahedron<27> >() , 3 );

    declare_cell_topology( shards::getCellTopologyData< shards::ShellTriangle<3> >() , 3 );
    declare_cell_topology( shards::getCellTopologyData< shards::ShellTriangle<6> >() , 3 );

    declare_cell_topology( shards::getCellTopologyData< shards::ShellQuadrilateral<4> >() , 3 );
    declare_cell_topology( shards::getCellTopologyData< shards::ShellQuadrilateral<8> >() , 3 );
    declare_cell_topology( shards::getCellTopologyData< shards::ShellQuadrilateral<9> >() , 3 );
  }
}

//----------------------------------------------------------------------------

EntityRank TopologicalMetaData::get_entity_rank( const CellTopologyData * top ) const
{
  typedef std::pair< const CellTopologyData * , EntityRank > ValueType ;

  std::vector< ValueType >::const_iterator i ;

  for ( i = m_top_rank.begin() ; i != m_top_rank.end() && top != i->first ; ++i );

  if (i == m_top_rank.end()) {
    throw std::runtime_error("stk::mesh::TopologicalMetaData::get_entity_rank ERROR, topology not found.");
  }

  return i->second ;
}

void TopologicalMetaData::internal_set_entity_rank( const CellTopologyData * top , EntityRank rank )
{
  typedef std::pair< const CellTopologyData * , EntityRank > ValueType ;

  m_top_rank.push_back( ValueType( top , rank ) );
}

//----------------------------------------------------------------------------

void TopologicalMetaData::declare_cell_topology(
  const CellTopologyData * top , EntityRank rank )
{
  static const char method[] = "stk::mesh::TopologicalMetaData::declare_cell_topology" ;

  typedef std::pair< const CellTopologyData * , EntityRank > ValueType ;

  std::vector< ValueType >::const_iterator i = m_top_rank.begin() ;
  for ( ; i != m_top_rank.end() && top != i->first ; ++i );

  const bool       duplicate     = i != m_top_rank.end();
  const EntityRank existing_rank = duplicate ? i->second : 0 ;

  const bool error_change = duplicate && existing_rank != rank ;
  const bool error_rank   = spatial_dimension < rank ;

  if ( error_rank || error_change ) {
    std::ostringstream msg ;
    msg << method << "( " << top->name
        << " , rank = " << rank << " ) ERROR " ;
    if ( error_rank ) {
      msg << ": rank exceeds maximum spatial_dimension = "
          << spatial_dimension ;
    }
    if ( error_change ) {
      msg << ": previously declared rank = " << existing_rank ;
    }
    throw std::runtime_error( msg.str() );
  }

  if ( ! duplicate ) {
    internal_set_entity_rank( top , rank );
  }
}

//----------------------------------------------------------------------------

Part & TopologicalMetaData::declare_part(
  const std::string & name , const CellTopologyData * top )
{
  //static const char method[] = "stk::mesh::TopologicalMetaData::declare_part" ;

  EntityRank entity_rank;
  try {
    entity_rank = get_entity_rank( top );
  }
  catch(std::runtime_error& /*x*/) {
    internal_set_entity_rank( top , top->dimension );
    entity_rank = top->dimension ;
  }

  Part & part = m_meta_data.declare_part( name , entity_rank );

  internal_set_cell_topology(part, entity_rank, top);
  return part ;
}

//----------------------------------------------------------------------------

void TopologicalMetaData::internal_set_cell_topology( 
    Part & part, 
    unsigned entity_rank, 
    const CellTopologyData * top
    ) 
{
  static const char method[] = "stk::mesh::TopologicalMetaData::internal_set_cell_topology" ;
  typedef std::pair< PartOrdinal , const CellTopologyData * > ValueType ;
  ValueType value( part.mesh_meta_data_ordinal() , top );

  std::vector< ValueType >::iterator
    i = std::lower_bound( m_part_top_map.begin() , m_part_top_map.end() , value );

  const bool duplicate  = i != m_part_top_map.end() && i->first == value.first ;
  const bool error_rank = part.primary_entity_rank() != entity_rank ;
  const bool error_top  = duplicate && i->second != value.second ;

  if ( error_rank || error_top ) {
    std::ostringstream msg ;
    msg << method << "( " << part.name()
        << " , " << top->name << " ) ERROR " ;
    if ( error_rank ) {
      msg << ": different entity_rank " << part.primary_entity_rank()
          << " != " << entity_rank ;
    }
    if ( error_top ) {
      msg << ": different topology " << i->second->name
          << " != " << top->name ;
    }
    throw msg.str();
  }

  if ( ! duplicate ) {
    m_part_top_map.insert( i , value );
  }

}

//----------------------------------------------------------------------------

const CellTopologyData *
TopologicalMetaData::internal_get_cell_topology( PartOrdinal part_ordinal ) const
{
  typedef std::pair< PartOrdinal , const CellTopologyData * > ValueType ;

  ValueType tmp( part_ordinal , NULL );

  std::vector< ValueType >::const_iterator
    i = std::lower_bound( m_part_top_map.begin() , m_part_top_map.end() , tmp );

  return i != m_part_top_map.end() && i->first == part_ordinal ? i->second : NULL ;
}

//----------------------------------------------------------------------------

void TopologicalMetaData::throw_ambiguous( const Part & part ) const
{
  static const char method[] = "stk::mesh::get_cell_topology" ;

  const PartVector & supersets = part.supersets();

  std::ostringstream msg ;

  msg << method << "( Part[" << part.name()
      << "] ) has ambiguous cell topology {" ;

  for ( PartVector::const_iterator
        j = supersets.begin(); j != supersets.end() ; ++j ) {

    const CellTopologyData * const top =
      internal_get_cell_topology( (*j)->mesh_meta_data_ordinal() );

    if ( top ) {
      msg << " ( Superset[" << (*j)->name()
          << "] -> " << top->name << " )" ;
    }
  }

  msg << " }" ;
  throw std::runtime_error( msg.str() );
}

void TopologicalMetaData::throw_ambiguous( const Bucket & bucket ) const
{
  static const char method[] = "stk::mesh::get_cell_topology" ;

  const PartVector & all_parts = m_meta_data.get_parts();

  std::ostringstream msg ;

  msg << method << "( Bucket ) has ambiguous cell topology {" ;

  const std::pair< const unsigned * , const unsigned * >
    supersets = bucket.superset_part_ordinals();

  for ( const unsigned * j = supersets.first ; j != supersets.second ; ++j ) {
    const CellTopologyData * const top =
      internal_get_cell_topology( *j );
    if ( top ) {
      msg << " ( Superset[" << all_parts[*j]->name()
          << "] -> " << top->name << " )" ;
    }
  }
  msg << " }" ;
  throw std::runtime_error( msg.str() );
}

void TopologicalMetaData::throw_required( const Part & part ,
                                          const char * required_by )
{
  static const char method[] = "stk::mesh::get_cell_topology" ;

  std::ostringstream msg ;

  msg << required_by << " Failed to obtain cell topology from "
      << method << "( Part[" << part.name() << "] );" ;

  throw std::runtime_error( msg.str() );
}

void TopologicalMetaData::throw_required( const Bucket & bucket ,
                                          const char * required_by )
{
  static const char method[] = "stk::mesh::get_cell_topology" ;

  std::ostringstream msg ;

  msg << required_by << " Failed to obtain cell topology from "
      << method << "( Bucket[" ;

  const BulkData   & bulk_data = bucket.mesh();
  const MetaData   & meta_data = bulk_data.mesh_meta_data();
  const PartVector & all_parts = meta_data.get_parts();

  const std::pair< const unsigned * , const unsigned * >
    supersets = bucket.superset_part_ordinals();

  for ( const unsigned * j = supersets.first ; j != supersets.second ; ++j ) {
    msg << " " << all_parts[*j]->name();
  }

  msg << " ] )" ;

  throw std::runtime_error( msg.str() );
}

//----------------------------------------------------------------------------

const CellTopologyData *
TopologicalMetaData::get_cell_topology( const Part & part ,
                                        const char * required_by )
{
  const TopologicalMetaData * const self = internal_get( part.mesh_meta_data() );

  const CellTopologyData * cell_top = NULL ;

  if ( self ) {
    cell_top = self->internal_get_cell_topology( part.mesh_meta_data_ordinal() );

    if ( ! cell_top ) {
      const PartVector & supersets = part.supersets();

      for ( PartVector::const_iterator
            j = supersets.begin(); j != supersets.end() ; ++j ) {

        const CellTopologyData * const top =
          self->internal_get_cell_topology( (*j)->mesh_meta_data_ordinal() );

        if ( ! cell_top ) {
          cell_top = top ;
        }
        else if ( top && top != cell_top ) {
          self->throw_ambiguous( part );
        }
      }
    }
  }

  if ( ! cell_top && required_by ) {
    throw_required( part , required_by );
  }

  return cell_top ;
}

//----------------------------------------------------------------------------

const CellTopologyData *
TopologicalMetaData::get_cell_topology( const Bucket & bucket ,
                                        const char * required_by )
{
  const BulkData   & bulk_data = bucket.mesh();
  const MetaData   & meta_data = bulk_data.mesh_meta_data();
  const PartVector & all_parts = meta_data.get_parts();

  const TopologicalMetaData * const self = internal_get( meta_data );

  const CellTopologyData * cell_top = NULL ;

  if ( self ) {
    const std::pair< const unsigned * , const unsigned * >
      supersets = bucket.superset_part_ordinals();

    for ( const unsigned * j = supersets.first ; j != supersets.second ; ++j ) {

      const Part & part = * all_parts[*j] ;

      if ( part.primary_entity_rank() == bucket.entity_rank() ) {

        const CellTopologyData * const top =
          self->internal_get_cell_topology( *j );

        if ( ! cell_top ) {
          cell_top = top ;
        }
        else if ( top && top != cell_top ) {
          self->throw_ambiguous( bucket );
        }
      }
    }
  }

  if ( ! cell_top && required_by ) {
    throw_required( bucket , required_by );
  }

  return cell_top ;
}

//----------------------------------------------------------------------------

const CellTopologyData *
TopologicalMetaData::get_cell_topology( const Entity & entity ,
                                        const char * required_by )
{
  return get_cell_topology( entity.bucket(), required_by );
}

}
}

