/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <sstream>
#include <stdexcept>

#include <Shards_BasicTopologies.hpp>

#include <stk_mesh/femImpl/PartCellTopologyMap.hpp>
#include <stk_mesh/femImpl/FiniteElementMeshImpl.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>

namespace stk {
namespace mesh {
namespace impl {

//----------------------------------------------------------------------------

PartCellTopologyMap::~PartCellTopologyMap() {}

PartCellTopologyMap::PartCellTopologyMap(
    MetaData & meta_data , unsigned spatial_dimension )
  : m_parts()
  , m_metaData( meta_data )
  , m_spatial_dimension( spatial_dimension )
{
  static const char method[] = "stk::mesh::impl::PartCellTopologyMap" ;

  verify_spatial_dimension( method , spatial_dimension );

  declare_part( * shards::getCellTopologyData< shards::Node >() , 0 );

  declare_part( * shards::getCellTopologyData< shards::Line<2> >() , 1 );
  declare_part( * shards::getCellTopologyData< shards::Line<3> >() , 1 );

  declare_part( * shards::getCellTopologyData< shards::Particle >() , spatial_dimension );

  if ( 1 < spatial_dimension ) {
    declare_part( * shards::getCellTopologyData< shards::Triangle<3> >() , 2 );
    declare_part( * shards::getCellTopologyData< shards::Triangle<6> >() , 2 );
    declare_part( * shards::getCellTopologyData< shards::Triangle<4> >() , 2 );

    declare_part( * shards::getCellTopologyData< shards::Quadrilateral<4> >() , 2 );
    declare_part( * shards::getCellTopologyData< shards::Quadrilateral<8> >() , 2 );
    declare_part( * shards::getCellTopologyData< shards::Quadrilateral<9> >() , 2 );

    declare_part( * shards::getCellTopologyData< shards::Beam<2> >() , spatial_dimension );
    declare_part( * shards::getCellTopologyData< shards::Beam<3> >() , spatial_dimension );
  }

  if ( 2 == spatial_dimension ) {
    declare_part( * shards::getCellTopologyData< shards::ShellLine<2> >() , 2 );
    declare_part( * shards::getCellTopologyData< shards::ShellLine<3> >() , 2 );
  }

  if ( 2 < spatial_dimension ) {
    declare_part( * shards::getCellTopologyData< shards::Tetrahedron<4> >() , 3 );
    declare_part( * shards::getCellTopologyData< shards::Tetrahedron<10> >() , 3 );
    declare_part( * shards::getCellTopologyData< shards::Tetrahedron<8> >() , 3 );

    declare_part( * shards::getCellTopologyData< shards::Pyramid<5> >() , 3 );
    declare_part( * shards::getCellTopologyData< shards::Pyramid<13> >() , 3 );
    declare_part( * shards::getCellTopologyData< shards::Pyramid<14> >() , 3 );

    declare_part( * shards::getCellTopologyData< shards::Wedge<6> >() , 3 );
    declare_part( * shards::getCellTopologyData< shards::Wedge<15> >() , 3 );
    declare_part( * shards::getCellTopologyData< shards::Wedge<18> >() , 3 );

    declare_part( * shards::getCellTopologyData< shards::Hexahedron<8> >() , 3 );
    declare_part( * shards::getCellTopologyData< shards::Hexahedron<20> >() , 3 );
    declare_part( * shards::getCellTopologyData< shards::Hexahedron<27> >() , 3 );

    declare_part( * shards::getCellTopologyData< shards::ShellTriangle<3> >() , 3 );
    declare_part( * shards::getCellTopologyData< shards::ShellTriangle<6> >() , 3 );

    declare_part( * shards::getCellTopologyData< shards::ShellQuadrilateral<4> >() , 3 );
    declare_part( * shards::getCellTopologyData< shards::ShellQuadrilateral<8> >() , 3 );
    declare_part( * shards::getCellTopologyData< shards::ShellQuadrilateral<9> >() , 3 );
  }
}

//----------------------------------------------------------------------------

inline
const CellTopologyData *
query_cell_topology( const Part & part )
{ return part.attribute<CellTopologyData>(); }

Part & PartCellTopologyMap::declare_part( const CellTopologyData & arg_top ,
                                          unsigned entity_rank )
{
  const bool bad_dimension = m_spatial_dimension < arg_top.dimension ;

  Part * part = m_metaData.get_part( arg_top.name , NULL );

  const CellTopologyData * const top =
    part ? query_cell_topology( *part ) : & arg_top ;

  const bool bad_existing = top != & arg_top ;
  const bool bad_rank = part ? entity_rank != part->primary_entity_rank() : false ;

  if ( bad_dimension || bad_existing || bad_rank ) {

    std::ostringstream msg ;
    msg << "stk::mesh::femImpl::PartCellTopologyMap::declare_part( "
        << arg_top.name << " , " << entity_rank << " ) ERROR" ;

    if ( bad_dimension ) {
      msg << ", spatial_dimension("
          << m_spatial_dimension << ") < topological_dimension("
          << arg_top.dimension << ")" ;
    }
    if ( bad_existing ) {
      if ( top ) { msg << ", has different topology " << top->name ; }
      else       { msg << ", has NO topology " ; }
    }
    if ( bad_rank ) {
      msg << ", has different entity_rank" ;
    }
    
    throw std::runtime_error( msg.str() );
  }
  else if ( ! part ) {
    part = & m_metaData.declare_part( top->name , entity_rank );

    m_metaData.declare_attribute_no_delete<CellTopologyData>(*part, &arg_top);

    m_parts.push_back( part );
  }

  return *part ;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

Part * PartCellTopologyMap::get_part( const CellTopologyData & top ,
                                      const char * const required_by ) const
{
  PartVector::const_iterator i = m_parts.begin();

  for ( ; i != m_parts.end() && query_cell_topology( **i ) != & top ; ++i );

  if ( required_by && i == m_parts.end() ) {
    std::ostringstream msg ;
    msg << required_by << " requires a Part for cell topology "
        << top.name << ", none was found" ;
    throw std::runtime_error( msg.str() );
  }

  return i != m_parts.end() ? *i : (Part *) NULL ;
}

//----------------------------------------------------------------------------

const CellTopologyData *
PartCellTopologyMap::get_cell_topology( const Part & part ,
                                        const char * const required_by )
{
  static const char method[] = "stk::mesh::femImpl::PartCellTopologyMap::get_cell_topology" ;

  const MetaData & meta_data = part.mesh_meta_data();

  const CellTopologyData * cell_top = NULL ;

  const bool has_entity_rank = part.primary_entity_rank() < meta_data.entity_rank_count();

  if ( has_entity_rank ) {

    cell_top = query_cell_topology( part );

    if ( 0 == cell_top ) {

      bool ok = true ;

      const PartVector & supersets = part.supersets();

      for ( PartVector::const_iterator
             i = supersets.begin() ; ok && i != supersets.end() ; ++i ) {

        const Part & super = **i ;

        if ( super.primary_entity_rank() == part.primary_entity_rank() ) {

          const CellTopologyData * const top = query_cell_topology( super );

          if ( 0 == cell_top ) { // Haven't found one
            cell_top = top ;
          }
          else { // Error to find another
            ok = 0 == top ;
          }
        }
      }

      if ( ! ok ) {

        std::ostringstream msg ;

        msg << method << "( " << part.name()
            << " ) ERROR, has multiple cell topologies: " ;

        for ( PartVector::const_iterator
               i = supersets.begin() ; i != supersets.end() ; ++i ) {

          const Part & super = **i ;

          if ( super.primary_entity_rank() == part.primary_entity_rank() &&
               0 != query_cell_topology( super ) ) {
            msg << " " << super.name();
          }
        }

        throw std::runtime_error( msg.str() );
      }
    }
  }

  if ( required_by && ! cell_top ) {
    std::ostringstream msg ;
    msg << required_by << " required Part["
        << part.name() << "] to have a cell topology, none was found" ;
    throw std::runtime_error( msg.str() );
  }

  return cell_top ;
}

const CellTopologyData *
PartCellTopologyMap::get_cell_topology( const Bucket & bucket ,
                                        const char * const required_by )
{
  static const char method[] = "stk::mesh::femImpl::PartCellTopologyMap::get_cell_topology" ;

  const BulkData   & bulk_data = bucket.mesh();
  const MetaData   & meta_data = bulk_data.mesh_meta_data();
  const PartVector & all_parts = meta_data.get_parts();
  const unsigned entity_rank = bucket.entity_rank();

  const CellTopologyData * cell_top = NULL ;

  bool ok = true ;

  // Iterate the bucket's superset parts looking for a part with
  // the same primary entity rank and a topology.

  std::pair<const unsigned * , const unsigned * >
    part_ordinals = bucket.superset_part_ordinals();

  for ( ; ok && part_ordinals.first != part_ordinals.second ; ++part_ordinals.first ) {
    const Part & part = *all_parts[ * part_ordinals.first ];

    if ( entity_rank == part.primary_entity_rank() ) {

      const CellTopologyData * const top = query_cell_topology( part );

      if ( 0 == cell_top ) { // Not yet found
        cell_top = top ;
      }
      else { // Has been found, cannot have another
        ok = 0 == top ;
      }
    }
  }

  if ( ! ok ) {
    part_ordinals = bucket.superset_part_ordinals();

    std::ostringstream msg ;
    msg << method << "( Bucket ) ERROR, has multiple cell topologies:" ;

    for ( ; part_ordinals.first != part_ordinals.second ; ++part_ordinals.first ) {
      const Part & part = *all_parts[ * part_ordinals.first ];
      if ( entity_rank == part.primary_entity_rank() &&
           0 != query_cell_topology( part ) ) {
        msg << " " << part.name() ;
      }
    }
    throw std::runtime_error( msg.str() );
  }

  if ( required_by && ! cell_top ) {
    std::ostringstream msg ;
    msg << required_by << " required a Bucket to have a cell topology, none was found" ;
    throw std::runtime_error( msg.str() );
  }

  return cell_top ;
}

//----------------------------------------------------------------------------

}
}
}

