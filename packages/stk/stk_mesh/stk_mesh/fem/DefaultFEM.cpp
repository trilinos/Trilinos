
#include <algorithm>

#include <Shards_CellTopology.hpp>

#include <stk_mesh/fem/DefaultFEM.hpp>

#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
#include <stk_mesh/fem/TopologyHelpersDeprecated.hpp>
#endif

namespace stk {
namespace mesh {
// namespace fem {


DefaultFEM::DefaultFEM(
  MetaData &            meta_data, 
  size_t                spatial_dimension)
  : FEMInterface(),
    m_spatialDimension(spatial_dimension),
    m_topEntityRank(),
    m_partCellTopologyMap()
{
  set_fem_interface(meta_data, this);

  initialize(spatial_dimension);
}


DefaultFEM::DefaultFEM(
  MetaData &            meta_data)
  : FEMInterface(),
    m_spatialDimension(fem::INVALID_RANK),
    m_topEntityRank(),
    m_partCellTopologyMap()
{
  set_fem_interface(meta_data, this);
}


void
DefaultFEM::set_spatial_dimension(
  size_t                spatial_dimension)
{
  m_spatialDimension = spatial_dimension;

  initialize(spatial_dimension);
}

  
void
DefaultFEM::initialize(
  size_t                spatial_dimension)
{
  // Load up appropriate standard cell topologies.
  set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Node >()), fem::NODE_RANK );

  if (spatial_dimension == 1) {
    const EntityRank element_rank = fem::element_rank(spatial_dimension);

    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Particle >()), element_rank );

    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Line<2> >()), element_rank ); // ???
    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Line<3> >()), element_rank ); // ???
  }

  else if (spatial_dimension == 2) {
    const EntityRank side_rank = fem::side_rank(spatial_dimension);
    const EntityRank element_rank = fem::element_rank(spatial_dimension);

    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Line<2> >()), side_rank );
    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Line<3> >()), side_rank );

    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Particle >()), element_rank );

    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Triangle<3> >()), element_rank);
    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Triangle<6> >()), element_rank);
    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Triangle<4> >()), element_rank);

    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()), element_rank);
    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<8> >()), element_rank);
    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<9> >()), element_rank);

    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Beam<2> >()), element_rank );
    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Beam<3> >()), element_rank );

    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::ShellLine<2> >()), element_rank);
    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::ShellLine<3> >()), element_rank);
  }

  else if (spatial_dimension == 3) {
    const EntityRank edge_rank = fem::edge_rank(spatial_dimension);
    const EntityRank side_rank = fem::side_rank(spatial_dimension);
    const EntityRank element_rank = fem::element_rank(spatial_dimension);

    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Line<2> >()), edge_rank );
    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Line<3> >()), edge_rank );

    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Triangle<3> >()), side_rank);
    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Triangle<6> >()), side_rank);
    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Triangle<4> >()), side_rank);

    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()), side_rank);
    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<8> >()), side_rank);
    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<9> >()), side_rank);

    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Particle >()), element_rank );

    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Beam<2> >()), element_rank );
    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Beam<3> >()), element_rank );

    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Tetrahedron<4> >()), element_rank);
    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Tetrahedron<10> >()), element_rank);
    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Tetrahedron<8> >()), element_rank);

    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Pyramid<5> >()), element_rank);
    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Pyramid<13> >()), element_rank);
    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Pyramid<14> >()), element_rank);

    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Wedge<6> >()), element_rank);
    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Wedge<15> >()), element_rank);
    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Wedge<18> >()), element_rank);

    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Hexahedron<8> >()), element_rank);
    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Hexahedron<20> >()), element_rank);
    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::Hexahedron<27> >()), element_rank);

    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::ShellTriangle<3> >()), element_rank);
    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::ShellTriangle<6> >()), element_rank);

    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::ShellQuadrilateral<4> >()), element_rank);
    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::ShellQuadrilateral<8> >()), element_rank);
    set_entity_rank( fem::CellTopology(shards::getCellTopologyData< shards::ShellQuadrilateral<9> >()), element_rank);
  }
}


EntityRank
DefaultFEM::get_entity_rank(
  fem::CellTopology             top) const
{
  typedef std::pair< fem::CellTopology, EntityRank > ValueType ;

  std::vector< ValueType >::const_iterator i ;

  for ( i = m_topEntityRank.begin() ; i != m_topEntityRank.end() && top != i->first ; ++i )
    ;

  if (i == m_topEntityRank.end()) {
    return fem::INVALID_RANK;
  }

  return i->second ;
}


void
DefaultFEM::set_cell_topology( 
  const Part &                  part, 
  fem::CellTopology             top) 
{
  static const char method[] = "stk::mesh::DefaultFEM::set_cell_topology" ;

  EntityRank entity_rank = get_entity_rank( top );
  if (entity_rank == fem::INVALID_RANK)
    entity_rank = top.getDimension();

  typedef std::pair< PartOrdinal, fem::CellTopology > ValueType ;
  ValueType value( part.mesh_meta_data_ordinal(), top );

  std::vector< ValueType >::iterator
    i = std::lower_bound( m_partCellTopologyMap.begin(), m_partCellTopologyMap.end(), value );

  const bool duplicate  = i != m_partCellTopologyMap.end() && i->first == value.first ;
  const bool error_rank = part.primary_entity_rank() != entity_rank ;
  const bool error_top  = duplicate && i->second != value.second ;

  if ( error_rank || error_top ) {
    std::ostringstream oss ;
    oss << method << "( " << part.name()
        << ", " << top.getName() << " ) ERROR " ;
    if ( error_rank ) {
      oss << ": different entity_rank " << part.primary_entity_rank()
          << " != " << entity_rank ;
    }
    if ( error_top ) {
      oss << ": different topology " << i->second.getName()
          << " != " << top.getName() ;
    }
    throw std::runtime_error(oss.str());
  }

  if ( ! duplicate ) {
    m_partCellTopologyMap.insert( i, value );
  }

#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  set_cell_topology_deprecated(const_cast<Part&>(part), top.getTopologyData());
#endif
}


void
DefaultFEM::set_entity_rank(
  const fem::CellTopology       top, 
  EntityRank                    rank)
{
  static const char method[] = "stk::mesh::DefaultFEM::set_entity_rank" ;

  typedef std::pair< fem::CellTopology, EntityRank > ValueType ;

  std::vector< ValueType >::const_iterator i = m_topEntityRank.begin() ;
  for ( ; i != m_topEntityRank.end() && top != i->first ; ++i )
    ;

  const bool       duplicate     = i != m_topEntityRank.end();
  const EntityRank existing_rank = duplicate ? i->second : 0 ;

  const bool error_change = duplicate && existing_rank != rank ;
  const bool error_rank   = m_spatialDimension < rank ;

  if ( error_rank || error_change ) {
    std::ostringstream oss ;
    oss << method << "( " << top.getName()
        << ", rank = " << rank << " ) ERROR " ;
    if ( error_rank ) {
      oss << ": rank exceeds maximum spatial_dimension = "
          << m_spatialDimension ;
    }
    if ( error_change ) {
      oss << ": previously declared rank = " << existing_rank ;
    }
    throw std::runtime_error( oss.str() );
  }
  
  if ( ! duplicate ) {
    typedef std::pair< const fem::CellTopology, EntityRank > ValueType ;

    m_topEntityRank.push_back( ValueType( top, rank ) );
  }
}


fem::CellTopology
DefaultFEM::get_cell_topology(
  const Part &          part) const
{
  PartOrdinal part_ordinal = part.mesh_meta_data_ordinal();
  
  typedef std::pair< PartOrdinal, fem::CellTopology > ValueType ;

  ValueType tmp( part_ordinal, fem::CellTopology());

  std::vector< ValueType >::const_iterator
    i = std::lower_bound( m_partCellTopologyMap.begin(), m_partCellTopologyMap.end(), tmp );

  return i != m_partCellTopologyMap.end() && i->first == part_ordinal ? i->second : NULL ;
}

// } // namespace fem
} // namespace mesh
} // namespace stk

