#include <algorithm>

#include <Shards_CellTopology.hpp>

#include <stk_mesh/fem/DefaultFEM.hpp>

#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
#include <stk_mesh/fem/TopologyHelpersDeprecated.hpp>
#endif

namespace stk {
namespace mesh {

namespace {

template<class T>
void
insert(
  std::vector<T> &      vector,
  size_t                index,
  const T &             value)
{
  if (index >= vector.size())
    vector.resize(index + 1);

  vector[index] = value;
}

} // namespace <unnamed>

// namespace fem {


DefaultFEM::DefaultFEM(
  MetaData &            meta_data,
  size_t                spatial_dimension)
  : FEMInterface(),
    m_metaData(meta_data),
    m_spatialDimension(spatial_dimension),
    m_cellTopologyPartEntityRankMap(),
    m_partCellTopologyVector()
{
  if (m_metaData.entity_rank_names().empty()) {
    m_metaData.set_entity_rank_names(fem::entity_rank_names(spatial_dimension));
  }

  set_fem_interface(meta_data, this);

  initialize(spatial_dimension);
}


DefaultFEM::DefaultFEM(
  MetaData &            meta_data)
  : FEMInterface(),
    m_metaData(meta_data),
    m_spatialDimension(fem::INVALID_RANK),
    m_cellTopologyPartEntityRankMap(),
    m_partCellTopologyVector()
{
  set_fem_interface(m_metaData, this);
}


void
DefaultFEM::set_spatial_dimension(
  size_t                spatial_dimension)
{
  if (m_spatialDimension != fem::INVALID_RANK && m_spatialDimension != spatial_dimension) {
    std::ostringstream oss;

    oss << "The spatial dimension has already been defined as " << m_spatialDimension << " and cannot be changed to " << spatial_dimension;
    
    throw std::runtime_error(oss.str());
  }
  else if (m_spatialDimension == fem::INVALID_RANK) { 
    m_spatialDimension = spatial_dimension;

    initialize(spatial_dimension);
  }
}


void
DefaultFEM::initialize(
  size_t                spatial_dimension)
{
  // Load up appropriate standard cell topologies.
  register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Node >()), fem::NODE_RANK);

  if (spatial_dimension == 1) {
    const EntityRank element_rank = fem::element_rank(spatial_dimension);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Particle >()), element_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Line<2> >()), element_rank); // ???
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Line<3> >()), element_rank); // ???
  }

  else if (spatial_dimension == 2) {
    const EntityRank side_rank = fem::side_rank(spatial_dimension);
    const EntityRank element_rank = fem::element_rank(spatial_dimension);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Line<2> >()), side_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Line<3> >()), side_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Particle >()), element_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Triangle<3> >()), element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Triangle<6> >()), element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Triangle<4> >()), element_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()), element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<8> >()), element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<9> >()), element_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Beam<2> >()), element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Beam<3> >()), element_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::ShellLine<2> >()), element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::ShellLine<3> >()), element_rank);
  }

  else if (spatial_dimension == 3) {
    const EntityRank edge_rank = fem::edge_rank(spatial_dimension);
    const EntityRank side_rank = fem::side_rank(spatial_dimension);
    const EntityRank element_rank = fem::element_rank(spatial_dimension);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Line<2> >()), edge_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Line<3> >()), edge_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Triangle<3> >()), side_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Triangle<6> >()), side_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Triangle<4> >()), side_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()), side_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<8> >()), side_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<9> >()), side_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Particle >()), element_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Beam<2> >()), element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Beam<3> >()), element_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Tetrahedron<4> >()), element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Tetrahedron<10> >()), element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Tetrahedron<8> >()), element_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Pyramid<5> >()), element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Pyramid<13> >()), element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Pyramid<14> >()), element_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Wedge<6> >()), element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Wedge<15> >()), element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Wedge<18> >()), element_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Hexahedron<8> >()), element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Hexahedron<20> >()), element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Hexahedron<27> >()), element_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::ShellTriangle<3> >()), element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::ShellTriangle<6> >()), element_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::ShellQuadrilateral<4> >()), element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::ShellQuadrilateral<8> >()), element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::ShellQuadrilateral<9> >()), element_rank);
  }
}


EntityRank
DefaultFEM::get_entity_rank(
  const fem::CellTopology       cell_topology) const
{
  CellTopologyPartEntityRankMap::const_iterator it = m_cellTopologyPartEntityRankMap.find(cell_topology);
  if (it == m_cellTopologyPartEntityRankMap.end())
    return fem::INVALID_RANK;
  else
    return (*it).second.second;
}


Part &
DefaultFEM::get_part(
  const fem::CellTopology       cell_topology) const
{
  CellTopologyPartEntityRankMap::const_iterator it = m_cellTopologyPartEntityRankMap.find(cell_topology);
  ThrowErrorMsgIf(it == m_cellTopologyPartEntityRankMap.end(),
                  "Cell topology " << cell_topology.getName() <<
                  " has not been registered");

  return *(*it).second.first;
}


void
DefaultFEM::set_cell_topology(
  Part &                        part,
  fem::CellTopology             cell_topology)
{
  EntityRank entity_rank = get_entity_rank(cell_topology);
  if (entity_rank == fem::INVALID_RANK)
    entity_rank = cell_topology.getDimension();

  PartOrdinal part_ordinal = part.mesh_meta_data_ordinal();
  fem::CellTopology existing_cell_topology = part_ordinal < m_partCellTopologyVector.size() ? m_partCellTopologyVector[part_ordinal] : fem::CellTopology();

  const bool duplicate  = existing_cell_topology.getCellTopologyData() != 0;

  ThrowErrorMsgIf( part.primary_entity_rank() != entity_rank,
    "For args: part " << part.name() << " and topology " << cell_topology.getName() << ", " <<
    "different entity_rank " << part.primary_entity_rank() << " != " << entity_rank );

  ThrowErrorMsgIf( duplicate && cell_topology != existing_cell_topology,
    "For args: part " << part.name() << " and topology " << cell_topology.getName() << ", " <<
    "different topology " << cell_topology.getName() << " != " << existing_cell_topology.getName() );

  if (!duplicate) {
    insert(m_partCellTopologyVector, part.mesh_meta_data_ordinal(), cell_topology);

    Part &root_part = get_part(cell_topology);
    m_metaData.declare_part_subset(root_part, part);
  }

#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  set_cell_topology_deprecated(const_cast<Part &>(part), cell_topology.getCellTopologyData());
#endif
}


void
DefaultFEM::register_cell_topology(
  const fem::CellTopology       cell_topology,
  EntityRank                    entity_rank)
{
  CellTopologyPartEntityRankMap::const_iterator it = m_cellTopologyPartEntityRankMap.find(cell_topology);

  const bool       duplicate     = it != m_cellTopologyPartEntityRankMap.end();
  const EntityRank existing_rank = duplicate ? (*it).second.second : 0;

  ThrowInvalidArgMsgIf(m_spatialDimension < entity_rank,
    "entity_rank " << entity_rank << ", " <<
    "exceeds maximum spatial_dimension = " << m_spatialDimension );

  ThrowErrorMsgIf(duplicate && existing_rank != entity_rank,
    "For args: cell_topolgy " << cell_topology.getName() << " and entity_rank " << entity_rank << ", " <<
    "previously declared rank = " << existing_rank );
  
  if (! duplicate) {
    Part &part = m_metaData.declare_part(cell_topology.getName(), entity_rank);
    m_cellTopologyPartEntityRankMap[cell_topology] = CellTopologyPartEntityRankMap::mapped_type(&part, entity_rank);

    insert(m_partCellTopologyVector, part.mesh_meta_data_ordinal(), cell_topology);
  }
}


fem::CellTopology
DefaultFEM::get_cell_topology(
  const Part &          part) const
{
  fem::CellTopology cell_topology;

  PartOrdinal part_ordinal = part.mesh_meta_data_ordinal();
  if (part_ordinal < m_partCellTopologyVector.size())
    cell_topology = m_partCellTopologyVector[part_ordinal];

  return cell_topology;
}


// } // namespace fem
} // namespace mesh
} // namespace stk

