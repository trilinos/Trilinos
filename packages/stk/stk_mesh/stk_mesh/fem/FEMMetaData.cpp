#include <set>
#include <stk_mesh/fem/FEMMetaData.hpp>

#include <Shards_CellTopology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Ghosting.hpp>

namespace stk {
namespace mesh {
namespace fem {

namespace {

void assign_cell_topology(
  FEMMetaData::PartCellTopologyVector & part_cell_topology_vector,
  size_t                   part_ordinal,
  const fem::CellTopology  cell_topology)
{
  if (part_ordinal >= part_cell_topology_vector.size())
    part_cell_topology_vector.resize(part_ordinal + 1);

  part_cell_topology_vector[part_ordinal] = cell_topology;

  if (!cell_topology.getCellTopologyData())
    {
      std::cout << "bad topology in FEMMetaData::assign_cell_topology" << std::endl;
    }

  ThrowRequireMsg(cell_topology.getCellTopologyData(), "bad topology in FEMMetaData::assign_cell_topology");
}

} // namespace

FEMMetaData::FEMMetaData()
:
    m_fem_initialized(false),
    m_spatial_dimension(0),
    m_side_rank(INVALID_RANK),
    m_element_rank(INVALID_RANK)
{
  // Attach FEMMetaData as attribute on MetaData to enable "get accessors" to FEMMetaData
  m_meta_data.declare_attribute_no_delete<FEMMetaData>(this);
}

FEMMetaData::FEMMetaData(size_t spatial_dimension,
                         const std::vector<std::string>& in_entity_rank_names)
  :
    m_fem_initialized(false),
    m_spatial_dimension(0),
    m_side_rank(INVALID_RANK),
    m_element_rank(INVALID_RANK)
{
  // Attach FEMMetaData as attribute on MetaData to enable "get accessors" to FEMMetaData
  m_meta_data.declare_attribute_no_delete<FEMMetaData>(this);

  FEM_initialize(spatial_dimension, in_entity_rank_names);
}

void FEMMetaData::FEM_initialize(size_t spatial_dimension, const std::vector<std::string>& rank_names)
{
  ThrowRequireMsg(!m_fem_initialized,"FEM functionality in FEMMetaData can only be initialized once.");
  if ( rank_names.empty() ) {
    m_entity_rank_names = fem::entity_rank_names(spatial_dimension);
  }
  else {
    ThrowRequireMsg(rank_names.size() >= spatial_dimension+1,
                    "Entity rank name vector must name every rank");
    m_entity_rank_names = rank_names;
  }
  internal_set_spatial_dimension_and_ranks(spatial_dimension);
  m_meta_data.set_entity_rank_names(m_entity_rank_names);
  m_fem_initialized = true;
  internal_declare_known_cell_topology_parts();
}

void FEMMetaData::internal_set_spatial_dimension_and_ranks(size_t spatial_dimension)
{
  ThrowRequireMsg( spatial_dimension != 0, "FEMMetaData::internal_set_spatial_dimension_and_ranks: spatial_dimension == 0!");
  m_spatial_dimension = spatial_dimension;

  // TODO:  Decide on correct terminology for the FEM Entity Ranks (consider topological vs spatial names and incompatibilities).
  // spatial_dimension = 1
  // node = 0, edge = 0, face = 0, side = 0, element = 1
  // spatial_dimension = 2
  // node = 0, edge = 1, face = 1, side = 1, element = 2
  // spatial_dimension = 3
  // node = 0, edge = 1, face = 2, side = 2, element = 3
  // spatial_dimension = 4
  // node = 0, edge = 1, face = 2, side = 3, element = 4
  m_side_rank = m_spatial_dimension - 1;
  m_element_rank = m_spatial_dimension;

}

void FEMMetaData::internal_declare_known_cell_topology_parts()
{
  // Load up appropriate standard cell topologies.
  register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Node >()), NODE_RANK);

  if (m_spatial_dimension == 1) {

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Particle >()), m_element_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Line<2> >()), m_element_rank); // ???
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Line<3> >()), m_element_rank); // ???

  }

  else if (m_spatial_dimension == 2) {

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Line<2> >()), m_side_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Line<3> >()), m_side_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Particle >()), m_element_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Triangle<3> >()), m_element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Triangle<6> >()), m_element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Triangle<4> >()), m_element_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()), m_element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<8> >()), m_element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<9> >()), m_element_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Beam<2> >()), m_element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Beam<3> >()), m_element_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::ShellLine<2> >()), m_element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::ShellLine<3> >()), m_element_rank);
  }

  else if (m_spatial_dimension == 3) {

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Line<2> >()), EDGE_RANK);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Line<3> >()), EDGE_RANK);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Triangle<3> >()), m_side_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Triangle<6> >()), m_side_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Triangle<4> >()), m_side_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()), m_side_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<8> >()), m_side_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<9> >()), m_side_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Particle >()), m_element_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Beam<2> >()), m_element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Beam<3> >()), m_element_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Tetrahedron<4> >()), m_element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Tetrahedron<10> >()), m_element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Tetrahedron<8> >()), m_element_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Pyramid<5> >()), m_element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Pyramid<13> >()), m_element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Pyramid<14> >()), m_element_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Wedge<6> >()), m_element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Wedge<15> >()), m_element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Wedge<18> >()), m_element_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Hexahedron<8> >()), m_element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Hexahedron<20> >()), m_element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::Hexahedron<27> >()), m_element_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::ShellTriangle<3> >()), m_element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::ShellTriangle<6> >()), m_element_rank);

    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::ShellQuadrilateral<4> >()), m_element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::ShellQuadrilateral<8> >()), m_element_rank);
    register_cell_topology(fem::CellTopology(shards::getCellTopologyData< shards::ShellQuadrilateral<9> >()), m_element_rank);
  }
}

void FEMMetaData::register_cell_topology(const fem::CellTopology cell_topology, EntityRank entity_rank)
{
  ThrowRequireMsg(is_FEM_initialized(),"FEMMetaData::register_cell_topology: FEM_initialize() must be called before this function");

  CellTopologyPartEntityRankMap::const_iterator it = m_cellTopologyPartEntityRankMap.find(cell_topology);

  const bool       duplicate     = it != m_cellTopologyPartEntityRankMap.end();
  const EntityRank existing_rank = duplicate ? (*it).second.second : 0;

  ThrowInvalidArgMsgIf(m_spatial_dimension < entity_rank,
    "entity_rank " << entity_rank << ", " <<
    "exceeds maximum spatial_dimension = " << m_spatial_dimension );

  ThrowErrorMsgIf(duplicate && existing_rank != entity_rank,
    "For args: cell_topolgy " << cell_topology.getName() << " and entity_rank " << entity_rank << ", " <<
    "previously declared rank = " << existing_rank );

  if (! duplicate) {
    std::string part_name = std::string("{FEM_ROOT_CELL_TOPOLOGY_PART_") + std::string(cell_topology.getName()) + std::string("}");

    ThrowErrorMsgIf(get_part(part_name) != 0, "Cannot register topology with same name as existing part '" << cell_topology.getName() << "'" );

    Part &part = declare_part(part_name, entity_rank);
    m_cellTopologyPartEntityRankMap[cell_topology] = CellTopologyPartEntityRankMap::mapped_type(&part, entity_rank);

    assign_cell_topology(m_partCellTopologyVector, part.mesh_meta_data_ordinal(), cell_topology);
  }
  //check_topo_db();
}


fem::CellTopology
FEMMetaData::get_cell_topology(
  const std::string &   topology_name) const
{
  std::string part_name = std::string("{FEM_ROOT_CELL_TOPOLOGY_PART_") + topology_name + std::string("}");

  Part *part = get_part(part_name);
  if (part)
    return get_cell_topology(*part);
  else
    return fem::CellTopology();
}


Part &FEMMetaData::get_cell_topology_root_part(const fem::CellTopology cell_topology) const
{
  ThrowRequireMsg(is_FEM_initialized(),"FEMMetaData::get_cell_topology_root_part: FEM_initialize() must be called before this function");
  CellTopologyPartEntityRankMap::const_iterator it = m_cellTopologyPartEntityRankMap.find(cell_topology);
  ThrowErrorMsgIf(it == m_cellTopologyPartEntityRankMap.end(),
                  "Cell topology " << cell_topology.getName() <<
                  " has not been registered");

  return *(*it).second.first;
}

/// Note:  This function only uses the PartCellTopologyVector to look up the
/// cell topology for a given part.
/// This depends on declare_part_subset to update this vector correctly.  If a
/// cell topology is not defined for the given part, then an invalid Cell
/// Topology object will be returned.
fem::CellTopology FEMMetaData::get_cell_topology( const Part & part) const
{
  ThrowRequireMsg(is_FEM_initialized(),"FEMMetaData::get_cell_topology: FEM_initialize() must be called before this function");
  fem::CellTopology cell_topology;

  PartOrdinal part_ordinal = part.mesh_meta_data_ordinal();
  if (part_ordinal < m_partCellTopologyVector.size())
    {
      cell_topology = m_partCellTopologyVector[part_ordinal];
    }

  return cell_topology;
}

#if 0
  void FEMMetaData::check_topo_db()
  {
    std::cout << "FEMMetaData::check_topo_db... m_partCellTopologyVector.size() = " << m_partCellTopologyVector.size() <<  std::endl;

  fem::CellTopology cell_topology;

  for (unsigned i = 0; i <  m_partCellTopologyVector.size(); i++)
    {
      cell_topology = m_partCellTopologyVector[i];
  if (!cell_topology.getCellTopologyData())
    {
      std::cout << "bad topology in FEMMetaData::check_topo_db" << std::endl;
    }
      ThrowRequireMsg(cell_topology.getCellTopologyData(), "bad topology in FEMMetaData::check_topo_db");

    }
    std::cout << "FEMMetaData::check_topo_db...done" << std::endl;

  }
#endif

namespace {

bool root_part_in_subset(stk::mesh::Part & part)
{
  if (is_cell_topology_root_part(part)) {
    return true;
  }
  const PartVector & subsets = part.subsets();
  for (PartVector::const_iterator it=subsets.begin() ; it != subsets.end() ; ++it) {
    if (is_cell_topology_root_part( **it )) {
      return true;
    }
  }
  return false;
}

void find_cell_topologies_in_part_and_subsets_of_same_rank(const Part & part, EntityRank rank, std::set<fem::CellTopology> & topologies_found)
{
  fem::FEMMetaData & fem_meta = fem::FEMMetaData::get(part);
  fem::CellTopology top = fem_meta.get_cell_topology(part);
  if ((top.isValid() && (part.primary_entity_rank() == rank))) {
    topologies_found.insert(top);
  }
  const PartVector & subsets = part.subsets();
  for (PartVector::const_iterator it=subsets.begin() ; it != subsets.end() ; ++it) {
    top = fem_meta.get_cell_topology(**it);
    if (top.isValid() && ( (**it).primary_entity_rank() == rank) ) {
      topologies_found.insert(top);
    }
  }
}

} // namespace

void FEMMetaData::declare_part_subset( Part & superset , Part & subset )
{
  ThrowRequireMsg(is_FEM_initialized(),"FEMMetaData::declare_part_subset: FEM_initialize() must be called before this function");
  fem::CellTopology superset_top = get_cell_topology(superset);

  const bool no_superset_topology = !superset_top.isValid();
  if ( no_superset_topology ) {
    m_meta_data.declare_part_subset(superset,subset);
    return;
  }
  // Check for cell topology root parts in subset or subset's subsets
  const bool subset_has_root_part = root_part_in_subset(subset);
  ThrowErrorMsgIf( subset_has_root_part, "FEMMetaData::declare_part_subset:  Error, root cell topology part found in subset or below." );

  std::set<fem::CellTopology> cell_topologies;
  find_cell_topologies_in_part_and_subsets_of_same_rank(subset,superset.primary_entity_rank(),cell_topologies);

  ThrowErrorMsgIf( cell_topologies.size() > 1,
      "FEMMetaData::declare_part_subset:  Error, multiple cell topologies of rank "
      << superset.primary_entity_rank()
      << " defined below subset"
      );
  const bool non_matching_cell_topology = ((cell_topologies.size() == 1) && (*cell_topologies.begin() != superset_top));
  ThrowErrorMsgIf( non_matching_cell_topology,
      "FEMMetaData::declare_part_subset:  Error, superset topology = "
      << superset_top.getName() << " does not match the topology = "
      << cell_topologies.begin()->getName()
      << " coming from the subset part"
      );
  // Everything is Okay!
  m_meta_data.declare_part_subset(superset,subset);
  // Update PartCellTopologyVector for "subset" and same-rank subsets, ad nauseum
  if (subset.primary_entity_rank() == superset.primary_entity_rank()) {
    assign_cell_topology(m_partCellTopologyVector, subset.mesh_meta_data_ordinal(), superset_top);
    const PartVector & subset_parts = subset.subsets();
    for (PartVector::const_iterator it=subset_parts.begin() ; it != subset_parts.end() ; ++it) {
      const Part & it_part = **it;
      if (it_part.primary_entity_rank() == superset.primary_entity_rank()) {
        assign_cell_topology(m_partCellTopologyVector, it_part.mesh_meta_data_ordinal(), superset_top);
      }
    }
  }
}

EntityRank FEMMetaData::get_entity_rank(
  const fem::CellTopology       cell_topology) const
{
  CellTopologyPartEntityRankMap::const_iterator it = m_cellTopologyPartEntityRankMap.find(cell_topology);
  if (it == m_cellTopologyPartEntityRankMap.end())
    return INVALID_RANK;
  else
    return (*it).second.second;
}

//--------------------------------------------------------------------------------
// Free Functions
//--------------------------------------------------------------------------------

bool is_cell_topology_root_part(const Part & part) {
  fem::FEMMetaData & fem_meta = fem::FEMMetaData::get(part);
  fem::CellTopology top = fem_meta.get_cell_topology(part);
  if (top.isValid()) {
    const Part & root_part = fem_meta.get_cell_topology_root_part(top);
    return (root_part == part);
  }
  return false;
}

/// This is a convenience function to get the root cell topology part and then
/// call declare_part_subset.
/// Note:  FEMMetaData::declare_part_subset is the function that actually
/// updates the PartCellTopologyVector in FEMMetaData for fast look-up of the
/// Cell Topology.
void set_cell_topology(
  Part &                        part,
  fem::CellTopology             cell_topology)
{
  FEMMetaData& fem_meta = FEMMetaData::get(part);

  ThrowRequireMsg(fem_meta.is_FEM_initialized(),"set_cell_topology: FEM_initialize() must be called before this function");

  Part &root_part = fem_meta.get_cell_topology_root_part(cell_topology);
  fem_meta.declare_part_subset(root_part, part);
}

std::vector<std::string>
entity_rank_names( size_t spatial_dimension )
{
  ThrowInvalidArgMsgIf( spatial_dimension < 1 || 3 < spatial_dimension,
                        "Invalid spatial dimension = " << spatial_dimension );

  std::vector< std::string > names ;

  names.reserve( spatial_dimension + 1 );

  names.push_back( std::string( "NODE" ) );

  if ( 1 < spatial_dimension ) { names.push_back( std::string("EDGE") ); }
  if ( 2 < spatial_dimension ) { names.push_back( std::string("FACE") ); }

  names.push_back( std::string("ELEMENT") );

  return names ;
}


CellTopology
get_cell_topology(
  const Bucket &                bucket)
{
  const BulkData   & bulk_data = BulkData::get(bucket);
  const MetaData   & meta_data = MetaData::get(bulk_data);
  const PartVector & all_parts = meta_data.get_parts();

  FEMMetaData &fem = FEMMetaData::get(meta_data);

  CellTopology cell_topology;

  const std::pair< const unsigned *, const unsigned * > supersets = bucket.superset_part_ordinals();

  if (supersets.first != supersets.second) {
    const Part *first_found_part = 0;

    for ( const unsigned * it = supersets.first ; it != supersets.second ; ++it ) {

      const Part & part = * all_parts[*it] ;

      if ( part.primary_entity_rank() == bucket.entity_rank() ) {

        CellTopology top = fem.get_cell_topology( part );

        if ( ! cell_topology.getCellTopologyData() ) {
          cell_topology = top ;

          if (!first_found_part)
            first_found_part = &part;
        }
        else {
          ThrowErrorMsgIf( top.getCellTopologyData() && top != cell_topology,
            "Cell topology is ambiguously defined. It is defined as " << cell_topology.getName() <<
            " on part " << first_found_part->name() << " and as " << top.getName() << " on its superset part " << part.name() );
        }
      }
    }
  }

  return cell_topology ;
}

} // namespace fem
} // namespace mesh
} // namespace stk
