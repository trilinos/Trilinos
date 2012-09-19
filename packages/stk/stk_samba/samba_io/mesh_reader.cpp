#include <Ioss_SubSystem.h>
#include <init/Ionit_Initializer.h>

#include <samba_io/mesh_reader.hpp>
#include <samba_io/ioss_topology_map.hpp>
#include <samba_io/sidesets_manager.hpp>

#include <algorithm>


#define MESH_READER_BABBLES

namespace samba {
namespace io {


bool mesh_reader::supports_sidesets(const mesh &mesh_arg)
{
  if (mesh_arg.spatial_dimension() == 2)
  {
    // Suport later...
    return false;
  }
  if (mesh_arg.spatial_dimension() == 3)
  {
    connectivity_map const &rel_map = mesh_arg.connectivity_map();
    if (!is_valid(rel_map(entity_rank::face(),entity_rank::node())))
    {
      std::cerr << "Mesh does not support sidesets!" << std::endl;
      return false;
    }
    return true;
  }
  return false;
}


bool mesh_reader::supports_direct_get_elements(const mesh &mesh_arg)
{
  if (mesh_arg.spatial_dimension() != 3)
  {
    return false;
  }
  connectivity_map const &rel_map = mesh_arg.connectivity_map();
  return is_valid(rel_map(entity_rank::face(),entity_rank::element()));
}


void mesh_reader::process_mesh()
{
  // For now.  Need to do something different if we decide to use this to pull data from
  // different Ioss::Regions into the mesh.
  BOOST_ASSERT(m_meshRankProcessed == -1);
  m_mesh.begin_modification();

  Ioss::NodeBlockContainer node_blocks = m_io_region->get_node_blocks();
  process_blocks(node_blocks, entity_rank::node());

  Ioss::ElementBlockContainer element_blocks = m_io_region->get_element_blocks();
  process_blocks(element_blocks, entity_rank::element());

  process_nodesets();

  m_mesh.end_modification();

  read_nodal_field(m_ioss_mapper_ptr->m_node_coordinates, "mesh_model_coordinates");
}


void mesh_reader::process_nodesets()
{
  const Ioss::NodeSetContainer& node_sets = m_io_region->get_nodesets();

  if (!node_sets.empty() && (m_ioss_mapper_rawptr->m_node_blocks.size() != 1))
  {
    throw std::runtime_error("samba::io::mesh_reader encounters node sets but does not have exactly one node block");
  }

  for (Ioss::NodeSetContainer::const_iterator nsets_iter=node_sets.begin(), nsets_end=node_sets.end();
       nsets_iter != nsets_end;
       ++nsets_iter)
  {
    process_nodeset(*nsets_iter);
  }
}


void mesh_reader::process_nodeset(const Ioss::NodeSet *nset)
{
  std::string nodeset_name = nset->name();
  std::vector<int> node_ids;

  nset->get_field_data("ids_raw", node_ids);

  size_t num_nodes = node_ids.size();

#ifdef MESH_READER_BABBLES
  std::cout << " Found a nodeset named \"" << nodeset_name << "\" with "
            << num_nodes << " members." << std::endl;
#endif

  for (size_t i = 0; i < num_nodes; ++i)
  {
    node_ids[i] -= 1;  // Change from one-based to zero-based.
  }

  entity_key_vector nodeset_node_keys(num_nodes);
  for (size_t i = 0; i < num_nodes; ++i)
  {
    nodeset_node_keys[i] = entity_key::create(entity_topology::node(),
                                              process_id::invalid(),
                                              entity_local_id::create(node_ids[i]));
  }

  entity_block_key nodeset_block = m_mesh.add_entity_block(nodeset_name, entity_rank::node());
  entity_block_key none = entity_block_key::invalid();
  m_mesh.move_entities(nodeset_node_keys.begin(), nodeset_node_keys.end(),
                       &nodeset_block, &nodeset_block + 1,
                       &none, &none );
}


std::vector<SidesetPtr> mesh_reader::process_sidesets()
{
  sidesets_manager ssm(m_mesh);

  m_ioss_mapper_rawptr->m_sidesets_info = ssm.process_region(m_io_region, *this);
  return m_ioss_mapper_rawptr->m_sidesets_info;
}


void mesh_reader::process_block(Ioss::EntityBlock *block, const entity_rank rank)
{
  BOOST_ASSERT(m_mesh.find_entity_block(block->name()) == entity_block_key::invalid());

  // std::cout << " Processing block named \"" << block->name() << "\" with rank "
  //           << rank << std::endl;

  // TO DO: ADAPT THIS TO HANDLE block WHOSE FIELD DATA VECTOR ELEMENTS ARE int64_t
  std::vector<int> global_ids;
  std::vector<int> connectivity;
  block->get_field_data("ids", global_ids);
  if (global_ids.empty())
  {
    return;
  }

  const Ioss::ElementTopology *topology = block->topology();
  const entity_topology block_topo = map_topology_ioss_to_cell(topology);

  const bool induce_membership = true;
  entity_block_key block_ebk = m_mesh.add_entity_block(block->name(), rank, induce_membership);
  BOOST_ASSERT(block_ebk != entity_block_key::invalid());
  m_entityBlocks[block_ebk] = block;

  const size_t block_entity_count = block->get_property("entity_count").get_int();
  const entity_key_interval block_entities =
    m_mesh.add_entities(block_topo, block_entity_count, &block_ebk, &block_ebk + 1);

  entity_key_vector &ordered_entity_keys = m_ordered_entity_keys[rank];
  size_t first_local = ordered_entity_keys.size();
  ordered_entity_keys.resize(first_local + block_entity_count);

  entity_key_to_ioss_id_field &to_idx_in_block = m_ioss_mapper_rawptr->m_to_ioss_idx_in_block;
  entity_key_to_ioss_id_field &to_local = m_ioss_mapper_rawptr->m_to_ioss_local_id;
  entity_key_to_ioss_id_field &to_global = m_ioss_mapper_rawptr->m_to_ioss_global_id;
  for (size_t i = 0; i < block_entity_count; ++i)
  {
    entity_key curr_entity = block_entities[i];
    to_idx_in_block[curr_entity] = i;

    int64_t zero_based_local_id = first_local + i;
    to_local[curr_entity] = zero_based_local_id;
    ordered_entity_keys[zero_based_local_id] = curr_entity;

    int64_t zero_based_global_id = global_ids[i] - 1;
    to_global[curr_entity] = zero_based_global_id;
  }

  // Need this at write time.
  m_ioss_mapper_rawptr->m_new_idx_in_block[block_ebk] = block_entity_count;

  BOOST_ASSERT(topology->number_nodes() >= 0);
  const size_t nodes_per_entity = topology->number_nodes();

  if (rank == entity_rank::node())
  {
    m_ioss_mapper_rawptr->m_node_blocks.push_back(block_ebk);
    return;
  }

  block->get_field_data("connectivity_raw", connectivity);

  for (size_t i = 0; i < block_entity_count; ++i)
  {
    const entity_key entity_from = block_entities[i];
    // Setting up connectivity connectivity from entity to nodes
    for (size_t j = 0; j < nodes_per_entity; j++)
    {
      // Find the node_id (in the mesh file) of the jth node that the ith entity in the block is
      // connected to.

      // Exodus begins array indices at 1, so we need to subtract 1 to get 0-based indexing.
      const int node_id = connectivity[i*nodes_per_entity + j] - 1;
      const entity_key node_to = entity_key::create(entity_topology::node(),
                                                    process_id::invalid(),
                                                    entity_local_id::create(node_id));
      // std::cout << "Adding connectivity between " << entity_from << " and " << node_to << std::endl;
      const connectivity_ordinal ro_j = {j};
      m_mesh.add_connectivity(entity_from, node_to, ro_j);
      m_mesh.add_connectivity(node_to, entity_from, ro_j);
    }
  }
}


}  // namespace io
}  // namespace samba
