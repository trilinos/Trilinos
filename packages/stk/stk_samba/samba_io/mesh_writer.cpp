#include <Ioss_SubSystem.h>
#include <init/Ionit_Initializer.h>

#include <samba_io/mesh_writer.hpp>
#include <samba_io/ioss_topology_map.hpp>
#include <samba_io/sidesets_manager.hpp>

namespace samba {
namespace io {

mesh_writer:: mesh_writer(MPI_Comm comm,
                          const std::string& file_name,
                          mesh mesh_arg,
                          ioss_mapper::Ptr ioss_mapper_arg)
  : m_io_region(NULL)
  , m_property_manager(NULL)
  , m_mesh(mesh_arg)
  , m_ioss_mapper_ptr(ioss_mapper_arg)
  , m_next_element_local_id(0)
{
  if (!m_ioss_mapper_ptr) {
    m_ioss_mapper_ptr = ioss_mapper::Ptr(new ioss_mapper(m_mesh));
  }

  m_property_manager = new Ioss::PropertyManager;
  Ioss::DatabaseIO *dbo = Ioss::IOFactory::create("exodusII", file_name, Ioss::WRITE_RESULTS,
                                                  comm, *m_property_manager);
  if (dbo == NULL || !dbo->ok()) {
    std::ostringstream oss;
    oss << "ERROR: Could not open results exodusII database '" << file_name;
    throw std::runtime_error(oss.str());
  }

  m_io_region = new Ioss::Region(dbo, "results_output");

  define_output_db(*m_io_region);
  write_output_db(*m_io_region);
}

mesh_writer::~mesh_writer()
{
  delete m_io_region;
  delete m_property_manager;
}

void mesh_writer::define_output_db(Ioss::Region & io_region)
{
  ioss_mapper &mapper = *m_ioss_mapper_ptr;

  io_region.begin_mode( Ioss::STATE_DEFINE_MODEL );

  m_io_region->property_add(Ioss::Property("spatial_dimension", m_mesh.spatial_dimension()() ));

  const std::vector<entity_block_key> &node_blocks = mapper.m_node_blocks;
  if (node_blocks.size() > 1)
  {
    throw std::runtime_error("ERROR:  samba_io::mesh_writer found more than one node block");
  }
  for (size_t nbi = 0; nbi < node_blocks.size(); ++nbi)
  {
    define_node_block(node_blocks[nbi], io_region);
  }

  std::vector<entity_block_key> element_blocks;
  m_mesh.get_entity_blocks(entity_rank::element(), element_blocks);
  for (size_t ebi = 0; ebi < element_blocks.size(); ++ebi)
  {
    define_element_block(element_blocks[ebi], io_region);
  }

  BOOST_ASSERT_MSG(m_mesh.spatial_dimension()() == 3,
                   "samba::io::mesh_writer only handles 3D meshes");
  std::vector<entity_block_key> face_sets;
  m_mesh.get_entity_blocks(entity_rank::face(), face_sets);
  for (size_t fsi = 0; fsi < face_sets.size(); ++fsi)
  {
    define_sideset(face_sets[fsi], io_region);
  }

  std::vector<entity_block_key> nodesets;
  m_mesh.get_entity_blocks(entity_rank::node(), nodesets);
  for (size_t nsi = 0; nsi < nodesets.size(); ++nsi)
  {
    define_nodeset(nodesets[nsi], io_region);
  }

  io_region.end_mode( Ioss::STATE_DEFINE_MODEL );
}

void mesh_writer::define_node_block(entity_block_key node_block_key, Ioss::Region& io_region)
{
  std::vector<entity_block_key> &node_blocks = m_ioss_mapper_ptr->m_node_blocks;
  if (std::find(node_blocks.begin(), node_blocks.end(), node_block_key) == node_blocks.end())
  {
    // Got a node set confused with a node block.  This problem will go away when samba
    // differentiates between sets and blocks.
    return;
  }

  std::string block_name = m_mesh.get_name(node_block_key);
  size_t num_nodes = m_mesh.num_entities(node_block_key);

  Ioss::NodeBlock* nb =
    new Ioss::NodeBlock(io_region.get_database(), block_name, num_nodes,
                        m_mesh.spatial_dimension()() );
  io_region.add(nb);
}

void mesh_writer::define_sideset(entity_block_key sideset_key, Ioss::Region &io_region)
{
  std::string sideset_name = m_mesh.get_name(sideset_key);
  Ioss::SideSet* ss = new Ioss::SideSet(io_region.get_database(), sideset_name);
  io_region.add(ss);

  size_t num_sides = m_mesh.num_entities(sideset_key & entity_rank::face());
  Ioss::SideBlock* sb = new Ioss::SideBlock(ss->get_database(), sideset_name,
                                            "unknown","unknown", num_sides );
  ss->add(sb);
}

void mesh_writer::define_element_block(entity_block_key elt_block_key, Ioss::Region &io_region)
{
  std::vector<partition_id> partitions;
  m_mesh.get_partitions(set_expression(elt_block_key) & set_expression(entity_rank::element()), partitions);
  partition_proxy some_partition = m_mesh[partitions.front()];
  entity_topology topology = some_partition.topology();

  // TO DO: CHECK THAT THE TOPOLOGY IS THE SAME FOR ALL THE PARTITIONS

  std::string block_name = m_mesh.get_name(elt_block_key);
  std::string topology_name = map_topology_to_ioss(topology);
  size_t num_elts = m_mesh.num_entities(elt_block_key & entity_rank::element());

  Ioss::ElementBlock* eb =
    new Ioss::ElementBlock(io_region.get_database(),block_name, topology_name, num_elts);
  io_region.add(eb);
}

void mesh_writer::define_nodeset(entity_block_key nodeset_key, Ioss::Region &io_region)
{
  std::vector<entity_block_key> &node_blocks = m_ioss_mapper_ptr->m_node_blocks;
  if (std::find(node_blocks.begin(), node_blocks.end(), nodeset_key) != node_blocks.end())
  {
    // Got a node set confused with a node block.  This problem will go away when samba
    // differentiates between sets and blocks.
    return;
  }

  std::string nodeset_name = m_mesh.get_name(nodeset_key);

  Ioss::NodeSet* ns =
    new Ioss::NodeSet(io_region.get_database(), m_mesh.get_name(nodeset_key),
                      m_mesh.num_entities(nodeset_key));
  io_region.add(ns);
}

void mesh_writer::write_output_db(Ioss::Region& io_region)
{
  // The order in which elements are written determines local id in an
  // Exodus file.
  m_next_element_local_id = 0;

  ioss_mapper &mapper = *m_ioss_mapper_ptr;

  io_region.begin_mode(Ioss::STATE_MODEL);

  // Write out node block.
  const std::vector<entity_block_key> &node_blocks = mapper.m_node_blocks;
  BOOST_ASSERT_MSG(node_blocks.size() == 1, "mesh can only be output if it has one node block");
  Ioss::NodeBlock& nb = *io_region.get_node_blocks()[0]; //TODO hard-coded that only 1 node-block exists?
  write_node_block(node_blocks[0], nb);

  write_node_coordinates(nb);

  // Write out element blocks
  const Ioss::ElementBlockContainer& elem_blocks = io_region.get_element_blocks();
  for (Ioss::ElementBlockContainer::const_iterator it = elem_blocks.begin(),end = elem_blocks.end();
       it != end;
       ++it)
  {
    entity_block_key element_block = m_mesh.find_entity_block((*it)->name());
    write_element_block(element_block, **it);
  }

  // Write out sidesets
  const Ioss::SideSetContainer& side_sets = io_region.get_sidesets();
  for (Ioss::SideSetContainer::const_iterator it = side_sets.begin();
       it != side_sets.end(); ++it) {
    Ioss::SideSet* ss = *it;
    write_sideset(*ss);
  }

  // SHOULD QUERY THE io_region FOR THE nodesets, WHICH WERE ALL DEFINED ABOVE.
  std::vector<entity_block_key> nodesets;
  m_mesh.get_entity_blocks(entity_rank::node(), nodesets);
  for (size_t nsi = 0; nsi < nodesets.size(); ++nsi)
  {
    write_nodeset(nodesets[nsi], io_region);
  }

  io_region.end_mode(Ioss::STATE_MODEL);
}

void mesh_writer::write_node_block(entity_block_key node_block_key,  Ioss::NodeBlock& nb)
{
  ioss_mapper &mapper = *m_ioss_mapper_ptr;

  // Because there is only one node block, the local id for a node will be the same
  // as its index into its node block
  int new_local_id;
  std::map<entity_block_key, int>::iterator new_idx_probe =
    mapper.m_new_idx_in_block.find(node_block_key);
  if (new_idx_probe != mapper.m_new_idx_in_block.end())
  {
    new_local_id = new_idx_probe->second;
  }
  else
  {
    new_local_id = 0;
  }

  entity_key_vector nodes;
  m_mesh.get_entities(node_block_key, nodes);
  size_t num_nodes = nodes.size();
  std::vector<int> node_ids(num_nodes);

  entity_key_to_ioss_id_field to_idx_in_block = mapper.m_to_ioss_idx_in_block;
  entity_key_to_ioss_id_field to_local_id = mapper.m_to_ioss_local_id;
  entity_key_to_ioss_id_field to_global_id = mapper.m_to_ioss_global_id;

  for (size_t i = 0; i < num_nodes; ++i)
  {
    entity_key node_key = nodes[i];

    // Write order of elements determines local id in an Exodus file.
    int local_id = to_local_id[node_key];
    if (local_id < 0)
    {
      local_id = new_local_id++;
      to_idx_in_block[node_key] = to_local_id[node_key] = local_id;
    }

    node_ids[local_id] = to_global_id[node_key] + 1;
  }

  // Be ready for new nodes to be added.
  mapper.m_new_idx_in_block[node_block_key] = new_local_id;

  size_t num_ids_written = nb.put_field_data("ids", node_ids);
  if (node_ids.size() != num_ids_written)
  {
    throw std::runtime_error("samba_io::mesh_writer FAILED in Ioss::NodeBlock::put_field_data");
  }
}

void mesh_writer::write_element_block(entity_block_key element_block_key,  Ioss::ElementBlock &eb)
{
  ioss_mapper &mapper = *m_ioss_mapper_ptr;

  int new_idx_in_block;
  std::map<entity_block_key, int>::iterator new_idx_probe =
    mapper.m_new_idx_in_block.find(element_block_key);
  if (new_idx_probe != mapper.m_new_idx_in_block.end())
  {
    new_idx_in_block = new_idx_probe->second;
  }
  else
  {
    new_idx_in_block = 0;
  }

  entity_key_vector elements;
  m_mesh.get_entities(element_block_key & entity_rank::element(), elements);
  size_t num_elements = elements.size();
  std::vector<int> elem_ids(num_elements);
  entity_key_vector elements_out(num_elements);

  //
  // Write out the elements in the block.
  //
  entity_key_to_ioss_id_field to_idx_in_block = mapper.m_to_ioss_idx_in_block;
  entity_key_to_ioss_id_field to_local_id = mapper.m_to_ioss_local_id;
  entity_key_to_ioss_id_field to_global_id = mapper.m_to_ioss_global_id;
  for (size_t i = 0; i < num_elements; ++i)
  {
    entity_key elem_key = elements[i];

    // Write order of elements determines local id in an Exodus file.
    to_local_id[elem_key] = m_next_element_local_id++;

    int idx = to_idx_in_block[elem_key];
    // If the element has not been written before, we need to give it a new
    // index in its block.
    if (idx < 0)
    {
      idx = new_idx_in_block++;
      to_idx_in_block[elem_key] = idx;
    }

    // Assume that the global exodus id has been set!
    elem_ids[idx] = to_global_id[elem_key] + 1;

    elements_out[idx] = elem_key;
  }

  // Be ready for new elements to be added.
  mapper.m_new_idx_in_block[element_block_key] = new_idx_in_block;

  size_t num_ids_written = eb.put_field_data("ids", elem_ids);

  if (elem_ids.size() != num_ids_written)
  {
    throw std::runtime_error("samba_io::mesh_writer FAILED in Ioss::ElementBlock::put_field_data");
  }

  // Now write out the connectivity for the block.  This must match the output order of the
  // elements.

  std::vector<int> connectivity;
  for (size_t i = 0; i < num_elements; ++i)
  {
    entity_proxy elem_proxy = m_mesh[elements_out[i]];
    partition_index_iterator nodes_iter, nodes_end;
    for (nodes_iter = elem_proxy.begin_nodes(), nodes_end = elem_proxy.end_nodes();
         nodes_iter != nodes_end;
         ++nodes_iter)
    {
      entity_key node_key = m_mesh.convert(*nodes_iter);
      connectivity.push_back(to_local_id[node_key] + 1);
    }
  }

  eb.put_field_data("connectivity_raw", connectivity);
}

void mesh_writer::write_sideset(Ioss::SideSet &io_sideset)
{
  entity_block_key sideset_key = m_mesh.find_entity_block(io_sideset.name());
  entity_key_vector side_keys;
  m_mesh.get_entities(sideset_key & entity_rank::face(), side_keys);

  entity_key_to_ioss_id_field to_local_id = m_ioss_mapper_ptr->m_to_ioss_local_id;

  std::vector<int> raw_side_data;
  size_t num_sides = side_keys.size();
  for (size_t i = 0; i < num_sides; ++i)
  {
    entity_key side_key = side_keys[i];
    entity_proxy side_proxy = m_mesh[side_key];
    BOOST_ASSERT_MSG(side_proxy.num_connectivity(entity_rank::element()) == 1,
                     (debug_message() << "Expected one element back-relation for " << side_key << ", found: " <<
                      side_proxy.num_connectivity(entity_rank::element())));
    entity_key_iterator element_connectivity_it = side_proxy.begin_elements<entity_key>();
    ordinal_iterator element_ord_it = side_proxy.begin_elements<connectivity_ordinal>();
    entity_key element_key = *element_connectivity_it;
    connectivity_ordinal side_ordinal = *element_ord_it;

    int local_element_id = to_local_id[element_key];
    raw_side_data.push_back(local_element_id + 1);
    raw_side_data.push_back(side_ordinal() + 1);
  }

  //we know there is only 1 side-block (we set up the sideset in the above define_output_db method)
  Ioss::SideBlock* sb = io_sideset.get_side_blocks()[0];

  const size_t num_sides_written = sb->put_field_data("element_side_raw",raw_side_data);
  if (num_sides_written != raw_side_data.size()/2) {
    throw std::runtime_error("samba::io::mesh_writer ERROR, failed in sideset output.");
  }
}

void mesh_writer::write_nodeset(entity_block_key nodeset_key, Ioss::Region &io_region)
{

}

void mesh_writer::write_node_coordinates( Ioss::NodeBlock& nb)
{
  ioss_mapper &mapper = *m_ioss_mapper_ptr;

  coordinates_field_type node_coordinates = mapper.m_node_coordinates;

  entity_key_vector nodes;
  m_mesh.get_entities(entity_rank::node(), nodes);
  size_t num_nodes = nodes.size();

  std::vector<double> coordinates_out(3 * num_nodes, 0);
  entity_key_to_ioss_id_field to_local_id = mapper.m_to_ioss_local_id;

  // Coordinates data to be written out must be consistent with the ordering of the nodes output
  // by write_node_block.
  size_t dimension = m_mesh.spatial_dimension()();
  for (size_t i = 0; i < num_nodes; ++i)
  {
    entity_key node_key = nodes[i];
    int local_id = to_local_id[node_key];
    double *xyz_begin = &coordinates_out[dimension * local_id];
    std::copy(node_coordinates[node_key], node_coordinates[node_key] + dimension, xyz_begin);
  }

  size_t num_coords_written = nb.put_field_data("mesh_model_coordinates",
                                                const_cast<std::vector<double>&>(coordinates_out));

  if (num_nodes != num_coords_written)
  {
    throw std::runtime_error("FAILED in Ioss::NodeBlock::put_field_data for coordinates.");
  }
}

} // namespace io
} // namespace samba
